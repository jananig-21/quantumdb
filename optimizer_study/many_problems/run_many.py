"""Large-scale arm: all 10 optimizers on every problem in problems.jsonl (random-start QAOA).

Per problem, exactly as in the P1-P3 study's random-initialization arm:
  * QUBO from the team's QUBOGenerator.generate_IBMQ_QUBO_for_left_deep_trees_v2 (unchanged);
  * QAOA depth p = 2, the team's fixed start [0.5] * 4;
  * every energy evaluation is a 10,240-shot estimate on the QUBO scale (fresh shot noise per call);
  * the optimizers from ../optimizers.py (make_optimizer: COBYLA/AQGD/SPSA built as in the team's
    get_optimizer, the rest library defaults + the documented FD step / Adam lr), budget 1,000
    evaluations; at the budget the best point seen is returned, like ../optimizers.BudgetGuard;
  * the returned point is scored EXACTLY (statevector) and decoded with the team's
    Postprocessing decoder (raw + fallback) into plan-quality metrics.

What differs from run_study.py, and why: running 10,880 optimizations through qiskit's VQE + Aer
path would take months on one machine (35 s per evaluation on 20 qubits). Instead the circuit is
simulated with fastqaoa.qaoa_probs (== ../fastsim.qaoa_probs to 1e-10, which is verified == qiskit
QAOAAnsatz statevector to 1e-10) and the 10,240 shots are drawn from that exact distribution
(fastqaoa.sample_energy) - the same distribution a noiseless QasmSimulator samples from. The optimizer, start, ansatz, QUBO, decoder and metrics are
unchanged. verify_against_driver.py checks this objective against the team's driver on P1-P3.

Usage:  python run_many.py [--workers 9] [--budget 1000] [--only ID ...]
Writes results/<problem id>.json (resumable: finished problems are skipped).
"""
import argparse, contextlib, io, itertools, json, os, sys, time, zlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.abspath(os.path.join(HERE, ".."))
BASE = os.path.abspath(os.path.join(STUDY, "..", "base"))
sys.path[:0] = [HERE, STUDY, BASE]

OPTIMIZERS = ["COBYLA", "AQGD", "SPSA", "NELDER_MEAD", "POWELL", "BOBYQA", "L_BFGS_B", "SLSQP", "ADAM", "NFT"]
REPS, SHOTS, ITERATIONS, X0 = 2, 10240, 10000, [0.5] * 4
CHECKPOINTS = [1, 10, 25, 50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000]


# ---- the problem, through the team's code ---------------------------------------------------
def build_problem(p):
    import Scripts.QUBOGenerator as QUBOGenerator
    import Scripts.Postprocessing as Postprocessing
    card, pred, pred_sel = p["card"], [tuple(x) for x in p["pred"]], p["pred_sel"]
    with contextlib.redirect_stdout(io.StringIO()):
        qubo, penalty = QUBOGenerator.generate_IBMQ_QUBO_for_left_deep_trees_v2(card, pred, pred_sel)
    n = qubo.get_num_vars()
    lin = qubo.objective.linear.to_array()
    quad = qubo.objective.quadratic.to_array(symmetric=False)
    idx = np.arange(2 ** n)
    # energy of every bitstring, term by term (an n x 2^n bit matrix would need ~170 MB per worker at 20 qubits)
    bit = lambda k: ((idx >> k) & 1).astype(np.float64)              # bit k of index = variable k
    E = np.full(2 ** n, float(qubo.objective.constant))
    for i in np.flatnonzero(lin):
        E += lin[i] * bit(i)
    for i, j in zip(*np.nonzero(quad)):
        E += quad[i, j] * (((idx >> i) & (idx >> j) & 1).astype(np.float64))
    for k in np.random.default_rng(0).integers(0, len(E), 20):      # vs qiskit's own evaluator
        xk = [int((k >> v) & 1) for v in range(n)]
        assert abs(E[k] - qubo.objective.evaluate(xk)) < 1e-6 * max(1, abs(E[k]))
    _, offset = qubo.to_ising()
    # decode: the plan depends only on the R*(R-2) v-variables, so decode each prefix once
    R = len(card)
    nv = R * (R - 2)
    w = np.arange(1, R - 1)[R - 3::-1]
    plan_cost, cd = {}, {}
    def cost(o):
        if tuple(o) not in plan_cost:
            plan_cost[tuple(o)] = Postprocessing.get_costs_for_leftdeep_tree(list(o), card, pred, pred_sel, cd)
        return plan_cost[tuple(o)]
    raw_tab, fb_tab = np.empty(2 ** nv), np.empty(2 ** nv)
    for v in range(2 ** nv):
        bits = np.array([(v >> k) & 1 for k in range(nv)])
        cv = bits.reshape(R, R - 2).dot(w)
        raw = Postprocessing.get_raw_join_order(cv)
        fb = Postprocessing.postprocess_join_order(raw, cv, R, pred)
        raw_tab[v], fb_tab[v] = cost(raw), cost(fb)
    prefix = idx & (2 ** nv - 1)
    all_costs = {o: cost(list(o)) for o in itertools.permutations(range(R))}
    opt = min(all_costs.values())
    distinct = sorted(set(np.round(list(all_costs.values()), 6)))
    gs = float(E.min())
    return dict(n=n, E=E, diag=E - float(offset), raw_cost=raw_tab[prefix], fb_cost=fb_tab[prefix],
                opt=opt, distinct=distinct, gs=gs, penalty=penalty,
                E_uniform=float(E.mean()), worst_cost=max(all_costs.values()))


def exact_metrics(P, theta):
    from fastqaoa import qaoa_probs
    pr = qaoa_probs(P["diag"], np.asarray(theta, float), REPS)
    second = P["distinct"][1] if len(P["distinct"]) > 1 else np.inf
    isopt = lambda c: np.isclose(c, P["opt"], rtol=1e-9)
    return dict(E_exact=float(pr @ P["E"]),
                P_opt_fb=float(pr @ isopt(P["fb_cost"])), P_opt_raw=float(pr @ isopt(P["raw_cost"])),
                P_top2_fb=float(pr @ (P["fb_cost"] <= second * (1 + 1e-9))),
                mean_cost_fb=float(pr @ P["fb_cost"]),
                P_ground=float(pr[np.isclose(P["E"], P["gs"], rtol=0, atol=1e-9)].sum()))


# ---- one optimization ---------------------------------------------------------------------------
class _Budget(Exception):
    pass


def run_one(P, name, budget, seed):
    from qiskit.utils import algorithm_globals
    from qiskit.algorithms.optimizers import BOBYQA
    from fastqaoa import qaoa_probs, sample_energy
    from optimizers import make_optimizer
    algorithm_globals.random_seed = seed            # SPSA / NFT internal randomness
    np.random.seed(seed % 2**32)
    rng = np.random.default_rng(seed)               # shot noise
    st = dict(n=0, best_f=np.inf, best_x=np.array(X0, float), log=[], inc=[])

    def one(x):
        if st["n"] >= budget:
            raise _Budget()
        st["n"] += 1
        f = sample_energy(qaoa_probs(P["diag"], x, REPS), P["E"], SHOTS, rng)
        if f < st["best_f"]:
            st["best_f"], st["best_x"] = f, np.array(x, float)
        st["log"].append(f)
        if st["n"] in CHECKPOINTS:
            st["inc"].append((st["n"], st["best_x"].tolist()))
        return f

    def fun(x):                                      # AQGD passes a flattened batch of points
        x = np.asarray(x, float).reshape(-1, 4)
        vals = [one(r) for r in x]
        return vals[0] if len(vals) == 1 else np.array(vals)

    opt = make_optimizer(name, ITERATIONS, budget)
    bounds = [(x - np.pi, x + np.pi) for x in X0] if isinstance(opt, BOBYQA) else None
    t0 = time.time()
    try:
        res = opt.minimize(fun=fun, x0=np.array(X0, float), bounds=bounds)
        x_ret, hit = np.asarray(res.x, float), False
    except _Budget:
        x_ret, hit = st["best_x"], True
    wall = time.time() - t0
    out = dict(optimizer=name, seed=seed, evaluations=st["n"], budget_hit=hit, wall_s=wall,
               x_returned=x_ret.tolist(), E_logged_best=st["best_f"], E_logged_first=st["log"][0],
               **exact_metrics(P, x_ret))
    # anytime curve: exact energy of the best-seen point after k evaluations (k = checkpoints)
    out["anytime"] = [(k, exact_metrics(P, x)["E_exact"]) for k, x in st["inc"]]
    return out


def solve_problem(args):
    p, budget, outdir = args
    path = os.path.join(outdir, p["id"] + ".json")
    if os.path.exists(path):
        return p["id"], "skip"
    t0 = time.time()
    P = build_problem(p)
    start = exact_metrics(P, X0)
    runs = []
    for name in OPTIMIZERS:
        seed = zlib.crc32(f"{p['id']}|{name}".encode())
        r = run_one(P, name, budget, seed)
        r["gap_closed"] = (start["E_exact"] - r["E_exact"]) / (start["E_exact"] - P["gs"])
        runs.append(r)
        print(f"  [progress] {p['id']} {P['n']}q {name}: {r['evaluations']} evals in {r['wall_s']:.0f}s "
              f"({r['wall_s'] / max(1, r['evaluations']):.3f} s/eval)", flush=True)
    rec = dict(problem=p, num_qubits=P["n"], ground_state_energy=P["gs"], E_uniform=P["E_uniform"],
               optimal_cost=P["opt"], worst_cost=P["worst_cost"], num_distinct_costs=len(P["distinct"]),
               penalty_weight=P["penalty"], start=start, budget=budget, runs=runs,
               wall_s=time.time() - t0)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(rec, f, default=float)
    os.replace(tmp, path)
    return p["id"], f"done {P['n']}q {time.time() - t0:.0f}s"


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 1))
    ap.add_argument("--budget", type=int, default=1000)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--source", choices=["all", "team", "generated"], default="all",
                    help="team = the 88 queries committed in the team's project; generated = the 1,000 new ones")
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    probs = [json.loads(l) for l in open(os.path.join(HERE, "problems.jsonl"))]
    if a.only:
        probs = [p for p in probs if p["id"] in a.only]
    if a.source == "team":
        probs = [p for p in probs if p["source"] != "generated"]
    elif a.source == "generated":
        probs = [p for p in probs if p["source"] == "generated"]
    probs.sort(key=lambda p: -p["num_qubits"])      # biggest first for load balance
    jobs = [(p, a.budget, a.out) for p in probs]
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    from multiprocessing import Pool
    done = 0
    with Pool(a.workers, maxtasksperchild=20) as pool:
        for pid, msg in pool.imap_unordered(solve_problem, jobs):
            done += 1
            print(f"[{done}/{len(jobs)}] {pid}: {msg}", flush=True)
