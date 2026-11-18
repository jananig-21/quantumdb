"""Do better optimizers help under the new "qaoa" formulation? Pilot on the team's saved SPIQ points.

Energy: H = lambda*(A+B) + C, written out from docs/problem-formulation.md, checked at start-up against the team's
saved exact ground state, its plan cost, and the SPIQ point's stim energy (which also checks the bit order).
Starts are read from saved files; the SPIQ search is not run. Every evaluation is 10,240 shots, scored exactly at
the best-seen point, as in ../many_problems/run_spiq_many.py.

Usage (QAOA env): python pilot.py <formulation worktree> <problem[,problem]> <optimizers> <trials> <budget> <workers>
"""
import itertools, json, os, sys, time, zlib
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.join(HERE, "..")
sys.path[:0] = [os.path.join(STUDY, "many_problems"), STUDY, os.path.join(STUDY, "..", "base")]
import Scripts.Postprocessing as PP
from fastsim import compile_circuit
from fastqaoa import qaoa_probs, sample_energy
from relaxed_sim import compile_layers, probs as relaxed_probs
from optimizers import make_optimizer
from qiskit.algorithms.optimizers import BOBYQA
from qiskit.utils import algorithm_globals
try:
    from qiskit import qpy
except ImportError:
    from qiskit.circuit import qpy_serialization as qpy

SHOTS, ITER = 10240, 10000
SI = "spiq_init_outputs"
FILES = {
    "P1": {"single": f"{SI}/single_Week105/spiq_initial_point_input0_reps2_3rel_2pred_qaoa.json", "multi": f"{SI}/diverse_Week102/spiq_initial_point_input0_reps2_3rel_2pred_qaoa.json"},
    "P2": {"single": f"{SI}/single_Week106/spiq_initial_point_input1_reps2_4rel_2pred_qaoa.json", "multi": f"{SI}/diverse_Week103/spiq_initial_point_input1_reps2_4rel_2pred_qaoa.json"},
    "P3": {"single": f"{SI}/single_Week107/spiq_initial_point_input2_reps2_4rel_6pred_qaoa.json", "multi": f"{SI}/diverse_Week104/spiq_initial_point_input2_reps2_4rel_6pred_qaoa.json"},
    "P4": {"single": f"{SI}/single_Week108/spiq_initial_point_input3_reps2_5rel_5pred_qaoa.json"},
    "P5": {"single": f"{SI}/single_Week109/spiq_initial_point_input4_reps2_6rel_7pred_qaoa.json"},
}
_C = {}


def plan_cost(o, card, pred, sel):
    """Sum of intermediate sizes, final join excluded; every predicate applies, repeated pairs included."""
    total, size = 0.0, float(card[o[0]])
    for j in range(1, len(o) - 1):
        size *= card[o[j]]
        for (a, b), s in zip(pred, sel):
            if (a == o[j] and b in o[:j]) or (b == o[j] and a in o[:j]):
                size *= s
        total += size
    return total


def build(pi):
    card, pred, sel = pi["card"], [tuple(p) for p in pi["pred"]], pi["pred_sel"]
    R, J = len(card), len(card) - 2
    n = R * J
    lam, c, s = pi["penalty_weight"], np.log10(card), np.log10(sel)
    idx = np.arange(2 ** n, dtype=np.int64)
    V = lambda t, j: ((idx >> (t * J + j)) & 1).astype(np.int8)
    E = np.zeros(2 ** n)
    for j in range(J):
        vs = [V(t, j) for t in range(R)]
        E += lam * (j + 2 - sum(v.astype(np.int16) for v in vs)) ** 2
        for t in range(R):
            E += c[t] * vs[t]
            if j:
                E += lam * (V(t, j - 1) * (1 - vs[t]))
        for (a, b), sp in zip(pred, s):
            E += sp * (vs[a] * vs[b])
    # decode: identical score vectors give identical plans, so decode each distinct one once
    w = np.arange(1, R - 1)[R - 3::-1]
    cv = np.stack([sum(int(w[j]) * ((idx >> (t * J + j)) & 1) for j in range(J)) for t in range(R)], 1).astype(np.int16)
    base = int(w.sum()) + 1
    key = (cv * (base ** np.arange(R))).sum(1)
    u, first, inv = np.unique(key, return_index=True, return_inverse=True)
    cost = {o: plan_cost(o, card, pred, sel) for o in itertools.permutations(range(R))}
    raw_u, fb_u = np.empty(len(u)), np.empty(len(u))
    for i, f in enumerate(first):
        ro = PP.get_raw_join_order(cv[f].astype(float))
        raw_u[i] = cost[tuple(ro)]
        fb_u[i] = cost[tuple(int(x) for x in PP.postprocess_join_order(ro, cv[f].astype(float), R, pred))]
    return dict(n=n, E=E, gs=float(E.min()), raw=raw_u[inv], fb=fb_u[inv], opt=min(cost.values()))


def problem(name, repo):
    if name in _C:
        return _C[name]
    arms = {}
    P = None
    for arm, rel in FILES[name].items():
        p = os.path.join(repo, rel)
        j = json.load(open(p))
        pc = qpy.load(open(os.path.join(os.path.dirname(p), os.path.basename(j["pcirc_qpy"])), "rb"))[0]
        P = P or build(j["problem_input"])
        n, ops = compile_circuit(pc)
        comp = compile_layers(n, [[a, b, None if c is None else list(c)] for a, b, c in ops])
        by = dict(zip(j["relaxed_param_names"], j["relaxed_initial_point"]))
        x0 = [by[q.name] for q in pc.parameters]
        e0 = float(relaxed_probs(comp, np.array(x0)) @ P["E"])
        assert abs(e0 - j["energy_best_qubo"]) < 1e-6, (name, arm, e0, j["energy_best_qubo"])
        gs = (j.get("exact_ground_state") or {}).get("exact_qubo_energy")
        assert gs is None or abs(P["gs"] - gs) < 1e-6, (name, P["gs"], gs)
        arms["spiq_" + arm] = (comp, x0)
    arms["random"] = (None, [0.5] * 4)
    offset = P["E"].mean()               # the mean is the identity term of the Ising form; a global phase in QAOA
    P["diag"] = P["E"] - offset
    _C[name] = (P, arms)
    return _C[name]


def metrics(P, pr):
    isopt = lambda c: np.isclose(c, P["opt"], rtol=1e-9)
    return dict(E=float(pr @ P["E"]), P_fb=float(pr @ isopt(P["fb"])), P_raw=float(pr @ isopt(P["raw"])),
                P_ground=float(pr[np.isclose(P["E"], P["gs"], rtol=0, atol=1e-9)].sum()))


def run(t):
    repo, name, arm, opt_name, trial, budget = t
    P, arms = problem(name, repo)
    comp, x0 = arms[arm.split("+k")[0]]
    state = (lambda x: relaxed_probs(comp, x)) if comp is not None else (lambda x: qaoa_probs(P["diag"], x, 2))
    seed = zlib.crc32(f"pilot|{name}|{arm}|{opt_name}|{trial}".encode())
    algorithm_globals.random_seed = seed
    np.random.seed(seed % 2 ** 32)
    rng = np.random.default_rng(seed)
    npar = len(x0)
    kick = float(arm.split("+k")[1]) if "+k" in arm else 0.0
    start = np.array(x0, float) + kick * rng.standard_normal(npar)      # "+k<sigma>": random kick off the start point
    st = dict(n=0, f=np.inf, x=start.copy(), first_low=0)

    def one(x):
        if st["n"] >= budget:
            return st["f"]
        st["n"] += 1
        f = sample_energy(state(x), P["E"], SHOTS, rng)
        if f < st["f"]:
            st.update(f=f, x=np.array(x, float), first_low=st["n"])
        return f

    def fun(x):
        x = np.asarray(x, float).reshape(-1, npar)
        v = [one(r) for r in x]
        return v[0] if len(v) == 1 else np.array(v)

    opt = make_optimizer(opt_name, ITER, budget)
    bounds = [(v - np.pi, v + np.pi) for v in x0] if isinstance(opt, BOBYQA) else None
    t0 = time.time()
    opt.minimize(fun=fun, x0=start, bounds=bounds)
    s, e = metrics(P, state(np.array(x0, float))), metrics(P, state(st["x"]))
    return dict(problem=name, arm=arm, optimizer=opt_name, trial=trial, budget=budget, evals=st["n"], best_eval=st["first_low"],
                n_angles=npar, E_ground=P["gs"], opt_cost=P["opt"], E_start=s["E"], P_fb_start=s["P_fb"], P_raw_start=s["P_raw"],
                E_logged_best=st["f"], E_end=e["E"], P_fb_end=e["P_fb"], P_raw_end=e["P_raw"], P_ground_end=e["P_ground"], wall_s=time.time() - t0)


if __name__ == "__main__":
    repo, names, opts, trials, budget, workers = sys.argv[1], sys.argv[2].split(","), sys.argv[3].split(","), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    T0 = int(os.environ.get("TRIAL_OFFSET", 0))      # different trial numbers give different seeds
    kicks = [k for k in os.environ.get("KICKS", "").split(",") if k]
    tasks = [(repo, n, a + s, o, t, budget) for n in names for a in ["spiq_" + k for k in FILES[n]] + ["random"]
             for s in ([""] if not kicks else [f"+k{k}" for k in kicks]) if not (s and a == "random") for o in opts for t in range(T0 + 1, T0 + trials + 1)]
    out = os.path.join(HERE, f"runs_{'_'.join(names)}_b{budget}{'_kicks' if kicks else ''}_t{T0}.jsonl")
    with Pool(workers) as pool, open(out, "a") as f:
        for r in pool.imap_unordered(run, tasks):
            f.write(json.dumps(r, default=float) + "\n")
            f.flush()
            print(r["problem"], r["arm"], r["optimizer"], r["trial"], f"E {r['E_start']:.3f}->{r['E_end']:.3f} ground {r['E_ground']:.3f} fb {r['P_fb_start']:.3f}->{r['P_fb_end']:.3f} raw {r['P_raw_end']:.3f} evals {r['evals']} {r['wall_s']:.0f}s", flush=True)
