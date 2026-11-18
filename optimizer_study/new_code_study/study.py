"""All 10 optimizers on the team's new code: Divya's "qaoa" formulation (P1-P5) and seeded multi-start SPIQ, Fangyi's
relaxed fixed-start baseline (PR #2: per-gate circuit started at the exact lift of [0.5]*4), and the old "q2o"
formulation on P1-P3 for comparison.

Starts (where they exist):
  std4         stock QAOA, 2*reps = 4 shared angles, all 0.5 (the team's old baseline)
  lift         per-gate circuit started at the state of std4 (Fangyi's PR #2 baseline): cost angles coefficient x 0.5,
               mixer angles 1 - pi (see problem()); checked equal to std4's state to 1e-9
  spiq_single  the team's saved single-start SPIQ point
  spiq_multi   the team's saved seeded multi-start SPIQ point (P1-P3 only)

Checks at start-up (any failure stops the run): exact ground-state energy and optimal plan cost against the team's
saved values, every saved SPIQ point's energy against its stim energy, and the lifted start's state against std4's.
Scoring as in ../many_problems: 10,240 shots per evaluation, decode the best-seen point, score it exactly.

Usage (QAOA env): python study.py <problems e.g. qaoa:P1,q2o:P2> <budget> <trials> <workers> [--check]
"""
import itertools, json, os, sys, time, zlib
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.join(HERE, "..")
REPO_NEW = os.path.expanduser("~/Downloads/umich-db-quantumdb-wt/formulation")
sys.path[:0] = [os.path.join(STUDY, "many_problems"), STUDY, os.path.join(STUDY, "..", "base")]
import run_many as RM
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

SHOTS, ITER, REPS = 10240, 10000, 2
OPTIMIZERS = ["COBYLA", "AQGD", "SPSA", "NELDER_MEAD", "POWELL", "BOBYQA", "L_BFGS_B", "SLSQP", "ADAM", "NFT"]
CHECKPOINTS = (50, 100, 250, 500, 1000, 2500, 5000, 10000)
SI = os.path.join(REPO_NEW, "spiq_init_outputs")
# saved SPIQ files; the q2o single-start files are the team's April-June artifacts the earlier study used
SPIQ = {
    ("qaoa", "P1"): {"spiq_single": f"{SI}/single_Week105/spiq_initial_point_input0_reps2_3rel_2pred_qaoa.json",
                     "spiq_multi": f"{SI}/diverse_Week102/spiq_initial_point_input0_reps2_3rel_2pred_qaoa.json"},
    ("qaoa", "P2"): {"spiq_single": f"{SI}/single_Week106/spiq_initial_point_input1_reps2_4rel_2pred_qaoa.json",
                     "spiq_multi": f"{SI}/diverse_Week103/spiq_initial_point_input1_reps2_4rel_2pred_qaoa.json"},
    ("qaoa", "P3"): {"spiq_single": f"{SI}/single_Week107/spiq_initial_point_input2_reps2_4rel_6pred_qaoa.json",
                     "spiq_multi": f"{SI}/diverse_Week104/spiq_initial_point_input2_reps2_4rel_6pred_qaoa.json"},
    ("qaoa", "P4"): {"spiq_single": f"{SI}/single_Week108/spiq_initial_point_input3_reps2_5rel_5pred_qaoa.json"},
    ("qaoa", "P5"): {"spiq_single": f"{SI}/single_Week109/spiq_initial_point_input4_reps2_6rel_7pred_qaoa.json"},
    ("q2o", "P1"): {"spiq_single": f"{SI}/spiq_initial_point_input0_reps2_3rel_2pred.json",
                    "spiq_multi": f"{SI}/diverse_Week93/spiq_initial_point_input0_reps2_3rel_2pred.json"},
    ("q2o", "P2"): {"spiq_single": os.path.join(REPO_NEW, "new_spiq_outputs/spiq_initial_point_input0_reps2_4rel_2pred.json"),
                    "spiq_multi": f"{SI}/diverse_Week94/spiq_initial_point_input1_reps2_4rel_2pred.json"},
    ("q2o", "P3"): {"spiq_single": f"{SI}/spiq_initial_point_input0_reps2_4rel2.json",
                    "spiq_multi": f"{SI}/diverse_Week95/spiq_initial_point_input2_reps2_4rel_6pred.json"},
}
# best plan cost with every predicate applied (Divya's fixed checker), from docs/problem-formulation.md section 2
OPT_COST = {"P1": 15, "P2": 315, "P3": 29.4, "P4": 165, "P5": 1605}
INPUT_IDX = {"P1": 0, "P2": 1, "P3": 2, "P4": 3, "P5": 4}
_C = {}


def plan_cost(o, card, pred, sel):
    """Sum of intermediate result sizes, final join excluded; every predicate applies (repeated pairs too)."""
    total, size = 0.0, float(card[o[0]])
    for j in range(1, len(o) - 1):
        size *= card[o[j]]
        for (a, b), s in zip(pred, sel):
            if (a == o[j] and b in o[:j]) or (b == o[j] and a in o[:j]):
                size *= s
        total += size
    return total


def qaoa_energies(card, pred, sel, lam):
    """Every bitstring's energy, H = lam*(A+B) + C (docs/problem-formulation.md, 'Correctness'). Variable t*J+j = bit."""
    R, J = len(card), len(card) - 2
    c, s = np.log10(card), np.log10(sel)
    idx = np.arange(2 ** (R * J), dtype=np.int64)
    V = lambda t, j: ((idx >> (t * J + j)) & 1).astype(np.int8)
    E = np.zeros(len(idx))
    for j in range(J):
        vs = [V(t, j) for t in range(R)]
        E += lam * (j + 2 - sum(v.astype(np.int16) for v in vs)) ** 2
        for t in range(R):
            E += c[t] * vs[t]
            if j:
                E += lam * (V(t, j - 1) * (1 - vs[t]))
        for (a, b), sp in zip(pred, s):
            E += sp * (vs[a] * vs[b])
    return E


def decode_tables(card, pred, sel, n):
    """Raw and fallback plan cost of every bitstring; the plan depends only on the R*(R-2) table bits."""
    R, J = len(card), len(card) - 2
    nv = R * J
    idx = np.arange(2 ** nv, dtype=np.int64)
    w = np.arange(1, R - 1)[R - 3::-1]
    cv = np.stack([sum(int(w[j]) * ((idx >> (t * J + j)) & 1) for j in range(J)) for t in range(R)], 1).astype(np.int16)
    key = (cv.astype(np.int64) * ((int(w.sum()) + 1) ** np.arange(R))).sum(1)
    u, first, inv = np.unique(key, return_index=True, return_inverse=True)
    cost = {o: plan_cost(o, card, pred, sel) for o in itertools.permutations(range(R))}
    raw_u, fb_u = np.empty(len(u)), np.empty(len(u))
    for i, f in enumerate(first):
        ro = PP.get_raw_join_order(cv[f].astype(float))
        raw_u[i] = cost[tuple(ro)]
        fb_u[i] = cost[tuple(int(x) for x in PP.postprocess_join_order(ro, cv[f].astype(float), R, pred))]
    raw, fb = raw_u[inv], fb_u[inv]
    if n > nv:                                   # q2o: predicate bits follow the table bits and don't change the plan
        prefix = np.arange(2 ** n, dtype=np.int64) & (2 ** nv - 1)
        raw, fb = raw[prefix], fb[prefix]
    return raw, fb, min(cost.values())


def load_pcirc(jpath):
    j = json.load(open(jpath))
    b = os.path.basename(j["pcirc_qpy"])
    for d in (os.path.dirname(jpath), os.path.join(REPO_NEW, "new_spiq_outputs")):
        if os.path.exists(os.path.join(d, b)):
            return j, qpy.load(open(os.path.join(d, b), "rb"))[0]
    raise FileNotFoundError(b)


def problem(key):
    if key in _C:
        return _C[key]
    form, name = key
    files = SPIQ[key]
    j0, _ = load_pcirc(files["spiq_single"])
    if form == "qaoa":
        pi = j0["problem_input"]
        card, pred, sel = pi["card"], [tuple(p) for p in pi["pred"]], pi["pred_sel"]
        assert pi["input_idx"] == INPUT_IDX[name] and pi.get("formulation") == "qaoa", (key, pi)
        E = qaoa_energies(card, pred, sel, pi["penalty_weight"])
        gs_ref = j0["exact_ground_state"]["exact_qubo_energy"]
    else:
        g = json.load(open(os.path.join(STUDY, "ground_truth.json")))[name]
        card, pred, sel = g["card"], [tuple(p) for p in g["pred"]], g["pred_sel"]
        E = RM.build_problem(dict(card=card, pred=pred, pred_sel=sel))["E"]      # the team's q2o QUBO, as in the study
        gs_ref = g["ground_state_energy"]
    n = int(np.log2(len(E)))
    raw, fb, opt = decode_tables(card, pred, sel, n)
    assert abs(E.min() - gs_ref) < 1e-6, (key, E.min(), gs_ref)
    assert np.isclose(opt, OPT_COST[name], rtol=2e-3), (key, opt)
    P = dict(key=key, n=n, E=E, diag=E - E.mean(), gs=float(E.min()), raw=raw, fb=fb, opt=opt, starts={})
    std_pr = qaoa_probs(P["diag"], np.full(2 * REPS, 0.5), REPS)
    P["starts"]["std4"] = (None, [0.5] * (2 * REPS))
    for arm, f in files.items():
        j, pc = load_pcirc(f)
        nq, ops = compile_circuit(pc)
        assert nq == n, (key, arm, nq, n)
        comp = compile_layers(nq, [[a, b, None if c is None else list(c)] for a, b, c in ops])
        by = dict(zip(j["relaxed_param_names"], j["relaxed_initial_point"]))
        x0 = [by[p.name] for p in pc.parameters]
        e0 = float(relaxed_probs(comp, np.array(x0)) @ E)
        ref = j.get("energy_best_qubo")
        if ref is None:                                # older artifacts store only the Ising energy
            ref = j["energy_best"] + (E.mean() if form == "q2o" else 0)
        assert abs(e0 - ref) < 1e-6 * max(1, abs(ref)), (key, arm, e0, ref)
        P["starts"][arm] = (comp, x0)
        if "lift" not in P["starts"]:
            # The saved circuits implement each mixer as s, sx, rz(theta), sx, rz(5pi/2) = RX(theta + pi), so the lift of
            # beta = 0.5 is 2*0.5 - pi on mixer angles; cost angles are coefficient x 0.5. Checked against std4 below.
            lift = [float(p.name.split("*")[0]) * 0.5 - (np.pi if "beta" in p.name else 0.0) for p in pc.parameters]
            dev = float(np.abs(relaxed_probs(comp, np.array(lift)) - std_pr).max())
            assert dev < 1e-9, (key, "lift differs from std4", dev)
            P["starts"]["lift"] = (comp, lift)
            P["lift_dev"] = dev
    _C[key] = P
    return P


def metrics(P, pr):
    isopt = lambda c: np.isclose(c, P["opt"], rtol=1e-9)
    return dict(E=float(pr @ P["E"]), P_fb=float(pr @ isopt(P["fb"])), P_raw=float(pr @ isopt(P["raw"])),
                P_ground=float(pr[np.isclose(P["E"], P["gs"], rtol=0, atol=1e-9)].sum()))


def run(t):
    key, arm, opt_name, trial, budget = t
    P = problem(key)
    comp, x0 = P["starts"][arm]
    state = (lambda x: relaxed_probs(comp, x)) if comp is not None else (lambda x: qaoa_probs(P["diag"], x, REPS))
    seed = zlib.crc32(f"newcode|{key[0]}|{key[1]}|{arm}|{opt_name}|{trial}".encode())
    algorithm_globals.random_seed = seed
    np.random.seed(seed % 2 ** 32)
    rng = np.random.default_rng(seed)
    npar = len(x0)
    st = dict(n=0, f=np.inf, x=np.array(x0, float), at=0, ck={})
    cks = set(range(25, budget + 1, 25)) | {1} if os.environ.get("TRACE") else set(CHECKPOINTS)

    def one(x):
        if st["n"] >= budget:
            return st["f"]
        st["n"] += 1
        f = sample_energy(state(x), P["E"], SHOTS, rng)
        if f < st["f"]:
            st.update(f=f, x=np.array(x, float), at=st["n"])
        if st["n"] in cks:
            st["ck"][st["n"]] = st["x"].copy()
        return f

    def fun(x):
        x = np.asarray(x, float).reshape(-1, npar)
        v = [one(r) for r in x]
        return v[0] if len(v) == 1 else np.array(v)

    opt = make_optimizer(opt_name, ITER, budget)
    bounds = [(v - np.pi, v + np.pi) for v in x0] if isinstance(opt, BOBYQA) else None
    t0 = time.time()
    opt.minimize(fun=fun, x0=np.array(x0, float), bounds=bounds)
    s, e = metrics(P, state(np.array(x0, float))), metrics(P, state(st["x"]))
    span = s["E"] - P["gs"]
    return dict(formulation=key[0], problem=key[1], qubits=P["n"], n_params=npar, start=arm, optimizer=opt_name, trial=trial, seed=seed,
                budget=budget, evals=st["n"], best_eval=st["at"], E_ground=P["gs"], opt_cost=P["opt"],
                E_start=s["E"], P_fb_start=s["P_fb"], P_raw_start=s["P_raw"],
                E_logged_best=st["f"], E_end=e["E"], P_fb_end=e["P_fb"], P_raw_end=e["P_raw"], P_ground_end=e["P_ground"],
                gap_closed=(s["E"] - e["E"]) / span if span > 1e-9 else None,
                anytime={k: float(state(v) @ P["E"]) for k, v in st["ck"].items()},
                anytime_P_fb=({k: metrics(P, state(v))["P_fb"] for k, v in st["ck"].items()} if os.environ.get("TRACE") else None),
                wall_s=time.time() - t0)


if __name__ == "__main__":
    keys = [tuple(k.split(":")) for k in sys.argv[1].split(",")]
    budget, trials, workers = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    if "--check" in sys.argv:
        for k in keys:
            t0 = time.time()
            P = problem(k)
            print(k, f"{P['n']} qubits, ground {P['gs']:.6f}, opt cost {P['opt']:g}, lift max |dp| {P['lift_dev']:.1e},",
                  ", ".join(f"{a}: {len(x)} params E {metrics(P, (relaxed_probs(c, np.array(x)) if c is not None else qaoa_probs(P['diag'], np.array(x), REPS)))['E']:.4f}"
                            for a, (c, x) in P["starts"].items()), f"({time.time() - t0:.0f}s)", flush=True)
        sys.exit()
    out = os.path.join(HERE, os.environ.get("OUT", f"runs_b{budget}.jsonl"))
    done = set()
    if os.path.exists(out):                      # resumable
        for l in open(out):
            r = json.loads(l)
            done.add((r["formulation"], r["problem"], r["start"], r["optimizer"], r["trial"]))
    tasks = [(k, a, o, t, budget) for k in keys for a in (["std4", "lift"] + list(SPIQ[k])) for o in OPTIMIZERS
             for t in range(1, trials + 1) if (k[0], k[1], a, o, t) not in done]
    only = [a for a in os.environ.get("STARTS", "").split(",") if a]     # e.g. STARTS=std4,spiq_single
    if only:
        tasks = [t for t in tasks if t[1] in only]
    if os.environ.get("TRIAL_MAJOR"):                                      # every cell gets trial 1 before any gets trial 2
        tasks.sort(key=lambda t: t[3])
    if os.environ.get("STD4_FIRST"):                                       # all std4 trials before remaining spiq trials
        tasks.sort(key=lambda t: (t[1] != "std4" and t[3] > 1, t[3]))      # (spiq trial 1 stays high priority)
    print(len(tasks), "runs to do", flush=True)
    with Pool(workers) as pool, open(out, "a") as f:
        for r in pool.imap_unordered(run, tasks):
            f.write(json.dumps(r, default=float) + "\n")
            f.flush()
            print(r["formulation"], r["problem"], r["start"], r["optimizer"], r["trial"],
                  f"E {r['E_start']:.3f}->{r['E_end']:.3f} (ground {r['E_ground']:.3f}) fb {r['P_fb_start']:.3f}->{r['P_fb_end']:.3f} evals {r['evals']} {r['wall_s']:.0f}s", flush=True)
    print("ALL DONE", flush=True)
