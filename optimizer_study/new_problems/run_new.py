"""The 2x2 comparison (old/new formulation x standard/SPIQ start) on 12 NEW 4-relation problems, plus kick + NFT.

Problems: new_problems.jsonl (3 CHAIN, 3 STAR, 3 CYCLE from the generated set; 3 brand-new CLIQUEs).
Formulations: q2o (old, team's v2 QUBO; 14-20 qubits) and qaoa (Divya's new QUBO; 8 qubits).
Starts: std4 (4 shared angles, all 0.5) and spiq (spiq_<form>/<id>.json, team's SPIQ code, 400 generations).
Pipelines:
  COBYLA, NFT  1,000 evaluations, 3 trials each (same scoring as ../new_code_study/study.py, TRACE every 25 evals)
  kick         SPIQ start only: SPIQ point + 5 x [kick sigma = 1.0 rad, NFT 3,000 evaluations], 3 trials
Plan costs use Divya's fixed checker (every predicate applies) in both formulations, so P(optimal plan) is comparable.
Start-up checks: the SPIQ start's exact energy equals SPIQ's own (stim) energy, and the new-formulation energy
equals the team's QUBO evaluated by qiskit on random bitstrings.

Usage (QAOA env): python run_new.py <workers>   -> runs_new.jsonl (resumable)
"""
import contextlib, importlib.util, io, json, os, sys, time, zlib
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [os.path.join(HERE, "..", "new_code_study"), os.path.join(HERE, "..", "many_problems")]
import study as S
import run_many as RM
from fastqaoa import qaoa_probs, sample_energy
from relaxed_sim import compile_layers, probs as relaxed_probs
from optimizers import make_optimizer
from qiskit.utils import algorithm_globals

PROBS = {p["id"]: p for p in map(json.loads, open(os.path.join(HERE, "new_problems.jsonl")))}
OUT = os.path.join(HERE, "runs_new.jsonl")
BUDGET, TRIALS, KICK, RESTARTS, FULL = 1000, 3, 1.0, 5, 3000
_C = {}


def problem(form, pid):
    if (form, pid) in _C:
        return _C[form, pid]
    p = PROBS[pid]
    card, pred, sel = p["card"], [tuple(x) for x in p["pred"]], p["pred_sel"]
    sp = json.load(open(os.path.join(HERE, f"spiq_{form}", pid + ".json")))
    if form == "qaoa":
        E = S.qaoa_energies(card, pred, sel, sp["penalty_weight"])
    else:
        with contextlib.redirect_stdout(io.StringIO()):
            E = RM.build_problem(p)["E"]
    n = int(np.log2(len(E)))
    raw, fb, opt = S.decode_tables(card, pred, sel, n)
    assert np.isclose(opt, p["opt_cost"]), (pid, opt)
    comp = compile_layers(sp["gates"]["num_qubits"], sp["gates"]["ops"])
    assert sp["gates"]["num_qubits"] == n
    x0 = np.array(sp["relaxed_initial_point"], float)
    e0 = float(relaxed_probs(comp, x0) @ E)
    assert abs(e0 - sp["energy_best_qubo"]) < 1e-6 * max(1, abs(e0)), (form, pid, e0, sp["energy_best_qubo"])
    gs = float(E.min())
    P = dict(key=(form, pid), n=n, E=E, diag=E - E.mean(), gs=gs, raw=raw, fb=fb, opt=opt,
             ground_is_opt=bool(np.isclose(fb[np.isclose(E, gs, rtol=0, atol=1e-9)], opt).all()),
             starts={"std4": (None, [0.5] * 4), "spiq": (comp, list(x0))})
    _C[form, pid] = P
    return P


def run(t):
    form, pid, arm, pipe, trial = t
    P = problem(form, pid)
    comp, x0 = P["starts"][arm]
    x0 = np.array(x0, float)
    state = (lambda x: relaxed_probs(comp, x)) if comp is not None else (lambda x: qaoa_probs(P["diag"], x, 2))
    seed = zlib.crc32(f"newprob|{form}|{pid}|{arm}|{pipe}|{trial}".encode())
    rng = np.random.default_rng(seed)
    ck = set(range(25, BUDGET + 1, 25)) | {1}
    G = dict(n=0, f=np.inf, x=x0.copy(), trace={})

    def measured(x):
        G["n"] += 1
        f = sample_energy(state(np.asarray(x, float)), P["E"], S.SHOTS, rng)
        if f < G["f"]:
            G.update(f=f, x=np.array(x, float))
        if pipe != "kick" and G["n"] in ck:
            m = S.metrics(P, state(G["x"]))
            G["trace"][G["n"]] = [m["E"], m["P_fb"]]
        elif pipe == "kick" and G["n"] % 250 == 0:
            m = S.metrics(P, state(G["x"]))
            G["trace"][G["n"]] = [m["E"], m["P_fb"]]
        return f

    def segment(opt_name, xs, budget):
        sub = int(rng.integers(2 ** 31))
        algorithm_globals.random_seed = sub
        np.random.seed(sub % 2 ** 32)
        st = dict(n=0)

        def fun(x):
            x = np.asarray(x, float).reshape(-1, len(x0))
            v = []
            for r in x:
                if st["n"] >= budget:
                    v.append(G["f"]); continue
                st["n"] += 1
                v.append(measured(r))
            return v[0] if len(v) == 1 else np.array(v)
        make_optimizer(opt_name, S.ITER, budget).minimize(fun=fun, x0=np.array(xs, float))

    t0 = time.time()
    if pipe == "kick":
        measured(x0)                                         # the SPIQ start is candidate 0
        for _ in range(RESTARTS):
            segment("NFT", x0 + KICK * rng.standard_normal(len(x0)), FULL)
    else:
        segment(pipe, x0, BUDGET)
    s, e = S.metrics(P, state(x0)), S.metrics(P, state(G["x"]))
    span = s["E"] - P["gs"]
    return dict(formulation=form, problem=pid, graph=PROBS[pid]["graph"], qubits=P["n"], n_params=len(x0), start=arm, pipeline=pipe,
                trial=trial, seed=seed, evals=G["n"], E_ground=P["gs"], ground_is_opt=P["ground_is_opt"], opt_cost=P["opt"],
                E_start=s["E"], P_fb_start=s["P_fb"], E_end=e["E"], P_fb_end=e["P_fb"], P_raw_end=e["P_raw"], P_ground_end=e["P_ground"],
                gap_closed=(s["E"] - e["E"]) / span if span > 1e-9 else None, trace=G["trace"], wall_s=time.time() - t0)


if __name__ == "__main__":
    workers = int(sys.argv[1])
    forms = sys.argv[2].split(",") if len(sys.argv) > 2 else ["qaoa", "q2o"]
    done = set()
    if os.path.exists(OUT):
        done = {(r["formulation"], r["problem"], r["start"], r["pipeline"], r["trial"]) for r in map(json.loads, open(OUT))}
    tasks = []
    for form in forms:
        for pid in PROBS:
            if not os.path.exists(os.path.join(HERE, f"spiq_{form}", pid + ".json")):
                continue
            for arm in ("std4", "spiq"):
                for pipe in ("COBYLA", "NFT") + (("kick",) if arm == "spiq" else ()):
                    for tr in range(1, TRIALS + 1):
                        if (form, pid, arm, pipe, tr) not in done:
                            tasks.append((form, pid, arm, pipe, tr))
    tasks.sort(key=lambda t: (t[3] == "kick", t[0] == "q2o"))       # quick runs first, kicks last
    print(len(tasks), "runs to do", flush=True)
    with Pool(workers) as pool, open(OUT, "a") as f:
        for r in pool.imap_unordered(run, tasks):
            f.write(json.dumps(r, default=float) + "\n")
            f.flush()
            print(r["formulation"], r["problem"], r["start"], r["pipeline"], r["trial"],
                  f"E {r['E_start']:.3f}->{r['E_end']:.3f} (ground {r['E_ground']:.3f}) fb {r['P_fb_start']:.3f}->{r['P_fb_end']:.3f} evals {r['evals']} {r['wall_s']:.0f}s", flush=True)
    print("ALL DONE", flush=True)
