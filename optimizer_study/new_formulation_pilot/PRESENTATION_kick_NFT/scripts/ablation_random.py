"""Control experiment: does the kick need to start from SPIQ, or would 5 fully random restarts + NFT do as well?

Same pipeline as kick 1.0 x5 + NFT (measure the SPIQ point once as candidate 0, 5 restarts of 3,000 NFT
evaluations, keep the lowest measured energy), except each restart starts from angles drawn uniformly in [-pi, pi)
instead of SPIQ + N(0, 1). New seeds (label "ablation|..."). Same circuits, shots, scoring.

Cases: P2 single-start and P2 multi-start (20 runs each, compare with round 5), and the 12 new problems in the new
formulation (3 runs each, compare with runs_new.jsonl kick). Output: ../data/runs_ablation_random.jsonl (resumable)
Usage (QAOA env): python ablation_random.py <workers>
"""
import json, os, sys, time, zlib, contextlib, io
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
OS = os.path.join(HERE, "..", "..")
sys.path[:0] = [OS, os.path.join(OS, "..", "base"), os.path.join(OS, "new_formulation_pilot"), os.path.join(OS, "new_problems"), os.path.join(OS, "new_code_study"), os.path.join(OS, "many_problems")]
REPO = os.path.expanduser("~/Downloads/umich-db-quantumdb-wt/formulation")
OUT = os.path.join(HERE, "..", "data", "runs_ablation_random.jsonl")
RESTARTS, FULL = 5, 3000


def load(case):
    import pilot as PL
    import run_new as RN
    kind, name, arm = case
    if kind == "pilot":
        P, arms = PL.problem(name, REPO)
        comp, x0 = arms[arm]
        return P, comp, np.array(x0, float), PL.metrics
    with contextlib.redirect_stdout(io.StringIO()):
        P = RN.problem("qaoa", name)
    comp, x0 = P["starts"]["spiq"]
    return P, comp, np.array(x0, float), RN.S.metrics


def run(t):
    from optimizers import make_optimizer
    from fastqaoa import sample_energy
    from relaxed_sim import probs as relaxed_probs
    from qiskit.utils import algorithm_globals
    case, rep = t
    P, comp, x0, metrics = load(case)
    seed = zlib.crc32(f"ablation|{case[1]}|{case[2]}|random_uniform|{rep}".encode())
    rng = np.random.default_rng(seed)
    t0 = time.time()
    G = dict(n=0)

    def measured(x):
        G["n"] += 1
        return sample_energy(relaxed_probs(comp, np.asarray(x, float)), P["E"], 10240, rng)

    def segment(xs, sub):
        algorithm_globals.random_seed = sub
        np.random.seed(sub % 2 ** 32)
        st = dict(n=0, f=np.inf, x=np.array(xs, float))

        def fun(x):
            if st["n"] >= FULL:
                return st["f"]
            st["n"] += 1
            f = measured(x)
            if f < st["f"]:
                st.update(f=f, x=np.array(x, float))
            return f
        make_optimizer("NFT", 10000, FULL).minimize(fun=fun, x0=np.array(xs, float))
        return st["f"], st["x"]

    cands = [(measured(x0), x0)]
    for _ in range(RESTARTS):
        cands.append(segment(rng.uniform(-np.pi, np.pi, len(x0)), int(rng.integers(2 ** 31))))
    f, x = min(cands, key=lambda c: c[0])
    s, e = metrics(P, relaxed_probs(comp, x0)), metrics(P, relaxed_probs(comp, x))
    picked_start = bool(np.array_equal(x, x0))
    return dict(kind=case[0], problem=case[1], arm=case[2], pipeline="random5", rep=rep, seed=seed, evals=G["n"], E_ground=float(P["gs"]),
                E_start=s["E"], P_fb_start=s["P_fb"], E_end=e["E"], P_fb_end=e["P_fb"], picked_spiq=picked_start, wall_s=time.time() - t0)


if __name__ == "__main__":
    workers = int(sys.argv[1])
    sys.path.insert(0, os.path.join(OS, "new_problems"))
    probs = [json.loads(l)["id"] for l in open(os.path.join(OS, "new_problems", "new_problems.jsonl"))]
    cases = [(("pilot", "P2", "spiq_single"), r) for r in range(1, 21)] + [(("pilot", "P2", "spiq_multi"), r) for r in range(1, 21)] + \
            [(("pilot", "P4", "spiq_single"), r) for r in range(1, 21)] + \
            [(("new", p, "spiq"), r) for p in probs for r in range(1, 4)]
    done = set()
    if os.path.exists(OUT):
        done = {(r["problem"], r["arm"], r["rep"]) for r in map(json.loads, open(OUT))}
    tasks = [(c, r) for c, r in cases if (c[1], c[2], r) not in done]
    print(len(tasks), "runs to do", flush=True)
    with Pool(workers) as pool, open(OUT, "a") as fo:
        for r in pool.imap_unordered(run, tasks):
            fo.write(json.dumps(r, default=float) + "\n")
            fo.flush()
            print(r["problem"], r["arm"], r["rep"], f"P {r['P_fb_start']:.3f}->{r['P_fb_end']:.3f} E {r['E_start']:.3f}->{r['E_end']:.3f} "
                  f"(ground {r['E_ground']:.3f}) picked_spiq {r['picked_spiq']} {r['wall_s']:.0f}s", flush=True)
    print("ALL DONE", flush=True)
