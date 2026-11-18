"""Kick 1.0 x 5 + NFT on P5 (24 qubits, 252 angles) from the team's single-start SPIQ point.

Same pipeline as round 5 ("fixed1.0" in restarts.py): measure the SPIQ start once (candidate 0), then 5 times
[add Gaussian noise sigma = 1.0 rad to every angle of the SPIQ point, run NFT for up to 3,000 evaluations], and
return the candidate with the lowest measured 10,240-shot energy. Selection never sees P(opt) or the ground state.

Extra logging for the graphs (does not affect the run):
  trace      every 100 evaluations: lowest measured energy so far (whole run)
  segments   after each restart: exact energy and exact P(opt) of that restart's best point, and of the
             best-so-far candidate

Seeds are new (label "p5kick|..."), so these runs share nothing with rounds 4-5.
Usage (QAOA env): python p5_kick.py <formulation worktree> <reps> <workers>
"""
import json, os, sys, time, zlib
import numpy as np
from multiprocessing import Pool
import pilot as PL
from optimizers import make_optimizer
from fastqaoa import sample_energy
from relaxed_sim import probs as relaxed_probs
from qiskit.utils import algorithm_globals

KICK, RESTARTS, FULL = 1.0, 5, 3000
OUT = os.path.join(PL.HERE, "runs_p5_kick.jsonl")


def run(t):
    repo, rep = t
    P, arms = PL.problem("P5", repo)
    comp, x0 = arms["spiq_single"]
    x0 = np.array(x0, float)
    seed = zlib.crc32(f"p5kick|P5|spiq_single|fixed1.0|{rep}".encode())
    rng = np.random.default_rng(seed)
    t0 = time.time()
    G = dict(n=0, best=np.inf, trace=[])

    def measured(x):
        f = sample_energy(relaxed_probs(comp, np.asarray(x, float)), P["E"], PL.SHOTS, rng)
        G["n"] += 1
        G["best"] = min(G["best"], f)
        if G["n"] % 100 == 0:
            G["trace"].append([G["n"], G["best"]])
        return f

    def segment(xs, budget, sub):
        algorithm_globals.random_seed = sub
        np.random.seed(sub % 2 ** 32)
        st = dict(n=0, f=np.inf, x=np.array(xs, float))

        def fun(x):
            if st["n"] >= budget:
                return st["f"]
            st["n"] += 1
            f = measured(x)
            if f < st["f"]:
                st.update(f=f, x=np.array(x, float))
            return f
        make_optimizer("NFT", PL.ITER, budget).minimize(fun=fun, x0=np.array(xs, float))
        return st["f"], st["x"]

    cands = [(measured(x0), x0)]
    s = PL.metrics(P, relaxed_probs(comp, x0))
    segs = [dict(restart=0, E_seg=s["E"], P_fb_seg=s["P_fb"], E_best=s["E"], P_fb_best=s["P_fb"], evals=G["n"], wall_s=time.time() - t0)]
    for k in range(1, RESTARTS + 1):
        f, x = segment(x0 + KICK * rng.standard_normal(len(x0)), FULL, int(rng.integers(2 ** 31)))
        cands.append((f, x))
        m = PL.metrics(P, relaxed_probs(comp, x))
        b = PL.metrics(P, relaxed_probs(comp, min(cands, key=lambda c: c[0])[1]))
        segs.append(dict(restart=k, E_logged=f, E_seg=m["E"], P_fb_seg=m["P_fb"], E_best=b["E"], P_fb_best=b["P_fb"], evals=G["n"], wall_s=time.time() - t0))
        print(f"rep {rep} restart {k}: E {m['E']:.3f} fb {m['P_fb']:.3f} | best E {b['E']:.3f} fb {b['P_fb']:.3f} | {G['n']} evals {time.time()-t0:.0f}s", flush=True)
    f, x = min(cands, key=lambda c: c[0])
    e = PL.metrics(P, relaxed_probs(comp, x))
    return dict(problem="P5", arm="spiq_single", pipeline="fixed1.0", rep=rep, seed=seed, evals=G["n"], E_logged=f, E_ground=P["gs"], opt_cost=P["opt"],
                E_start=s["E"], P_fb_start=s["P_fb"], P_raw_start=s["P_raw"], E_end=e["E"], P_fb_end=e["P_fb"], P_raw_end=e["P_raw"],
                P_ground_end=e["P_ground"], segments=segs, trace=G["trace"], wall_s=time.time() - t0)


if __name__ == "__main__":
    repo, reps, workers = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    done = set()
    if os.path.exists(OUT):
        done = {json.loads(l)["rep"] for l in open(OUT)}
    tasks = [(repo, r) for r in range(1, reps + 1) if r not in done]
    print(len(tasks), "runs to do", flush=True)
    with Pool(workers) as pool, open(OUT, "a") as fo:
        for r in pool.imap_unordered(run, tasks):
            fo.write(json.dumps(r, default=float) + "\n")
            fo.flush()
            print("DONE rep", r["rep"], f"E {r['E_start']:.3f}->{r['E_end']:.3f} (ground {r['E_ground']:.3f}) fb {r['P_fb_start']:.3f}->{r['P_fb_end']:.3f} {r['wall_s']/3600:.1f}h", flush=True)
    print("ALL DONE", flush=True)
