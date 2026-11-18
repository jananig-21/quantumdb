"""Round 4: real restarts and a kick-size rule that needs no tuning per problem.

Plan fixed before any result was seen. Each pipeline is run independently 20 times, with new seeds, and spends
the same total of 15,000 evaluations (plus 1 to measure the start, where it is a candidate). Each pipeline
returns the candidate with the lowest measured (10,240-shot) energy; it never sees P(opt) or the ground state.

  cobyla       COBYLA from the SPIQ point, no kick (the team's optimizer; it may stop early by itself)
  nft          NFT from the SPIQ point, no kick
  fixed0.5     5 restarts of 3,000 evaluations, kick 0.5
  fixed1.0     5 restarts of 3,000 evaluations, kick 1.0
  portfolio    5 restarts of 3,000 evaluations, one each at kick 0.2, 0.35, 0.5, 0.75, 1.0
  halving      15 short runs of 400 (3 at each of those 5 kicks), then the 3 lowest continue for 3,000 each

The primary pipeline is `halving`; it was chosen in advance because it needs no kick size per problem. Kicked
pipelines also keep the unkicked start as a candidate, so a kick cannot replace a start that is already better.

Usage (QAOA env): python restarts.py <formulation worktree> <workers>
"""
import json, os, sys, time, zlib
import numpy as np
from multiprocessing import Pool
import pilot as PL
from optimizers import make_optimizer
from fastqaoa import sample_energy
from relaxed_sim import probs as relaxed_probs
from qiskit.utils import algorithm_globals

KICKS = [0.2, 0.35, 0.5, 0.75, 1.0]
TOTAL, FULL, SHORT = 15000, 3000, 400
CASES = [("P1", "spiq_multi"), ("P3", "spiq_multi"), ("P2", "spiq_multi"), ("P2", "spiq_single"), ("P4", "spiq_single")]
PIPES = ["cobyla", "nft", "fixed0.5", "fixed1.0", "portfolio", "halving"]
REPS = 20
OUT = os.path.join(PL.HERE, "runs_round4.jsonl")
if os.environ.get("CONFIRM"):      # round 5: fixed1.0 was picked after seeing round 4, so it is re-tested on seeds it has not seen
    PIPES, OUT = ["cobyla", "fixed1.0"], os.path.join(PL.HERE, "runs_round5_confirm.jsonl")
FIRST_REP = 101 if os.environ.get("CONFIRM") else 1


def segment(P, comp, opt_name, x0, budget, rng, seed):
    """One optimizer run; returns (lowest measured energy, its angles, evaluations used)."""
    algorithm_globals.random_seed = seed
    np.random.seed(seed % 2 ** 32)
    st = dict(n=0, f=np.inf, x=np.array(x0, float))

    def fun(x):
        if st["n"] >= budget:
            return st["f"]
        st["n"] += 1
        f = sample_energy(relaxed_probs(comp, np.asarray(x, float)), P["E"], PL.SHOTS, rng)
        if f < st["f"]:
            st.update(f=f, x=np.array(x, float))
        return f

    make_optimizer(opt_name, PL.ITER, budget).minimize(fun=fun, x0=np.array(x0, float))
    return st["f"], st["x"], st["n"]


def run(t):
    repo, name, arm, pipe, rep = t
    P, arms = PL.problem(name, repo)
    comp, x0 = arms[arm]
    x0 = np.array(x0, float)
    label = f"round4|{name}|{arm}|{pipe}|{rep}"
    seed = zlib.crc32(label.encode())
    rng = np.random.default_rng(seed)
    t0 = time.time()
    used, cands = 0, []
    kicked = lambda k: x0 + k * rng.standard_normal(len(x0))
    sub = lambda: int(rng.integers(2 ** 31))
    if pipe in ("cobyla", "nft"):
        f, x, n = segment(P, comp, "COBYLA" if pipe == "cobyla" else "NFT", x0, TOTAL, rng, sub())
        cands, used, picked_kick = [(f, x, 0.0)], n, 0.0
    else:
        cands.append((sample_energy(relaxed_probs(comp, x0), P["E"], PL.SHOTS, rng), x0, 0.0))     # the unkicked start
        used = 1
        if pipe == "halving":
            short = []
            for k in KICKS:
                for _ in range(3):
                    f, x, n = segment(P, comp, "NFT", kicked(k), SHORT, rng, sub())
                    short.append((f, x, k))
                    used += n
            cands += short
            for f, x, k in sorted(short, key=lambda c: c[0])[:3]:
                f2, x2, n = segment(P, comp, "NFT", x, FULL, rng, sub())
                cands.append((f2, x2, k))
                used += n
        else:
            ks = KICKS if pipe == "portfolio" else [float(pipe[5:])] * 5
            for k in ks:
                f, x, n = segment(P, comp, "NFT", kicked(k), FULL, rng, sub())
                cands.append((f, x, k))
                used += n
    f, x, picked_kick = min(cands, key=lambda c: c[0])
    s, e = PL.metrics(P, relaxed_probs(comp, x0)), PL.metrics(P, relaxed_probs(comp, x))
    return dict(problem=name, arm=arm, pipeline=pipe, rep=rep, evals=used, picked_kick=picked_kick, E_logged=f, E_ground=P["gs"], opt_cost=P["opt"],
                E_start=s["E"], P_fb_start=s["P_fb"], P_raw_start=s["P_raw"], E_end=e["E"], P_fb_end=e["P_fb"], P_raw_end=e["P_raw"], wall_s=time.time() - t0)


if __name__ == "__main__":
    repo, workers = sys.argv[1], int(sys.argv[2])
    done = set()
    if os.path.exists(OUT):                      # resumable: a reboot loses only the runs in progress
        for l in open(OUT):
            r = json.loads(l)
            done.add((r["problem"], r["arm"], r["pipeline"], r["rep"]))
    tasks = [(repo, n, a, p, k) for n, a in CASES for p in PIPES for k in range(FIRST_REP, FIRST_REP + REPS) if (n, a, p, k) not in done]
    print(len(tasks), "runs to do", flush=True)
    with Pool(workers) as pool, open(OUT, "a") as f:
        for r in pool.imap_unordered(run, tasks):
            f.write(json.dumps(r, default=float) + "\n")
            f.flush()
            print(r["problem"], r["arm"], r["pipeline"], r["rep"], f"E {r['E_start']:.3f}->{r['E_end']:.3f} fb {r['P_fb_start']:.3f}->{r['P_fb_end']:.3f} raw {r['P_raw_end']:.3f} kick {r['picked_kick']} evals {r['evals']} {r['wall_s']:.0f}s", flush=True)
    print("ALL DONE", flush=True)
