"""Why kick + NFT: the energy landscape around the team's SPIQ points (new formulation), measured exactly.

For every SPIQ point (P1-P4 single/multi, and the 12 new problems):
  1. Check that energy is an exact sinusoid in each angle (5-point fit on 4 random angles). This is what makes
     NFT's 1-D fits exact and the parameter-shift rule exact.
  2. Exact gradient (parameter shift) and single-angle curvature at the SPIQ point:
       g_i = [E(x+pi/2 e_i) - E(x-pi/2 e_i)]/2,  c_i = [E(x+pi/2 e_i) + E(x-pi/2 e_i) - 2E(x)]/2
     c_i > 0: uphill both ways, c_i ~ 0: flat, c_i < 0: downhill.
  3. Kick test: 200 random kicks of size sigma in {0.2, 0.5, 1.0, 2.0} (Gaussian on every angle), exact E and P(opt)
     of the kicked point BEFORE any optimizer runs.
  4. 1-D cuts: E and P(opt) along x0 + t*d, t in [0, 3], d ~ N(0, I) (so t is the kick size), 12 directions.
All energies are exact statevector values (no shots). Output: ../data/landscape.json
Usage (QAOA env): python landscape.py
"""
import json, os, sys, contextlib, io
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OS = os.path.join(HERE, "..", "..")
sys.path[:0] = [os.path.join(OS, "new_formulation_pilot"), os.path.join(OS, "new_problems"), os.path.join(OS, "new_code_study"), os.path.join(OS, "many_problems")]
import pilot as PL
import run_new as RN
from relaxed_sim import probs as relaxed_probs

REPO = os.path.expanduser("~/Downloads/umich-db-quantumdb-wt/formulation")
OUT = os.path.join(HERE, "..", "data", "landscape.json")
rng = np.random.default_rng(20261118)
SIGMAS = [0.2, 0.5, 1.0, 2.0]


def analyse(label, comp, x0, P, metrics):
    x0 = np.array(x0, float)
    n = len(x0)
    E = lambda x: float(relaxed_probs(comp, x) @ P["E"])
    M = lambda x: metrics(P, relaxed_probs(comp, x))
    f0 = E(x0)
    # 1. sinusoid check
    res = []
    for i in rng.choice(n, size=min(4, n), replace=False):
        ts = np.array([0, 0.7, 1.9, 3.1, 4.4])
        ys = np.array([E(x0 + t * np.eye(n)[i]) for t in ts])
        A = np.c_[np.ones(5), np.cos(ts), np.sin(ts)]
        coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
        res.append(float(np.abs(A @ coef - ys).max()))
    # 2. gradient and curvature
    up = np.array([E(x0 + np.pi / 2 * np.eye(n)[i]) for i in range(n)])
    dn = np.array([E(x0 - np.pi / 2 * np.eye(n)[i]) for i in range(n)])
    g, c = (up - dn) / 2, (up + dn - 2 * f0) / 2
    tol = 1e-9 * max(1, abs(f0))
    xr = x0 + rng.normal(0, 0.3, n)
    gr = np.array([(E(xr + np.pi / 2 * np.eye(n)[i]) - E(xr - np.pi / 2 * np.eye(n)[i])) / 2 for i in range(n)])
    # 3. kicks
    m0 = M(x0)
    kick = {}
    for s in SIGMAS:
        ms = [M(x0 + s * rng.standard_normal(n)) for _ in range(200)]
        kick[str(s)] = dict(E=[m["E"] for m in ms], P=[m["P_fb"] for m in ms])
    # 4. cuts
    ts = np.linspace(0, 3, 31)
    cuts = []
    for _ in range(12):
        d = rng.standard_normal(n)
        ms = [M(x0 + t * d) for t in ts]
        cuts.append(dict(E=[m["E"] for m in ms], P=[m["P_fb"] for m in ms]))
    r = dict(label=label, angles=n, qubits=int(np.log2(len(P["E"]))), E_spiq=f0, E_ground=float(P["gs"]), P_spiq=m0["P_fb"],
             sinusoid_max_residual=max(res), grad_norm=float(np.linalg.norm(g)), grad_norm_random=float(np.linalg.norm(gr)),
             uphill=int((c > tol).sum()), flat=int((np.abs(c) <= tol).sum()), downhill=int((c < -tol).sum()),
             kick=kick, cut_t=list(ts), cuts=cuts)
    print(f"{label:28s} {n:4d} angles  |grad| {r['grad_norm']:.1e} (random pt {r['grad_norm_random']:.2f})  up/flat/down {r['uphill']}/{r['flat']}/{r['downhill']}"
          f"  sinusoid resid {r['sinusoid_max_residual']:.1e}  E {f0:.3f} ground {P['gs']:.3f} P {m0['P_fb']:.3f}", flush=True)
    return r


out = {}
for name, arm in [("P1", "spiq_single"), ("P2", "spiq_single"), ("P2", "spiq_multi"), ("P3", "spiq_single"), ("P4", "spiq_single")]:
    P, arms = PL.problem(name, REPO)
    comp, x0 = arms[arm]
    out[f"{name}|{arm}"] = analyse(f"{name} {arm}", comp, x0, P, PL.metrics)
for pid in RN.PROBS:
    with contextlib.redirect_stdout(io.StringIO()):
        P = RN.problem("qaoa", pid)
    comp, x0 = P["starts"]["spiq"]
    out[f"new|{pid}"] = analyse(pid, comp, x0, P, RN.S.metrics)
json.dump(out, open(OUT, "w"))
print("wrote", OUT)
