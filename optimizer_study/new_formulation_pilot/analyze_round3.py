"""Held-out test (trials 101-120, equal 3,000-evaluation budget). Restart pipeline picks by logged (noisy) energy only."""
import glob, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
pd.set_option("display.width", 220, "display.max_rows", 200, "display.max_columns", 30)
r = pd.concat([pd.read_json(f, lines=True) for f in glob.glob("pilot_round3/runs_*_t100.jsonl")])
g = r.groupby(["problem", "arm", "optimizer"]).agg(n=("E_end", "size"), evals=("evals", "median"), E_start=("E_start", "mean"), E_mean=("E_end", "mean"), E_min=("E_end", "min"),
    ground=("E_ground", "first"), fb_start=("P_fb_start", "mean"), fb_mean=("P_fb_end", "mean"), fb_med=("P_fb_end", "median"), fb_min=("P_fb_end", "min"), fb_max=("P_fb_end", "max"),
    raw_mean=("P_raw_end", "mean"), over90=("P_fb_end", lambda s: int((s > .9).sum())), wall=("wall_s", "median")).round(3).reset_index()
print(g.to_string(index=False))
print("\nkick NFT vs same-start no-kick NFT and COBYLA (fallback P(opt), one-sided Mann-Whitney)")
for P, arm in (("P2", "spiq_multi"), ("P2", "spiq_single"), ("P4", "spiq_single")):
    for k in ("0.5", "1.0"):
        a = r[(r.problem == P) & (r.arm == f"{arm}+k{k}") & (r.optimizer == "NFT")].P_fb_end
        for bname, b in (("NFT", r[(r.problem == P) & (r.arm == arm) & (r.optimizer == "NFT")].P_fb_end), ("COBYLA", r[(r.problem == P) & (r.arm == arm) & (r.optimizer == "COBYLA")].P_fb_end)):
            print(P, arm, "k" + k, "vs", bname, f"mean {a.mean():.3f} vs {b.mean():.3f}, p={mannwhitneyu(a, b, alternative='greater').pvalue:.2e}")
print("\nbest of N restarts, picked by logged energy; 2000 bootstrap draws from the 20 runs (budget = N x 3000)")
rng = np.random.default_rng(0)
for P, arm in (("P2", "spiq_multi+k0.5"), ("P2", "spiq_multi+k1.0"), ("P2", "spiq_single+k1.0"), ("P4", "spiq_single+k0.5"), ("P4", "spiq_single+k1.0")):
    x = r[(r.problem == P) & (r.arm == arm) & (r.optimizer == "NFT")]
    for N in (1, 3, 5):
        idx = rng.integers(0, len(x), (2000, N))
        pick = idx[np.arange(2000), x.E_logged_best.values[idx].argmin(1)]
        v = x.P_fb_end.values[pick]
        print(P, arm, "N", N, f"mean P(opt) {v.mean():.3f}, share over 0.9: {(v > .9).mean():.3f}")
print("\nenergy vs fallback P(opt), rank correlation over all runs at this budget")
for P in ("P2", "P4"):
    for name, x in (("all SPIQ arms", r[(r.problem == P) & r.arm.str.startswith("spiq")]), ("kicked NFT only", r[(r.problem == P) & r.arm.str.contains(r"\+k") & (r.optimizer == "NFT")])):
        print(P, name, len(x), f"rho {spearmanr(x.E_end, x.P_fb_end)[0]:.3f}")
