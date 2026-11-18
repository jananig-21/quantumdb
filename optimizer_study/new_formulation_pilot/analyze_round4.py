"""Round 4 tables. Every number is over 20 independent pipeline runs."""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
pd.set_option("display.width", 220, "display.max_rows", 200, "display.max_columns", 30)
r = pd.read_json("runs_round4.jsonl", lines=True)
print("rows", len(r), "cells with n!=20:", (r.groupby(["problem", "arm", "pipeline"]).size() != 20).sum())
order = ["cobyla", "nft", "fixed0.5", "fixed1.0", "portfolio", "halving"]
r["pipeline"] = pd.Categorical(r.pipeline, order)
g = r.groupby(["problem", "arm", "pipeline"], observed=True).agg(evals=("evals", "median"), E_start=("E_start", "mean"), E_mean=("E_end", "mean"), E_sd=("E_end", "std"), E_min=("E_end", "min"), ground=("E_ground", "first"),
    fb_start=("P_fb_start", "mean"), fb_mean=("P_fb_end", "mean"), fb_sd=("P_fb_end", "std"), fb_med=("P_fb_end", "median"), fb_min=("P_fb_end", "min"), raw_mean=("P_raw_end", "mean"),
    over90=("P_fb_end", lambda s: int((s > .9).sum())), at_ground=("E_end", lambda s: 0), wall=("wall_s", "median")).round(3).reset_index().drop(columns="at_ground")
print(g.to_string(index=False))
print("\nhalving vs cobyla and nft (one-sided Mann-Whitney)")
for (p, a), x in r.groupby(["problem", "arm"]):
    h = x[x.pipeline == "halving"]
    for b in ("cobyla", "nft"):
        y = x[x.pipeline == b]
        same_p = np.allclose(h.P_fb_end.values, h.P_fb_end.values[0]) and np.allclose(y.P_fb_end.values, h.P_fb_end.values[0])
        pv_p = 1.0 if same_p else mannwhitneyu(h.P_fb_end, y.P_fb_end, alternative="greater").pvalue
        pv_e = mannwhitneyu(h.E_end, y.E_end, alternative="less").pvalue
        print(p, a, "vs", b, f"P(opt) {h.P_fb_end.mean():.3f} vs {y.P_fb_end.mean():.3f} p={pv_p:.2e} | energy {h.E_end.mean():.3f} vs {y.E_end.mean():.3f} p={pv_e:.2e}")
print("\nkick picked by halving")
print(r[r.pipeline == "halving"].groupby(["problem", "arm"]).picked_kick.value_counts().unstack(fill_value=0))
print("\nwinner: does selection pick the higher-P(opt) outcome? rank correlation energy vs P(opt), kicked pipelines")
for (p, a), x in r[r.pipeline.isin(["fixed0.5", "fixed1.0", "portfolio", "halving"])].groupby(["problem", "arm"]):
    if x.P_fb_end.std() > 1e-6 and x.E_end.std() > 1e-6:
        print(p, a, f"rho {spearmanr(x.E_end, x.P_fb_end)[0]:.3f}")
print("\nP4 halving, each run"); print(r[(r.problem == "P4") & (r.pipeline == "halving")].sort_values("E_end")[["rep", "E_end", "P_fb_end", "P_raw_end", "picked_kick"]].round(3).to_string(index=False))
