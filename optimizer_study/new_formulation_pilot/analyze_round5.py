"""Round 5: fixed1.0 vs COBYLA on seeds not used to choose it (reps 101-120), plus exact 95% intervals on the success counts."""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, beta
pd.set_option("display.width", 200, "display.max_columns", 30)
r = pd.read_json("runs_round5_confirm.jsonl", lines=True)
print("rows", len(r), "reps", r.rep.min(), r.rep.max(), "cells n!=20:", (r.groupby(["problem", "arm", "pipeline"]).size() != 20).sum())
for (p, a), x in r.groupby(["problem", "arm"]):
    f, c = x[x.pipeline == "fixed1.0"], x[x.pipeline == "cobyla"]
    k = int((f.P_fb_end > .9).sum()); n = len(f)
    lo = beta.ppf(.025, k, n - k + 1) if k else 0.0; hi = beta.ppf(.975, k + 1, n - k) if k < n else 1.0
    same = np.allclose(f.P_fb_end.mean(), c.P_fb_end.mean()) and f.P_fb_end.std() < 1e-6
    pv = 1.0 if same else mannwhitneyu(f.P_fb_end, c.P_fb_end, alternative="greater").pvalue
    pe = 1.0 if f.E_end.std() < 1e-6 else mannwhitneyu(f.E_end, c.E_end, alternative="less").pvalue
    g = f.E_ground.iloc[0]
    print(f"{p} {a}: start P {c.P_fb_start.mean():.3f} E {c.E_start.mean():.3f} ground {g:.3f} | COBYLA P {c.P_fb_end.mean():.3f} E {c.E_end.mean():.3f} evals {c.evals.median():.0f} | kick P mean {f.P_fb_end.mean():.3f} med {f.P_fb_end.median():.3f} sd {f.P_fb_end.std():.3f} min {f.P_fb_end.min():.3f} raw {f.P_raw_end.mean():.3f} E {f.E_end.mean():.3f} min {f.E_end.min():.3f} | over 0.9: {k}/{n} [{lo:.2f},{hi:.2f}] | p(P) {pv:.1e} p(E) {pe:.1e}")
x = r[(r.problem == "P4") & (r.pipeline == "fixed1.0")].sort_values("E_end")
print(x[["rep", "E_end", "P_fb_end", "P_raw_end", "picked_kick"]].round(3).to_string(index=False))
