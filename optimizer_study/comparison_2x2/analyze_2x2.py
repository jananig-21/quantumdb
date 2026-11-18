"""Analysis and figures for the 2x2 comparison (old/new formulation x standard/SPIQ start), the 12 new problems,
and kick + NFT on P5. Reads only result files; safe to re-run while runs are still going.

Inputs
  ../new_code_study/runs_trace_b1000.jsonl   10 optimizers, P1-P3 (old q2o) and P1-P4 (new qaoa), 5 trials, trend every 25 evals
  ../new_code_study/runs_b1000.jsonl         Study 2 (new qaoa P5, 3 trials; also the originals the trace re-run must reproduce)
  ../new_problems/runs_new.jsonl             12 new problems x 2 formulations x 2 starts x {COBYLA, NFT, kick}, 3 trials
  ../new_formulation_pilot/runs_p5_kick.jsonl  kick 1.0 x 5 + NFT on P5
Outputs: figures/*.png, tables/*.csv, summary.json
"""
import json, os, sys
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
ST = os.path.join(HERE, "..")
FIG, TAB = os.path.join(HERE, "figures"), os.path.join(HERE, "tables")
os.makedirs(FIG, exist_ok=True); os.makedirs(TAB, exist_ok=True)

OPTS = ["COBYLA", "AQGD", "SPSA", "NELDER_MEAD", "POWELL", "BOBYQA", "L_BFGS_B", "SLSQP", "ADAM", "NFT"]
LABEL = {"NELDER_MEAD": "Nelder–Mead", "L_BFGS_B": "L-BFGS-B", "POWELL": "Powell"}
lab = lambda o: LABEL.get(o, o)
FORM = {"q2o": "Old formulation", "qaoa": "New formulation"}
STARTL = {"std4": "standard start", "spiq_single": "SPIQ start", "spiq": "SPIQ start", "spiq_multi": "multi-search SPIQ"}
CONDS = [("q2o", "std4"), ("q2o", "spiq_single"), ("qaoa", "std4"), ("qaoa", "spiq_single")]
cond_name = lambda f, s: f"{'Old' if f == 'q2o' else 'New'} · {'standard' if s == 'std4' else 'SPIQ'}"
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1", "#fcfcfb"
BLUE, ORANGE, AQUA, YELLOW, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#b9b8b3"
SEQ = LinearSegmentedColormap.from_list("seq", ["#f4f7fb", "#bcd3f0", "#5b9be6", "#2a78d6", "#0f4c99"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlecolor": INK, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False})


def load(p):
    p = os.path.join(ST, p)
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []


def save(fig, name):
    fig.savefig(os.path.join(FIG, name), dpi=170, bbox_inches="tight")
    plt.close(fig)


trace = pd.DataFrame(load("new_code_study/runs_trace_b1000.jsonl"))
study2 = pd.DataFrame(load("new_code_study/runs_b1000.jsonl"))
newp = pd.DataFrame(load("new_problems/runs_new.jsonl"))
p5k = load("new_formulation_pilot/runs_p5_kick.jsonl")
summary = {}

# ---------------------------------------------------------------------------------------------------------------
# 0. Reproducibility: the traced re-run uses the same seeds as Study 2, so its results must be identical.
if len(trace) and len(study2):
    k = ["formulation", "problem", "start", "optimizer", "trial"]
    m = trace.merge(study2[study2.formulation == "qaoa"], on=k, suffixes=("", "_s2"))
    if len(m):
        d = (m.E_end - m.E_end_s2).abs()
        summary["reproducibility"] = dict(matched_runs=int(len(m)), max_abs_diff_E=float(d.max()),
                                          identical_runs=int((d < 1e-9).sum()))

# ---------------------------------------------------------------------------------------------------------------
# 1. The 2x2 table: one row per (formulation, start, problem, optimizer), mean over trials
rows = []
if len(trace):
    allp = pd.concat([trace, study2[(study2.formulation == "qaoa") & (study2.problem == "P5") & study2.start.isin(["std4", "spiq_single"])]], ignore_index=True)
    for (f, s, p, o), g in allp.groupby(["formulation", "start", "problem", "optimizer"]):
        rows.append(dict(formulation=f, start=s, problem=p, optimizer=o, trials=len(g), qubits=int(g.qubits.iloc[0]),
                         n_params=int(g.n_params.iloc[0]), E_start=g.E_start.mean(), E_end=g.E_end.mean(), E_end_sd=g.E_end.std(),
                         E_ground=g.E_ground.iloc[0], gap_closed=g.gap_closed.mean() if g.gap_closed.notna().any() else np.nan,
                         P_fb_start=g.P_fb_start.mean(), P_fb_end=g.P_fb_end.mean(), P_fb_end_min=g.P_fb_end.min(),
                         P_fb_end_max=g.P_fb_end.max(), evals=g.evals.mean()))
T = pd.DataFrame(rows)
if len(T):
    T.to_csv(os.path.join(TAB, "table_2x2_all.csv"), index=False)


def heat(ax, M, rowlabels, collabels, title, vmin=0, vmax=1, fmt="{:.2f}", note=None):
    im = ax.imshow(np.ma.masked_invalid(M), cmap=SEQ, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(collabels))); ax.set_xticklabels(collabels, fontsize=8)
    ax.set_yticks(range(len(rowlabels))); ax.set_yticklabels(rowlabels, fontsize=8)
    ax.grid(False)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            txt = "—" if np.isnan(v) else fmt.format(v)
            if note is not None and isinstance(note[i][j], str):
                txt = note[i][j]
            dark = (not np.isnan(v)) and (v - vmin) / (vmax - vmin + 1e-12) > 0.55
            ax.text(j, i, txt, ha="center", va="center", fontsize=7.5, color="#ffffff" if dark else INK)
    ax.set_title(title, fontsize=10, loc="left", fontweight="bold")
    for sp in ax.spines.values():
        sp.set_visible(False)
    return im


def probs_for(f):
    return ["P1", "P2", "P3"] if f == "q2o" else ["P1", "P2", "P3", "P4", "P5"]


if len(T):
    # Figure 1 + 2: heatmaps, 10 optimizers x problems, for each of the four conditions
    for metric, fname, ttl in [("gap_closed", "fig1_heatmap_gap_closed.png", "Energy: share of the gap to the ground state closed (mean over trials; 1 = reached ground state)"),
                               ("P_fb_end", "fig2_heatmap_P_optimal_plan.png", "Optimal plan: P(optimal plan) at the end (mean over trials)")]:
        fig, axes = plt.subplots(2, 2, figsize=(11.5, 9), gridspec_kw=dict(width_ratios=[3, 5]))
        for ax, (f, s) in zip([axes[0, 0], axes[1, 0], axes[0, 1], axes[1, 1]], CONDS):
            ps = probs_for(f)
            M = np.full((len(OPTS), len(ps)), np.nan)
            note = [[None] * len(ps) for _ in OPTS]
            for j, p in enumerate(ps):
                for i, o in enumerate(OPTS):
                    r = T[(T.formulation == f) & (T.start == s) & (T.problem == p) & (T.optimizer == o)]
                    if len(r):
                        v = r[metric].iloc[0]
                        M[i, j] = v
                        if metric == "gap_closed" and np.isnan(v) and abs(r.E_start.iloc[0] - r.E_ground.iloc[0]) < 1e-6:
                            note[i][j] = "at ground"
            cl = []
            for p in ps:
                r = T[(T.formulation == f) & (T.start == s) & (T.problem == p)]
                if metric == "P_fb_end" and len(r):
                    cl.append(f"{p}\n(start {r.P_fb_start.iloc[0]:.2f})")
                elif len(r):
                    cl.append(f"{p}\n{int(r.qubits.iloc[0])} qubits")
                else:
                    cl.append(p)
            vmin = -0.5 if metric == "gap_closed" else 0
            heat(ax, np.clip(M, vmin, 1) if metric == "gap_closed" else M, [lab(o) for o in OPTS], cl, cond_name(f, s), vmin=vmin, vmax=1, note=note)
            if metric == "gap_closed":
                for i in range(M.shape[0]):
                    for j in range(M.shape[1]):
                        if not np.isnan(M[i, j]) and M[i, j] < vmin:
                            ax.texts[i * M.shape[1] + j].set_text(f"{M[i, j]:.1f}")
        fig.suptitle(ttl, x=0.01, ha="left", fontsize=11, fontweight="bold", color=INK)
        fig.text(0.01, 0.005, "Old formulation exists only for P1–P3 (P4 and P5 would need 30 and 52 qubits). SPIQ = the team's single-search SPIQ point. "
                 "P5: 3 trials, from Study 2." + (" Negative gap closed = ended worse than the start; 'at ground' = SPIQ already starts at the ground state." if metric == "gap_closed" else
                 " Start = P(optimal plan) of the starting point, before any optimizer runs."), fontsize=7.5, color=MUTED)
        fig.tight_layout(rect=(0, 0.02, 1, 0.97))
        save(fig, fname)

    # Figure 3: best P(opt) any optimizer reached, and COBYLA, per condition and problem (old vs new side by side)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, s in zip(axes, ["std4", "spiq_single"]):
        ps = ["P1", "P2", "P3", "P4", "P5"]
        x = np.arange(len(ps))
        for k, (f, col) in enumerate([("q2o", ORANGE), ("qaoa", BLUE)]):
            best, cob, st0 = [], [], []
            for p in ps:
                r = T[(T.formulation == f) & (T.start == s) & (T.problem == p)]
                best.append(r.P_fb_end.max() if len(r) else np.nan)
                cob.append(r[r.optimizer == "COBYLA"].P_fb_end.mean() if len(r) else np.nan)
                st0.append(r.P_fb_start.iloc[0] if len(r) else np.nan)
            off = -0.2 if k == 0 else 0.2
            ax.bar(x + off, best, width=0.36, color=col, label=f"{FORM[f]}: best optimizer", zorder=2)
            ax.scatter(x + off, cob, marker="D", s=26, color=INK, zorder=3, label="COBYLA (team)" if k == 1 else None)
            ax.scatter(x + off, st0, marker="_", s=260, color=MUTED, linewidths=2, zorder=3, label="start" if k == 1 else None)
            for xi, v in zip(x + off, best):
                if np.isnan(v):
                    ax.text(xi, 0.02, "n/a", ha="center", fontsize=7, color=MUTED)
        ax.set_xticks(x); ax.set_xticklabels(ps)
        ax.set_ylim(0, 1.05); ax.set_title(f"From the {STARTL[s]}", loc="left", fontsize=10, fontweight="bold")
    axes[0].set_ylabel("P(optimal plan) at the end")
    axes[1].legend(frameon=False, fontsize=8, loc="upper right")
    fig.suptitle("Old vs new formulation: how often the circuit returns the optimal plan", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig3_old_vs_new_P_optimal.png")

    # Figure 4 + 5: trends (mean over trials) of energy gap closed and P(opt) of the best-so-far point
    def curves(df, field):
        out = {}
        for o, g in df.groupby("optimizer"):
            xs = list(range(25, 1001, 25))
            Y = []
            for _, r in g.iterrows():
                d = {int(k): v for k, v in (r[field] or {}).items()}
                if not d:
                    continue
                last, y = None, []
                for x in xs:
                    if x in d:
                        last = d[x]
                    elif last is None:
                        last = d.get(1, np.nan)
                    y.append(last)
                if field == "anytime":
                    span = r.E_start - r.E_ground
                    y = [(r.E_start - v) / span if span > 1e-9 else np.nan for v in y]
                Y.append(y)
            if Y:
                out[o] = (xs, np.nanmean(np.array(Y, float), 0))
        return out

    for field, fname, ylab, ttl in [("anytime", "fig4_trend_energy.png", "gap closed (best-so-far point)", "Energy trend over the 1,000 evaluations"),
                                    ("anytime_P_fb", "fig5_trend_P_optimal.png", "P(optimal plan) of best-so-far point", "Optimal-plan trend over the 1,000 evaluations")]:
        probs = ["P1", "P2", "P3", "P4"]
        fig, axes = plt.subplots(len(probs), 4, figsize=(13, 2.5 * len(probs)), sharex=True)
        for j, (f, s) in enumerate(CONDS):
            for i, p in enumerate(probs):
                ax = axes[i, j]
                g = trace[(trace.formulation == f) & (trace.start == s) & (trace.problem == p)]
                if not len(g):
                    ax.text(0.5, 0.5, "n/a (old formulation:\nP4 needs 30 qubits)" if f == "q2o" else "no data", ha="center", va="center",
                            transform=ax.transAxes, color=MUTED, fontsize=8)
                    ax.set_yticks([]); ax.grid(False)
                elif field == "anytime" and abs(g.E_start.iloc[0] - g.E_ground.iloc[0]) < 1e-6:
                    ax.text(0.5, 0.5, "SPIQ starts at the ground state:\nno energy gap to close\n(every optimizer stays there)", ha="center", va="center",
                            transform=ax.transAxes, color=MUTED, fontsize=8)
                    ax.set_yticks([]); ax.grid(False)
                else:
                    cv = curves(g, field)
                    for o, (xs, y) in cv.items():
                        if o in ("COBYLA", "NFT"):
                            continue
                        ax.plot(xs, y, color=GRAY, lw=1, zorder=1)
                    for o, col in [("COBYLA", ORANGE), ("NFT", BLUE)]:
                        if o in cv:
                            ax.plot(*cv[o], color=col, lw=2, zorder=3)
                    best = max(cv, key=lambda o: cv[o][1][-1]) if cv else None
                    if best and best not in ("COBYLA", "NFT"):
                        ax.plot(*cv[best], color=AQUA, lw=2, zorder=2)
                        ax.text(1000, cv[best][1][-1], f" {lab(best)}", color=INK2, fontsize=7, va="center")
                    if field == "anytime_P_fb":
                        ax.set_ylim(-0.02, 1.02)
                if i == 0:
                    ax.set_title(cond_name(f, s), fontsize=10, fontweight="bold", loc="left")
                if j == 0:
                    ax.set_ylabel(f"{p}\n{ylab}", fontsize=8)
                if i == len(probs) - 1:
                    ax.set_xlabel("evaluations")
        handles = [plt.Line2D([], [], color=ORANGE, lw=2), plt.Line2D([], [], color=BLUE, lw=2), plt.Line2D([], [], color=AQUA, lw=2), plt.Line2D([], [], color=GRAY, lw=1)]
        fig.legend(handles, ["COBYLA (team)", "NFT", "best other optimizer (labelled)", "other optimizers"], loc="upper right", ncol=4, frameon=False, fontsize=8)
        fig.suptitle(ttl + " (mean of 5 trials; a flat line after an early stop means the optimizer stopped)", x=0.01, ha="left", fontsize=11, fontweight="bold")
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        save(fig, fname)

    # Summary numbers for the report
    for (f, s), g in T.groupby(["formulation", "start"]):
        sub = g[g.problem.isin(["P1", "P2", "P3"])]
        rk = sub.groupby("problem").apply(lambda d: d.set_index("optimizer").E_end.rank()).reset_index()
        mean_rank = rk.groupby("optimizer").E_end.mean().sort_values() if "E_end" in rk else None
        summary[f"{f}|{s}"] = dict(
            mean_P_fb_end_P1_P3=sub.groupby("optimizer").P_fb_end.mean().round(3).to_dict(),
            mean_gap_closed_P1_P3=sub.groupby("optimizer").gap_closed.mean().round(3).to_dict(),
            mean_rank_E_P1_P3=None if mean_rank is None else mean_rank.round(2).to_dict())

# ---------------------------------------------------------------------------------------------------------------
# 6. New problems
if len(newp):
    newp["cond"] = [cond_name(f, "std4" if s == "std4" else "spiq") for f, s in zip(newp.formulation, newp.start)]
    G = newp.groupby(["problem", "graph", "formulation", "start", "pipeline"]).agg(
        trials=("trial", "count"), P_fb_start=("P_fb_start", "mean"), P_fb_end=("P_fb_end", "mean"), gap=("gap_closed", "mean"),
        E_start=("E_start", "mean"), E_end=("E_end", "mean"), E_ground=("E_ground", "first"), ground_is_opt=("ground_is_opt", "first"),
        qubits=("qubits", "first")).reset_index()
    G.to_csv(os.path.join(TAB, "table_new_problems.csv"), index=False)
    pids = sorted(G.problem.unique(), key=lambda p: (["CHAIN", "STAR", "CYCLE", "CLIQUE"].index(p.split("_")[1]), p))
    cols = [(f, s, pp) for f in ("q2o", "qaoa") for s, pps in (("std4", ("COBYLA", "NFT")), ("spiq", ("COBYLA", "NFT", "kick"))) for pp in pps]
    M = np.full((len(pids), len(cols)), np.nan)
    for i, p in enumerate(pids):
        for j, (f, s, pp) in enumerate(cols):
            r = G[(G.problem == p) & (G.formulation == f) & (G.start == s) & (G.pipeline == pp)]
            if len(r):
                M[i, j] = r.P_fb_end.iloc[0]
    fig, ax = plt.subplots(figsize=(12, 0.42 * len(pids) + 2.2))
    cl = [f"{'Old' if f == 'q2o' else 'New'}\n{'std' if s == 'std4' else 'SPIQ'}\n{'kick+NFT' if pp == 'kick' else pp}" for f, s, pp in cols]
    rl = []
    for p in pids:
        q = G[(G.problem == p) & (G.formulation == "q2o")].qubits
        rl.append(f"{p.replace('Gen_', '').replace('_4r_s', ' #')} ({int(q.iloc[0]) if len(q) else '?'}q old / 8q new)")
    heat(ax, M, rl, cl, "12 new problems: P(optimal plan) at the end (mean of 3 trials)")
    ax.axvline(4.5, color=INK, lw=1.5)
    fig.tight_layout()
    save(fig, "fig6_new_problems_heatmap.png")

    # Figure 7: averaged over the 12 problems, start vs end, per condition and pipeline
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, metric, yl in [(axes[0], "P_fb_end", "P(optimal plan)"), (axes[1], "gap", "gap closed")]:
        x, labels = 0, []
        for f, col in (("q2o", ORANGE), ("qaoa", BLUE)):
            for s in ("std4", "spiq"):
                for pp in ("COBYLA", "NFT", "kick"):
                    r = G[(G.formulation == f) & (G.start == s) & (G.pipeline == pp)]
                    if not len(r):
                        continue
                    v = r[metric].mean()
                    ax.bar(x, v, color=col, alpha=1.0 if pp != "COBYLA" else 0.55, width=0.8, zorder=2)
                    if metric == "P_fb_end":
                        ax.scatter(x, r.P_fb_start.mean(), marker="_", s=200, color=INK, linewidths=2, zorder=3)
                    ax.text(x, (v if not np.isnan(v) else 0) + 0.02, f"{v:.2f}", ha="center", fontsize=7, color=INK2)
                    labels.append((x, f"{'Old' if f == 'q2o' else 'New'}\n{'std' if s == 'std4' else 'SPIQ'}\n{'kick' if pp == 'kick' else pp}"))
                    x += 1
            x += 0.6
        ax.set_xticks([a for a, _ in labels]); ax.set_xticklabels([b for _, b in labels], fontsize=7)
        ax.set_title(yl + (" (bar = end, — = start)" if metric == "P_fb_end" else ""), loc="left", fontsize=10, fontweight="bold")
    axes[0].set_ylim(0, 1.08)
    fig.suptitle("12 new problems, averaged: orange = old formulation, blue = new formulation", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig7_new_problems_summary.png")

    # Figure 8: does the lowest-energy state decode to the optimal plan? (new problems + P1-P3)
    gi = G.groupby(["problem", "formulation"]).ground_is_opt.first().unstack()
    summary["new_problems_ground_is_opt"] = {f: int(gi[f].sum()) for f in gi.columns} | {"n": int(len(gi))}
    summary["new_problems_mean"] = G.groupby(["formulation", "start", "pipeline"]).agg(P_start=("P_fb_start", "mean"), P_end=("P_fb_end", "mean"),
                                                                                          gap=("gap", "mean")).round(3).reset_index().to_dict("records")

    # Figure 9: trends on new problems (COBYLA vs NFT), averaged over problems and trials
    fig, axes = plt.subplots(1, 5, figsize=(15, 3.2), sharey=True, gridspec_kw=dict(width_ratios=[1, 1, 1, 1, 1.4]))
    for ax, (f, s) in zip(axes[:4], [("q2o", "std4"), ("q2o", "spiq"), ("qaoa", "std4"), ("qaoa", "spiq")]):
        for pp, col in (("COBYLA", ORANGE), ("NFT", BLUE)):
            g = newp[(newp.formulation == f) & (newp.start == s) & (newp.pipeline == pp)]
            Y = []
            for _, r in g.iterrows():
                d = {int(k): v[1] for k, v in r.trace.items()}
                last, y = d.get(1, r.P_fb_start), []
                for x in range(25, 1001, 25):
                    last = d.get(x, last); y.append(last)
                Y.append(y)
            if Y:
                ax.plot(range(25, 1001, 25), np.mean(Y, 0), color=col, lw=2, label=pp)
        ax.set_title(cond_name(f, "std4" if s == "std4" else "spiq"), loc="left", fontsize=10, fontweight="bold")
        ax.set_xlabel("evaluations"); ax.set_ylim(-0.02, 1.02)
    ax = axes[4]
    for f, col in (("q2o", ORANGE), ("qaoa", BLUE)):
        g = newp[(newp.formulation == f) & (newp.pipeline == "kick")]
        xs = list(range(0, 15001, 250))
        Y = []
        for _, r in g.iterrows():
            d = {int(k): v[1] for k, v in r.trace.items()}
            last, y = r.P_fb_start, []
            for x in xs:
                last = d.get(x, last); y.append(last)
            Y.append(y)
        if Y:
            ax.plot(xs, np.mean(Y, 0), color=col, lw=2, label=f"{FORM[f]} ({g.problem.nunique()} problems)")
    for k in range(1, 5):
        ax.axvline(1 + 3000 * k, color=GRAY, lw=0.8, ls=":")
    ax.set_title("Kick + NFT from SPIQ (15,000 evals)", loc="left", fontsize=10, fontweight="bold")
    ax.set_xlabel("evaluations (dotted = next kick)"); ax.legend(frameon=False, fontsize=8, loc="lower right")
    axes[0].set_ylabel("P(optimal plan)\nbest-so-far point"); axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("12 new problems: optimal-plan trend, COBYLA vs NFT (mean over problems and trials)", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig8_new_problems_trend.png")

    # Figure 10: per kick run, P(ground state) vs P(optimal plan): is the lowest-energy state the optimal plan?
    k = newp[newp.pipeline == "kick"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, f in zip(axes, ("q2o", "qaoa")):
        g = k[k.formulation == f]
        ok, bad = g[g.ground_is_opt], g[~g.ground_is_opt]
        ax.plot([0, 1], [0, 1], color=GRAY, lw=1, ls="--", zorder=1)
        ax.scatter(ok.P_ground_end, ok.P_fb_end, s=46, color=BLUE, edgecolor="white", lw=1.2, zorder=3, label="ground state = optimal plan")
        ax.scatter(bad.P_ground_end, bad.P_fb_end, s=46, color=ORANGE, marker="X", edgecolor="white", lw=0.8, zorder=4, label="ground state ≠ optimal plan")
        for _, r in bad.sort_values("P_ground_end").groupby("problem").tail(1).iterrows():
            if r.P_ground_end > 0.5:
                ax.annotate(r.problem.replace("Gen_", "").replace("_4r_s", " #"), (r.P_ground_end, r.P_fb_end), xytext=(-8, 10),
                            textcoords="offset points", fontsize=7.5, color=INK2, ha="right")
        nb = int((~G[G.formulation == f].groupby("problem").ground_is_opt.first()).sum())
        ax.set_title(f"{FORM[f]}\n{nb} of 12 problems: lowest-energy state ≠ optimal plan", loc="left", fontsize=9.5, fontweight="bold")
        ax.set_xlabel("P(lowest-energy state) at the end of kick + NFT"); ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03)
        ax.legend(frameon=False, fontsize=8, loc="upper left")
    axes[0].set_ylabel("P(optimal plan) at the end")
    fig.suptitle("Kick + NFT runs on the 12 new problems: reaching the ground state only helps if the ground state IS the optimal plan",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig10_ground_vs_optimal.png")

# ---------------------------------------------------------------------------------------------------------------
# 10. P5 kick + NFT
if p5k:
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    for r in p5k:
        ks = [s["restart"] for s in r["segments"]]
        axes[0].plot(ks, [s["E_best"] for s in r["segments"]], marker="o", ms=5, lw=2, label=f"run {r['rep']}")
        axes[1].plot(ks, [s["P_fb_best"] for s in r["segments"]], marker="o", ms=5, lw=2, label=f"run {r['rep']}")
    E0, Eg = p5k[0]["E_start"], p5k[0]["E_ground"]
    axes[0].axhline(Eg, color=INK, ls="--", lw=1); axes[0].text(0, Eg, " ground state 8.34", va="bottom", fontsize=8)
    c = study2[(study2.problem == "P5") & (study2.start == "spiq_single") & (study2.optimizer == "COBYLA")]
    if len(c):
        axes[0].axhline(c.E_end.mean(), color=ORANGE, lw=1.5); axes[0].text(5, c.E_end.mean(), " COBYLA from SPIQ", color=ORANGE, va="bottom", ha="right", fontsize=8)
    axes[0].set_xlabel("restart (0 = SPIQ start)"); axes[0].set_ylabel("exact energy of best candidate")
    axes[1].set_xlabel("restart (0 = SPIQ start)"); axes[1].set_ylabel("P(optimal plan) of best candidate"); axes[1].set_ylim(-0.02, 1.02)
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("P5 (24 qubits, 252 angles): kick 1.0 × 5 + NFT, progress per restart", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig9_P5_kick.png")
    summary["P5_kick"] = [dict(rep=r["rep"], E_start=r["E_start"], E_end=r["E_end"], P_fb_start=r["P_fb_start"], P_fb_end=r["P_fb_end"],
                               P_raw_end=r["P_raw_end"], hours=r["wall_s"] / 3600) for r in p5k]

json.dump(summary, open(os.path.join(HERE, "summary.json"), "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in summary.items() if k in ("reproducibility", "new_problems_ground_is_opt", "P5_kick")}, indent=1, default=float))
print("figures:", sorted(os.listdir(FIG)))
