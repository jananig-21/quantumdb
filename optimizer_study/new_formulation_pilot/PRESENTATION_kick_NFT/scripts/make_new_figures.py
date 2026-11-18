"""New figures for the kick + NFT presentation (the rest are copied from MEETING_2026-11-18, see copy_figures.sh).

  01_all10_optimizers_stuck_at_SPIQ.png   10 optimizers from the team's SPIQ point, new formulation, P1-P5 (Study 2)
  06_round1_which_optimizer_to_kick.png   round 1: 5 optimizers x 3 kick sizes on P2 (why NFT)
  07_round2_kick_size_sweep.png           round 2: kick size vs P(optimal plan), NFT vs COBYLA, P2 and P4
  12_new_problems_kick_heatmap.png        12 new problems x pipelines, incl. random-restart control and old formulation

Inputs (read only):
  ../../../new_code_study/runs_b1000.jsonl                       Study 2 (10 optimizers, P1-P5, new + old formulation)
  ../../pilot_round1/runs_P1_P2_P3_b1000_kicks.jsonl             round 1 kicks
  ../../pilot_round2/*.jsonl                                     round 2
  ../../../new_problems/runs_new.jsonl                           12 new problems
  ../../../MEETING_2026-11-18/data/runs_ablation_random.jsonl    random-restart control
Output: ../figures/, ../data/new_figures_numbers.json
Usage (QAOA env): python make_new_figures.py
"""
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PILOT = os.path.join(HERE, "..", "..")
OS = os.path.join(PILOT, "..")
FIG = os.path.join(HERE, "..", "figures")
DATA = os.path.join(HERE, "..", "data")
os.makedirs(FIG, exist_ok=True); os.makedirs(DATA, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#c9ced6", "axes.labelcolor": "#3b4350",
                     "xtick.color": "#3b4350", "ytick.color": "#3b4350", "axes.grid": False,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#fcfcfa", "axes.facecolor": "#fcfcfa", "savefig.dpi": 170})
RJ = lambda p: pd.read_json(p, lines=True)
NUM = {}
OPTS = ["COBYLA", "NFT", "BOBYQA", "POWELL", "NELDER_MEAD", "SLSQP", "L_BFGS_B", "ADAM", "AQGD", "SPSA"]


def heat(ax, M, rows, cols, txt, cmap, vmin, vmax, center=None):
    im = ax.imshow(M, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, fontsize=8.5)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=8.5)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if np.isnan(v):
                ax.text(j, i, "—", ha="center", va="center", fontsize=8, color="#8a919c"); continue
            norm = (v - vmin) / (vmax - vmin)
            dark = norm > 0.62 if center is None else abs(v) > 0.6 * max(abs(vmin), abs(vmax))
            ax.text(j, i, txt[i][j], ha="center", va="center", fontsize=7.6, color="white" if dark else "#14181f")
    for s in ax.spines.values():
        s.set_visible(False)
    return im


# ------------------------------------------------------------------ 01: all 10 optimizers from SPIQ, new formulation
S2 = RJ(os.path.join(OS, "new_code_study", "runs_b1000.jsonl"))
S2 = S2[(S2.formulation == "qaoa") & (S2.start == "spiq_single")]
probs = ["P1", "P2", "P3", "P4", "P5"]
P_end = np.full((10, 5), np.nan); dE = np.full((10, 5), np.nan); txtP = [[""] * 5 for _ in range(10)]; txtE = [[""] * 5 for _ in range(10)]
start = {}
for j, p in enumerate(probs):
    d = S2[S2.problem == p]
    start[p] = dict(E=float(d.E_start.iloc[0]), P=float(d.P_fb_start.iloc[0]), gs=float(d.E_ground.iloc[0]), n=int(d.groupby("optimizer").size().iloc[0]),
                    angles=int(d.n_params.iloc[0]), qubits=int(d.qubits.iloc[0]))
    for i, o in enumerate(OPTS):
        e = d[d.optimizer == o]
        if len(e):
            P_end[i, j] = e.P_fb_end.mean(); dE[i, j] = (e.E_end - e.E_start).mean()
            txtP[i][j] = f"{e.P_fb_end.mean():.2f}"; txtE[i][j] = f"{(e.E_end - e.E_start).mean():+.2f}"
            NUM[f"01|{p}|{o}"] = dict(P_end=float(e.P_fb_end.mean()), dE=float((e.E_end - e.E_start).mean()), n=len(e))
NUM["01|start"] = start
cols = [f"{p}\n{start[p]['qubits']} qubits\n{start[p]['angles']} angles\nstart {start[p]['P']:.2f}" for p in probs]
fig, ax = plt.subplots(1, 2, figsize=(13.5, 6.2))
im = heat(ax[0], P_end, OPTS, cols, txtP, "Blues", 0, 1)
ax[0].set_title("P(optimal plan) after 1,000 evaluations\n(column label: start = value at the SPIQ point)", loc="left", fontsize=10, fontweight="bold")
fig.colorbar(im, ax=ax[0], fraction=0.035)
lim = 1.2
im = heat(ax[1], np.clip(dE, -lim, lim), OPTS, cols, txtE, "RdBu_r", -lim, lim, center=0)
ax[1].set_title("Energy change from the SPIQ start (end − start)\nblue = lower (better), red = higher (worse); 0.00 = did not move", loc="left", fontsize=10, fontweight="bold")
fig.colorbar(im, ax=ax[1], fraction=0.035)
fig.suptitle("Why kick + NFT was needed: from the team's SPIQ point, all 10 optimizers stay where they started (new formulation)",
             x=0.01, ha="left", fontsize=11.5, fontweight="bold")
fig.text(0.01, 0.005, "Study 2 (new_code_study/runs_b1000.jsonl): single-start SPIQ point, 1,000 evaluations, 10,240 shots, depth 2; mean of 5 trials "
         "(P5: 3). P1 and P3 start at the optimum already. SPSA is the only one that moves much, and it moves uphill.", fontsize=8, color="#5b6370")
fig.tight_layout(rect=(0, 0.03, 1, 0.94))
fig.savefig(os.path.join(FIG, "01_all10_optimizers_stuck_at_SPIQ.png"), bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ 06: round 1, which optimizer to kick
K = RJ(os.path.join(PILOT, "pilot_round1", "runs_P1_P2_P3_b1000_kicks.jsonl"))
N0 = RJ(os.path.join(PILOT, "pilot_round1", "runs_P1_P2_P3_b1000.jsonl"))
k5 = ["NFT", "COBYLA", "POWELL", "ADAM", "BOBYQA"]
fig, axes = plt.subplots(2, 2, figsize=(13, 7.4))
for c, arm in enumerate(["spiq_multi", "spiq_single"]):
    e0 = float(N0[(N0.problem == "P2") & (N0.arm == arm)].E_start.iloc[0]); gs = float(N0[N0.problem == "P2"].E_ground.iloc[0])
    p0 = float(N0[(N0.problem == "P2") & (N0.arm == arm)].P_fb_start.iloc[0])
    colsk = ["no kick", "kick 0.05", "kick 0.2", "kick 0.5"]
    ME = np.full((5, 4), np.nan); MP = np.full((5, 4), np.nan); tE = [[""] * 4 for _ in range(5)]; tP = [[""] * 4 for _ in range(5)]
    for i, o in enumerate(k5):
        for j, kk in enumerate([None, "0.05", "0.2", "0.5"]):
            d = N0[(N0.problem == "P2") & (N0.arm == arm) & (N0.optimizer == o)] if kk is None else K[(K.problem == "P2") & (K.arm == f"{arm}+k{kk}") & (K.optimizer == o)]
            ME[i, j] = (d.E_end.mean() - e0); MP[i, j] = d.P_fb_end.mean()
            n9 = int((d.P_fb_end > 0.9).sum())
            tE[i][j] = f"{d.E_end.mean() - e0:+.2f}"; tP[i][j] = f"{d.P_fb_end.mean():.2f}" + (f"\n({n9}/5 >0.9)" if n9 else "")
            NUM[f"06|P2|{arm}|{o}|{kk}"] = dict(E_end=float(d.E_end.mean()), dE=float(d.E_end.mean() - e0), P=float(d.P_fb_end.mean()), n90=n9, n=len(d))
    lim = 2.0
    im = heat(axes[0, c], np.clip(ME, -lim, lim), k5, colsk, tE, "RdBu_r", -lim, lim, center=0)
    nm = "multi-start" if arm == "spiq_multi" else "single-start"
    axes[0, c].set_title(f"P2, {nm} SPIQ point: energy change vs SPIQ ({e0:.3f}; ground {gs:.3f})\nblue = ended below SPIQ, red = ended above", loc="left", fontsize=9.5, fontweight="bold")
    im2 = heat(axes[1, c], MP, k5, colsk, tP, "Blues", 0, 1)
    axes[1, c].set_title(f"P2, {nm}: P(optimal plan) at the end (SPIQ start: {p0:.2f})", loc="left", fontsize=9.5, fontweight="bold")
fig.suptitle("Round 1 (5 Nov): kick tried on 5 optimizers. Only NFT ends BELOW the SPIQ energy after a kick, and only NFT produces runs above 0.9",
             x=0.01, ha="left", fontsize=11, fontweight="bold")
fig.text(0.01, 0.005, "pilot_round1/runs_P1_P2_P3_b1000(_kicks).jsonl: 1,000 evaluations, 5 trials per cell. Kick = Gaussian noise of that σ (radians) "
         "on every angle of the SPIQ point before optimizing. P1 and P3 not shown: SPIQ is already exact there.", fontsize=8, color="#5b6370")
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
fig.savefig(os.path.join(FIG, "06_round1_which_optimizer_to_kick.png"), bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ 07: round 2 kick-size sweep
R2 = pd.concat([RJ(os.path.join(PILOT, "pilot_round2", f)) for f in ["runs_P2_P4_b1000.jsonl", "runs_P2_b1000_kicks.jsonl", "runs_P2_b3000_kicks.jsonl",
                                                                       "runs_P4_b1000_kicks.jsonl", "runs_P4_b3000_kicks.jsonl"]], ignore_index=True)
R2["base"] = R2.arm.apply(lambda a: a.split("+k")[0])
R2["kick"] = R2.arm.apply(lambda a: float(a.split("+k")[1]) if "+k" in a else 0.0)
cases = [("P2", "spiq_multi", "P2 multi-start SPIQ"), ("P2", "spiq_single", "P2 single-start SPIQ"), ("P4", "spiq_single", "P4 single-start SPIQ")]
fig, axes = plt.subplots(1, 3, figsize=(14, 4.6), sharey=True)
for a, (p, arm, nm) in zip(axes, cases):
    for (opt, b), col, ls, lab in [(("COBYLA", 1000), "#e8743b", "-", "COBYLA, 1,000 evals"), (("NFT", 1000), "#1b5fb5", "-", "NFT, 1,000 evals"),
                                    (("NFT", 3000), "#1b5fb5", "--", "NFT, 3,000 evals")]:
        d = R2[(R2.problem == p) & (R2.base == arm) & (R2.optimizer == opt) & ((R2.budget == b) | ((R2.kick == 0) & (R2.budget == 1000)))]
        if b == 3000:
            d = pd.concat([d[d.kick > 0][d[d.kick > 0].budget == 3000]])
        g = d.groupby("kick").agg(P=("P_fb_end", "mean"), n=("P_fb_end", "size"), n9=("P_fb_end", lambda x: int((x > 0.9).sum()))).reset_index()
        if not len(g):
            continue
        a.plot(g.kick, g.P, marker="o", color=col, ls=ls, lw=2, label=lab)
        for _, r in g.iterrows():
            if r.n9:
                a.annotate(f"{int(r.n9)}/{int(r.n)}", (r.kick, r.P), textcoords="offset points", xytext=(0, -13 if opt == "COBYLA" else 7), ha="center", fontsize=7.5, color=col)
            NUM[f"07|{p}|{arm}|{opt}|{b}|{r.kick}"] = dict(P=float(r.P), n=int(r.n), n90=int(r.n9))
    a.set_title(nm, loc="left", fontsize=10.5, fontweight="bold"); a.set_xlabel("kick size σ (radians; 0 = no kick)")
    a.set_ylim(-0.03, 1.05); a.grid(True, color="#e9ecef")
axes[0].set_ylabel("mean P(optimal plan)"); fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(0.01, 0.925), ncol=3)
fig.suptitle("Round 2 (5 Nov): kick size sweep. Kicked NFT climbs, kicked COBYLA does not; the best σ differs by start. "
             "(Used to choose settings — NOT claimed as evidence)", x=0.01, ha="left", fontsize=10.5, fontweight="bold")
fig.text(0.01, -0.02, "pilot_round2/*.jsonl. P2: 20 trials per point; P4 kicks: 10 trials. Labels = runs above 0.9. "
         "Round 2's P4 σ=1.0 NFT 3,000 result (0.574) did NOT replicate in round 3 (0.276).", fontsize=8, color="#5b6370")
fig.tight_layout(rect=(0, 0, 1, 0.86))
fig.savefig(os.path.join(FIG, "07_round2_kick_size_sweep.png"), bbox_inches="tight"); plt.close(fig)

# ------------------------------------------------------------------ 12: 12 new problems heatmap
NP = RJ(os.path.join(OS, "new_problems", "runs_new.jsonl"))
AB = RJ(os.path.join(OS, "MEETING_2026-11-18", "data", "runs_ablation_random.jsonl"))
AB = AB[AB.kind == "new"]
ids = list(dict.fromkeys(NP.sort_values(["graph", "problem"], key=lambda s: s.map({"CHAIN": 0, "STAR": 1, "CYCLE": 2, "CLIQUE": 3}) if s.name == "graph" else s).problem))
new = NP[NP.formulation == "qaoa"]; old = NP[NP.formulation == "q2o"]
colsN = ["SPIQ start\n(no optimizer)", "SPIQ +\nCOBYLA\n1,000", "SPIQ +\nNFT\n1,000", "SPIQ + kick 1.0\n×5 + NFT\n15,001", "random ×5\n+ NFT (control)\n15,001", "OLD formulation\nkick ×5 + NFT\n15,001"]
M = np.full((len(ids), 6), np.nan); T = [[""] * 6 for _ in ids]
for i, pid in enumerate(ids):
    sp = new[(new.problem == pid) & (new.start == "spiq")]
    vals = [sp.P_fb_start.mean(), sp[sp.pipeline == "COBYLA"].P_fb_end.mean(), sp[sp.pipeline == "NFT"].P_fb_end.mean(),
            sp[sp.pipeline == "kick"].P_fb_end.mean(), AB[AB.problem == pid].P_fb_end.mean(), old[(old.problem == pid) & (old.pipeline == "kick")].P_fb_end.mean()]
    runs = [None, None, None, sp[sp.pipeline == "kick"].P_fb_end.tolist(), AB[AB.problem == pid].P_fb_end.tolist(), old[(old.problem == pid) & (old.pipeline == "kick")].P_fb_end.tolist()]
    for j, v in enumerate(vals):
        M[i, j] = v
        T[i][j] = f"{v:.2f}" + ("" if runs[j] is None else "\n" + " ".join(f"{x:.2f}" for x in runs[j]))
    NUM[f"12|{pid}"] = dict(zip(["spiq", "cobyla", "nft", "kick", "random5", "old_kick"], [float(v) for v in vals]),
                             kick_runs=runs[3], random_runs=runs[4], old_kick_runs=runs[5])
rowsN = [f"{pid.replace('Gen_', '').replace('_4r_s', ' #')}" for pid in ids]
fig, ax = plt.subplots(figsize=(13, 9.2))
im = heat(ax, M, rowsN, colsN, T, "Blues", 0, 1)
ax.axvline(4.5, color="#14181f", lw=1.6)
mean = np.nanmean(M, 0)
ax.set_xticklabels([f"{c}\nmean {m:.3f}" for c, m in zip(colsN, mean)], fontsize=8.3)
NUM["12|mean"] = dict(zip(["spiq", "cobyla", "nft", "kick", "random5", "old_kick"], [float(m) for m in mean]))
fig.colorbar(im, ax=ax, fraction=0.03, label="P(optimal plan)")
ax.set_title("12 new problems: P(optimal plan) per problem and pipeline. Big number = mean; small numbers = each run (3 runs per cell)\n"
             "Left of the black line: new formulation (8 qubits, 64 angles). Right: old formulation (14–20 qubits).", loc="left", fontsize=10.5, fontweight="bold")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "12_new_problems_kick_heatmap.png"), bbox_inches="tight"); plt.close(fig)

json.dump(NUM, open(os.path.join(DATA, "new_figures_numbers.json"), "w"), indent=1, default=float)
print("01 means:", {p: float(np.nanmean(P_end[:, j])) for j, p in enumerate(probs)})
print("12 means:", NUM["12|mean"])
for k in ["06|P2|spiq_multi|NFT|0.5", "06|P2|spiq_single|NFT|0.2", "06|P2|spiq_multi|COBYLA|0.5"]:
    print(k, NUM[k])
