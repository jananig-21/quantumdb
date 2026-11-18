"""Figures for Part 2 (kick + NFT). Reads saved results only; runs no optimizer.

Inputs
  ../data/landscape.json                                   landscape.py (gradient, curvature, kicks, 1-D cuts)
  new_formulation_pilot/pilot_round3/*.jsonl              round 3 (kick size, kicked COBYLA vs kicked NFT, 20 runs)
  new_formulation_pilot/runs_round4.jsonl                 round 4 (6 pipelines, 20 runs)
  new_formulation_pilot/runs_round5_confirm.jsonl         round 5 (kick 1.0 x5 vs COBYLA on unseen seeds)
  new_problems/runs_new.jsonl                             12 new problems
  new_formulation_pilot/runs_p5_kick.jsonl                P5 (if finished)
Output: ../Part2_Kick_plus_NFT/figures/*.png and ../data/part2_numbers.json
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
OS = os.path.join(HERE, "..", "..")
FIG = os.path.join(HERE, "..", "Part2_Kick_plus_NFT", "figures")
os.makedirs(FIG, exist_ok=True)
ORANGE, BLUE, AQUA, GRAY, INK, INK2, MUTED, RED = "#e8673a", "#2a78d6", "#1aa37a", "#b9bec5", "#14181f", "#3b4350", "#7a828c", "#c23b3b"
NAVY, MAIZE = "#00274c", "#ffcb05"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": "#c9ced6", "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": "#e9ecef", "grid.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#fcfcfa", "axes.facecolor": "#fcfcfa",
                     "savefig.dpi": 170})
L = lambda f: [json.loads(l) for l in open(f)]
NUM = {}


def save(fig, name):
    fig.savefig(os.path.join(FIG, name), bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def short(pid):
    return pid.replace("Gen_", "").replace("_4r_s", " #")


# ---------------------------------------------------------------------------------------------------------------
# K1  the pipeline, drawn
fig, ax = plt.subplots(figsize=(13, 4.3))
ax.set_xlim(0, 13); ax.set_ylim(0, 4.3); ax.axis("off")


def box(x, y, w, h, text, fc, ec=None, tc=INK, fs=9, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12", fc=fc, ec=ec or fc, lw=1.4))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc, fontweight="bold" if bold else "normal", linespacing=1.35)


def arrow(x1, y1, x2, y2):
    ax.annotate("", (x2, y2), (x1, y1), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.3, shrinkA=0, shrinkB=0))


box(0.1, 1.55, 1.75, 1.2, "SPIQ point\n(team's code)\nangles x₀", "#e8eef7", NAVY, bold=True)
box(2.3, 3.05, 1.9, 0.95, "Measure x₀ once\n(10,240 shots)\n= candidate 0", "#eef0f3", GRAY)
arrow(1.85, 2.4, 2.3, 3.45)
for k in range(5):
    y = 2.35 - k * 0.52
    box(2.3, y, 1.9, 0.42, f"kick {k + 1}:  x₀ + 1.0·N(0, I)", "#fff6d6", "#e0b400", fs=8.2)
    arrow(1.85, 2.0, 2.3, y + 0.21)
    box(4.6, y, 3.0, 0.42, f"NFT, 3,000 evaluations → candidate {k + 1}", "#e3effc", BLUE, fs=8.2)
    arrow(4.2, y + 0.21, 4.6, y + 0.21)
    arrow(7.6, y + 0.21, 8.0, 1.97)
arrow(4.2, 3.5, 9.1, 2.55)
box(8.0, 1.4, 2.2, 1.15, "Keep the candidate\nwith the LOWEST\nmeasured energy", "#e9f6f0", AQUA, bold=True)
arrow(10.2, 1.97, 10.55, 1.97)
box(10.55, 1.4, 2.35, 1.15, "Score it exactly:\nenergy and\nP(optimal plan)", "#f3f4f6", GRAY)
ax.text(0.1, 0.12, "Total = 1 + 5 × 3,000 = 15,001 evaluations.  The pipeline never sees P(optimal plan) or the ground state while choosing; "
        "it only compares measured energies.", fontsize=8.5, color=MUTED)
ax.text(0.1, 4.15, "Kick + NFT: how one run is orchestrated", fontsize=12, fontweight="bold", color=INK)
save(fig, "K1_pipeline_diagram.png")

# ---------------------------------------------------------------------------------------------------------------
# K2  the SPIQ point is a flat spot: gradient 0, no single angle goes downhill
LS = json.load(open(os.path.join(HERE, "..", "data", "landscape.json")))
keys = list(LS)
lab = [(k.split("|")[0] + (" (team, multi-start)" if "multi" in k else " (team)")) if not k.startswith("new") else short(k.split("|")[1]) for k in keys]
up = np.array([LS[k]["uphill"] / LS[k]["angles"] for k in keys])
fl = np.array([LS[k]["flat"] / LS[k]["angles"] for k in keys])
dn = np.array([LS[k]["downhill"] / LS[k]["angles"] for k in keys])
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.4), gridspec_kw=dict(width_ratios=[1.35, 1]))
y = np.arange(len(keys))[::-1]
a1.barh(y, up, color=ORANGE, height=0.66, label="uphill both ways", edgecolor="white", lw=1)
a1.barh(y, fl, left=up, color=GRAY, height=0.66, label="flat", edgecolor="white", lw=1)
a1.barh(y, dn, left=up + fl, color=AQUA, height=0.66, label="downhill (0 everywhere)", edgecolor="white", lw=1)
for yi, k in zip(y, keys):
    a1.text(1.01, yi, f"{LS[k]['uphill']}/{LS[k]['flat']}/{LS[k]['downhill']}  of {LS[k]['angles']}", va="center", fontsize=7.5, color=INK2)
a1.set_yticks(y); a1.set_yticklabels(lab, fontsize=8); a1.set_xlim(0, 1.28); a1.set_xticks([0, .25, .5, .75, 1]); a1.grid(axis="y", visible=False)
a1.set_xlabel("share of angles: what one angle moving by ±π/2 does to the energy")
a1.legend(frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=3)
g0 = [LS[k]["grad_norm"] for k in keys]; gr = [LS[k]["grad_norm_random"] for k in keys]
a2.scatter(gr, y, color=BLUE, s=34, zorder=3, label="nearby random point (SPIQ + small noise)")
a2.scatter(np.maximum(g0, 1e-16), y, color=RED, marker="D", s=30, zorder=3, label="at the SPIQ point")
for yi, a, b in zip(y, g0, gr):
    a2.plot([max(a, 1e-16), b], [yi, yi], color=GRAY, lw=1, zorder=1)
a2.set_xscale("log"); a2.set_xlim(1e-16, 50); a2.set_yticks(y); a2.set_yticklabels([]); a2.grid(axis="y", visible=False)
a2.set_xlabel("size of the exact gradient |∇E| (log scale)")
a2.legend(frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1)
fig.suptitle("Why COBYLA and NFT stall at SPIQ: the gradient is zero (≈10⁻¹⁴, rounding error) and no single angle goes downhill — on all 17 SPIQ points",
             x=0.01, ha="left", fontsize=11, fontweight="bold", y=1.04)
fig.tight_layout()
save(fig, "K2_spiq_flat_spot.png")
NUM["landscape"] = {k: dict(angles=LS[k]["angles"], grad=LS[k]["grad_norm"], grad_random=LS[k]["grad_norm_random"], up=LS[k]["uphill"],
                            flat=LS[k]["flat"], down=LS[k]["downhill"], sinusoid_resid=LS[k]["sinusoid_max_residual"]) for k in keys}

# ---------------------------------------------------------------------------------------------------------------
# K3  1-D cuts: walk away from SPIQ in random directions; t = kick size
picks = ["P2|spiq_single", "P4|spiq_single", "new|Gen_CLIQUE_4r_s43001"]
fig, axes = plt.subplots(2, 3, figsize=(13, 6.2), sharex=True)
for j, k in enumerate(picks):
    r = LS[k]; t = np.array(r["cut_t"])
    for c in r["cuts"]:
        axes[0, j].plot(t, c["E"], color=BLUE, alpha=0.28, lw=1)
        axes[1, j].plot(t, c["P"], color=BLUE, alpha=0.28, lw=1)
    axes[0, j].plot(t, np.mean([c["E"] for c in r["cuts"]], 0), color=NAVY, lw=2.2)
    axes[0, j].axhline(r["E_ground"], color=INK, ls="--", lw=1); axes[0, j].text(3, r["E_ground"], "ground state ", ha="right", va="bottom", fontsize=7.5)
    axes[0, j].axhline(r["E_spiq"], color=ORANGE, lw=1.2); axes[0, j].text(3, r["E_spiq"], "SPIQ energy ", ha="right", va="bottom", fontsize=7.5, color=ORANGE)
    axes[1, j].axhline(r["P_spiq"], color=ORANGE, lw=1.2)
    for a in axes[:, j]:
        a.axvline(1.0, color=MAIZE, lw=6, alpha=0.35, zorder=0)
    nm = k.split("|")[0] + (" (single SPIQ)" if "single" in k else "") if not k.startswith("new") else short(k.split("|")[1]) + " (new problem)"
    axes[0, j].set_title(f"{nm}: {r['angles']} angles", loc="left", fontsize=10, fontweight="bold")
    axes[1, j].set_xlabel("distance t along a random direction  (x₀ + t·d,  d ~ N(0, I))")
    axes[1, j].set_ylim(-0.03, 1.03)
axes[0, 0].set_ylabel("exact energy"); axes[1, 0].set_ylabel("P(optimal plan)")
fig.suptitle("Walking away from the SPIQ point in 12 random directions (thin lines; thick = mean), no optimizer yet. Yellow band = kick size 1.0",
             x=0.01, ha="left", fontsize=11, fontweight="bold")
fig.tight_layout()
save(fig, "K3_landscape_cuts.png")

# ---------------------------------------------------------------------------------------------------------------
# K4  what a kick alone does (before NFT): energy and P(opt) at kicked points
fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
for ax, k in zip(axes, picks):
    r = LS[k]
    data = [r["kick"][s]["E"] for s in ("0.2", "0.5", "1.0", "2.0")]
    bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False, medianprops=dict(color=INK))
    for p, c in zip(bp["boxes"], ["#cfe0f6", "#a8c8f0", MAIZE, "#f2b8a2"]):
        p.set_facecolor(c); p.set_edgecolor(INK2)
    ax.axhline(r["E_spiq"], color=ORANGE, lw=1.2); ax.axhline(r["E_ground"], color=INK, ls="--", lw=1)
    ax.set_xticklabels(["σ 0.2", "σ 0.5", "σ 1.0", "σ 2.0"])
    nm = k.split("|")[0] if not k.startswith("new") else short(k.split("|")[1])
    ax.set_title(f"{nm}", loc="left", fontsize=10, fontweight="bold")
    NUM.setdefault("kick_alone", {})[k] = {s: dict(E_median=float(np.median(r["kick"][s]["E"])), P_mean=float(np.mean(r["kick"][s]["P"])),
                                                    share_below_spiq=float(np.mean(np.array(r["kick"][s]["E"]) < r["E_spiq"])))
                                            for s in r["kick"]}
axes[0].set_ylabel("exact energy right after the kick")
fig.suptitle("A kick alone makes the energy WORSE (orange = SPIQ, dashed = ground). It is a jump to a new starting point, not an improvement; NFT does the improving.",
             x=0.01, ha="left", fontsize=10.5, fontweight="bold")
fig.tight_layout()
save(fig, "K4_kick_alone.png")

# ---------------------------------------------------------------------------------------------------------------
# K5  round 3: kick size, and kick needs NFT (kicked COBYLA does not work)
R3 = L(os.path.join(OS, "new_formulation_pilot", "pilot_round3", "runs_P2_P4_b3000_t100.jsonl")) + \
     L(os.path.join(OS, "new_formulation_pilot", "pilot_round3", "runs_P2_P4_b3000_kicks_t100.jsonl"))
cases = [("P2", "spiq_multi"), ("P2", "spiq_single"), ("P4", "spiq_single")]
fig, axes = plt.subplots(1, 3, figsize=(13, 3.9), sharey=True)
r3 = {}
for ax, (p, a) in zip(axes, cases):
    x = np.arange(3)
    for k, (opt, col) in enumerate([("COBYLA", ORANGE), ("NFT", BLUE)]):
        vals, hits = [], []
        for suf in ("", "+k0.5", "+k1.0"):
            rs = [r for r in R3 if r["problem"] == p and r["arm"] == a + suf and r["optimizer"] == opt]
            v = [r["P_fb_end"] for r in rs]
            vals.append(np.mean(v)); hits.append(sum(np.array(v) > 0.9))
            r3[f"{p}|{a}{suf}|{opt}"] = dict(P=float(np.mean(v)), over09=int(hits[-1]), n=len(v), E=float(np.mean([r["E_end"] for r in rs])))
            xs = x[len(vals) - 1] + (-0.19 if k == 0 else 0.19) + np.random.default_rng(len(vals) + 7 * k).uniform(-0.07, 0.07, len(v))
            ax.scatter(xs, v, s=9, color=col, alpha=0.45, zorder=3, lw=0)
        ax.bar(x + (-0.19 if k == 0 else 0.19), vals, width=0.36, color=col, alpha=0.35, zorder=2, label=opt)
        for xi, v_, h in zip(x + (-0.19 if k == 0 else 0.19), vals, hits):
            ax.text(xi, min(v_ + 0.04, 1.02), f"{v_:.2f}", ha="center", fontsize=7.5, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(["no kick", "kick 0.5", "kick 1.0"]); ax.set_ylim(-0.03, 1.1)
    ax.set_title(f"{p} · {'multi-start' if 'multi' in a else 'single-start'} SPIQ", loc="left", fontsize=10, fontweight="bold")
axes[0].set_ylabel("P(optimal plan) at the end\n(bar = mean, dots = 20 runs)"); axes[0].legend(frameon=False, fontsize=8, loc="upper left")
fig.suptitle("Round 3 (one kick, 3,000 evaluations, 20 runs): a kick only pays off when NFT follows it: kicked COBYLA never gets a run above 0.9 (0 of 20 in every cell)",
             x=0.01, ha="left", fontsize=11, fontweight="bold")
fig.tight_layout()
save(fig, "K5_round3_kick_needs_NFT.png")
NUM["round3"] = r3

# ---------------------------------------------------------------------------------------------------------------
# K6  round 4: six pipelines with the same 15,000 budget
R4 = L(os.path.join(OS, "new_formulation_pilot", "runs_round4.jsonl"))
pipes = [("cobyla", "COBYLA"), ("nft", "NFT"), ("fixed0.5", "kick 0.5 ×5"), ("fixed1.0", "kick 1.0 ×5"), ("portfolio", "mixed kicks ×5"), ("halving", "halving")]
fig, axes = plt.subplots(1, 3, figsize=(13, 4.1), sharey=True)
r4 = {}
for ax, (p, a) in zip(axes, cases):
    for i, (pp, nm) in enumerate(pipes):
        v = np.array([r["P_fb_end"] for r in R4 if r["problem"] == p and r["arm"] == a and r["pipeline"] == pp])
        col = AQUA if pp == "fixed1.0" else (ORANGE if pp == "cobyla" else (RED if pp == "halving" else BLUE))
        ax.bar(i, v.mean(), color=col, alpha=0.35 if pp != "fixed1.0" else 0.6, width=0.7, zorder=2)
        ax.scatter(i + np.random.default_rng(i).uniform(-0.2, 0.2, len(v)), v, s=9, color=col, alpha=0.6, lw=0, zorder=3)
        ax.text(i, min(v.mean() + 0.04, 1.03), f"{v.mean():.2f}", ha="center", fontsize=7.5, color=INK2)
        r4[f"{p}|{a}|{pp}"] = dict(P=float(v.mean()), over09=int((v > 0.9).sum()), n=len(v))
    ax.set_xticks(range(len(pipes))); ax.set_xticklabels([n for _, n in pipes], rotation=30, ha="right", fontsize=8)
    ax.set_title(f"{p} · {'multi-start' if 'multi' in a else 'single-start'} SPIQ", loc="left", fontsize=10, fontweight="bold")
    ax.set_ylim(-0.03, 1.1)
axes[0].set_ylabel("P(optimal plan) at the end\n(bar = mean, dots = 20 runs)")
fig.suptitle("Round 4: six ways to spend the same 15,000 evaluations (plan fixed in advance). Kick 1.0 ×5 (green) is the only one at the top on all three",
             x=0.01, ha="left", fontsize=11, fontweight="bold")
fig.tight_layout()
save(fig, "K6_round4_six_pipelines.png")
NUM["round4"] = r4

# ---------------------------------------------------------------------------------------------------------------
# K7  why halving failed: which kick it picked
fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
KS = [0.2, 0.35, 0.5, 0.75, 1.0]
hv = {}
for ax, (p, a) in zip(axes, [("P2", "spiq_single"), ("P4", "spiq_single")]):
    rs = [r for r in R4 if r["problem"] == p and r["arm"] == a and r["pipeline"] == "halving"]
    cnt = [sum(abs(r["picked_kick"] - k) < 1e-9 for r in rs) for k in KS]
    ax.bar(range(5), cnt, color=[RED if k == 0.2 else GRAY for k in KS], width=0.65, zorder=2)
    for i, c in enumerate(cnt):
        ax.text(i, c + 0.3, str(c), ha="center", fontsize=8.5, color=INK2)
    ax.set_xticks(range(5)); ax.set_xticklabels([f"kick {k}" for k in KS]); ax.set_ylim(0, 21)
    ax.set_title(f"{p} single-start: kick size behind the final answer (20 runs)", loc="left", fontsize=9.5, fontweight="bold")
    hv[p] = dict(zip(map(str, KS), map(int, cnt)))
axes[0].set_ylabel("runs")
fig.suptitle("Halving screened 15 kicks for 400 evaluations each and kept the lowest: small kicks look best early, so it picked kick 0.2 — and lost on P4",
             x=0.01, ha="left", fontsize=10.5, fontweight="bold")
fig.tight_layout()
save(fig, "K7_halving_failure.png")
NUM["halving_picked"] = hv

# ---------------------------------------------------------------------------------------------------------------
# K8  round 5: confirmation on unseen seeds, per run, P(opt) and energy
R5 = L(os.path.join(OS, "new_formulation_pilot", "runs_round5_confirm.jsonl"))
fig, axes = plt.subplots(2, 3, figsize=(13, 6.4))
r5 = {}
for j, (p, a) in enumerate(cases):
    for i, (pp, col, nm) in enumerate([("cobyla", ORANGE, "COBYLA (team)"), ("fixed1.0", AQUA, "kick 1.0 ×5 + NFT")]):
        rs = [r for r in R5 if r["problem"] == p and r["arm"] == a and r["pipeline"] == pp]
        P = np.array([r["P_fb_end"] for r in rs]); E = np.array([r["E_end"] for r in rs])
        jit = np.random.default_rng(i + 3 * j).uniform(-0.18, 0.18, len(rs))
        axes[0, j].scatter(i + jit, P, color=col, s=22, alpha=0.8, lw=0, zorder=3)
        axes[0, j].hlines(P.mean(), i - 0.3, i + 0.3, color=INK, lw=2, zorder=4)
        axes[0, j].text(i + 0.33, P.mean(), f"{P.mean():.3f}", va="center", fontsize=8, color=INK)
        axes[1, j].scatter(i + jit, E, color=col, s=22, alpha=0.8, lw=0, zorder=3)
        axes[1, j].hlines(E.mean(), i - 0.3, i + 0.3, color=INK, lw=2, zorder=4)
        r5[f"{p}|{a}|{pp}"] = dict(P=float(P.mean()), over09=int((P > 0.9).sum()), n=len(P), E=float(E.mean()),
                                   evals=float(np.mean([r["evals"] for r in rs])), E_ground=rs[0]["E_ground"], E_start=rs[0]["E_start"], P_start=rs[0]["P_fb_start"])
    eg = [r["E_ground"] for r in R5 if r["problem"] == p][0]
    axes[1, j].axhline(eg, color=INK, ls="--", lw=1); axes[1, j].text(1.45, eg, "ground", va="bottom", ha="right", fontsize=7.5)
    for a_ in axes[:, j]:
        a_.set_xticks([0, 1]); a_.set_xticklabels(["COBYLA (team)", "kick 1.0 ×5 + NFT"]); a_.set_xlim(-0.5, 1.6)
    axes[0, j].set_ylim(-0.03, 1.05)
    axes[0, j].set_title(f"{p} · {'multi-start' if 'multi' in a else 'single-start'} SPIQ", loc="left", fontsize=10, fontweight="bold")
axes[0, 0].set_ylabel("P(optimal plan), 20 runs"); axes[1, 0].set_ylabel("exact energy, 20 runs")
fig.suptitle("Round 5, the claim: kick 1.0 ×5 + NFT vs COBYLA on 20 unseen seeds (black bar = mean)", x=0.01, ha="left", fontsize=11, fontweight="bold")
fig.tight_layout()
save(fig, "K8_round5_confirmation.png")
NUM["round5"] = r5

# ---------------------------------------------------------------------------------------------------------------
# K9  12 new problems: every kick + NFT run, P(opt) over the 15,000 evaluations (new formulation)
NP = L(os.path.join(OS, "new_problems", "runs_new.jsonl"))
pids = sorted({r["problem"] for r in NP}, key=lambda p: (["CHAIN", "STAR", "CYCLE", "CLIQUE"].index(p.split("_")[1]), p))
fig, axes = plt.subplots(3, 4, figsize=(13, 7.4), sharex=True, sharey=True)
newk = {}
for ax, pid in zip(axes.flat, pids):
    rs = [r for r in NP if r["problem"] == pid and r["formulation"] == "qaoa" and r["pipeline"] == "kick"]
    for r in rs:
        xs = [0] + sorted(int(k) for k in r["trace"]); d = {int(k): v[1] for k, v in r["trace"].items()}
        ax.plot(xs, [r["P_fb_start"]] + [d[x] for x in xs[1:]], color=AQUA, lw=1.6, alpha=0.85)
    for pp, col in (("COBYLA", ORANGE), ("NFT", BLUE)):
        v = [r["P_fb_end"] for r in NP if r["problem"] == pid and r["formulation"] == "qaoa" and r["start"] == "spiq" and r["pipeline"] == pp]
        ax.scatter([1000], [np.mean(v)], color=col, marker="D", s=26, zorder=4)
    for k in range(1, 5):
        ax.axvline(1 + 3000 * k, color=GRAY, lw=0.7, ls=":")
    ax.set_title(short(pid), loc="left", fontsize=9, fontweight="bold"); ax.set_ylim(-0.03, 1.05)
    newk[pid] = [round(r["P_fb_end"], 3) for r in rs]
for a in axes[-1]:
    a.set_xlabel("evaluations")
for a in axes[:, 0]:
    a.set_ylabel("P(optimal plan)")
h = [plt.Line2D([], [], color=AQUA, lw=2), plt.Line2D([], [], color=ORANGE, marker="D", ls=""), plt.Line2D([], [], color=BLUE, marker="D", ls="")]
fig.legend(h, ["kick + NFT, each of 3 runs (best-so-far point)", "COBYLA from SPIQ at 1,000 evals (mean)", "NFT from SPIQ at 1,000 evals (mean)"],
           loc="upper left", bbox_to_anchor=(0.005, 0.975), ncol=3, frameon=False, fontsize=8.5)
fig.suptitle("12 new problems, new formulation: every kick + NFT run (dotted = next kick)", x=0.01, ha="left", fontsize=11, fontweight="bold", y=1.0)
fig.tight_layout(rect=(0, 0, 1, 0.95))
save(fig, "K9_new_problems_every_run.png")
NUM["new_problems_kick_runs"] = newk

# ---------------------------------------------------------------------------------------------------------------
# K10  P5 (24 qubits) progress per restart, if any run has finished
P5 = os.path.join(OS, "new_formulation_pilot", "runs_p5_kick.jsonl")
p5 = L(P5) if os.path.exists(P5) else []
if p5:
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    for r in p5:
        ks = [s["restart"] for s in r["segments"]]
        axes[0].plot(ks, [s["E_best"] for s in r["segments"]], marker="o", lw=2, label=f"run {r['rep']}")
        axes[1].plot(ks, [s["P_fb_best"] for s in r["segments"]], marker="o", lw=2, label=f"run {r['rep']}")
    axes[0].axhline(p5[0]["E_ground"], color=INK, ls="--", lw=1)
    axes[0].set_xlabel("restart (0 = SPIQ start)"); axes[0].set_ylabel("exact energy, best candidate")
    axes[1].set_xlabel("restart (0 = SPIQ start)"); axes[1].set_ylabel("P(optimal plan), best candidate"); axes[1].set_ylim(-0.03, 1.03)
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle("P5 (24 qubits, 252 angles): kick 1.0 ×5 + NFT, progress per restart", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    save(fig, "K10_P5_progress.png")
    NUM["P5"] = [dict(rep=r["rep"], E_start=r["E_start"], E_end=r["E_end"], E_ground=r["E_ground"], P_start=r["P_fb_start"], P_end=r["P_fb_end"],
                      hours=r["wall_s"] / 3600) for r in p5]

# ---------------------------------------------------------------------------------------------------------------
# K10b  P5 interim, parsed from the log (used until runs_p5_kick.jsonl has finished runs)
import re, csv
if not p5:
    lines = open(os.path.join(OS, "new_formulation_pilot", "log_p5_kick.txt")).read().splitlines()
    pat = re.compile(r"rep (\d+) restart (\d+): E ([\d.]+) fb ([\d.]+) \| best E ([\d.]+) fb ([\d.]+)")
    T2 = list(csv.DictReader(open(os.path.join(OS, "comparison_2x2", "tables", "table_2x2_all.csv"))))
    row = [t for t in T2 if t["formulation"] == "qaoa" and t["start"] == "spiq_single" and t["problem"] == "P5"][0]
    E0, P0, EG = float(row["E_start"]), float(row["P_fb_start"]), float(row["E_ground"])
    runs = {}
    for l in lines:
        m = pat.search(l)
        if m:
            rep, k = int(m.group(1)), int(m.group(2))
            runs.setdefault(rep, [(0, E0, P0)]).append((k, float(m.group(5)), float(m.group(6))))
    if runs:
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
        for rep in sorted(runs):
            ks, Es, Ps = zip(*runs[rep])
            axes[0].plot(ks, Es, marker="o", lw=2, label=f"run {rep}")
            axes[1].plot(ks, Ps, marker="o", lw=2, label=f"run {rep}")
        axes[0].axhline(EG, color=INK, ls="--", lw=1); axes[0].text(5, EG, "ground state ", ha="right", va="bottom", fontsize=8)
        for a in axes:
            a.set_xlim(-0.3, 5.3); a.set_xticks(range(6)); a.set_xlabel("restart (0 = SPIQ start; 1–5 = kicks)")
        axes[0].set_ylabel("exact energy, best candidate so far"); axes[1].set_ylabel("P(optimal plan), best candidate so far")
        axes[1].set_ylim(-0.03, 1.03); axes[1].legend(frameon=False, fontsize=8)
        axes[1].axhline(0, color=ORANGE, lw=1.2); axes[1].text(5.2, 0.02, "every optimizer from SPIQ (Part 1)", color=ORANGE, ha="right", fontsize=7.5)
        done_k = max(k for r in runs.values() for k, _, _ in r)
        fig.suptitle(f"P5 (24 qubits, 252 angles), INTERIM: kick 1.0 + NFT after {done_k} of 5 kicks", x=0.01, ha="left", fontsize=11, fontweight="bold")
        fig.tight_layout()
        save(fig, "K10_P5_interim.png")
        NUM["P5_interim"] = {str(k): v for k, v in runs.items()}

# ---------------------------------------------------------------------------------------------------------------
# K11  control: 5 random restarts + NFT vs 5 kicks (1.0) from SPIQ + NFT, same budget, same selection
AB = os.path.join(HERE, "..", "data", "runs_ablation_random.jsonl")
if os.path.exists(AB):
    ab = L(AB)
    rows = []
    for p, a, nm in [("P2", "spiq_multi", "P2 multi-start"), ("P2", "spiq_single", "P2 single-start"), ("P4", "spiq_single", "P4 single-start")]:
        rnd = [r["P_fb_end"] for r in ab if r["problem"] == p and r["arm"] == a]
        kk = [r["P_fb_end"] for r in R5 if r["problem"] == p and r["arm"] == a and r["pipeline"] == "fixed1.0"]
        cb = [r["P_fb_end"] for r in R5 if r["problem"] == p and r["arm"] == a and r["pipeline"] == "cobyla"]
        if rnd:
            rows.append((nm, cb, kk, rnd))
    rnd = [r["P_fb_end"] for r in ab if r["kind"] == "new"]
    kk = [r["P_fb_end"] for r in NP if r["formulation"] == "qaoa" and r["pipeline"] == "kick"]
    cb = [r["P_fb_end"] for r in NP if r["formulation"] == "qaoa" and r["pipeline"] == "COBYLA" and r["start"] == "spiq"]
    if rnd:
        rows.append(("12 new problems", cb, kk, rnd))
    fig, axes = plt.subplots(1, len(rows), figsize=(3.4 * len(rows), 3.9), sharey=True, squeeze=False)
    ctl = {}
    for ax, (nm, cb, kk, rn) in zip(axes[0], rows):
        for i, (v, col, lab_) in enumerate([(cb, ORANGE, "COBYLA\nfrom SPIQ"), (kk, AQUA, "kick 1.0 ×5\n+ NFT"), (rn, "#7d5ba6", "random ×5\n+ NFT")]):
            v = np.array(v)
            ax.bar(i, v.mean(), color=col, alpha=0.45, width=0.7, zorder=2)
            ax.scatter(i + np.random.default_rng(i).uniform(-0.2, 0.2, len(v)), v, s=10, color=col, lw=0, alpha=0.8, zorder=3)
            ax.text(i, min(v.mean() + 0.04, 1.03), f"{v.mean():.3f}\n({(v > 0.9).sum()}/{len(v)} > 0.9)", ha="center", fontsize=7, color=INK2)
        ax.set_xticks(range(3)); ax.set_xticklabels(["COBYLA\nfrom SPIQ", "kick 1.0 ×5\n+ NFT", "random ×5\n+ NFT"], fontsize=8)
        ax.set_title(nm, loc="left", fontsize=10, fontweight="bold"); ax.set_ylim(-0.03, 1.15)
        ctl[nm] = dict(cobyla=float(np.mean(cb)), kick=float(np.mean(kk)), random=float(np.mean(rn)), n_random=len(rn), n_kick=len(kk),
                       kick_over09=int((np.array(kk) > 0.9).sum()), random_over09=int((np.array(rn) > 0.9).sum()))
    axes[0, 0].set_ylabel("P(optimal plan) at the end")
    fig.suptitle("Control: is it the kick from SPIQ, or just restarts + NFT? Same budget (15,001), same selection (lowest measured energy)",
                 x=0.01, ha="left", fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    save(fig, "K11_random_restart_control.png")
    NUM["control_random"] = ctl

json.dump(NUM, open(os.path.join(HERE, "..", "data", "part2_numbers.json"), "w"), indent=1, default=float)
print("numbers -> data/part2_numbers.json")
