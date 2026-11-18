# Part 2: Kick + NFT — how it is orchestrated, why it has to be built that way, and what it achieved

**The question this part answers.** Part 1 showed that from the team's SPIQ starting point, COBYLA never improves anything and NFT improves only slowly. This part explains **why** the optimizers get stuck there, using exact measurements of the energy landscape. It then explains **how** kick + NFT gets around that, **why each piece of the pipeline is there**, how the design was chosen and tested over five rounds, and **what it achieved** on P2, P4, 12 new problems and (in progress) P5.

**Short answer.** At every SPIQ point we measured, the energy's slope is exactly zero, and turning any single angle either raises the energy or leaves it flat. So an optimizer that follows slopes, or tries one direction at a time, has nowhere to go. Kick + NFT gets out of this by repeatedly *jumping* away from the SPIQ point (a random "kick" of 1 radian on every angle), then letting NFT, an optimizer built for exactly this kind of circuit, descend from the new spot. Five jumps, and the lowest measured energy wins. On unseen seeds, this raised P(optimal plan) from 0.50 → 0.977 (P2 multi-start), 0.00 → 0.973 (P2 single-start), 0.25 → 0.693 (P4), and on 12 new problems from 0.365 → 0.904 on average. A control experiment (Section 2.7) then showed that **fully random restarts + NFT do just as well** (0.999 / 0.975 / 0.712 / 0.963). The things that matter are leaving the flat spot, NFT, and several restarts chosen by measured energy; starting the restarts near SPIQ's angles does not.

All figure files named here are in `figures/` next to this document.

---

## 2.1 The pipeline, step by step

**Figure: `K1_pipeline_diagram.png`.** One run of the pipeline, as a flow chart.

One run, start to finish:

1. **Start from the team's SPIQ point** x₀: its angles, one per gate (64 for P2, 140 for P4, 252 for P5).
2. **Measure x₀ once** with 10,240 shots. This measured energy is *candidate 0*. Keeping it means the pipeline can never return something with a worse measured energy than SPIQ itself.
3. **Repeat 5 times:**
   a. **Kick:** make a new starting point x₀ + 1.0 · N(0, I). Every angle gets its own random nudge drawn from a bell curve with standard deviation 1.0 radian (about 57°).
   b. **NFT for 3,000 evaluations** from that kicked point. Remember the lowest measured energy it saw, and where. That is candidate 1, 2, 3, 4 or 5.
4. **Return the candidate with the lowest measured energy.** The pipeline chooses using only the 10,240-shot energy it measured. It **never sees P(optimal plan) or the ground state**, so the choice is fair and could be done on real hardware.
5. **Score the returned angles exactly:** exact energy and exact P(optimal plan). This scoring is only for our report; the pipeline itself never uses it.

Total cost: 1 + 5 × 3,000 = **15,001 evaluations** per run.

📘 **New word: kick (perturbation).** Adding random noise to every angle at once. Here: x_kicked = x₀ + σ · z, where each zᵢ is an independent standard normal number and **σ = 1.0 radian**. About 68% of angles move by less than 1 radian, about 95% by less than 2.

📘 **New word: restart.** Running the optimizer again from a different starting point and keeping the best result. Five kicks = five restarts.

📘 **New word: measured energy vs exact energy.** *Measured* = the average energy of 10,240 simulated shots: what a real quantum computer would report, with shot noise. *Exact* = Σ probability × energy, computed from the full state with no noise. The pipeline chooses with measured energy; we report exact energy.

---

## 2.2 NFT: the optimizer built for these circuits

📘 **New word: NFT (Nakanishi–Fujii–Todo, 2019, "Sequential minimal optimization for quantum-classical hybrid algorithms").** An optimizer that changes **one angle at a time** and, for that angle, jumps straight to the **exact best value** instead of taking small steps.

### Why one angle can be solved exactly

In SPIQ's circuit, every angle drives exactly one rotation gate. For such a circuit, if you hold all other angles fixed and turn one angle θ, the energy traces out a perfect sine wave:

  E(θ) = a + b·cos θ + c·sin θ

Three numbers (a, b, c) describe the whole curve, so **three energy values pin it down**. Once you know the curve, its lowest point is known exactly.

**We checked this directly** on all 17 SPIQ points (P1–P4 and the 12 new problems). For 4 random angles per point, we evaluated the energy at 5 different values of the angle and fitted a + b cos θ + c sin θ. The largest error of the fit was **3 × 10⁻¹⁵ to 2.5 × 10⁻¹⁴**, which is floating-point rounding. The energy is an exact sine wave in each angle. (`data/landscape.json`, field `sinusoid_max_residual`.)

### What NFT does at each step (qiskit's implementation)

For angle i (cycling through all angles in order):

1. z0 = energy at the current angles (reused from the previous step; re-measured every 32 steps to stop errors piling up)
2. z1 = energy with angle i turned by **+π/2**
3. z3 = energy with angle i turned by **−π/2**
4. From z0, z1, z3, compute the sine wave's amplitude and phase, and **set angle i to the bottom of the wave**.

So each step costs **2 evaluations** (3 every 32nd step), and it moves the angle to its exact best value given the other angles.

> **Worked example of budget.** 3,000 evaluations ≈ 1,480 single-angle steps.
> P2 / new problems (64 angles): ≈ 23 full sweeps over all angles.
> P4 (140 angles): ≈ 10.5 sweeps. P5 (252 angles): ≈ 6 sweeps.

### Why NFT and not the others

- **No step size to tune.** It jumps to the exact 1-D minimum, so it can make big moves when the landscape allows and never overshoots along that angle.
- **Robust to shot noise.** It uses evaluations ±π/2 apart (big differences in energy), not tiny finite-difference steps that shot noise swamps. Gradient methods (L-BFGS-B, SLSQP, ADAM) estimate slopes from nearby points and suffer from noise.
- **It fits the circuit.** The sine-wave property is exactly the structure of SPIQ's per-gate circuit.
- **Evidence:** in Part 1, NFT is the only one of 10 optimizers that kept lowering the energy from SPIQ (B1, B4), and it had the best average P(optimal plan) from the standard start in both formulations.

---

## 2.3 Why every optimizer stalls at the SPIQ point

**Figure: `K2_spiq_flat_spot.png`.** Left: for each of the 17 SPIQ points, the share of angles where turning that one angle by ±π/2 goes uphill both ways (orange), stays flat (gray), or goes downhill (green). Right: the size of the exact gradient at the SPIQ point (red diamond) vs at a nearby point with small random noise added (blue dot), on a log scale.

📘 **New word: gradient.** The slope of the energy in every direction at once: how fast the energy changes if you nudge each angle. An optimizer that follows slopes (gradient descent, COBYLA's linear model) needs a nonzero slope to know where to go.

📘 **New word: parameter-shift rule.** Because the energy is a sine wave in each angle (Section 2.2), its exact slope is

  slope_i = [E(θ + π/2 in angle i) − E(θ − π/2 in angle i)] / 2

No approximation, no tiny steps.

📘 **New word: curvature (single-angle).** How the energy bends along one angle:

  curvature_i = [E(+π/2) + E(−π/2) − 2·E(θ)] / 2

Positive = both directions go up (a valley floor along that angle). Zero = flat. Negative = it goes down in at least one direction.

### What we measured (exact, no shots)

| SPIQ point | Angles | Slope size at SPIQ | Slope size nearby | Uphill / flat / **downhill** angles |
|---|---|---|---|---|
| P1 (team) | 18 | 1.4 × 10⁻¹⁵ | 1.26 | 15 / 3 / **0** |
| P2 (team, single-start) | 64 | 4.6 × 10⁻¹⁵ | 2.98 | 51 / 13 / **0** |
| P2 (team, multi-start) | 64 | 4.5 × 10⁻¹⁵ | 2.49 | 56 / 8 / **0** |
| P3 (team) | 64 | 6.0 × 10⁻¹⁵ | 5.03 | 59 / 5 / **0** |
| P4 (team) | 140 | 1.1 × 10⁻¹⁴ | 5.06 | 119 / 21 / **0** |
| 12 new problems | 64 each | 8.7 × 10⁻¹⁵ – 1.8 × 10⁻¹⁴ | 6.65 – 12.33 | 48–60 / 4–16 / **0** |

**On all 17 SPIQ points, the slope is zero (to rounding error) and not one angle goes downhill.** Move a little away from SPIQ and the slope is suddenly 1–12, completely normal.

### What this means for each optimizer

- **COBYLA** builds a linear model of the energy from nearby points. The model's slope is zero (plus shot noise), so it has no direction to follow. It shrinks its trust region and stops where it started. This is exactly Part 1's "COBYLA ends where SPIQ started, every time".
- **Gradient methods** see a zero gradient (plus noise). Same outcome.
- **NFT** looks along one angle at a time. For "uphill" angles, the bottom of the sine wave *is* the current value, so it stays put. For "flat" angles the wave has (almost) no height, so where it lands is decided by shot noise. That is why NFT drifts slowly away from SPIQ and sometimes finds a way down (old P3: 40% of the gap in 1,000 evaluations), but slowly.
- **SPSA** moves all angles at once by a random step. From here, almost every random direction goes uphill, so it ends worse (Part 1: gap closed −0.5 to −1.4).

📘 **New word: saddle point / flat spot.** A point where the slope is zero but it is not the lowest point around. Being uphill along every *single* angle does not mean uphill along every *combination* of angles. The SPIQ point can still have downhill directions that need several angles to move together, which one-angle-at-a-time and slope-following methods cannot find from here.

**Why SPIQ points look like this.** SPIQ only searches Clifford angle settings (multiples of π/2) and keeps the best one it finds. Turning any one angle by ±π/2 gives another Clifford setting that SPIQ could have chosen, so it is no surprise that those moves are uphill or flat. The measured slope of exactly zero at all 17 points is an empirical finding, not something we assumed.

---

## 2.4 What the landscape looks like further away

**Figure: `K3_landscape_cuts.png`.** Starting at the SPIQ point, walk in 12 random directions (x₀ + t·d, where d is a fresh N(0, I) vector, so t matches the kick size σ). Top row: exact energy. Bottom row: P(optimal plan). The yellow band marks t = 1.0, the kick size used. Shown for P2, P4 and a new clique problem.

- **Energy rises in every direction.** The SPIQ point sits at the bottom of a narrow dip. By t ≈ 1 the energy reaches a high plateau: around 10 on P2 (SPIQ 5.1, ground 3.65), 24 on P4 (SPIQ 7.4, ground 4.9), 32 on the clique (SPIQ 13.8, ground 8.9).
- **P(optimal plan) along these straight lines stays low (≈ 0.05–0.2).** No straight path out of SPIQ leads directly to a good plan. The better regions have to be *found* by an optimizer, not stumbled on.

**Figure: `K4_kick_alone.png`.** The exact energy right after a kick, before NFT runs (200 kicks per size). Orange line = SPIQ energy, dashed = ground state.

| Kick size σ | P2: median energy after kick (SPIQ 5.12) | P4 (SPIQ 7.42) | Clique #43001 (SPIQ 13.78) | Share of kicks below SPIQ |
|---|---|---|---|---|
| 0.2 | 6.23 | 11.74 | 18.07 | 0% |
| 0.5 | 8.88 | 20.87 | 27.84 | 0% |
| 1.0 | 10.39 | 24.00 | 32.53 | 0% |
| 2.0 | 10.47 | 24.18 | 32.80 | 0% |

**A kick on its own always makes things worse: 0% of 200 kicks land below SPIQ, at every size.** The kick is not an improvement step. It is a jump to a new starting point, and NFT does all the improving afterwards.

**Why σ = 1.0 matters here.** The energy after a 1.0 kick is almost the same as after a 2.0 kick: by σ = 1 the circuit has left SPIQ's dip completely. A smaller kick (0.2, 0.5) leaves the circuit on the slope of the dip, where NFT tends to slide straight back to SPIQ's point. Section 2.5 shows this in the results.

---

## 2.5 How the design was chosen: five rounds, with one pre-registered failure

The pipeline was not designed in one go. It was built in rounds, and the final claim is only made on seeds that were never used to choose anything.

All rounds: new formulation, the team's saved SPIQ points (`formulation` branch), 10,240 shots per evaluation, depth 2, exact scoring of the best-seen point.

### Rounds 1–2: exploring (not used as evidence)

- **Round 1:** the 10 optimizers from SPIQ on P1–P3, 5 trials, 1,000 evaluations. Everyone stuck, as in Part 1.
- **Round 2:** tried kicks of several sizes on P2 and P4 to choose candidates. **Because these runs were used to pick settings, their numbers are not claimed.**

### Round 3: one kick, then 3,000 evaluations, 20 new seeds each

**Figure: `K5_round3_kick_needs_NFT.png`.** Bars = mean P(optimal plan), dots = each of 20 runs. COBYLA (orange) vs NFT (blue), with no kick, kick 0.5, kick 1.0.

| Start | COBYLA | NFT | COBYLA + kick 0.5 | NFT + kick 0.5 | COBYLA + kick 1.0 | NFT + kick 1.0 |
|---|---|---|---|---|---|---|
| P2 multi-start | 0.500 (0/20 > 0.9) | 0.501 (0) | 0.334 (0) | **0.950 (18)** | 0.275 (0) | 0.531 (8) |
| P2 single-start | 0.000 (0) | 0.000 (0) | 0.065 (0) | 0.142 (0) | 0.187 (0) | **0.653 (12)** |
| P4 single-start | 0.250 (0) | 0.296 (0) | 0.164 (0) | 0.229 (2) | 0.188 (0) | 0.276 (2) |

What round 3 taught:

1. **A kick needs NFT after it.** Kicked COBYLA never produced a single run above 0.9 (0 of 20 in every cell), and on P2 multi-start a kick made COBYLA *worse* (0.50 → 0.33 / 0.28). Kicked NFT produced 18 of 20 and 12 of 20.
2. **The best kick size depended on the start:** 0.5 for P2 multi-start, 1.0 for P2 single-start.
3. **One kick was not enough for P4.** No setting helped (statistical test: p = 0.98 and 0.86, no gain). Some kicks land well and some land badly, which suggests trying several.
4. **The P4 kick 1.0 result from round 2 (0.574) did not replicate** (0.276 here). We report that openly.

### Round 4: six ways to spend the same 15,000 evaluations (plan written down before running)

To fix both problems (kick size depends on the start; one kick is a gamble), round 4 compared six pipelines. Each gets the same total budget, 15,000 evaluations, and each picks its answer by measured energy only. The plan, including which pipeline was the "primary" one, was fixed in `restarts.py` before any result was seen.

| Pipeline | What it does |
|---|---|
| COBYLA | team's optimizer from SPIQ, no kick (may stop early by itself) |
| NFT | NFT from SPIQ, no kick, all 15,000 evaluations |
| kick 0.5 ×5 | 5 restarts × 3,000 evaluations, kick 0.5 |
| kick 1.0 ×5 | 5 restarts × 3,000 evaluations, kick 1.0 |
| mixed kicks ×5 | 5 restarts, one each at kick 0.2, 0.35, 0.5, 0.75, 1.0 |
| **halving (primary)** | 15 short runs of 400 evaluations (3 at each kick size), then the 3 lowest continue for 3,000 each |

"Halving" was chosen as the primary pipeline in advance because it needs no kick size: it is supposed to discover the right one.

**Figure: `K6_round4_six_pipelines.png`.** Mean P(optimal plan) (bars) and every run (dots), 20 runs per pipeline.

| Start | COBYLA | NFT | kick 0.5 ×5 | **kick 1.0 ×5** | mixed ×5 | halving |
|---|---|---|---|---|---|---|
| P2 multi | 0.500 (0) | 0.904 (16) | 0.999 (20) | **0.998 (20)** | 0.999 (20) | 0.998 (20) |
| P2 single | 0.000 (0) | 0.050 (1) | 0.533 (9) | **0.949 (19)** | 0.768 (13) | 0.786 (15) |
| P4 single | 0.250 (0) | 0.398 (1) | 0.337 (3) | **0.761 (13)** | 0.603 (9) | 0.196 (1) |

(numbers in brackets = runs out of 20 ending above 0.9)

**Kick 1.0 ×5 is the only pipeline at or near the top on all three starts.**

### The primary pipeline failed, and why

**Figure: `K7_halving_failure.png`.** Which kick size halving's final answer came from, out of 20 runs.

On P4, halving picked **kick 0.2 in 19 of 20 runs** and ended at 0.196, *worse than COBYLA*. Its 400-evaluation screen kept the runs with the lowest energy after 400 evaluations. Small kicks start close to SPIQ, so they look best early: they slide back down SPIQ's dip quickly. But they stay stuck near SPIQ, while kick-1.0 runs start high and end much lower after 3,000 evaluations. **A short screen rewards the wrong thing.** We report this as a failed pre-registered pipeline. Kick 1.0 ×5 was picked *after* seeing round 4, which is why round 5 was needed.

### Round 5: the claim, on seeds never seen before

Kick 1.0 ×5 vs COBYLA, 20 runs each, with brand-new seeds (reps 101–120).

**Figure: `K8_round5_confirmation.png`.** Top: P(optimal plan) of every run. Bottom: exact energy of every run. Black bar = mean.

| Start | COBYLA P(opt) | **Kick 1.0 ×5 + NFT P(opt)** | Runs > 0.9 (95% interval for the share) | COBYLA energy | Kick + NFT energy | Ground | Test p-value |
|---|---|---|---|---|---|---|---|
| P2 multi-start | 0.500 | **0.977** | 19/20 (0.75–1.00) | 3.742 | 3.660 | 3.653 | 3.4 × 10⁻⁸ |
| P2 single-start | 0.000 | **0.973** | 19/20 (0.75–1.00) | 5.127 | 3.661 | 3.653 | 3.4 × 10⁻⁸ |
| P4 single-start | 0.250 | **0.693** | 11/20 (0.32–0.77) | 7.425 | 5.190 | 4.908 | 6.1 × 10⁻⁴ |
| P1, P3 | 1.000 | 1.000 | 20/20 | at ground | at ground | | no change |

(one-sided Mann–Whitney test on P(optimal plan), kick + NFT vs COBYLA)

- **It replicated:** P2 0.949 → 0.973 and P4 0.761 → 0.693 between round 4 and round 5.
- **P4 is not solved.** 11 of 20 runs succeed, and all 11 have exact energy ≤ 4.943 (ground 4.908). The other 9 end on mixed or wrong plans at 5.095–6.219. When the energy gets within about 0.035 of the ground state, the plan is right.
- **It never hurts P1 and P3,** where SPIQ is already exact, because SPIQ's point is candidate 0 and nothing beats its energy.

---

## 2.6 Twelve new problems: does it carry over?

Same pipeline, nothing re-tuned, on the 12 new problems from Part 1 (3 runs each, new formulation, 8 qubits, 64 angles).

**Figure: `K9_new_problems_every_run.png`.** Every kick + NFT run: P(optimal plan) of the best point so far, over 15,000 evaluations. Dotted lines = each new kick. Diamonds = COBYLA (orange) and NFT (blue) from SPIQ at 1,000 evaluations.

| Pipeline (new formulation, 12 problems) | Mean P(optimal plan) | Problems ≥ 0.99 |
|---|---|---|
| SPIQ start (no optimizer) | 0.365 | 1 of 12 |
| SPIQ + COBYLA (team) | 0.365 | 1 of 12 |
| SPIQ + NFT, 1,000 evaluations | 0.465 | 2 of 12 |
| **SPIQ + kick 1.0 ×5 + NFT** | **0.904** | **9 of 12** |

- **Solved (≥ 0.99, mean of 3 runs):** all 3 chains, all 3 cycles, STAR #41000, CLIQUE #43001 and #43002.
- **CLIQUE #43000 (0.83):** two runs reach 0.999; one stays at 0.50.
- **STAR #41003 (0.52):** runs end at 0.55, 0.996 and 0.00. In the third run, a kick found a point with *lower energy than SPIQ* (14.81 vs 15.74; ground 14.57) but on a wrong plan. Lowest-energy selection then preferred it. This is visible in K9 as a line dropping from 0.5 to 0.
- **STAR #41004 (0.51):** all three runs end at energy 11.778 (ground 11.761): half on the optimal plan, half on one competing plan only **0.017** above it in energy. That competitor is a near-tie the optimizer could not separate.

### Is it just the bigger budget?

Kick + NFT spends 15,001 evaluations, against 1,000 for COBYLA or NFT. Compare at **equal budget**:

| At 1,000 evaluations | Mean P(optimal plan), 12 new problems |
|---|---|
| COBYLA from SPIQ | 0.365 |
| NFT from SPIQ | 0.465 |
| Kick + NFT (first kick, 999 evaluations of NFT so far) | **0.660** |

**At the same 1,000 evaluations, kick + NFT is already well ahead.** Most of the rest comes by the second kick (0.765 at 3,000 evaluations), and five kicks reach 0.904.

The old formulation on the same pipeline (all 12 problems, 36 runs) reaches 0.628 on average (0.383 at SPIQ), vs 0.904 for the new formulation. It reaches ≥ 0.99 on 4 of 12 problems (new: 9 of 12). On the three cliques it ends at 0.20, 0.03 and 0.02. The old formulation cannot go higher on problems where its ground state is the wrong plan (Part 1, Section 1.7).

---

## 2.7 Control experiment: does the kick need to start from SPIQ?

Section 2.4 showed that a kick of σ = 1.0 throws the circuit about as far from SPIQ as a kick of 2.0. So a professor will fairly ask: **is staying near SPIQ doing anything, or would fully random restarts work just as well?** We tested it directly.

**Control pipeline ("random ×5 + NFT").** Everything is identical to kick + NFT except where each restart starts:

- the same SPIQ circuit (one angle per gate);
- SPIQ measured once as candidate 0;
- 5 restarts × 3,000 NFT evaluations;
- the same lowest-measured-energy selection;
- the same 15,001 evaluations and exact scoring;
- new seeds.

The one change: each restart begins from angles drawn **uniformly at random in [−π, π)** instead of SPIQ + 1.0·N(0, I).

**Figure: `K11_random_restart_control.png`.** COBYLA from SPIQ vs kick 1.0 ×5 + NFT vs random ×5 + NFT. Bars = mean, dots = every run.

| Case (runs) | COBYLA from SPIQ | Kick 1.0 ×5 + NFT | **Random ×5 + NFT** |
|---|---|---|---|
| P2 multi-start (20) | 0.500 (0 > 0.9) | 0.977 (19) | **0.999 (20)** |
| P2 single-start (20) | 0.000 (0) | 0.973 (19) | **0.975 (19)** |
| P4 single-start (20) | 0.250 (0) | 0.693 (11) | **0.712 (11)** |
| 12 new problems (36) | 0.365 (3) | 0.904 (30) | **0.963 (33)** |

(numbers in brackets = runs ending above 0.9)

**Result: random restarts do as well as kicks from SPIQ, and slightly better on the new problems.** The differences (0.02–0.06) are within run-to-run noise at these sample sizes. Random restarts had fewer of the "stuck at 0.5" endings: CLIQUE #43000 1.00 / 1.00 / 1.00 vs 1.00 / 1.00 / 0.50, and STAR #41003 0.70 / 1.00 / 1.00 vs 0.55 / 1.00 / 0.00. STAR #41004, the near-tie problem, stayed hard for both. In only 1 of 76 control runs did SPIQ's own point remain the best candidate.

**What this means: we now know which parts of the pipeline do the work.**

| Ingredient | Needed? | Evidence |
|---|---|---|
| **Leaving SPIQ's flat spot** | **Yes** | From SPIQ without leaving it, COBYLA stays put and NFT stays low (round 4 "NFT": P2 single 0.050, P4 0.398) |
| **NFT as the optimizer** | **Yes** | Kicked COBYLA: 0 of 20 runs above 0.9 in every round-3 cell |
| **Several restarts, chosen by measured energy** | **Yes** | One kick: P4 0.28 (round 3); five: 0.69–0.76 |
| **Restarts centred on SPIQ's angles** | **No (on these problems)** | Uniform random restarts match or beat them on P2, P4 and all 12 new problems |
| **SPIQ's circuit (one angle per gate)** | Not tested separately | Every run here uses it; the 4-angle standard circuit with NFT reaches only 0.22 (Part 1) — but that also used a single start, so the circuit alone is not isolated |
| **SPIQ's point as candidate 0** | A safety net | It guarantees the result is never worse than SPIQ by measured energy (P1, P3 stay at 1.00) |

**So the honest description of the method is "multi-start NFT": several NFT runs from fresh starting points, keeping the lowest measured energy.** The kick from SPIQ is one good way to make those fresh starts, and a uniformly random start is an equally good one. SPIQ's value in this pipeline is the circuit layout it provides and the safe candidate 0, not the neighbourhood of its angles. This also matches Section 2.4: at σ = 1.0 the kicked point has already lost SPIQ's energy advantage completely.

---

## 2.8 P5 (24 qubits, 252 angles): in progress

P5 is the scale test. From SPIQ, every one of the 10 optimizers leaves P5 at P(optimal plan) = 0.00 (Part 1, B2).

**Figure: `K10_P5_interim.png`.** Energy and P(optimal plan) of each run's best candidate so far, after each kick. Restart 0 is the SPIQ start.

Four kick + NFT runs with new seeds are running. P5's statevector is 16.8 million amplitudes, and each evaluation takes about **10 seconds** on this laptop with the other jobs running. The restarts took **8.4, 5.9, 5.6 and 4.0 hours**. Results after four of five kicks (each cell: where that restart ended, energy / P(optimal plan); **bold** = the run's kept best after that kick):

| Run | Kick 1 | Kick 2 | Kick 3 | Kick 4 | Kept best after kick 4 |
|---|---|---|---|---|---|
| 1 | **12.596 / 0.003** | 15.131 / 0.002 | 14.536 / 0.021 | 14.426 / 0.000 | 12.596 / 0.003 |
| 2 | **11.086 / 0.380** | 17.551 / 0.160 | 17.597 / 0.120 | **10.316 / 0.039** | 10.316 / **0.039** |
| 3 | 17.765 / 0.000 (SPIQ kept) | **9.768 / 0.976** | 15.831 / 0.024 | 13.167 / 0.293 | 9.768 / **0.976** |
| 4 | **12.000 / 0.612** | 17.445 / 0.144 | 12.683 / 0.415 | 12.477 / 0.057 | 12.000 / **0.612** |

SPIQ start: energy 17.558, P(optimal plan) 0.000. Ground state: 8.340.

**What this shows so far:**

- **Kick + NFT escapes SPIQ on P5.** Every optimizer from SPIQ stays at 0.00; here run 3 reaches **0.976** and run 4 **0.612**. Every run's best energy is well below SPIQ's 17.558.
- **Several restarts matter.** Run 3's first kick found nothing better than SPIQ; its second solved the problem.
- **Selection by measured energy keeps the best point,** so a worse restart never overwrites a better one (kicks 2–4 in runs 1, 3, 4).
- **But lower energy is not always a better plan probability (run 2, kick 4).** Kick 4 found a lower energy (10.316 vs 11.086), so the selection rule kept it, and P(optimal plan) fell from **0.380 to 0.039**. Energy is an average over all plans. A point can lower that average by putting most of its weight on the second- and third-best plans instead of the best one. Neither point is at the ground state (8.340), and only the ground state guarantees the optimal plan. We have not yet checked which plans run 2's new point favours; that needs one more 24-qubit evaluation and will be added after the runs finish.
- **Honest consequence:** the pipeline can only see energy, not P(optimal plan), so on a problem this large an unconverged run can trade one for the other. The fix is more budget per restart (to get closer to the ground state, where the two agree), not a different selection rule, because P(optimal plan) is not observable without knowing the answer.

The last kick per run will take about another 4 hours. Results are interim and may change.

---

## 2.9 What it costs

| | COBYLA from SPIQ | Kick 1.0 ×5 + NFT |
|---|---|---|
| Evaluations per run | ~1,240 (P2), ~2,900 (P4) — stops by itself | 15,001 |
| Shots per run (× 10,240) | 12.7 million (P2) | **154 million** |
| Ratio | | ~12× (P2), ~5× (P4) |
| Wall time, 8 qubits (simulated) | seconds | ~100 s |
| Wall time, P5 24 qubits (simulated) | | ~40 h per run on a laptop |

On real hardware, 154 million shots per run is expensive. But the equal-budget comparison (Section 2.6) shows the gain is not only from spending more: at 1,000 evaluations the kick already beats both alternatives.

---

## 2.10 Conclusions of Part 2

1. **Why optimizers stall at SPIQ (measured):** at all 17 SPIQ points, the gradient is exactly zero and no single angle goes downhill. Slope-following and one-direction methods have nowhere to go.
2. **Why a kick:** to leave SPIQ's dip completely. σ = 1.0 does that; small kicks slide back. A kick alone always makes things worse (0% of kicks improve), so it only works with a strong optimizer after it.
3. **Why NFT:** the energy is an exact sine wave in each angle (checked to 10⁻¹⁴), so NFT's one-angle exact jumps are the natural optimizer. Kicked COBYLA never succeeded.
4. **Why five restarts and lowest-measured-energy selection:** single kicks are a gamble (round 3, P4). Five tries with honest selection make it reliable, and keeping SPIQ as candidate 0 means it can never return something worse than SPIQ by measured energy.
5. **Result:** P2 0.50 → 0.977 and 0.00 → 0.973; P4 0.25 → 0.693; 12 new problems 0.365 → 0.904 (9 of 12 solved); P5 0.00 → 0.976 and 0.612 on two of four runs after four of five kicks (interim).
6. **Control experiment:** fully random restarts + NFT match or beat kicks from SPIQ (P2 0.999 / 0.975, P4 0.712, new problems 0.963). The working ingredients are escaping the flat spot, NFT, and multiple restarts selected by measured energy. The recommendation is therefore **multi-start NFT on SPIQ's circuit, keeping SPIQ's point as candidate 0**. Kick or random starts both work.
7. **Honesty:** the pre-registered "halving" pipeline failed. Kick 1.0 was confirmed only on seeds it had not seen. The control showed that the kick's link to SPIQ is not what drives the gain. P4 remains partly unsolved, and some problems have near-tie competing plans.

## 2.11 Limits

- Noise-free simulation; no hardware noise. 15,001 × 10,240 shots per run is costly on hardware.
- σ = 1.0 and 5 × 3,000 were chosen on P2 and P4 and then confirmed on unseen seeds and 12 new problems. It is not proven to be the best setting for every problem size.
- P5 is unfinished.
- Selection by lowest measured energy can prefer a wrong plan that has lower energy than SPIQ (STAR #41003). In the new formulation this happens only when the optimizer stops short of the ground state.

---

## Likely questions and answers

**Q1. Why not just use a better starting point than SPIQ?**
SPIQ's point is very good in energy: it starts lower than anything the standard start reaches after 1,000 evaluations (Part 1, A1). The problem is not where it is, but that it is a flat spot. The pipeline keeps SPIQ's point as candidate 0, so it never does worse than SPIQ, and escapes the flat spot with restarts.

**Q2. Isn't a kick of 1 radian basically a random restart?**
Yes, and we tested it. Five uniformly random restarts + NFT gave 0.999 / 0.975 (P2), 0.712 (P4) and 0.963 (12 new problems), against 0.977 / 0.973 / 0.693 / 0.904 for kicks (K11). The kick's link to SPIQ's angles is not what makes it work. NFT plus several fresh starts is.

**Q2b. Then what is SPIQ still good for?**
Three things. It provides the circuit with one angle per gate that every successful run here uses. Its point is the safety-net candidate 0 (P1 and P3 stay at 1.00). And on its own it already beats 1,000 evaluations from the standard start (Part 1). What it does not do is give the optimizer a good place to *continue from*, because its point is a flat spot.

**Q3. How do you know the pipeline isn't "cheating" by looking at the answer?**
It picks only by measured 10,240-shot energy. P(optimal plan) and the ground state are computed only afterwards, for the report. The selection could be done on real hardware exactly as written.

**Q4. Why is the gradient exactly zero at SPIQ?**
We measured it to be zero (10⁻¹⁵–10⁻¹⁴) at all 17 points. SPIQ chooses the best Clifford setting it finds, and turning any one angle by ±π/2 lands on another Clifford setting that is no better. That explains why no single angle goes downhill. The exact zero slope is a measured fact we do not claim to have proven in general.

**Q5. Why does NFT work here when COBYLA doesn't?**
Away from SPIQ's point the landscape has normal slopes (1–12). NFT solves each angle exactly with big ±π/2 probes that shot noise does not swamp, and it never steps uphill along the angle it is fixing. COBYLA, given a kick, still stops early on a local fit (0 of 20 runs above 0.9).

**Q6. Why 5 restarts of 3,000, not 15 of 1,000 or 1 of 15,000?**
One long NFT run (round 4, "NFT" column) reached only 0.050 on P2 single-start and 0.398 on P4, so it stays stuck. Many short runs are the halving idea, and short screens picked the wrong kicks (K7). 5 × 3,000 gives each restart about 23 sweeps on 64 angles: long enough to converge, with enough restarts to survive bad kicks.

**Q7. P4 is only 0.693. Why not higher?**
9 of 20 runs end at energies 5.095–6.219, on mixed or wrong plans. All successes are within 0.035 of the ground state. P4 has 140 angles, so 3,000 evaluations is only ~10.5 sweeps. More restarts or longer segments might help, but that has not been tested yet.

**Q8. What about the failed halving result: doesn't it weaken the claim?**
It strengthens the method of the claim. Halving was the pre-registered primary pipeline and it failed on P4 for a reason we can explain (short screens favour small kicks). The pipeline we now recommend was then tested on 20 new seeds per case and 12 new problems it had never seen.

**Q9. How many evaluations does this need compared with COBYLA?**
15,001 vs ~1,240 (P2) and ~2,900 (P4): about 12× and 5×. At equal budget (1,000 evaluations) it is still ahead on the new problems: 0.660 vs 0.365.

**Q10. Does it work in the old formulation?**
It improves it (0.383 → 0.628 over all 12 new problems, vs 0.904 for the new formulation), but it cannot fix problems where the old formulation's ground state is the wrong plan. Lowering the energy there moves *away* from the optimal plan (Part 1, Section 1.7).

**Q11. On P5 run 2, the energy went down but P(optimal plan) went from 0.38 to 0.04. Doesn't that break your method?**
It shows its limit, not a bug. The pipeline can only measure energy; P(optimal plan) needs the answer, so it cannot be used to choose. Energy is an average over plans, so above the ground state a lower average can come from favouring the second-best plans over the best one. Both of run 2's points are far from the ground state (10.3 and 11.1 vs 8.34). Where a run gets close to the ground state, as run 3 did (9.77), energy and P(optimal plan) agree (0.976). On P2, P4 and the 12 new problems this trade-off was rare because runs got much closer to the ground state. The remedy is a larger budget per restart on big problems.

---

## Test yourself

1. Write the 5 steps of one kick + NFT run.
2. Why does NFT need only 2 evaluations per angle step?
3. At a SPIQ point, how many single angles go downhill? What is the gradient size?
4. A kick of σ = 1.0 is applied. Does the energy go up or down before NFT runs?
5. Why was kick 1.0 tested again in round 5?
6. Why does the pipeline keep the SPIQ point as candidate 0?

**Answers.**
1. Measure SPIQ (candidate 0) → 5 × [kick x₀ + 1.0·N(0, I), NFT 3,000 evaluations] → keep lowest measured energy → score exactly.
2. The energy along one angle is a sine wave with 3 unknowns. The current value is known, so 2 more evaluations (at ±π/2) determine the whole curve and its minimum.
3. Zero downhill angles at all 17 points; gradient ≈ 10⁻¹⁵–10⁻¹⁴, i.e. zero.
4. Up, every time (0 of 200 kicks went below SPIQ). NFT does the improving afterwards.
5. Because it was chosen after seeing round 4's results. A fair claim needs seeds that played no part in the choice.
6. So the pipeline can never return something worse than SPIQ by measured energy: on P1 and P3, where SPIQ is exact, it keeps 1.00.

---

## Slide plan for Part 2 (about 11 slides)

| # | Title | Figure | Say this in one sentence |
|---|---|---|---|
| 1 | The problem: stuck at SPIQ | Part 1 `B2` (SPIQ row) | From SPIQ, COBYLA ends exactly where it started on every problem. |
| 2 | Why: a flat spot | `K2_spiq_flat_spot.png` | At all 17 SPIQ points the gradient is zero and no single angle goes downhill. |
| 3 | The landscape around SPIQ | `K3_landscape_cuts.png` | SPIQ sits in a narrow dip; by distance 1 the energy is on a high plateau. |
| 4 | NFT in one slide | sine-wave formula + 4 steps | Energy is an exact sine wave per angle, so NFT jumps to each angle's exact best. |
| 5 | The pipeline | `K1_pipeline_diagram.png` | Measure SPIQ, 5 × (kick 1.0 + NFT 3,000), keep lowest measured energy. |
| 6 | A kick alone makes it worse | `K4_kick_alone.png` | The kick is a jump, not an improvement; NFT does the work. |
| 7 | Round 3: kick needs NFT | `K5_round3_kick_needs_NFT.png` | Kicked COBYLA never succeeds; kicked NFT does; one kick is a gamble on P4. |
| 8 | Round 4: six pipelines | `K6…` + `K7_halving_failure.png` | Kick 1.0 ×5 wins on all three; the pre-registered halving failed on P4, for a clear reason. |
| 9 | Round 5: confirmed on unseen seeds | `K8_round5_confirmation.png` | P2 0.977 / 0.973, P4 0.693, vs COBYLA 0.50 / 0.00 / 0.25. |
| 10 | 12 new problems | `K9_new_problems_every_run.png` + equal-budget table | 0.365 → 0.904; ahead even at equal budget (0.660 vs 0.365). |
| 11 | Control: what really drives it | `K11_random_restart_control.png` + ingredients table | Random restarts do as well as kicks: the gain is NFT + leaving the flat spot + several restarts. |
| 12 | P5 and cost | `K10_P5_interim.png` + cost table | P5: 0.00 → 0.976 / 0.612 on 2 of 4 runs after 4 of 5 kicks (interim); run 2 shows lower energy can lower P(opt); costs ~5–12× COBYLA's evaluations. |
