# Kick + NFT: the complete guide

**What this file is.** Everything about the kick + NFT work, in one place, written so you can present it and answer any question about it:

- why all 10 optimizers fail from the team's SPIQ point;
- why we came up with a "kick", and why we put it on NFT and not another optimizer;
- what seeds are and how they were used;
- every problem and every round of runs;
- every figure: what is plotted, how to read it, what it proves;
- the 12 new problems, the control experiment and P5;
- conclusions, what they imply, the limits, likely questions with answers.

Every number below comes from a result file in this project. The file is named next to each claim, so you can always point to the evidence.

Formulation throughout: the **new** ("qaoa", Divya's) formulation, unless a section says "old formulation".

---

## Contents

0. [What is in this folder (file map)](#0-what-is-in-this-folder-file-map)
1. [The whole story in 12 lines](#1-the-whole-story-in-12-lines)
2. [Concepts you must know (glossary)](#2-concepts-you-must-know-glossary)
3. [Seeds: what they are and how they were used](#3-seeds-what-they-are-and-how-they-were-used)
4. [The problems we tested on](#4-the-problems-we-tested-on)
5. [Step 1. All 10 optimizers fail from SPIQ](#5-step-1-all-10-optimizers-fail-from-spiq)
6. [Step 2. Why they fail: SPIQ's point is a flat spot](#6-step-2-why-they-fail-spiqs-point-is-a-flat-spot)
7. [Step 3. Where the "kick" idea came from, and what a kick does](#7-step-3-where-the-kick-idea-came-from-and-what-a-kick-does)
8. [Step 4. Why NFT, out of all 10 optimizers](#8-step-4-why-nft-out-of-all-10-optimizers)
   - [8b. The intuition: what the kick means, and why NFT is the one that works](#8b-the-intuition-what-the-kick-means-and-why-nft-is-the-one-that-works)
9. [Step 5. Building the pipeline: rounds 2–5](#9-step-5-building-the-pipeline-rounds-25)
10. [Step 6. The 12 new problems](#10-step-6-the-12-new-problems)
11. [Step 7. Control experiment: is it the kick, or just restarting?](#11-step-7-control-experiment-is-it-the-kick-or-just-restarting)
12. [Step 8. P5, the 24-qubit scale test (interim)](#12-step-8-p5-the-24-qubit-scale-test-interim)
13. [What it costs](#13-what-it-costs)
14. [Final conclusions and what they imply](#14-final-conclusions-and-what-they-imply)
15. [Limits (say these before you are asked)](#15-limits-say-these-before-you-are-asked)
16. [Figure-by-figure guide](#16-figure-by-figure-guide)
17. [Likely questions and answers](#17-likely-questions-and-answers)
18. [Slide plan](#18-slide-plan)
19. [Test yourself](#19-test-yourself)

---

## 0. What is in this folder (file map)

### This presentation folder: `new_formulation_pilot/PRESENTATION_kick_NFT/`

```
PRESENTATION_kick_NFT/
├── KICK_NFT_COMPLETE_GUIDE.md        this file
├── figures/                          17 figures, numbered in the order you present them
│   ├── 01_all10_optimizers_stuck_at_SPIQ.png       the problem: 10 optimizers × P1–P5 from SPIQ
│   ├── 02_spiq_flat_spot.png                       the cause: zero gradient, no downhill angle
│   ├── 03_landscape_around_SPIQ.png                energy and P(opt) walking away from SPIQ
│   ├── 04_kick_alone_always_worse.png              a kick by itself always raises the energy
│   ├── 05_pipeline_diagram.png                     the final pipeline as a flow chart
│   ├── 06_round1_which_optimizer_to_kick.png       why NFT: 5 optimizers kicked, only NFT improves
│   ├── 07_round2_kick_size_sweep.png               choosing the kick size (exploratory)
│   ├── 08_round3_kick_needs_NFT.png                round 3: kicked COBYLA vs kicked NFT, new seeds
│   ├── 09_round4_six_pipelines.png                 round 4: six pipelines, same budget
│   ├── 10_round4_halving_failure.png               why the pre-registered pipeline failed
│   ├── 11_round5_confirmation_unseen_seeds.png     round 5: the main claim
│   ├── 12_new_problems_kick_heatmap.png            12 new problems × every pipeline (heatmap)
│   ├── 13_new_problems_every_run.png               12 new problems, every kick + NFT run over time
│   ├── 14_new_problems_trend_vs_COBYLA_NFT.png     12 new problems, averaged trends
│   ├── 15_random_restart_control.png               control: random restarts vs kicks
│   ├── 16_P5_interim.png                           P5 (24 qubits), progress per kick (still running)
│   └── 17_old_formulation_ground_state_wrong_plan.png  why kick + NFT cannot fix the old formulation
├── data/
│   ├── runs_new.jsonl                one line per run on the 12 new problems (both formulations; COBYLA, NFT, kick)
│   ├── new_problems.jsonl            the 12 problems: tables, sizes, predicates, selectivities, optimal cost
│   ├── runs_ablation_random.jsonl    the control experiment (random restarts + NFT), 96 runs
│   ├── landscape.json                gradient, curvature, kick-alone and line-cut measurements at 17 SPIQ points
│   ├── part2_numbers.json            every number shown in figures 02–05, 08–11, 13, 15, 16
│   └── new_figures_numbers.json      every number shown in figures 01, 06, 07, 12
└── scripts/
    ├── make_new_figures.py           makes figures 01, 06, 07, 12 (run it from this folder; it reads the original files)
    ├── landscape.py                  measures the landscape at SPIQ points (copy, for reference)
    ├── part2_figures.py              makes figures 02–05, 08–11, 13, 15, 16 (copy, for reference)
    └── ablation_random.py            runs the control experiment (copy, for reference)
```

The three copied scripts use paths relative to `optimizer_study/MEETING_2026-11-18/scripts/`. To re-run them, run the originals there. Figures 14 and 17 come from `optimizer_study/comparison_2x2/analyze_2x2.py`.

### The pilot folder itself: `new_formulation_pilot/` (all the original runs, 5–6 Nov, plus P5)

| File / folder | What it contains | What the information means |
|---|---|---|
| `REPORT.md` | the original write-up of rounds 1–5 (6 Nov) | the main claim table is at the top |
| `pilot.py` | the engine for rounds 1–3: builds each problem's energy table, loads the team's SPIQ points, runs one optimizer with an optional kick | every round 1–3 result line was made by this file's `run()` |
| `summ.py` | a one-screen summary table of any results file | for quick checks only |
| `pilot_round1/runs_P1_P2_P3_b1000.jsonl` | **Round 1, no kick**: 10 optimizers × P1–P3 × 3 starts (random, SPIQ single, SPIQ multi) × 5 trials = 450 runs, 1,000 evaluations | shows every optimizer stays at SPIQ |
| `pilot_round1/runs_P1_P2_P3_b1000_kicks.jsonl` | **Round 1, kicked**: 5 optimizers (NFT, COBYLA, POWELL, ADAM, BOBYQA) × kick 0.05 / 0.2 / 0.5 × P1–P3 × 2 SPIQ starts × 5 trials = 450 runs | shows only NFT improves after a kick (figure 06) |
| `pilot_round1/log_p123.txt`, `log_p123k.txt` | the console output of those runs | one line per run |
| `pilot_round2/runs_P2_P4_b1000.jsonl` | **Round 2, no kick**: NFT and COBYLA, P2 and P4, 20 trials | baseline for round 2 |
| `pilot_round2/runs_P2_b1000_kicks.jsonl`, `runs_P2_b3000_kicks.jsonl` | **Round 2, P2 kicks** σ = 0.2, 0.35, 0.5, 0.75, 1.0; NFT and COBYLA at 1,000; NFT at 3,000; 20 trials | kick-size sweep (figure 07) |
| `pilot_round2/runs_P4_b1000_kicks.jsonl`, `runs_P4_b3000_kicks.jsonl` | **Round 2, P4 kicks** σ = 0.2, 0.5, 1.0; 10 trials | kick-size sweep on P4 |
| `pilot_round2/launch.sh`, `log_*.txt`, `round2.done` | the exact commands, console output, a "finished" marker | |
| `pilot_round3/runs_P2_P4_b3000_t100.jsonl` | **Round 3, no kick**, NFT and COBYLA, 3,000 evaluations, 20 trials, **new seeds (trials 101–120)** | |
| `pilot_round3/runs_P2_P4_b3000_kicks_t100.jsonl` | **Round 3, kicked** 0.5 and 1.0, NFT and COBYLA, same budget and seeds | figure 08 |
| `pilot_round3/launch_round3.sh`, logs | commands and output | |
| `analyze_round3.py` | round 3 tables, Mann–Whitney tests | |
| `restarts.py` | the engine for rounds 4 and 5: six multi-restart pipelines; with `CONFIRM=1` it runs round 5 | the plan (incl. "halving is primary") is written in its first lines, before any run |
| `runs_round4.jsonl` | **Round 4**: 6 pipelines × 5 problem-starts × 20 runs = 600 runs, 15,000 evaluations each | figures 09, 10 |
| `analyze_round4.py`, `launch_round4.sh`, `log_round4.txt` | tables, command, output | |
| `runs_round5_confirm.jsonl` | **Round 5**: COBYLA vs kick 1.0 ×5 + NFT, 5 problem-starts × 20 runs, **brand-new seeds (reps 101–120)** = 200 runs | figure 11; **the main claim** |
| `analyze_round5.py`, `launch_round5.sh`, `log_round5.txt` | tables, command, output | |
| `p5_kick.py`, `launch_p5_kick.sh` | P5 kick 1.0 ×5 + NFT, 4 runs, new seeds | |
| `log_p5_kick.txt` | P5 progress, one line per finished kick per run | figure 16 is drawn from this until the runs finish |
| `runs_p5_kick.jsonl` | P5 final results (empty until the runs finish) | |

### Related folders elsewhere in `optimizer_study/`

| Folder / file | What it is |
|---|---|
| `new_code_study/runs_b1000.jsonl` | Study 2: 10 optimizers on P1–P5 from every start (source of figure 01) |
| `new_problems/` | the 12 new problems: `pick_problems.py` (how they were chosen), `spiq_new.py` (SPIQ run on them), `run_new.py` (COBYLA, NFT, kick runs), `runs_new.jsonl` |
| `many_problems/SUMMARY.md` | the earlier 1,088-problem study; its "Future work" item 4 is where the kick idea came from |
| `MEETING_2026-11-18/` | the meeting folder: Part 1 (formulation × start) and Part 2 (kick + NFT) |
| `comparison_2x2/optimizer_2x2_report.html` | the full 2×2 report (HTML) |

---

## 1. The whole story in 12 lines

1. The team's SPIQ method gives a starting point for QAOA. From that point, **all 10 optimizers we tried stay exactly where they started** on every problem (figure 01).
2. **Why:** at all 17 SPIQ points we measured, the slope (gradient) is exactly zero and **not one angle goes downhill** when turned on its own (figure 02). It is a flat spot.
3. The idea of a **kick** (random noise on every angle, to leave the flat spot) came from the future-work list of our earlier 1,088-problem study.
4. A kick on its own always makes the energy worse (0 of 13,600 kicks improve: 200 per size × 4 sizes × 17 SPIQ points; figure 04). So the optimizer after the kick has to do all the work.
5. **Round 1** tried kicks on 5 different optimizers. **Only NFT** ended below SPIQ's energy after a kick, and only NFT produced any run above 0.9 (figure 06). That is why NFT was chosen.
6. NFT fits these circuits exactly: the energy is a perfect sine wave in each angle (checked to 10⁻¹⁴), and NFT jumps to the exact bottom of that sine wave, one angle at a time.
7. **Rounds 2–3** chose a kick size and showed one kick is a gamble (figures 07, 08). **Round 4** compared six pipelines; the one picked in advance ("halving") **failed** (figures 09, 10). Kick 1.0 ×5 + NFT did best.
8. **Round 5**, on seeds never seen before: P2 **0.50 → 0.977** and **0.00 → 0.973**, P4 **0.25 → 0.693** chance of the optimal plan (figure 11).
9. **12 new problems**, nothing re-tuned: **0.365 → 0.904** on average, 9 of 12 solved (figures 12–14). It is ahead even at equal budget (0.660 vs 0.365 at 1,000 evaluations).
10. **Control:** fully random restarts + NFT do just as well (0.999 / 0.975 / 0.712 / 0.963; figure 15). What matters is **leaving the flat spot + NFT + several restarts chosen by measured energy**: "multi-start NFT".
11. **P5** (24 qubits), interim after 4 of 5 kicks: two of four runs reach **0.976** and **0.612**, where every optimizer from SPIQ stays at 0.00 (figure 16). One run shows the limit: lower energy, lower chance of the optimal plan.
12. In the **old formulation**, kick + NFT helps less (0.628) and cannot fix problems where that formulation's lowest-energy state is the wrong plan (figure 17).

---

## 2. Concepts you must know (glossary)

**QAOA (Quantum Approximate Optimization Algorithm).** A quantum circuit with adjustable angles. You run it, measure bitstrings, and each bitstring is a candidate answer (here: a join order). A classical optimizer changes the angles to make good answers more likely.

**Join-order problem.** A database query joins several tables; the order of joins changes the cost hugely. The goal is the cheapest order, the *optimal plan*. P1–P5 have 3, 4, 4, 5 and 6 tables (relations).

**Formulation (QUBO).** The rule that turns the join-order problem into an energy over bitstrings, so that low energy = good plan.
- **Old formulation ("q2o")**: the team's earlier QUBO, with table bits and predicate bits. P1–P3 need 5 / 12 / 20 qubits; P4 and P5 would need 30 and 52 (too many to simulate).
- **New formulation ("qaoa", Divya's)**: R·(R−2) bits for R tables. P1–P5 need **3 / 8 / 8 / 15 / 24 qubits**. All kick + NFT work in rounds 1–5 uses this one.

**Energy.** The QUBO's value, averaged over the bitstrings the circuit produces. Lower is better.
- **Exact energy** = Σ (probability of each bitstring × its energy), from the full simulated state, no noise. Used only for *reporting*.
- **Measured energy** = the average energy of **10,240 simulated shots** (samples): what a real quantum computer would report, with shot noise. The optimizers and the pipeline only ever *see* measured energy.

**Ground state / ground energy.** The bitstring with the lowest energy, and that energy. In the new formulation, the ground state *is* the optimal plan on all 16 problems we checked. In the old formulation it is the wrong plan on 5 of 15 (P2, STAR #41003, cliques #43000/#43001/#43002).

**P(optimal plan).** The probability that a measured bitstring, decoded into a join order, gives the cheapest plan. Decoding uses the team's fallback decoder (`get_raw_join_order` + `postprocess_join_order`) and Divya's fixed cost checker. **This is the number that matters for the database**, and the one we report as the main result. 1.00 = the circuit always outputs the optimal plan.

**Gap closed.** (E_start − E_end) / (E_start − E_ground). 1 = reached the ground state; 0 = did not move; negative = got worse.

**Evaluation.** One run of the circuit with 10,240 shots, returning one measured energy. **Budget** = how many evaluations an optimizer is allowed.

**Angles (parameters).** The numbers the optimizer changes.
- **Standard circuit**: 4 angles (2 layers × γ, β), the textbook QAOA.
- **SPIQ circuit**: **one angle per gate**: 18 (P1), 64 (P2, P3, every new problem), 140 (P4), 252 (P5).

**SPIQ.** The team's method for choosing a starting point. It searches only **Clifford** angle settings (multiples of π/2), which can be simulated cheaply, and keeps the one with the lowest energy. Then a normal optimizer is supposed to continue from it.
- **Single-start SPIQ point**: one SPIQ search.
- **Multi-start ("diverse") SPIQ point**: the best of several SPIQ searches. Available for P1–P3 only.
- Because the angles are Clifford, a SPIQ point's P(optimal plan) is always a fraction m/2ᵏ (0, 0.25, 0.5, 1.00…), which is why so many cells show exactly 0.50 or 0.25.

**Optimizer.** The classical algorithm that changes the angles. The 10 we tested:

| Optimizer | Family | How it moves |
|---|---|---|
| COBYLA | trust region, linear model | the **team's default**; fits a linear model to nearby points, steps, shrinks its region |
| BOBYQA | trust region, quadratic model | like COBYLA with a curved model |
| NFT | sequential, one angle at a time | solves each angle exactly using the sine-wave property (Section 8) |
| POWELL | direction set | line searches along a set of directions |
| NELDER_MEAD | simplex | moves a shape of n+1 points |
| SLSQP | gradient, quasi-Newton | estimates slopes by finite differences |
| L_BFGS_B | gradient, quasi-Newton | same family |
| ADAM | gradient, stochastic | gradient steps with momentum |
| AQGD | gradient, quantum | analytic gradient descent |
| SPSA | stochastic perturbation | moves all angles at once by a random step |

**Kick (perturbation).** x_kicked = x₀ + σ·z, where x₀ = SPIQ's angles, each zᵢ is an independent standard normal random number, and **σ** (sigma, in radians) is the kick size. With σ = 1.0, about 68% of angles move by less than 1 radian (57°), about 95% by less than 2.

**Restart.** Running the optimizer again from a different starting point. **Candidate** = the best point one restart found. **Selection rule** = return the candidate with the **lowest measured energy**. It never looks at P(optimal plan) or the ground state, so it could be used on real hardware.

**Candidate 0.** SPIQ's own point, measured once and kept as a candidate. It guarantees the pipeline never returns something with a worse measured energy than SPIQ.

**Sweep.** NFT updating every angle once. 3,000 evaluations ≈ 1,480 single-angle steps ≈ 23 sweeps over 64 angles, 10.5 over 140, 6 over 252.

**Gradient (slope).** How fast the energy changes when each angle is nudged. Measured exactly with the **parameter-shift rule**: slope_i = [E(θᵢ + π/2) − E(θᵢ − π/2)] / 2.

**Curvature (single-angle).** [E(+π/2) + E(−π/2) − 2·E(θ)] / 2. Positive = uphill both ways along that angle, 0 = flat, negative = downhill one way.

**Flat spot / saddle point.** A point with zero slope that is not the lowest point. One-angle-at-a-time and slope-following methods cannot see a way out.

**Trial / run / rep.** One independent run of a method on a problem with its own seed. "20 runs" = 20 different seeds.

**Pre-registration.** Writing down the plan, including which pipeline is the main one, **before** running. Round 4's plan is the first lines of `restarts.py`.

**Held-out (unseen) seeds.** Seeds that played no part in choosing the method. Results on them are a fair test.

**Mann–Whitney test.** A test that asks "are the values from method A typically higher than those from method B?" without assuming a bell curve. Our p-values are one-sided (is kick + NFT higher than COBYLA?). Small p (e.g. 3.4 × 10⁻⁸) = very unlikely to be luck.

**95% interval for "19/20 runs above 0.9".** The range of the true success rate consistent with 19 of 20 (0.75–1.00).

**Control experiment.** The same pipeline with exactly one thing changed, to test whether that one thing matters.

**Statevector simulation.** We simulate the circuit exactly on a laptop (2ⁿ amplitudes: 256 for 8 qubits, 16.8 million for 24). No hardware noise; shot noise is simulated by sampling 10,240 shots.

---

## 3. Seeds: what they are and how they were used

### What a seed is

Every run uses random numbers: the shot sampling (which bitstrings come out of 10,240 shots), the kick (the random z), and some optimizers' internal randomness. A computer makes these with a **pseudo-random generator**, which produces the same sequence every time it is started with the same **seed** (a starting number).

- **Same seed → exactly the same run**, bit for bit. That makes every result reproducible: anyone can re-run it and get the identical number.
- **Different seed → a different, independent run** of the same method. "20 runs" means 20 different seeds, so we see the spread of outcomes, not one lucky result.

### How our seeds are made

Each run's seed is computed from a **text label** that describes the run, using the CRC-32 checksum (`zlib.crc32`):

```
seed = crc32("round4|P2|spiq_single|fixed1.0|7")
```

So a seed is never typed in by hand or chosen to look good: it follows from what the run *is*. Changing any part of the label (round, problem, start, pipeline, run number) gives a completely different seed.

Inside a run, that one seed starts a generator (`numpy.random.default_rng(seed)`), and every random draw comes from it: the kick vectors, the shot samples, and a fresh sub-seed for each NFT restart.

### The labels used in each experiment

| Experiment | Seed label | Run numbers | New or reused? |
|---|---|---|---|
| Round 1 | `pilot|problem|start(+kick)|optimizer|trial` | 1–5 | new |
| Round 2 | same format | 1–20 | **trials 1–5 share round 1's seeds where the label matches** (e.g. P2 NFT kick 0.2), noted in `pilot_round2/launch.sh` |
| Round 3 | same format | **101–120** | new (trial numbers never used before) |
| Round 4 | `round4|problem|start|pipeline|rep` | 1–20 | new |
| Round 5 | `round4|problem|start|pipeline|rep` | **101–120** | new: same format, but run numbers 101–120 were never used, so every seed is new |
| 12 new problems | `newprob|formulation|problem|start|pipeline|trial` | 1–3 | new |
| Control experiment | `ablation|problem|start|random_uniform|rep` | 1–20 (pilot), 1–3 (new problems) | new |
| P5 | `p5kick|P5|spiq_single|fixed1.0|rep` | 1–4 | new |

### Why seeds matter for the argument

- **Rounds 1, 2 and 4 were used to choose settings.** Their seeds "saw" the decisions, so their numbers are not used as the final claim.
- **Round 3 (trials 101–120) and round 5 (reps 101–120) use seeds that played no part in choosing anything.** That is why round 5 is the main claim: kick 1.0 ×5 was chosen *after* seeing round 4, so it had to be re-tested on new seeds.
- The **12 new problems** and the **control** are new problems *and* new seeds, so they test generalisation.

**One-line answer if asked "what are seeds?":** "A seed is the starting number of the random generator. The same seed reproduces a run exactly; different seeds give independent repeats. Our seeds are computed from each run's label, and the final claim uses only seeds that were never used to choose the method."

---

## 4. The problems we tested on

### The team's five problems (new formulation)

| Problem | Tables (relations) | Predicates | Qubits | SPIQ angles | SPIQ start energy | Ground energy | SPIQ start P(optimal plan) |
|---|---|---|---|---|---|---|---|
| P1 | 3 | 2 | 3 | 18 | 1.176 | 1.176 | **1.00** (already optimal) |
| P2 | 4 | 2 | 8 | 64 | single 5.125 · multi 3.741 | 3.653 | single **0.00** · multi **0.50** |
| P3 | 4 | 6 | 8 | 64 | 2.334 | 2.334 | **1.00** (already optimal) |
| P4 | 5 | 5 | 15 | 140 | 7.425 | 4.908 | **0.25** |
| P5 | 6 | 7 | 24 | 252 | 17.558 | 8.340 | **0.00** |

So the **interesting problems are P2, P4 and P5**: SPIQ leaves room to improve. On P1 and P3 SPIQ is already exact, and the only requirement is "don't break it".

### The 12 new problems

Chosen by `new_problems/pick_problems.py`, deterministically, without looking at any optimizer result:

- **3 chains, 3 stars, 3 cycles**: the first three 4-table queries of each shape (by id) from the team's generated query set that the team had not used, **and** whose optimal plan the team's decoder can actually output (otherwise P(optimal plan) is 0 for every possible circuit).
- **3 cliques**: the generated set had no 4-table cliques, so three were made with the team's own Steinbrunn query generator (seeds 43000+), by the same rule.

| Shape | Problems | Old formulation qubits | New formulation qubits / angles |
|---|---|---|---|
| Chain | #40000, #40001, #40003 | 14 | 8 / 64 |
| Star | #41000, #41003, #41004 | 14 | 8 / 64 |
| Cycle | #42000, #42001, #42002 | 16 | 8 / 64 |
| Clique | #43000, #43001, #43002 | 20 | 8 / 64 |

SPIQ was run on each of them with the team's own SPIQ code (400 generations, 4 starts), for both formulations.

### Which experiment used which problems

| Experiment | Problems | Runs per cell |
|---|---|---|
| Study 2 (figure 01) | P1–P5, 10 optimizers | 5 (P5: 3) |
| Round 1 | P1–P3, 10 optimizers; kicks on 5 | 5 |
| Round 2 | P2, P4 | 20 (P4 kicks: 10) |
| Round 3 | P2 (both starts), P4 | 20 |
| Round 4 | P1, P3 (multi-start), P2 (both), P4 | 20 |
| Round 5 | P1, P3 (multi-start), P2 (both), P4 | 20 |
| 12 new problems | all 12, both formulations | 3 |
| Control | P2 (both), P4, all 12 new | 20 / 3 |
| P5 | P5 | 4 (running) |

---

## 5. Step 1. All 10 optimizers fail from SPIQ

**Figure: `01_all10_optimizers_stuck_at_SPIQ.png`** (data: `new_code_study/runs_b1000.jsonl`, Study 2).

From the team's single-start SPIQ point, each of the 10 optimizers had 1,000 evaluations on each of P1–P5.

| | P1 | P2 | P3 | P4 | P5 |
|---|---|---|---|---|---|
| SPIQ start P(optimal plan) | 1.00 | 0.00 | 1.00 | 0.25 | 0.00 |
| After COBYLA, NFT, BOBYQA, POWELL, NELDER_MEAD, SLSQP, L_BFGS_B, ADAM, AQGD | 1.00 | 0.00 | 1.00 | 0.25 | 0.00 |
| After SPSA | 0.97 | 0.03 | 0.64 | 0.19 | 0.00 |
| Energy change, 9 optimizers | 0.00 | 0.00 to −0.18 (NFT) | 0.00 | 0.00 to +0.01 | −0.01 to +0.02 |
| Energy change, SPSA | +0.07 | +0.76 | +1.34 | +3.55 | +7.75 |

**What this shows.**
- Nine optimizers end **exactly where they started**: same P(optimal plan) to two decimals, energy change of at most 0.18.
- The only real move is NFT on P2 (energy −0.18), and even that does not change P(optimal plan).
- SPSA is the only one that moves a lot, and **it moves uphill** (worse on every problem).
- So the choice of optimizer does not matter from SPIQ: **they are all stuck.** This is the problem kick + NFT solves.

---

## 6. Step 2. Why they fail: SPIQ's point is a flat spot

### The energy is a sine wave in each angle

In SPIQ's circuit every angle drives exactly one rotation gate. If you hold all other angles fixed and turn one angle θ, the energy is

  E(θ) = a + b·cos θ + c·sin θ

an exact sine wave. We checked this on all 17 SPIQ points (P1–P4 and the 12 new problems): fitting that formula to 5 energies along an angle leaves an error of 3 × 10⁻¹⁵ to 2.5 × 10⁻¹⁴, which is rounding error (`data/landscape.json`, `sinusoid_max_residual`).

### What we measured at SPIQ points

**Figure: `02_spiq_flat_spot.png`.**

| SPIQ point | Angles | Slope size at SPIQ | Slope size nearby | Uphill / flat / **downhill** angles |
|---|---|---|---|---|
| P1 | 18 | 1.4 × 10⁻¹⁵ | 1.26 | 15 / 3 / **0** |
| P2 single-start | 64 | 4.6 × 10⁻¹⁵ | 2.98 | 51 / 13 / **0** |
| P2 multi-start | 64 | 4.5 × 10⁻¹⁵ | 2.49 | 56 / 8 / **0** |
| P3 | 64 | 6.0 × 10⁻¹⁵ | 5.03 | 59 / 5 / **0** |
| P4 | 140 | 1.1 × 10⁻¹⁴ | 5.06 | 119 / 21 / **0** |
| 12 new problems | 64 each | 8.7 × 10⁻¹⁵ – 1.8 × 10⁻¹⁴ | 6.65 – 12.33 | 48–60 / 4–16 / **0** |

**At all 17 SPIQ points the slope is zero (to rounding error), and not one angle goes downhill on its own.** A little away from SPIQ, the slope is a normal 1–12.

### What that does to each optimizer

- **COBYLA / BOBYQA** fit a model from nearby points. The model's slope is zero plus shot noise, so there is no direction to follow. They shrink their region and stop at the start.
- **Gradient methods (SLSQP, L_BFGS_B, ADAM, AQGD)** see a zero gradient plus noise. Same result.
- **NFT** looks at one angle at a time. For an uphill angle the bottom of its sine wave *is* the current value, so it stays. For a flat angle, shot noise decides, so it drifts slowly. That is the −0.18 on P2.
- **SPSA** moves all angles at once in a random direction. From here nearly every direction goes uphill, so it gets worse.

### Why SPIQ points look like this

SPIQ keeps the best **Clifford** setting it finds. Turning any one angle by ±π/2 lands on another Clifford setting, which SPIQ could have picked, so it is no better. That explains "no single angle goes downhill". The exact zero slope is something we **measured**, not something we assumed or proved in general.

Being uphill along every *single* angle does not mean uphill along every *combination* of angles. A better region can exist that needs many angles to move together, which none of these optimizers can find from here.

### The landscape further away

**Figure: `03_landscape_around_SPIQ.png`.** Walking from SPIQ along 12 random directions (x₀ + t·d), for P2, P4 and clique #43001:

- **Energy rises in every direction.** SPIQ sits at the bottom of a narrow dip. By t ≈ 1 the energy is on a high plateau: about 10 on P2 (SPIQ 5.1, ground 3.65), 24 on P4 (SPIQ 7.4, ground 4.9), 32 on the clique (SPIQ 13.8, ground 8.9).
- **P(optimal plan) along these straight lines stays low (about 0.05–0.2).** No straight path out of SPIQ leads to a good plan. Good regions have to be *found* by an optimizer.

---

## 7. Step 3. Where the "kick" idea came from, and what a kick does

### Origin

Our earlier study on 1,088 generated problems (`many_problems/SUMMARY.md`) already found that optimizers stall at SPIQ points. Its **Future work, item 4** reads: "Escape the SPIQ stationary point: SPIQ followed by a small random kick, or SPIQ → NFT → BOBYQA." The pilot in `new_formulation_pilot/` started on 5 Nov by testing that idea on the new formulation.

**The reasoning in one sentence:** if the optimizer cannot see a way out because every slope is zero, move it somewhere the slopes are not zero, and let it optimize from there.

### A kick alone always makes things worse

**Figure: `04_kick_alone_always_worse.png`.** 200 kicks at each size, exact energy *right after* the kick (no optimizer yet):

| Kick size σ | P2 median energy (SPIQ 5.12) | P4 (SPIQ 7.42) | Clique #43001 (SPIQ 13.78) | Kicks that land below SPIQ |
|---|---|---|---|---|
| 0.2 | 6.23 | 11.74 | 18.07 | **0%** |
| 0.5 | 8.88 | 20.87 | 27.84 | **0%** |
| 1.0 | 10.39 | 24.00 | 32.53 | **0%** |
| 2.0 | 10.47 | 24.18 | 32.80 | **0%** |

- **Not one of the 13,600 kicks (200 per size × 4 sizes × all 17 SPIQ points; the table shows 3 of them) lowered the energy.** The kick is not an improvement step; it is a jump to a new starting point. **All the improvement has to come from the optimizer afterwards.** That is why the choice of optimizer after the kick is the key decision (Section 8).
- By σ = 1.0 the energy is already as high as at σ = 2.0: the circuit has completely left SPIQ's dip. A smaller kick leaves it on the side of the dip, where an optimizer tends to slide straight back to SPIQ.

---

## 8. Step 4. Why NFT, out of all 10 optimizers

The answer has three parts: **evidence (round 1)**, **mechanism (how NFT works)**, and **confirmation (rounds 2–3)**.

### Part 1, evidence: round 1 kicked five optimizers, and only NFT improved

**Figure: `06_round1_which_optimizer_to_kick.png`** (data: `pilot_round1/runs_P1_P2_P3_b1000_kicks.jsonl`).

Round 1 (5 Nov) first ran all 10 optimizers from SPIQ on P1–P3 (all stuck, as in figure 01). It then kicked **five** of them, NFT, COBYLA, POWELL, ADAM and BOBYQA, at σ = 0.05, 0.2 and 0.5, 5 runs each, 1,000 evaluations. Results on P2 (P1 and P3 are already exact):

**Energy change vs SPIQ (negative = ended lower than SPIQ = better):**

| P2 multi-start (SPIQ 3.741) | no kick | 0.05 | 0.2 | 0.5 |
|---|---|---|---|---|
| **NFT** | +0.00 | +0.00 | **−0.01** | **−0.04** |
| COBYLA | +0.00 | +0.07 | +0.10 | +0.60 |
| POWELL | +0.00 | +0.04 | +0.32 | +1.69 |
| ADAM | +0.00 | +0.07 | +0.20 | +2.24 |
| BOBYQA | +0.00 | +0.09 | +1.11 | +3.85 |

| P2 single-start (SPIQ 5.125) | no kick | 0.05 | 0.2 | 0.5 |
|---|---|---|---|---|
| **NFT** | −0.12 | **−0.45** | **−1.08** | **−0.89** |
| COBYLA | +0.00 | +0.06 | +0.13 | +0.04 |
| POWELL | +0.00 | +0.05 | +0.15 | +0.95 |
| ADAM | −0.00 | +0.07 | +0.14 | +1.16 |
| BOBYQA | −0.00 | +0.07 | +0.79 | +3.19 |

**P(optimal plan), P2 multi-start (SPIQ 0.50):** NFT reaches **0.57** (kick 0.2) and **0.76** (kick 0.5, with 2 of 5 runs above 0.9). Every other kicked optimizer ends at or below 0.50 (down to 0.14 for BOBYQA). On P2 single-start nobody passes 0.09 at this budget.

Across all kicked P2 runs in round 1: **NFT 2 runs above 0.9; COBYLA, POWELL, ADAM, BOBYQA 0** (30 runs each).

**What this shows:**
- After a kick, **NFT is the only optimizer that ends *below* SPIQ's energy**. The other four end *above* it, and the bigger the kick, the worse they end. They cannot recover from the energy a kick adds (Section 7).
- NFT is the only one that turns a kick into a better plan.
- That is the data-driven reason NFT was chosen.

**Why these five, and not all ten?** The run log does not record a written reason for the choice of five. The five cover the main optimizer families: the team's COBYLA, two other model/direction methods (BOBYQA, POWELL), one gradient method (ADAM), and NFT. BOBYQA was also named in the future-work list. The remaining five were already shown to be stuck or harmful from SPIQ in figure 01 (SPSA makes things worse from SPIQ; SLSQP and L_BFGS_B stop early). If asked, say plainly: "we kicked one representative of each optimizer family; NFT was the only one that improved."

### Part 2, mechanism: why NFT is the right tool for these circuits

**NFT = Nakanishi–Fujii–Todo (2019), "Sequential minimal optimization for quantum-classical hybrid algorithms".**

Because the energy is an exact sine wave in each angle (Section 6), three energy values pin down that whole curve, and its lowest point is known exactly. NFT uses this. For angle i (cycling through all angles):

1. z0 = energy at the current angles (reused from the previous step; re-measured every 32 steps so errors don't build up),
2. z1 = energy with angle i turned **+π/2**,
3. z3 = energy with angle i turned **−π/2**,
4. compute the sine wave from z0, z1, z3 and **set angle i to its exact minimum**.

So each step costs **2 evaluations** (3 every 32nd step) and puts that angle at its exact best value, given the others.

Why this beats the others after a kick:
- **No step size.** It jumps to the exact 1-D minimum, so it can make big moves when needed and never overshoots along that angle. COBYLA, BOBYQA and POWELL take cautious local steps and stop early on a local model.
- **Robust to shot noise.** It compares energies ±π/2 apart (big differences). Gradient methods compare nearby points, whose small differences shot noise swamps.
- **It is built for exactly this circuit structure** (one rotation gate per angle).
- **It never steps uphill along the angle it is fixing.** After a kick, every NFT step can only lower (or keep) the energy along that angle, so it walks down from the high plateau the kick landed on.

### Part 3, confirmation: rounds 2 and 3

- **Round 2 (P2, 1,000 evaluations, 200 kicked runs each):** kicked NFT gave 34 runs above 0.9, kicked COBYLA 8. At 3,000 evaluations kicked NFT gave 89 of 200. On P4 at 1,000 evaluations neither succeeded (0 of 30 each).
- **Round 3 (new seeds, 3,000 evaluations, 120 kicked runs each): kicked NFT 42 runs above 0.9, kicked COBYLA 0.** Kicked COBYLA stopped by itself after 1,475–3,000 evaluations without escaping. On P2 multi-start a kick even made COBYLA *worse* (0.50 → 0.33 / 0.28).

**So the claim is precise:** kicked COBYLA occasionally succeeded in round 2 (8 of 200), never in round 3 (0 of 120). Kicked NFT succeeded far more often in both.

---

## 8b. The intuition: what the kick means, and why NFT is the one that works

This section explains the *why* in plain words. **Measured facts** (from the data) and **explanations** (our reasoning about the mechanism) are marked separately, so you never claim more than the data shows.

### The picture to keep in mind

Think of the energy as a landscape and the angles as your position in it. The optimizer is a hiker trying to walk downhill, but it can only feel the ground near its feet (it only sees measured energies).

- **SPIQ puts the hiker on a flat ledge.** Measured: the slope there is zero and every single-angle step goes up or stays flat (Section 6). A hiker who only walks downhill has no step to take, so every optimizer stays put (figure 01).
- **The real valley is elsewhere.** Measured: from SPIQ the energy rises in every straight direction, and P(optimal plan) along those lines stays low (figure 03). The good region is not "just next door".

### What the kick means

**Definition:** kick = add a random number to every angle at once: x_kicked = x₀ + σ·z, with z drawn from a bell curve and σ = 1.0 radian.

**In the picture:** pick the hiker up and drop them at a random spot about one step-size away in every direction at once.

**What the kick is NOT:** an improvement. Measured: 0 of 13,600 kicks landed lower than SPIQ (figure 04). The hiker always lands *higher* up the hill.

**What the kick IS:** a way to put the hiker on **sloped ground**. Measured: slopes near SPIQ are zero; slopes at nearby random points are 1–12 (Section 6). On sloped ground, an optimizer has a direction to walk.

**Why σ = 1.0 and not smaller:** measured, the energy after a σ = 1.0 kick is the same as after σ = 2.0 (figure 04): the hiker has completely left SPIQ's ledge. Explanation: a small kick (0.2) leaves the hiker on the side of SPIQ's dip, and walking downhill from there leads straight back to the same ledge. That is exactly what the failed halving pipeline picked (figure 10).

**Why we first kicked "from SPIQ" and not from anywhere:** the original idea (future-work list) was to stay near SPIQ's good energy and just get off the ledge. The control experiment later showed that the starting neighbourhood does not matter: uniformly random starts work just as well (Section 11). So the kick's real job is just "start somewhere with a slope".

### Why NFT works when the others don't: it is the nature of the optimizer

Yes: the reason is **how NFT is built, and how well that matches these circuits**. Five properties, each with what we measured and what we infer from it.

**1. Its model of the landscape is exact, not approximate.**
- *Measured:* along any single angle, the energy is an exact sine wave, E(θ) = a + b cos θ + c sin θ (fit error 10⁻¹⁵–10⁻¹⁴ at all 17 SPIQ points).
- *How NFT uses it:* three energy values pin down that sine wave completely, so NFT knows the whole curve along that angle, not just the slope at one point.
- *The others:* COBYLA fits a straight-line (linear) model, BOBYQA a curved (quadratic) one, gradient methods a local slope. These are approximations that are only good very close to the current point.

**2. It jumps to the best value along each angle, anywhere in the full circle.**
- *How NFT works:* the minimum of a sine wave can be anywhere from −π to π, and NFT moves the angle straight there, even by a large amount, in one step.
- *Why that matters after a kick:* the kicked point is high on the hill. NFT can make big, exact moves downhill one angle at a time. Local methods take small, careful steps.

**3. It never stops early: it keeps sweeping until the budget is used.**
- *Measured:* kicked COBYLA stopped by itself after 1,475–3,000 evaluations in round 3 (its trust region, started at radius 2.0 here, shrinks until it decides it has converged). From SPIQ without a kick, COBYLA stops after about 1,240 (P2) and 2,900 (P4) evaluations.
- *How NFT works:* it has no "converged, stop" rule of that kind. It cycles through all angles again and again (about 23 full sweeps in 3,000 evaluations on 64 angles), each sweep lowering the energy a little more.

**4. It is robust to shot noise.**
- *How NFT works:* it compares energies measured at ±π/2 apart, which differ a lot wherever the slope is non-zero. Shot noise (10,240 shots) is small compared with those differences.
- *The others:* gradient methods estimate slopes from points very close together. Their energy differences are tiny and shot noise can swamp them.

**5. Along the angle it is fixing, it never goes uphill.**
- *How NFT works:* each step puts the angle at the bottom of its sine wave, so (up to shot noise) the energy along that angle can only go down or stay.

**What the round 1 data shows these properties doing** (figure 06):
- After a kick, NFT ended **below** SPIQ's energy (as much as 1.08 lower on P2 single-start).
- The four others ended **above** SPIQ, and the bigger the kick, the worse (BOBYQA +3.85 at kick 0.5).
- *Explanation:* a kick adds a lot of energy (figure 04). NFT (properties 1–3) can take that back off and then go lower. The others take small steps or stop early, so they only take back part of what the kick added.

### Why NFT alone (without a kick) is not enough

The same properties explain why NFT from SPIQ does almost nothing (figure 01: P2 only −0.18; round 4 "NFT only": P2 single 0.050, P4 0.398):
- At SPIQ, for each angle the bottom of its sine wave **is** the current value (uphill) or the wave is flat. So NFT's exact jump is a jump of zero.
- **The kick and NFT complement each other:** the kick provides a starting point where the sine waves have real height. NFT then exploits that height exactly.

### Why several restarts, not one

- *Measured:* one kick on P4 gave 0.28 (round 3); five kicks with lowest-energy selection gave 0.69–0.76.
- *Explanation:* NFT is a one-angle-at-a-time method. It walks down to the bottom of whichever valley it starts in, but cannot climb out to a deeper one. Different starts land in different valleys. Five tries make it likely that one lands in a good one, and lowest measured energy picks it.

### What to say in one breath

"SPIQ leaves every optimizer on a flat ledge with zero slope. A kick moves the circuit to sloped ground. It always makes the energy worse at first, so everything depends on the optimizer after it. NFT is the one that works because its model is exact for these circuits: the energy is a perfect sine wave in each angle. So it jumps to the exact best value of each angle, makes big moves, is robust to shot noise, and doesn't stop early. Local methods like COBYLA take small steps, stop early, and end higher than SPIQ. One kick is a gamble, so we do five and keep the lowest measured energy. And the control showed the kick doesn't need to start near SPIQ: what matters is NFT plus several fresh starts."

### Be honest about what is measured vs explained

| Statement | Status |
|---|---|
| Slope is zero at SPIQ; no single angle goes downhill | **Measured** (17 points) |
| Energy is an exact sine wave in each angle | **Measured** (fit error ≤ 2.5 × 10⁻¹⁴) |
| A kick always raises the energy | **Measured** (0 of 13,600) |
| After a kick, only NFT ends below SPIQ | **Measured** (round 1, 5 optimizers) |
| Kicked COBYLA stops early | **Measured** (1,475–3,000 evaluations in round 3) |
| NFT works *because* its model is exact and it jumps to each angle's minimum | **Explanation**, consistent with all the measurements above, but we did not run an experiment that isolates each property |
| Gradient methods fail because of shot noise | **Explanation**; not tested separately with a kick |
| Restarts work because each lands in a different valley | **Explanation**; consistent with run-to-run spread (figures 08, 13) |

---

## 9. Step 5. Building the pipeline: rounds 2–5

All rounds: new formulation, the team's saved SPIQ points, 10,240 shots per evaluation, depth 2 (two QAOA layers), best-seen point scored exactly.

### Round 2 (5 Nov): kick size sweep (exploratory)

**Figure: `07_round2_kick_size_sweep.png`.**

- **P2 multi-start:** NFT is best at **σ = 0.5** (0.83 at 1,000 evaluations; **0.974** with 18/20 above 0.9 at 3,000).
- **P2 single-start:** NFT keeps improving up to **σ = 1.0** (0.60 at 3,000 evaluations, 10/20).
- **P4:** only NFT with σ = 1.0 at 3,000 evaluations helped (0.574, 5/10).
- **Lesson: the best kick size depends on the start.** A kick that is too small does not leave the dip.
- **These numbers are not claimed**, because they were used to choose settings. And the P4 0.574 **did not replicate** in round 3 (0.276). Say this openly.

### Round 3 (5–6 Nov): one kick, then 3,000 evaluations, 20 new seeds

**Figure: `08_round3_kick_needs_NFT.png`** (bars = mean, dots = each run).

| Start | COBYLA | NFT | COBYLA + kick 0.5 | NFT + kick 0.5 | COBYLA + kick 1.0 | NFT + kick 1.0 |
|---|---|---|---|---|---|---|
| P2 multi-start | 0.500 (0/20 > 0.9) | 0.501 (0) | 0.334 (0) | **0.950 (18)** | 0.275 (0) | 0.531 (8) |
| P2 single-start | 0.000 (0) | 0.000 (0) | 0.065 (0) | 0.142 (0) | 0.187 (0) | **0.653 (12)** |
| P4 single-start | 0.250 (0) | 0.296 (0) | 0.164 (0) | 0.229 (2) | 0.188 (0) | 0.276 (2) |

Statistical test (one-sided Mann–Whitney, kicked NFT vs COBYLA): P2 multi kick 0.5 p = 3.4 × 10⁻⁸; P2 single kick 1.0 p = 5.8 × 10⁻⁵; **P4: p = 0.98 and 0.86, no gain.**

What round 3 taught:
1. **A kick needs NFT after it.**
2. **The best kick size depends on the start** (0.5 vs 1.0).
3. **One kick is a gamble.** On P4 some kicks land well and many badly, so no single-kick setting helped. Hence: try **several kicks** and keep the best.

### Round 4 (6 Nov): six ways to spend the same 15,000 evaluations (planned in advance)

Written in `restarts.py` **before** any round 4 result was seen. Every pipeline gets the same budget and chooses by measured energy only.

| Pipeline | What it does |
|---|---|
| COBYLA | team's optimizer from SPIQ, no kick (stops by itself) |
| NFT | NFT from SPIQ, no kick, all 15,000 evaluations |
| kick 0.5 ×5 | 5 restarts × 3,000, kick 0.5 |
| kick 1.0 ×5 | 5 restarts × 3,000, kick 1.0 |
| mixed ×5 | 5 restarts, one each at 0.2, 0.35, 0.5, 0.75, 1.0 |
| **halving (pre-registered primary)** | 15 short runs of 400 (3 per kick size), then the 3 lowest continue for 3,000 each |

"Halving" was the primary because it needs no kick size: it was supposed to discover the right one.

**Figure: `09_round4_six_pipelines.png`.**

| Start | COBYLA | NFT | kick 0.5 ×5 | **kick 1.0 ×5** | mixed ×5 | halving |
|---|---|---|---|---|---|---|
| P2 multi | 0.500 (0) | 0.904 (16) | 0.999 (20) | **0.998 (20)** | 0.999 (20) | 0.998 (20) |
| P2 single | 0.000 (0) | 0.050 (1) | 0.533 (9) | **0.949 (19)** | 0.768 (13) | 0.786 (15) |
| P4 single | 0.250 (0) | 0.398 (1) | 0.337 (3) | **0.761 (13)** | 0.603 (9) | 0.196 (1) |

(brackets = runs out of 20 above 0.9.) **Kick 1.0 ×5 is the only pipeline at or near the top on all three.** One long NFT run without a kick (column "NFT") stays stuck on P2 single (0.050) and P4 (0.398).

### Why the pre-registered pipeline failed

**Figure: `10_round4_halving_failure.png`.** On P4, halving picked **kick 0.2 in 19 of 20 runs** (0.35 in the other) and ended at **0.196, worse than COBYLA**. Its 400-evaluation screen kept the runs with the lowest energy after 400 evaluations. Small kicks start close to SPIQ, so after 400 evaluations they look best: they slide back into SPIQ's dip. But they stay stuck there, while kick-1.0 runs start high and end much lower after 3,000. **A short screen rewards the wrong thing.** We report this as a failed pre-registered pipeline.

### Round 5 (6 Nov): the main claim, on seeds never seen before

Kick 1.0 ×5 was picked *after* seeing round 4, so it had to be re-tested on fresh seeds: reps 101–120, 20 runs per start.

**The final pipeline (figure `05_pipeline_diagram.png`):**
1. Measure SPIQ's point x₀ once (candidate 0).
2. Five times: kick x₀ + 1.0·N(0, I), run NFT for 3,000 evaluations, keep its lowest measured point (candidates 1–5).
3. Return the candidate with the lowest measured energy (15,001 evaluations in all).
4. Score it exactly (for the report only).

**Figure: `11_round5_confirmation_unseen_seeds.png`.**

| Start | COBYLA P(opt) | **Kick 1.0 ×5 + NFT P(opt)** | Runs > 0.9 (95% interval) | COBYLA energy | Kick + NFT energy | Ground | p-value |
|---|---|---|---|---|---|---|---|
| P2 multi-start | 0.500 | **0.977** | 19/20 (0.75–1.00) | 3.742 | 3.660 | 3.653 | 3.4 × 10⁻⁸ |
| P2 single-start | 0.000 | **0.973** | 19/20 (0.75–1.00) | 5.127 | 3.661 | 3.653 | 3.4 × 10⁻⁸ |
| P4 single-start | 0.250 | **0.693** | 11/20 (0.32–0.77) | 7.425 | 5.190 | 4.908 | 6.1 × 10⁻⁴ |
| P1, P3 | 1.000 | 1.000 | 20/20 | at ground | at ground | | no change |

- **It replicated:** P2 single 0.949 → 0.973, P4 0.761 → 0.693 from round 4 to round 5.
- **P4 is not solved.** 11 of 20 runs succeed and all 11 have exact energy ≤ 4.943 (ground 4.908). The other 9 end at 5.095–6.219 on mixed or wrong plans. **When the energy gets within about 0.035 of the ground state, the plan is right.**
- **It never hurts P1 or P3.** SPIQ's point (candidate 0) was the chosen answer in all 20 runs on each, because no restart beat its energy. On P2 and P4 SPIQ's point was never chosen (0 of 20): a restart always beat it.

---

## 10. Step 6. The 12 new problems

Same pipeline, **nothing re-tuned**, 3 runs per problem with new seeds, both formulations (`data/runs_new.jsonl`).

### Per problem

**Figure: `12_new_problems_kick_heatmap.png`.** Rows = problems; columns = pipelines; big number = mean P(optimal plan); small numbers = each of the 3 runs.

| Problem | SPIQ start | SPIQ + COBYLA | SPIQ + NFT | **Kick ×5 + NFT** | Random ×5 + NFT (control) | Old formulation kick ×5 + NFT |
|---|---|---|---|---|---|---|
| CHAIN #40000 | 0.50 | 0.50 | 1.00 | **1.00** | 1.00 | 1.00 |
| CHAIN #40001 | 1.00 | 1.00 | 1.00 | **1.00** | 1.00 | 1.00 |
| CHAIN #40003 | 0.00 | 0.00 | 0.48 | **1.00** | 1.00 | 0.79 |
| STAR #41000 | 0.25 | 0.25 | 0.50 | **1.00** | 0.84 | 1.00 |
| STAR #41003 | 0.50 | 0.50 | 0.50 | **0.52** (0.55, 1.00, 0.00) | 0.90 | 0.00 |
| STAR #41004 | 0.50 | 0.50 | 0.50 | **0.51** (0.50, 0.50, 0.52) | 0.83 | 0.82 |
| CYCLE #42000 | 0.00 | 0.00 | 0.00 | **1.00** | 1.00 | 1.00 |
| CYCLE #42001 | 0.50 | 0.50 | 0.50 | **1.00** | 1.00 | 0.75 |
| CYCLE #42002 | 0.12 | 0.13 | 0.10 | **1.00** | 1.00 | 0.93 |
| CLIQUE #43000 | 0.50 | 0.50 | 0.50 | **0.83** (1.00, 1.00, 0.50) | 1.00 | 0.20 |
| CLIQUE #43001 | 0.00 | 0.00 | 0.00 | **1.00** | 1.00 | 0.02 |
| CLIQUE #43002 | 0.50 | 0.50 | 0.50 | **1.00** | 1.00 | 0.03 |
| **Mean** | **0.365** | **0.365** | **0.465** | **0.904** | **0.963** | **0.628** |

COBYLA and NFT here use 1,000 evaluations; kick and random use 15,001.

**Reading it:**
- **COBYLA adds nothing** (0.365 → 0.365): the same flat-spot story on 12 problems it has never seen.
- **Kick + NFT solves 9 of 12** (≥ 0.99): all chains, all cycles, STAR #41000, CLIQUE #43001 and #43002.
- **The three misses, explained:**
  - **CLIQUE #43000 (0.83):** two runs reach 0.999; one stays at 0.50.
  - **STAR #41003 (0.52):** runs end at 0.55, 0.996 and 0.00. In the third run a restart found a point with *lower energy than SPIQ* (14.81 vs 15.74; ground 14.57) but on a wrong plan, and lowest-energy selection kept it. Below the ground state's neighbourhood, lower energy does not guarantee the right plan.
  - **STAR #41004 (0.51):** all three runs end at energy 11.778 (ground 11.761), half on the optimal plan and half on one competing plan only **0.017** higher in energy. A near-tie the optimizer could not separate.

### Every run over time

**Figure: `13_new_problems_every_run.png`.** One panel per problem; green lines = each of 3 kick + NFT runs (P(optimal plan) of the best point so far), dotted vertical lines = each new kick, diamonds = COBYLA (orange) and NFT (blue) from SPIQ at 1,000 evaluations. Most runs jump to 1.0 within the first one or two kicks. Lines that drop (STAR #41000, #41003, #41004) are the cases where a lower-energy point on a different plan was found and kept.

### Averaged trends and the equal-budget question

**Figure: `14_new_problems_trend_vs_COBYLA_NFT.png`.** Four panels: old/new formulation × standard/SPIQ start, COBYLA (orange) vs NFT (blue), over 1,000 evaluations. COBYLA's line is flat in every panel. The fifth panel shows kick + NFT over 15,000 evaluations: new formulation (blue) rises to 0.904, old (orange) to 0.628.

**Is it just the bigger budget?** Compare at the same 1,000 evaluations:

| At 1,000 evaluations | Mean P(optimal plan), 12 new problems |
|---|---|
| COBYLA from SPIQ | 0.365 |
| NFT from SPIQ | 0.465 |
| Kick + NFT (first kick, 999 NFT evaluations so far) | **0.660** |

**At the same budget kick + NFT is already far ahead.** It reaches 0.765 by 3,000 evaluations, and 0.904 after five kicks.

### The old formulation

**Figure: `17_old_formulation_ground_state_wrong_plan.png`.** One dot per kick + NFT run. x = probability on the lowest-energy state; y = P(optimal plan). Blue dots = problems where the lowest-energy state is the optimal plan; orange crosses = problems where it is not.

- **New formulation (right):** the lowest-energy state is the optimal plan on all 12, so dots lie on the diagonal: reaching the ground state = getting the plan.
- **Old formulation (left):** on 4 of the 12 (STAR #41003 and the three cliques) the lowest-energy state is a **wrong plan**. Orange crosses sit at the bottom: the better the optimizer does at lowering the energy, the further it moves from the optimal plan. STAR #41003 reached 99.5% on the lowest-energy state with **0.000** on the optimal plan.
- So kick + NFT reaches 0.628 in the old formulation (≥ 0.99 on 4 of 12) vs 0.904 in the new one (9 of 12). On the cliques it ends at 0.20, 0.03 and 0.02, vs 0.83, 0.999 and 0.999 with the new formulation. **No optimizer can fix a formulation whose minimum is the wrong answer.**

---

## 11. Step 7. Control experiment: is it the kick, or just restarting?

### Why we ran it

Section 7 showed that a σ = 1.0 kick throws the circuit about as far as σ = 2.0, so it has completely left SPIQ's dip. A fair question follows: **is starting near SPIQ doing anything, or would fully random restarts work just as well?** We tested it.

### What was changed (exactly one thing)

Same SPIQ circuit, SPIQ measured once as candidate 0, 5 restarts × 3,000 NFT evaluations, same lowest-measured-energy selection, same 15,001 evaluations, same scoring, new seeds. **Only change:** each restart starts from angles drawn **uniformly at random in [−π, π)** instead of SPIQ + 1.0·N(0, I). (`scripts/ablation_random.py`, `data/runs_ablation_random.jsonl`.)

### Result

**Figure: `15_random_restart_control.png`.**

| Case (runs) | COBYLA from SPIQ | Kick 1.0 ×5 + NFT | **Random ×5 + NFT** |
|---|---|---|---|
| P2 multi-start (20) | 0.500 (0 > 0.9) | 0.977 (19) | **0.999 (20)** |
| P2 single-start (20) | 0.000 (0) | 0.973 (19) | **0.975 (19)** |
| P4 single-start (20) | 0.250 (0) | 0.693 (11) | **0.712 (11)** |
| 12 new problems (36) | 0.365 (3) | 0.904 (30) | **0.963 (33)** |

- **Random restarts do as well as kicks, slightly better on the new problems.** Differences of 0.02–0.06 are within run-to-run noise at these sample sizes.
- In only **1 of 76** control runs did SPIQ's own point stay the best candidate.

### What drives the gain

| Ingredient | Needed? | Evidence |
|---|---|---|
| **Leaving SPIQ's flat spot** | **Yes** | From SPIQ without leaving, nothing moves (figure 01; round 4 "NFT" column) |
| **NFT as the optimizer** | **Yes** | Round 1: only kicked NFT improved; round 3: kicked COBYLA 0 of 120 above 0.9 |
| **Several restarts, chosen by measured energy** | **Yes** | One kick on P4: 0.28 (round 3); five: 0.69–0.76 |
| **Restarts centred on SPIQ's angles** | **No (on these problems)** | Random restarts match or beat kicks |
| **SPIQ's circuit (one angle per gate)** | Not tested separately | Every run uses it |
| **SPIQ's point as candidate 0** | A safety net | P1 and P3 stay at 1.00 |

**So the honest name for the method is "multi-start NFT": several NFT runs from fresh starting points, keeping the lowest measured energy.** A kick from SPIQ is one good way to make the fresh starts; a uniformly random start is an equally good one. SPIQ's value in this pipeline is the circuit layout and the safe candidate 0, not the neighbourhood of its angles.

**Say this yourself before you are asked.** It turns the obvious objection ("isn't a 1-radian kick just a random restart?") into evidence that you tested your own method.

---

## 12. Step 8. P5, the 24-qubit scale test (interim)

From SPIQ, every one of the 10 optimizers leaves P5 at P(optimal plan) = **0.00** (figure 01). Four kick + NFT runs with new seeds are running (`log_p5_kick.txt`). P5's state has 16.8 million amplitudes; each evaluation took about 5–10 seconds on the laptop. The kicks took 8.4, 5.9, 5.6 and 4.0 hours.

**Figure: `16_P5_interim.png`.** Left: exact energy of each run's best candidate so far. Right: its P(optimal plan). Restart 0 = SPIQ.

Results after **4 of 5 kicks** (each cell: where that restart ended, energy / P(optimal plan); **bold** = kept as the run's best):

| Run | Kick 1 | Kick 2 | Kick 3 | Kick 4 | Kept best |
|---|---|---|---|---|---|
| 1 | **12.596 / 0.003** | 15.131 / 0.002 | 14.536 / 0.021 | 14.426 / 0.000 | 12.596 / 0.003 |
| 2 | **11.086 / 0.380** | 17.551 / 0.160 | 17.597 / 0.120 | **10.316 / 0.039** | 10.316 / **0.039** |
| 3 | 17.765 / 0.000 (SPIQ kept) | **9.768 / 0.976** | 15.831 / 0.024 | 13.167 / 0.293 | 9.768 / **0.976** |
| 4 | **12.000 / 0.612** | 17.445 / 0.144 | 12.683 / 0.415 | 12.477 / 0.057 | 12.000 / **0.612** |

SPIQ start: 17.558 / 0.000. Ground: 8.340.

- **It escapes on P5:** run 3 reaches **0.976**, run 4 **0.612**, where every optimizer from SPIQ gives 0.00. Every run's best energy is far below SPIQ's.
- **Several restarts matter:** run 3's first kick found nothing better than SPIQ; its second solved it.
- **Worse restarts never overwrite better ones** (runs 1, 3, 4 after kick 2).
- **The limit (run 2, kick 4):** a lower energy (10.316 vs 11.086) was kept, and P(optimal plan) fell from 0.380 to 0.039. Energy is an average over plans; above the ground state, a lower average can come from favouring the 2nd- and 3rd-best plans. Both points are far from the ground state (8.34). Run 3 got close (9.77) and there energy and plan agree. **The pipeline can only see energy, so on a problem this large an unconverged run can trade one for the other. The fix is more budget per restart, not a different rule, because P(optimal plan) cannot be measured without knowing the answer.** We have not yet checked which plans run 2's point favours.

**Status:** kick 5 is running. Present P5 as "after 4 of 5 kicks, interim".

---

## 13. What it costs

| | COBYLA from SPIQ | Kick 1.0 ×5 + NFT |
|---|---|---|
| Evaluations per run (round 5, median) | 1,237–1,242 (P2), 2,907 (P4); stops by itself | 15,001 |
| Shots per run (× 10,240) | about 12.7 million (P2) | about **154 million** |
| Ratio | | about **12× (P2), 5× (P4)** |
| Laptop time per run (round 5, median) | seconds | 23–24 s (P2), 172 s (P4); P5 about 28 h (4 of 5 kicks took 24 h) |

On real hardware 154 million shots per run is expensive. But the equal-budget comparison (Section 10) shows the gain is not only from spending more: at 1,000 evaluations it already beats both COBYLA and NFT from SPIQ.

---

## 14. Final conclusions and what they imply

### Conclusions

1. **From SPIQ, optimizer choice does not matter: all 10 are stuck** (figure 01). The team's COBYLA never improves on SPIQ.
2. **The cause is measured, not guessed:** at all 17 SPIQ points the gradient is zero to 10⁻¹⁴ and no single angle goes downhill (figure 02).
3. **A kick is needed to leave the flat spot, but a kick alone always makes things worse** (0 of 13,600 kicks improve). Everything depends on the optimizer after it.
4. **NFT is that optimizer.** Of five kicked optimizers, only NFT ended below SPIQ's energy (round 1). Kicked COBYLA got 0 of 120 runs above 0.9 in round 3; kicked NFT got 42. NFT matches the circuit's exact sine-wave structure.
5. **Several restarts are needed:** one kick is a gamble (P4: 0.28). Five, with lowest-measured-energy selection, make it reliable.
6. **The pre-registered pipeline (halving) failed**, for a clear reason: short screens favour small kicks.
7. **Main result (round 5, unseen seeds):** P2 0.50 → **0.977** and 0.00 → **0.973**, P4 0.25 → **0.693**, all highly significant; P1 and P3 kept at 1.00.
8. **It generalises:** 12 new problems 0.365 → **0.904** (9 of 12 solved), ahead even at equal budget.
9. **Control experiment:** fully random restarts do as well. The method is **multi-start NFT**; SPIQ contributes the circuit and a safe candidate 0, not a good neighbourhood.
10. **Scale (interim):** on P5 two of four runs reach 0.976 and 0.612 from 0.00.
11. **Formulation matters more than the optimizer:** in the old formulation the same pipeline reaches only 0.628 and fails where that formulation's minimum is the wrong plan.

### What the results imply

- **For the team's pipeline:** replace "SPIQ + COBYLA" with "SPIQ + multi-start NFT" (keep SPIQ's point as a candidate). That is the change with a measured effect.
- **For how to think about SPIQ:** SPIQ gives a low-energy, cheap-to-find point, but it is a flat spot. It is a good *answer candidate* and circuit layout, not a good place to *continue optimizing from*.
- **For optimizer studies:** comparing optimizers from SPIQ points measures nothing, since all of them stay put. Comparisons need restarts or starts away from the flat spot.
- **For the formulation question:** use the new formulation. Its lowest-energy state is the optimal plan on 16 of 16 problems, so lowering the energy helps. In the old one it can actively hurt.
- **For scaling:** budget per restart has to grow with the number of angles (6 sweeps on P5 vs 23 on P2). The P5 run 2 case shows what happens when it doesn't.
- **For hardware:** the selection uses only measured energy, so the method can run unchanged on a real device. The cost (about 154 million shots per run) is the obstacle.

---

## 15. Limits (say these before you are asked)

- **Noiseless simulation.** Shot noise is included; hardware noise is not.
- **Settings chosen on P2 and P4** (σ = 1.0, 5 × 3,000), then confirmed on unseen seeds and 12 new problems. Not shown to be the best setting for every size.
- **P4 is only partly solved** (11 of 20 runs).
- **Near-ties** (STAR #41004: a competing plan 0.017 above the ground state) are not separated.
- **Lowest-energy selection can pick a wrong plan when runs stop short of the ground state** (STAR #41003; P5 run 2).
- **Cost:** about 5–12× COBYLA's evaluations.
- **P5 is interim** (4 runs, 4 of 5 kicks).
- **12 new problems × 3 runs** is enough for averages, but per-problem results with 3 runs are noisy.
- **Only 5 of 10 optimizers were kicked** in round 1; the other 5 were not tested with a kick.
- **SPIQ's circuit was not tested separately** from the restarts.

---

## 16. Figure-by-figure guide

Format for each: **what is plotted → how to read it → key numbers → the one-sentence takeaway.**

**01 `01_all10_optimizers_stuck_at_SPIQ.png`**
- *Plotted:* two heatmaps, rows = 10 optimizers, columns = P1–P5 (with qubits, angles and the SPIQ start's P(optimal plan) under each). Left: P(optimal plan) after 1,000 evaluations (white 0 → dark blue 1). Right: energy change from SPIQ (blue = lower, red = higher, grey 0 = did not move).
- *Read:* compare each left column with its "start" value: they are identical. On the right, everything is grey except SPSA (red).
- *Numbers:* nine optimizers keep 1.00 / 0.00 / 1.00 / 0.25 / 0.00 exactly; SPSA +0.07 to +7.75 energy.
- *Takeaway:* "From SPIQ, no optimizer moves, so the optimizer is not the issue: the starting point is."

**02 `02_spiq_flat_spot.png`**
- *Plotted:* left, for each of 17 SPIQ points, a stacked bar of angles that go uphill both ways (orange), flat (grey), downhill (green) when turned ±π/2. Right, gradient size (log scale) at SPIQ (red diamond) vs a nearby point (blue dot).
- *Read:* no green anywhere; red diamonds at 10⁻¹⁵–10⁻¹⁴, blue dots at 1–12.
- *Takeaway:* "SPIQ's point is a flat spot: zero slope and no single angle goes downhill."

**03 `03_landscape_around_SPIQ.png`**
- *Plotted:* top row exact energy, bottom row P(optimal plan), along 12 random straight lines from SPIQ (x-axis t = distance, same scale as the kick σ). Panels P2, P4, a clique. Yellow band = t = 1.0.
- *Read:* energy rises steeply then flattens on a plateau; P(opt) stays low.
- *Takeaway:* "SPIQ sits in a narrow dip; no straight line out of it leads to a good plan, so an optimizer has to find one."

**04 `04_kick_alone_always_worse.png`**
- *Plotted:* distribution of exact energy right after 200 kicks per size σ, for P2, P4, a clique. Orange line = SPIQ energy, dashed = ground.
- *Read:* every distribution sits above the orange line; σ = 1 and σ = 2 look the same.
- *Takeaway:* "A kick never improves anything by itself (0%); it is a jump, and the optimizer after it does the work."

**05 `05_pipeline_diagram.png`**
- *Plotted:* the flow chart: SPIQ point → measure (candidate 0) → 5 × [kick → NFT 3,000] → candidates → pick lowest measured energy → exact scoring.
- *Takeaway:* "Five fresh NFT runs plus SPIQ's own point; the lowest measured energy wins; it never sees the answer."

**06 `06_round1_which_optimizer_to_kick.png`**
- *Plotted:* top row, energy change vs SPIQ for 5 optimizers × (no kick, kick 0.05, 0.2, 0.5) on P2 (left multi-start, right single-start); blue = below SPIQ. Bottom row, P(optimal plan), same layout.
- *Read:* the NFT row is the only blue row; other rows turn deeper red as the kick grows. Bottom left: NFT reaches 0.76 at kick 0.5, others fall.
- *Numbers:* NFT −1.08 (single, 0.2), −0.04 (multi, 0.5); BOBYQA +3.85; NFT 2 runs above 0.9, others 0.
- *Takeaway:* "We kicked five optimizers; only NFT could turn a kick into an improvement, so we chose NFT."

**07 `07_round2_kick_size_sweep.png`**
- *Plotted:* mean P(optimal plan) vs kick size σ, three panels (P2 multi, P2 single, P4). Orange = COBYLA 1,000 evals, solid blue = NFT 1,000, dashed blue = NFT 3,000. Labels = runs above 0.9.
- *Read:* dashed blue is highest; its peak is at 0.5 for P2 multi, 1.0 for P2 single and P4.
- *Takeaway:* "Kicked NFT climbs, kicked COBYLA doesn't, and the best kick size depends on the start. Exploratory only, not claimed."

**08 `08_round3_kick_needs_NFT.png`**
- *Plotted:* bars = mean P(optimal plan), dots = each of 20 runs, for COBYLA and NFT with no kick / kick 0.5 / kick 1.0, on P2 multi, P2 single, P4. New seeds.
- *Read:* orange dots never reach the top; blue dots cluster near 1 for P2 at the right kick; P4 dots spread.
- *Takeaway:* "On new seeds, a kick only works with NFT, and one kick is a gamble on P4."

**09 `09_round4_six_pipelines.png`**
- *Plotted:* bars = mean P(optimal plan), dots = 20 runs, six pipelines, three starts, same 15,000-evaluation budget.
- *Read:* the "kick 1.0 ×5" bar is at or near the top in every panel; halving collapses on P4.
- *Takeaway:* "With equal budgets, five kicks of 1.0 is the only pipeline that wins everywhere."

**10 `10_round4_halving_failure.png`**
- *Plotted:* for halving's 20 runs per start, which kick size the returned answer came from.
- *Read:* on P4, almost all bars at kick 0.2 (19 of 20).
- *Takeaway:* "The pre-registered pipeline failed because a short screen rewards small kicks that slide back to SPIQ."

**11 `11_round5_confirmation_unseen_seeds.png`**
- *Plotted:* top, P(optimal plan) of every run (dots) for COBYLA vs kick 1.0 ×5 + NFT; bottom, exact energy of every run, with the ground state. Black bar = mean.
- *Read:* P2: 19 of 20 kick dots at the top vs COBYLA at 0.5 / 0.0; P4: 11 of 20 at the top.
- *Takeaway:* "On seeds never used to design it: P2 0.977 / 0.973, P4 0.693, vs COBYLA 0.50 / 0.00 / 0.25."

**12 `12_new_problems_kick_heatmap.png`**
- *Plotted:* rows = 12 new problems; columns = SPIQ start, SPIQ + COBYLA, SPIQ + NFT, kick + NFT, random restarts + NFT (control), old-formulation kick + NFT. Big number = mean; small numbers = each of 3 runs. Column labels give the mean over the 12.
- *Read:* the first three columns are pale (stuck); the kick and random columns are mostly dark (solved); the old-formulation column is pale on the cliques.
- *Takeaway:* "On 12 unseen problems: 0.365 → 0.904; random restarts 0.963; old formulation 0.628 and fails on the cliques."

**13 `13_new_problems_every_run.png`**
- *Plotted:* one panel per problem; each green line = one kick + NFT run's best-so-far P(optimal plan) over 15,000 evaluations; dotted verticals = new kick; diamonds = COBYLA / NFT from SPIQ at 1,000.
- *Read:* most lines jump to 1.0 within one or two kicks; falling lines are lower-energy wrong-plan picks.
- *Takeaway:* "Most problems are solved within the first two kicks; the misses are near-ties or wrong-plan picks."

**14 `14_new_problems_trend_vs_COBYLA_NFT.png`**
- *Plotted:* average P(optimal plan) over evaluations for COBYLA (orange) and NFT (blue) in four conditions (old/new × standard/SPIQ); fifth panel kick + NFT for both formulations.
- *Read:* COBYLA is flat everywhere; kick + NFT rises to 0.90 (new) vs 0.63 (old).
- *Takeaway:* "COBYLA never improves; kick + NFT does, and much more in the new formulation."

**15 `15_random_restart_control.png`**
- *Plotted:* bars = mean, dots = every run; COBYLA from SPIQ vs kick ×5 + NFT vs random ×5 + NFT, for P2 multi, P2 single, P4, 12 new problems.
- *Read:* the kick and random bars are about equal; both far above COBYLA.
- *Takeaway:* "Random restarts do as well: the gain comes from NFT plus several restarts, not from staying near SPIQ."

**16 `16_P5_interim.png`**
- *Plotted:* each of 4 runs, best-so-far exact energy (left) and P(optimal plan) (right) after each kick; dashed = ground; orange line on the right = every optimizer from SPIQ (0.00).
- *Read:* energies drop far below SPIQ; run 3 jumps to 0.976 at kick 2; run 2 drops to 0.039 at kick 4.
- *Takeaway:* "At 24 qubits it escapes too (0.976, 0.612), but an unconverged run can lower energy while losing the plan."

**17 `17_old_formulation_ground_state_wrong_plan.png`**
- *Plotted:* one point per kick + NFT run on the new problems; x = probability on the lowest-energy state, y = P(optimal plan); blue = lowest-energy state is the optimal plan, orange × = it is not. Left old, right new formulation.
- *Read:* right: all on the diagonal; left: orange crosses at the bottom, one at x ≈ 1, y = 0 (STAR #41003).
- *Takeaway:* "In the old formulation, getting to the lowest energy can mean getting the wrong plan, so no optimizer can fix it."

---

## 17. Likely questions and answers

**Q1. Why did all 10 optimizers fail from SPIQ?**
SPIQ's point is a flat spot: we measured the gradient as zero (10⁻¹⁵–10⁻¹⁴) at all 17 SPIQ points, and turning any single angle goes uphill or stays flat. Slope-followers see no slope; one-angle methods see no downhill angle; SPSA's random steps almost all go uphill.

**Q2. Why is SPIQ's point a flat spot?**
SPIQ searches Clifford settings (multiples of π/2) and keeps the best. Turning one angle by ±π/2 lands on another Clifford setting that SPIQ could have chosen, so it is no better. The exact zero gradient we measured; we don't claim it is proven in general.

**Q3. Where did the kick idea come from?**
From the future-work list of our earlier 1,088-problem study ("SPIQ followed by a small random kick"). The logic: if the optimizer sees no slope, move it to where there is one.

**Q4. Why NFT and not one of the other 9?**
Evidence first: round 1 kicked five optimizers (NFT, COBYLA, POWELL, ADAM, BOBYQA). Only NFT ended below SPIQ's energy; the other four ended higher, more so with bigger kicks, and only NFT had runs above 0.9. Mechanism: the energy is an exact sine wave in each angle, and NFT jumps to the exact minimum of that wave with ±π/2 probes, so it can walk down from the high-energy point a kick creates. Confirmation: round 3, new seeds, kicked NFT 42 of 120 runs above 0.9, kicked COBYLA 0 of 120.

**Q4b. Is NFT's success because of the nature of the optimizer?**
Yes, as far as we can tell. Its model is exact for these circuits (the energy is a perfect sine wave in each angle, measured), so it jumps to each angle's exact best value. It makes big moves, is robust to shot noise and does not stop early. COBYLA's local model and shrinking trust region make it stop early (1,475–3,000 evaluations), so after a kick it ends higher than SPIQ. The mechanism is our explanation; the measured facts behind it are listed in Section 8b.

**Q5. Why weren't all 10 kicked?**
We kicked one representative per family. The other five were already stuck or harmful from SPIQ (SPSA moves uphill; SLSQP and L-BFGS-B stop early). It is a limitation that they were not kicked; NFT's advantage is explained by its mechanism, which none of the others share.

**Q6. What is a seed?**
The starting number of the random generator. Same seed → identical run; different seeds → independent repeats. Ours are computed from each run's label (CRC-32), and the main claim (round 5) uses run numbers 101–120 that were never used to choose anything.

**Q7. Why σ = 1.0?**
Smaller kicks leave the circuit on the side of SPIQ's dip, where NFT slides back to SPIQ (round 2; halving failure). At σ = 1.0 the circuit has fully left the dip (energy after the kick equals that at σ = 2.0). In round 4, kick 1.0 ×5 was the only pipeline at or near the top on all three starts.

**Q8. Why 5 restarts of 3,000, not 1 of 15,000 or 15 of 1,000?**
One long NFT run stays stuck (round 4: P2 single 0.050, P4 0.398). Many short runs (halving) pick the wrong kicks. 5 × 3,000 gives about 23 sweeps over 64 angles per restart, enough to converge, and 5 tries to survive bad kicks.

**Q9. Isn't a 1-radian kick just a random restart?**
Yes, and we tested it: uniformly random restarts + NFT gave 0.999 / 0.975 (P2), 0.712 (P4), 0.963 (new problems), the same as kicks. The method is multi-start NFT.

**Q10. Then what is SPIQ still good for?**
Its circuit (one angle per gate), which every successful run uses; its point as the safe candidate 0 (P1, P3 stay at 1.00); and on its own it beats 1,000 evaluations from the standard start. It is not a good place to continue optimizing from.

**Q11. How do you know the pipeline doesn't cheat?**
It chooses only by measured 10,240-shot energy. P(optimal plan) and the ground state are computed afterwards for the report. It could run unchanged on hardware.

**Q12. Is the improvement just from a bigger budget?**
No. At equal budget (1,000 evaluations) on the 12 new problems: kick + NFT 0.660 vs NFT 0.465 vs COBYLA 0.365.

**Q13. Why did the pre-registered halving pipeline fail?**
Its 400-evaluation screen kept the lowest early energies, which come from small kicks sliding back to SPIQ. Those stay stuck. On P4 it picked kick 0.2 in 19 of 20 runs and ended at 0.196.

**Q14. Doesn't a failed pre-registered pipeline weaken the claim?**
It makes the process more credible. We report the failure, explain it, and test the replacement on unseen seeds (round 5) and 12 unseen problems.

**Q15. Why is P4 only 0.693?**
9 of 20 runs end at 5.095–6.219, on mixed or wrong plans; all 11 successes are within 0.035 of the ground state. P4 has 140 angles, so 3,000 evaluations is only about 10.5 sweeps. Longer restarts might help; not tested yet.

**Q16. Why did STAR #41003 get a run at 0.00?**
A restart found a point with lower energy than SPIQ (14.81 vs 15.74) but on a wrong plan; ground is 14.57. Lowest-energy selection kept it. Away from the ground state, lower energy doesn't guarantee the right plan.

**Q17. What about STAR #41004?**
All runs end at 11.778 vs ground 11.761, split between the optimal plan and a competitor only 0.017 higher. A near-tie.

**Q18. P5 run 2: lower energy but P(optimal plan) fell from 0.38 to 0.04. Is the method broken?**
No, it shows its limit. The pipeline can only see energy. Above the ground state, a lower average energy can come from favouring 2nd- and 3rd-best plans. Run 3 got close to the ground state (9.77 vs 8.34) and there energy and plan agree (0.976). The remedy is more budget per restart on large problems.

**Q19. Does it work in the old formulation?**
Partly: 0.383 → 0.628 on the 12 new problems (vs 0.904 new). It cannot fix problems where the old formulation's lowest-energy state is the wrong plan (cliques end at 0.20 / 0.03 / 0.02).

**Q20. How significant are the results?**
Round 5, one-sided Mann–Whitney vs COBYLA: P2 p = 3.4 × 10⁻⁸ (both starts), P4 p = 6.1 × 10⁻⁴. 19/20 runs above 0.9 on P2 (95% interval 0.75–1.00).

**Q21. Does it ever make things worse than SPIQ?**
Not by measured energy: SPIQ's point is candidate 0. On P1 and P3 it kept 1.00 in all 20 runs. P(optimal plan) can fall below SPIQ's only if a lower-energy wrong-plan point is found (STAR #41003, one run).

**Q22. How expensive is it?**
15,001 evaluations (about 154 million shots) vs about 1,240 (P2) and 2,900 (P4) for COBYLA: about 12× and 5×. On the laptop: 24 s (P2), 172 s (P4), about 28 h (P5) per run.

**Q23. Why is P(optimal plan) so often exactly 0.50 or 0.25 at SPIQ?**
SPIQ points are Clifford circuits, whose output probabilities are multiples of 1/2ᵏ. So the probability of the optimal plan is m/2ᵏ: 0, 0.25, 0.5, 1.

**Q24. Why report P(optimal plan) and not just energy?**
The database needs the plan. Energy is only a proxy; it is a reliable proxy near the ground state in the new formulation, and an unreliable one in the old formulation.

**Q25. What would you do next?**
Longer restarts on large problems (P4, P5); measure which plans P5 run 2 favours; test on noisy hardware models; test whether a cheaper restart count (e.g. 3) keeps most of the gain.

---

## 18. Slide plan

| # | Title | Figure | Say this |
|---|---|---|---|
| 1 | The problem: everyone is stuck at SPIQ | 01 | From SPIQ, all 10 optimizers end where they started on P1–P5. |
| 2 | Why: a flat spot | 02 | At all 17 SPIQ points the gradient is zero and no single angle goes downhill. |
| 3 | The landscape around SPIQ | 03 | A narrow dip; by distance 1 the energy is on a plateau, and no straight path leads to the plan. |
| 4 | The idea: kick, and what a kick does | 04 | From our future-work list; a kick alone always raises the energy, so the optimizer must do the work. |
| 5 | Why NFT | 06 + sine-wave formula | Five optimizers kicked; only NFT ended below SPIQ; it solves each angle exactly. |
| 6 | Choosing the kick size | 07 | Exploratory: best σ depends on the start; not claimed. |
| 7 | Round 3: one kick needs NFT, and is a gamble | 08 | Kicked COBYLA 0/120 above 0.9; kicked NFT 42/120; P4 no gain from one kick. |
| 8 | Round 4: six pipelines, one failure | 09 + 10 | Kick 1.0 ×5 wins everywhere; pre-registered halving failed on P4, for a clear reason. |
| 9 | The pipeline | 05 | Measure SPIQ, 5 × (kick 1.0 + NFT 3,000), keep the lowest measured energy. |
| 10 | Round 5: the main result | 11 | Unseen seeds: P2 0.977 / 0.973, P4 0.693 vs COBYLA 0.50 / 0.00 / 0.25. |
| 11 | 12 new problems | 12 (+13, 14 as backup) | 0.365 → 0.904, 9 of 12 solved; ahead even at equal budget. |
| 12 | Control: what really drives it | 15 + ingredients table | Random restarts do as well: it is multi-start NFT. |
| 13 | Scale: P5 (interim) | 16 | 0.00 → 0.976 and 0.612 on 2 of 4 runs; run 2 shows the energy-vs-plan limit. |
| 14 | Formulation matters | 17 | Old formulation 0.628; its minimum is the wrong plan on the cliques. |
| 15 | Conclusions, cost and limits | tables in 13–15 | Use multi-start NFT on the new formulation; costs 5–12×; P4 partly solved; noiseless. |

---

## 19. Test yourself

1. From SPIQ, what P(optimal plan) does COBYLA reach on P2 single-start, and why?
2. What is the measured gradient size at SPIQ points, and how many angles go downhill?
3. Which five optimizers were kicked in round 1, and which one ended below SPIQ's energy?
4. Why does NFT need only 2 evaluations per step?
5. What is a seed, and why are round 5's runs numbered 101–120?
6. Why did the halving pipeline fail?
7. What is the round 5 result on P4, and how many of 20 runs succeeded?
8. What did the control experiment change, and what did it show?
9. Why can the old formulation not be fixed by a better optimizer?
10. In P5 run 2, why did P(optimal plan) fall when the energy went down?

**Answers.**
1. 0.00. It sees a zero slope at SPIQ's flat spot and never moves.
2. 10⁻¹⁵ to 10⁻¹⁴ (zero); 0 downhill angles at all 17 points.
3. NFT, COBYLA, POWELL, ADAM, BOBYQA; only NFT.
4. The energy along one angle is a sine wave with 3 unknowns; the current value is known, so 2 more (at ±π/2) give the whole curve and its exact minimum.
5. The starting number of the random generator; same seed = same run. 101–120 were never used before, so they played no part in choosing the method.
6. Its 400-evaluation screen favoured small kicks, which look good early by sliding back to SPIQ but stay stuck; on P4 it picked kick 0.2 in 19 of 20 runs.
7. 0.693 (vs COBYLA 0.250); 11 of 20 runs above 0.9.
8. Only where restarts start: uniformly random instead of SPIQ + kick. Random restarts did as well, so the gain is from NFT + several restarts (multi-start NFT).
9. On some problems its lowest-energy state is the wrong plan, so lowering the energy moves away from the optimal plan.
10. It is far from the ground state; there a lower average energy can come from favouring 2nd- and 3rd-best plans. The pipeline can only see energy.
