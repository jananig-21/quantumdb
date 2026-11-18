# Optimizers under the new "qaoa" formulation: pilot

> **Presentation material (18 Nov 2026):** `PRESENTATION_kick_NFT/KICK_NFT_COMPLETE_GUIDE.md` explains every concept, run, file and figure, and adds the 12 new problems, the random-restart control and P5. Figures are in `PRESENTATION_kick_NFT/figures/`, numbered in talk order.

**Current result: rounds 4-5 below.** Five restarts of kick 1.0 + NFT, keeping the lowest measured energy, confirmed on unseen seeds:

| Problem, start | COBYLA P(opt) | Kick + NFT P(opt) | Runs over 0.9 (95% interval) | COBYLA energy | Kick + NFT energy | Ground | p (P(opt)) |
|---|---|---|---|---|---|---|---|
| P2 multi-start point | 0.500 | **0.977** | 19/20 (0.75-1.00) | 3.742 | 3.660 | 3.653 | 3.4e-8 |
| P2 single-start point | 0.000 | **0.973** | 19/20 (0.75-1.00) | 5.127 | 3.661 | 3.653 | 3.4e-8 |
| P4 single-start point | 0.250 | **0.693** | 11/20 (0.32-0.77) | 7.425 | 5.190 | 4.908 | 6.1e-4 |
| P1, P3 | 1.000 | 1.000 | 20/20 | exact | exact | | no change |

Question: can the choice of optimizer improve energy or join plans once the formula is fixed?

Starts are the team's saved SPIQ points (`formulation` 79daf329). The SPIQ search was not run. Energy is written from `docs/problem-formulation.md` and checked against the saved ground states and SPIQ energies. Noiseless simulation, 10,240 shots per evaluation, reps 2, scored exactly at the best-seen point. `+k` = Gaussian kick of that sigma (radians) on every angle before optimizing.

- Round 1 (`pilot_round1/`): 10 optimizers, 5 trials, P1-P3. Budget 1,000.
- Round 2 (`pilot_round2/`): kick sizes on P2 and P4. Used to choose settings, so its numbers are not the claim.
- **Round 3 (this folder, the claim):** new seeds, 20 trials per cell, allowance of 3,000 evaluations for every method. COBYLA stops by itself when it converges; a kicked COBYLA used 1,640-3,000.

## Round 3: mean fallback P(opt), and mean exact energy

| Problem, start | COBYLA | NFT | NFT + k0.5 | NFT + k1.0 | Ground |
|---|---|---|---|---|---|
| P2 multi-start point: P(opt) | 0.500 | 0.501 | **0.950** (18/20 over 0.9) | 0.531 (8/20) | |
| P2 multi-start point: energy | 3.742 | 3.743 | **3.665** | 4.058 | 3.653 |
| P2 single-start point: P(opt) | 0.000 | 0.000 | 0.142 (0/20) | **0.653** (12/20) | |
| P2 single-start point: energy | 5.127 | 4.396 | 3.909 | 3.921 | 3.653 |
| P4 single-start point: P(opt) | 0.250 | 0.296 | 0.229 (2/20) | 0.276 (2/20) | |
| P4 single-start point: energy | 7.426 | 6.969 | **6.014** | 6.640 | 4.908 |

One-sided Mann-Whitney on P(opt), kicked NFT vs COBYLA from the same start: P2 multi k0.5 p = 3.4e-8; P2 single k1.0 p = 5.8e-5; P4 k0.5 p = 0.98, P4 k1.0 p = 0.86 (no gain).

Random start `[0.5]*4`, energy: P2 COBYLA 9.373 vs NFT 7.645; P4 COBYLA 21.662 vs NFT 18.754. P(opt): P2 0.093 vs 0.136; P4 0.166 vs 0.140. COBYLA stops at about 115 evaluations.

Rank correlation of energy and P(opt), kicked NFT runs: P2 -0.906, P4 -0.300.

P1 and P3: the SPIQ start is already the exact optimum, so there is nothing to gain. In round 1 (5 trials) kicks did not improve them and larger kicks can hurt.

## Did it replicate?

| Setting | Round 2 | Round 3 (new seeds) |
|---|---|---|
| P2 multi k0.5, NFT | 0.974 (18/20) | **0.950 (18/20)** replicated |
| P4 k1.0, NFT | 0.574 (5/10) | **0.276 (2/20)** did not replicate |

## Limits

- Two problems only, and one saved SPIQ point per start. P5 and larger query sets are untested.
- The best kick size differs by start (0.5 vs 1.0), and no rule for choosing it in advance exists yet.
- `analyze_round3.py` estimates best-of-N restarts by resampling the 20 runs. That is an estimate, not N real runs per pipeline, and it costs N x 3,000 evaluations.
- BOBYQA at this budget stalled for 30 minutes and was dropped. At 1,000 evaluations (round 1) it stayed stuck from the SPIQ start and kicks made it worse.
- Noiseless simulation; 3,000 x 10,240 shots per run is expensive on real hardware.

## Round 4: real restarts, 6 pipelines, 20 runs each, 15,000-evaluation allowance

Plan fixed in advance (`restarts.py`). Each pipeline picks by measured energy only and keeps the unkicked start as a candidate. Mean fallback P(opt), with runs over 0.9 in brackets:

| Problem, start | COBYLA | NFT | Kick 0.5 x5 | Kick 1.0 x5 | Mixed kicks x5 | Halving (primary) |
|---|---|---|---|---|---|---|
| P2 multi | 0.500 (0) | 0.904 (16) | 0.999 (20) | 0.998 (20) | 0.999 (20) | 0.998 (20) |
| P2 single | 0.000 (0) | 0.050 (1) | 0.533 (9) | 0.949 (19) | 0.768 (13) | 0.786 (15) |
| P4 single | 0.250 (0) | 0.398 (1) | 0.337 (3) | 0.761 (13) | 0.603 (9) | **0.196 (1)** |
| P1, P3 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

**The primary pipeline, halving, failed on P4.** It picked kick 0.2 in 19 of 20 runs: its 400-evaluation screen selected the smallest kicks, which then ended at worse energies than kick 1.0 reached. So it is not the kick-size rule. Kick 1.0 x5 was chosen after seeing this table, hence round 5.

## Round 5: kick 1.0 x5 vs COBYLA on unseen seeds (the claim)

Reps 101-120, table at the top of this file. Round 4 to round 5 on P4: 0.761 (13/20) to 0.693 (11/20).

On P4, 9 of 20 runs still end below 0.9, and all 11 successes have energy at or below 4.943. Other runs end on mixed or wrong plans (5.095-6.219).

## Limits of rounds 4-5

- Two problems show a gain, each from saved SPIQ points. P5 and larger query sets are untested. Kick 1.0 is confirmed only here, and is not shown to be a general rule.
- Cost: the pipeline uses 15,001 evaluations. COBYLA stops by itself after about 1,240 on P2 and 2,900 on P4, so the gain costs about 12x and 5x the evaluations.
- Noiseless simulation; 15,001 x 10,240 shots is about 154 million shots per pipeline run.
- Divya's run-time figures were not re-measured here.
