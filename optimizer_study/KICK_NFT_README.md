# Kick + NFT: escaping the SPIQ flat spot

This branch adds the kick + NFT optimizer study, built on `spiq-joo`. Every file is under `optimizer_study/`, with the same paths the code was run from, so the scripts' relative imports work.

**Start here:** [`new_formulation_pilot/PRESENTATION_kick_NFT/KICK_NFT_COMPLETE_GUIDE.md`](new_formulation_pilot/PRESENTATION_kick_NFT/KICK_NFT_COMPLETE_GUIDE.md). It explains every concept, run, file and figure, with the evidence behind each number.

## Main results

| Problem, start | COBYLA from SPIQ | Kick 1.0 ×5 + NFT | Random ×5 + NFT (control) |
|---|---|---|---|
| P2, multi-start SPIQ point | 0.500 | **0.977** | 0.999 |
| P2, single-start SPIQ point | 0.000 | **0.973** | 0.975 |
| P4, single-start SPIQ point | 0.250 | **0.693** | 0.712 |
| 12 new problems (mean) | 0.365 | **0.904** | 0.963 |

The values are the mean probability of the optimal join plan:
- **P2 and P4:** 20 runs each (round 5), on seeds never used to choose the method.
- **12 new problems:** 3 runs each.

All results use noiseless simulation, QAOA depth 2, and 10,240 shots per evaluation.

**P5 (24 qubits) is interim.** The run is in progress, and `log_p5_kick.txt` shows the state after 4 of 5 kicks.

## Contents

| Path | What it is |
|---|---|
| `new_formulation_pilot/PRESENTATION_kick_NFT/` | the complete guide, 17 figures in presentation order, the data behind them, and the figure scripts |
| `new_formulation_pilot/pilot.py` | runs optimizers from the team's saved SPIQ points on P1–P5 (new "qaoa" formulation), with optional kicks |
| `new_formulation_pilot/restarts.py` | round 4 (six pipelines, same 15,000-evaluation budget) and round 5 (`CONFIRM=1`, unseen seeds) |
| `new_formulation_pilot/p5_kick.py` | kick 1.0 ×5 + NFT on P5 |
| `new_formulation_pilot/pilot_round1/` | round 1: all 10 optimizers from SPIQ, then 5 optimizers kicked at σ 0.05 / 0.2 / 0.5 (P1–P3) |
| `new_formulation_pilot/pilot_round2/` | round 2: kick-size sweep, 1,000 and 3,000 evaluations (P2, P4) |
| `new_formulation_pilot/pilot_round3/` | round 3: one kick + 3,000 evaluations, 20 new seeds |
| `new_formulation_pilot/runs_round4.jsonl`, `runs_round5_confirm.jsonl` | round 4 and round 5 results |
| `new_formulation_pilot/analyze_round3.py`, `analyze_round4.py`, `analyze_round5.py`, `summ.py` | the analysis behind the numbers above |
| `new_formulation_pilot/launch_*.sh`, `log_*.txt` | the exact commands and run logs |
| `new_formulation_pilot/REPORT.md` | the new-formulation pilot report |
| `new_problems/` | the 12 new problems (chain, star, cycle, clique, 4 tables): how they were picked, their SPIQ points in both formulations, and kick + NFT runs in both formulations (`runs_new.jsonl`, 360 runs) |
| `MEETING_2026-11-18/Part2_Kick_plus_NFT/` | the kick + NFT meeting write-up and its figures |
| `MEETING_2026-11-18/scripts/` | `landscape.py` (slopes at SPIQ), `ablation_random.py` (random-restart control) and `part2_figures.py` |
| `MEETING_2026-11-18/data/` | landscape measurements, the control experiment's runs and logs, and the figure numbers |
| `comparison_2x2/analyze_2x2.py` | the formulation × start analysis, including the kick columns |
| `optimizers.py`, `fastsim.py`, `many_problems/{fastqaoa,relaxed_sim,run_many,spiq_many,fast_clapton,SteinbrunnQueryGenerator}.py`, `new_code_study/study.py` | the shared simulation and optimizer code the scripts above import |

## Running

The scripts import the team's code from `base/` and take the path of a checkout of the `formulation` branch as their first argument (see the `launch_*.sh` files).

Seeds are derived from each run's label with `zlib.crc32`, so any run can be reproduced exactly.
