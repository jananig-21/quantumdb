# How the join-order problem is turned into a QAOA problem

This file explains, in plain terms, how we turn "pick the best join order" into something QAOA can solve.
It covers what we had before, what was wrong with it, and what we changed.
Keep it updated whenever the formulation changes or an experiment tests it. Code changes also go in `CHANGELOG.md`.

---

## 1. The idea in one paragraph

QAOA cannot "search for the cheapest join order" directly. It can only find bitstrings with **low energy**.
So we write an energy formula over bits with two goals:
- every bitstring that is **not a real join order** gets a **penalty**, so its energy is high;
- among real join orders, **cheaper plans get lower energy**.

If both goals hold, the lowest-energy bitstring is the best join order, and QAOA, by looking for low energy, finds good plans.
If either goal fails, QAOA can be doing its job perfectly and still return a bad plan.

**Formula vs Hamiltonian.** The energy formula is called a QUBO. The code (`qubo.to_ising()`) turns it into a Hamiltonian, which is the same formula rewritten for qubits. That conversion is exact and was never the problem. Everything below is about the formula itself.

## 2. Our test problems

| | Tables (rows in each) | Join conditions (selectivity) | Best plan's cost |
|---|---|---|---|
| P1 | 10, 15, 20 | 0–1 (0.1), 1–2 (0.1) | 15 |
| P2 | 10, 15, 20, 30 | 0–1 (0.1), 2–3 (0.1) | 315 |
| P3 | 10, 15, 20, 30 | 0–1 (0.1), 1–2 (0.2), 2–3 (0.3), 0–2 (0.4), 0–3 (0.5), 1–2 again (0.6) | 29.4 |
| P4 (5 tables, 15 qubits) | 10, 15, 20, 30, 50 | 0–1 (0.1), 1–2 (0.2), 2–3 (0.1), 3–4 (0.3), 1–3 (0.5) | 165 |
| P5 (6 tables, 24 qubits) | 10, 15, 20, 30, 50, 100 | 0–1 (0.1), 1–2 (0.2), 2–3 (0.1), 3–4 (0.3), 4–5 (0.1), 1–4 (0.5), 0–5 (0.4) | 1605 |

P4 and P5 (input_idx 3 and 4) only fit with the new formulation: it needs R × (R − 2) qubits, 15 and 24. The old one would need (5 + 5) × 3 = 30 and (6 + 7) × 4 = 52.

A plan's **cost** is the total size of the intermediate results. Example for P2 with order 0, 1, 2, 3:
- joining 0 and 1 gives 10 × 15 × 0.1 = **15** rows;
- adding table 2 gives 15 × 20 = **300** rows;
- so the cost is 15 + 300 = **315**. The final join always contains every table, so it is the same for every plan and isn't counted.

## 3. What we had before (the "Q2O" formulation, built for digital annealers)

### The bits

**What a join condition is.** A join condition is a predicate in the query's WHERE clause that links two tables, such as `T1.id = T2.id`. Its **selectivity** is the fraction of row pairs that survive it. The input files list the conditions in `pred.txt` (which two tables) and `pred_sel.txt` (the selectivity).

**Why P3 has 6.** P3 was defined with six entries: 0–1, 1–2, 2–3, 0–2, 0–3, and 1–2 again. Four tables have 6 possible pairs, and P3 covers 5 of them (1–3 is missing). The pair 1–2 is listed twice (0.2 and 0.6). That means two separate conditions between the same two tables, for example `T1.a = T2.a AND T1.b = T2.b`. Together they filter by 0.2 × 0.6 = 0.12.

**Which join steps get bits.** With 4 tables there are 3 joins. After the 1st join the result holds 2 tables, after the 2nd it holds 3, and after the 3rd it holds all 4. The last result is the same for every plan, so only the first **2 join steps** (R − 2 in general) need bits.

- One bit per **(table, join step)**: "is this table in the result after this step?" With 4 tables and 2 steps that is 4 × 2 = **8 bits**.
- One extra bit per **(join condition, join step)**: "is this condition's selectivity applied to this step's result?" With P3's 6 conditions and 2 steps that is 6 × 2 = **12 bits**.

That makes (4 + 6) × 2 = **20 qubits** for P3. P1 has (3 + 2) × 1 = 5 and P2 has (4 + 2) × 2 = 12.

**These extra bits are redundant.** A condition applies exactly when both of its tables are in the result, so the table bits already decide it. For the order T0, T1, T2, T3:
- after step 1 the result is {T0, T1}, so only 0–1 should be on;
- after step 2 the result is {T0, T1, T2}, so 0–1, 1–2, 0–2 and the second 1–2 should be on.

In the old formulation those 12 bits were separate variables. Penalty C (below) only stopped a condition being on when one of its tables was missing. Nothing forced a condition on when both tables were present. As §4 shows, the optimizer could exploit this gap. The new formulation (§5) removes these bits: it writes "condition applies" as the product of the two table bits, which is 1 only when both tables are in. That takes P3 from 20 qubits to 8.

### The energy
- **Penalty A:** each join step must contain the right number of tables (2 after the first join, 3 after the second, and so on).
- **Penalty B:** once a table is in, it stays in.
- **Penalty C:** a join condition may only be "on" if both of its tables are in.
- **Cost:** for each join step, take the log of the result size and **square it**. Logs turn the multiplications into additions.
- All penalties are multiplied by a weight λ, which was set to **the number of join steps** (1 for P1, 2 for P2 and P3).

**The old Hamiltonian, in full.** Notation:
- v_{t,j} = "table t is in the result after step j", and w_{p,j} = "condition p is applied at step j";
- c_t = log10 of table t's size;
- s_p = log10 of condition p's selectivity, which links tables a_p and b_p;
- J = R − 2 steps, and λ = J.

$$H_{\text{q2o}} = \lambda\Big[\underbrace{\sum_j\big(j+2-\textstyle\sum_t v_{t,j}\big)^2}_{\text{A}} + \underbrace{\sum_t\sum_{j\ge1}\big(v_{t,j-1}-v_{t,j-1}v_{t,j}\big)}_{\text{B}} + \underbrace{\sum_{p,j}\big(2w_{p,j}-w_{p,j}v_{a_p,j}-w_{p,j}v_{b_p,j}\big)}_{\text{C}}\Big] + \underbrace{\sum_j\Big(\sum_t c_t\,v_{t,j}+\sum_p s_p\,w_{p,j}\Big)^2}_{\text{cost}}$$

Its Pauli (Z/ZZ) form and every coefficient are in §9.0.

## 4. What was wrong with it

### Problem 1: the penalty was too weak, so the lowest energy belonged to an invalid bitstring
Real-plan energies were ranked correctly; the cheapest plan had the lowest energy among real plans. But a **fake plan** that breaks a rule could get an even lower energy, because breaking a rule cost only λ = 2 while saving more than that.

**P2 example** (all logs are base 10):

| Bitstring | What it means | Energy |
|---|---|---|
| Best real plan | step 1 = {0,1}, step 2 = {0,1,2} | 1.18² + 2.48² = **7.52** |
| Fake plan | step 1 = {0,1}, step 2 = {0,1}: table 2 forgotten, so Penalty A = 1 | 2 × 1 + 1.18² + 1.18² = **4.77** |

The fake plan wins. When decoded, it turns into the order 1, 0, 3, 2, which costs **465**. That is exactly what we saw: 99.99% of P2 samples had cost 465.
For the real plan to win, λ needed to be above 4.75. It was 2.

**P1 example.** The lowest energies belonged to bitstrings that switched on a join condition whose tables weren't both present, which Penalty C allows at the price of only 1. One of them decodes to the cost-**30** plan. That explains P1's 50/50 split between 15 and 30.

### Problem 2: squaring the cost rewards "making the result look small"
With a squared cost, anything that pushes the log-size toward 0 lowers the energy, including cheating: dropping a table, or switching on a join condition that doesn't apply. And the bigger the join, the more one cheat saves:
- squared cost: dropping a table of log-size ℓ from a join of log-size L saves 2·L·ℓ − ℓ². In P2 that is **4.75**;
- plain (unsquared) cost: it saves only ℓ. In P2 that is **1.30**.

For a 10-table query with tables of about 100,000 rows, the squared version would need a penalty around **175**, compared with about **5** for the plain one.

### Problem 3: too many qubits
The join-condition bits add (number of conditions) × (number of join steps) qubits: 2, 4 and 12 extra for P1, P2 and P3. They exist only because the squared cost can't be written without them. For larger queries this grows quickly: see "How the qubit count grows with problem size" in §5.

### Problem 4: the cost checker disagreed with the energy on P3
P3 lists the pair 1–2 twice. That means two separate conditions between the same tables, for example `A.x = B.x AND A.y = B.y`, so both selectivities should apply: 0.2 × 0.6. The energy formula applied both, but the cost checker (`Postprocessing.py`) applied only the first. So the checker said P3's best cost was 39, while the true best is 29.4. It was the same plan either way.

### What was fine
- Penalties A and B together are exactly right. A bitstring passes both if and only if it is a real join order. Nothing is missing and nothing is extra.
- Turning the formula into the Hamiltonian, and the QAOA circuit built from it, are correct. They faithfully encode the formula; the formula was the problem.

## 5. What we changed (the "qaoa" formulation)

### The intuition, in plain terms

Think of the energy as a **price tag** on every bitstring, and of QAOA as a shopper who mostly picks cheap items. The price has two parts: **fines** for breaking the rules of a real plan (λ × penalties), and the **plan's own cost**. Each change below makes the price tag more honest.

**1. "Does this join condition apply?" is no longer a separate choice.**
- *In database terms:* a query optimizer never *decides* whether to apply `WHERE T1.id = T2.id`. The predicate applies automatically as soon as T1 and T2 are both in the intermediate result. Planners derive the applicable predicates from the set of joined relations.
- *The old formulation* gave every (condition, step) its own on/off switch, and a rule (Penalty C) to catch switches that were on illegally. The switch could still be on when it shouldn't be, or off when it should be on. Either way the optimizer could make a result *look* smaller than it really is.
- *The new formulation* does what a planner does: "condition applies" = "table a is in" AND "table b is in", the product of two bits. There is nothing left to cheat with, no Penalty C, and many fewer qubits.

**2. The cost is the log of the result size, not its square.**
- *In database terms:* cardinality estimates multiply (rows × rows × selectivity), so taking logs turns the estimate into a sum. Each table in the result adds its log row count, and each applicable predicate subtracts the log of its selectivity. This is the usual log-space trick in cost models.
- *Why the square was there:* the real cost (total rows produced) is dominated by the biggest intermediate result, and squaring makes big results count more.
- *Why it backfired:* squaring rewards *anything* that pushes the estimated size toward 1 row. That includes impossible plans that forget a table, or switch on a predicate that doesn't apply. It is like a cost model that loves plans with an estimated 0 rows, even when the estimate is wrong because the plan is broken. With the plain log, forgetting a table saves exactly that table's log row count, no more. That makes the fine needed to stop it small and easy to compute.

**3. The fine (λ) is set from the data.**
- *In everyday terms:* a parking fine must be higher than the parking fee, or people stop paying. But a fine of a million dollars makes every other price look like zero.
- *Here:* λ must exceed what one cheat can save (one table's log size, change 2), so we set λ = 1.1 × the largest log table size. If λ is too large, the price differences between real plans become a tiny fraction of the total price range, and QAOA, which only sees the overall landscape, can't tell the best plan from the second best.

**What it costs us.** The energy now ranks real plans by the **product** of their intermediate sizes, the sum of their logs. The real cost is the **sum** of the sizes. In database terms, we rank plans by the geometric mean of the work instead of the total work. The best plan comes out on top in every problem we checked. Among mediocre plans the order disagrees a little more often than with the square (§9).

| | Before ("q2o") | Now ("qaoa") |
|---|---|---|
| Bits | tables × steps **+ conditions × steps** | tables × steps only |
| Qubits (P1 / P2 / P3) | 5 / 12 / 20 | **3 / 8 / 8** |
| How a join condition is counted | its own bit, kept honest by Penalty C | "table a is in **and** table b is in", written as the product of their two bits, so it **cannot cheat** |
| Cost per join step | (log size)² | **log size** (not squared) |
| Penalty weight λ | number of join steps (1, 2, 2) | **1.1 × the largest log table size** (1.43, 1.63, 1.63), computed from the data |
| Penalties A and B | yes | yes, unchanged |
| Penalty C | yes | **gone**: there are no condition bits left to police |
| Repeated pair (P3 1–2) | energy: both applied; checker: first only | both applied in the energy **and** in the checker |

### The new Hamiltonian, in full

The notation is the same as in §3, but there are no w bits, and λ = 1.1 · max_t c_t.

$$H_{\text{qaoa}} = \lambda\Big[\underbrace{\sum_j\big(j+2-\textstyle\sum_t v_{t,j}\big)^2}_{\text{A}} + \underbrace{\sum_t\sum_{j\ge1}\big(v_{t,j-1}-v_{t,j-1}v_{t,j}\big)}_{\text{B}}\Big] + \underbrace{\sum_j\Big[\sum_t c_t\,v_{t,j} + \sum_p s_p\,v_{a_p,j}\,v_{b_p,j}\Big]}_{\text{cost}}$$

Compared with §3:
- Penalty C is gone.
- w_{p,j} is replaced by v_{a_p,j} v_{b_p,j}.
- The cost is no longer squared.
- λ comes from the data.

Its Pauli (Z/ZZ) form and every coefficient are in §9.0.

### Term by term: what each piece means, and how old and new differ

Both Hamiltonians have the same shape: **λ × (rule-breaking fines) + (plan cost)**. A real join order pays no fines, so its energy is just its cost. Below, each piece is shown with what it measures, why it has that form, and what it becomes on the qubits. To get qubit terms, substitute x = (1 − Z)/2. A product of two bits then becomes a ZZ coupling between their qubits, and a single bit becomes a Z field.

**The bits.**
- **Both:** v_{t,j} = 1 means "table t is in the intermediate result after join step j". There are J = R − 2 steps: the last join always holds every table, so it needs no bits.
- **Old only:** w_{p,j} = 1 means "join condition p is applied at step j". It is a separate yes/no choice that the optimizer can set either way.
- **New:** there is no w. "Condition p applies at step j" is computed as v_{a_p,j} · v_{b_p,j}, which is 1 exactly when both of its tables are in. A planner works the same way: a predicate applies once both its relations have been joined.

**Penalty A, the right number of tables at each step: (j + 2 − Σ_t v_{t,j})².** *Identical in old and new.*
- *Meaning:* after the first join the result holds 2 tables, after the second 3, and so on. The bracket counts how far the bitstring is from that number. It is 0 when the count is right, and 1, 4, 9, … when it is 1, 2, 3 off.
- *Why squared:* it has to be zero only at the right count and positive on both sides. Too many tables and too few must both be fined.
- *On the qubits:* expanding the square couples **every pair of tables at the same step**, with ZZ coefficient 2λ/4. It also gives each table a field that depends on the step. This is the densest part of the new Hamiltonian, with J·R(R−1)/2 couplings.

**Penalty B, once in, always in: v_{t,j−1} − v_{t,j−1}·v_{t,j} = v_{t,j−1}(1 − v_{t,j}).** *Identical in old and new.*
- *Meaning:* it is 1 only if table t is in the result at step j−1 and gone at step j. A left-deep plan only ever adds tables, so a table can't leave the result.
- *On the qubits:* each table is coupled to itself at the next step (ZZ −λ/4), a chain through the steps, R(J−1) couplings in total. This is what makes the steps describe one join order, not unrelated sets.

**Penalty C, a condition may be on only if both its tables are in: w_{p,j}(2 − v_{a_p,j} − v_{b_p,j}).** *Old only.*
- *Meaning:* it is 0 if the condition is off. If the condition is on, the fine is 0, 1 or 2 when 0, 1 or 2 of its tables are missing.
- *The gap:* nothing fines a condition that is *off* while both its tables are in. Since a condition's selectivity lowers the cost, leaving it off never pays, but switching it on illegally can. C has to catch that, with only λ = J behind it. This is how the old ground state on P1 and P2 came to be a fake plan (§4).
- *On the qubits:* each condition bit is coupled to its two tables (ZZ −λ/4) and has its own field, which adds P·J qubits.
- *New:* no w bits, so C is not needed. A condition can't be on wrongly or off wrongly, because it is computed from the table bits.

**The cost: the size of each intermediate result, in logs.**
- *What it measures:* the estimated size of step j's result is (product of its tables' row counts) × (product of the selectivities of the conditions that apply). In logs that becomes a sum: log size_j = Σ_t c_t v_{t,j} + Σ (s_p over applied conditions). Each table adds its log row count. Each applied condition adds its log selectivity, which is negative, so it lowers the size.
- *Old:* Σ_j (Σ_t c_t v_{t,j} + Σ_p s_p w_{p,j})², the **square** of each step's log size, with conditions counted through the free w bits.
- *New:* Σ_j (Σ_t c_t v_{t,j} + Σ_p s_p v_{a_p,j} v_{b_p,j}), the log size itself, with conditions counted through the table bits.
- *Why the old one squared:* the real cost is the **sum** of the sizes, which the biggest result dominates. Squaring the log makes big results count more, which is closer to that.
- *What the square does on the qubits:* (a + b + …)² produces a cross term for **every pair** of bits in the step. That couples each table with each other table (c_t c_u / 2), each table with each condition (c_t s_p / 2, even one that shares no table with it), and each condition with each other condition (s_p s_q / 2). That is why the old Hamiltonian has 2–6× more ZZ terms (3 vs 10 up to 78 vs 330), and why the old fields all depend on the size of the *whole* join (§9.0).
- *What the square does to the ranking:* it rewards pulling a step's log size toward 0, from either side. "Forget a table" and "switch on a condition that doesn't apply" both shrink it, so both are rewarded (§4, Problem 2).
- *What the linear cost does on the qubits:* each table gets a field equal to its own log size. Each condition adds s_p/4 to the coupling of **its own two tables only**, a coupling that Penalty A already creates. So conditions add no new qubits and no new couplings, only weight.
- *The trade-off:* the linear cost ranks plans by the **sum of the logs** of the sizes, which is the log of their **product**, not by the sum of the sizes. The best plan is still ranked first on P1–P5. A few mediocre plans swap places (§9.4).

**The penalty weight λ: how big the fines are compared with the cost.**
- *Old:* λ = J, the number of steps (1, 2, 2, 3, 4 for P1–P5), whatever the table sizes are. Nothing ties it to how much a cheat can save, and on P1 and P2 it was too small (§4).
- *New:* λ = 1.1 × the largest c_t. With a linear cost, leaving a table out lowers the cost by at most that table's c_t (by less when it also drops a condition, since that removes a negative s_p). So a fine slightly above the largest c_t makes that cheat a loss. It also stays small, which keeps the gaps between real plans a visible fraction of the energy range (§9.4).

**The constant (offset).** In both, λ Σ_j (j+2)², plus whatever the conversion to Z adds. It shifts every energy by the same amount, so it changes nothing about which bitstring wins. It's only kept so that QUBO and Hamiltonian energies can be compared directly.

**How each one is built in the code (`base/Scripts/QUBOGenerator.py`).** Both write the formula with docplex, convert it with `from_docplex_mp` to a qiskit `QuadraticProgram` (the QUBO), and QAOA turns that into the Hamiltonian with `qubo.to_ising()`. Both take log10 of the sizes and selectivities, without rounding. `qubo-formulation` in `base/config.py` (or the environment variable `QUBO_FORMULATION`) picks the formulation.

| Step | Old: `generate_IBMQ_QUBO_for_left_deep_trees_v2` | New: `generate_qaoa_qubo_left_deep` |
|---|---|---|
| Variables | `v` (tables × steps) **and** `pred_vars` (conditions × steps) | `v` only |
| Penalties | `H_A + H_B + H_pred` (C = `H_pred_a + H_pred_b`, one per table of the condition) | `H_A + H_B`, written the same way |
| Cost | a loop over steps that adds `(Σ card·v + Σ sel·pred_vars)**2` | one sum of `card·v + sel·v[a]·v[b]` over steps, not squared |
| λ | a counter that goes up by 1 in the cost loop, giving J | `1.1 * max(card)`, computed before the model is built |
| Repeated condition (P3 lists tables 1–2 twice) | two separate w bits | both selectivities land on the same v_1·v_2 coupling (they add up) |

### How the qubit count grows with problem size

With R tables and P join conditions there are R − 2 join steps that need bits:
- **New:** R × (R − 2) qubits. This depends only on the number of tables, and grows roughly as R².
- **Old:** (R + P) × (R − 2) qubits. It grows with the number of conditions too. A query has at least R − 1 conditions (a chain) and at most R(R − 1)/2 (every pair joined), so the old count grows roughly as 2R² to R³/2.

| Tables (R) | New | Old, chain (P = R − 1) | Old, every pair joined (P = R(R − 1)/2) |
|---|---|---|---|
| 3 | 3 | 5 | 6 |
| 4 | 8 | 14 | 20 |
| 5 | 15 | 27 | 45 |
| 6 | 24 | 44 | 84 |
| 7 | 35 | 65 | 140 |
| 8 | 48 | 90 | 216 |
| 10 | 80 | 152 | 440 |

Our problems: P1 3 vs 5, P2 8 vs 12, P3 8 vs 20, P4 15 vs 30, and P5 24 vs 52 (new vs old).

**What fits in a simulation.** An exact simulation stores 2^n numbers, at 16 bytes each: 256 MB at 24 qubits, 512 MB at 25, 16 GB at 30, and 16 TB at 40. With our 25-qubit limit (the container has about 7.8 GB of RAM), the new formulation reaches **6 tables** (24 qubits). The old one reaches **4 tables**, and only with up to 8 conditions. Every 5-table query is at least 27 qubits in the old formulation.

**The circuit grows too.** SPIQ gives every gate its own angle, so the number of angles follows the number of terms in the energy. Measured: P4 has 140 angles (new) vs 410 (old), and P5 has 252 vs 868. The SPIQ search took 3–4× longer with the old formulation (see `docs/RESULTS_INDEX.md`, "Findings").

### Why each change helps
- **Removing the condition bits** saves qubits, makes illegal conditions impossible, and removes a whole penalty.
- **Not squaring** removes the reward for cheating. It still ranks the plans just as well; in all three problems the cheapest plan stays cheapest. With a single join step, the energy of a real plan is exactly log(cost): P1's plans get 1.18 = log 15, 1.48 = log 30 and 2.30 = log 200.
- **λ from the data.** Cheating by dropping a table saves at most that table's log size, so λ only has to beat the largest one. 1.1 × the largest is just enough. Keeping λ small matters for QAOA: a big penalty stretches the energy range, and the small difference between the best and second-best plan becomes hard to tell apart.

### Checked (exact energies of every bitstring, no QAOA run)

| | λ | Lowest-energy bitstring | Its decoded cost | Best cost |
|---|---|---|---|---|
| P1 | 1.43 | a real plan | 15 | 15 |
| P2 | 1.63 | a real plan | 315 | 315 |
| P3 | 1.63 | a real plan | 29 (29.4, stored as a whole number) | 29.4 |

### Expected effect on QAOA (estimates, not yet measured)

| | Two-qubit gates per layer (before → now) | Share of the starting state that is a real plan (before → now) |
|---|---|---|
| P1 | 10 → 3 | 15.6% → 37.5% |
| P2 | 34 → 16 | 0.78% → 4.7% |
| P3 | 94 → 16 | 0.024% → 4.7% |

### Correctness

The same style of argument the original papers use: check each part of the energy separately, then choose λ so that breaking a rule never pays. Notation:
- R tables, so there are J = R − 2 join steps with bits, numbered j = 0 … J − 1.
- Bit v[t, j] = 1 means "table t is in the result after step j". S_j is the set of tables with v[t, j] = 1.
- c_t = log(rows of t), which is ≥ 0. s_p = log(selectivity of condition p), which is ≤ 0. Condition p joins tables a_p and b_p.
- The energy is E = λ · (A + B) + C, where:
  - A = Σ_j (j + 2 − |S_j|)²
  - B = Σ_t Σ_{j ≥ 1} v[t, j−1] · (1 − v[t, j])
  - C = Σ_j ( Σ_t c_t · v[t, j] + Σ_p s_p · v[a_p, j] · v[b_p, j] )

**Claim 1: A + B = 0 exactly on the valid plans, and A + B ≥ 1 on everything else.** *Proven.*
- A is a sum of squares of whole numbers. It is 0 exactly when every step j holds j + 2 tables.
- Each term of B is 0 or 1. It is 1 exactly when a table is in at step j − 1 and out at step j. So B = 0 means S_0 ⊆ S_1 ⊆ … ⊆ S_{J−1}.
- Both 0 together means S_0 is a pair, each later step adds exactly one table, and one table is left over for the final join. That is a left-deep join order (t1, t2, …, tR) with S_j = {t1, …, t_{j+2}}. Conversely, every join order gives such bits.
- Two orders give the same bits only when they differ by swapping the first two tables. That doesn't change the cost, because joining 1 with 2 is the same as joining 2 with 1. So there are R!/2 valid bitstrings.
- A and B are whole numbers ≥ 0, so a bitstring that breaks any rule has A + B ≥ 1.

**Claim 2: on a valid plan, C is exactly the sum of the log result sizes.** *Proven.*
- For a valid bitstring, v[a_p, j] · v[b_p, j] = 1 exactly when both tables of condition p are in S_j. So step j contributes Σ_{t ∈ S_j} c_t + Σ_{p inside S_j} s_p.
- That is log( Π rows × Π selectivities ) = log |I_j|, the estimated size of the j-th intermediate result. A condition listed twice contributes both selectivities, as in the cost checker.
- So C = Σ_j log |I_j| = log Π_j |I_j|, with no approximation.
- The old formulation only reached this once Penalty C held. A condition bit could be on with a missing table, or off with both tables present (§4).

**Claim 3: with a large enough λ, no invalid bitstring has lower energy than the best plan.** *Proven for the safe λ below. Our default λ is checked per problem, not proven.*
- Let C* be the cost term of the best plan, and let C_min be the lowest C over all 2^n bitstrings.
- An invalid bitstring has E ≥ λ · 1 + C_min (Claim 1), and the best plan has E = C*. So any λ > C* − C_min is safe.
- Neither number needs a search:
  - C* ≤ C(any valid plan you write down).
  - C_min ≥ J × Σ_p s_p, because the c_t terms are ≥ 0, and at most every condition is switched on at every step.
- So this weight is always safe:

  λ_safe > C(any valid plan) − J × Σ_p s_p

  (The "Known limit" formula below is this, plus a margin of 1.)
- Our default λ = 1.1 × max c_t is smaller. It covers "drop one table", which saves at most max c_t. It does not cover "swap a table between steps" (breaking B), or several breaks at once that each save a selective condition. For example, on P1 λ_safe is about 3.2, but we use 1.43.
- So for the default λ, Claim 3 is checked instance by instance. `base/test_qubo_formulation.py` brute-forces all 2^n bitstrings and confirms that the lowest-energy one is a valid plan, for P1–P5. A new problem needs this check, or an exact classical solver, before trusting the default λ. Otherwise use λ_safe, at the price of smaller visible gaps between plans.

**Claim 4: the lowest energy belongs to the plan with the lowest true cost.** *Not provable in general, for either formulation. Checked on P1–P5.*
- By Claims 1–3, the ground state is the valid plan with the smallest Σ_j log |I_j|, which is the smallest product Π_j |I_j|. The cost we report is the sum Σ_j |I_j|.
- With a single step (R = 3, P1) the two agree, because log is increasing. With two or more steps they can disagree. Take plan X with intermediate sizes (10, 1000): product 10 000, sum 1010. Plan Y has (100, 200): product 20 000, sum 300. Our energy prefers X, while the true cost prefers Y.
- The old formulation's Σ_j (log |I_j|)² is also a stand-in, with the same kind of gap. It weights large results more, so it may track the sum more closely, but it is no more exact.
- What we know: on P1–P5 the ground state decodes to the true optimal cost (brute force over all join orders, in `test_qubo_formulation.py`). How often the two disagree on bigger or random queries is not yet measured. Comparing the product-best plan with the sum-best plan over many random queries would measure it, using only classical enumeration.

**Summary.** The bits, penalties and cost term are exact (Claims 1–2). Validity of the ground state is guaranteed with λ_safe, and checked per problem with our smaller default λ (Claim 3). Optimality for the true cost is checked, not guaranteed (Claim 4). The original formulations have the same gap: they prove the encoding and the penalty argument, and rely on experiments for the rest.

### Known limit
The λ rule is argued only for the "forgot a table" kind of cheat. It is not proven for the "swapped a table between steps" kind (see Claim 3 in "Correctness"). It holds on P1–P5; `base/test_qubo_formulation.py` checks this by brute force. If a bigger problem ever ends up with a fake plan at the lowest energy, switch to the guaranteed-safe weight:

  λ = (energy of any real plan) − (number of join steps × sum of the negative log selectivities) + 1

The safe weight is larger, so the best-vs-second-best difference gets harder to see. A comment in `QUBOGenerator.py` points to this.

## 6. What about the "Ready to Leap" formulation?
It was built for gate-based machines, but it isn't a better fit for QAOA.
- It approximates each result size with "is the size above threshold 1? above threshold 2? …" bits, plus extra helper bits to make the inequalities work. That needs many more qubits; even small problems run past our 23-qubit limit.
- Plans whose sizes fall between the same two thresholds get the same energy, so QAOA can't tell them apart.
- It still depends on a penalty weight, with the same risk as Problem 1.

(This description is from memory of the paper; re-check it before citing.)

## 7. How to use it
- `base/config.py`: `"qubo-formulation": "qaoa"` (the default) or `"q2o"` (the old one). `IBMQExperiments.py`, `spiq_initialization.py` and `postprocess_results.py` all read this setting.
- SPIQ output files now end in `_<formulation>`, e.g. `spiq_initial_point_input1_reps2_4rel_2pred_qaoa.json`, and record the formulation inside the file. `IBMQExperiments.py` refuses to load a SPIQ file made for the other formulation. Older files, which don't record one, count as `"q2o"`.
- To re-process **old** results with `postprocess_results.py`, set the config to `"q2o"` first.
- Old CSVs (Week93–95) were scored with the old checker, so P3's target in `graphs/seeded-multi-spiq.rmd` stays at **39** for them. New P3 runs score the best plan as **29**.
- The decoder needed no change; it already reads only the table bits.

## 8. To do
- [x] Run `python3 test_qubo_formulation.py` in the container (from `base/`). It passes.
- [x] Random-start QAOA with the new formulation (Week96–98). Results are in `docs/RESULTS_INDEX.md`, "Findings".
- [x] SPIQ-started QAOA with the new formulation (Week102–104, same settings as Week93–95). Results are in `docs/RESULTS_INDEX.md`, "Findings".
- [ ] P2 still splits 50/50 between the best plan (315) and the second-best (465). Ideas: more SPIQ starts, or a slightly stronger cost term so the two plans are further apart in energy.
- [ ] Later, only if needed: an exact-cost version, or special QAOA mixers that never leave the set of real plans (these need no penalty, but cost more qubits).

## 9. Verification: the Hamiltonians and the exact eigensolver, old vs new

Everything here is classical: no QAOA. `base/hamiltonian_report.py`, run in the container from `base/`, regenerates it. The full Hamiltonians of every problem are in `docs/hamiltonians/P<k>_<formulation>.txt`, with the numbers in `docs/hamiltonians/summary.json`.

**Reading a Hamiltonian.** `qubo.to_ising()` replaces each bit x by (1 − Z)/2, so the QUBO becomes a constant (the *offset*) plus terms on one qubit (`Z`) and on two qubits (`ZZ`). The labels are qiskit's: qubit 0 is the **rightmost** letter. The conversion is exact: QUBO energy = Ising energy + offset, for every bitstring.

### 9.0 The general form of both Hamiltonians

**Notation.**
- R tables and P join conditions. J = R − 2 join steps get bits, numbered j = 0 … J−1.
- c_t = log10(card of table t), and s_p = log10(selectivity of condition p), which is negative.
- Condition p links tables a_p and b_p.
- v_{t,j} = 1 means "table t is in the result after join step j".
- w_{p,j} = 1 means "condition p is applied at step j". Only the old formulation has these bits.

**The new formulation ("qaoa"): R·J qubits, λ = 1.1 · max_t c_t**

$$H = \lambda\Big[\sum_j \big(j+2-\textstyle\sum_t v_{t,j}\big)^2 + \sum_t\sum_{j\ge1}\big(v_{t,j-1}-v_{t,j-1}v_{t,j}\big)\Big] + \sum_j\Big[\sum_t c_t\,v_{t,j} + \sum_p s_p\,v_{a_p,j}\,v_{b_p,j}\Big]$$

**The old formulation ("q2o"): (R+P)·J qubits, λ = J**

$$H = \lambda\Big[H_A + H_B + \sum_{p,j}\big(2w_{p,j} - w_{p,j}v_{a_p,j} - w_{p,j}v_{b_p,j}\big)\Big] + \sum_j\Big(\sum_t c_t\,v_{t,j} + \sum_p s_p\,w_{p,j}\Big)^2$$

H_A and H_B are the same two penalty sums as in the new formulation.

**Written out bit by bit (the QUBO).** Since x² = x for a bit, each formula is a constant, plus a weight h on each bit, plus a weight K on each pair of bits. Pairs not listed have weight 0. [j < J−1] is 1 when step j is not the last one and 0 otherwise.

| Term | New ("qaoa") | Old ("q2o") |
|---|---|---|
| h on v_{t,j} | λ(1 − 2(j+2)) + λ·[j < J−1] + **c_t** | λ(1 − 2(j+2)) + λ·[j < J−1] + **c_t²** |
| h on w_{p,j} | (no such bit) | 2λ + s_p² |
| K between two tables at the same step, v_{t,j}–v_{u,j} | 2λ + **Σ s_p over conditions linking t and u** | 2λ + **2 c_t c_u** |
| K between a table at two consecutive steps, v_{t,j}–v_{t,j+1} | −λ | −λ |
| K between table and condition at the same step, v_{t,j}–w_{p,j} | (no such bit) | 2 c_t s_p − λ·[t is a_p or b_p] |
| K between two conditions at the same step, w_{p,j}–w_{q,j} | (no such bit) | 2 s_p s_q |
| Constant | λ Σ_j (j+2)² | λ Σ_j (j+2)² |

**As a Hamiltonian.** Replace every bit x_i by (1 − Z_i)/2. Then:
- each pair weight K_ik becomes a **ZZ term with coefficient K_ik / 4**;
- each bit becomes a **Z term with coefficient −h_i/2 − Σ_k K_ik/4**;
- the **offset** is constant + Σ h_i/2 + Σ K_ik/4.

**Both Hamiltonians in Pauli form.** Doing that substitution gives the following. Here
- C = Σ_t c_t and Ŝ = Σ_p s_p, so C + Ŝ is the log of the full join's size;
- S_t = the sum of s_p over the conditions that touch table t;
- s_{tu} = the sum of s_p over the conditions that link tables t and u (0 if none);
- d_t = the number of conditions that touch table t;
- δ_j = ½([j < J−1] − [j > 0]), which is +½ at the first step, −½ at the last, and 0 in between (0 when J = 1).

New ("qaoa"):

$$H = \text{offset} \;-\; \tfrac12\sum_{t,j}\Big[c_t + \tfrac12 S_t + \lambda(R-2j-4+\delta_j)\Big]Z_{t,j} \;+\; \sum_j\sum_{t<u}\tfrac{2\lambda+s_{tu}}{4}\,Z_{t,j}Z_{u,j} \;-\; \tfrac{\lambda}{4}\sum_t\sum_{j<J-1} Z_{t,j}Z_{t,j+1}$$

Old ("q2o"):

$$\begin{aligned}H = \text{offset} &- \tfrac12\sum_{t,j}\Big[c_t(C+\hat S) + \lambda\big(R-2j-4-\tfrac{d_t}{2}+\delta_j\big)\Big]Z_{t,j} \;-\; \tfrac12\sum_{p,j}\Big[\lambda + s_p(C+\hat S)\Big]W_{p,j}\\ &+ \sum_j\sum_{t<u}\tfrac{\lambda+c_tc_u}{2}\,Z_{t,j}Z_{u,j} \;-\; \tfrac{\lambda}{4}\sum_t\sum_{j<J-1} Z_{t,j}Z_{t,j+1}\\ &+ \sum_{j}\sum_{t,p}\tfrac{2c_ts_p-\lambda[t\in p]}{4}\,Z_{t,j}W_{p,j} \;+\; \sum_j\sum_{p<q}\tfrac{s_ps_q}{2}\,W_{p,j}W_{q,j}\end{aligned}$$

Z_{t,j} is the Pauli Z on the qubit of v_{t,j}, and W_{p,j} is the Pauli Z on the qubit of w_{p,j}. The offset is constant + Σ h_i/2 + Σ K_ik/4, using the values from the table above.

How to read the difference:
- **Table–table couplings.** In the new formulation every pair of tables at a step has nearly the same coupling, 2λ/4, moved only by the conditions that link those two tables. In the old one the coupling is (λ + c_t c_u)/2, so it grows with the product of the two tables' sizes. Big tables dominate.
- **Single-qubit fields.** In the old formulation every field depends on C + Ŝ, the size of the *whole* join. So each qubit's field depends on every table and every condition at once. In the new one, a table's field depends only on its own size and its own conditions.
- **Extra couplings.** Only the old formulation has table–condition and condition–condition couplings. They exist only because the cost is squared, and they are the terms that make it large: J·P·R table–condition pairs plus J·P(P−1)/2 condition pairs.

**What the tables say about size.** Count the pairs with a non-zero K:
- **New:** J·R(R−1)/2 pairs inside the steps, plus R(J−1) pairs between steps. **The number of conditions P doesn't appear**: a condition only changes the weight of a table pair that is already there.
- **Old:** J·(R+P)(R+P−1)/2 + R(J−1) pairs. Squaring the cost connects **every** bit of a step to every other bit of that step, including condition bits that share no table.

That gives exactly the ZZ counts in 9.2 (3/10, 16/34, 16/94, 40/145, 78/330). It also shows why P3, with 6 conditions, has the same size as P2 in the new formulation but is 3× larger in the old one.

**Checked:** `base/test_hamiltonian_form.py` builds every coefficient from these formulas and compares it with `qubo.to_ising()`, for P1–P5 and both formulations. It also checks the closed-form Z coefficients of both formulations. Every coefficient and offset matches to within 3·10⁻¹³, and the set of ZZ pairs is identical.

### 9.1 The full Hamiltonians for P1

New (3 qubits: tables 0, 1, 2 at the single join step; λ = 1.431, offset = 2.670):
```
+0.315051 ZII    +0.627521 IZI    +0.465566 IIZ
+0.465566 ZZI    +0.715566 ZIZ    +0.465566 IZZ
```
Old (5 qubits: tables 0–2, then conditions 0–1 and 1–2; λ = 1, offset = 4.064):
```
+0.238561 ZIIII  +0.238561 IZIII  -0.210890 IIZII  +0.131385 IIIZI  +0.011439 IIIIZ
+0.500000 ZZIII  -0.900515 ZIZII  -0.650515 IZZII  -0.838046 ZIIZI  -0.838046 IZIZI
+1.265065 IIZZI  -0.500000 ZIIIZ  -0.750000 IZIIZ  +1.150515 IIZIZ  +1.088046 IIIZZ
```
P2–P5 are too long to print here; see the files.

### 9.2 Size and constants

| | Qubits (new / old) | One-qubit Z terms | Two-qubit ZZ terms (new / old) | Higher-order terms | λ (new / old) | Offset (new / old) | Largest \|coefficient\| (new / old) | Largest ÷ smallest \|coefficient\| (new / old) |
|---|---|---|---|---|---|---|---|---|
| P1 | 3 / 5 | one per qubit | 3 / 10 | none | 1.43 / 1 | 2.67 / 4.06 | 0.72 / 1.27 | **2.3 / 111** |
| P2 | 8 / 12 | one per qubit | 16 / 34 | none | 1.62 / 2 | 10.45 / 20.49 | 1.78 / 2.18 | 4.4 / 4.6 |
| P3 | 8 / 20 | one per qubit | 16 / 94 | none | 1.62 / 2 | 9.88 / 25.80 | 1.96 / 3.32 | **4.8 / 99** |
| P4 | 15 / 30 | one per qubit | 40 / 145 | none | 1.87 / 3 | 24.16 / 65.85 | 3.18 / 5.66 | 15 / 87 |
| P5 | 24 / 52 | one per qubit | 78 / 330 | none | 2.20 / 4 | 48.69 / 153.09 | 4.86 / 9.81 | 55 / 164 |

**What differs in the terms.**
- **Kept in both:**
  - a Z on every (table, step) bit;
  - ZZ between two tables at the same step (Penalty A, plus, in the new formulation, the join conditions);
  - ZZ between the same table at consecutive steps (Penalty B).
- **Only in the old formulation:**
  - a Z on every (condition, step) bit;
  - ZZ between a condition bit and each of its two tables (Penalty C);
  - ZZ between condition bits, and between condition and table bits, at the same step. These come from squaring the cost.

  That is where the extra 7 to 252 two-qubit terms come from. Each ZZ term is one two-qubit rotation per QAOA layer, so the old circuits have 2–6× more two-qubit gates.
- **Constants.**
  - The new formulation's coefficients all sit within a factor of 2–55 of each other. The old ones span up to a factor of 164. In P1 and P3 the old formulation has terms of 0.011 and 0.033 next to terms above 1.
  - The QUBO coefficients show the same thing. For P1, the new linear terms are all in [−4.72, −4.42] and the pairwise terms in [1.43, 2.86]. The old ones are in [−4, 2] and [−3.60, 5.06].
  - Small, widely spread coefficients are the first thing to be lost to noise or limited angle precision on hardware.
  - The offset only shifts every energy by the same amount. It doesn't change which bitstring is lowest, or the QAOA dynamics. It explains why the two formulations' energies can't be compared directly.

### 9.3 The exact eigensolver

`NumPyMinimumEigensolver` gives the lowest energy of the Hamiltonian. The script cross-checks it against a brute-force evaluation of all 2^n bitstrings, and the two agree to within 10⁻¹³ in every case. "Best plan's energy" is the energy of the cheapest join order, with condition bits set correctly for the old formulation.

| | Lowest energy (new / old) | Best plan's energy (new / old) | Is the lowest-energy bitstring a real plan? (new / old) | Cost it decodes to (new / old) | Bitstrings below the best plan (new / old) |
|---|---|---|---|---|---|
| P1 (best 15) | 1.176 / 1.031 | 1.176 / 1.383 | yes / **no**: condition 1–2 is on without table 2 | 15 / 15 (luck) | 0 / **2** |
| P2 (best 315) | 3.653 / 4.766 | 3.653 / 7.519 | yes / **no**: table 2 is forgotten at step 2 | 315 / **465** | 0 / **10** |
| P3 (best 29.4) | 2.334 / 2.725 | 2.334 / 2.725 | yes / yes | 29 / 29 | 0 / 0 |
| P4 (best 165) | 4.908 / too big (30 qubits) | 4.908 / 8.364 | yes / unknown | 165 / – | 0 / – |
| P5 (best 1605) | 8.340 / too big (52 qubits) | 8.340 / 19.430 | yes / unknown | 1605 / – | 0 / – |

**Should the eigensolver give the same number for both formulations?** No, and that is expected. They are two different energy functions: one squares the cost and the other doesn't, they use different λ, and the old one has extra bits. So their lowest energies are on different scales and differ (1.176 vs 1.031 on P1). What *should* be the same is the **plan** the lowest energy decodes to: the cheapest join order.
- In the new formulation this holds on all five problems.
- In the old one it holds on P3 only. On P1 and P2 the lowest energy is a fake plan, the cheats described in §4. On P1 the fake plan happens to decode to the best cost. On P2 it decodes to 465.

This is the main fix, now confirmed by an exact solver rather than by QAOA samples.

### 9.4 Did anything get worse?

Three measures of how clearly the Hamiltonian separates the best plan from the rest:

| | Gap: second-best plan − best plan (new / old) | Gap ÷ full energy range over all bitstrings (new / old) | Gap ÷ Σ\|coefficients\| (new / old) | Plan pairs ranked opposite to their cost (new / old) |
|---|---|---|---|---|
| P1 | 0.301 / 0.799 | 0.066 / 0.066 | 0.099 / 0.086 | 0 of 3 / 0 of 3 |
| P2 | 0.176 / 0.903 | 0.0065 / 0.0137 | 0.0088 / 0.0174 | 9 of 65 (13.8%) / 8 of 65 (12.3%) |
| P3 | 0.380 / 1.039 | 0.0140 / 0.0114 | 0.0194 / 0.0141 | 2 of 66 (3.0%) / 2 of 66 (3.0%) |
| P4 | 0.574 / 1.779 | 0.0081 / – | 0.0102 / 0.0092 | 177 of 1767 (10.0%) / 140 (7.9%) |
| P5 | 0.602 / 1.779 | 0.0038 / – | 0.0046 / 0.0036 | 8415 of 64521 (13.0%) / 6705 (10.4%) |

- **The raw gap is smaller in the new formulation** because the cost is no longer squared. That alone means nothing: QAOA sees the gap relative to the whole energy landscape, so only the two normalized columns matter. On those the formulations are about equal. The new one is ahead on P1, P3, P4 and P5, and behind on P2, where the old formulation's gap is moot because 10 fake bitstrings sit below its best plan.
- **One real, small degradation: ranking the non-best plans.** Both formulations put some pairs of plans in the opposite order to their true cost, because both are stand-ins for the sum of result sizes (Claim 4 in §5). The linear version does it slightly more often: 13.8% vs 12.3% on P2, 10.0% vs 7.9% on P4, 13.0% vs 10.4% on P5. The squared version weights big results more, which tracks the sum a little better. The **best plan is still ranked first by both** in all five problems. The degradation is in the ordering of mediocre plans, which affects how good the *second-best* answers are when QAOA misses the best one.

### 9.5 Pros and cons of each Hamiltonian

| | New ("qaoa") | Old ("q2o") |
|---|---|---|
| **Lowest energy = a real plan** | Yes on P1–P5 (exact solver); provable with λ_safe (§5 Correctness) | No on P1 and P2: λ too small for the squared cost |
| **Qubits** | R(R − 2), independent of the number of conditions | (R + P)(R − 2); 2–3.4× more on P4/P5 |
| **Two-qubit terms (circuit depth)** | 3–78 | 10–330 (2–6× more) |
| **Coefficient spread** | 2–55× | up to 164×, with some tiny terms (harder on noisy hardware) |
| **Penalty weight** | Derived from the data; one cheat saves at most one table's log size | Number of steps; not tied to how much a cheat saves |
| **Ranking of non-best plans** | Slightly worse (product of sizes) | Slightly better (squared logs weight big results) |
| **Raw best-vs-second gap** | Smaller (no square) | Larger, but on P1/P2 fake states sit below it |
| **Fits in a 25-qubit simulation** | Up to 6 tables | Up to 4 tables |

**Bottom line:** the new Hamiltonian is smaller, shallower, has better-conditioned coefficients, and its lowest energy is always the best plan. The price is a slightly less faithful ordering of mediocre plans. No measure shows a loss for the best plan.

### 9.6 Circuit depth and how the Pauli terms are spread

`base/circuit_depth_report.py`, run in the container, builds the QAOA circuit for each problem and formulation at 1–3 layers and decomposes it into CNOT, RZ, RX and H gates without reordering them. It only builds the circuits; it doesn't simulate them. The numbers are in `docs/hamiltonians/circuit_depth.json`, and the term counts come from the `P<k>_<f>.txt` files.

**Circuit depth at 2 layers (the setting all our runs use)**

| Problem | Qubits (new / old) | Depth (new / old) | CNOT depth (new / old) | CNOTs (new / old) |
|---|---|---|---|---|
| P1 | 3 / 5 | 23 / 41 | 12 / 24 | 12 / 40 |
| P2 | 8 / 12 | 38 / 53 | 22 / 32 | 64 / 136 |
| P3 | 8 / 20 | 38 / 89 | 22 / 56 | 64 / 376 |
| P4 | 15 / 30 | 65 / 107 | 40 / 68 | 160 / 580 |
| P5 | 24 / 52 | 80 / 140 | 50 / 90 | 312 / 1320 |

- **The depth changes with every workload.** It grows with the number of tables and, in the old formulation, with the number of conditions too. P2 and P3 give identical circuits in the new formulation (depth 38) but very different ones in the old (53 vs 89), because P3 has 6 conditions.
- **The depth grows linearly with the number of layers.** For example, at 1 → 2 → 3 layers, P5 goes 57 → 80 → 115 (new) and 96 → 140 → 196 (old). The CNOT count is exactly proportional to the layers.
- **The new formulation is shallower on every workload:** 1.4–2.4× shallower, with 2.1–5.9× fewer CNOTs. In our noiseless simulations this only changes run time. On real hardware every CNOT adds error, and the old formulation's denser connectivity would also need more SWAP gates to route on a chip with limited connectivity. The routing hasn't been measured.

**How depth relates to the number of terms**
- **The CNOT count is fixed exactly by the ZZ terms:** CNOTs = 2 × (ZZ terms) × (layers). Each ZZ term becomes CNOT, RZ, CNOT, and Z terms need no CNOTs. For example, P5 old: 2 × 330 × 2 = 1320.
- **The depth depends on the busiest qubit, not on the total.** Gates on different qubits run in parallel. So one layer needs at least 2 × (the most couplings on any single qubit) CNOT steps. That maximum is 2 / 4 / 4 / 6 / 7 for the new formulation and 4 / 6 / 10 / 11 / 14 for the old. That's why the depth gap (1.4–2.4×) is smaller than the CNOT gap (2–6×).
- Qiskit's default term order comes out at about 1.5–2.5× that lower bound in both formulations, so reordering the terms could shorten both circuits. That hasn't been tried.

**How the Pauli terms are spread**

| Problem | Qubits (new / old) | Z | ZZ | Share of qubit pairs coupled | Most couplings on one qubit | ZZ coefficient size range |
|---|---|---|---|---|---|---|
| P1 | 3 / 5 | 3 / 5 | 3 / 10 | 100% / 100% | 2 / 4 | 0.47–0.72 / 0.50–1.27 |
| P2 | 8 / 12 | 8 / 12 | 16 / 34 | 57% / 52% | 4 / 6 | 0.41–0.81 / 0.50–1.96 |
| P3 | 8 / 20 | 8 / 20 | 16 / 94 | 57% / 49% | 4 / 10 | 0.41–0.81 / 0.03–1.96 |
| P4 | 15 / 30 | 15 / 30 | 40 / 145 | 38% / 33% | 6 / 11 | 0.47–0.93 / 0.08–2.76 |
| P5 | 24 / 52 | 24 / 52 | 78 / 330 | 28% / 25% | 7 / 14 | 0.55–1.10 / 0.06–3.70 |

- **Only Z and ZZ terms.** Both formulations have exactly one Z term per qubit plus the ZZ couplings, and no terms on three or more qubits.
- **The same shape, at different sizes.** In both, the qubits of one join step are all coupled to each other, and each table is chained to itself across steps. The old formulation has more qubits per step (R + P instead of R), so its blocks are larger.
- **Coefficient sizes.** The new coefficients stay within about 2× of each other. The old ones span up to 60× (P3): the square adds many small condition–condition and table–condition couplings (as small as 0.03) next to large table–table couplings (up to 3.7). A small coupling costs as many gates as a large one, but changes the energy very little.

**Is the number of terms tied to the number of qubits? Yes, by the same rule in both formulations.** With J join steps and b qubits per step (so n = b·J qubits in total):

ZZ terms = J · b(b − 1)/2 + R(J − 1) ≈ n² / (2J)

Here b = R in the new formulation and R + P in the old. The rule reproduces every count in the table: for P3 old, 2 · 10 · 9/2 + 4 = 94.
- **For a fixed number of join steps, the ZZ count grows with the square of the qubit count, in both formulations.** The old one has more terms only because it has more qubits per step. The ratio is about ((R + P)/R)², for example (10/4)² ≈ 6 for P3.
- **The share of qubit pairs that are coupled is about 1/J** in both: 100% at J = 1, about 50% at J = 2, 33–38% at J = 3, and 25–28% at J = 4. Larger queries give sparser coupling graphs relative to their size, though the number of couplings still grows.

**How the term count scales**

| Grows with | New ("qaoa") | Old ("q2o") |
|---|---|---|
| Qubits n, for fixed J | ≈ n² / (2J) | ≈ n² / (2J): the same rule, but more qubits |
| Workload (R tables, P conditions) | ≈ R³/2, **independent of P** | ≈ R(R + P)²/2: grows with the conditions too |
| Circuit depth (p layers) | The Hamiltonian doesn't change; the circuit applies every term once per layer: CNOTs = 2·ZZ·p, rotations = (Z + ZZ + qubits)·p, SPIQ angles = the same count | same |

For example, the SPIQ circuit for P5 has 126 angles per layer in the new formulation and 434 in the old.
