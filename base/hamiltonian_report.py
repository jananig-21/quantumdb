"""Hamiltonian + eigensolver comparison, old (q2o) vs new (qaoa) formulation, P1-P5. Classical only, no QAOA.

Run from base/ in the container (needs qiskit 0.19). Writes the full Ising Hamiltonians to
../docs/hamiltonians/P<k>_<f>.txt, a summary.json, and prints one line per (problem, formulation).
Up to 24 qubits it also runs NumPyMinimumEigensolver and evaluates all 2^n energies.
"""
import itertools, os, json
import numpy as np
from qiskit.algorithms import NumPyMinimumEigensolver
import Scripts.Postprocessing as P, Scripts.ProblemGenerator as G, Scripts.QUBOGenerator as Q

MAX_BRUTE = 24
OUT = "../docs/hamiltonians"
os.makedirs(OUT, exist_ok=True)
rows = []


def all_energies(qubo, chunk=1 << 20):
    n = qubo.get_num_vars()
    lin = qubo.objective.linear.to_array(); quad = qubo.objective.quadratic.to_array(); c = qubo.objective.constant
    shifts = np.arange(n)
    out = np.empty(1 << n)
    for s in range(0, 1 << n, chunk):
        idx = np.arange(s, min(s + chunk, 1 << n))
        x = ((idx[:, None] >> shifts) & 1).astype(np.float64)
        out[s:s + len(idx)] = c + x @ lin + np.einsum("ij,ij->i", x @ quad, x)
    return out


for i in range(5):
    card, pred, sel = G.get_join_ordering_problem(f"ExperimentalAnalysis/IBMQ/QPUPerformance/Problems/JSON/{i}_predicates")
    R = len(card); J = R - 2
    orders = list(itertools.permutations(range(R)))
    costs = {o: P.get_costs_for_leftdeep_tree(list(o), card, pred, sel, {}) for o in orders}
    for f in ("qaoa", "q2o"):
        qubo, lam = Q.FORMULATIONS[f](card, pred, sel)
        n = qubo.get_num_vars()
        op, offset = qubo.to_ising()
        labels = op.primitive.paulis.to_labels(); coeffs = op.primitive.coeffs.real
        weight = [l.count("Z") for l in labels]
        with open(f"{OUT}/P{i+1}_{f}.txt", "w") as fh:
            fh.write(f"# P{i+1} {f}: {n} qubits, lambda={lam:.6f}, Ising offset={offset:.6f}\n")
            fh.write("# H = offset + sum coeff * Pauli (qiskit labels, qubit 0 rightmost)\n")
            for l, cf in sorted(zip(labels, coeffs), key=lambda t: (t[0].count('Z'), t[0][::-1])):
                fh.write(f"{cf:+.6f} {l}\n")
        lin = qubo.objective.linear.to_array(); quad = qubo.objective.quadratic.to_array()
        r = dict(problem=f"P{i+1}", formulation=f, qubits=n, lam=lam, offset=float(offset),
                 n_terms=len(labels), n_Z=weight.count(1), n_ZZ=weight.count(2), n_other=sum(w > 2 for w in weight),
                 max_abs_coeff=float(np.abs(coeffs).max()), min_abs_coeff=float(np.abs(coeffs[np.abs(coeffs) > 1e-12]).min()),
                 sum_abs_coeff=float(np.abs(coeffs).sum()),
                 qubo_lin_range=[float(lin.min()), float(lin.max())],
                 qubo_quad_range=[float(quad[quad != 0].min()), float(quad[quad != 0].max())] if (quad != 0).any() else None)

        # Real plan: valid table bits and, for q2o, condition bits set honestly. Energy of every real plan (table bits from the order, condition bits honest for q2o)
        def plan_bits(o):
            sets = [set(o[:j + 2]) for j in range(J)]
            x = [int(t in sets[j]) for t in range(R) for j in range(J)]
            if f == "q2o":
                x += [int(a in sets[j] and b in sets[j]) for (a, b) in pred for j in range(J)]
            return x
        plan_e = {o: qubo.objective.evaluate(plan_bits(o)) for o in orders}
        best_cost = min(costs.values())
        e_best = min(e for o, e in plan_e.items() if costs[o] == best_cost)
        e_second = min((e for o, e in plan_e.items() if costs[o] > best_cost), default=None)
        r.update(best_cost=best_cost, E_best_plan=e_best, E_second_plan=e_second,
                 plan_gap=None if e_second is None else e_second - e_best,
                 rank_ok=all((costs[a] < costs[b]) <= (plan_e[a] < plan_e[b] + 1e-9) for a in orders for b in orders))

        if n <= MAX_BRUTE:
            ev = NumPyMinimumEigensolver().compute_minimum_eigenvalue(op)
            r["eigensolver_qubo"] = float(ev.eigenvalue.real) + float(offset)
            E = all_energies(qubo)
            k = int(np.argmin(E))
            x = [(k >> b) & 1 for b in range(n)]
            v = np.array(x[:R * J]).reshape(R, J)
            valid = all(v[:, j].sum() == j + 2 for j in range(J)) and (v[:, :-1] <= v[:, 1:]).all()
            _, sol = P.readout([[[x, 1, 0.0, "", 1.0]], 0.0], card, pred, sel, {})
            r.update(brute_min=float(E.min()), E_max=float(E.max()), ground_bits="".join(map(str, x)),
                     ground_tables_valid=bool(valid), ground_is_real_plan=bool(abs(E.min() - e_best) < 1e-9), ground_cost_raw=sol[0][2], ground_cost_fallback=sol[1][2],
                     n_below_best_plan=int((E < e_best - 1e-9).sum()),
                     gap_over_range=None if e_second is None else (e_second - e_best) / (E.max() - E.min()))
        rows.append(r)
        print(json.dumps(r))

json.dump(rows, open(f"{OUT}/summary.json", "w"), indent=1)
