# Circuit depth, CX count and ZZ-coupling degree of the QAOA ansatz, both formulations, P1-P5, reps 1-3.
# Builds and transpiles the circuits only (no simulation). Run in the container from base/:
#   python3 circuit_depth_report.py   -> ../docs/hamiltonians/circuit_depth.json
import json, collections
from qiskit import transpile
from qiskit.circuit.library import QAOAAnsatz
import Scripts.ProblemGenerator as G, Scripts.QUBOGenerator as Q
rows = []
for i in range(5):
    card, pred, sel = G.get_join_ordering_problem(f"ExperimentalAnalysis/IBMQ/QPUPerformance/Problems/JSON/{i}_predicates")
    for f in ("qaoa", "q2o"):
        qubo, _ = Q.FORMULATIONS[f](card, pred, sel)
        op, _ = qubo.to_ising()
        n = qubo.get_num_vars()
        labels = op.primitive.paulis.to_labels()
        deg = collections.Counter()
        nzz = 0
        for l in labels:
            q = [k for k, ch in enumerate(l) if ch == "Z"]
            if len(q) == 2:
                nzz += 1; deg[q[0]] += 1; deg[q[1]] += 1
        maxdeg = max(deg.values())
        for reps in (1, 2, 3):
            qc = transpile(QAOAAnsatz(op, reps=reps), basis_gates=["cx", "rz", "rx", "h"], optimization_level=0)
            ops = qc.count_ops()
            d2 = qc.depth(lambda g: g[0].num_qubits == 2)
            r = dict(problem=f"P{i+1}", f=f, n=n, z=len(labels) - nzz, zz=nzz, maxdeg=maxdeg,
                     avgdeg=round(2 * nzz / n, 2), reps=reps, depth=qc.depth(), cx=ops.get("cx", 0),
                     cx_depth=d2, params=qc.num_parameters)
            rows.append(r); print(json.dumps(r))
json.dump(rows, open("../docs/hamiltonians/circuit_depth.json", "w"), indent=1)
