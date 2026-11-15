"""Gate inventory of the QAOA circuits, old ("q2o") vs new ("qaoa") formulation, P1-P5. No simulation.

Three views of each circuit (reps 1-3 for the first two, reps 2 for the third):
  plain    QAOAAnsatz in {h, rx, rz, cx}: the textbook circuit.
  spiq     the relaxed circuit SPIQ hands to COBYLA/NFT (build_vanilla_spiq_objects): what our runs execute.
  hardware the spiq circuit compiled for FakeBrooklyn (65-qubit heavy-hex, basis cx/rz/sx/x), optimization
           level 3, best of 5 layout seeds. Adds routing SWAPs (as extra CNOTs) and an estimated success
           probability = product of (1 - error) over every gate and readout.

Run in the container (needs stim/pygad from spiq-site):
  PYTHONPATH=/home/repro/spiq-site:/home/repro/sigmod-repro python3 circuit_gate_report.py
Writes ../docs/hamiltonians/circuit_gates.json.
"""
import contextlib, io, json
from qiskit import transpile
from qiskit.circuit.library import QAOAAnsatz
from qiskit.test.mock import FakeBrooklyn
import Scripts.ProblemGenerator as G, Scripts.QUBOGenerator as Q
from spiq_initialization import build_vanilla_spiq_objects

backend = FakeBrooklyn()
props = backend.properties()


def summary(qc):
    ops = dict(qc.count_ops())
    ops.pop("barrier", None)
    return dict(ops=ops, total=sum(v for k, v in ops.items() if k != "measure"),
                two_q=sum(v for k, v in ops.items() if k in ("cx", "cz", "ecr", "swap")),
                depth=qc.depth(), two_q_depth=qc.depth(lambda g: g[0].num_qubits == 2),
                params=qc.num_parameters)


def esp(qc):
    p = 1.0
    for inst, qargs, _ in qc.data:
        q = [x.index for x in qargs]
        if inst.name in ("cx", "sx", "x"):
            p *= 1 - props.gate_error(inst.name, q)
        elif inst.name == "measure":
            p *= 1 - props.readout_error(q[0])
    return p


rows = []
for i in range(5):
    card, pred, sel = G.get_join_ordering_problem(f"ExperimentalAnalysis/IBMQ/QPUPerformance/Problems/JSON/{i}_predicates")
    for f in ("qaoa", "q2o"):
        with contextlib.redirect_stdout(io.StringIO()):
            qubo, _ = Q.FORMULATIONS[f](card, pred, sel)
        op, _ = qubo.to_ising()
        labels = op.primitive.paulis.to_labels()
        zz = sum(l.count("Z") == 2 for l in labels)
        base = dict(problem=f"P{i+1}", f=f, R=len(card), P=len(pred), n=qubo.get_num_vars(), z=len(labels) - zz, zz=zz)
        for reps in (1, 2, 3):
            plain = transpile(QAOAAnsatz(op, reps=reps), basis_gates=["cx", "rz", "rx", "h"], optimization_level=0)
            with contextlib.redirect_stdout(io.StringIO()):
                spiq = build_vanilla_spiq_objects(qubo, reps)[3]
            r = dict(base, reps=reps, plain=summary(plain), spiq=summary(spiq))
            if reps == 2:
                m = spiq.copy(); m.measure_all()
                best = min((transpile(m, backend, optimization_level=3, seed_transpiler=s) for s in range(5)),
                           key=lambda c: c.count_ops().get("cx", 0))
                r["hardware"] = dict(summary(best), swaps=(best.count_ops().get("cx", 0) - r["spiq"]["two_q"]) / 3,
                                     esp=esp(best))
            rows.append(r)
            print(json.dumps(r), flush=True)
json.dump(rows, open("../docs/hamiltonians/circuit_gates.json", "w"), indent=1)
