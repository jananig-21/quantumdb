"""SPIQ initialization for every problem in problems.jsonl, with the team's own code.

Runs in the SPIQ environment (requirements-spiq.txt: qiskit < 1.0, qiskit-aer 0.13.3,
qiskit-algorithms 0.3.1, stim, pygad 3.3.1). For each problem it calls the team's
base/spiq_initialization.run_spiq_initialization UNMODIFIED (relaxed per-gate QAOA ansatz,
clapton genetic search over Clifford angles k*pi/2, best point mapped back to angles), exactly
as its CLI does, except that the QUBO is built from the problem's card/pred/pred_sel instead of
the IBMQ problem folders.

Output per problem, spiq/<id>.json:
  relaxed_initial_point, relaxed_param_names, energy_best (Ising), energy_best_qubo, ising_offset,
  ks_best_raw, n_gens / n_starts / n_proc, wall time, and `gates`: the relaxed circuit as a plain
  gate list (../fastsim.compile_circuit format) so the optimizer phase - which runs in the older
  qiskit-terra 0.19 environment that cannot read QPY written by qiskit 0.46 - can simulate it.
  The team's own self-check (qiskit Statevector at the SPIQ point == stim energy) runs inside
  run_spiq_initialization; we additionally check the gate list against Statevector here.

Usage (SPIQ env):  python spiq_many.py [--n_gens 400] [--n_proc 8] [--source team] [--only ID ...]
"""
import argparse, contextlib, io, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.abspath(os.path.join(HERE, ".."))
ROOT = os.path.abspath(os.path.join(STUDY, ".."))
BASE = os.path.join(ROOT, "base")
sys.path[:0] = [BASE, ROOT, STUDY]
OUT = os.environ.get("SPIQ_OUT", os.path.join(HERE, "spiq"))


def gate_list(qc):
    """../fastsim.compile_circuit, returned in a JSON-friendly form."""
    from fastsim import compile_circuit
    n, ops = compile_circuit(qc)
    return dict(num_qubits=n, ops=[[name, qs, None if arg is None else list(arg)] for name, qs, arg in ops],
                param_names=[p.name for p in qc.parameters])


def check_gate_list(qc, gl, theta):
    from fastsim import circuit_probs
    from qiskit.quantum_info import Statevector
    ops = [(o[0], o[1], None if o[2] is None else tuple(o[2])) for o in gl["ops"]]
    p_fast = circuit_probs((gl["num_qubits"], ops), np.asarray(theta, float))
    p_ref = Statevector.from_instruction(qc.assign_parameters(dict(zip(qc.parameters, theta)))).probabilities()
    return float(np.max(np.abs(p_fast - p_ref)))


def spiq_one(p, n_gens, n_proc, n_starts):
    import Scripts.QUBOGenerator as QUBOGenerator
    import spiq_initialization as SP                      # the team's module, unmodified
    import fast_clapton                                   # identical stim circuits, ~500x faster to build
    fast_clapton.apply()
    fast_clapton.apply_single_core_eval()                # no per-generation Manager process (n_proc/n_starts == 1)
    path = os.path.join(OUT, p["id"] + ".json")
    if os.path.exists(path):
        return "skip"
    t0 = time.time()
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        qubo, penalty = QUBOGenerator.generate_IBMQ_QUBO_for_left_deep_trees_v2(
            p["card"], [tuple(x) for x in p["pred"]], p["pred_sel"])
        os.makedirs(os.path.join(OUT, "traces"), exist_ok=True)
        objs = SP.build_vanilla_spiq_objects(qubo, reps=2)
        fast_clapton.verify(objs[4], trials=20)          # patched builder == clapton's, on this circuit
        res = SP.run_spiq_initialization(qubo=qubo, reps=2, n_gens=n_gens, n_proc=n_proc, n_starts=n_starts,
                                         n_rounds=1, out_file=os.path.join(OUT, "traces", p["id"] + ".txt"))
    pc = res["pcirc"]
    gl = gate_list(pc)
    by_name = dict(zip(res["relaxed_param_names"], res["relaxed_initial_point"]))
    x0 = [by_name[q.name] for q in pc.parameters]
    err0 = check_gate_list(pc, gl, x0)
    err1 = check_gate_list(pc, gl, np.random.default_rng(0).uniform(-np.pi, np.pi, len(x0)))
    assert max(err0, err1) < 1e-8, (p["id"], err0, err1)
    diag = [l for l in log.getvalue().splitlines() if "rel_err" in l]
    rec = dict(problem=p["id"], penalty_weight=penalty, n_gens=n_gens, n_starts=n_starts, n_proc=n_proc, reps=2,
               relaxed_initial_point=x0, relaxed_param_names=gl["param_names"],
               energy_best=res["energy_best"], energy_best_qubo=res["energy_best_qubo"],
               ising_offset=res["ising_offset"], ks_best_raw=res["ks_best_raw"],
               team_selfcheck=diag[0] if diag else None, gatelist_max_dp=max(err0, err1),
               num_params=len(x0), gates=gl, wall_s=time.time() - t0)
    tmp = path + ".tmp"
    json.dump(rec, open(tmp, "w"))
    os.replace(tmp, path)
    return f"done {gl['num_qubits']}q {len(x0)} angles E_spiq={res['energy_best_qubo']:.3f} {time.time() - t0:.0f}s"


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.set_start_method("fork")   # clapton's worker pool relies on fork (the team ran on Linux)
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_gens", type=int, default=400)
    ap.add_argument("--n_proc", type=int, default=8)
    ap.add_argument("--n_starts", type=int, default=4)
    ap.add_argument("--source", choices=["all", "team", "generated"], default="all")
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--shard", default="0/1", help="i/k: this process handles problems i, i+k, i+2k, ...")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    os.chdir(BASE)
    probs = [json.loads(l) for l in open(os.path.join(HERE, "problems.jsonl"))]
    if a.only:
        probs = [p for p in probs if p["id"] in a.only]
    if a.source != "all":
        probs = [p for p in probs if (p["source"] == "generated") == (a.source == "generated")]
    probs.sort(key=lambda p: (p["source"] == "generated", p["num_qubits"]))   # team queries first
    si, sk = map(int, a.shard.split("/"))
    probs = probs[si::sk]
    for i, p in enumerate(probs, 1):
        print(f"[{i}/{len(probs)}] {p['id']}: {spiq_one(p, a.n_gens, a.n_proc, a.n_starts)}", flush=True)
