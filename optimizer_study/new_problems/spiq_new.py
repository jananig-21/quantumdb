"""SPIQ points for the 12 new problems, in both formulations, with the team's own SPIQ code (unmodified).

Same call as ../many_problems/spiq_many.py (400 generations, 4 starts, fast stim builder that produces identical
circuits), except the QUBO:
  q2o   the old formulation, base/Scripts/QUBOGenerator.generate_IBMQ_QUBO_for_left_deep_trees_v2 (study repo)
  qaoa  Divya's new formulation, generate_qaoa_qubo_left_deep from the team's `formulation` branch
For q2o, problems that already have a Study 1 SPIQ point (many_problems/spiq/<id>.json, same settings) reuse it.

Usage (SPIQ env): python spiq_new.py <q2o|qaoa> [--n_proc 2]
Output: spiq_<form>/<id>.json (resumable)
"""
import argparse, contextlib, importlib.util, io, json, os, shutil, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MP = os.path.join(HERE, "..", "many_problems")
sys.path.insert(0, MP)
import spiq_many as SM                     # gate_list, check_gate_list, BASE on sys.path

FORMULATION_QG = os.path.expanduser("~/Downloads/umich-db-quantumdb-wt/formulation/base/Scripts/QUBOGenerator.py")


def qubo_for(form, p):
    if form == "q2o":
        import Scripts.QUBOGenerator as QG
        return QG.generate_IBMQ_QUBO_for_left_deep_trees_v2(p["card"], [tuple(x) for x in p["pred"]], p["pred_sel"])
    spec = importlib.util.spec_from_file_location("QG_formulation", FORMULATION_QG)
    QG = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(QG)
    return QG.generate_qaoa_qubo_left_deep(p["card"], [tuple(x) for x in p["pred"]], p["pred_sel"])


def one(form, p, n_gens, n_proc, n_starts, out):
    path = os.path.join(out, p["id"] + ".json")
    if os.path.exists(path):
        return "skip"
    old = os.path.join(MP, "spiq", p["id"] + ".json")
    if form == "q2o" and os.path.exists(old):
        shutil.copy(old, path)
        return "reused Study 1 point"
    import spiq_initialization as SP
    import fast_clapton
    fast_clapton.apply()
    fast_clapton.apply_single_core_eval()
    t0 = time.time()
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        qubo, penalty = qubo_for(form, p)
        os.makedirs(os.path.join(out, "traces"), exist_ok=True)
        objs = SP.build_vanilla_spiq_objects(qubo, reps=2)
        fast_clapton.verify(objs[4], trials=20)
        res = SP.run_spiq_initialization(qubo=qubo, reps=2, n_gens=n_gens, n_proc=n_proc, n_starts=n_starts,
                                         n_rounds=1, out_file=os.path.join(out, "traces", p["id"] + ".txt"))
    pc = res["pcirc"]
    gl = SM.gate_list(pc)
    by = dict(zip(res["relaxed_param_names"], res["relaxed_initial_point"]))
    x0 = [by[q.name] for q in pc.parameters]
    err = max(SM.check_gate_list(pc, gl, x0), SM.check_gate_list(pc, gl, np.random.default_rng(0).uniform(-np.pi, np.pi, len(x0))))
    assert err < 1e-8, (p["id"], err)
    rec = dict(problem=p["id"], formulation=form, penalty_weight=penalty, n_gens=n_gens, n_starts=n_starts, reps=2,
               relaxed_initial_point=x0, relaxed_param_names=gl["param_names"], energy_best=res["energy_best"],
               energy_best_qubo=res["energy_best_qubo"], ising_offset=res["ising_offset"], num_params=len(x0),
               gates=gl, gatelist_max_dp=err, wall_s=time.time() - t0)
    json.dump(rec, open(path + ".tmp", "w"))
    os.replace(path + ".tmp", path)
    return f"done {gl['num_qubits']}q {len(x0)} angles E_spiq={res['energy_best_qubo']:.3f} {time.time() - t0:.0f}s"


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.set_start_method("fork")
    ap = argparse.ArgumentParser()
    ap.add_argument("form", choices=["q2o", "qaoa"])
    ap.add_argument("--n_gens", type=int, default=400)
    ap.add_argument("--n_proc", type=int, default=1)
    ap.add_argument("--n_starts", type=int, default=4)
    a = ap.parse_args()
    out = os.path.join(HERE, f"spiq_{a.form}")
    os.makedirs(out, exist_ok=True)
    probs = [json.loads(l) for l in open(os.path.join(HERE, "new_problems.jsonl"))]
    os.chdir(SM.BASE)
    for i, p in enumerate(probs, 1):
        print(f"[{i}/{len(probs)}] {a.form} {p['id']}: {one(a.form, p, a.n_gens, a.n_proc, a.n_starts, out)}", flush=True)
