"""Pick 12 generated 4-relation queries (3 per graph shape) not used by the team, whose optimal plan is reachable
with the team's fallback decoder (otherwise P(optimal plan) is 0 for every state). Deterministic: first eligible by id."""
import json, os, sys
sys.path[:0] = [os.path.join(os.path.dirname(__file__), "..", "new_code_study")]
import study as S
probs = [json.loads(l) for l in open(os.path.join(os.path.dirname(__file__), "..", "many_problems", "problems.jsonl"))]
sel = []
for g in ["CHAIN", "STAR", "CYCLE", "CLIQUE"]:
    k = 0
    for p in sorted([p for p in probs if p["source"] == "generated" and p["graph"] == g and p["relations"] == 4], key=lambda p: p["id"]):
        card, pred, s = p["card"], [tuple(x) for x in p["pred"]], p["pred_sel"]
        raw, fb, opt = S.decode_tables(card, pred, s, 8)
        if (abs(fb - opt) < 1e-9 * max(1, opt)).any():
            sel.append(dict(p, opt_cost=opt)); k += 1
        if k == 3: break
json.dump(sel, open("selected.json", "w"), indent=1)
for p in sel: print(p["id"], p["num_qubits"], "old qubits;", p["card"], p["pred"], p["pred_sel"], "opt", p["opt_cost"])

# The 1,088 set has no generated 4-relation CLIQUEs, so make new ones exactly as make_problems.py does (seed 43000+).
import itertools, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "many_problems"))
import SteinbrunnQueryGenerator as SQG
k, seed = 0, 43000
while k < 3:
    np.random.seed(seed)
    card, sl = SQG.produce_steinbrunn_query(4, "CLIQUE", "Random", max_cardinality=100000)
    pred, ps = [], []
    for t1, t2 in itertools.combinations(range(len(card)), 2):
        if sl[t1][t2] < 1.0:
            pred.append([t1, t2]); ps.append(sl[t1][t2])
    p = dict(id=f"Gen_CLIQUE_4r_s{seed}", source="generated_new", graph="CLIQUE", relations=4, card=[float(c) for c in card],
             pred=pred, pred_sel=ps, seed=seed, num_qubits=(4 + len(pred)) * 2)
    raw, fb, opt = S.decode_tables(card, [tuple(x) for x in pred], ps, 8)
    if (abs(fb - opt) < 1e-9 * max(1, opt)).any():
        sel.append(dict(p, opt_cost=opt)); k += 1
        print(p["id"], p["num_qubits"], "old qubits;", p["card"], pred, "opt", opt)
    seed += 1
json.dump(sel, open("selected.json", "w"), indent=1)
with open("new_problems.jsonl", "w") as f:
    for p in sel: f.write(json.dumps(p) + "\n")
print(len(sel), "problems")
