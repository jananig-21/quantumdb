# Checks the closed-form QUBO/Ising coefficients in docs/problem-formulation.md §9.0 against QUBOGenerator, P1-P5.
# Classical only. Run from base/: python3 test_hamiltonian_form.py
import itertools
import numpy as np
import Scripts.ProblemGenerator as G, Scripts.QUBOGenerator as Q


def closed_form(f, card, pred, sel):
    c = np.log10(card); s = np.log10(sel)
    R, P, J = len(c), len(s), len(c) - 2
    lam = 1.1 * c.max() if f == "qaoa" else J
    nv = R * J; n = nv + (P * J if f == "q2o" else 0)
    v = lambda t, j: t * J + j
    w = lambda p, j: nv + p * J + j
    h = np.zeros(n); Jm = {}
    add = lambda a, b, x: Jm.__setitem__((min(a, b), max(a, b)), Jm.get((min(a, b), max(a, b)), 0) + x)
    const = lam * sum((j + 2) ** 2 for j in range(J))
    for t in range(R):
        for j in range(J):
            h[v(t, j)] = lam * (1 - 2 * (j + 2)) + lam * (j < J - 1) + (c[t] if f == "qaoa" else c[t] ** 2)
            if j < J - 1: add(v(t, j), v(t, j + 1), -lam)
    for j in range(J):
        for t, u in itertools.combinations(range(R), 2):
            add(v(t, j), v(u, j), 2 * lam + (0 if f == "qaoa" else 2 * c[t] * c[u]))
        for p, (a, b) in enumerate(pred):
            if f == "qaoa":
                add(v(a, j), v(b, j), s[p])
            else:
                h[w(p, j)] = 2 * lam + s[p] ** 2
                for t in range(R):
                    add(v(t, j), w(p, j), 2 * c[t] * s[p] - lam * (t in (a, b)))
                for q in range(p + 1, P):
                    add(w(p, j), w(q, j), 2 * s[p] * s[q])
    # Ising: x = (1 - Z)/2
    zz = {k: x / 4 for k, x in Jm.items()}
    z = -h / 2
    for (a, b), x in Jm.items():
        z[a] -= x / 4; z[b] -= x / 4
    off = const + h.sum() / 2 + sum(Jm.values()) / 4
    if f == "qaoa":  # explicit Z formula from the doc
        for t in range(R):
            St = sum(s[p] for p, ab in enumerate(pred) if t in ab)
            for j in range(J):
                zf = -0.5 * (c[t] + St / 2 + lam * (R - 2 * j - 4) + lam / 2 * ((j < J - 1) - (j > 0)))
                assert abs(zf - z[v(t, j)]) < 1e-9, (t, j, zf, z[v(t, j)])
    else:  # explicit q2o Z formulas from the doc; C + S_all = log10 of the full join
        C, S_all = c.sum(), s.sum()
        for t in range(R):
            d = sum(t in ab for ab in pred)
            for j in range(J):
                zf = -0.5 * (c[t] * (C + S_all) + lam * (R - 2 * j - 4 - d / 2) + lam / 2 * ((j < J - 1) - (j > 0)))
                assert abs(zf - z[v(t, j)]) < 1e-9, (t, j, zf, z[v(t, j)])
        for p in range(P):
            for j in range(J):
                assert abs(-0.5 * (lam + s[p] * (C + S_all)) - z[w(p, j)]) < 1e-9, (p, j)
    return n, zz, z, off


for i in range(5):
    card, pred, sel = G.get_join_ordering_problem(f"ExperimentalAnalysis/IBMQ/QPUPerformance/Problems/JSON/{i}_predicates")
    for f in ("qaoa", "q2o"):
        qubo, _ = Q.FORMULATIONS[f](card, pred, sel)
        op, offset = qubo.to_ising()
        n, zz, z, off = closed_form(f, card, pred, sel)
        got_z, got_zz = np.zeros(n), {}
        for lab, cf in zip(op.primitive.paulis.to_labels(), op.primitive.coeffs.real):
            q = [n - 1 - k for k, ch in enumerate(lab) if ch == "Z"]
            if len(q) == 1: got_z[q[0]] += cf
            else: got_zz[tuple(sorted(q))] = got_zz.get(tuple(sorted(q)), 0) + cf
        zz = {k: x for k, x in zz.items() if abs(x) > 1e-12}
        assert set(zz) == {k for k, x in got_zz.items() if abs(x) > 1e-12}, (f, i)
        err = max([abs(zz[k] - got_zz[k]) for k in zz] + list(np.abs(z - got_z)) + [abs(off - offset)])
        assert err < 1e-9, (f, i, err)
        print(f"P{i+1} {f}: n={n} ZZ={len(zz)} offset={off:.4f} max err {err:.1e} OK")
