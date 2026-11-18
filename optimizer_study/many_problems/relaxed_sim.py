"""Fast exact simulator for SPIQ's relaxed QAOA circuits (one angle per gate).

Input: the gate list of ../fastsim.compile_circuit (gates cx / rz / s / sx / ...; parameter
arguments ('param', index, coeff) or ('const', value)). ../fastsim.circuit_probs applies these
one gate at a time - ~900 full passes over the 2^n state for a 20-qubit problem. This module
groups the gates into two kinds of layers and applies each layer at once:

  * single-qubit layer  - a run of single-qubit gates. Gates on the same qubit are multiplied into
                          one 2x2 matrix; the layer is applied in blocks of 4 qubits (~5 passes).
  * phase layer         - a run of CX and diagonal gates (rz, s, sdg, z) starting with a CX.
                          Tracking each qubit's parity mask through the CX network turns every
                          diagonal gate into a term w * (-1)^{popcount(x & mask)}; the whole phase
                          is then a Walsh-Hadamard transform of the sparse coefficient vector
                          (n passes), followed by the net CX permutation (identity for QAOA's
                          CX-RZ-CX ladders, but handled in general).

Same state, same parameter order as fastsim.circuit_probs; verify() checks this to 1e-10.
"""
import numpy as np

_1Q = {
    "s":   lambda th: np.array([[1, 0], [0, 1j]]),
    "sdg": lambda th: np.array([[1, 0], [0, -1j]]),
    "sx":  lambda th: 0.5 * np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]]),
    "h":   lambda th: np.array([[1, 1], [1, -1]]) / np.sqrt(2),
    "x":   lambda th: np.array([[0, 1], [1, 0]]),
    "z":   lambda th: np.array([[1, 0], [0, -1]]),
    "rz":  lambda th: np.array([[np.exp(-0.5j * th), 0], [0, np.exp(0.5j * th)]]),
    "rx":  lambda th: np.array([[np.cos(th / 2), -1j * np.sin(th / 2)], [-1j * np.sin(th / 2), np.cos(th / 2)]]),
}
# diagonal gates as a phase: phase(b) = g + c * (-1)^b, with c either constant or -coeff*theta/2
_DIAG_CONST = {"s": np.pi / 4, "sdg": -np.pi / 4, "z": np.pi / 2}     # phase = g - c(-1)^b  (s: pi/2 * b)
BLOCK = 4


def _arg(arg, theta):
    if arg is None:
        return None
    return arg[1] if arg[0] == "const" else arg[2] * theta[arg[1]]


_PHASE = ("cx", "rz", "s", "sdg", "z")


def _schedule(n, ops):
    """Reorder gates into alternating runs without changing the circuit: a gate may only move past
    gates acting on DIFFERENT qubits (its per-qubit predecessors always stay before it). Greedily
    takes every available phase-type gate (cx / diagonal), then every available single-qubit gate."""
    from collections import deque
    queues = [deque() for _ in range(n)]
    for i, (_, qs, _) in enumerate(ops):
        for q in qs:
            queues[q].append(i)
    done = [False] * len(ops)
    out, mode = [], "phase" if ops and ops[0][0] == "cx" else "1q"
    remaining = len(ops)
    while remaining:
        progressed = True
        taken_any = False
        while progressed:
            progressed = False
            for q in range(n):
                while queues[q]:
                    i = queues[q][0]
                    name, qs, _ = ops[i]
                    ok = (name in _PHASE) if mode == "phase" else (len(qs) == 1)
                    if not ok or any(queues[r][0] != i for r in qs):
                        break
                    for r in qs:
                        queues[r].popleft()
                    out.append(ops[i]); done[i] = True; remaining -= 1
                    progressed = taken_any = True
        mode = "1q" if mode == "phase" else "phase"
    return out


def compile_layers(n, ops):
    """Split the gate list into single-qubit and phase layers (done once per circuit)."""
    ops = [(o[0], list(o[1]), None if o[2] is None else tuple(o[2])) for o in ops]
    ops = _schedule(n, ops)
    layers, i = [], 0
    while i < len(ops):
        if ops[i][0] == "cx":                                      # phase layer
            masks = [1 << q for q in range(n)]
            terms = []                                             # (mask, kind, arg) ; kind: 'rz' or const
            while i < len(ops) and (ops[i][0] == "cx" or ops[i][0] in ("rz", "s", "sdg", "z")):
                name, qs, arg = ops[i]
                if name == "cx":
                    c, t = qs
                    masks[t] ^= masks[c]
                else:
                    terms.append((masks[qs[0]], name, arg))
                i += 1
            uniq = sorted({m for m, _, _ in terms})
            pos = {m: k for k, m in enumerate(uniq)}
            perm = None
            if masks != [1 << q for q in range(n)]:               # net CX permutation x -> y
                idx = np.arange(2 ** n)
                y = np.zeros(2 ** n, dtype=np.int64)
                for q in range(n):
                    par = np.zeros(2 ** n, dtype=np.int64)
                    m, b = masks[q], 0
                    while m:
                        if m & 1:
                            par ^= (idx >> b) & 1
                        m >>= 1; b += 1
                    y |= par << q
                perm = y
            layers.append(("phase", dict(masks=np.array(uniq, dtype=np.int64),
                                         terms=[(pos[m], name, arg) for m, name, arg in terms], perm=perm)))
        else:                                                      # single-qubit layer
            gates = []
            while i < len(ops) and ops[i][0] != "cx":
                gates.append(ops[i]); i += 1
            layers.append(("1q", gates))
    return n, layers


def _fwht(v):
    n = int(np.log2(len(v)))
    for q in range(n):
        w = v.reshape(-1, 2, 2 ** q)
        a, b = w[:, 0, :].copy(), w[:, 1, :]
        w[:, 0, :] += b
        w[:, 1, :] = a - b
    return v


def probs(compiled, theta):
    n, layers = compiled
    theta = np.asarray(theta, float)
    psi = np.zeros(2 ** n, complex); psi[0] = 1
    for kind, L in layers:
        if kind == "1q":
            U = [np.eye(2, dtype=complex) for _ in range(n)]
            for name, qs, arg in L:
                U[qs[0]] = _1Q[name](_arg(arg, theta)) @ U[qs[0]]
            lo = 0
            while lo < n:
                k = min(BLOCK, n - lo)
                M = np.array([[1.0 + 0j]])
                for q in range(lo + k - 1, lo - 1, -1):           # kron(U_{lo+k-1}, ..., U_lo)
                    M = np.kron(M, U[q])
                psi = np.matmul(M, psi.reshape(2 ** (n - lo - k), 2 ** k, 2 ** lo)).reshape(-1)
                lo += k
        else:
            c = np.zeros(2 ** n)
            for k, name, arg in L["terms"]:
                if name == "rz":
                    c[L["masks"][k]] += -0.5 * _arg(arg, theta)
                else:
                    c[L["masks"][k]] += -_DIAG_CONST[name]       # global phase dropped
            phase = _fwht(c)                                       # sum_k c_k (-1)^{popcount(x & m_k)}
            psi = psi * np.exp(1j * phase)
            if L["perm"] is not None:
                out = np.empty_like(psi); out[L["perm"]] = psi; psi = out
    return psi.real ** 2 + psi.imag ** 2


def verify(gl, trials=3, seed=0, x0=None):
    """Compare with ../fastsim.circuit_probs on random angles (and x0 if given)."""
    import os, sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    from fastsim import circuit_probs
    n = gl["num_qubits"]
    ops = [(o[0], o[1], None if o[2] is None else tuple(o[2])) for o in gl["ops"]]
    comp = compile_layers(n, gl["ops"])
    npar = len(gl["param_names"])
    rng = np.random.default_rng(seed)
    pts = [rng.uniform(-np.pi, np.pi, npar) for _ in range(trials)] + ([np.asarray(x0)] if x0 is not None else [])
    err = max(np.max(np.abs(probs(comp, t) - circuit_probs((n, ops), t))) for t in pts)
    assert err < 1e-10, err
    return err


if __name__ == "__main__":
    import os, sys, time
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    from spiq_artifacts import load_spiq
    from fastsim import compile_circuit, circuit_probs
    for P in ("P1", "P2", "P3"):
        pc, x0, _ = load_spiq(P)
        n, ops = compile_circuit(pc)
        gl = dict(num_qubits=n, ops=[[a, b, None if c is None else list(c)] for a, b, c in ops],
                  param_names=[p.name for p in pc.parameters])
        err = verify(gl, trials=2, x0=x0)
        comp = compile_layers(n, gl["ops"])
        t = time.time(); probs(comp, x0); t1 = time.time() - t
        t = time.time(); circuit_probs((n, ops), np.asarray(x0)); t2 = time.time() - t
        print(f"{P}: {n} qubits, {len(ops)} gates -> {len(comp[1])} layers; max |dp| vs fastsim = {err:.1e}; "
              f"{1000 * t1:.1f} ms vs {1000 * t2:.1f} ms per evaluation")
