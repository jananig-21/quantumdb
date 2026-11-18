"""Speed patch for clapton (the team's SPIQ search), producing IDENTICAL stim circuits.

clapton's ParametrizedCliffordCircuit.stim_circuit() builds the circuit with one
stim.Circuit.append() call per gate. With stim 1.16 on this machine each append costs
~60-600 us, so building a 180-gate circuit takes ~0.14 s and dominates the genetic search
(~99.9 % of the loss-function time; the tableau simulation itself takes ~1 ms). For the
noiseless case used here (no depolarization model, no Pauli measurement - the team's runs
use err=None), this patch writes the same instructions as text and parses them once.

`verify(pcirc, trials)` asserts str(fast) == str(original) for random parameter assignments.
Nothing in the team's files is modified; the patch is applied at runtime by spiq_many.py.
"""
import numpy as np
import stim

_original = None


def _fast_stim_circuit(self, pauli=None):
    if pauli is not None or self.depolarization_model is not None:
        return _original(self, pauli)
    lines = []
    for gate in self.gates:
        gate_id = gate.get_stim_id()
        if gate_id is None:
            continue
        ids = gate_id if isinstance(gate_id, (list, tuple)) else [gate_id]
        qbs = " ".join(str(int(q)) for q in gate.qbs)
        for gid in ids:
            lines.append(f"{gid} {qbs}")
    return stim.Circuit("\n".join(lines))


def apply():
    global _original
    from clapton.clifford import ParametrizedCliffordCircuit
    if _original is None:
        _original = ParametrizedCliffordCircuit.stim_circuit
        ParametrizedCliffordCircuit.stim_circuit = _fast_stim_circuit


def verify(pcirc, trials=50, seed=0):
    """Fast and original builders must give the same circuit for random Clifford assignments."""
    rng = np.random.default_rng(seed)
    space = pcirc.parameter_space()
    for _ in range(trials):
        x = [int(rng.choice(s)) for s in space]
        pcirc.assign(x)
        a, b = str(_original(pcirc)), str(_fast_stim_circuit(pcirc))
        assert a == b, "fast stim circuit differs from clapton's"
    return trials


# ---- second patch: no per-generation Manager process when a GA runs on a single core -----------
# clapton.eval_xs_terms_mp creates mp.Manager().Queue() on EVERY generation, which starts a
# server process (~0.1 s) even when n_proc == 1 and no worker process will ever use the queue.
# With n_proc == 1 the queue stays empty, so an in-process queue gives identical losses.
class _LocalManager:
    def Queue(self):
        import queue
        return queue.Queue()


class _MpShim:
    """Stands in for the `multiprocessing` module inside clapton.clapton during a single-core call."""
    def __init__(self, real):
        self._real = real

    def Manager(self):
        return _LocalManager()

    def __getattr__(self, name):
        return getattr(self._real, name)


def apply_single_core_eval():
    import clapton.clapton as C
    if getattr(C, "_fast_eval_patched", False):
        return
    original = C.eval_xs_terms_mp

    def eval_xs_terms_mp(xs, paulis, coeffs, vqe_pcirc, trans_pcirc=None, n_proc=1, *args, **kw):
        if n_proc != 1:
            return original(xs, paulis, coeffs, vqe_pcirc, trans_pcirc, n_proc, *args, **kw)
        real = C.mp
        C.mp = _MpShim(real)
        try:
            return original(xs, paulis, coeffs, vqe_pcirc, trans_pcirc, n_proc, *args, **kw)
        finally:
            C.mp = real

    C.eval_xs_terms_mp = eval_xs_terms_mp
    C._fast_eval_patched = True
