"""Faster drop-in for ../fastsim.qaoa_probs (same state, same parameter order), plus shot sampling.

fastsim applies the mixer RX(2*beta) one qubit at a time: ~10 full passes over the 2^n state per
qubit, which on 20 qubits is memory-bandwidth bound (~0.2 s alone, ~1.6 s with 8 processes sharing
memory). Every qubit gets the SAME 2x2 rotation, so a block of k qubits gets RX^{(x)k}: applying it as
one 2^k x 2^k matrix to blocks of 4 qubits needs ~5 passes per layer instead of ~200 (4-5x faster). verify() checks the
result against fastsim.qaoa_probs to 1e-10.

sample_energy draws the 10,240 shots by inverse-CDF sampling (cumulative sum + searchsorted): the
same multinomial distribution a noiseless QasmSimulator samples from, without numpy's O(2^n)
sequential binomial draws.
"""
import functools
import numpy as np

BLOCK = 4


@functools.lru_cache(maxsize=64)
def _rx_block(beta, k):
    c, s = np.cos(beta), -1j * np.sin(beta)
    m = np.array([[c, s], [s, c]])
    out = np.array([[1.0 + 0j]])
    for _ in range(k):
        out = np.kron(out, m)
    return out


def qaoa_probs(ising_diag, params, reps):
    n = int(np.log2(len(ising_diag)))
    betas, gammas = params[:reps], params[reps:]
    psi = np.full(2 ** n, 2 ** (-n / 2), dtype=complex)
    for l in range(reps):
        psi *= np.exp(-1j * gammas[l] * ising_diag)
        lo = 0
        while lo < n:
            k = min(BLOCK, n - lo)
            v = psi.reshape(2 ** (n - lo - k), 2 ** k, 2 ** lo)
            psi = np.matmul(_rx_block(float(betas[l]), k), v).reshape(-1)
            lo += k
    p = psi.real ** 2 + psi.imag ** 2
    return p


def sample_energy(probs, energies, shots, rng):
    cdf = np.cumsum(probs)
    cdf /= cdf[-1]
    idx = np.searchsorted(cdf, rng.random(shots), side="right")
    return float(energies[np.minimum(idx, len(cdf) - 1)].mean())


def verify():
    import os, sys, time
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
    from fastsim import qaoa_probs as ref
    rng = np.random.default_rng(1)
    for n in (5, 6, 7, 12, 14, 16, 20):
        d = rng.normal(size=2 ** n) * 10
        for _ in range(3):
            th = rng.uniform(-np.pi, np.pi, 4)
            err = np.max(np.abs(qaoa_probs(d, th, 2) - ref(d, th, 2)))
            assert err < 1e-10, (n, err)
        t = time.time()
        for _ in range(5):
            p = qaoa_probs(d, th, 2); sample_energy(p, d, 10240, rng)
        print(f"n={n:2d}: matches fastsim (max |dp| < 1e-10); {1000 * (time.time() - t) / 5:.1f} ms per evaluation")
    # sampling: mean of many shot estimates == exact expectation
    d = rng.normal(size=2 ** 6); p = qaoa_probs(d, rng.uniform(-1, 1, 4), 2)
    est = np.mean([sample_energy(p, d, 10240, rng) for _ in range(2000)])
    assert abs(est - p @ d) < 5 * np.sqrt(p @ d ** 2 - (p @ d) ** 2) / np.sqrt(10240 * 2000), (est, p @ d)
    print("sample_energy: unbiased (2,000 x 10,240-shot estimates match the exact expectation)")


if __name__ == "__main__":
    verify()
