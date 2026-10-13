"""SYNTHETIC regression tests for BQ-1 (numerical correctness, tie policy, allocation lifetimes).
Run: ../wpd_local_retrieval/build/venv/bin/python -B -m pytest -p no:cacheprovider test_b2opt.py
"""
import gc
import sys
import tracemalloc
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import b2opt_lib as B  # noqa: E402

RNG = np.random.default_rng(7)
ARMS = ('V0', 'V1', 'V2b')                                   # V2 (frozen) is tested separately below


def chain(n, L=20, step=3.8):
    """SYNTHETIC random-walk Cα chains (n, L, 3), Å."""
    v = RNG.normal(size=(n, L, 3))
    v *= step / np.linalg.norm(v, axis=-1, keepdims=True)
    return np.cumsum(v, axis=1)


def arm_dist(arm, X, q, di, s):
    if arm == 'V0':
        return B.v0_window(X, q, di, s)
    P, G = B.prepare(X, di, [s])
    return B.prepared_window(arm, X[q, s:s + B.W], P[0], G[0])[0]


def rand_rot():
    Q, _ = np.linalg.qr(RNG.normal(size=(3, 3)))
    return Q * np.sign(np.linalg.det(Q))


# --- tie policy and ranking --------------------------------------------------------------------------------------------

def test_topk_lex_matches_bruteforce_with_exact_ties():
    for _ in range(200):
        d = RNG.integers(0, 6, 40).astype(float)                # many exact ties
        ref = np.lexsort((np.arange(40), d))[:10]
        assert np.array_equal(B.topk_lex(d), ref)
    assert np.array_equal(B.topk_lex(np.array([2.0, 1.0, 1.0])), [1, 2, 0])


def test_rank_check_tolerance_and_strictness():
    d = np.arange(20, dtype=float)
    d[10] = d[9] + 1e-6                                         # near tie at the boundary
    ok = B.topk_lex(d)
    swapped = ok.copy()
    swapped[9] = 10
    assert B.rank_check(ok, d) == (True, True, True)
    assert B.rank_check(swapped, d)[:2] == (True, False)        # within τ_rank: equivalent, not strict
    far = ok.copy()
    far[9] = 15
    assert B.rank_check(far, d)[0] is False
    dup = ok.copy()
    dup[9] = dup[8]
    assert B.rank_check(dup, d)[0] is False


# --- polynomial algebra ------------------------------------------------------------------------------------------------

def test_charpoly_and_det4_against_numpy():
    M = RNG.normal(size=(200, 3, 3)) * 30
    Km = B.key_matrix(M)
    assert np.allclose(np.trace(Km, axis1=1, axis2=2), 0, atol=1e-9)
    c2, c1, c0 = B.charpoly(M)
    for i in range(200):
        p = np.poly(Km[i])
        scale = np.abs(p).max()
        assert np.allclose([p[1], p[2], p[3], p[4]], [0, c2[i], c1[i], c0[i]], atol=1e-9 * scale)
    assert np.allclose(B.det4(Km), np.linalg.det(Km), rtol=1e-9, atol=1e-6)


def test_qcp_lambda_equals_eigvalsh_and_forced_fallback():
    X = chain(50)
    a = X[0, :8] - X[0, :8].mean(0)
    b = X[1:, :8] - X[1:, :8].mean(1, keepdims=True)
    M = np.einsum('wi,dwj->dij', a, b)
    E0 = 0.5 * (np.sum(a * a) + np.sum(b * b, axis=(1, 2)))
    lam, bad = B.lam_max_qcp(M, E0)
    ref = np.linalg.eigvalsh(B.key_matrix(M))[:, -1]
    assert not bad.any() and np.allclose(lam, ref, rtol=1e-12, atol=1e-9)
    lam0, bad0 = B.lam_max_qcp(M, E0, max_iter=0)               # no iterations → every row falls back
    assert bad0.all() and np.array_equal(lam0, ref)


# --- distances against independent references ---------------------------------------------------------------------------

@pytest.mark.parametrize('arm', ARMS)
def test_random_windows_match_three_references(arm):
    X = chain(60)
    di = np.arange(1, 60)
    for s in (0, 5, 12):
        d = arm_dist(arm, X, 0, di, s)
        r = B.ref_svd(X[0, s:s + 8], X[di, s:s + 8])
        assert np.max(np.abs(d - r)) <= B.TOL_D
        for j in (0, 17, 58):
            a, b = X[0, s:s + 8], X[di[j], s:s + 8]
            assert abs(B.ref_kabsch_pair(a, b) - r[j]) <= B.TOL_D
            assert abs(B.ref_qcprot_pair(a, b) - r[j]) <= B.TOL_D


def degenerate_db():
    base = chain(1)[0]
    t = np.linspace(0, 20, 20)[:, None]
    line = t * np.array([[1.0, 2.0, 2.0]]) / 3.0
    plane = np.c_[RNG.normal(size=(20, 2)) * 5, np.zeros(20)]
    mirror = base * np.array([1.0, 1.0, -1.0])
    return np.stack([base, line, plane, mirror, base @ rand_rot().T + 7.0, base.copy(), base.copy()])


@pytest.mark.parametrize('arm', ARMS)
def test_degenerate_identical_mirror_and_exact_ties(arm):
    X = degenerate_db()
    for q in range(4):                                          # base, collinear, planar, mirror queries
        di = np.arange(len(X))
        for s in (0, 6):
            d = arm_dist(arm, X, q, di, s)
            r = B.ref_svd(X[q, s:s + 8], X[di, s:s + 8])
            assert np.max(np.abs(d - r)) <= B.TOL_D
    d = arm_dist(arm, X, 0, np.arange(len(X)), 0)
    assert d[0] <= B.TOL_D and d[4] <= 1e-5 and d[5] <= B.TOL_D
    assert d[5] == d[6] and d[0] == d[5]                        # exact duplicates: bitwise-equal distances
    assert np.array_equal(B.topk_lex(d, 4), np.lexsort((np.arange(len(d)), d))[:4])
    assert d[3] > 0.1                                           # mirror image is not superposable by a proper rotation


@pytest.mark.parametrize('arm', ARMS)
def test_permutation_rigid_transform_and_near_ties(arm):
    X = chain(80)
    di = np.arange(1, 80)
    d = arm_dist(arm, X, 0, di, 3)
    perm = RNG.permutation(len(di))
    dp = arm_dist(arm, X, 0, di[perm], 3)
    assert np.max(np.abs(dp - d[perm])) <= 1e-9
    assert B.rank_check(perm[B.topk_lex(dp)], d)[0]
    Xr = X.copy()
    Xr[0] = X[0] @ rand_rot().T + np.array([3.0, -40.0, 11.0])
    assert np.max(np.abs(arm_dist(arm, Xr, 0, di, 3) - d)) <= 1e-8
    Xn = np.concatenate([X, X[1:11] + 1e-9 * RNG.normal(size=(10, 20, 3))])   # near-tie twins
    dn = arm_dist(arm, Xn, 0, np.arange(1, len(Xn)), 3)
    rn = B.ref_svd(Xn[0, 3:11], Xn[1:, 3:11])
    assert B.rank_check(B.topk_lex(dn), rn)[0]


# --- allocation lifetimes ----------------------------------------------------------------------------------------------

def test_prepared_arrays_owned_exact_and_query_does_not_copy_them():
    X = chain(1501, L=60)
    di = np.arange(1, 1501)
    peaks = {}
    for nwin in (5, 40):
        wins = list(range(nwin))
        P, G = B.prepare(X, di, wins)
        assert P.flags['OWNDATA'] and G.flags['OWNDATA'] and P.flags['C_CONTIGUOUS']
        assert P.nbytes + G.nbytes == nwin * len(di) * (B.W * 3 + 1) * 8
        for arm in ('V1', 'V2'):
            B.search(arm, X, 0, di, wins, (P, G))
            gc.collect()
            tracemalloc.start()
            out, _, _ = B.search(arm, X, 0, di, wins, (P, G))
            peaks[(arm, nwin)] = tracemalloc.get_traced_memory()[1] - out.nbytes
            tracemalloc.stop()
            del out
        del P, G
    for arm in ('V1', 'V2'):
        # per-window transients only: the query peak does not grow with the number of prepared windows
        assert peaks[(arm, 40)] < 1.2 * peaks[(arm, 5)] + 64 * 1024
        assert peaks[(arm, 40)] < 40 * len(di) * (B.W * 3 + 1) * 8 / 4


def test_repeated_queries_release_outputs():
    X = chain(301, L=60)
    di = np.arange(1, 301)
    wins = list(range(52))
    P, G = B.prepare(X, di, wins)
    tracemalloc.start()
    out, _, _ = B.search('V2', X, 0, di, wins, (P, G))
    one = tracemalloc.get_traced_memory()[1]
    tracemalloc.reset_peak()
    for q in range(3):
        del out
        gc.collect()
        out, _, _ = B.search('V2', X, q, di, wins, (P, G), keep_dist=True)
    many = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert many < 1.3 * one + 52 * 300 * 8 * 1.1                 # one output (+ its distance array) alive at a time


def test_frozen_v2_precision_floor_on_identical_collinear_windows_and_v2b_fallback():
    """Documents the frozen-V2 defect found before measurement (Amendment 1): for identical collinear windows λ_max
    is a multiple root and Newton on P(λ) stops at a √ε-level error, beyond the 1e-6 Å tolerance; V2b falls back."""
    X = degenerate_db()
    di = np.arange(len(X))
    err2, nfb_total = 0.0, 0
    for s in range(0, 13):
        r = B.ref_svd(X[1, s:s + 8], X[di, s:s + 8])
        err2 = max(err2, abs(arm_dist('V2', X, 1, di, s)[1] - r[1]))
        P, G = B.prepare(X, di, [s])
        d2b, nfb = B.prepared_window('V2b', X[1, s:s + 8], P[0], G[0])
        assert np.max(np.abs(d2b - r)) <= B.TOL_D
        nfb_total += nfb
    assert err2 > B.TOL_D and nfb_total >= 1


def test_v2b_no_fallback_on_generic_chains():
    X = chain(200)
    P, G = B.prepare(X, np.arange(1, 200), [4])
    _, nfb = B.prepared_window('V2b', X[0, 4:12], P[0], G[0])
    assert nfb == 0
