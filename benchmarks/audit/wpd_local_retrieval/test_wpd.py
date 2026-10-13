"""WP-D harness tests on SYNTHETIC fixtures (protocol §10 A3). No AdK data used.

Run: build/venv/bin/python -B -m pytest -p no:cacheprovider test_wpd.py
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from wpd_lib import (circular_mean_R, d_chord, d_circ, d_scalar, dihedral, horn_rmsd,  # noqa: E402
                     kabsch_rmsd, membership, query_metrics, sin_angle, spearman)
from run_adk import split_ranges  # noqa: E402


def rng(s=0):
    return np.random.default_rng(s)


# --- metrics on SYNTHETIC distance matrices -----------------------------------------------

def test_identical_rankings():
    T = rng(1).random((20, 50))
    m = query_metrics(T.copy(), T, 10)
    assert np.allclose(m['overlap'], 1.0) and np.all(m['regret'] == 0) and np.allclose(m['norm_regret'], 0)


def test_complete_ties_equal_random_expectation():
    T = rng(2).random((20, 50))
    m = query_metrics(np.zeros((20, 50)), T, 10)
    assert np.allclose(m['overlap'], 10 / 50)
    assert np.allclose(m['norm_regret'], 1.0) and m['ties'].all()


def test_boundary_tie_hand_computed():
    # one query, 5 database items; method distances: [1, 2, 2, 2, 9]; target best two: items 0, 3
    R = np.array([[1.0, 2.0, 2.0, 2.0, 9.0]])
    T = np.array([[0.0, 5.0, 5.0, 1.0, 9.0]])
    P, ok, tie = membership(R, 2)
    assert np.allclose(P, [[1, 1 / 3, 1 / 3, 1 / 3, 0]]) and tie[0]
    m = query_metrics(R, T, 2)
    # overlap = (1·1 + 1/3·1) / 2 = 2/3 ; retrieved mean T = (0 + (5+5+1)/3)/2 ; oracle = (0+1)/2
    assert abs(m['overlap'][0] - 2 / 3) < 1e-15
    assert abs(m['regret'][0] - ((0 + 11 / 3) / 2 - 0.5)) < 1e-12


def test_random_retrieval_near_reference():
    r = rng(3)
    T = r.random((400, 300))
    m = query_metrics(r.random((400, 300)), T, 10)
    assert abs(m['overlap'].mean() - 10 / 300) < 0.01
    assert abs(np.nanmean(m['norm_regret']) - 1.0) < 0.05


def test_constant_target_normalized_regret_undefined():
    m = query_metrics(rng(4).random((10, 30)), np.full((10, 30), 0.7), 5)
    assert not m['norm_defined'].any() and np.all(np.isnan(m['norm_regret']))
    assert np.allclose(m['overlap'], 5 / 30)


def test_reversed_retrieval_is_worst():
    T = rng(5).random((30, 100))
    m = query_metrics(-T, T, 10)
    assert np.allclose(m['overlap'], 0.0) and np.all(m['norm_regret'] > 1.0)


def test_ill_conditioned_database_ranked_last_and_undefined_query():
    mq = np.array([0.1, np.nan])
    md = np.array([0.1, np.nan, 0.2, 0.3])
    D = d_scalar(mq, md)
    assert np.isinf(D[0, 1]) and np.all(np.isinf(D[1]))
    P, ok, _ = membership(D, 2)
    assert ok[0] and not ok[1] and P[0, 1] == 0.0


# --- geometry references ---------------------------------------------------------------------

def test_dihedral_known_values():
    p0, p1, p2 = np.array([1.0, 0, 0]), np.zeros(3), np.array([0, 0, 1.0])
    for ang in (-2.5, -0.3, 0.0, 1.2, math.pi - 1e-9):
        p3 = np.array([math.cos(ang), math.sin(ang), 1.0])
        assert abs(math.remainder(float(dihedral(p0, p1, p2, p3)) - ang, 2 * math.pi)) < 1e-12


def test_sin_angle_detects_collinear():
    assert sin_angle(np.array([0, 0, 0.0]), np.array([1, 0, 0.0]), np.array([2, 1e-6, 0.0])) < 1e-3
    assert sin_angle(np.array([0, 0, 0.0]), np.array([1, 0, 0.0]), np.array([1, 1, 0.0])) > 0.99


def test_horn_matches_kabsch_and_rejects_reflection():
    r = rng(6)
    A, B = r.normal(size=(5, 8, 3)) * 4, r.normal(size=(7, 8, 3)) * 4
    H = horn_rmsd(A, B)
    for i in range(5):
        for j in range(7):
            assert abs(H[i, j] - kabsch_rmsd(A[i], B[j])) < 1e-6
    mirror = A[:1] * np.array([1, 1, -1])
    assert horn_rmsd(A[:1], mirror)[0, 0] > 0.1


def test_chord_and_circ_relation():
    r = rng(7)
    a, b = r.uniform(-np.pi, np.pi, (20, 5)), r.uniform(-np.pi, np.pi, (30, 5))
    ch, ci = d_chord(a, b), d_circ(a, b)
    assert np.all(ch <= ci + 1e-12) and np.all(ch >= 2 / np.pi * ci - 1e-12)


def test_circular_mean_resultant():
    m, R = circular_mean_R(np.array([[0.1, np.pi + 0.1]]))
    assert R[0] < 1e-12


def test_spearman_matches_scipy_with_ties():
    stats = pytest.importorskip('scipy.stats')
    r = rng(8)
    a = np.round(r.random(500), 2)
    b = a + r.normal(scale=0.3, size=500)
    assert abs(spearman(a, b) - stats.spearmanr(a, b).statistic) < 1e-12


def test_split_ranges_exact_and_nonoverlapping():
    assert split_ranges(2253) == ((0, 1013), (1013, 1238), (1238, 2253))
    assert split_ranges(1678) == ((0, 755), (755, 922), (922, 1678))


def _helix(n, rise=1.5, r=2.3, step=1.745):
    t = np.arange(n) * step
    return np.stack([r * np.cos(t), r * np.sin(t), rise * np.arange(n)], axis=1)


def test_binding_window_indexing_matches_numpy_reference():
    tf = pytest.importorskip('topofold')
    from wpd_lib import turning_angle
    X = _helix(12) + rng(9).normal(scale=0.2, size=(12, 3))
    X = X.astype(np.float32).astype(np.float64)
    kap, tau, _ = tf.compute_invariants(X.astype(np.float32))
    W, s = 5, 4
    ref_tau = dihedral(X[s:s + W - 3], X[s + 1:s + W - 2], X[s + 2:s + W - 1], X[s + 3:s + W])
    ref_kap = turning_angle(X[s:s + W - 2], X[s + 1:s + W - 1], X[s + 2:s + W])
    assert np.max(np.abs(np.asarray(tau)[s:s + W - 3] - ref_tau)) < 1e-9
    assert np.max(np.abs(np.asarray(kap)[s:s + W - 2] - ref_kap)) < 1e-9


def test_degenerate_window_excluded_for_all_methods():
    from run_adk import window_validity
    F, N = 3, 10
    sinCA = np.full((F, N), 0.8)
    sinCA[1, 4] = 1e-4                       # frame 1: collinear Cα triple centred on residue index 4
    ones = np.full((F, N), 0.9)
    geo = {'sinCA': sinCA, 's_cna': ones, 's_nac': ones, 's_acn': ones,
           'b_nca': np.full((F, N), 1.46), 'b_cac': np.full((F, N), 1.52), 'b_cn': np.full((F, N - 1), 1.33),
           'binfail': np.array([False, False, True])}
    valid, why = window_validity(geo, 2, 5)  # interior indices 3..5 include 4
    assert valid.tolist() == [True, False, False]
    assert why['degenerate_ca'].tolist() == [False, True, False] and why['binding_failure'][2]
    valid, _ = window_validity(geo, 5, 5)    # interior 6..8: unaffected except binding failure
    assert valid.tolist() == [True, True, False]
