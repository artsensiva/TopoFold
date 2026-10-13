"""SYNTHETIC regression tests for the amended LE-1 analysis (Amendment 1 §1, §4, §5). No AdK data.

Run: ../wpd_local_retrieval/build/venv/bin/python -B -m pytest -p no:cacheprovider test_analyse_a1.py
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyse_dev as AD0  # noqa: E402  (original, preserved)
import analyse_dev_a1 as AD  # noqa: E402


def unequal():
    """Windows with unequal valid query counts, including one fully undefined window."""
    r = np.random.default_rng(0)
    X = r.random((45, 310))
    X[0, 50:] = np.nan                  # window 0: 50 valid queries, high values
    X[0, :50] += 5.0
    X[1, :] = np.nan                    # window 1: undefined
    X[2:10, ::3] = np.nan               # partial missingness
    return X


def test_identity_resample_equals_point_estimate_with_unequal_counts():
    X = unequal()
    theta = AD.window_mean_stat(X)
    for qb in (100, 250):
        S, C = AD.window_block_sums(X, qb)
        one_q = np.ones((1, S.shape[1]))
        one_w = np.ones((1, X.shape[0]))
        assert AD.boot_stat(S, C, one_q, one_w)[0] == pytest.approx(theta, abs=1e-12)


def test_original_pooled_bootstrap_estimand_differs_here():
    X = unequal()
    S, C = AD0.block_sums(X, qb=100, wb=20)
    pooled = S.sum() / C.sum()
    assert abs(pooled - AD.window_mean_stat(X)) > 0.05     # the defect the amendment removes


def test_bootstrap_centred_on_point_estimate_and_paired():
    X = unequal()
    Y = X + 0.1
    CQ, CW = AD.resample_weights(*X.shape, 100, 500, 3)
    bx = AD.boot_stat(*AD.window_block_sums(X, 100), CQ, CW)
    by = AD.boot_stat(*AD.window_block_sums(Y, 100), CQ, CW)
    assert np.allclose(by - bx, 0.1)                        # constant shift exact in every resample
    assert abs(np.median(bx) - AD.window_mean_stat(X)) < 0.2


def test_paired_delta_on_units_vs_difference_of_thetas_with_different_missingness():
    X = unequal()
    Y = X.copy()
    Y[3, :] = np.nan
    d = X - Y
    assert AD.window_mean_stat(d) == pytest.approx(0.0, abs=1e-12)
    assert AD.window_mean_stat(X) - AD.window_mean_stat(Y) != pytest.approx(0.0, abs=1e-6)


def test_governing_rule_and_outcomes():
    assert AD.governing(40.0) == 100 and AD.governing(180.0) == 250
    assert AD.governing(300.0) is None and AD.governing(None) is None
    assert AD.outcome(['within_margin'] * 4) == 'retain_as_heldout_candidate'
    assert AD.outcome(['within_margin', 'beyond_margin', 'within_margin', 'inconclusive']) == 'restrict'
    assert AD.outcome(['not_assessable'] + ['within_margin'] * 3) == 'insufficient_evidence'
    (lo, hi), f = AD.interval(np.r_[np.full(30, np.nan), np.linspace(0, 1, 970)])
    assert lo == -np.inf and hi == np.inf and f == pytest.approx(0.03)
