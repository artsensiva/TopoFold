"""SYNTHETIC tests for the LE-1 block-bootstrap analysis helpers (no AdK data).

Run: ../wpd_local_retrieval/build/venv/bin/python -B -m pytest -p no:cacheprovider test_analyse.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyse_dev as AD  # noqa: E402


def test_block_sums_match_direct_sums_with_nan_and_partial_blocks():
    r = np.random.default_rng(0)
    X = r.random((47, 233))
    X[3, 5] = np.nan
    S, C = AD.block_sums(X, qb=100, wb=20)
    assert S.shape == (3, 3) and C.sum() == X.size - 1
    assert abs(S.sum() - np.nansum(X)) < 1e-9
    assert abs(S[2, 2] - np.nansum(X[40:, 200:])) < 1e-12


def test_bootstrap_paired_identical_methods_have_zero_difference_and_unit_weights_recover_mean():
    r = np.random.default_rng(1)
    base = r.random((60, 300))
    arr = {m: base.copy() for m in AD.METHODS}
    arr['A'] = base + 0.05
    AD.N_BOOT = 200
    b = AD.bootstrap_cell(arr, 100, 7)
    assert np.allclose(b['M1'] - b['M2'], 0.0)
    assert np.allclose(b['A'] - b['M1'], 0.05)            # constant shift survives every resample exactly
    S, C = AD.block_sums(base, 100, 20)
    assert abs(S.sum() / C.sum() - base.mean()) < 1e-12


def test_interval_conservative_invalid_handling_and_classes():
    v = np.linspace(0, 0.05, 1000)
    inv = np.zeros(1000, bool)
    inv[:30] = True
    assert AD.interval(v, inv) == [-np.inf, np.inf]
    assert AD.classify(0.0, 0.05, True) == 'within_margin'
    assert AD.classify(0.2, 0.3, True) == 'beyond_margin'
    assert AD.classify(0.05, 0.2, True) == 'inconclusive'
    assert AD.classify(0.0, 0.05, False) == 'not_assessable'
