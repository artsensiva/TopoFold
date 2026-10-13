"""SYNTHETIC tests for held-out dependence diagnostics and the Amendment-2 candidate set (no mdCATH/AdK data).

Run: build/venv/bin/python -B -m pytest -p no:cacheprovider test_heldout_dependence.py
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import heldout_decision as HD  # noqa: E402
import heldout_decision_a2 as HD2  # noqa: E402
from heldout_dependence import (classify_series, comparator_for_dependence,  # noqa: E402
                                is_structural_identity)


def ar1(phi, n, seed):
    r = np.random.default_rng(seed)
    x = np.empty(n)
    x[0] = r.normal()
    for i in range(1, n):
        x[i] = phi * x[i - 1] + r.normal()
    return x


def test_structural_identity_self_comparison_and_violation():
    assert is_structural_identity('A', 'A', 8) and is_structural_identity('F', 'F_lock', 5)
    assert not is_structural_identity('F', 'F_lock', 8) and not is_structural_identity('M1', 'M2', 5)
    assert classify_series(np.zeros(800), identity=True)['class'] == 'structural_identity'
    with pytest.raises(ValueError):
        classify_series(np.r_[np.zeros(799), 1e-6], identity=True)


def test_constant_observed_is_not_identity_and_not_failure():
    c = classify_series(np.full(800, 0.031))
    assert c['class'] == 'constant_observed' and c['value'] == pytest.approx(0.031)
    assert classify_series(np.zeros(800))['class'] == 'constant_observed'   # zero but not declared


def test_estimation_failures_are_distinguished():
    assert classify_series(np.r_[np.full(30, np.nan), np.arange(20.0)])['reason'] == 'too_short_or_nonfinite'
    slow = classify_series(ar1(0.999, 400, 1))
    assert slow['class'] == 'estimation_failure' and slow['reason'] == 'tau_not_resolvable_from_series_length'


def test_estimated_tau_matches_ar1_theory():
    phi = 0.8                                   # τ_int = (1 + φ)/(1 − φ) = 9
    taus = [classify_series(ar1(phi, 20000, s))['tau_int_frames'] for s in range(5)]
    assert abs(np.mean(taus) - 9.0) < 1.5
    assert classify_series(ar1(0.0, 2000, 3))['tau_int_frames'] < 1.5


def test_runner_up_used_when_method_is_best():
    assert comparator_for_dependence('A', best='A', runner_up='B2') == 'B2'
    assert comparator_for_dependence('M1', best='A', runner_up='B2') == 'A'


def test_candidate_set_adds_A_without_changing_original_module():
    assert HD.BEST_ORDER == ('B2', 'B1', 'B3', 'M1', 'M2')
    assert HD2.BEST_ORDER == ('B2', 'B1', 'B3', 'M1', 'M2', 'A')
    med = {'B2': 0.70, 'B1': 0.72, 'B3': 0.75, 'M1': 0.80, 'M2': 0.80, 'A': 0.65}
    assert HD2.select_best(med) == 'A' and HD.select_best(med) == 'B2'


def test_a2_estimand_identical_to_amendment1_when_A_not_best():
    r = np.random.default_rng(4)
    n = 30
    nr = {m: r.uniform(0.5, 0.9, n) for m in ('B2', 'B1', 'B3', 'M1', 'M2', 'C1')}
    nr['A'] = np.full(n, 0.99)
    nr_cm = {m: nr[m].copy() for m in nr}
    s1 = HD.cell_statistics(nr, nr_cm, np.arange(n))
    s2 = HD2.cell_statistics(nr, nr_cm, np.arange(n))
    assert s1['best'] == s2['best']
    for m in HD.BEST_ORDER:
        assert s1['L'][m] == s2['L'][m]
    assert s2['L']['A'] == pytest.approx((np.median(nr['A']) - np.median(nr[s2['best']])) / s2['gain_best'])
