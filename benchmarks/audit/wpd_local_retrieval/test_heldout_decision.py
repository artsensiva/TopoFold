"""SYNTHETIC tests for the amended held-out decision statistics (amendment 1). No mdCATH or AdK data.

Run: build/venv/bin/python -B -m pytest -p no:cacheprovider test_heldout_decision.py
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from heldout_decision import (BEST_ORDER, C1, EPS_G, analyse_cell, cell_statistics,  # noqa: E402
                              classify_loss, conservative_interval, overall_outcome, post_burn_frames,
                              select_best)


def make(base, n, others=None, c1=None, c1_cov=None):
    """SYNTHETIC cell: every method equals `base` unless overridden; coverage-matched = primary on covered domains."""
    nr = {m: np.array(base, float) for m in BEST_ORDER}
    nr.update({k: np.array(v, float) for k, v in (others or {}).items()})
    nr[C1] = np.array(c1 if c1 is not None else base, float)
    cov = np.ones(n, bool) if c1_cov is None else np.asarray(c1_cov, bool)
    nr_cm = {m: np.where(cov, nr[m], np.nan) for m in BEST_ORDER + (C1,)}
    return nr, nr_cm


# 1. tied medians with unequal paired differences -------------------------------------------------

def test_tied_medians_unequal_paired_differences_give_zero_loss():
    a = [0.1, 0.2, 0.6, 0.7, 0.8]
    b = [0.5, 0.55, 0.6, 0.9, 0.95]
    assert np.median(np.subtract(b, a)) == pytest.approx(0.2)      # old estimand would be 0.2 / 0.4 = 0.5
    nr, nr_cm = make([0.9] * 5, 5, others={'B2': a, 'M1': b})
    s = cell_statistics(nr, nr_cm, np.arange(5))
    assert s['best'] == 'B2' and s['L']['M1'] == 0.0 and s['L']['B2'] == 0.0


def test_loss_is_difference_of_medians_and_nonnegative():
    r = np.random.default_rng(1)
    nr, nr_cm = make(r.uniform(0.5, 0.9, 30), 30,
                     others={m: r.uniform(0.4, 0.95, 30) for m in BEST_ORDER})
    s = cell_statistics(nr, nr_cm, np.arange(30))
    med = {m: np.median(nr[m]) for m in BEST_ORDER}
    b = min(BEST_ORDER, key=lambda m: med[m])
    for m in BEST_ORDER:
        assert s['L'][m] == pytest.approx((med[m] - med[b]) / (1 - med[b]), abs=1e-15)
        assert s['L'][m] >= 0.0


# 2. identical methods -------------------------------------------------------------------------------

def test_identical_methods_zero_loss_in_every_resample():
    r = np.random.default_rng(2)
    nr, nr_cm = make(r.uniform(0.3, 0.7, 30), 30)
    res = analyse_cell(nr, nr_cm, n_boot=300, seed=5)
    assert res['observed']['best'] == BEST_ORDER[0]
    for m in BEST_ORDER:
        assert res['methods'][m]['L'] == 0.0 and res['methods'][m]['L_interval'] == [0.0, 0.0]
        assert res['methods'][m]['loss_class'] == 'within_margin'
    assert select_best({m: 0.5 for m in BEST_ORDER}) == BEST_ORDER[0]


# 3. zero, negative and negligible best gain in the observed sample --------------------------------

@pytest.mark.parametrize('level', [1.0, 1.2, 1.0 - EPS_G / 2])
def test_observed_nonpositive_or_negligible_best_gain_not_assessable(level):
    nr, nr_cm = make([level] * 30, 30, others={'M1': [level + 0.05] * 30})
    res = analyse_cell(nr, nr_cm, n_boot=200, seed=3)
    assert not res['observed']['loss_valid']
    for m in BEST_ORDER:
        assert res['methods'][m]['L'] is None
        assert res['methods'][m]['loss_class'] == 'not_assessable'
        assert res['methods'][m]['gain_class'] == 'uninformative'
        assert res['methods'][m]['invalid_fraction'] == 1.0
    assert overall_outcome([(res['methods']['M1']['loss_class'], res['methods']['M1']['gain_class'])] * 4) \
        == 'fails_criterion_in_all_cells'
    assert overall_outcome([('not_assessable', 'inconclusive')] * 4) == 'inconclusive'


# 4. invalid bootstrap denominators ------------------------------------------------------------------

def test_conservative_interval_places_invalid_at_both_tails():
    v = np.linspace(0.0, 0.05, 1000)
    inv = np.zeros(1000, bool)
    inv[:10] = True                                    # 1 % invalid: bounds stay finite but widen
    lo, hi, f = conservative_interval(v, inv)
    lo0, hi0, _ = conservative_interval(v, np.zeros(1000, bool))
    assert f == 0.01 and np.isfinite(lo) and np.isfinite(hi) and lo <= lo0 and hi >= hi0
    inv[:30] = True                                    # 3 % > alpha/2: both bounds infinite
    lo, hi, f = conservative_interval(v, inv)
    assert lo == -np.inf and hi == np.inf and f == 0.03
    assert classify_loss(True, lo, hi) == 'inconclusive'


def test_resamples_with_invalid_best_gain_are_counted_not_dropped():
    # SYNTHETIC: best method NR is 0.5 on 16 domains and 1.0 on 14; observed median 0.5 (valid),
    # but many resamples have median 1.0 (gain 0).
    best = np.array([0.5] * 16 + [1.0] * 14)
    nr, nr_cm = make(np.full(30, 1.0), 30, others={'B2': best, 'M1': best + 0.01})
    res = analyse_cell(nr, nr_cm, n_boot=1000, seed=7)
    assert res['observed']['loss_valid']
    f = res['methods']['M1']['invalid_fraction']
    assert 0.025 < f < 1.0
    assert res['methods']['M1']['L_interval'] == [-np.inf, np.inf]
    assert res['methods']['M1']['loss_class'] == 'inconclusive'


# 5. C1 coverage matching ----------------------------------------------------------------------------

def test_c1_uses_common_population_and_is_never_best():
    n = 30
    r = np.random.default_rng(4)
    b2 = r.uniform(0.5, 0.7, n)
    cov = np.array([True] * 26 + [False] * 4)
    nr, nr_cm = make(np.full(n, 0.9), n, others={'B2': b2}, c1=np.full(n, 0.1), c1_cov=cov)
    # primary B2 differs from its coverage-matched value on covered domains
    nr_cm['B2'] = np.where(cov, b2 + 0.1, np.nan)
    s = cell_statistics(nr, nr_cm, np.arange(n))
    assert s['best'] == 'B2'                                            # C1 lowest, but not a candidate
    mb = np.median(nr_cm['B2'][cov])
    mc = np.median(nr_cm[C1][cov])
    assert s['c1_n_domains'] == 26
    assert s['gain_best_cm'] == pytest.approx(1 - mb)
    assert s['L_C1'] == pytest.approx((mc - mb) / (1 - mb))           # numerator and denominator matched
    assert s['L_C1'] < 0 and s['G_C1'] == pytest.approx(0.9)


def test_c1_insufficient_coverage_not_assessable_and_mismatched_population_rejected():
    n = 30
    cov = np.array([True] * 20 + [False] * 10)                         # 20 < 24 domains
    nr, nr_cm = make(np.full(n, 0.6), n, c1_cov=cov)
    res = analyse_cell(nr, nr_cm, n_boot=100, seed=1)
    assert res['methods'][C1]['loss_class'] == 'not_assessable'
    assert res['methods'][C1]['gain_class'] == 'not_assessable'
    nr_cm['M1'] = nr['M1'].copy()                                      # breaks the common population
    with pytest.raises(ValueError):
        analyse_cell(nr, nr_cm, n_boot=10, seed=1)


# 6. frame-count arithmetic (amendment §5) -------------------------------------------------------------

def test_post_burn_frames_minimum_and_k_fraction():
    assert post_burn_frames(200) == 180
    vals = [post_burn_frames(n) for n in range(200, 2001)]
    assert min(vals) == 180 and all(b >= a for a, b in zip(vals, vals[1:]))
    assert 10 / 180 == pytest.approx(0.0556, abs=1e-4)


def test_overall_outcome_precedence():
    W, I, B, U = 'within_margin', 'informative', 'beyond_margin', 'inconclusive'
    assert overall_outcome([(W, I)] * 4) == 'meets_criterion_in_all_cells'
    assert overall_outcome([(B, I)] * 4) == 'fails_criterion_in_all_cells'
    assert overall_outcome([(W, I), (W, I), (B, I), (B, I)]) == 'cell_dependent'
    assert overall_outcome([(W, I), (W, I), (W, I), (U, I)]) == 'inconclusive'
    assert overall_outcome([(W, U)] * 4) == 'inconclusive'
