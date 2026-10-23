"""Tests for topofold_inputs.compare (histogram JS distance, binning, windows, blocks). SYNTHETIC inputs only.

NumPy only; the SciPy cross-check runs where SciPy is installed (e.g. the WP-D audit venv).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from topofold_inputs.compare import (  # noqa: E402
    circular_hist, compare_channel, equal_duration_blocks, equalize, js_distance, linear_hist, spearman,
    window_distances,
)


def test_js_distance_analytical_cases():
    p = np.array([3, 1, 0, 6])
    assert js_distance(p, p) == 0.0
    assert js_distance(p, 2 * p) == 0.0                      # counts are normalised
    assert js_distance([1, 0], [0, 1]) == pytest.approx(1.0)  # disjoint support: maximum in base 2
    assert js_distance([5, 0, 0], [0, 2, 7]) == pytest.approx(1.0)
    q = np.array([1, 4, 2, 0])
    assert js_distance(p, q) == pytest.approx(js_distance(q, p), abs=0)
    # two-bin closed form: P = (1, 0), Q = (1/2, 1/2) -> JSD = 1.5 - 0.75 log2 3
    expected = np.sqrt(0.5 * np.log2(1 / 0.75) + 0.5 * (0.5 * np.log2(0.5 / 0.75) + 0.5 * np.log2(0.5 / 0.25)))
    assert js_distance([1, 0], [1, 1]) == pytest.approx(expected, abs=1e-15)
    assert js_distance([1, 0], [1, 1]) == pytest.approx(np.sqrt(1.5 - 0.75 * np.log2(3)), abs=1e-15)
    for bad in (([1, 2], [1, 2, 3]), ([0, 0], [1, 1]), ([-1, 2], [1, 1])):
        with pytest.raises(ValueError):
            js_distance(*bad)


def test_js_distance_matches_scipy():
    sp = pytest.importorskip("scipy.spatial.distance")
    rng = np.random.default_rng(0)
    for _ in range(50):
        p, q = rng.integers(0, 20, 30), rng.integers(0, 20, 30)
        p[0] += 1
        q[1] += 1
        assert js_distance(p, q) == pytest.approx(sp.jensenshannon(p, q, base=2), abs=1e-12)


def test_linear_hist_counts_out_of_range():
    c, oor = linear_hist(np.array([2.9, 3.0, 5.5, 8.0, 8.1, 9.0]), 3.0, 8.0, 5)
    assert c.sum() == 6 and oor == {"underflow": 1, "overflow": 2}
    assert c[:5].tolist() == [1, 0, 1, 0, 1] and c[5:].tolist() == [1, 2]
    with pytest.raises(ValueError):
        linear_hist(np.array([np.nan]), 0, 1, 2)


def test_circular_hist_wraps_and_origin():
    eps = 1e-9
    near_pi = np.array([np.pi - eps, -np.pi + eps, np.pi])
    c = circular_hist(near_pi, 36, origin=-np.pi)
    assert c.sum() == 3 and c[0] == 2 and c[-1] == 1         # π ≡ −π lands in the first arc
    shifted = circular_hist(near_pi, 36, origin=-np.pi + np.pi / 36)
    assert shifted.sum() == 3 and shifted[-1] == 3           # all three share the arc across ±π
    # angles that differ by 2π are binned identically
    t = np.random.default_rng(1).uniform(-np.pi, np.pi, 200)
    assert (circular_hist(t, 18) == circular_hist(t + 2 * np.pi, 18)).all()
    # distributions straddling ±π are identical in JS distance after a 2π relabelling
    a, b = np.array([3.1, -3.1, 3.0]), np.array([3.1 - 2 * np.pi, -3.1 + 2 * np.pi, 3.0 - 2 * np.pi])
    assert js_distance(circular_hist(a, 36), circular_hist(b, 36)) == 0.0


def test_window_distances_mapping():
    X = np.zeros((2, 5, 3))
    X[:, :, 0] = np.arange(5) * 3.8
    X[1, 4, 0] += 1.0
    d2, d3 = window_distances(X, 2), window_distances(X, 3)
    assert d2.shape == (2, 5) and d3.shape == (2, 4)
    assert np.isnan(d2[:, [0, 4]]).all() and d2[0, 1] == pytest.approx(7.6)   # κ index 1 ↔ residues 0, 2
    assert np.isnan(d3[:, [0, 3]]).all() and d3[0, 1] == pytest.approx(11.4)  # τ hinge 1 ↔ residues 0, 3
    assert d3[1, 2] == pytest.approx(12.4) and d2[1, 3] == pytest.approx(8.6)  # residue 4 moved in frame 1


def test_equal_duration_blocks_and_equalize():
    t = np.arange(0, 335400 + 1, 800.0)[:420]               # 0 … 335 200 ps, every 800 ps
    blocks = equal_duration_blocks(t, 0.0, 335400.0, 3)
    assert [len(b) for b in blocks] == [140, 140, 140]
    assert np.concatenate(blocks).tolist() == list(range(420))  # contiguous, disjoint, complete
    eq, n = equalize([np.arange(5), np.arange(5, 9), np.arange(9, 15)])
    assert n == 4 and [b.tolist() for b in eq] == [[0, 1, 2, 3], [5, 6, 7, 8], [9, 10, 11, 12]]
    with pytest.raises(ValueError):
        equalize([np.arange(3), np.array([], dtype=int)])


def test_compare_channel_identical_shifted_and_excluded():
    rng = np.random.default_rng(2)
    n = 300
    a = rng.normal(5.0, 0.2, size=(n, 3))
    b = a.copy()
    b[:, 1] += 1.0                                          # window 1 shifted
    va = np.ones_like(a, bool)
    vb = va.copy()
    vb[10, 2] = False                                       # window 2 has an invalid sample
    blocks = [np.arange(0, 100), np.arange(100, 200), np.arange(200, 300)]
    out = compare_channel(a, b, va, vb, blocks, blocks, lambda x: linear_hist(x, 3.0, 8.0, 36)[0])
    assert out[0]["full_between"] == 0.0 and out[0]["fraction_block_between_above_max_within"] == 0.0
    assert out[1]["full_between"] > 0.9 and out[1]["fraction_block_between_above_max_within"] == 1.0
    assert out[2]["excluded"] and out[2]["invalid_samples"] == {"a": 0, "b": 1}
    assert len(out[1]["block_between"]) == 9 and len(out[1]["block_within"]) == 6


def test_spearman():
    assert spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    assert spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)
    assert spearman([1, 1, 2], [3, 3, 4]) == pytest.approx(1.0)


def _h(center, rng, n=400):
    return linear_hist(rng.normal(center, 0.3, n), 3.0, 8.0, 36)[0]


def test_replica_comparison_categories():
    from topofold_inputs.compare import replica_comparison
    rng = np.random.default_rng(4)
    same = replica_comparison([_h(5.0, rng) for _ in range(5)], [_h(5.0, rng) for _ in range(5)])
    assert same["category"] == "inconclusive"
    shifted = replica_comparison([_h(4.5, rng) for _ in range(5)], [_h(6.0, rng) for _ in range(5)])
    assert shifted["category"] == "consistent" and shifted["fraction_between_above_max_within"] == 1.0
    assert len(shifted["pairs"]["between"]) == 25 and len(shifted["pairs"]["within_a"]) == 10
    assert len(shifted["loro"]) == 10 and all(r["fraction"] == 1.0 for r in shifted["loro"])
    # one outlier replica in condition A overlaps condition B: separation depends on leaving it out
    a = [_h(4.5, rng) for _ in range(4)] + [_h(6.0, rng)]
    dep = replica_comparison(a, [_h(6.0, rng) for _ in range(5)])
    assert dep["category"] == "replica_dependent" and dep["separating_when_left_out"] == ["a4"]
    lo, hi = dep["pooled_loro_range"]
    assert 0.0 <= lo <= hi <= 1.0
    with pytest.raises(ValueError):
        replica_comparison([_h(5.0, rng)], [_h(5.0, rng), _h(5.0, rng)])
