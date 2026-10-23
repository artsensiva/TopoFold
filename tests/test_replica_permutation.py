"""Tests for the prospective replica permutation method: hand-calculated example, ties, label symmetry and an
independent brute-force implementation. SYNTHETIC inputs only."""

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from topofold_inputs.compare import linear_hist  # noqa: E402
from topofold_inputs.replica_permutation import allocations, statistic  # noqa: E402
import topofold_inputs.replica_permutation as rp  # noqa: E402


def test_hand_calculated_two_vs_two():
    # replicas at positions 0, 0, 1, 1 with D_ij = |x_i - x_j|
    x = np.array([0.0, 0.0, 1.0, 1.0])
    D = np.abs(x[:, None] - x[None, :])
    masks = allocations(4, 2)
    assert len(masks) == 6
    t = {tuple(np.flatnonzero(m)): float(v) for m, v in zip(masks, statistic(D, masks))}
    # {0,1} vs {2,3}: between 1, within 0 and 0 -> T = 1; same for the label swap
    assert t[(0, 1)] == pytest.approx(1.0) and t[(2, 3)] == pytest.approx(1.0)
    # mixed allocations: between mean (0+1+1+0)/4 = 0.5, within 1 and 1 -> T = -0.5
    for k in [(0, 2), (0, 3), (1, 2), (1, 3)]:
        assert t[k] == pytest.approx(-0.5)
    # observed allocation {0,1}: 2 of 6 allocations reach T_obs -> p = 1/3
    obs = np.array([True, True, False, False])
    t_obs = statistic(D, obs)
    assert np.sum(statistic(D, masks) >= t_obs - 1e-12) / 6 == pytest.approx(1 / 3)


def test_ties_count_toward_p():
    D = np.ones((6, 6)) - np.eye(6)                      # every distance equal: every T equals 0
    masks = allocations(6, 3)
    t = statistic(D, masks)
    assert np.allclose(t, 0.0)
    obs = np.array([1, 1, 1, 0, 0, 0], bool)
    t_obs = float(statistic(D, obs))
    tol = rp.TIE_RTOL * max(1.0, abs(t_obs))
    assert np.sum(t >= t_obs - tol) == len(masks)          # all tied -> p = 1


def _hists(centers, seed):
    rng = np.random.default_rng(seed)
    return [linear_hist(rng.normal(c, 0.2, 387), 3.0, 8.0, 36)[0] for c in centers]


def brute_force(ha, hb):
    """Independent implementation: explicit loops over allocations and pairs."""
    hs = list(ha) + list(hb)
    n, na = len(hs), len(ha)
    D = [[rp.js_distance(hs[i], hs[j]) if i != j else 0.0 for j in range(n)] for i in range(n)]

    def T(A):
        Bs = [i for i in range(n) if i not in A]
        wa = [D[i][j] for i, j in combinations(A, 2)]
        wb = [D[i][j] for i, j in combinations(Bs, 2)]
        bt = [D[i][j] for i in A for j in Bs]
        return sum(bt) / len(bt) - 0.5 * (sum(wa) / len(wa) + sum(wb) / len(wb))

    t_obs = T(tuple(range(na)))
    vals = [T(A) for A in combinations(range(n), na)]
    return t_obs, sum(v >= t_obs - 1e-12 * max(1.0, abs(t_obs)) for v in vals) / len(vals)


@pytest.mark.parametrize("shift,seed", [(0.0, 1), (0.1, 2), (1.0, 3)])
def test_matches_brute_force_and_label_symmetry(shift, seed):
    ha = _hists([5.5] * 5, seed)
    hb = _hists([5.5 + shift * 0.2] * 5, seed + 100)
    r = rp.permutation_test(ha, hb)
    t_bf, p_bf = brute_force(ha, hb)
    assert r["n_allocations"] == 252
    assert r["t_obs"] == pytest.approx(t_bf, abs=1e-12) and r["p"] == pytest.approx(p_bf, abs=0)
    s = rp.permutation_test(hb, ha)                      # swapping condition labels
    assert s["t_obs"] == pytest.approx(r["t_obs"], abs=1e-12) and s["p"] == r["p"]
    assert r["n_ge"] % 2 == 0                             # each partition counted twice (label symmetry)


def test_separated_groups_reach_minimum_p():
    r = rp.permutation_test(_hists([4.5] * 5, 7), _hists([6.5] * 5, 8))
    assert r["p"] == pytest.approx(2 / 252)


def test_group_size_validation():
    with pytest.raises(ValueError):
        statistic(np.zeros((3, 3)), np.array([True, False, False]))
