"""Tests for the dependence-aware bootstrap prototype (SYNTHETIC inputs; analytical and published references)."""

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from topofold_inputs import dependent_bootstrap as db  # noqa: E402


@pytest.mark.parametrize("rho,n,expected", [
    (0.7, 200, 11.47), (0.7, 800, 18.20), (0.1, 200, 2.01), (0.1, 800, 3.20), (-0.4, 200, 5.66), (-0.4, 800, 8.99),
])
def test_theoretical_block_length_matches_patton_politis_white_table1(rho, n, expected):
    # exact AR(1) autocovariances (unit innovations), all lags with weight 1 (window wide enough to be exact)
    K = 2000
    R = np.array([rho ** k / (1 - rho * rho) for k in range(K + 1)])
    b = db.sb_block_length_from_acov(R, n, M=1e9)
    assert b == pytest.approx(expected, abs=0.01)
    assert b == pytest.approx((abs(2 * rho / (1 - rho * rho))) ** (2 / 3) * n ** (1 / 3), rel=1e-9)


def test_estimated_block_length_tracks_theory_on_long_ar1():
    rng = np.random.default_rng(0)
    for rho in (0.0, 0.8):
        n = 20000
        e = rng.standard_normal(n)
        x = np.empty(n)
        x[0] = e[0] / math.sqrt(1 - rho * rho) if rho else e[0]
        for t in range(1, n):
            x[t] = rho * x[t - 1] + e[t]
        b = db.politis_white_sb(x)
        theory = (2 * rho / (1 - rho * rho)) ** (2 / 3) * n ** (1 / 3) if rho else 1.0
        if rho:
            assert 0.6 * theory < b < 1.4 * theory
        else:
            assert b < 3.0
        assert 1.0 <= b <= math.ceil(min(3 * math.sqrt(n), n / 3))


def test_stationary_bootstrap_indices():
    rng = np.random.default_rng(1)
    n, L, B = 387, 20.0, 400
    idx = db.stationary_bootstrap_indices(rng, n, L, B)
    assert idx.shape == (B, n) and idx.min() >= 0 and idx.max() < n
    cont = (idx[:, 1:] == (idx[:, :-1] + 1) % n)
    # a new block starts with probability 1/L; a "new" start can coincide with continuation with probability 1/n
    p_new = 1 - cont.mean()
    assert abs(p_new - (1 / L) * (1 - 1 / n)) < 0.004
    iid = db.stationary_bootstrap_indices(rng, n, 1.0, 50)
    assert abs((iid[:, 1:] == (iid[:, :-1] + 1) % n).mean() - 1 / n) < 0.003


def test_ecdf_on_grid_matches_brute_force():
    rng = np.random.default_rng(2)
    x = rng.normal(5.5, 0.3, 200)
    grid = np.linspace(3.0, 8.0, 64)
    F = db.ecdf_from_positions(db.grid_positions(x, grid), len(grid))[0]
    brute = np.array([(x <= g).mean() for g in grid])
    np.testing.assert_allclose(F, brute, atol=0)
    idx = rng.integers(0, 200, (5, 200))
    F2 = db.ecdf_from_positions(db.grid_positions(x, grid)[idx], len(grid))
    for b in range(5):
        np.testing.assert_allclose(F2[b], [(x[idx[b]] <= g).mean() for g in grid], atol=0)


def test_kuiper_analytical_and_rotation_invariance():
    grid = np.linspace(-math.pi, math.pi, 73)[1:]                  # 72 arcs of 5 degrees, origin -pi
    G = len(grid)
    a, b = np.full(10, 0.3), np.full(10, 2.0)
    Fa = db.ecdf_from_positions(db.grid_positions(a, grid), G)[0]
    Fb = db.ecdf_from_positions(db.grid_positions(b, grid), G)[0]
    assert db.kuiper(Fa - Fb) == pytest.approx(1.0)                 # disjoint point masses
    assert db.kuiper(Fa - Fa) == 0.0
    rng = np.random.default_rng(3)
    x, y = rng.vonmises(0.5, 2.0, 300), rng.vonmises(1.5, 2.0, 300)
    step = 2 * math.pi / 72

    def V(u, v):
        return db.kuiper(db.ecdf_from_positions(db.grid_positions(u, grid), G)[0]
                         - db.ecdf_from_positions(db.grid_positions(v, grid), G)[0])

    wrap = lambda t: math.pi - np.mod(math.pi - t, 2 * math.pi)  # noqa: E731
    for k in (7, 30, 71):                                            # rotations by whole grid steps
        assert V(wrap(x + k * step + 1e-9), wrap(y + k * step + 1e-9)) == pytest.approx(V(x + 1e-9, y + 1e-9), abs=1e-12)


def test_identical_samples_give_p_one_and_replicas_stay_separate():
    rng = np.random.default_rng(4)
    reps = [rng.normal(0, 1, 120) for _ in range(3)]
    grid = np.linspace(-4, 4, 50)
    r = db.two_sample_test(reps, [r.copy() for r in reps], grid, 99, np.random.default_rng(5), circular=False)
    assert r["V"] == 0.0 and r["p"] == 1.0
    # resampling never mixes replicas: bootstrap ECDF of a replica uses only its own values
    a = np.zeros(100)
    b = np.ones(100)
    pos_a, pos_b = db.grid_positions(a, grid), db.grid_positions(b, grid)
    idx = db.stationary_bootstrap_indices(np.random.default_rng(6), 100, 10.0, 20)
    assert np.all(db.ecdf_from_positions(pos_a[idx], len(grid)) == db.ecdf_from_positions(pos_a, len(grid))[0])
    assert np.all(db.ecdf_from_positions(pos_b[idx], len(grid)) == db.ecdf_from_positions(pos_b, len(grid))[0])


def test_iid_bootstrap_variance_of_ecdf_point():
    rng = np.random.default_rng(7)
    x = rng.normal(0, 1, 400)
    grid = np.array([0.0])
    pos = db.grid_positions(x, grid)
    idx = db.stationary_bootstrap_indices(rng, 400, 1.0, 4000)
    Fb = db.ecdf_from_positions(pos[idx], 1)[:, 0]
    p = (x <= 0).mean()
    assert Fb.var() == pytest.approx(p * (1 - p) / 400, rel=0.1)


def test_shifted_samples_reject():
    rng = np.random.default_rng(8)
    ra = [rng.normal(0, 1, 200) for _ in range(3)]
    rb = [rng.normal(1.0, 1, 200) for _ in range(3)]
    r = db.two_sample_test(ra, rb, np.linspace(-4, 5, 64), 199, np.random.default_rng(9), circular=False)
    assert r["p"] <= 0.01
