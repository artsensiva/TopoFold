"""Tests for the SYNTHETIC REPCHAR-1 generators (marginals, autocorrelation, mixture weights and dwell, wrapping)."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/audit/synth_rep_char"))
from synthetic_generators import (  # noqa: E402
    ar1_latent, circular_series, ess_gaussian_mean, linear_series, markov_states, mixture_series, wrap,
)


def lag1(x):
    x = x - x.mean()
    return float((x[:-1] * x[1:]).mean() / (x * x).mean())


@pytest.mark.parametrize("phi", [0.0, 0.8, 0.95])
def test_ar1_marginal_and_autocorrelation(phi):
    rng = np.random.default_rng(1)
    z = np.concatenate([ar1_latent(rng, 2000, phi) for _ in range(100)])          # 200k samples
    assert abs(z.mean()) < 0.05 and abs(z.std() - 1.0) < 0.03                     # stationary N(0, 1) marginal
    first = np.array([ar1_latent(rng, 5, phi)[0] for _ in range(20000)])
    assert abs(first.std() - 1.0) < 0.03                                           # stationary initialization
    r = np.mean([lag1(ar1_latent(rng, 2000, phi)) for _ in range(50)])
    assert abs(r - phi) < 0.02


def test_linear_and_circular_marginals():
    rng = np.random.default_rng(2)
    x = np.concatenate([linear_series(rng, 1000, 5.5, 0.2, 0.8) for _ in range(100)])
    assert abs(x.mean() - 5.5) < 0.01 and abs(x.std() - 0.2) < 0.005
    th = np.concatenate([circular_series(rng, 1000, np.pi, 0.3, 0.8) for _ in range(100)])
    assert th.min() > -np.pi and th.max() <= np.pi                                 # wrapped to (-pi, pi]
    assert (th > 3.0).any() and (th < -3.0).any()                                  # straddles the branch cut
    R = abs(np.mean(np.exp(1j * th)))
    assert abs(np.sqrt(-2 * np.log(R)) - 0.3) < 0.01                                # circular SD = s
    assert abs(np.angle(np.mean(np.exp(1j * th))) - np.pi) < 0.01 or abs(np.angle(np.mean(np.exp(1j * th))) + np.pi) < 0.01
    assert wrap(np.pi) == pytest.approx(np.pi) and wrap(-np.pi) == pytest.approx(np.pi)


@pytest.mark.parametrize("w,D", [(0.5, 1), (0.3, 20), (0.1, 100), (0.4, 20)])
def test_markov_weight_and_dwell(w, D):
    rng = np.random.default_rng(3)
    s = np.concatenate([markov_states(rng, 5000, w, D) for _ in range(40)])
    assert abs(s.mean() - w) < (0.02 if D <= 20 else 0.04)                        # stationary weight (slow chains vary more)
    runs1 = np.diff(np.flatnonzero(np.diff(np.concatenate([[0], s, [0]]) != 0)))[::2]
    assert abs(runs1.mean() - D / (1 - w)) / (D / (1 - w)) < 0.15                  # mean dwell in state 1
    first = np.array([markov_states(rng, 2, w, D)[0] for _ in range(20000)])
    assert abs(first.mean() - w) < 0.015                                           # stationary initialization


def test_mixture_and_validation():
    rng = np.random.default_rng(4)
    x, s = mixture_series(rng, 20000, 5.1, 5.9, 0.2, 0.3, 1)
    assert abs(x[s == 1].mean() - 5.9) < 0.01 and abs(x[s == 0].std() - 0.2) < 0.01
    with pytest.raises(ValueError):
        markov_states(rng, 10, 0.5, 0.5)
    with pytest.raises(ValueError):
        ar1_latent(rng, 10, 1.0)
    assert ess_gaussian_mean(387, 0.8) == pytest.approx(387 * 0.2 / 1.8)
