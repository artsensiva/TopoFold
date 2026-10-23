"""Prototype: dependence-aware two-sample comparison of marginal distributions (evaluated in SYNTHETIC_DABOOT-1 only).

Method (established components, no new estimator):

* **Statistic:** two-sample Kuiper statistic V = max_g D(g) − min_g D(g) with D = F_A − F_B, where F_c is the empirical
  CDF of all frames of condition c (its replicas pooled), evaluated on a fixed grid. Kuiper (1960) is invariant to the
  choice of origin for circular data, so linear and circular features share one statistic.
* **Null distribution:** centered stationary bootstrap (Politis & Romano 1994), applied to each condition separately,
  so the two conditions may have different temporal dependence. For each bootstrap draw and each replica, a
  stationary-bootstrap resample of that replica's own frames is taken (geometric block lengths with mean L_c,
  wrapping circularly **within** the replica). Replicas are never joined to each other.
  D*(g) = (F*_A − F_A)(g) − (F*_B − F_B)(g), V* = max D* − min D*, p = (1 + #{V* ≥ V}) / (B + 1).
* **Block length:** Politis & White (2004) automatic selection for the stationary bootstrap, with the correction of
  Patton, Politis & White (2009) (D_SB = 2 ĝ(0)²). Applied to each replica's scalar series (circular: the larger of
  the values for cos θ and sin θ). The condition's L_c is the median over its replicas.

Assumptions: within each replica the series is stationary and weakly dependent, with correlation length short
compared with the replica; replicas are independent given their start and identically distributed within a
condition. Validity is asymptotic; finite-sample calibration must be measured (SYNTHETIC_DABOOT-1).
"""

from __future__ import annotations

import math
from typing import Any, Sequence

import numpy as np

TWO_PI = 2.0 * math.pi


# --- block length ---------------------------------------------------------------------------------------------------

def _flat_top(t: np.ndarray) -> np.ndarray:
    """Trapezoidal flat-top lag window: 1 on |t| <= 1/2, 2(1 − |t|) on 1/2 < |t| <= 1, 0 beyond."""
    a = np.abs(t)
    return np.where(a <= 0.5, 1.0, np.where(a <= 1.0, 2.0 * (1.0 - a), 0.0))


def autocovariances(x: np.ndarray, max_lag: int) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    n = len(x)
    xc = x - x.mean()
    return np.array([np.dot(xc[: n - k], xc[k:]) / n for k in range(max_lag + 1)])


def sb_block_length_from_acov(R: np.ndarray, n: int, M: float) -> float:
    """b_opt,SB = (2 Ĝ² / D̂_SB)^(1/3) n^(1/3), D̂_SB = 2 ĝ(0)², with flat-top weights λ(k/M) on lags |k| <= M."""
    k = np.arange(len(R))
    w = _flat_top(k / M) if M > 0 else (k == 0).astype(float)
    g0 = R[0] * w[0] + 2.0 * np.sum(w[1:] * R[1:])
    G = 2.0 * np.sum(w[1:] * k[1:] * R[1:])
    if g0 <= 0:
        return 1.0
    return float((2.0 * G * G / (2.0 * g0 * g0)) ** (1.0 / 3.0) * n ** (1.0 / 3.0))


def politis_white_sb(x: np.ndarray) -> float:
    """Automatic expected block length for the stationary bootstrap of one scalar series.

    Settings follow Politis & White (2004) with the Patton, Politis & White (2009) correction and the defaults of the
    R package ``blocklength``: K_N = max(5, ln n), c = 1.96, threshold c·sqrt(log10(n)/n) on the sample
    autocorrelation, m̂ = smallest m with |ρ̂(m + k)| below threshold for k = 1..K_N, M = 2 m̂,
    result clipped to [1, b_max], b_max = ceil(min(3 sqrt(n), n / 3)).
    """
    x = np.asarray(x, dtype=np.float64)
    n = len(x)
    K = max(5, int(math.ceil(math.log(n))))
    m_max = int(math.ceil(math.sqrt(n))) + K
    R = autocovariances(x, min(n - 1, 2 * m_max + K))
    if R[0] <= 0:
        return 1.0
    rho = R / R[0]
    thr = 1.96 * math.sqrt(math.log10(n) / n)
    m_hat = m_max
    for m in range(0, m_max + 1):
        if np.all(np.abs(rho[m + 1: m + 1 + K]) < thr):
            m_hat = m
            break
    M = 2 * m_hat
    b = sb_block_length_from_acov(R[: M + 1], n, M) if M > 0 else 1.0
    b_max = math.ceil(min(3 * math.sqrt(n), n / 3))
    return float(min(max(b, 1.0), b_max))


def block_length(x: np.ndarray, circular: bool) -> float:
    if circular:
        return max(politis_white_sb(np.cos(x)), politis_white_sb(np.sin(x)))
    return politis_white_sb(x)


# --- stationary bootstrap -------------------------------------------------------------------------------------------

def stationary_bootstrap_indices(rng: np.random.Generator, n: int, mean_block: float, B: int) -> np.ndarray:
    """(B, n) indices into one series: each position starts a new block with probability 1/L (uniform start),
    otherwise continues the previous index + 1, wrapping circularly within the series."""
    p = 1.0 / max(mean_block, 1.0)
    idx = np.empty((B, n), dtype=np.int64)
    idx[:, 0] = rng.integers(0, n, B)
    new = rng.random((B, n)) < p
    starts = rng.integers(0, n, (B, n))
    for t in range(1, n):
        idx[:, t] = np.where(new[:, t], starts[:, t], (idx[:, t - 1] + 1) % n)
    return idx


# --- ECDF on a grid and Kuiper ----------------------------------------------------------------------------------------

def grid_positions(x: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """For each value, the first grid index j with x <= grid[j] (len(grid) if none); F(grid[j]) counts positions <= j."""
    return np.searchsorted(grid, np.asarray(x, dtype=np.float64), side="left")


def ecdf_from_positions(pos: np.ndarray, G: int) -> np.ndarray:
    """ECDF on the grid from positions; pos may be 1-D (one sample) or 2-D (rows = bootstrap draws)."""
    pos = np.atleast_2d(pos)
    B, n = pos.shape
    flat = (pos + (G + 1) * np.arange(B)[:, None]).ravel()
    counts = np.bincount(flat, minlength=B * (G + 1)).reshape(B, G + 1)[:, :G]
    return np.cumsum(counts, axis=1) / n


def kuiper(d: np.ndarray) -> np.ndarray:
    return np.max(d, axis=-1) - np.min(d, axis=-1)


def two_sample_test(reps_a: Sequence[np.ndarray], reps_b: Sequence[np.ndarray], grid: np.ndarray, B: int,
                    rng: np.random.Generator, circular: bool) -> dict[str, Any]:
    """Centered stationary-bootstrap Kuiper test of equal marginal distributions (replicas resampled separately)."""
    G = len(grid)
    out = {}
    F, Fstar = {}, {}
    for lab, reps in (("a", reps_a), ("b", reps_b)):
        lens = [block_length(np.asarray(r), circular) for r in reps]
        L = float(np.median(lens))
        pos = [grid_positions(r, grid) for r in reps]
        n_tot = sum(len(r) for r in reps)
        F[lab] = sum(ecdf_from_positions(p_, G)[0] * len(p_) for p_ in pos) / n_tot
        acc = np.zeros((B, G))
        for p_ in pos:
            idx = stationary_bootstrap_indices(rng, len(p_), L, B)
            acc += ecdf_from_positions(p_[idx], G) * len(p_)
        Fstar[lab] = acc / n_tot
        out[f"block_length_{lab}"] = L
        out[f"block_lengths_{lab}"] = lens
    v = float(kuiper(F["a"] - F["b"]))
    vstar = kuiper((Fstar["a"] - F["a"]) - (Fstar["b"] - F["b"]))
    out.update(V=v, p=float((1 + np.sum(vstar >= v - 1e-12)) / (B + 1)), B=B)
    return out
