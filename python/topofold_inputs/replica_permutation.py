"""Prospective replica-level permutation method (evaluated on SYNTHETIC data in SYNTHETIC_PERM-1 only).

Unit: one replica histogram. With Jensen–Shannon distances D between all replicas of both conditions, the statistic
for an allocation of replicas into groups A and B is

    T = mean_{i in A, j in B} D_ij − ½ (mean_{i<j in A} D_ij + mean_{i<j in B} D_ij).

The exact upper-tail permutation probability enumerates every labelled allocation of the n_A + n_B replicas into a
group of size n_A and a group of size n_B (252 for 5 vs 5), including the observed one:

    p = #{allocations with T_perm >= T_obs − tol} / #allocations,   tol = 1e-12 · max(1, |T_obs|)   (tie rule).

Because T is symmetric in the two labels, every partition appears twice, so for 5 vs 5 p takes values 2k/252; p ≤ 0.05
then corresponds to at most 12 of the 252 allocations (nominal size 12/252 ≈ 0.048 when the null is exact and T has
no ties).

**Exchangeability assumption:** under the null, the joint distribution of the replica histograms is unchanged by any
permutation of replica labels across conditions. Replica heterogeneity common to both conditions keeps this; a
component shared only within one condition, or an unusual replica confined to one condition, violates it. The
probability is then not a valid null probability, and in no case does it attribute a difference to its cause.
"""

from __future__ import annotations

from itertools import combinations
from typing import Any, Sequence

import numpy as np

from .compare import js_distance

TIE_RTOL = 1e-12


def distance_matrix(hists: Sequence[np.ndarray]) -> np.ndarray:
    n = len(hists)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = js_distance(hists[i], hists[j])
    return D


def allocations(n: int, k: int) -> np.ndarray:
    """All labelled allocations of n replicas with group A of size k, as a (C(n, k), n) boolean array (A = True)."""
    out = np.zeros((0, n), dtype=bool)
    rows = []
    for comb in combinations(range(n), k):
        m = np.zeros(n, dtype=bool)
        m[list(comb)] = True
        rows.append(m)
    return np.array(rows) if rows else out


def statistic(D: np.ndarray, mask_a: np.ndarray) -> np.ndarray:
    """T for one allocation (1-D mask) or many (2-D masks, one row each)."""
    M = np.atleast_2d(mask_a).astype(np.float64)
    B = 1.0 - M
    ka, kb = M.sum(axis=1), B.sum(axis=1)
    if np.any(ka < 2) or np.any(kb < 2):
        raise ValueError("each group needs at least two replicas")
    within_a = 0.5 * np.einsum("ri,ij,rj->r", M, D, M) / (ka * (ka - 1) / 2)
    within_b = 0.5 * np.einsum("ri,ij,rj->r", B, D, B) / (kb * (kb - 1) / 2)
    between = np.einsum("ri,ij,rj->r", M, D, B) / (ka * kb)
    t = between - 0.5 * (within_a + within_b)
    return t if np.ndim(mask_a) == 2 else t[0]


def permutation_test(hists_a: Sequence[np.ndarray], hists_b: Sequence[np.ndarray]) -> dict[str, Any]:
    """Exact upper-tail permutation probability of T over all labelled allocations (observed one included)."""
    hists = list(hists_a) + list(hists_b)
    na, n = len(hists_a), len(hists)
    D = distance_matrix(hists)
    obs = np.zeros(n, dtype=bool)
    obs[:na] = True
    t_obs = float(statistic(D, obs))
    masks = allocations(n, na)
    t_all = statistic(D, masks)
    tol = TIE_RTOL * max(1.0, abs(t_obs))
    n_ge = int(np.sum(t_all >= t_obs - tol))
    return {"t_obs": t_obs, "p": n_ge / len(masks), "n_allocations": int(len(masks)), "n_ge": n_ge}
