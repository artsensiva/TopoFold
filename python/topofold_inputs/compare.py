"""Descriptive histogram comparison of per-window geometric features between two runs (development tool).

Histogram-based comparison in the spirit of established per-feature ensemble comparison (e.g. PENSA's
Jensen–Shannon distances); it is not a reproduction or validation of PENSA. Nothing here is a statistical test.

* **Jensen–Shannon distance** (not divergence): ``sqrt(JSD)`` with base-2 logarithms, so it lies in [0, 1]. Empty
  bins contribute 0 (0·log 0 = 0); no pseudocounts.
* **Linear features** use fixed edges ``[lo, hi]`` plus explicit underflow (< lo) and overflow (> hi) bins, so every
  observation is counted; the out-of-range counts are reported.
* **Circular features** (radians) use ``bins`` equal arcs starting at ``origin``; every angle falls in one arc.
* **Windows.** κ at residue i and the Cα(i−1)–Cα(i+1) distance share the three-residue window (i−1, i, i+1);
  τ at hinge i and the Cα(i−1)–Cα(i+2) distance share the four-residue window (i−1, …, i+2).
"""

from __future__ import annotations

from typing import Any, Sequence

import numpy as np

TWO_PI = 2.0 * np.pi


def js_distance(p_counts: np.ndarray, q_counts: np.ndarray) -> float:
    """Jensen–Shannon distance (base 2) between two histograms given as non-negative counts."""
    p = np.asarray(p_counts, dtype=np.float64)
    q = np.asarray(q_counts, dtype=np.float64)
    if p.shape != q.shape or p.ndim != 1:
        raise ValueError("histograms must be 1-D with identical bins")
    if (p < 0).any() or (q < 0).any() or p.sum() <= 0 or q.sum() <= 0:
        raise ValueError("histograms need non-negative counts and at least one observation each")
    p, q = p / p.sum(), q / q.sum()
    m = 0.5 * (p + q)

    def kl(a):
        nz = a > 0
        return float(np.sum(a[nz] * np.log2(a[nz] / m[nz])))

    jsd = 0.5 * kl(p) + 0.5 * kl(q)
    return float(np.sqrt(min(max(jsd, 0.0), 1.0)))


def linear_hist(x: np.ndarray, lo: float, hi: float, bins: int) -> tuple[np.ndarray, dict[str, int]]:
    """Counts in ``bins`` equal bins over [lo, hi] followed by underflow and overflow bins."""
    x = np.asarray(x, dtype=np.float64)
    if not np.isfinite(x).all():
        raise ValueError("non-finite values")
    inside = (x >= lo) & (x <= hi)
    counts, _ = np.histogram(x[inside], bins=bins, range=(lo, hi))
    under, over = int((x < lo).sum()), int((x > hi).sum())
    return np.concatenate([counts, [under, over]]), {"underflow": under, "overflow": over}


def circular_hist(theta: np.ndarray, bins: int, origin: float = -np.pi) -> np.ndarray:
    """Counts in ``bins`` equal arcs of the circle starting at ``origin`` (radians); all angles are counted."""
    t = np.asarray(theta, dtype=np.float64)
    if not np.isfinite(t).all():
        raise ValueError("non-finite values")
    w = TWO_PI / bins
    idx = np.floor(np.mod(t - origin, TWO_PI) / w).astype(np.int64)
    idx = np.clip(idx, 0, bins - 1)  # mod can return exactly 2π after rounding
    return np.bincount(idx, minlength=bins)


def window_distances(coords: np.ndarray, offset: int) -> np.ndarray:
    """Cα(i−1)–Cα(i−1+offset) distances aligned to the feature index i used by κ (offset 2) or τ hinges (offset 3).

    Returns (frames, n_index) with NaN where the window leaves the chain; n_index is N for offset 2, N−1 for 3.
    """
    F, N, _ = coords.shape
    n_index = N if offset == 2 else N - 1
    out = np.full((F, n_index), np.nan)
    for i in range(n_index):
        a, b = i - 1, i - 1 + offset
        if a >= 0 and b < N:
            out[:, i] = np.linalg.norm(coords[:, b] - coords[:, a], axis=-1)
    return out


def equal_duration_blocks(times: np.ndarray, t_start: float, t_end: float, n_blocks: int) -> list[np.ndarray]:
    """Sample indices of ``n_blocks`` contiguous equal-duration blocks of [t_start, t_end] (last edge inclusive)."""
    t = np.asarray(times, dtype=np.float64)
    edges = t_start + (t_end - t_start) * np.arange(n_blocks + 1) / n_blocks
    blocks = []
    for b in range(n_blocks):
        hi_ok = t <= edges[b + 1] if b == n_blocks - 1 else t < edges[b + 1]
        blocks.append(np.flatnonzero((t >= edges[b]) & hi_ok))
    return blocks


def equalize(blocks: Sequence[np.ndarray]) -> tuple[list[np.ndarray], int]:
    """Truncate every block to the smallest block size, keeping each block's earliest samples."""
    n = min(len(b) for b in blocks)
    if n == 0:
        raise ValueError("an empty block")
    return [np.asarray(b)[:n] for b in blocks], n


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman rank correlation with average ranks for ties (descriptive)."""
    def ranks(x):
        x = np.asarray(x, dtype=np.float64)
        order = np.argsort(x, kind="mergesort")
        r = np.empty(len(x))
        r[order] = np.arange(len(x), dtype=np.float64)
        for v in np.unique(x):
            m = x == v
            if m.sum() > 1:
                r[m] = r[m].mean()
        return r
    ra, rb = ranks(a), ranks(b)
    ra, rb = ra - ra.mean(), rb - rb.mean()
    den = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / den) if den > 0 else float("nan")


def compare_channel(x_a: np.ndarray, x_b: np.ndarray, valid_a: np.ndarray, valid_b: np.ndarray,
                    blocks_a: Sequence[np.ndarray], blocks_b: Sequence[np.ndarray], hist) -> list[dict[str, Any]]:
    """Per-window full-run and block-level JS distances for one feature channel.

    ``x_*`` are (frames, windows); ``valid_*`` masks select matched valid observations; ``hist`` maps a 1-D array to
    counts with identical bins for every call. Windows with any invalid observation in either run are excluded
    (reported), so that every distance uses the same samples.
    """
    out = []
    for w in range(x_a.shape[1]):
        va, vb = valid_a[:, w], valid_b[:, w]
        if not (va.all() and vb.all()):
            out.append({"window": w, "excluded": True,
                        "invalid_samples": {"a": int((~va).sum()), "b": int((~vb).sum())}})
            continue
        a, b = x_a[:, w], x_b[:, w]
        ha = [hist(a[idx]) for idx in blocks_a]
        hb = [hist(b[idx]) for idx in blocks_b]
        between = [js_distance(p, q) for p in ha for q in hb]
        within = [js_distance(ha[i], ha[j]) for i in range(len(ha)) for j in range(i + 1, len(ha))] + \
                 [js_distance(hb[i], hb[j]) for i in range(len(hb)) for j in range(i + 1, len(hb))]
        max_within = max(within)
        out.append({
            "window": w, "excluded": False,
            "full_between": js_distance(hist(a), hist(b)),
            "block_between": between, "block_within": within,
            "matched_block_between": [js_distance(ha[k], hb[k]) for k in range(len(ha))],
            "median_block_between": float(np.median(between)), "max_block_within": max_within,
            "fraction_block_between_above_max_within": float(np.mean(np.array(between) > max_within)),
        })
    return out
