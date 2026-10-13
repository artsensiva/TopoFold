"""BQ-2 helpers (docs/B2_SCALING_CONFIRMATION_PROTOCOL.md). Harness code, not production.

The V2b solver, fallback, distance objective and ranking policy are imported unchanged from b2opt_lib (hash-checked
by the runner). This module adds only:
- `tiled_rows`: SYNTHETIC/performance-only L = 1000 chains built by concatenating real D1+ frames end to end
  (repeated data; junction windows are artificial). Used only for scaling, never for biological claims.
- `prepare_reps`: repeated preparation with each prepared array released before the next repetition.
- `predicted_peak_bytes` / `memory_ok`: the frozen memory-headroom rule.
"""
import gc
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'b2opt'))
import b2opt_lib as B  # noqa: E402  (V2b unchanged)

L_STRESS = 1000
SEG = 214                          # D1+ chain length
JUNCTION_A = 3.8                   # Å, gap placed between consecutive segments
QUERY_OFFSET = 100_000             # tiled query rows use a frame sequence disjoint in index from database rows
BASE_RSS = 60 * 1024 ** 2          # interpreter + NumPy baseline allowance (bytes)
TRANSIENT = 64 * 1024 ** 2         # per-window transient allowance (bytes)
HEADROOM = 600 * 1024 ** 2         # MemAvailable that must remain after the predicted peak


def tiled_rows(src, rows, L=L_STRESS):
    """SYNTHETIC: row r = segments of real frames (5r + j) mod F, j = 0…⌈L/214⌉−1, each translated so it starts
    3.8 Å along +x from the previous segment's last atom; truncated to L residues. src: (F, 214, 3) float64."""
    F = src.shape[0]
    nseg = -(-L // SEG)
    out = np.empty((len(rows), nseg * SEG, 3))
    for i, r in enumerate(rows):
        end = np.zeros(3)
        for j in range(nseg):
            seg = src[(5 * r + j) % F]
            start = end + np.array([JUNCTION_A, 0.0, 0.0]) if j else seg[0]
            seg = seg - seg[0] + start
            out[i, j * SEG:(j + 1) * SEG] = seg
            end = seg[-1]
    return np.ascontiguousarray(out[:, :L])


def predicted_peak_bytes(n_rows, L, n_d, n_win, W=8):
    return BASE_RSS + TRANSIENT + n_rows * L * 3 * 8 + n_d * n_win * (W * 3 + 1) * 8


def memory_ok(pred_bytes, mem_available_kb):
    return mem_available_kb * 1024 - pred_bytes >= HEADROOM


def prep_reps_for(pred_s):
    """Frozen repetition rule for preparation timing (predicted from the BQ-1 V2b coefficients)."""
    return 9 if pred_s < 0.05 else (5 if pred_s < 1.0 else 3)


def prepare_reps(X, di, wins, reps):
    """Time `reps` preparations, releasing each prepared array before the next; returns the last (kept) and timings."""
    walls, cpus, prep = [], [], None
    for r in range(reps):
        if prep is not None:
            del prep
            prep = None
            gc.collect()
        t0, c0 = time.perf_counter(), time.process_time()
        prep = B.prepare(X, di, wins)
        walls.append(time.perf_counter() - t0)
        cpus.append(time.process_time() - c0)
    return prep, walls, cpus
