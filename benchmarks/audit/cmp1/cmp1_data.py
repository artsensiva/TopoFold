"""CMP-1 canonical input (protocol §5): one read-only float64 array in Å, float32-representable.

Real workloads: S2 rounded arrays (already float32-representable; verified). SYNTHETIC_tiled workloads: BQ-2 rows,
recast once through float32 (documented change from BQ-2/BQ-3). Loader: bq2_run.load (unchanged, hash-checked).
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'bq2'))
sys.path.insert(0, str(HERE.parent / 'eb1'))
import bq2_run as B2R  # noqa: E402


def load(w):
    """w: {'kind': 'real'|'SYNTHETIC_tiled', 'd', 'nd', 'nwin', 'queries'}. Returns (X read-only, di, windows, info)."""
    X, di, wins = B2R.load(w)
    rep = X.astype(np.float32).astype(np.float64)
    info = {'float32_representable_before_recast': bool(np.array_equal(rep, X))}
    if w['kind'] != 'real':
        X = rep
    elif not info['float32_representable_before_recast']:
        raise ValueError('real canonical input is not float32-representable')
    X = np.ascontiguousarray(X)
    X.setflags(write=False)
    return X, di, wins, info
