"""EB-1 Amendment 1 library (docs/EB1_AMENDMENT_1.md). Harness code, not production; no optimization implemented.

- Method-specific operational paths: A = S1 → S2 → S3 (κ, τ) → S4_A; B2 = S1 → S2 → S4_B2 (no S3).
- Work and memory accounting with effective batch b = min(n_Q, BATCH_LIMIT), retained top-k outputs over all
  windows, resident inputs owned by the process (copies, never slices of a larger backing array).
- `search_instrumented` performs exactly the EB-1 harness computation (eb1_lib.search) while timing its parts:
  per-call slice copies (harness overhead), distance kernel, top-k selection, and the remaining Python loop.
"""
import math
import time

import numpy as np

import eb1_lib as E

BATCH_LIMIT = 128
K = 10
STAGES = {'A': ('S1', 'S2', 'S3', 'S4_A'), 'B2': ('S1', 'S2', 'S4_B2')}


# --- work counts ---------------------------------------------------------------------------------------------------

def effective_batch(n_q, batch_limit=BATCH_LIMIT):
    return max(1, min(n_q, batch_limit))


def work(n_q, n_d, n_win, batch_limit=BATCH_LIMIT):
    """Comparisons and kernel calls for exhaustive search over n_win windows."""
    b = effective_batch(n_q, batch_limit)
    calls = n_win * math.ceil(n_q / b)
    return {'comparisons': n_q * n_d * n_win, 'calls': calls, 'call_db_rows': calls * n_d, 'b': b}


# --- memory accounting (bytes) -------------------------------------------------------------------------------------

def resident_bytes(method, n_frames, L):
    """Operational resident inputs: A holds κ (L−2) and τ (L−3) float64 per frame; B2 holds Cα float64 (L, 3)."""
    return n_frames * ((2 * L - 5) * 8 if method == 'A' else L * 3 * 8)


def retained_bytes(n_q, n_win, k=K):
    """Top-k indices (int64) retained for every query and window."""
    return n_q * n_win * k * 8


def temp_bytes(method, b, n_d, W, coef):
    """Peak transient allocation of one kernel call: coef['alpha'][method] bytes per (query, database) pair plus
    coef['beta'][method] bytes per database row (per-call slice copies); calibrated before measurement."""
    return coef['alpha'][method] * b * n_d + coef['beta'][method] * n_d


def predicted_bytes(method, n_frames, L, n_q, n_d, n_win, W, coef, batch_limit=BATCH_LIMIT):
    b = effective_batch(n_q, batch_limit)
    return resident_bytes(method, n_frames, L) + retained_bytes(n_q, n_win) + temp_bytes(method, b, n_d, W, coef)


# --- instrumented search (identical computation to eb1_lib.search) -------------------------------------------------

def search_instrumented(method, W, windows, feat, qi, di, batch=BATCH_LIMIT):
    """feat: {'X'} for B2 or {'kap', 'tau'} for A. Returns (results {s: (n_q, k)}, timing parts in s)."""
    t = {'slice_copy': 0.0, 'kernel': 0.0, 'topk': 0.0}
    t_all = time.perf_counter()
    out = {}
    for s in windows:
        res = []
        for b0 in range(0, len(qi), batch):
            q = qi[b0:b0 + batch]
            t0 = time.perf_counter()
            if method == 'A':
                kq, tq = feat['kap'][q, s:s + W - 2], feat['tau'][q, s:s + W - 3]
                kd, td = feat['kap'][di, s:s + W - 2], feat['tau'][di, s:s + W - 3]
                t1 = time.perf_counter()
                D = E.d_aligned(kq, tq, kd, td)
            else:
                xq, xd = feat['X'][q, s:s + W], feat['X'][di, s:s + W]
                t1 = time.perf_counter()
                D = E.horn_rmsd(xq, xd)
            t2 = time.perf_counter()
            res.append(E.topk(D))
            t3 = time.perf_counter()
            t['slice_copy'] += t1 - t0
            t['kernel'] += t2 - t1
            t['topk'] += t3 - t2
        out[s] = np.vstack(res)
    total = time.perf_counter() - t_all
    t['loop_other'] = total - sum(t.values())
    t['total'] = total
    return out, t


def owned_subset(arr, idx):
    """Copy rows so that the result owns its data (no retained backing allocation)."""
    out = np.ascontiguousarray(arr[idx])
    assert out.base is None or out.base is not arr
    return out
