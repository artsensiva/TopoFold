#!/usr/bin/env python3
"""LE-1 analysis (protocol §7–§8): cell means, paired contrasts, two-way block bootstrap, τ_int diagnostic,
development classification of A and F. Reads the per-(window, query) NR@10 arrays written by run_dev.py.

Third-party MD derivative. Run: ../wpd_local_retrieval/build/venv/bin/python -B analyse_dev.py
"""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RES = HERE / 'results'
SUMMARY = RES / 'le_adk_THIRD_PARTY_MD_DERIVATIVE_summary.json'
OUT = RES / 'le_adk_THIRD_PARTY_MD_DERIVATIVE_analysis.json'
LOCAL = ('M1', 'M2', 'K', 'A', 'F', 'F_lock', 'B1', 'B2')
METHODS = LOCAL + ('G_RMSD', 'G_F')
BEST_SET = ('M1', 'M2', 'K', 'A', 'F', 'B1', 'B2')       # protocol §8
TARGETS = ('phipsi', 'pp')
CONTRASTS = {'H1_A_minus_M2': ('A', 'M2'), 'H1_A_minus_K': ('A', 'K'), 'H2_F_minus_A': ('F', 'A'),
             'H2_F_minus_Flock': ('F', 'F_lock'), 'H2_Flock_minus_A': ('F_lock', 'A'),
             'H4_GF_minus_F': ('G_F', 'F'), 'H4_GRMSD_minus_B2': ('G_RMSD', 'B2'), 'H4_F_minus_B2': ('F', 'B2'),
             'H5_B2_minus_A': ('B2', 'A'), 'H5_B1_minus_A': ('B1', 'A'), 'M1_minus_M2': ('M1', 'M2')}
QBLOCK, QBLOCK_SENS, WBLOCK = 100, 250, 20
N_BOOT, SEED, ALPHA, RHO, EPS_G = 2000, 20261012, 0.05, 0.10, 1e-6


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def tau_int_fn():
    spec = importlib.util.spec_from_file_location('adk_lib', HERE.parent / 'adk_feasibility' / 'adk_lib.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.tau_int


def block_sums(X, qb, wb):
    """X (n_windows, n_queries) with NaN = undefined → sums and counts over (window block, query block)."""
    nw, nq = X.shape
    wi = np.arange(nw) // wb
    qi = np.arange(nq) // qb
    S = np.zeros((wi.max() + 1, qi.max() + 1))
    C = np.zeros_like(S)
    ok = np.isfinite(X)
    np.add.at(S, (wi[:, None].repeat(nq, 1)[ok], qi[None, :].repeat(nw, 0)[ok]), X[ok])
    np.add.at(C, (wi[:, None].repeat(nq, 1)[ok], qi[None, :].repeat(nw, 0)[ok]), 1.0)
    return S, C


def interval(v, invalid=None):
    v = np.asarray(v, float)
    inv = np.zeros(len(v), bool) if invalid is None else np.asarray(invalid, bool)
    lo = np.quantile(np.where(inv, -np.inf, v), ALPHA / 2, method='inverted_cdf')
    hi = np.quantile(np.where(inv, np.inf, v), 1 - ALPHA / 2, method='inverted_cdf')
    return [float(lo), float(hi)]


def bootstrap_cell(arr, qb, rng_seed):
    """Two-way block bootstrap of cell means for every method (paired: one resample for all)."""
    sums = {m: block_sums(arr[m], qb, WBLOCK) for m in METHODS}
    nwb, nqb = sums[METHODS[0]][0].shape
    rng = np.random.default_rng(rng_seed)
    means = {m: np.empty(N_BOOT) for m in METHODS}
    for r in range(N_BOOT):
        cw = np.bincount(rng.integers(0, nwb, nwb), minlength=nwb).astype(float)
        cq = np.bincount(rng.integers(0, nqb, nqb), minlength=nqb).astype(float)
        wgt = cw[:, None] * cq[None, :]
        for m in METHODS:
            S, C = sums[m]
            means[m][r] = (S * wgt).sum() / (C * wgt).sum()
    return means


def classify(lo, hi, assessable):
    if not assessable:
        return 'not_assessable'
    if hi <= RHO:
        return 'within_margin'
    if lo > RHO:
        return 'beyond_margin'
    return 'inconclusive'


def main():
    summ = json.loads(SUMMARY.read_text())
    npz_path = REPO / summ['_meta']['per_window_query_nr10_npz']
    if sha256(npz_path) != summ['_meta']['npz_sha256']:
        raise SystemExit('STOP: per-window array file does not match the summary')
    Z = np.load(npz_path)
    tau_int = tau_int_fn()
    out = {'_meta': {'producer': 'analyse_dev.py', 'label': 'THIRD-PARTY MD DERIVATIVE (not SYNTHETIC)',
                     'summary_sha256': sha256(SUMMARY), 'npz_sha256': summ['_meta']['npz_sha256'],
                     'bootstrap': {'query_block_frames': QBLOCK, 'query_block_sensitivity': QBLOCK_SENS,
                                   'window_block_starts': WBLOCK, 'n_boot': N_BOOT, 'seed': SEED, 'alpha': ALPHA,
                                   'interval': 'percentile, numpy inverted_cdf'},
                     'rho': RHO, 'eps_G': EPS_G}, 'cells': {}}
    cell_i = 0
    for label in ('open_start', 'closed_start'):
        for W in (5, 8):
            for t in TARGETS:
                cell_i += 1
                arr = {m: Z[f'{label}__W{W}__{m}__{t}'].astype(np.float64) for m in METHODS}
                cell = {'n_windows': int(arr['A'].shape[0]), 'n_queries': int(arr['A'].shape[1]),
                        'undefined_nr_units': {m: int(np.sum(~np.isfinite(arr[m]))) for m in METHODS}}
                point = {m: float(np.nanmean(np.nanmean(arr[m], axis=1))) for m in METHODS}
                cell['mean_nr10'] = point
                cell['median_of_window_means'] = {m: float(np.nanmedian(np.nanmean(arr[m], axis=1)))
                                                  for m in METHODS}
                cell['gain_over_random'] = {m: 1.0 - point[m] for m in METHODS}
                boots = {}
                for qb in (QBLOCK, QBLOCK_SENS):
                    boots[qb] = bootstrap_cell(arr, qb, SEED + 1000 * cell_i + qb)
                cell['method_intervals'] = {m: {f'qblock{qb}': interval(boots[qb][m]) for qb in boots}
                                            for m in METHODS}
                con = {}
                for name, (x, y) in CONTRASTS.items():
                    d_point = point[x] - point[y]
                    unit = np.nanmean(arr[x] - arr[y], axis=0)            # per-query, window-averaged
                    e = {'delta': d_point, 'tau_int_frames_per_query_series': float(tau_int(unit[np.isfinite(unit)]))}
                    for qb in boots:
                        e[f'interval_qblock{qb}'] = interval(boots[qb][x] - boots[qb][y])
                    gov = QBLOCK_SENS if e['tau_int_frames_per_query_series'] > QBLOCK else QBLOCK
                    lo, hi = e[f'interval_qblock{gov}']
                    e['governing_qblock'] = gov
                    e['statement'] = 'below_0' if hi < 0 else ('above_0' if lo > 0 else 'unresolved')
                    con[name] = e
                cell['contrasts'] = con
                # development classification (primary target only), protocol §8, Amendment-1 estimand
                if t == 'phipsi':
                    cls = {}
                    for qb in boots:
                        b = boots[qb]
                        best_r = np.array([min(BEST_SET, key=lambda m: b[m][r]) for r in range(N_BOOT)])
                        gb = 1.0 - np.array([b[best_r[r]][r] for r in range(N_BOOT)])
                        inv = gb <= EPS_G
                        best_obs = min(BEST_SET, key=lambda m: point[m])
                        g_obs = 1.0 - point[best_obs]
                        for m in ('A', 'F'):
                            Lr = (b[m] - (1.0 - gb)) / np.where(inv, 1.0, gb)
                            lo, hi = interval(Lr, inv)
                            cls.setdefault(m, {})[f'qblock{qb}'] = {
                                'L': (point[m] - point[best_obs]) / g_obs if g_obs > EPS_G else None,
                                'L_interval': [lo, hi], 'invalid_fraction': float(inv.mean()),
                                'class': classify(lo, hi, g_obs > EPS_G)}
                        cls['best_observed'] = best_obs
                        cls[f'best_counts_qblock{qb}'] = {m: int(np.sum(best_r == m)) for m in BEST_SET}
                    cell['development_classification'] = cls
                out['cells'][f'{label}/W{W}/{t}'] = cell
                print(label, W, t, {m: round(point[m], 3) for m in METHODS})
    # overall development outcome per method (primary cells, governing block = 100 unless tau_int rule)
    overall = {}
    for m in ('A', 'F'):
        per = {}
        for qb in (QBLOCK, QBLOCK_SENS):
            cl = [out['cells'][f'{lab}/W{W}/phipsi']['development_classification'][m][f'qblock{qb}']['class']
                  for lab in ('open_start', 'closed_start') for W in (5, 8)]
            if any(c == 'not_assessable' for c in cl):
                o = 'insufficient_evidence'
            elif all(c == 'within_margin' for c in cl):
                o = 'retain_as_heldout_candidate'
            elif all(c == 'beyond_margin' for c in cl):
                o = 'reject_for_local_same_sequence_task'
            elif any(c == 'within_margin' for c in cl) and any(c == 'beyond_margin' for c in cl):
                o = 'restrict'
            else:
                o = 'inconclusive'
            per[f'qblock{qb}'] = {'cells': cl, 'outcome': o}
        overall[m] = per
    out['development_outcome'] = overall
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(overall, indent=1))


if __name__ == '__main__':
    sys.exit(main())
