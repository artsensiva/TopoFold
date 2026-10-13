#!/usr/bin/env python3
"""LE-1 analysis as amended (docs/LEGACY_ENGINE_EVALUATION_AMENDMENT_1.md). Supersedes analyse_dev.py, which is
preserved unmodified and was never run on AdK data.

Equal-window-weighted point estimates and a two-way block bootstrap of the *same* statistic (§1); paired
contrasts on per-unit differences; both block lengths with the τ_int statement rule (§4); H4 wording (§2);
verification gates (§5). Third-party MD derivative.
Run: ../wpd_local_retrieval/build/venv/bin/python -B analyse_dev_a1.py
"""
import hashlib
import importlib.util
import json
import sys
import warnings
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RES = HERE / 'results'
SUMMARY = RES / 'le_adk_THIRD_PARTY_MD_DERIVATIVE_summary.json'
OUT = RES / 'le_adk_THIRD_PARTY_MD_DERIVATIVE_analysis_a1.json'
LOCAL = ('M1', 'M2', 'K', 'A', 'F', 'F_lock', 'B1', 'B2')
METHODS = LOCAL + ('G_RMSD', 'G_F')
BEST_SET = ('M1', 'M2', 'K', 'A', 'F', 'B1', 'B2')
TARGETS = ('phipsi', 'pp')
CONTRASTS = {'H1_A_minus_M2': ('A', 'M2'), 'H1_A_minus_K': ('A', 'K'),
             'overall_F_minus_A': ('F', 'A'), 'coupling_F_minus_Flock': ('F', 'F_lock'),
             'aggregation_boundary_Flock_minus_A': ('F_lock', 'A'),
             'H4_GF_minus_F': ('G_F', 'F'), 'H4_GRMSD_minus_B2': ('G_RMSD', 'B2'),
             'descriptive_F_minus_B2': ('F', 'B2'), 'H5_B2_minus_A': ('B2', 'A'), 'H5_B1_minus_A': ('B1', 'A'),
             'M1_minus_M2': ('M1', 'M2')}
QBLOCKS, WBLOCK = (100, 250), 20
N_BOOT, SEED, ALPHA, RHO, EPS_G = 2000, 20261012, 0.05, 0.10, 1e-6
TAU_SHORT, TAU_LONG = 100, 250
TOL_FRECHET, TOL_REF, MAX_UNDEF_FRAC = 1e-9, 1e-9, 0.01


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def tau_int_fn():
    spec = importlib.util.spec_from_file_location('adk_lib', HERE.parent / 'adk_feasibility' / 'adk_lib.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.tau_int


def window_mean_stat(X):
    """θ = mean over windows with ≥1 defined unit of the window's mean over defined units (§1)."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)      # windows with no defined unit → NaN, excluded
        wm = np.nanmean(np.where(np.isfinite(X), X, np.nan), axis=1)
    return float(np.mean(wm[np.isfinite(wm)])) if np.isfinite(wm).any() else float('nan')


def window_block_sums(X, qb):
    """S[w, b], C[w, b]: sums and counts of defined units of window w in query block b."""
    nw, nq = X.shape
    bi = np.arange(nq) // qb
    ok = np.isfinite(X)
    S = np.zeros((nw, bi.max() + 1))
    C = np.zeros_like(S)
    for b in range(bi.max() + 1):
        cols = bi == b
        S[:, b] = np.where(ok[:, cols], X[:, cols], 0.0).sum(1)
        C[:, b] = ok[:, cols].sum(1)
    return S, C


def resample_weights(nw, nq, qb, n_boot, seed):
    """Query-block multiplicities CQ (B, n_qblocks) and per-window multiplicities CW (B, nw)."""
    rng = np.random.default_rng(seed)
    nqb = (nq - 1) // qb + 1
    nwb = (nw - 1) // WBLOCK + 1
    CQ = np.stack([np.bincount(rng.integers(0, nqb, nqb), minlength=nqb) for _ in range(n_boot)]).astype(float)
    CWB = np.stack([np.bincount(rng.integers(0, nwb, nwb), minlength=nwb) for _ in range(n_boot)]).astype(float)
    CW = CWB[:, np.arange(nw) // WBLOCK]
    return CQ, CW


def boot_stat(S, C, CQ, CW):
    """θ* for every resample: Σ c_w m*_w / Σ c_w over windows with defined m*_w."""
    num = CQ @ S.T
    den = CQ @ C.T
    with np.errstate(invalid='ignore', divide='ignore'):
        m = np.where(den > 0, num / np.where(den > 0, den, 1.0), np.nan)
    ok = np.isfinite(m)
    wsum = np.where(ok, CW, 0.0).sum(1)
    with np.errstate(invalid='ignore'):
        return np.where(wsum > 0, np.where(ok, CW * m, 0.0).sum(1) / np.where(wsum > 0, wsum, 1.0), np.nan)


def interval(v, invalid=None):
    v = np.asarray(v, float)
    inv = ~np.isfinite(v) if invalid is None else (np.asarray(invalid, bool) | ~np.isfinite(v))
    lo = np.quantile(np.where(inv, -np.inf, v), ALPHA / 2, method='inverted_cdf')
    hi = np.quantile(np.where(inv, np.inf, v), 1 - ALPHA / 2, method='inverted_cdf')
    return [float(lo), float(hi)], float(inv.mean())


def classify(lo, hi, assessable):
    if not assessable:
        return 'not_assessable'
    if hi <= RHO:
        return 'within_margin'
    if lo > RHO:
        return 'beyond_margin'
    return 'inconclusive'


def statement(lo, hi):
    return 'below_0' if hi < 0 else ('above_0' if lo > 0 else 'unresolved')


def governing(tau):
    if tau is None or not np.isfinite(tau) or tau > TAU_LONG:
        return None
    return TAU_SHORT if tau <= TAU_SHORT else TAU_LONG


def outcome(classes):
    if any(c == 'not_assessable' for c in classes):
        return 'insufficient_evidence'
    if all(c == 'within_margin' for c in classes):
        return 'retain_as_heldout_candidate'
    if all(c == 'beyond_margin' for c in classes):
        return 'reject_for_local_same_sequence_task'
    if any(c == 'within_margin' for c in classes) and any(c == 'beyond_margin' for c in classes):
        return 'restrict'
    return 'inconclusive'


def per_query_series(diff):
    """Window-averaged per-query series (time order) of a per-unit difference."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        s = np.nanmean(np.where(np.isfinite(diff), diff, np.nan), axis=0)
    return s[np.isfinite(s)]


def main():
    summ = json.loads(SUMMARY.read_text())
    npz_path = REPO / summ['_meta']['per_window_query_nr10_npz']
    if sha256(npz_path) != summ['_meta']['npz_sha256']:
        raise SystemExit('STOP: per-window array file does not match the summary')
    Z = np.load(npz_path)
    tau_int = tau_int_fn()
    out = {'_meta': {'producer': 'analyse_dev_a1.py', 'amendment': 'docs/LEGACY_ENGINE_EVALUATION_AMENDMENT_1.md',
                     'label': 'THIRD-PARTY MD DERIVATIVE (not SYNTHETIC)', 'summary_sha256': sha256(SUMMARY),
                     'npz_sha256': summ['_meta']['npz_sha256'],
                     'bootstrap': {'query_blocks': QBLOCKS, 'window_block_starts': WBLOCK, 'n_boot': N_BOOT,
                                   'seed': SEED, 'alpha': ALPHA, 'statistic': 'equal-window-weighted mean'},
                     'rho': RHO, 'eps_G': EPS_G}, 'cells': {}}
    gates = {'reference_check': all(t['reference_check']['frames_checked'] == t['n_frames']
                                    and t['reference_check']['max_tau_circ_diff'] <= TOL_REF
                                    and t['reference_check']['max_kappa_diff'] <= TOL_REF
                                    for t in summ['trajectories'].values()),
             'frechet_rust_numpy': max(x['frechet_rust_numpy_max_abs_diff'] for x in summ['windows']) <= TOL_FRECHET,
             'binding_failures_zero': all(t['binding_failure_frames'] == 0 for t in summ['trajectories'].values())}
    nonfinite = {m: sum(x['methods'][m][t]['nonfinite_distances'] for x in summ['windows'] if 'methods' in x
                        for t in ('phipsi',)) for m in METHODS}
    failures = {m: sum(x['methods'][m][t]['failures'] for x in summ['windows'] if 'methods' in x for t in TARGETS)
                for m in METHODS}
    gates['nonfinite_zero_best_set'] = all(nonfinite[m] == 0 for m in BEST_SET)
    gates['regret_failures_zero'] = all(v == 0 for v in failures.values())
    out['verification'] = {'gates': gates, 'nonfinite_distances_phipsi': nonfinite, 'regret_failures': failures}
    cell_i = 0
    undef_ok = True
    for label in ('open_start', 'closed_start'):
        ws = [x for x in summ['windows'] if x['trajectory'] == label]
        for W in (5, 8):
            wW = [x for x in ws if x['W'] == W]
            nvq = [x['n_query_valid'] for x in wW]
            for t in TARGETS:
                cell_i += 1
                arr = {m: Z[f'{label}__W{W}__{m}__{t}'].astype(np.float64) for m in METHODS}
                nw, nq = arr['A'].shape
                defined = {m: np.isfinite(arr[m]).sum(1) for m in METHODS}
                cell = {'n_windows': int(nw), 'n_query_frames': int(nq),
                        'valid_query_counts': {'min': int(min(nvq)), 'max': int(max(nvq)),
                                               'n_distinct': len(set(nvq))},
                        'defined_units_per_window': {m: {'min': int(defined[m].min()), 'max': int(defined[m].max())}
                                                     for m in METHODS},
                        'undefined_units': {m: int((~np.isfinite(arr[m])).sum()) for m in METHODS},
                        'windows_without_defined_units': {m: int((defined[m] == 0).sum()) for m in METHODS}}
                if t == 'phipsi':
                    undef_ok &= all(cell['undefined_units'][m] < MAX_UNDEF_FRAC * nw * nq for m in BEST_SET)
                point = {m: window_mean_stat(arr[m]) for m in METHODS}
                cell['theta'] = point
                cell['gain_over_random'] = {m: 1.0 - point[m] for m in METHODS}
                W8 = {}
                for qb in QBLOCKS:
                    CQ, CW = resample_weights(nw, nq, qb, N_BOOT, SEED + 1000 * cell_i + qb)
                    W8[qb] = (CQ, CW)
                bm = {qb: {m: boot_stat(*window_block_sums(arr[m], qb), *W8[qb]) for m in METHODS} for qb in QBLOCKS}
                cell['method_intervals'] = {m: {f'qblock{qb}': interval(bm[qb][m])[0] for qb in QBLOCKS}
                                            for m in METHODS}
                con = {}
                for name, (x, y) in CONTRASTS.items():
                    diff = arr[x] - arr[y]
                    series = per_query_series(diff)
                    tau = tau_int(series) if series.size > 2 and np.ptp(series) > 0 else None
                    e = {'delta_paired': window_mean_stat(diff), 'delta_of_thetas': point[x] - point[y],
                         'tau_int_frames': tau}
                    for qb in QBLOCKS:
                        e[f'interval_qblock{qb}'] = interval(boot_stat(*window_block_sums(diff, qb), *W8[qb]))[0]
                    g = governing(tau)
                    e['governing_qblock'] = g
                    e['statement'] = ('dependence_unresolved' if g is None
                                      else statement(*e[f'interval_qblock{g}']))
                    e['statement_qblock100'] = statement(*e['interval_qblock100'])
                    e['statement_qblock250'] = statement(*e['interval_qblock250'])
                    con[name] = e
                cell['contrasts'] = con
                if t == 'phipsi':
                    best_obs = min(BEST_SET, key=lambda m: point[m])
                    g_obs = 1.0 - point[best_obs]
                    cls = {'best_observed': best_obs}
                    for qb in QBLOCKS:
                        b = bm[qb]
                        B = np.stack([b[m] for m in BEST_SET])
                        bi = np.argmin(B, axis=0)
                        gb = 1.0 - B[bi, np.arange(N_BOOT)]
                        inv = gb <= EPS_G
                        cls[f'best_counts_qblock{qb}'] = {m: int(np.sum(bi == j)) for j, m in enumerate(BEST_SET)}
                        for m in ('A', 'F'):
                            Lr = (b[m] - (1.0 - gb)) / np.where(inv, 1.0, gb)
                            (lo, hi), fi = interval(Lr, inv)
                            cls.setdefault(m, {})[f'qblock{qb}'] = {
                                'L': (point[m] - point[best_obs]) / g_obs if g_obs > EPS_G else None,
                                'L_interval': [lo, hi], 'invalid_fraction': fi,
                                'class': classify(lo, hi, g_obs > EPS_G)}
                    for m in ('A', 'F'):
                        s = per_query_series(arr[m] - arr[best_obs])
                        cls[m]['tau_int_vs_best_frames'] = (tau_int(s) if s.size > 2 and np.ptp(s) > 0 else None)
                    cell['development_classification'] = cls
                out['cells'][f'{label}/W{W}/{t}'] = cell
                print(label, W, t, {m: round(point[m], 3) for m in METHODS}, flush=True)
    gates['undefined_below_1pct'] = bool(undef_ok)
    overall = {}
    for m in ('A', 'F'):
        per = {}
        for qb in QBLOCKS:
            cl = [out['cells'][f'{lab}/W{W}/phipsi']['development_classification'][m][f'qblock{qb}']['class']
                  for lab in ('open_start', 'closed_start') for W in (5, 8)]
            per[f'qblock{qb}'] = {'cells': cl, 'outcome': outcome(cl)}
        taus = [out['cells'][f'{lab}/W{W}/phipsi']['development_classification'][m]['tau_int_vs_best_frames']
                for lab in ('open_start', 'closed_start') for W in (5, 8)]
        reasons = [k for k, v in gates.items() if not v]
        if per['qblock100']['outcome'] != per['qblock250']['outcome']:
            reasons.append('block_lengths_disagree')
        if any(tt is None or tt > TAU_LONG for tt in taus):
            reasons.append('dependence_unresolved')
        per['tau_int_vs_best_frames'] = taus
        per['gate_failures'] = reasons
        per['final'] = (per['qblock100']['outcome'] if not reasons
                        else f"insufficient_evidence ({', '.join(reasons)})")
        overall[m] = per
    out['development_outcome'] = overall
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({m: overall[m]['final'] for m in overall}, indent=1))


if __name__ == '__main__':
    sys.exit(main())
