#!/usr/bin/env python3
"""
M3-A2g data verification: format, bookkeeping and derivative-freshness checks.

SYNTHETIC: every input is an output of the synthetic M3-A2g experiments
(simulated product-von-Mises / wrapped-normal draws). Nothing here is MD data.

Inputs (preserved originals, never modified by this script):
  results/m3_a2g_*.csv (10 files) and results/m3_a2g_trajectories.txt
Derivative (regenerated only on request):
  results/m3_a2g_data_verification.json

Default mode is READ-ONLY. It recomputes every check from the input bytes,
compares the stored derivative JSON against the fresh result, and exits
nonzero if an integrity check fails or the stored JSON is stale.
`--write` regenerates the JSON, and only after the integrity checks pass;
the comparison against the previously stored file is reported first.

What is checked (expected designs derived from the generator scripts):
  - headers: required columns present, no unexpected or duplicate columns
    (column order is not enforced; no consumer depends on it)
  - every field parses as its type; empty required values are errors;
    non-finite values are errors except documented `nan` sentinels, which
    must appear exactly where the generator writes them
  - value domains (selected_K in {1, 2}, counts within budgets, ...)
  - configuration/replicate coverage against the generator loops:
    missing, duplicate and unexpected keys
  - selected_K against the implemented Delta-BIC rule (K=2 iff dBIC < 0;
    equality selects K=1), for every selection pair a file records
  - trajectory blocks: declarations, per-block counts, iteration sequence,
    field parsing, malformed lines

What is NOT checked:
  - scientific validity of the optimizer or of any diagnostic. A consistent
    selection bookkeeping entry is not a validated selection; non-converged
    fits are reported separately and are not validated by these checks.
  - random seeds (they are not recorded in the CSVs)
  - any correspondence between trajectory blocks and CSV rows (none is
    recorded, so none is asserted)
"""
import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
BASE_DIR = SCRIPT_DIR.parent
RESULTS_DIR = BASE_DIR / "results"
JSON_NAME = "m3_a2g_data_verification.json"
TRAJECTORY_NAME = "m3_a2g_trajectories.txt"
BASELINE_REV = "394ce0ccf8c46d3f37a2bf0e9a81e255886de5c0"

MAX_LISTED_ERRORS = 50


# ============================================================================
# DECISION RULE
# ============================================================================

def expected_selected_k(delta_bic):
    """Delta-BIC rule as implemented by the generators.

    product_vm_fast_accurate.select_K_via_bic_corrected: selects_K2 = (delta_bic < 0)
    m3_a2g_stopping_criterion_validation.continuation_test:
        selected_K = 2 if delta_bic < 0 else 1
    delta_bic = BIC_K2 - BIC_K1; equality selects K=1.
    """
    return 2 if delta_bic < 0 else 1


# ============================================================================
# DATASET SPECIFICATIONS (derived from the generator scripts)
# ============================================================================

def _nonneg(v):
    return v >= 0


def _positive(v):
    return v > 0


def _in(*allowed):
    return lambda v: v in allowed


def _between(lo, hi):
    return lambda v: lo <= v <= hi


def _close(a, b):
    return math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-9)


def _row_ninit(r):
    errs = []
    if r['n_converged'] > r['n_init']:
        errs.append(f"n_converged={r['n_converged']} exceeds n_init={r['n_init']}")
    if not _close(r['convergence_rate'], r['n_converged'] / r['n_init']):
        errs.append("convergence_rate != n_converged / n_init")
    if r['min_loglik'] > r['max_loglik']:
        errs.append("min_loglik > max_loglik")
    if not _close(r['spread_range'], r['max_loglik'] - r['min_loglik']):
        errs.append("spread_range != max_loglik - min_loglik")
    return errs


def _row_stopping(r):
    errs = []
    if r['K_changed'] != (r['selected_K_cont'] != r['selected_K_conv']):
        errs.append("K_changed inconsistent with selected_K_conv/selected_K_cont")
    if not _close(r['loglik_gain'], r['loglik_cont'] - r['loglik_conv']):
        errs.append("loglik_gain != loglik_cont - loglik_conv")
    if not _close(r['bic_change'], r['delta_bic_cont'] - r['delta_bic_conv']):
        errs.append("bic_change != delta_bic_cont - delta_bic_conv")
    # continuation_test: convergence_iter = 2000 if not converged else 2000 - 100
    expected_iter = 1900 if r['converged'] else 2000
    if r['convergence_iter'] != expected_iter:
        errs.append(f"convergence_iter={r['convergence_iter']} != {expected_iter} "
                    f"(generator approximation for converged={r['converged']})")
    return errs


def _row_duplicate(r):
    errs = []
    if not _close(r['final_separation_deg'], r['final_separation_rad'] * 180 / math.pi):
        errs.append("final_separation_deg != final_separation_rad * 180/pi")
    return errs


def _row_iters_within(budget_col, iters_col):
    def check(r):
        if r[iters_col] > r[budget_col]:
            return [f"{iters_col}={r[iters_col]} exceeds {budget_col}={r[budget_col]}"]
        return []
    return check


def _keys_tolerance():
    keys = set()
    for case in ['exact_easy', 'wrapped_0.9']:
        for tol in [1e-4, 1e-5, 1e-6, 1e-7]:
            keys.add((case, 'absolute', tol))
        for tol in [1e-5, 1e-6, 1e-7, 1e-8]:
            keys.add((case, 'relative', tol))
    return keys


DATASETS = {
    'exact_family': {
        'file': 'm3_a2g_exact_family_convergence.csv',
        'generator': 'm3_a2g_exact_k2_convergence.py::audit_exact_family_convergence',
        'fields': {'q': 'int', 'kappa': 'int', 'F': 'int', 'rep': 'int', 'converged': 'bool',
                   'selected_K': 'int', 'delta_bic': 'float', 'loglik': 'float'},
        'domains': {'q': _positive, 'kappa': _nonneg, 'F': _positive, 'rep': _nonneg,
                    'selected_K': _in(1, 2)},
        'key': ('q', 'kappa', 'F', 'rep'),
        'expected_keys': set(product([1, 2, 5], [1, 5, 20], [100, 500, 1000, 2500], range(10))),
        'selection_pairs': [('selected_K', 'delta_bic')],
    },
    'k2_convergence': {
        'file': 'm3_a2g_k2_convergence.csv',
        'generator': 'm3_a2g_exact_k2_convergence.py::audit_k2_convergence',
        'fields': {'separation_deg': 'int', 'occupancy_1': 'float', 'occupancy_2': 'float',
                   'kappa': 'int', 'F': 'int', 'rep': 'int', 'converged': 'bool',
                   'selected_K': 'int', 'delta_bic': 'float', 'loglik': 'float'},
        'domains': {'separation_deg': _nonneg, 'occupancy_1': _between(0, 1),
                    'occupancy_2': _between(0, 1), 'kappa': _nonneg, 'F': _positive,
                    'rep': _nonneg, 'selected_K': _in(1, 2)},
        'key': ('separation_deg', 'occupancy_1', 'occupancy_2', 'kappa', 'F', 'rep'),
        'expected_keys': {(s, o1, o2, k, F, rep)
                          for s, (o1, o2), k, F, rep in product(
                              [10, 30, 60, 90], [(0.5, 0.5), (0.8, 0.2)], [5, 20], [250, 1000],
                              range(10))},
        'selection_pairs': [('selected_K', 'delta_bic')],
    },
    'stress_regime': {
        'file': 'm3_a2g_stress_regime.csv',
        'generator': 'm3_a2g_convergence_grid.py::test_wrapped_normal_stress',
        'fields': {'sigma': 'float', 'F': 'int', 'rep': 'int', 'converged': 'bool',
                   'selected_K': 'int', 'loglik': 'float', 'delta_bic': 'float',
                   'max_kappa': 'float'},
        'domains': {'sigma': _positive, 'F': _positive, 'rep': _nonneg,
                    'selected_K': _in(1, 2), 'max_kappa': _nonneg},
        'key': ('sigma', 'F', 'rep'),
        'expected_keys': set(product([0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2],
                                     [250, 500, 1000, 2500], range(10))),
        'selection_pairs': [('selected_K', 'delta_bic')],
    },
    'extended_maxiter': {
        'file': 'm3_a2g_extended_maxiter.csv',
        'generator': 'm3_a2g_extended_maxiter.py::test_extended_maxiter',
        'fields': {'sigma': 'float', 'F': 'int', 'max_iter': 'int', 'rep': 'int',
                   'converged': 'bool', 'selected_K': 'int', 'loglik': 'float',
                   'delta_bic': 'float', 'max_kappa': 'float', 'runtime_sec': 'float'},
        'domains': {'sigma': _positive, 'F': _positive, 'max_iter': _positive,
                    'rep': _nonneg, 'selected_K': _in(1, 2), 'max_kappa': _nonneg,
                    'runtime_sec': _nonneg},
        'key': ('sigma', 'F', 'max_iter', 'rep'),
        'expected_keys': {(s, F, m, rep)
                          for (s, F), m, rep in product(
                              [(0.6, 1000), (0.6, 2500), (0.7, 1000), (0.7, 2500),
                               (1.0, 1000), (1.0, 2500)], [500, 1000, 2000], range(10))},
        'selection_pairs': [('selected_K', 'delta_bic')],
    },
    'ninit_sensitivity': {
        'file': 'm3_a2g_ninit_sensitivity.csv',
        'generator': 'm3_a2g_ninit_sensitivity.py (INVALID diagnostic: inner reseeding defect)',
        'fields': {'case': 'str', 'n_init': 'int', 'best_loglik': 'float', 'min_loglik': 'float',
                   'max_loglik': 'float', 'spread_range': 'float', 'spread_std': 'float',
                   'n_converged': 'int', 'convergence_rate': 'float', 'best_n_iters': 'int',
                   'best_converged': 'bool', 'runtime_sec': 'float'},
        'domains': {'n_init': _positive, 'spread_range': _nonneg, 'spread_std': _nonneg,
                    'n_converged': _nonneg, 'convergence_rate': _between(0, 1),
                    'best_n_iters': _between(1, 1000), 'runtime_sec': _nonneg},
        'row_checks': [_row_ninit],
        'key': ('case', 'n_init'),
        'expected_keys': set(product(['exact_k1', 'wrapped_1.0', 'k2_close', 'k2_clear'],
                                     [1, 3, 5, 10, 20])),
        'selection_pairs': [],
    },
    'stopping_criterion_validation': {
        'file': 'm3_a2g_stopping_criterion_validation.csv',
        'generator': 'm3_a2g_stopping_criterion_validation.py (INVALID diagnostic)',
        'fields': {'converged': 'bool', 'convergence_iter': 'int', 'loglik_conv': 'float',
                   'loglik_cont': 'float', 'loglik_gain': 'float', 'weight_change': 'float',
                   'max_center_change_rad': 'float', 'kappa_change': 'float',
                   'delta_bic_conv': 'float', 'delta_bic_cont': 'float', 'bic_change': 'float',
                   'selected_K_conv': 'int', 'selected_K_cont': 'int', 'K_changed': 'bool',
                   'case': 'str', 'tolerance': 'str'},
        'domains': {'weight_change': _nonneg, 'max_center_change_rad': _between(0, math.pi),
                    'kappa_change': _nonneg, 'selected_K_conv': _in(1, 2),
                    'selected_K_cont': _in(1, 2)},
        'row_checks': [_row_stopping],
        'key': ('case', 'tolerance'),
        'expected_keys': set(product(['exact_k1', 'wrapped_0.9', 'wrapped_1.0', 'k2_close',
                                      'k2_clear'],
                                     ['abs_1e-6', 'rel_1e-6', 'rel_1e-7', 'rel_1e-8'])),
        'selection_pairs': [('selected_K_conv', 'delta_bic_conv'),
                            ('selected_K_cont', 'delta_bic_cont')],
    },
    'tolerance_comparison': {
        'file': 'm3_a2g_tolerance_comparison.csv',
        'generator': 'm3_a2g_convergence_grid.py::test_tolerance_criteria',
        'fields': {'case': 'str', 'criterion': 'str', 'tolerance': 'float', 'converged': 'bool',
                   'n_iters': 'int', 'final_loglik': 'float'},
        'domains': {'criterion': _in('absolute', 'relative'), 'tolerance': _positive,
                    'n_iters': _between(1, 1000)},
        'key': ('case', 'criterion', 'tolerance'),
        'expected_keys': _keys_tolerance(),
        'selection_pairs': [],
    },
    'max_iter_sensitivity': {
        'file': 'm3_a2g_max_iter_sensitivity.csv',
        'generator': 'm3_a2g_convergence_grid.py::test_max_iter_sensitivity',
        'fields': {'case': 'str', 'max_iter': 'int', 'converged': 'bool', 'n_iters': 'int',
                   'final_loglik': 'float', 'improvement_after_100': 'float',
                   'final_delta_abs': 'float', 'final_delta_rel': 'float', 'max_kappa': 'float'},
        'domains': {'max_iter': _positive, 'n_iters': _positive, 'max_kappa': _nonneg},
        # Generator writes nan exactly when the history is too short.
        'nan_sentinels': {
            'improvement_after_100': ('n_iters <= 100', lambda r: r['n_iters'] <= 100),
            'final_delta_abs': ('n_iters <= 1', lambda r: r['n_iters'] <= 1),
            'final_delta_rel': ('n_iters <= 1', lambda r: r['n_iters'] <= 1),
        },
        'row_checks': [_row_iters_within('max_iter', 'n_iters')],
        'key': ('case', 'max_iter'),
        'expected_keys': set(product(['exact_easy', 'wrapped_0.9', 'wrapped_1.0', 'k2_clear'],
                                     [100, 250, 500, 1000])),
        'selection_pairs': [],
    },
    'initialization_comparison': {
        'file': 'm3_a2g_initialization_comparison.csv',
        'generator': 'm3_a2g_convergence_grid.py::test_initialization_strategies',
        'fields': {'case': 'str', 'strategy': 'str', 'best_loglik': 'float',
                   'loglik_spread': 'float', 'n_converged': 'int', 'best_n_iters': 'int',
                   'best_converged': 'bool'},
        'domains': {'loglik_spread': _nonneg, 'n_converged': _between(0, 5),
                    'best_n_iters': _between(1, 500)},
        'key': ('case', 'strategy'),
        'expected_keys': set(product(['exact_easy', 'wrapped_1.0', 'k2_close'],
                                     ['split', 'random', 'farthest', 'embedding',
                                      'perturb_k1'])),
        'selection_pairs': [],
    },
    'duplicate_components': {
        'file': 'm3_a2g_duplicate_components.csv',
        'generator': 'm3_a2g_duplicate_component_investigation.py',
        'fields': {'case': 'str', 'converged': 'bool', 'n_iters': 'int', 'final_loglik': 'float',
                   'final_separation_rad': 'float', 'final_separation_deg': 'float',
                   'final_kappa_diff': 'float', 'final_weight_ratio': 'float',
                   'final_min_weight': 'float', 'final_max_kappa': 'float',
                   'effective_n_components': 'float'},
        'domains': {'n_iters': _between(1, 1000), 'final_separation_rad': _nonneg,
                    'final_kappa_diff': _nonneg, 'final_weight_ratio': _between(0, 1),
                    'final_min_weight': _between(0, 0.5), 'final_max_kappa': _nonneg,
                    # 1 / sum(pi_k^2) for K=2 lies in [1, 2]
                    'effective_n_components': _between(1, 2)},
        'row_checks': [_row_duplicate],
        'key': ('case',),
        'expected_keys': {('k1_kappa1',), ('k1_kappa5',), ('k1_kappa20',), ('wrapped_0.6',),
                          ('wrapped_0.8',), ('wrapped_0.9',), ('wrapped_1.0',),
                          ('wrapped_1.2',)},
        'selection_pairs': [],
    },
}

TRAJECTORY_HEADER = "iter,loglik,delta_abs,delta_rel,max_kappa,min_weight,eff_n_comp"
# m3_a2g_convergence_audit.py writes these cases in this order.
EXPECTED_TRAJECTORY_CASES = ['exact_easy', 'wrapped_0.9', 'wrapped_1.0', 'k2_clear', 'k2_close']
# loglik is written with 6 decimals; decreases smaller than this are not resolvable.
LOGLIK_DECREASE_TOL = 1e-6


# ============================================================================
# HELPERS
# ============================================================================

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def compute_input_hashes(results_dir=RESULTS_DIR):
    """SHA-256 of the 10 CSV inputs and the trajectory file. Missing files raise."""
    results_dir = Path(results_dir)
    names = sorted([spec['file'] for spec in DATASETS.values()] + [TRAJECTORY_NAME])
    return {name: sha256_file(results_dir / name) for name in names}


def _cap(errors):
    return {'count': len(errors), 'items': errors[:MAX_LISTED_ERRORS]}


def parse_value(raw, ftype):
    """Parse one CSV field. Returns (value, error). Non-finite floats are returned
    as values; the caller decides whether they are permitted."""
    if raw is None or raw.strip() == '':
        return None, 'missing required value'
    s = raw.strip()
    if ftype == 'str':
        return s, None
    if ftype == 'bool':
        if s in ('True', 'False'):
            return s == 'True', None
        return None, f"invalid bool {raw!r}"
    if ftype == 'int':
        try:
            return int(s), None
        except ValueError:
            return None, f"invalid int {raw!r}"
    if ftype == 'float':
        try:
            return float(s), None
        except ValueError:
            return None, f"invalid float {raw!r}"
    raise ValueError(f"unknown field type {ftype}")


# ============================================================================
# CSV VALIDATION
# ============================================================================

def validate_dataset(name, path=None, spec=None):
    """Validate one CSV against its generator-derived specification."""
    spec = spec or DATASETS[name]
    path = Path(path) if path else RESULTS_DIR / spec['file']
    fields = spec['fields']
    nan_sentinels = spec.get('nan_sentinels', {})

    with open(path, newline='') as f:
        raw_rows = list(csv.reader(f))

    if not raw_rows:
        header = []
        body = []
    else:
        header, body = raw_rows[0], raw_rows[1:]

    counts = Counter(header)
    header_report = {
        'columns': header,
        'missing_columns': sorted(set(fields) - set(header)),
        'unexpected_columns': sorted(set(header) - set(fields)),
        'duplicate_columns': sorted(c for c, n in counts.items() if n > 1),
    }
    header_report['schema_valid'] = not (header_report['missing_columns']
                                         or header_report['unexpected_columns']
                                         or header_report['duplicate_columns'])

    field_errors, row_errors, selection_errors = [], [], []
    parsed_rows = []
    key_counts = Counter()

    if header_report['schema_valid']:
        for i, raw in enumerate(body, start=2):  # line number in file
            if len(raw) != len(header):
                field_errors.append(f"line {i}: {len(raw)} fields, header has {len(header)}")
                continue
            raw_map = dict(zip(header, raw))
            row, ok = {}, True
            for col, ftype in fields.items():
                value, err = parse_value(raw_map[col], ftype)
                if err:
                    field_errors.append(f"line {i}: {col}: {err}")
                    ok = False
                row[col] = value
            if not ok:
                continue

            for col, ftype in fields.items():
                v = row[col]
                if ftype == 'float':
                    if col in nan_sentinels:
                        desc, cond = nan_sentinels[col]
                        if cond(row):
                            if not math.isnan(v):
                                field_errors.append(
                                    f"line {i}: {col}={v} but generator writes nan when {desc}")
                                ok = False
                            continue
                    if not math.isfinite(v):
                        field_errors.append(f"line {i}: {col}: non-finite value {v}")
                        ok = False
                        continue
                check = spec.get('domains', {}).get(col)
                if check is not None and not check(v):
                    field_errors.append(f"line {i}: {col}={v} outside allowed domain")
                    ok = False

            key_counts[tuple(row[k] for k in spec['key'])] += 1
            if not ok:
                continue

            for check in spec.get('row_checks', []):
                for err in check(row):
                    row_errors.append(f"line {i}: {err}")

            for k_col, d_col in spec['selection_pairs']:
                expected = expected_selected_k(row[d_col])
                if row[k_col] != expected:
                    selection_errors.append(
                        f"line {i}: {k_col}={row[k_col]} but {d_col}={row[d_col]!r} "
                        f"implies K={expected}")
            parsed_rows.append(row)

    present = set(key_counts)
    expected_keys = spec['expected_keys']
    coverage = {
        'expected': len(expected_keys),
        'present_unique': len(present),
        'missing': sorted(map(list, expected_keys - present), key=repr),
        'unexpected': sorted(map(list, present - expected_keys), key=repr),
        'duplicates': sorted(([list(k), n] for k, n in key_counts.items() if n > 1), key=repr),
    }
    coverage['complete'] = not (coverage['missing'] or coverage['unexpected']
                                or coverage['duplicates'])

    records_valid = header_report['schema_valid'] and not field_errors and not row_errors
    result = {
        'file': spec['file'],
        'generator': spec['generator'],
        'row_count': len(body),
        'header': header_report,
        'field_errors': _cap(field_errors),
        'row_consistency_errors': _cap(row_errors),
        'selection_bookkeeping_errors': _cap(selection_errors),
        'selection_pairs_checked': [list(p) for p in spec['selection_pairs']],
        'coverage': coverage,
        'records_valid': records_valid,
        'valid': records_valid and coverage['complete'] and not selection_errors,
    }
    return result, parsed_rows


# ============================================================================
# TRAJECTORY VALIDATION
# ============================================================================

def parse_trajectory_file(path=None, expected_cases=EXPECTED_TRAJECTORY_CASES,
                          decrease_tol=LOGLIK_DECREASE_TOL):
    """Strictly parse the trajectory file. Malformed lines are errors, never skipped."""
    path = Path(path) if path else RESULTS_DIR / TRAJECTORY_NAME
    errors = []
    blocks = []
    block = None
    state = 'preamble'  # preamble -> declarations -> records -> closed

    with open(path) as f:
        lines = f.read().split('\n')

    for ln, line in enumerate(lines, start=1):
        s = line.strip()
        if s.startswith('Case:'):
            block = {'case': s[len('Case:'):].strip(), 'line': ln, 'declared': {},
                     'header_seen': False, 'records': [], 'iter_lines': []}
            blocks.append(block)
            state = 'declarations'
            continue
        if state == 'preamble':
            if s and s != 'M3-A2g ITERATION TRAJECTORIES' and set(s) != {'='}:
                errors.append(f"line {ln}: unexpected content before first block: {s[:60]!r}")
            continue
        if state == 'declarations':
            if not s or set(s) == {'-'}:
                continue
            if s == TRAJECTORY_HEADER:
                block['header_seen'] = True
                state = 'records'
                continue
            if ':' in s:
                label, _, value = s.partition(':')
                label, value = label.strip(), value.strip()
                parsers = {'Converged': 'bool', 'Iterations': 'int', 'Final loglik': 'float'}
                if label in parsers:
                    if label in block['declared']:
                        errors.append(f"line {ln}: duplicate declaration {label!r}")
                    v, err = parse_value(value, parsers[label])
                    if err:
                        errors.append(f"line {ln}: {label}: {err}")
                    else:
                        block['declared'][label] = v
                    continue
            if s.startswith('iter,'):
                errors.append(f"line {ln}: unexpected record header {s!r}")
            else:
                errors.append(f"line {ln}: unexpected declaration line {s[:60]!r}")
            continue
        if state == 'records':
            if not s:
                state = 'closed'
                continue
            parts = s.split(',')
            if len(parts) != 7:
                errors.append(f"line {ln}: malformed record ({len(parts)} fields, expected 7)")
                block['records'].append(None)
                block['iter_lines'].append((parse_value(parts[0], 'int')[0], ln))
                continue
            vals = {}
            bad = False
            for col, raw, ftype in zip(
                    ['iter', 'loglik', 'delta_abs', 'delta_rel', 'max_kappa', 'min_weight',
                     'eff_n_comp'], parts,
                    ['int', 'float', 'float', 'float', 'float', 'float', 'float']):
                v, err = parse_value(raw, ftype)
                if err:
                    errors.append(f"line {ln}: {col}: {err}")
                    bad = True
                vals[col] = v
            block['iter_lines'].append((vals['iter'], ln))
            if bad:
                block['records'].append(None)
                continue
            for col in ['loglik', 'max_kappa', 'min_weight', 'eff_n_comp']:
                if not math.isfinite(vals[col]):
                    errors.append(f"line {ln}: {col}: non-finite value {vals[col]}")
                    bad = True
            # Generator writes nan deltas exactly at iteration 0.
            for col in ['delta_abs', 'delta_rel']:
                if vals['iter'] == 0 and not math.isnan(vals[col]):
                    errors.append(f"line {ln}: {col} should be nan at iteration 0")
                    bad = True
                elif vals['iter'] != 0 and not math.isfinite(vals[col]):
                    errors.append(f"line {ln}: {col}: non-finite value {vals[col]}")
                    bad = True
            if not bad:
                if vals['max_kappa'] < 0:
                    errors.append(f"line {ln}: max_kappa < 0")
                if not -1e-6 <= vals['min_weight'] <= 0.5 + 1e-6:
                    errors.append(f"line {ln}: min_weight outside [0, 0.5]")
                if not 1 - 1e-4 <= vals['eff_n_comp'] <= 2 + 1e-4:
                    errors.append(f"line {ln}: eff_n_comp outside [1, 2]")
            vals['line'] = ln
            block['records'].append(None if bad else vals)
            continue
        if state == 'closed' and s:
            errors.append(f"line {ln}: content after block records: {s[:60]!r}")

    decreases = []
    block_reports = []
    for b in blocks:
        recs = b['records']
        good = [r for r in recs if r is not None]
        name = b['case']
        for label in ['Converged', 'Iterations', 'Final loglik']:
            if label not in b['declared']:
                errors.append(f"block {name!r} (line {b['line']}): missing declaration {label!r}")
        if not b['header_seen']:
            errors.append(f"block {name!r} (line {b['line']}): missing record header")
        declared_n = b['declared'].get('Iterations')
        if declared_n is not None and declared_n != len(recs):
            errors.append(f"block {name!r}: declared Iterations={declared_n}, "
                          f"found {len(recs)} records")
        if not recs:
            errors.append(f"block {name!r}: no iteration records")
        iters = [it for it, _ in b['iter_lines']]
        # Writer emits history entries with iter = 0, 1, 2, ... in order.
        for idx, (it, ln) in enumerate(b['iter_lines']):
            if it != idx:
                errors.append(f"block {name!r}: iteration sequence broken at record {idx} "
                              f"(line {ln}: iter={it}, expected {idx})")
                break
        final = b['declared'].get('Final loglik')
        # Writer: Final loglik = best['loglik'], the loglik of the last history entry.
        if final is not None and good and len(good) == len(recs):
            if abs(final - good[-1]['loglik']) > 1e-6:
                errors.append(f"block {name!r}: Final loglik {final} != last record loglik "
                              f"{good[-1]['loglik']}")
        for prev, cur in zip(good, good[1:]):
            if cur['loglik'] < prev['loglik'] - decrease_tol:
                decreases.append(f"block {name!r} line {cur['line']}: "
                                 f"{prev['loglik']:.6f} -> {cur['loglik']:.6f}")
        block_reports.append({
            'case': name,
            'declared_converged': b['declared'].get('Converged'),
            'declared_iterations': declared_n,
            'records': len(recs),
            'first_iter': iters[0] if iters else None,
            'last_iter': iters[-1] if iters else None,
        })

    names = [b['case'] for b in blocks]
    dup_names = sorted(n for n, c in Counter(names).items() if c > 1)
    if dup_names:
        errors.append(f"duplicate blocks: {dup_names}")
    if expected_cases is not None and names != list(expected_cases):
        errors.append(f"block cases {names} != expected {list(expected_cases)}")

    return {
        'num_blocks': len(blocks),
        'iteration_records_total': sum(len(b['records']) for b in blocks),
        'blocks': block_reports,
        'errors': _cap(errors),
        'valid': not errors,
        'loglik_decreases': {
            'tolerance': decrease_tol,
            'scope': 'consecutive records within each block only; 5 blocks of single fits',
            'count': len(decreases),
            'items': decreases[:MAX_LISTED_ERRORS],
        },
    }


# ============================================================================
# SUMMARIES (bookkeeping counts from the CSVs; not scientific validation)
# ============================================================================

def _conv_counts(rows, key_fn, conv_col='converged'):
    acc = defaultdict(lambda: {'total': 0, 'converged': 0})
    for r in rows:
        k = key_fn(r)
        acc[k]['total'] += 1
        acc[k]['converged'] += int(r[conv_col])
    return {k: dict(v, rate=v['converged'] / v['total']) for k, v in sorted(acc.items())}


def _selection_counts(rows):
    out = {'converged': {'K1': 0, 'K2': 0}, 'not_converged': {'K1': 0, 'K2': 0}}
    for r in rows:
        out['converged' if r['converged'] else 'not_converged'][f"K{r['selected_K']}"] += 1
    return out


def _total(rows, conv_col='converged'):
    n = len(rows)
    c = sum(int(r[conv_col]) for r in rows)
    return {'total_fits': n, 'converged': c, 'rate': c / n if n else None}


def summarize(name, rows):
    if name == 'exact_family':
        configs = _conv_counts(rows, lambda r: f"q={r['q']},kappa={r['kappa']},F={r['F']}")
        return dict(_total(rows),
                    by_F=_conv_counts(rows, lambda r: f"F={r['F']}"),
                    low_conv_configs_below_20pct=[k for k, v in configs.items() if v['rate'] < 0.2],
                    zero_conv_configs=[k for k, v in configs.items() if v['converged'] == 0],
                    selected_K_counts=_selection_counts(rows))
    if name == 'k2_convergence':
        configs = _conv_counts(rows, lambda r: (
            f"sep={r['separation_deg']},occ={r['occupancy_1']}/{r['occupancy_2']},"
            f"kappa={r['kappa']},F={r['F']}"))
        return dict(_total(rows),
                    by_separation_deg=_conv_counts(rows, lambda r: f"sep={r['separation_deg']}"),
                    low_conv_configs_below_50pct=[k for k, v in configs.items() if v['rate'] < 0.5],
                    selected_K_counts=_selection_counts(rows))
    if name == 'stress_regime':
        return dict(_total(rows),
                    by_sigma_F=_conv_counts(rows, lambda r: f"sigma={r['sigma']},F={r['F']}"),
                    selected_K_counts=_selection_counts(rows))
    if name == 'extended_maxiter':
        return dict(_total(rows),
                    by_sigma_F_max_iter=_conv_counts(
                        rows, lambda r: f"sigma={r['sigma']},F={r['F']},max_iter={r['max_iter']}"),
                    design_note='UNPAIRED: seed = 1000 + row index, so draws differ per max_iter')
    if name == 'ninit_sensitivity':
        by_case = defaultdict(dict)
        for r in rows:
            by_case[r['case']][f"n_init={r['n_init']}"] = {
                'best_loglik': r['best_loglik'], 'n_converged': r['n_converged'],
                'best_converged': r['best_converged'], 'spread_range': r['spread_range']}
        w = {r['n_init']: r['best_loglik'] for r in rows if r['case'] == 'wrapped_1.0'}
        return dict(_total(rows, 'best_converged'),
                    status_column='best_converged',
                    by_case=dict(sorted(by_case.items())),
                    wrapped_1_0_best_loglik_diff_n_init_1_vs_3=(
                        abs(w[1] - w[3]) if 1 in w and 3 in w else None),
                    diagnostic_status='INVALID for robustness claims (inner reseeding defect); '
                                      'counts describe recorded output only')
    if name == 'stopping_criterion_validation':
        return {'rows': len(rows),
                'converged': sum(r['converged'] for r in rows),
                'K_changed': sum(r['K_changed'] for r in rows),
                'diagnostic_status': 'INVALID (relative tolerance unimplemented; approximate '
                                     'convergence_iter); counts describe recorded output only'}
    if name == 'initialization_comparison':
        return dict(_total(rows, 'best_converged'), status_column='best_converged',
                    by_case=_conv_counts(rows, lambda r: r['case'], 'best_converged'))
    if name in ('tolerance_comparison', 'max_iter_sensitivity', 'duplicate_components'):
        return _total(rows)
    raise KeyError(name)


# ============================================================================
# REPORT
# ============================================================================

def build_report(results_dir=RESULTS_DIR):
    """Recompute everything from the input bytes. Pure: writes nothing."""
    results_dir = Path(results_dir)
    datasets, summaries = {}, {}
    for name, spec in DATASETS.items():
        result, rows = validate_dataset(name, results_dir / spec['file'], spec)
        datasets[name] = result
        summaries[name] = summarize(name, rows) if result['valid'] else None
    trajectories = parse_trajectory_file(results_dir / TRAJECTORY_NAME)

    failures = [f"{n}: {k}" for n, r in datasets.items()
                for k, bad in [('header', not r['header']['schema_valid']),
                               ('field values', r['field_errors']['count'] > 0),
                               ('row consistency', r['row_consistency_errors']['count'] > 0),
                               ('selection bookkeeping',
                                r['selection_bookkeeping_errors']['count'] > 0),
                               ('coverage', not r['coverage']['complete'])] if bad]
    if not trajectories['valid']:
        failures.append('trajectories: structure')

    return {
        '_meta': {
            'type': 'SYNTHETIC',
            'artifact_role': 'Regenerated derivative summary of SYNTHETIC M3-A2g experiment '
                             'outputs; not an original experiment artifact',
            'generated_by': 'scripts/m3_a2g_data_verification.py --write',
            'baseline_rev': BASELINE_REV,
            'decision_rule': 'selected_K = 2 if delta_bic < 0 else 1 (equality selects K=1)',
            'provenance': (
                'Trajectory record counting (1634 records in 5 blocks) was already corrected '
                'in commit bf33346 and reported at HEAD 394ce0c. This revision adds '
                'generator-derived coverage, strict field parsing, both stopping-diagnostic '
                'selection pairs, strict trajectory validation and read-only verification.'),
            'limitations': [
                'Format, bookkeeping and provenance integrity only; no scientific validation',
                'Consistent selection bookkeeping is not a validated selection; '
                'non-converged fits are counted separately and are not validated',
                'n_init and stopping diagnostics are INVALID; their counts describe recorded '
                'output only',
                'extended_maxiter draws are UNPAIRED across max_iter',
                'Random seeds are not recorded in the CSVs and are not verified',
                'No correspondence between trajectory blocks and CSV rows is asserted',
                f'Likelihood decreases reported with tolerance {LOGLIK_DECREASE_TOL} '
                '(6-decimal resolution) over 5 single-fit trajectories only',
            ],
        },
        'input_hashes': compute_input_hashes(results_dir),
        'csv_row_counts': {n: r['row_count'] for n, r in datasets.items()},
        'total_csv_rows': sum(r['row_count'] for r in datasets.values()),
        'integrity': {'passed': not failures, 'failures': failures},
        'datasets': datasets,
        'summaries': summaries,
        'trajectories': trajectories,
    }


def _normalize(obj):
    return json.loads(json.dumps(obj))


def compare_with_stored(report, json_path):
    """Compare a fresh report with the stored derivative JSON (read-only)."""
    json_path = Path(json_path)
    if not json_path.exists():
        return {'exists': False, 'fresh': False, 'stale_input_hashes': [],
                'differing_sections': ['<file missing>']}
    with open(json_path) as f:
        stored = json.load(f)
    fresh = _normalize(report)
    stored_hashes = stored.get('input_hashes', {}) if isinstance(stored, dict) else {}
    stale = sorted(name for name in set(fresh['input_hashes']) | set(stored_hashes)
                   if stored_hashes.get(name) != fresh['input_hashes'].get(name))
    keys = set(fresh) | (set(stored) if isinstance(stored, dict) else set())
    differing = sorted(k for k in keys
                       if not isinstance(stored, dict) or stored.get(k) != fresh.get(k))
    return {'exists': True, 'fresh': not differing, 'stale_input_hashes': stale,
            'differing_sections': differing}


def write_report(report, json_path):
    with open(json_path, 'w') as f:
        json.dump(report, f, indent=2)
        f.write('\n')


def print_report(report):
    print("Datasets (header / fields / row checks / selection / coverage):")
    for name, r in report['datasets'].items():
        c = r['coverage']
        print(f"  {'OK  ' if r['valid'] else 'FAIL'} {name:30s} rows={r['row_count']:4d} "
              f"fields_err={r['field_errors']['count']} row_err={r['row_consistency_errors']['count']} "
              f"sel_err={r['selection_bookkeeping_errors']['count']} "
              f"keys={c['present_unique']}/{c['expected']} missing={len(c['missing'])} "
              f"unexpected={len(c['unexpected'])} dup={len(c['duplicates'])}")
        for section in ['field_errors', 'row_consistency_errors', 'selection_bookkeeping_errors']:
            for item in r[section]['items'][:5]:
                print(f"       {item}")
        if not r['header']['schema_valid']:
            print(f"       header: {r['header']}")
        for label in ['missing', 'unexpected', 'duplicates']:
            if c[label]:
                print(f"       {label} keys: {c[label][:5]}")
    print(f"  total CSV rows: {report['total_csv_rows']}")
    t = report['trajectories']
    print(f"Trajectories: {'OK' if t['valid'] else 'FAIL'} blocks={t['num_blocks']} "
          f"records={t['iteration_records_total']} "
          f"loglik decreases (tol {t['loglik_decreases']['tolerance']}): "
          f"{t['loglik_decreases']['count']}")
    for item in t['errors']['items'][:10]:
        print(f"       {item}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    parser.add_argument('--results-dir', default=str(RESULTS_DIR))
    parser.add_argument('--json', default=None, help='derivative JSON path')
    parser.add_argument('--write', action='store_true',
                        help='regenerate the derivative JSON (only if integrity passes)')
    args = parser.parse_args(argv)
    results_dir = Path(args.results_dir)
    json_path = Path(args.json) if args.json else results_dir / JSON_NAME

    print("M3-A2g DATA VERIFICATION (SYNTHETIC inputs; read-only unless --write)")
    report = build_report(results_dir)
    print_report(report)

    cmp = compare_with_stored(report, json_path)
    print(f"Stored derivative {json_path.name}: "
          f"{'matches fresh result' if cmp['fresh'] else 'STALE'}")
    if cmp['stale_input_hashes']:
        print(f"  input hashes differing from stored: {cmp['stale_input_hashes']}")
    if not cmp['fresh']:
        print(f"  differing sections: {cmp['differing_sections']}")

    if not report['integrity']['passed']:
        print(f"INTEGRITY FAILED: {report['integrity']['failures']}")
        if args.write:
            print("Refusing to write the derivative JSON.")
        return 1
    if args.write:
        write_report(report, json_path)
        print(f"Wrote {json_path}")
        return 0
    if not cmp['fresh']:
        print("FAIL: stored derivative JSON does not match the inputs "
              "(regenerate explicitly with --write after review)")
        return 1
    print("PASS: integrity checks passed and stored derivative is current "
          "(format/bookkeeping only; not scientific validation)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
