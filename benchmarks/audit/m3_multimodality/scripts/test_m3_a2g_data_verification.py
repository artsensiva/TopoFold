#!/usr/bin/env python3
"""
Regression tests for the M3-A2g verifier, acceptance checklist and manifest checker.

Every corruption is applied to a copy of the inputs in a pytest tmp_path;
repository artifacts are only read, never modified.

Run: python -B -m pytest -p no:cacheprovider scripts/test_m3_a2g_data_verification.py
"""
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))
import compute_manifest_hashes as manifest  # noqa: E402
from m3_a2g_data_verification import (  # noqa: E402
    DATASETS, JSON_NAME, RESULTS_DIR, TRAJECTORY_NAME, build_report, compare_with_stored,
    expected_selected_k, parse_trajectory_file, validate_dataset, write_report)

INPUTS = [spec['file'] for spec in DATASETS.values()] + [TRAJECTORY_NAME]
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')


@pytest.fixture
def results(tmp_path):
    """Isolated copy of the 11 original inputs (no derivative JSON)."""
    d = tmp_path / 'results'
    d.mkdir()
    for name in INPUTS:
        shutil.copy2(RESULTS_DIR / name, d / name)
    return d


def edit_line(path, line_index, fn):
    lines = path.read_text().split('\n')
    lines[line_index] = fn(lines[line_index])
    path.write_text('\n'.join(lines))


def set_field(path, line_index, column, value):
    header = path.read_text().split('\n')[0].split(',')
    col = header.index(column)

    def fn(line):
        parts = line.split(',')
        parts[col] = value
        return ','.join(parts)
    edit_line(path, line_index, fn)


def run_cli(script, *args):
    return subprocess.run([sys.executable, '-B', str(SCRIPTS / script), *map(str, args)],
                          capture_output=True, text=True, env=ENV)


def dataset(results, name):
    return validate_dataset(name, results / DATASETS[name]['file'])[0]


# ---------------------------------------------------------------------------
# Real inputs (read-only)
# ---------------------------------------------------------------------------

def test_real_inputs_pass_integrity():
    report = build_report(RESULTS_DIR)
    assert report['integrity'] == {'passed': True, 'failures': []}
    assert report['total_csv_rows'] == 1235
    assert report['trajectories']['iteration_records_total'] == 1634
    assert [b['records'] for b in report['trajectories']['blocks']] == [500, 500, 500, 9, 125]


def test_ninit_summary_uses_best_converged():
    s = build_report(RESULTS_DIR)['summaries']['ninit_sensitivity']
    assert s['status_column'] == 'best_converged'
    assert (s['converged'], s['total_fits']) == (15, 20)
    assert 'INVALID' in s['diagnostic_status']


def test_build_report_writes_nothing(results):
    before = sorted(p.name for p in results.iterdir())
    build_report(results)
    assert sorted(p.name for p in results.iterdir()) == before


# ---------------------------------------------------------------------------
# Decision rule
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('delta_bic, k', [(-1e-12, 2), (-5.0, 2), (0.0, 1), (-0.0, 1), (3.0, 1)])
def test_expected_selected_k_equality_selects_k1(delta_bic, k):
    assert expected_selected_k(delta_bic) == k


def test_selected_k2_at_zero_delta_bic_fails(results):
    f = results / DATASETS['exact_family']['file']
    set_field(f, 1, 'delta_bic', '0.0')
    set_field(f, 1, 'selected_K', '2')
    r = dataset(results, 'exact_family')
    assert r['selection_bookkeeping_errors']['count'] == 1
    assert not r['valid']


# ---------------------------------------------------------------------------
# Independently demonstrated corruptions (review findings)
# ---------------------------------------------------------------------------

def test_out_of_range_replicate_99_fails(results):
    f = results / DATASETS['exact_family']['file']
    set_field(f, 1, 'rep', '99')  # replaces (q=1, kappa=1, F=100, rep=0)
    r = dataset(results, 'exact_family')
    assert r['coverage']['unexpected'] == [[1, 1, 100, 99]]
    assert r['coverage']['missing'] == [[1, 1, 100, 0]]
    assert not r['valid']
    assert run_cli('m3_a2g_data_verification.py', '--results-dir', results).returncode == 1
    assert run_cli('acceptance_checklist.py', '--results-dir', results,
                   '--skip-provenance').returncode == 1


def test_inconsistent_continuation_selection_fails(results):
    f = results / DATASETS['stopping_criterion_validation']['file']
    # Row 1: delta_bic_cont > 0 and selected_K_cont = 1. Flip the selection and keep
    # K_changed consistent, so only the continuation selection pair is wrong.
    set_field(f, 1, 'selected_K_cont', '2')
    set_field(f, 1, 'K_changed', 'True')
    r = dataset(results, 'stopping_criterion_validation')
    errs = r['selection_bookkeeping_errors']
    assert errs['count'] == 1 and 'selected_K_cont' in errs['items'][0]
    assert r['row_consistency_errors']['count'] == 0
    assert not r['valid']
    assert run_cli('acceptance_checklist.py', '--results-dir', results,
                   '--skip-provenance').returncode == 1


def test_conv_selection_pair_also_checked(results):
    f = results / DATASETS['stopping_criterion_validation']['file']
    set_field(f, 1, 'selected_K_conv', '2')
    set_field(f, 1, 'K_changed', 'True')
    errs = dataset(results, 'stopping_criterion_validation')['selection_bookkeeping_errors']
    assert errs['count'] == 1 and 'selected_K_conv' in errs['items'][0]


def test_broken_trajectory_iteration_numbering_fails(results):
    f = results / TRAJECTORY_NAME
    lines = f.read_text().split('\n')
    first = next(i for i, line in enumerate(lines) if line.startswith('0,'))
    edit_line(f, first, lambda line: '99999' + line[1:])
    t = parse_trajectory_file(f)
    assert not t['valid']
    assert any('iteration sequence broken' in e for e in t['errors']['items'])
    assert run_cli('m3_a2g_data_verification.py', '--results-dir', results).returncode == 1
    assert run_cli('acceptance_checklist.py', '--results-dir', results,
                   '--skip-provenance').returncode == 1


# ---------------------------------------------------------------------------
# Coverage, headers, field values
# ---------------------------------------------------------------------------

def test_missing_replicate_fails(results):
    f = results / DATASETS['stress_regime']['file']
    lines = f.read_text().split('\n')
    del lines[5]
    f.write_text('\n'.join(lines))
    r = dataset(results, 'stress_regime')
    assert len(r['coverage']['missing']) == 1 and not r['valid']


def test_duplicate_replicate_fails(results):
    f = results / DATASETS['k2_convergence']['file']
    lines = f.read_text().split('\n')
    lines.insert(3, lines[2])
    f.write_text('\n'.join(lines))
    r = dataset(results, 'k2_convergence')
    assert len(r['coverage']['duplicates']) == 1 and r['row_count'] == 321
    assert not r['valid']


def test_duplicate_header_fails(results):
    f = results / DATASETS['exact_family']['file']
    edit_line(f, 0, lambda line: line.replace('loglik', 'delta_bic'))
    r = dataset(results, 'exact_family')
    assert r['header']['duplicate_columns'] == ['delta_bic']
    assert r['header']['missing_columns'] == ['loglik']
    assert not r['header']['schema_valid'] and not r['valid']


@pytest.mark.parametrize('column, value', [
    ('delta_bic', ''),        # empty required value
    ('F', 'abc'),             # unparseable int
    ('loglik', 'inf'),        # non-finite
    ('delta_bic', 'nan'),     # non-finite, no sentinel documented
    ('converged', 'yes'),     # invalid bool
    ('selected_K', '3'),      # outside domain
    ('q', '1.0'),             # int column with float text
])
def test_malformed_required_field_fails(results, column, value):
    f = results / DATASETS['exact_family']['file']
    set_field(f, 1, column, value)
    r = dataset(results, 'exact_family')
    assert r['field_errors']['count'] >= 1
    assert not r['records_valid'] and not r['valid']


def test_short_row_fails(results):
    f = results / DATASETS['tolerance_comparison']['file']
    edit_line(f, 1, lambda line: line.rsplit(',', 1)[0])
    assert dataset(results, 'tolerance_comparison')['field_errors']['count'] == 1


def test_nan_sentinel_must_match_generator_condition(results):
    f = results / DATASETS['max_iter_sensitivity']['file']
    set_field(f, 1, 'improvement_after_100', '0.1')   # n_iters=100: generator writes nan
    set_field(f, 2, 'improvement_after_100', 'nan')   # n_iters=250: must be finite
    r = dataset(results, 'max_iter_sensitivity')
    assert r['field_errors']['count'] == 2


def test_kappa_zero_is_permitted(results):
    f = results / DATASETS['stress_regime']['file']
    set_field(f, 1, 'max_kappa', '0.0')
    assert dataset(results, 'stress_regime')['valid']


def test_row_consistency_ninit_rate(results):
    f = results / DATASETS['ninit_sensitivity']['file']
    set_field(f, 1, 'convergence_rate', '0.5')  # n_converged=0, n_init=1
    r = dataset(results, 'ninit_sensitivity')
    assert r['row_consistency_errors']['count'] == 1


# ---------------------------------------------------------------------------
# Trajectory structure
# ---------------------------------------------------------------------------

def _traj(results):
    return results / TRAJECTORY_NAME


def test_trajectory_declared_count_mismatch_fails(results):
    f = _traj(results)
    lines = f.read_text().split('\n')
    last = max(i for i, line in enumerate(lines) if line[:1].isdigit())
    del lines[last]
    f.write_text('\n'.join(lines))
    t = parse_trajectory_file(f)
    assert any('declared Iterations=125, found 124' in e for e in t['errors']['items'])


def test_trajectory_malformed_record_not_skipped(results):
    f = _traj(results)
    lines = f.read_text().split('\n')
    idx = next(i for i, line in enumerate(lines) if line.startswith('5,'))
    edit_line(f, idx, lambda line: line.rsplit(',', 1)[0])
    t = parse_trajectory_file(f)
    assert not t['valid']
    assert t['iteration_records_total'] == 1634  # counted, flagged, not dropped
    assert any('malformed record' in e for e in t['errors']['items'])


def test_trajectory_missing_declaration_fails(results):
    f = _traj(results)
    f.write_text(f.read_text().replace('Iterations: 9\n', '', 1))
    t = parse_trajectory_file(f)
    assert any("missing declaration 'Iterations'" in e for e in t['errors']['items'])


def test_trajectory_likelihood_decrease_is_reported_not_failed(tmp_path):
    f = tmp_path / 't.txt'
    f.write_text('Case: exact_easy\nConverged: False\nIterations: 3\nFinal loglik: -100.5\n\n'
                 'iter,loglik,delta_abs,delta_rel,max_kappa,min_weight,eff_n_comp\n'
                 '0,-100.0,nan,nan,5.0,0.5,2.0\n1,-101.0,-1.0e+00,-1.0e-02,5.0,0.5,2.0\n'
                 '2,-100.5,5.0e-01,5.0e-03,5.0,0.5,2.0\n\n')
    t = parse_trajectory_file(f, expected_cases=['exact_easy'])
    assert t['valid']
    assert t['loglik_decreases']['count'] == 1


# ---------------------------------------------------------------------------
# Derivative freshness and command behaviour
# ---------------------------------------------------------------------------

def test_stale_checksum_detected(results):
    json_path = results / JSON_NAME
    write_report(build_report(results), json_path)
    assert compare_with_stored(build_report(results), json_path)['fresh']
    f = results / DATASETS['duplicate_components']['file']
    set_field(f, 1, 'final_loglik', '-3246.9')  # bytes change, record still valid
    assert build_report(results)['integrity']['passed']
    cmp = compare_with_stored(build_report(results), json_path)
    assert cmp['stale_input_hashes'] == [DATASETS['duplicate_components']['file']]
    assert not cmp['fresh']
    assert run_cli('m3_a2g_data_verification.py', '--results-dir', results).returncode == 1


def test_verifier_default_does_not_write(results):
    proc = run_cli('m3_a2g_data_verification.py', '--results-dir', results)
    assert proc.returncode == 1  # no stored JSON -> stale
    assert not (results / JSON_NAME).exists()


def test_write_refused_on_integrity_failure(results):
    set_field(results / DATASETS['exact_family']['file'], 1, 'rep', '99')
    proc = run_cli('m3_a2g_data_verification.py', '--results-dir', results, '--write')
    assert proc.returncode == 1
    assert not (results / JSON_NAME).exists()


def test_write_then_verify_and_accept(results):
    assert run_cli('m3_a2g_data_verification.py', '--results-dir', results,
                   '--write').returncode == 0
    stored = json.loads((results / JSON_NAME).read_text())
    assert stored['_meta']['type'] == 'SYNTHETIC'
    assert run_cli('m3_a2g_data_verification.py', '--results-dir', results).returncode == 0
    proc = run_cli('acceptance_checklist.py', '--results-dir', results, '--skip-provenance')
    assert proc.returncode == 0, proc.stdout[-2000:]


# ---------------------------------------------------------------------------
# Manifest checker (isolated directory)
# ---------------------------------------------------------------------------

def _manifest_dir(tmp_path, rows):
    d = tmp_path / 'm3'
    d.mkdir()
    (d / 'a.csv').write_text('x\n1\n')
    (d / 'b.md').write_text('doc\n')
    body = '\n'.join(f"| `{p}` | note | `{h}` |" for p, h in rows)
    (d / 'MANIFEST.md').write_text(f"# M\n\n| Path | Note | SHA-256 |\n|---|---|---|\n{body}\n")
    return d


def test_manifest_detects_stale_omitted_duplicate(tmp_path):
    good = manifest.hashlib.sha256(b'x\n1\n').hexdigest()
    d = _manifest_dir(tmp_path, [('a.csv', good), ('a.csv', good), ('gone.md', good)])
    classes = {'a.csv': ('c', ''), 'b.md': ('c', '')}
    m = manifest.check_manifest(d, classification=classes)
    assert m['duplicate_entries'] and m['omitted_from_manifest'] == ['b.md']
    assert m['listed_but_missing'] == ['gone.md'] and not m['passed']
    (d / 'a.csv').write_text('x\n2\n')
    m = manifest.check_manifest(d, classification=classes)
    assert len(m['checksum_mismatches']) == 2  # both duplicate rows are stale


def test_manifest_detects_unclassified_file(tmp_path):
    d = _manifest_dir(tmp_path, [])
    m = manifest.check_manifest(d, classification={'a.csv': ('c', '')})
    assert m['unclassified_files'] == ['b.md']
