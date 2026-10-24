#!/usr/bin/env python3
"""AI tester 1 independent, standard-library-only trial checks.

Prepared before trial output inspection. No TopoFold modules are imported.
Numerical tolerance is abs(got-ref) <= 1e-8 + 1e-7*abs(ref).
CSV missing/non-numeric values and JSON non-numeric leaves compare exactly.
Machine-specific provenance, software versions, and timestamps may differ.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys

ATOL = 1e-8
RTOL = 1e-7
NUMERIC = re.compile(r'^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$')
GROUPS = {
    3: ['kappa_A_rad', 'kappa_B_rad', 'd2_A_angstrom', 'd2_B_angstrom', 'd_kappa_rad', 'd_i_i2_angstrom', 'rmsd3_angstrom'],
    4: ['tau_A_rad', 'tau_B_rad', 'd3_A_angstrom', 'd3_B_angstrom', 'd_tau_rad', 'd_i_i3_angstrom', 'rmsd4_angstrom'],
}
SHARED = {3: ['d_kappa_rad', 'd_i_i2_angstrom', 'rmsd3_angstrom'], 4: ['d_tau_rad', 'd_i_i3_angstrom', 'rmsd4_angstrom']}
OWN = {3: ['kappa_{side}_rad', 'd2_{side}_angstrom'], 4: ['tau_{side}_rad', 'd3_{side}_angstrom']}

class Checks:
    def __init__(self):
        self.counts = Counter()
        self.failures = []
        self.max_abs_difference = 0.0
        self.max_tolerance_ratio = 0.0
    def check(self, ok, label, detail=None):
        self.counts['checks'] += 1
        self.counts['passed' if ok else 'failed'] += 1
        if not ok:
            self.failures.append({'check': label, 'detail': detail})
    def compare(self, got, ref, label):
        if isinstance(ref, dict) and isinstance(got, dict):
            self.check(set(got) == set(ref), label + ': keys', {'got': sorted(got), 'reference': sorted(ref)})
            for key in sorted(set(got) & set(ref)):
                self.compare(got[key], ref[key], label + '/' + str(key))
        elif isinstance(ref, list) and isinstance(got, list):
            self.check(len(got) == len(ref), label + ': length', {'got': len(got), 'reference': len(ref)})
            for i, (a, b) in enumerate(zip(got, ref)):
                self.compare(a, b, label + '/' + str(i))
        else:
            numeric = (isinstance(ref, (float, int)) and not isinstance(ref, bool) and isinstance(got, (float, int)) and not isinstance(got, bool))
            numeric = numeric or (isinstance(ref, str) and isinstance(got, str) and NUMERIC.fullmatch(ref) and NUMERIC.fullmatch(got))
            if numeric:
                a, b = float(got), float(ref)
                diff, tol = abs(a-b), ATOL + RTOL*abs(b)
                self.counts['numeric_values'] += 1
                if math.isfinite(diff):
                    self.max_abs_difference = max(self.max_abs_difference, diff)
                    self.max_tolerance_ratio = max(self.max_tolerance_ratio, diff/tol)
                self.check(math.isfinite(a) and math.isfinite(b) and diff <= tol, label, {'got': got, 'reference': ref, 'difference': diff, 'tolerance': tol})
            else:
                self.counts['exact_values'] += 1
                self.check(type(got) is type(ref) and got == ref, label, {'got': got, 'reference': ref})
    def report(self):
        return {'counts': dict(self.counts), 'max_absolute_difference': self.max_abs_difference, 'max_tolerance_ratio': self.max_tolerance_ratio, 'failures': self.failures}

def read_csv(path):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return reader.fieldnames, rows

def load_json(path):
    return json.loads(path.read_text())

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def ca_records(path, chain):
    """Parse deposited first-model fixed-column C-alpha records independently."""
    result = {}
    for line in path.read_text().splitlines():
        record = line[:6].strip()
        if record == 'ENDMDL':
            break
        if record not in ('ATOM', 'HETATM') or line[12:16].strip() != 'CA' or line[21:22] != chain:
            continue
        if line[17:20].strip() == 'CA':  # a calcium ion is not a protein C-alpha
            continue
        key = str(int(line[22:26])) + line[26:27].strip()
        result.setdefault(key, []).append({'resname': line[17:20].strip(), 'altloc': line[16:17].strip(), 'xyz': [float(line[30:38]), float(line[38:46]), float(line[46:54])]})
    return result

def compare_reference(trial, reference):
    checks = Checks()
    table_sizes = {}
    for filename in ('windows.csv', 'residues.csv', 'exclusions.csv'):
        ghead, got = read_csv(trial/filename)
        rhead, ref = read_csv(reference/filename)
        checks.compare(ghead, rhead, filename + '/header')
        checks.compare(got, ref, filename + '/rows')
        table_sizes[filename] = {'trial_rows': len(got), 'reference_rows': len(ref), 'columns': len(ghead)}
    checks.compare(load_json(trial/'summary.json'), load_json(reference/'summary.json'), 'summary.json')
    gm, rm = load_json(trial/'manifest.json'), load_json(reference/'manifest.json')
    for key in ('schema', 'name', 'statement', 'correspondence', 'policies', 'breaks'):
        checks.compare(gm[key], rm[key], 'manifest/' + key)
    checks.compare(gm['spec']['sha256'], rm['spec']['sha256'], 'manifest/spec/sha256')
    for side in ('A', 'B'):
        for key in ('label', 'sha256', 'chain', 'altloc_policy', 'units', 'frame_identity'):
            checks.compare(gm['inputs'][side][key], rm['inputs'][side][key], 'manifest/inputs/' + side + '/' + key)
    return dict(checks.report(), table_sizes=table_sizes, byte_identical_files={name: sha(trial/name) == sha(reference/name) for name in ('windows.csv', 'residues.csv', 'exclusions.csv', 'summary.json')}, scientific_summary=load_json(trial/'summary.json'))

def verify_manifest(directory, invocation_cwd):
    checks = Checks()
    manifest = load_json(directory/'manifest.json')
    outputs = manifest['outputs_sha256']
    for name in ('summary.json', 'windows.csv', 'residues.csv', 'exclusions.csv'):
        checks.check(name in outputs, 'manifest lists required output ' + name)
    for name, expected in outputs.items():
        path = directory/name
        checks.check(path.is_file(), 'output exists ' + name)
        if path.is_file():
            checks.compare(sha(path), expected, 'output SHA256 ' + name)
    for side, entry in manifest['inputs'].items():
        path = Path(entry['path'])
        checks.check(path.is_file(), 'input exists ' + side)
        if path.is_file():
            checks.compare(sha(path), entry['sha256'], 'input SHA256 ' + side)
    spec = Path(manifest['spec']['path'])
    if not spec.is_absolute():
        spec = invocation_cwd/spec
    checks.check(spec.is_file(), 'spec exists relative to recorded invocation cwd')
    if spec.is_file():
        checks.compare(sha(spec), manifest['spec']['sha256'], 'spec SHA256')
    return checks.report()

def verify_extra(directory, primary, spec_path):
    checks = Checks()
    spec = load_json(spec_path)
    pairs = [[str(a), str(b)] for a, b in spec['correspondence']['pairs']]
    checks.compare(pairs, [[str(i), str(i)] for i in range(55, 73)], 'preselected explicit correspondence')
    present = {}
    for side in ('A', 'B'):
        info = spec['structures'][side]
        path = (spec_path.parent/info['file']).resolve()
        present[side] = ca_records(path, info['chain'])
    _, residues = read_csv(directory/'residues.csv')
    _, windows = read_csv(directory/'windows.csv')
    _, exclusions = read_csv(directory/'exclusions.csv')
    _, primary_windows = read_csv(primary/'windows.csv')
    summary = load_json(directory/'summary.json')
    manifest = load_json(directory/'manifest.json')
    checks.compare(len(residues), 18, 'extra mapping row count')
    checks.compare(len(windows), 18, 'extra window position row count')
    checks.compare(manifest['correspondence']['pairs'], pairs, 'extra manifest correspondence')
    missing = {side: [] for side in ('A', 'B')}
    for i, row in enumerate(residues):
        checks.compare(row['position'], str(i), 'residue position ' + str(i))
        for j, side in enumerate(('A', 'B')):
            resid = pairs[i][j]
            label = spec['structures'][side]['chain'] + ':' + resid
            checks.compare(row[side + '_residue'], label, 'residue label ' + side + '/' + resid)
            exists = resid in present[side]
            if exists:
                records = present[side][resid]
                checks.compare(row[side + '_status'], 'present', 'present status ' + side + '/' + resid)
                checks.compare(row[side + '_resname'], records[0]['resname'], 'present residue name ' + side + '/' + resid)
                checks.check(len(records) == 1 and records[0]['altloc'] == '', 'unambiguous C-alpha ' + side + '/' + resid, records)
            else:
                missing[side].append(resid)
                checks.compare(row[side + '_status'], 'residue_not_found', 'missing status ' + side + '/' + resid)
                checks.compare(row[side + '_resname'], '', 'missing residue name ' + side + '/' + resid)
    expected_excluded = set()
    expected_missing_occurrences = set()
    nonoverlapping_valid_windows = []
    expected_valid = Counter()
    matched_metrics = 0
    primary_by_window = {size: {row['window' + str(size) + '_A']: row for row in primary_windows if row['window' + str(size) + '_A']} for size in (3, 4)}
    for i, row in enumerate(windows):
        checks.compare(row['position'], str(i), 'window position ' + str(i))
        for size in (3, 4):
            indices = list(range(i-1, i-1+size))
            complete = min(indices) >= 0 and max(indices) < len(pairs)
            if not complete:
                for column in GROUPS[size]:
                    checks.compare(row[column], '', 'boundary missing value ' + str(i) + '/' + column)
                continue
            label = 'A:' + pairs[indices[0]][0] + '..A:' + pairs[indices[-1]][0]
            checks.compare(row['window' + str(size) + '_A'], label, 'window label ' + str(i) + '/' + str(size))
            missing_sides = []
            for j, side in enumerate(('A', 'B')):
                absent = [pairs[k][j] for k in indices if pairs[k][j] not in present[side]]
                if absent:
                    missing_sides.append(side)
                    expected_excluded.add((size, i, side))
                    expected_missing_occurrences.update((size, i, side, resid) for resid in absent)
                    for template in OWN[size]:
                        column = template.format(side=side)
                        checks.compare(row[column], '', 'absent own geometry ' + str(i) + '/' + column)
            if missing_sides:
                for column in SHARED[size]:
                    checks.compare(row[column], '', 'excluded comparison value ' + str(i) + '/' + column)
            else:
                expected_valid[size] += 1
                prior = primary_by_window[size].get(label)
                if prior is None:
                    nonoverlapping_valid_windows.append({'size': size, 'window_A': label})
                for column in GROUPS[size]:
                    checks.check(bool(NUMERIC.fullmatch(row[column])), 'valid numeric window ' + label + '/' + column, row[column])
                    if prior is not None:
                        checks.compare(row[column], prior[column], 'matching example ' + label + '/' + column)
                        matched_metrics += 1
    observed_excluded = set()
    observed_missing_occurrences = []
    cause_counts = Counter()
    for row in exclusions:
        key = (int(row['window_size']), int(row['position']), row['structure'])
        observed_excluded.add(key)
        cause_counts[row['cause']] += 1
        checks.check(key in expected_excluded, 'exclusion predicted by deposited missing C-alpha', row)
        checks.check('missing' in row['cause'].lower(), 'exclusion names missing cause', row)
        checks.check(bool(row['detail']) and bool(row['location']), 'exclusion has detail and location', row)
        size, position, side = key
        side_index = ('A', 'B').index(side)
        absent = {pairs[k][side_index] for k in range(position-1, position-1+size) if pairs[k][side_index] not in present[side]}
        locations = {value.split(':')[-1] for value in row['location'].split('→')}
        named_missing = absent & locations
        checks.check(len(named_missing) == 1, 'cause location identifies a deposited missing C-alpha in this window', row)
        if len(named_missing) == 1:
            observed_missing_occurrences.append((size, position, side, next(iter(named_missing))))
    checks.check(observed_excluded == expected_excluded, 'all and only missing windows excluded', {'missing': sorted(expected_excluded-observed_excluded), 'unexpected': sorted(observed_excluded-expected_excluded)})
    checks.compare(sorted(observed_missing_occurrences), sorted(expected_missing_occurrences), 'one explanation per missing C-alpha occurrence in eligible window')
    checks.compare(len(exclusions), len(expected_missing_occurrences), 'independently expected cause-row count')
    checks.compare(summary['n_exclusion_rows'], len(exclusions), 'summary exclusion count')
    for metric, size in [('d_kappa', 3), ('d_i_i2', 3), ('rmsd3', 3), ('d_tau', 4), ('d_i_i3', 4), ('rmsd4', 4)]:
        checks.compare(summary['n_windows'][metric], expected_valid[size], 'summary valid ' + metric)
    common = sum(pair[0] in present['A'] and pair[1] in present['B'] for pair in pairs)
    checks.compare(summary['global']['n_ca'], common, 'global common C-alpha count')
    return dict(checks.report(), deposited_missing_CA=missing, expected_valid_windows=dict(expected_valid), expected_excluded_window_structure_pairs=len(expected_excluded), expected_missing_CA_occurrences=len(expected_missing_occurrences), nonoverlapping_valid_windows=nonoverlapping_valid_windows, actual_exclusion_rows=len(exclusions), exclusion_causes=dict(cause_counts), overlapping_numeric_metrics_checked=matched_metrics, scientific_summary=summary)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trial', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = args.trial.resolve()
    output = args.output or root/'evidence/validation_results.json'
    report = {'tester': 'AI tester 1', 'validation_script_sha256': sha(Path(__file__)), 'started_utc': datetime.now(timezone.utc).isoformat(), 'tolerance': {'absolute': ATOL, 'relative': RTOL, 'formula': 'abs(got-ref) <= 1e-8 + 1e-7*abs(ref)'}, 'scope': 'Reference numerical tables and summary; selected scientific manifest fields; current file hashes; independently parsed deposited first-model PDB C-alpha presence and preselected explicit mapping. No independent geometry implementation or biological validation.'}
    exit_code = 0
    try:
        report['reference_comparison'] = compare_reference(root/'results/kras2', root/'checkout/benchmarks/audit/structpair_demo/results/kras2')
        report['primary_manifest_integrity'] = verify_manifest(root/'results/kras2', root/'checkout')
        report['additional_manifest_integrity'] = verify_manifest(root/'results/explicit_55_72', root/'checkout')
        report['additional_case'] = verify_extra(root/'results/explicit_55_72', root/'results/kras2', root/'specs/explicit_55_72.json')
        exit_code = int(any(value.get('failures') for value in report.values() if isinstance(value, dict)))
    except Exception as exc:
        report['error'] = {'type': type(exc).__name__, 'message': str(exc)}
        exit_code = 2
    report['completed_utc'] = datetime.now(timezone.utc).isoformat()
    report['exit_code'] = exit_code
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'report': str(output), 'exit_code': exit_code, 'sections': {key: {'counts': val.get('counts'), 'failure_count': len(val.get('failures', []))} for key, val in report.items() if isinstance(val, dict) and 'counts' in val}, 'error': report.get('error')}, indent=2))
    return exit_code

if __name__ == '__main__':
    sys.exit(main())
