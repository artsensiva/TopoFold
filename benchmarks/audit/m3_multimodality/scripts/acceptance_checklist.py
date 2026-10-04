#!/usr/bin/env python3
"""
M3-A2g acceptance checklist (fresh, non-writing).

Every check is recomputed from the current files; nothing is written.
Exit status is nonzero if any required check fails.

Sections:
  A. Data integrity          format, coverage and selection bookkeeping (fresh)
  B. Documented figures      values cited in M3 documents vs fresh counts
  C. Derivative freshness    stored JSON vs fresh result
  D. Manifest / preservation inventory, checksums, baseline byte comparison

Passing this checklist establishes data-format and provenance integrity only.
It does not validate the optimizer, does not change the provisional O2
classification, and does not resolve KNOWN_ISSUES M2/M3 or production readiness.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from m3_a2g_data_verification import (  # noqa: E402
    BASELINE_REV, JSON_NAME, RESULTS_DIR, build_report, compare_with_stored)
import compute_manifest_hashes as manifest  # noqa: E402

# Reconciliation targets: the row total given in the correction directive and
# figures cited in M3_STATUS.md / README.md. They are compared, never trusted.
DOCUMENTED = {
    'total_csv_rows': 1235,
    'trajectory_blocks': 5,
    'trajectory_records': 1634,
    'exact_family': (111, 360),
    'k2_convergence': (258, 320),
    'stress_regime': (102, 280),
    'stress sigma=1.0,F=1000': (3, 10),
    'extended_maxiter max_iter=2000': (56, 60),
    'ninit wrapped_1.0 loglik diff (2 d.p.)': 2.63,
}


def _items(items, limit=5):
    return str(items[:limit]) if items else ''


class Checklist:
    def __init__(self):
        self.results = []

    def check(self, section, name, ok, detail=''):
        self.results.append((section, name, bool(ok), detail))
        print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  [{detail}]" if detail else ''))
        return ok

    @property
    def failed(self):
        return [r for r in self.results if not r[2]]


def run(results_dir=RESULTS_DIR, json_path=None, provenance=True, rev=BASELINE_REV):
    results_dir = Path(results_dir)
    json_path = Path(json_path) if json_path else results_dir / JSON_NAME
    cl = Checklist()
    report = build_report(results_dir)

    print("\nA. Data integrity (fresh)")
    for name, r in report['datasets'].items():
        c = r['coverage']
        cl.check('A', f"{name}: header", r['header']['schema_valid'],
                 '' if r['header']['schema_valid'] else str(r['header']))
        cl.check('A', f"{name}: field values", r['field_errors']['count'] == 0,
                 '; '.join(r['field_errors']['items'][:3]))
        cl.check('A', f"{name}: row consistency", r['row_consistency_errors']['count'] == 0,
                 '; '.join(r['row_consistency_errors']['items'][:3]))
        cl.check('A', f"{name}: coverage {c['present_unique']}/{c['expected']} keys",
                 c['complete'],
                 f"missing={c['missing'][:3]} unexpected={c['unexpected'][:3]} "
                 f"duplicates={c['duplicates'][:3]}" if not c['complete'] else '')
        if r['selection_pairs_checked']:
            cl.check('A', f"{name}: selection bookkeeping "
                          f"{'+'.join(k for k, _ in r['selection_pairs_checked'])}",
                     r['selection_bookkeeping_errors']['count'] == 0,
                     '; '.join(r['selection_bookkeeping_errors']['items'][:3]))
    t = report['trajectories']
    cl.check('A', 'trajectories: declarations, counts, sequences, records', t['valid'],
             '; '.join(t['errors']['items'][:3]))

    print("\nB. Documented figures vs fresh data")
    s = report['summaries']
    cl.check('B', f"total CSV rows = {DOCUMENTED['total_csv_rows']}",
             report['total_csv_rows'] == DOCUMENTED['total_csv_rows'],
             f"fresh {report['total_csv_rows']}")
    cl.check('B', f"trajectory blocks/records = {DOCUMENTED['trajectory_blocks']}/"
                  f"{DOCUMENTED['trajectory_records']}",
             (t['num_blocks'], t['iteration_records_total'])
             == (DOCUMENTED['trajectory_blocks'], DOCUMENTED['trajectory_records']),
             f"fresh {t['num_blocks']}/{t['iteration_records_total']}")

    def frac(name, got):
        want = DOCUMENTED[name]
        cl.check('B', f"{name} = {want[0]}/{want[1]}", got == want,
                 f"fresh {got[0]}/{got[1]}" if got else 'unavailable (dataset invalid)')

    def total(name):
        return (s[name]['converged'], s[name]['total_fits']) if s.get(name) else None

    for name in ['exact_family', 'k2_convergence', 'stress_regime']:
        frac(name, total(name))
    ex = s['stress_regime']['by_sigma_F'].get('sigma=1.0,F=1000') if s.get('stress_regime') else None
    frac('stress sigma=1.0,F=1000', (ex['converged'], ex['total']) if ex else None)
    em = None
    if s.get('extended_maxiter'):
        cells = [v for k, v in s['extended_maxiter']['by_sigma_F_max_iter'].items()
                 if k.endswith('max_iter=2000')]
        em = (sum(v['converged'] for v in cells), sum(v['total'] for v in cells))
    frac('extended_maxiter max_iter=2000', em)
    diff = s['ninit_sensitivity']['wrapped_1_0_best_loglik_diff_n_init_1_vs_3'] \
        if s.get('ninit_sensitivity') else None
    want = DOCUMENTED['ninit wrapped_1.0 loglik diff (2 d.p.)']
    cl.check('B', f"ninit wrapped_1.0 loglik diff = {want}",
             diff is not None and round(diff, 2) == want, f"fresh {diff}")

    print("\nC. Derivative freshness (read-only comparison)")
    cmp = compare_with_stored(report, json_path)
    cl.check('C', f"stored {json_path.name} exists", cmp['exists'])
    cl.check('C', 'stored input hashes match current input bytes',
             cmp['exists'] and not cmp['stale_input_hashes'], _items(cmp['stale_input_hashes']))
    cl.check('C', 'stored derivative equals fresh result', cmp['fresh'],
             f"differing sections: {cmp['differing_sections']}" if not cmp['fresh'] else '')

    if provenance:
        print("\nD. Manifest and preservation")
        m = manifest.check_manifest()
        for key in ['duplicate_entries', 'omitted_from_manifest', 'listed_but_missing',
                    'checksum_mismatches', 'unclassified_files', 'classified_but_missing']:
            cl.check('D', f"manifest: no {key.replace('_', ' ')}", not m[key], _items(m[key]))
        missing_tracked = sorted(p for p in manifest.tracked_files()
                                 if not (manifest.BASE_DIR / p).is_file())
        cl.check('D', 'all tracked M3 files present', not missing_tracked, _items(missing_tracked))
        p = manifest.compare_protected(rev=rev)
        cl.check('D', f"protected artifacts identical to {rev[:12]} "
                      f"({p['identical']}/{p['protected']})",
                 p['passed'], '' if p['passed'] else f"changed={p['changed']} missing={p['missing']}")

    return cl, report


def main(argv=None):
    parser = argparse.ArgumentParser(description='M3-A2g acceptance checklist')
    parser.add_argument('--results-dir', default=str(RESULTS_DIR))
    parser.add_argument('--json', default=None)
    parser.add_argument('--skip-provenance', action='store_true',
                        help='skip section D (manifest and git baseline comparison)')
    parser.add_argument('--rev', default=BASELINE_REV)
    args = parser.parse_args(argv)

    print("M3-A2g ACCEPTANCE CHECKLIST (fresh, non-writing; SYNTHETIC inputs)")
    cl, _ = run(args.results_dir, args.json, not args.skip_provenance, args.rev)

    print("\nScientific status (information only; not changed by this checklist)")
    print("  Provisional O2 - functional but fragile. n_init and stopping diagnostics INVALID;")
    print("  extended max_iter UNPAIRED; monotonicity evidence limited; dependent-toroidal")
    print("  validation not performed. Not production-ready. KNOWN_ISSUES M3 remains open.")

    print(f"\nSummary: {len(cl.results) - len(cl.failed)}/{len(cl.results)} checks passed")
    if cl.failed:
        for section, name, _, detail in cl.failed:
            print(f"  FAIL [{section}] {name}" + (f"  [{detail}]" if detail else ''))
        print("RESULT: FAIL")
        return 1
    print("RESULT: PASS (format, bookkeeping and provenance integrity only)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
