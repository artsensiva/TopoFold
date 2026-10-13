#!/usr/bin/env python3
"""EB-1 Amendment 1 report (REPORT_A1.md) from results/eb1a1_measurements.json and results/eb1a1_analysis.json.
Every number is read from those files. The original REPORT.md is not modified.
Run: ../wpd_local_retrieval/build/venv/bin/python -B make_report_a1.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
M = json.loads((HERE / 'results' / 'eb1a1_measurements.json').read_text())
A = json.loads((HERE / 'results' / 'eb1a1_analysis.json').read_text())
F = json.loads((HERE / 'results' / 'eb1a1_freeze.json').read_text())
OUT = HERE / 'REPORT_A1.md'


def f(x, d=3):
    if x is None:
        return '—'
    return f'{x:.{d}g}' if isinstance(x, float) else str(x)


def main():
    L = []
    w = L.append
    meta = M['_meta']
    w('# EB-1 Amendment 1: corrected accounting and targeted confirmation')
    w('')
    w(f"Amendment `docs/EB1_AMENDMENT_1.md` (SHA-256 `{F['sha256']['docs/EB1_AMENDMENT_1.md'][:16]}…`, frozen "
      f"{F['frozen_utc']}, after the original results and before any amended measurement). Original EB-1 protocol, "
      'measurements, analysis, report and outcome are unchanged. Development data only.')
    w('')
    w(f"Cells: {meta.get('n_cells_ok')} of {meta.get('n_cells_planned')} completed. CPU: amendment "
      f"{f(meta.get('amendment_cpu_h'))} h; aggregate including the original run {f(meta.get('aggregate_cpu_h_including_original'))} h "
      '(cap 2 h).')
    w('')
    d = A['decision']
    w(f"## Decision under the amended rules: **{d['outcome']}**" + (f" — order {d['order']}" if d.get('order') else '')
      + (f" — {d['reason']}" if d.get('reason') else ''))
    w('')
    c = M['correctness']
    w('## 1. Correctness (separate subprocess)')
    w('')
    w(f"- C1′ closed_start reconstruction vs WP-B derived: max |Δ| {f(c['C1p']['max_abs_diff_A'])} Å.")
    w(f"- C2′ features: κ {f(c['C2p']['max_kappa_diff'])}, τ {f(c['C2p']['max_tau_circ_diff'])} rad.")
    w('- C3′ D1+ search: ' + '; '.join(f"{k}: {v['topk_agree']}/{v['queries']} top-10 agree, near ties {v['near_ties']}, "
                                      f"max |Δd| {f(v['max_abs_dist_diff'])}" for k, v in c['C3p'].items()) + '.')
    w(f"- C4′ instrumented search identical to the harness on CW4: {c['C4p']['instrumented_identical_to_harness_CW4']}.")
    w(f"- Pass: **{c['pass']}**.")
    w('')
    w('## 2. S4 models and confirmation (W = 8; per method and mode)')
    w('')
    w('| Model | coefficients | fit residual | fit ok | confirmation: workload measured → predicted (rel. err.) | supported | P2 value | P2 extrapolation × |')
    w('|---|---|---|---|---|---|---|---|')
    for k, v in A['fits'].items():
        if not k.startswith('S4_'):
            continue
        conf = '; '.join(f"{e['workload']} {f(e['measured'])} → {f(e['predicted'])} ({f(e['rel_error'], 2)})"
                         for e in v.get('confirmation', []))
        p2 = v.get('P2', {})
        w(f"| {k} | {', '.join(f'{a}={f(b)}' for a, b in v['coef'].items())} | {f(v['max_fit_rel_residual'], 2)} | "
          f"{v['fit_supported']} | {conf} | **{v['supported']}** | {f(p2.get('value'))} | "
          f"{f(p2.get('extrapolation_factor_vs_largest_measured'))} |")
    w('')
    w('Units: batch time = CPU s; single time = wall s; memory = absolute peak RSS bytes.')
    w('')
    w('## 3. Other stage fits (S1 reused from the original run with corrected memory predictor; S2/S3 re-measured)')
    w('')
    w('| Fit | supported | c0 | c1 | max rel. residual | P2 projection |')
    w('|---|---|---|---|---|---|')
    for k, v in A['fits'].items():
        if k.startswith('S4_'):
            continue
        p = A['projections_simple'].get(k)
        w(f"| {k} | {v['supported']} | {f(v.get('c0'))} | {f(v.get('c1'))} | {f(v.get('max_rel_residual'))} | "
          f"{f(p['value']) if p else '—'} |")
    w('')
    w('## 4. Requirements (method-specific composition; thresholds unchanged)')
    w('')
    w('| Requirement | component | projected | limit | unit | ratio |')
    w('|---|---|---|---|---|---|')
    for k, v in A['requirements'].items():
        w(f"| {k} | {v['component']} | {f(v['value'])} | {f(v['limit'])} | {v['unit']} | {f(A['ratios_descriptive'][k])} |")
    w('')
    w('## 5. S4 time decomposition (share of total; incidental harness slice copies vs required kernel/top-k work)')
    w('')
    w('| Cell | slice copies | kernel | top-k | other loop |')
    w('|---|---|---|---|---|')
    for k, v in A['S4_time_share_decomposition'].items():
        if '/fit/' in k and 'nD1015' not in k:
            continue
        w(f"| {k} | {f(v['slice_copy'], 2)} | {f(v['kernel'], 2)} | {f(v['topk'], 2)} | {f(v['loop_other'], 2)} |")
    w('')
    w('## 6. What remains valid from EB-1, what was corrected')
    w('')
    w('- **Valid:** original correctness checks C0–C3; S1 ingestion timings and peak RSS (process unaffected; predictor '
      'corrected only); the original outcome as a record of the original protocol.')
    w('- **Corrected:** S4 memory accounting (effective batch, n_D, method-specific transients, retained outputs, '
      'per method and mode); S2/S3 memory (owned inputs; no reference in the S3 process); method-specific paths (B2 '
      'without S3); confirmation workloads; time decomposition.')
    w('')
    w('## 7. Limitations')
    w('')
    w('- Shared desktop under memory pressure; page caches not dropped; timings conditional on contention.')
    w('- P2 values remain extrapolations beyond the largest confirmed workload (factors in §2–§3); confirmation tests '
      'prediction within the measured range plus moderate extrapolation, not at P2.')
    w('- Transient-byte coefficients come from a SYNTHETIC tracemalloc calibration at W = 8; RSS includes allocator '
      'retention that tracemalloc does not see.')
    w('- D1+ adds real closed_start frames for capacity only; it makes no retrieval-quality claim.')
    w('- Only existing implementations are measured; no saving is claimed for an implementation that does not exist.')
    w('')
    OUT.write_text('\n'.join(L) + '\n')
    print('wrote', OUT.name)


if __name__ == '__main__':
    main()
