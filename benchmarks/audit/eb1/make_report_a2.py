#!/usr/bin/env python3
"""EB-1 Amendment 2 report (REPORT_A2.md) from results/eb1a2_*.json. Every number is read from those files.
EB-1 and Amendment 1 reports are not modified. Run: ../wpd_local_retrieval/build/venv/bin/python -B make_report_a2.py
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
M = json.loads((RES / 'eb1a2_measurements.json').read_text())
A = json.loads((RES / 'eb1a2_analysis.json').read_text())
F = json.loads((RES / 'eb1a2_freeze.json').read_text())
P = json.loads((RES / 'eb1a2_preflight.json').read_text())
A1 = json.loads((RES / 'eb1a1_analysis.json').read_text())
OUT = HERE / 'REPORT_A2.md'


def f(x, d=3):
    if x is None:
        return '—'
    return f'{x:.{d}g}' if isinstance(x, float) else str(x)


def main():
    L = []
    w = L.append
    meta = M['_meta']
    d = A['decision']
    w('# EB-1 Amendment 2: remaining defects')
    w('')
    w(f"Amendment `docs/EB1_AMENDMENT_2.md` (SHA-256 `{F['sha256']['docs/EB1_AMENDMENT_2.md'][:16]}…`, frozen "
      f"{F['frozen_utc']}, after the read-only preflight and before any Amendment 2 measurement). EB-1 and Amendment 1 "
      f"files and outcomes are unchanged (Amendment 1: `{A1['decision']['outcome']}`); Amendment 1 did not select "
      'engine work. Development data only.')
    w('')
    w(f"Cells: {meta['n_cells_ok']} of {meta['n_cells_planned']} completed. CPU: Amendment 2 {f(meta['amendment_cpu_h'])} h; "
      f"aggregate including EB-1 and Amendment 1 {f(meta['aggregate_cpu_h_including_prior'])} h (cap 2 h).")
    w('')
    w(f"## Decision under the frozen rule: **{d['outcome']}**" + (f" — order {d['order']}" if d.get('order') else '')
      + (f" — {d['reason']}" if d.get('reason') else ''))
    w('')
    w('This is the Amendment 2 outcome. It is not a retrospective outcome of EB-1 or Amendment 1. Under the frozen '
      'rule the engine flag comes from the violated engine requirements listed in §1; requirements within limits are '
      'projections, not guarantees.')
    w('')
    w('## 1. Validity by group')
    w('')
    w('| Group | valid | status | violated | requirements (projected / limit, ratio) |')
    w('|---|---|---|---|---|')
    for g, v in A['groups'].items():
        reqs = '; '.join(f"{k} {f(r['value'])}/{f(r['limit'])} {r['unit']} ({f(r['ratio'], 2)})"
                         for k, r in v['requirements'].items())
        w(f"| {g} | {v['valid']} | {v['status']} | {', '.join(v['violated']) or '—'} | {reqs} |")
    w('')
    w('## 2. Preflight (read-only; not a measurement)')
    w('')
    fr = P['frames']
    w(f"- D1 stored frames {fr['D1_stored_before_burn']} (no burn-in). D2 stored per trajectory "
      f"{fr['D2_stored_per_traj_before_burn']} = {fr['D2_stored_total_before_burn']} before burn-in; "
      f"{fr['D2_post_burn_per_traj']} = {fr['D2_post_burn_total']} after (array: {fr['D2_round_all_actual']}). "
      f"D1+ {fr['D1p_round']} rows.")
    w('- Invalid Amendment 1 requests: ' + '; '.join(f"{r['block']} {r['args']} requested {r.get('requested_frames')} "
                                                    f"of {r.get('available_frames')}" for r in P['amendment1_invalid'])
      + '.')
    cls = P['amendment1_memory_classification']
    w(f"- Amendment 1 memory results affected: {sum(x['memory_affected'] for x in cls)} of {len(cls)}; replaced by "
      f"re-measurement ({A['_meta']['memory_sources']['A2_remeasured']}); unaffected reused "
      f"{A['_meta']['memory_sources']['A1_unaffected_reused']}.")
    w('')
    w('## 3. Correctness')
    w('')
    lc = M['loadcheck']
    w(f"- New direct-read loaders vs Amendment 1 load paths: {sum(r['identical'] for r in lc['loads'])} of "
      f"{len(lc['loads'])} distinct inputs identical (dtype, shape, values); owned: "
      f"{all(r['owned'] for r in lc['loads'])}. FC1 search through `measure_reps` identical to the EB-1 harness: "
      f"{lc['FC1_search_identical_to_harness']}.")
    w(f"- Every measured cell's inputs owned (no backing array): {all(c.get('inputs_owned') for c in M['cells'])}.")
    c = A['correctness']
    w(f"- Reused: EB-1 correctness {c['original']}, Amendment 1 correctness {c['amendment1']}. Overall pass: **{c['pass']}**.")
    w('')
    w('## 4. A batch-time model')
    w('')
    v = A['fits']['S4_A_batch_time']
    ref = v['amendment1_model_for_reference']
    w(f"Amendment 1 model (with intercept), refitted for reference: coefficients "
      f"{', '.join(f'{k}={f(x)}' for k, x in ref['coef'].items())}; max fit residual {f(ref['max_fit_rel_residual'], 2)} "
      '(criterion 0.25).')
    w('')
    w(f"Replacement (frozen, chosen after seeing the Amendment 1 residuals): "
      f"{', '.join(f'{k}={f(x)}' for k, x in v['coef'].items())}; terms dropped: {v['dropped'] or 'none'}; max fit "
      f"residual **{f(v['max_fit_rel_residual'], 3)}** (criterion 0.25). Supported: **{v['supported']}**.")
    w('')
    w('| Fit cell n_Q / n_D / n_win | measured CPU s | predicted | abs. error s | rel. error |')
    w('|---|---|---|---|---|')
    for r in v['fit_residuals']:
        w(f"| {r['nq']} / {r['nd']} / {r['n_win']} | {f(r['measured'])} | {f(r['predicted'])} | "
          f"{f(r['abs_error'], 2)} | {f(r['rel_error'], 2)} |")
    w('')
    w('| Confirmation | role | measured CPU s | predicted | abs. error s | rel. error |')
    w('|---|---|---|---|---|---|')
    for e in v['confirmation']:
        w(f"| {e['workload']} ({e['data']}) | fresh, deciding | {f(e['measured'])} | {f(e['predicted'])} | "
          f"{f(e['abs_error'], 2)} | {f(e['rel_error'], 2)} |")
    for e in v['amendment1_confirmations_not_deciding']:
        w(f"| {e['workload']} ({e['data']}) | seen, not deciding | {f(e['measured'])} | {f(e['predicted'])} | "
          f"{f(e['abs_error'], 2)} | {f(e['rel_error'], 2)} |")
    w('')
    w(f"P2 value {f(v['P2']['value'])} CPU s, extrapolation × {f(v['P2']['extrapolation_factor_vs_largest_measured'])} "
      'beyond the largest measured workload.')
    w('')
    w('## 5. All S4 models (W = 8)')
    w('')
    w('| Model | fit residual | confirmation max rel. error | supported | P2 value | extrapolation × |')
    w('|---|---|---|---|---|---|')
    for k, v in A['fits'].items():
        if not k.startswith('S4_'):
            continue
        ce = max((e['rel_error'] for e in v['confirmation']), default=None)
        w(f"| {k} | {f(v['max_fit_rel_residual'], 2)} | {f(ce, 2)} | {v['supported']} | {f(v['P2']['value'])} | "
          f"{f(v['P2']['extrapolation_factor_vs_largest_measured'])} |")
    w('')
    w('Units: batch time = CPU s; single time = wall s; memory = absolute peak RSS bytes.')
    w('')
    w('## 6. Re-measured memory and reused times (descriptive)')
    w('')
    for b in ('S2', 'S3', 'S4'):
        z = [x['delta_kB'] for x in A['memory_change_A1_to_A2'] if x['block'] == b]
        if z:
            w(f"- {b}: peak RSS change A2 − A1 over {len(z)} matched cells: min {f(float(min(z)))} kB, median "
              f"{f(float(np.median(z)))} kB, max {f(float(max(z)))} kB.")
    r = [t['ratio_A2_over_A1'] for t in A['time_reuse_check_descriptive']]
    w(f"- Times in re-measured cells vs reused Amendment 1 times (ratio A2/A1, {len(r)} cells): min {f(float(min(r)))}, "
      f"median {f(float(np.median(r)))}, max {f(float(max(r)))}. Re-measured times do not replace reused ones.")
    w('')
    w('## 7. Requirements (thresholds unchanged)')
    w('')
    w('| Requirement | component | projected | limit | unit | ratio |')
    w('|---|---|---|---|---|---|')
    for k, v in A['requirements'].items():
        w(f"| {k} | {v['component']} | {f(v['value'])} | {f(v['limit'])} | {v['unit']} | "
          f"{f(A['ratios_descriptive'][k])} |")
    w('')
    w('## 8. Limitations')
    w('')
    w('- Shared desktop; page caches not dropped; timings conditional on contention.')
    w('- P2 values are extrapolations beyond the largest confirmed workload (factors above). Confirmation tests '
      'prediction within and moderately beyond the measured range, not at P2.')
    w('- The replacement A batch-time model was chosen after seeing the Amendment 1 residuals; only FC1–FC4 test it '
      'prospectively. Its fit residual is close to the 0.25 criterion.')
    w('- Memory transients use the SYNTHETIC Amendment 1 calibration; peak RSS includes allocator retention.')
    w('- The engine outcome rests on R2_B2 (single-query B2 at P2). No optimization exists or is measured here; the '
      'B2 single-query proposal is in `docs/EB1_B2_SINGLE_QUERY_OPTIMIZATION_PROPOSAL.md` (not implemented).')
    w('')
    OUT.write_text('\n'.join(L) + '\n')
    print('wrote', OUT.name)


if __name__ == '__main__':
    main()
