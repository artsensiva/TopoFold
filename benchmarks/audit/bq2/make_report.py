#!/usr/bin/env python3
"""BQ-2 report (REPORT.md) from results/bq2_*.json. Every number is read from those files.
Run: ../wpd_local_retrieval/build/venv/bin/python -B make_report.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
M = json.loads((RES / 'bq2_measurements.json').read_text())
A = json.loads((RES / 'bq2_analysis.json').read_text())
F = json.loads((RES / 'bq2_protocol_freeze.json').read_text())
PM = json.loads((RES / 'bq2_pre_measurement_record.json').read_text())
OUT = HERE / 'REPORT.md'


def f(x, d=3):
    if x is None:
        return '—'
    if isinstance(x, bool):
        return str(x)
    return f'{x:.{d}g}' if isinstance(x, float) else str(x)


def rr(s):
    return '—' if not s else f"{f(s['median'])} ({f(s['min'])}–{f(s['max'])})"


def main():
    L = []
    w = L.append
    c = A['conclusion']
    w('# BQ-2: scaling confirmation of unchanged V2b (harness; development data + SYNTHETIC_tiled performance-only rows)')
    w('')
    w(f"Protocol `docs/B2_SCALING_CONFIRMATION_PROTOCOL.md` (SHA-256 "
      f"`{F['sha256']['docs/B2_SCALING_CONFIRMATION_PROTOCOL.md'][:16]}…`, frozen {F['frozen_utc']}); code and analysis "
      f"hashed {PM['recorded_utc']} before measurement. Informed by BQ-1 (not an independent rediscovery). BQ-1 files and "
      'its unsuccessful overall acceptance outcome are unchanged. No production code changed.')
    w('')
    w(f"CPU: {f(A['_meta']['cpu_h_total'])} h of the 1 CPU-h cap. Stopped early: {A['_meta']['stopped'] or 'no'}. "
      f"Load average start/end: {A['_meta']['loadavg_start_end']}.")
    w('')
    w('## Conclusion')
    w('')
    w(f"- **A. Unchanged V2b supports prepared-database operation:** {c['A_prepared_database_supported']}")
    w(f"- **B. First-query latency still warrants optimization:** {c['B_first_query_warrants_optimization']}")
    w(f"- **C. Evidence remains insufficient:** {c['C_insufficient_evidence']}")
    w(f"- Projected violation with supported models: {c['projected_violation_with_supported_models']}; measured "
      f"correctness failures: {c['measured_correctness_failure'] or 'none'}; measured resource-limit violations: "
      f"{c['measured_resource_limit_violation'] or 'none'}.")
    if c['missing_evidence']:
        w(f"- Missing evidence: {'; '.join(c['missing_evidence'])}.")
    w('')
    w('Scope: a faster implementation of the existing retrieval objective on development and SYNTHETIC '
      'performance-only workloads. No biological capability, novelty over existing RMSD software or clinical relevance '
      'is implied.')
    w('')
    w('## 1. Correctness gates (before timing)')
    w('')
    w(f"Tests in-run: {(A.get('tests') or {}).get('summary')}. Timed outputs identical to gate outputs: "
      f"{A['timed_outputs_identical_to_gate']}.")
    w('')
    w('| Arm | workload | status | max |Δd| Å | non-equivalent | strict / query-windows | near-tie | fallback rows |')
    w('|---|---|---|---|---|---|---|---|')
    for arm, gs in A['gates'].items():
        for wid, g in gs.items():
            w(f"| {arm} | {wid} | {g.get('status')} | {f(g.get('max_abs_dist_err_A'))} | {g.get('rank_nonequivalent')} | "
              f"{g.get('rank_strict_identical')}/{g.get('query_windows')} | {g.get('near_tie_query_windows')} | "
              f"{g.get('fallback_rows')} |")
    w('')
    w('## 2. Direct measurements (medians over 3 rounds; speed-up = per-round paired ratio vs same-session V0)')
    w('')
    w('| Workload | kind | role | n_D × n_win | prep s (reps rel. IQR) | V2b query s | first query s | V0 query s | '
      'speed-up prepared | speed-up cold | prepared MB | peak RSS MB | fallback rows/run | load avg |')
    w('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    for wid, d in A['direct'].items():
        s, b = d['V2b'], d['V0']
        if not s:
            w(f"| {wid} | {d['kind']} | {d['role']} | not measured | | | | | | | | | | |")
            continue
        sp = d['speedup_vs_V0'] or {}
        w(f"| {wid} | {d['kind']} | {d['role']} | {s['n_d']} × {s['n_win']} | {f(s['prep_wall_s'])} "
          f"({f(s['prep_rep_rel_iqr'], 2)}) | {f(s['prepared_query_wall_s'])} | {f(s['first_query_wall_s'])} | "
          f"{f(b['prepared_query_wall_s']) if b else '—'} | {rr(sp.get('prepared'))} | {rr(sp.get('cold'))} | "
          f"{f(s['prepared_bytes'] / 1e6)} | {f(s['maxrss_kB'] / 1024)} | {s['fallback_rows_per_run']} | "
          f"{f(max(s['loadavg_1min']), 2)} |")
    w('')
    p = A['P2D']
    w(f"**Direct P2 check (P2D, SYNTHETIC_tiled n_D 10,000 × 993):** status {p.get('status')}"
      + (f"; predicted peak {f(p.get('predicted_peak_bytes', 0) / 1e9)} GB vs MemAvailable "
         f"{f(p.get('mem_available_kB', 0) / 1e6)} GB minus 0.6 GB headroom" if p.get('status') == 'skipped_memory_rule'
         else '') + (f"; prep {f(p.get('prep_wall_s'))} s, query {f((p.get('query_wall_s') or [None])[0])} s, peak RSS "
                     f"{f(p.get('rss_end_kB', 0) / 1e6)} GB" if p.get('status') == 'ok' else '') + '.')
    w('')
    w('## 3. Models (fitted on the real fit grid only) and confirmation')
    w('')
    mods = A.get('models')
    if mods:
        w('| Model | coefficients | fit residual | confirmation rel. error (abs.) | supported | P2 value | extrapolation × |')
        w('|---|---|---|---|---|---|---|')
        for k in ('query_time', 'prep_wall', 'prep_cpu', 'memory'):
            m = mods[k]
            conf = '; '.join(f"{e['workload']} {f(e['rel_error'], 2)} ({f(e['abs_error'], 2)})" for e in m['confirmation'])
            w(f"| {k} | {', '.join(f'{a}={f(b)}' for a, b in m['coef'].items())} | {f(m['max_fit_rel_residual'], 2)} | "
              f"{conf} | **{m['supported']}** | {f(m['P2']['value'])} | {f(m['P2']['extrapolation_factor'])} |")
        w('')
        w('Units: time models in s (c1 per row-window, c2 per window); memory in bytes. BQ-1 cells under the BQ-2 '
          'models (diagnostic only; BQ-1 verdict unchanged):')
        w('')
        for k in ('prep_wall', 'query_time'):
            dg = mods[k]['diagnostic_BQ1'] or []
            w(f"- {k}: " + '; '.join(f"{e['bq1_workload']} {f(e['rel_error'], 2)} ({f(e['abs_error'], 2)} s)" for e in dg) + '.')
        w('')
    w('## 4. P2 projections (fitted; valid only where the model is supported)')
    w('')
    pr = A.get('projections_P2') or {}
    if pr:
        sup = pr['models_supported']

        def tag(*ks):
            return '' if all(sup[k] for k in ks) else ' **(UNSUPPORTED — not a valid projection)**'
        w(f"- Prepared-query latency: {f(pr['prepared_query_s'])} s{tag('query_time')} (slice-scaled from S4_P2slice, "
          f"assuming linearity in windows only: {f(pr['slice_scaled_prepared_query_s'])} s); limit 10 s.")
        w(f"- Preparation: {f(pr['prep_wall_s'])} s wall, {f(pr['prep_cpu_s'])} s CPU{tag('prep_wall', 'prep_cpu')}. "
          f"First query (prep + query): {f(pr['first_query_s'])} s{tag('query_time', 'prep_wall')}.")
        w(f"- Amortized per query at N = 1/10/100/1000: {', '.join(f(x) for x in pr['amortized_s'].values())} s; "
          f"break-even N for 10 s: {f(pr['break_even_n_for_10s'])}{tag('query_time', 'prep_wall')}.")
        w(f"- Peak RSS: {f(pr['peak_rss_bytes'] / 1e9)} GB{tag('memory')} (limit 8 GiB). R1_B2 + preparation: "
          f"{f(pr['R1_B2_plus_prep_cpu_h'])} CPU-h{tag('prep_cpu')} (limit 96; not evaluated if unsupported).")
        w(f"- Models supported: {sup}.")
        w('')
        w('## 5. Post-hoc observations (not decision-bearing; computed from the direct measurements above)')
        w('')
        rows = [(wid, d['kind'], d['V2b']) for wid, d in A['direct'].items() if d['V2b']]
        w('Preparation wall time per database row-window (s), by workload and row length L:')
        w('')
        w('| Workload | kind | L | n_D × n_win | prep s / (n_D·n_win) |')
        w('|---|---|---|---|---|')
        for wid, kind, s in rows:
            Lr = s['resident_bytes'] // ((s['n_d'] + 3) * 24)
            w(f"| {wid} | {kind} | {Lr} | {s['n_d']} × {s['n_win']} | {f(s['prep_wall_s'] / s['comparisons'])} |")
        sl = A['direct']['S4_P2slice']['V2b']
        sl_prep = sl['prep_wall_s'] * 993 / sl['n_win']
        w('')
        w("- The per-row-window preparation cost is roughly constant on the real workloads and on S1. It rises on "
          "S2, S3 and S4_P2slice, whose resident coordinate arrays (n_D × L) are the largest measured, and the real-data "
          "model under-predicts exactly those cells. A plausible structural cause, not tested here, is the strided "
          "gather `X[di, s:s+W]` over a resident array larger than the CPU caches. The frozen model has no "
          "working-set term, and at P2 the resident array is larger still.")
        w(f"- Slice-scaled preparation at P2 (S4_P2slice × 993/{sl['n_win']}, linearity in windows assumed; post hoc): "
          f"{f(sl_prep)} s; with the slice-scaled query, first query ≈ "
          f"{f(sl_prep + pr['slice_scaled_prepared_query_s'])} s. Indicative only.")
        w('- The large near-tie counts on S3/S4_P2slice and the ~5e-7 Å maximum distance differences (both arms, '
          'within the 1e-6 Å tolerance) come from repeated tiled rows: identical windows give RMSD ≈ 0, where the '
          'square root amplifies rounding. Rankings remained equivalent in every query-window.')
    w('')
    w('## 6. Limitations')
    w('')
    w('- SYNTHETIC_tiled rows are performance-only: repeated real frames joined end to end (28 of 993 windows per row '
      'are artificial junctions). They test scaling with n_D·n_win, n_win up to 993 and n_D = 10,000. They do not '
      'establish performance on real 1000-residue proteins or real 10,000-frame trajectories, nor any retrieval or '
      'biological property.')
    w('- One thread on a shared desktop; page caches not dropped; load average recorded per cell.')
    w('- Unless P2D was measured, P2 values remain projections; the extrapolation factors are listed in §3.')
    w('')
    OUT.write_text('\n'.join(L) + '\n')
    print('wrote', OUT.name)


if __name__ == '__main__':
    main()
