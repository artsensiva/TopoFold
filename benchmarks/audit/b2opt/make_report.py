#!/usr/bin/env python3
"""BQ-1 report (REPORT.md) from results/b2opt_*.json. Every number is read from those files.
Run: ../wpd_local_retrieval/build/venv/bin/python -B make_report.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
M = json.loads((RES / 'b2opt_measurements.json').read_text())
A = json.loads((RES / 'b2opt_analysis.json').read_text())
F = json.loads((RES / 'b2opt_protocol_freeze.json').read_text())
PM = json.loads((RES / 'b2opt_pre_measurement_record.json').read_text())
WL = json.loads((RES / 'b2opt_workloads.json').read_text())
OUT = HERE / 'REPORT.md'
ROLE = {w['id']: w['role'] for w in WL['workloads']}


def f(x, d=3):
    if x is None:
        return '—'
    if isinstance(x, bool):
        return str(x)
    return f'{x:.{d}g}' if isinstance(x, float) else str(x)


def rng_(s):
    return '—' if s is None else f"{f(s['median'])} ({f(s['min'])}–{f(s['max'])})"


def main():
    L = []
    w = L.append
    w('# BQ-1: B2 single-query optimization experiment (harness only, development data)')
    w('')
    w(f"Protocol `docs/B2_SINGLE_QUERY_OPTIMIZATION_PROTOCOL.md` (SHA-256 "
      f"`{F['sha256']['docs/B2_SINGLE_QUERY_OPTIMIZATION_PROTOCOL.md'][:16]}…`, frozen {F['frozen_utc']}), Amendment 1 "
      f"(`…AMENDMENT_1.md`, `{PM['sha256']['docs/B2_SINGLE_QUERY_OPTIMIZATION_AMENDMENT_1.md'][:16]}…`), code hashed "
      f"{PM['recorded_utc']} before any real-data computation. EB-1 versions and outcomes are unchanged. No production "
      'code was changed.')
    w('')
    w(f"CPU used: {f(A['_meta']['cpu_h_total'])} h of the 1 CPU-h cap (tests, references, gates, calibration, "
      f"preparation and timing). Stopped early: {A['_meta']['stopped'] or 'no'}.")
    w('')
    w(f"## Outcome: **{A['outcome']}**")
    w('')
    w(f"Strict (cold first query, preparation included) qualifying arms: {A['qualifying_cold'] or 'none'}. "
      f"Amendment 2 R2_B2 for reference: {f(A['amendment2_R2_B2_s'])} s.")
    w('')
    w('| Arm | status | gate | (a) gate + identity | (b) query ≤ 10 s | (c) prep models | (d) memory ≤ 8 GB | '
      '(e) R1 + prep ≤ 96 h | R2 (prepared DB) | R2 (cold) |')
    w('|---|---|---|---|---|---|---|---|---|---|')
    for arm, v in A['arms'].items():
        c = v.get('criteria', {})
        w(f"| {arm} | {v['status']} | {f((v.get('gate') or {}).get('pass'))} | "
          + ' | '.join(f(c.get(k)) for k in ('a_gate_and_identity', 'b_query_supported_and_le_10s',
                                              'c_prep_models_supported', 'd_memory_supported_and_le_8GB',
                                              'e_R1_B2_plus_prep_le_96h'))
          + f" | **{f(v.get('R2_supported_prepared_database'))}** | {f(v.get('R2_supported_cold'))} |")
    w('')
    w('Failed criteria per timed arm (frozen rules; nothing relaxed):')
    w('')
    for arm, v in A['arms'].items():
        if 'criteria' not in v:
            if v['status'] != 'timed':
                w(f"- {arm}: {v['status']}.")
            continue
        why = []
        for name in ('query_time', 'prep_wall', 'prep_cpu', 'memory'):
            mm = v['models'][name]
            sup = mm.get('supported_full', mm.get('supported'))
            if not sup:
                res = mm.get('max_fit_rel_residual', mm.get('max_rel_residual'))
                cell = ''
                if 'x' in mm and 'c0' in mm:
                    errs = [abs(y - (mm['c0'] + mm['c1'] * x)) / y for x, y in zip(mm['x'], mm['y'])]
                    i = max(range(len(errs)), key=errs.__getitem__)
                    cell = (f"; worst fit cell n_D·n_win = {int(mm['x'][i])}, measured {f(mm['y'][i])} s, "
                            f"residual {f(errs[i], 2)}")
                conf = max((e['rel_error'] for e in mm.get('confirmation', []) if e['deciding']), default=None)
                why.append(f"{name} unsupported (fit residual {f(res, 2)} > 0.25{cell}; max fresh confirmation error "
                           f"{f(conf, 2)})")
        p = v['projections_P2']
        if p['prepared_query_s'] is not None and p['prepared_query_s'] > 10:
            why.append(f"projected prepared-query latency {f(p['prepared_query_s'])} s > 10 s")
        failed = [k for k, ok in v['criteria'].items() if not ok]
        w(f"- {arm}: failed {', '.join(failed) or 'none'}; " + ('; '.join(why) or '—') + '.')
    w('')
    w('## 1. Checks before timing')
    w('')
    t = M.get('tests') or {}
    w(f"- Regression tests (in-run): {t.get('summary')} (pass {t.get('pass')}).")
    rc = M.get('refcheck') or {}
    if rc.get('workloads'):
        e1 = max(v['max_vs_kabsch_explicit'] for v in rc['workloads'].values())
        e2 = max(v['max_vs_mdanalysis_qcprot'] for v in rc['workloads'].values())
        w(f"- Reference R_svd vs explicit-rotation Kabsch: max {f(e1)} Å; vs MDAnalysis qcprot: max {f(e2)} Å "
          f"(2,000 sampled pairs per workload). Pass {rc.get('pass')}.")
    g0 = M.get('gate_V0') or {}
    w(f"- V0 (harness) gate: {g0.get('pass')}.")
    w('')
    w('## 2. Correctness gates per arm (all workloads, 3 queries each, every window)')
    w('')
    w('| Arm | workload | role | max |Δd| Å | non-equivalent | strict identical / query-windows | near-tie | fallback rows (fraction) |')
    w('|---|---|---|---|---|---|---|---|')
    gates = {'V0': g0, **{a: (M['arms'].get(a) or {}).get('gate') or {} for a in ('V1', 'V2', 'V2b')}}
    for arm, g in gates.items():
        for wid, v in (g.get('workloads') or {}).items():
            w(f"| {arm} | {wid} | {ROLE[wid]} | {f(v['max_abs_dist_err_A'])} | {v['rank_nonequivalent']} | "
              f"{v['rank_strict_identical']}/{v['query_windows']} | {v['near_tie_query_windows']} | "
              f"{v['fallback_rows']} ({f(v['fallback_fraction'])}) |")
    w('')
    w('SYNTHETIC adversarial suite (protocol §5.3):')
    w('')
    for arm, g in gates.items():
        s = g.get('synthetic')
        if s:
            w(f"- {arm}: pass {s['pass']}; " + '; '.join(
                f"{k} {v['pass']}" + (f" (max {f(v['max_err'])})" if 'max_err' in v else '')
                + (f", fallback rows {v['fallback_rows']}" if 'fallback_rows' in v else '')
                for k, v in s.items() if k != 'pass'))
    w('')
    w('## 3. Direct measurements (same-session V0; medians over rounds; speed-up = per-round paired ratio, median (min–max))')
    w('')
    for arm, v in A['arms'].items():
        if 'direct' not in v:
            continue
        w(f"### {arm} (timed outputs identical to gate outputs: {v['timed_outputs_identical_to_gate']})")
        w('')
        w('| Workload | role | n_D × n_win | V0 query s | arm query s | prep s | first query s | speed-up prepared | '
          'speed-up cold | prepared MB | peak RSS MB arm / V0 | fallback rows/run |')
        w('|---|---|---|---|---|---|---|---|---|---|---|---|')
        for wid, d in v['direct'].items():
            a, b = d['arm'], d['V0']
            w(f"| {wid} | {ROLE[wid]} | {a['n_d']} × {a['n_win']} | {f(b['prepared_query_wall_s'])} | "
              f"{f(a['prepared_query_wall_s'])} | {f(a['prep_wall_s'])} | {f(a['first_query_wall_s'])} | "
              f"{rng_(d['speedup_vs_V0']['prepared'])} | {rng_(d['speedup_vs_V0']['cold'])} | "
              f"{f(a['prepared_bytes'] / 1e6)} | {f(a['maxrss_kB'] / 1024)} / {f(b['maxrss_kB'] / 1024)} | "
              f"{a['fallback_rows_per_run']} |")
        cw = [d['arm']['cpu_over_wall_query'] for d in v['direct'].values()]
        w('')
        w(f"Query CPU/wall ratio (one thread expected ≈ 1): {f(min(cw))}–{f(max(cw))}.")
        w('')
    w('## 4. Models, confirmation and P2 projections (fitted; not measurements)')
    w('')
    for arm, v in A['arms'].items():
        if 'models' not in v:
            continue
        w(f"### {arm}")
        w('')
        w('| Model | fit residual | supported | confirmation (fresh = deciding): workload rel. error (abs. error) | P2 value | extrapolation × |')
        w('|---|---|---|---|---|---|---|')
        for name, mm in v['models'].items():
            res = mm.get('max_fit_rel_residual', mm.get('max_rel_residual'))
            sup = mm.get('supported_full', mm.get('supported'))
            conf = '; '.join(f"{e['workload']}{'' if e['deciding'] else ' (reused)'} {f(e['rel_error'], 2)} "
                             f"({f(e['abs_error'], 2)})" for e in mm.get('confirmation', []))
            p2 = mm.get('P2', {})
            val = p2.get('value', p2.get('value_bytes'))
            w(f"| {name} | {f(res, 2)} | {sup} | {conf} | {f(val)} | {f(p2.get('extrapolation_factor'))} |")
        w('')
        p = v['projections_P2']
        w(f"P2: prepared-query {f(p['prepared_query_s'])} s; preparation {f(p['prep_wall_s'])} s wall, "
          f"{f(p['prep_cpu_s'])} s CPU; cold first query {f(p['first_query_cold_s'])} s; peak RSS "
          f"{f(None if p['peak_rss_bytes'] is None else p['peak_rss_bytes'] / 1e9)} GB; R1_B2 + preparation "
          f"{f(p['R1_B2_plus_prep_cpu_h'])} CPU-h; same-session V0 {f(p['V0_same_session_query_s'])} s. "
          f"Amortized per query at N = 1/10/100/1000: "
          + (', '.join(f(x) for x in p['amortized_s'].values()) if p['amortized_s'] else '—')
          + f" s; break-even N for 10 s: {f(p['break_even_n_for_10s'])}.")
        w('')
    w('Units: time = wall s; memory = absolute peak RSS bytes. Time models: c0 + c1·(n_D·n_win) + c2·n_win; '
      'preparation: p0 + p1·(n_D·n_win); memory: m0 + m1·P.')
    w('')
    w('## 5. Limitations')
    w('')
    w('- Shared desktop, one thread, page caches not dropped; timings are conditional on contention.')
    w('- P2 values are extrapolations (factors above), confirmed only on the development workloads.')
    w('- Arms are NumPy harness implementations. The result says nothing about compiled or other B2 '
      'implementations, and does not generalize to all B2 implementations.')
    w('- V2b was added before measurement (Amendment 1) after a SYNTHETIC test exposed frozen V2 on degenerate '
      'windows. Its fallback constant was fixed from an error estimate, not tuned.')
    w('- Preparation is treated as a database-side stage under the primary verdict. The strict cold verdict counts it '
      'in the query latency.')
    w('')
    OUT.write_text('\n'.join(L) + '\n')
    print('wrote', OUT.name)


if __name__ == '__main__':
    main()
