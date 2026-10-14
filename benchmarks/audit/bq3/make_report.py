#!/usr/bin/env python3
"""BQ-3 report (REPORT_SYNTHETIC_P2.md) from results/bq3_SYNTHETIC_*.json; every number is read from those files.
Run: ../wpd_local_retrieval/build/venv/bin/python -B make_report.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
A = json.loads((RES / 'bq3_SYNTHETIC_analysis.json').read_text())
F = json.loads((RES / 'bq3_protocol_freeze.json').read_text())
PM = json.loads((RES / 'bq3_pre_measurement_record.json').read_text())
OUT = HERE / 'REPORT_SYNTHETIC_P2.md'


def f(x, d=3):
    if x is None:
        return '—'
    return f'{x:.{d}g}' if isinstance(x, float) else str(x)


def s(v, scale=1.0, unit=''):
    if not v:
        return '—'
    return f"median {f(v['median'] * scale)}{unit} (min {f(v['min'] * scale)}, max {f(v['max'] * scale)}; n = {v['n']})"


def main():
    d, g, ds = A['decision'], A['gate'], A['direct_summary']
    L = []
    w = L.append
    w('# BQ-3: direct P2 qualification of unchanged V2b (SYNTHETIC_tiled, performance-only)')
    w('')
    w(f"Protocol `docs/B2_DIRECT_P2_QUALIFICATION_PROTOCOL.md` (SHA-256 "
      f"`{F['sha256']['docs/B2_DIRECT_P2_QUALIFICATION_PROTOCOL.md'][:16]}…`, frozen {F['frozen_utc']}, starting "
      f"checkpoint `{F['starting_checkpoint'][:7]}`); code hashed {PM['recorded_utc']} before measurement. Harness "
      'only. BQ-1, BQ-2 (conclusion C) and EB-1 records are unchanged. Not pushed.')
    w('')
    w('**Workload:** n_D = 10,000, L = 1,000, W = 8, 993 windows, 3 queries. Rows are real D1+ frames joined end to end '
      '(SYNTHETIC_tiled, repeated after 3,706 rows). This tests computational scale at P2 dimensions only. It says '
      'nothing about real 1,000-residue proteins or 10,000-frame trajectories, retrieval quality or biology.')
    w('')
    w(f"CPU used: {f(A['_meta']['cpu_h_total'])} h of the 1 CPU-h cap. Stopped early: {A['_meta']['stopped'] or 'no'}.")
    w('')
    w('## Verdicts (frozen rule: every repetition must pass)')
    w('')
    w(f"- **Prepared-database operation:** {d['prepared_database']}")
    w(f"- **First-query operation (preparation + first query):** {d['first_query']}")
    w(f"- Verification: {d['verification_ok']} {d['verification_reasons'] or ''}; complete: {d['complete']} "
      f"{d['completeness_reasons'] or ''}; peak RSS within 8 GiB in every repetition: {d['memory_within_R4_all_reps']}.")
    w('')
    w('## 1. Verification (outside timing)')
    w('')
    w(f"- Failure-path tests: {(A.get('tests') or {}).get('summary')}.")
    w(f"- Reference: status {A['reference'].get('status')}, finite {A['reference'].get('finite')}.")
    w(f"- Gate: status {g.get('status')}; candidate finite {g.get('candidate_finite')}; max |Δd| "
      f"{f(g.get('max_abs_dist_err_A'))} Å (tolerance 1e-6); rank non-equivalent {g.get('rank_nonequivalent')} of "
      f"{g.get('query_windows')} query-windows; strict identical {g.get('rank_strict_identical')}; near-tie "
      f"{g.get('near_tie_query_windows')}; fallback rows {g.get('fallback_rows')}; gate peak RSS "
      f"{f((g.get('peak_rss_kB') or 0) / 1024 ** 2)} GiB.")
    w('')
    w('## 2. Direct measurements (all repetitions)')
    w('')
    w('| Rep | start (UTC) | P1 s | Q1 s | first-query span s | Q2 s | Q3 s | repeated prep s | peak RSS GiB | '
      'fallback rows | load avg start | MemAvailable start GB | SwapFree start MB | outputs = gate |')
    w('|---|---|---|---|---|---|---|---|---|---|---|---|---|---|')
    gh = g.get('topk_sha256')
    for r in A['repetitions']:
        if r.get('status') != 'ok':
            w(f"| {r['rep']} | {r.get('started_utc', '—')} | {r.get('status')} | | | | | | | | | | | |")
            continue
        q = r['query_wall_s']
        w(f"| {r['rep']} | {r['started_utc']} | {f(r['P1_wall_s'])} | {f(q[0])} | {f(r['first_query_span_s'])} | "
          f"{f(q[1])} | {f(q[2])} | {', '.join(f(x) for x in r['P_repeat_wall_s'])} | "
          f"{f(r['peak_rss_kB'] / 1024 ** 2)} | {r['fallback_rows']} | {r['start']['loadavg'][0]} | "
          f"{f(r['start']['MemAvailable_kB'] / 1e6)} | {f(r['start']['SwapFree_kB'] / 1e3)} | {r['topk_sha256'] == gh} |")
    w('')
    w(f"- Initial preparation P1: {s(ds['P1_initial_prep_wall_s'], unit=' s')}; CPU {s(ds['P1_initial_prep_cpu_s'], unit=' s')}.")
    w(f"- First query Q1: {s(ds['Q1_first_query_wall_s'], unit=' s')}. First-query span (measured end to end): "
      f"{s(ds['first_query_span_s'], unit=' s')}.")
    w(f"- Subsequent queries Q2–Q3: {s(ds['subsequent_query_wall_s'], unit=' s')}. All queries: "
      f"{s(ds['all_query_wall_s'], unit=' s')}.")
    w(f"- Repeated preparations: {s(ds['repeated_prep_wall_s'], unit=' s')}; CPU {s(ds['repeated_prep_cpu_s'], unit=' s')}.")
    w(f"- Peak RSS: {s(ds['peak_rss_kB'], 1 / 1024 ** 2, ' GiB')}; prepared arrays {f(ds['prepared_bytes'] / 1e9)} GB, "
      f"resident rows {f(ds['resident_bytes'] / 1e9)} GB.")
    w(f"- Query CPU/wall: {s(ds['query_cpu_over_wall'])}. OS threads per process: {ds['os_threads']}; threadpools: "
      f"{ds['threadpool']}.")
    w('')
    m = A['R1_mixed_descriptive']
    w('## 3. R1 (not measured here)')
    w('')
    w(f"R1_B2 remains the EB-1 Amendment 2 batch **projection** ({f(m['R1_B2_EB1_A2_projection_cpu_h'])} CPU-h). "
      f"Adding the directly measured P2 preparation CPU (median {f(m['direct_P2_prep_cpu_s_median'])} s) gives "
      f"{f(m['sum_cpu_h'])} CPU-h against the 96 CPU-h limit. This mixes a projection with a measurement, is "
      'descriptive only, and does not decide.')
    w('')
    w('## 4. Limitations')
    w('')
    w('- SYNTHETIC_tiled rows are performance-only: repeated real frames, with artificial junction windows. Fallback '
      'frequency, cache behaviour and timings may differ on real P2-scale data.')
    w(f"- Ranking verification is weak on this workload: {g.get('near_tie_query_windows')} of {g.get('query_windows')} "
      "query-windows contain near ties, because repeated tiled rows create exact duplicate windows. Rankings were "
      f"nevertheless strictly identical to the reference in {g.get('rank_strict_identical')} query-windows. Distance "
      "differences up to the reported maximum sit at the square-root rounding floor of RMSD ≈ 0 duplicates.")
    w('- One machine, one thread, shared desktop with swap nearly exhausted (recorded per repetition). Page caches '
      'were not dropped.')
    w('- Five repetitions characterize this machine and session only.')
    w('- Disclosure: after the measurement, `make_report.py` gained the ranking-verification limitation above (report '
      'text only). The acceptance code (`bq3_lib.py`, `bq3_analyse.py`) is unchanged from its pre-measurement hash.')
    w('- A qualification here is an engineering result for a harness implementation. It establishes no biological '
      'capability, no novelty over existing RMSD software and no production readiness.')
    w('')
    OUT.write_text('\n'.join(L) + '\n')
    print('wrote', OUT.name)


if __name__ == '__main__':
    main()
