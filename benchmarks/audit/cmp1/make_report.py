#!/usr/bin/env python3
"""CMP-1 report (REPORT.md) from results/cmp1_*.json; every number is read from those files.
Run: build/venv/bin/python -B make_report.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
A = json.loads((RES / 'cmp1_analysis.json').read_text())
F = json.loads((RES / 'cmp1_protocol_freeze.json').read_text())
PM = json.loads((RES / 'cmp1_pre_measurement_record.json').read_text())
ENV = json.loads((RES / 'cmp1_environment.json').read_text())
DEV = json.loads((RES / 'cmp1_adapter_dev_checks.json').read_text())
OUT = HERE / 'REPORT.md'
IMPLS = ('T-V2b', 'T-V0', 'E-MDT', 'E-MDA')


def f(x, d=3):
    if x is None:
        return '—'
    if isinstance(x, bool):
        return str(x)
    return f'{x:.{d}g}' if isinstance(x, float) else str(x)


def main():
    L = []
    w = L.append
    w('# CMP-1: bounded competitive evaluation of window Cα RMSD retrieval (harness; development + SYNTHETIC data)')
    w('')
    w(f"Protocol `docs/CMP1_COMPETITIVE_EVALUATION_PROTOCOL.md` (SHA-256 "
      f"`{F['sha256']['docs/CMP1_COMPETITIVE_EVALUATION_PROTOCOL.md'][:16]}…`, frozen {F['frozen_utc']}); code hashed "
      f"{PM['recorded_utc']} before deciding measurements. Environment: Python {ENV['python']}; "
      + ', '.join(p for p in ENV['pip_freeze'] if p.split('==')[0] in ('mdtraj', 'MDAnalysis', 'numpy', 'scipy'))
      + '; isolated venv; wheel hashes in `results/cmp1_environment.json`. CPPTRAJ excluded (protocol §2). '
      'This is a bounded comparison of the listed implementations on the listed workloads, not a field-wide '
      'best-tool claim.')
    w('')
    w(f"CPU used: {f(A['_meta']['cpu_h_total'])} h of the 2 CPU-h cap. Stopped early: {A['_meta']['stopped'] or 'no'}. "
      f"Timed cells completed: {A['coverage']['timed_cells_ok']} of {A['coverage']['timed_cells_total']}.")
    w('')
    w('## Conclusion')
    w('')
    w(f"- **Primary (common F32 accuracy track): {A['conclusion_primary_F32']}** — per workload: "
      + '; '.join(f'{k} {v}' for k, v in A['outcomes_primary_F32'].items()) + '.')
    w(f"- Secondary (strict F64 track, float64-accurate implementations only): {A['conclusion_secondary_F64']} — "
      + '; '.join(f'{k} {v}' for k, v in A['outcomes_secondary_F64'].items()) + '. The secondary track cannot by '
      'itself justify integration (protocol §7).')
    w(f"- V2b degeneracy suite (F64): {A['V2b_degeneracy_pass_F64']}.")
    w('')
    w('## 1. Accuracy (verification outside timing; identical requirements per track)')
    w('')
    w('F64: |d − d_ref| ≤ 1e-6 Å, rank τ = 2e-6 Å. F32: |d − d_ref| ≤ 1e-3 Å, rank τ = 2e-3 Å. Reference: float64 SVD '
      'Kabsch on the canonical input. Native ratio = max |d² − d_ref²| / (u·(G_q + G_j)/W), descriptive only.')
    w('')
    w('| Workload | impl | F64 pass | F32 pass | max |Δd| Å | native ratio | rank non-eq. (F32) | strict identical (F32) | near-tie (F32) | input unchanged |')
    w('|---|---|---|---|---|---|---|---|---|---|')
    for wid, p in A['per_workload'].items():
        for i in IMPLS:
            x = p['impl'][i]
            t = x['tracks'] or {}
            if not t:
                w(f"| {wid} | {i} | {x['gate_status']} | | | | | | | |")
                continue
            w(f"| {wid} | {i} | {t['F64']['pass']} | {t['F32']['pass']} | {f(t['F64']['max_abs_err_A'])} | "
              f"{f(t['F64']['native_msd_ratio_descriptive'])} | {t['F32']['rank_nonequivalent']} | "
              f"{t['F32']['strict_identical']}/{t['F32']['windows']} | {t['F32']['near_tie']} | {x['input_unchanged_gate']} |")
    w('')
    w('**Competitor- and candidate-specific numerical limitations (listed separately; never sufficient for integration):**')
    w('')
    for lim in A['competitor_and_candidate_numerical_limitations']:
        extra = (f" failed {lim['failed']}; errors {', '.join(f'{k}={f(v)}' for k, v in lim['errors_A'].items())}"
                 if 'failed' in lim else f" max |Δd| {f(lim['max_abs_err_A'])} Å, rank non-eq. {lim['rank_nonequivalent']}")
        w(f"- {lim['impl']} / {lim['track']} / {lim['where']}:{extra}")
    if not A['competitor_and_candidate_numerical_limitations']:
        w('- none')
    w('')
    ms = DEV['mdtraj_semantics']
    w(f"MDTraj semantics (pinned source and adapter-development checks): distances are {ms['distances_dtype']}; "
      f"`precentered` with `atom_indices` is ignored with a warning ({ms['precentered_with_atom_indices_emits_warning']}); "
      f"the `atom_indices=None` path centres the input in place ({ms['atom_indices_None_path_centres_input_in_place']}); "
      f"pre-sliced precentered windows agree with subset RMSDs to {f(ms['pre_sliced_precentered_vs_subset_max_abs_nm'])} nm. "
      'MDTraj unit conversion (Å → nm float32 → Å) representation error on development cases: '
      + ', '.join(f"{k} {f(v['E-MDT_representation_error_A'])} Å" for k, v in DEV['cases'].items()) + '.')
    w('')
    w('## 2. Performance (fresh process per cell, one thread, 5 seeded rounds; medians, with ranges)')
    w('')
    w('L2 = complete retrieval workflow per query (all windows, shared top-k). L1 = solver-kernel time within the same '
      'queries (distance computation only). First-query span = preparation + query 1, excluding imports, warm-up and '
      'canonical input loading. Incremental peak RSS = lifetime peak − RSS after imports and warm-up.')
    w('')
    w('| Workload | impl | setup s | prepare s | first-query span s | repeated query L2 s (range) | kernel L1 s | peak RSS MB (incr.) | prepared MB | reps ok | outputs = gate |')
    w('|---|---|---|---|---|---|---|---|---|---|---|')
    for wid, p in A['per_workload'].items():
        for i in IMPLS:
            x = p['impl'][i]
            rr = x['repeated_query_range']
            w(f"| {wid} ({p['kind']}) | {i} | {f(x['setup_import_s'])} | {f(x['prepare_wall_s'])} | "
              f"{f(x['first_query_span_s'])} | {f(x['repeated_query_s'])}"
              + (f" ({f(rr[0])}–{f(rr[1])})" if rr else '') + f" | {f(x['kernel_L1_repeated_s'])} | "
              f"{f(x['peak_rss_kB'] / 1024 if x['peak_rss_kB'] else None)} "
              f"({f(x['incremental_peak_rss_kB'] / 1024 if x['incremental_peak_rss_kB'] else None)}) | "
              f"{f(x['prepared_bytes'] / 1e6 if x['prepared_bytes'] else None)} | {x['reps_ok']}/{x['reps_total']} | "
              f"{x['outputs_equal_gate']} |")
    w('')
    w('## 3. Per-workload comparison (primary F32 track)')
    w('')
    w('| Workload | E* | outcome | V2b/E* median repeated query | paired E*/V2b per round | memory criterion | ratios to V2b: L2, L1, first query, prepare, incr. memory (descriptive) |')
    w('|---|---|---|---|---|---|---|')
    for wid, p in A['per_workload'].items():
        o = p['outcome_F32_primary']
        rat = '; '.join(f"{i}: {f(v['repeated_query_s'], 2)}, {f(v['kernel_L1_repeated_s'], 2)}, "
                        f"{f(v['first_query_span_s'], 2)}, {f(v['prepare_wall_s'], 2)}, {f(v['incremental_peak_rss_kB'], 2)}"
                        for i, v in p['ratios_vs_V2b_descriptive'].items())
        w(f"| {wid} | {o.get('E_star') or '—'} | {o['outcome']}{' (' + o['reason'] + ')' if o.get('reason') else ''} | "
          f"{f(o.get('median_ratio_V2b_over_Estar'), 3)} | "
          f"{', '.join(f(x, 3) for x in o.get('paired_ratios_Estar_over_V2b', [])) or '—'} | {f(o.get('memory_ok'))} | {rat} |")
    w('')
    w('Ratios > 1 mean the other implementation takes longer (or uses more memory) than V2b.')
    w('')
    w('## 4. Limitations')
    w('')
    w('- Bounded comparison: four implementations, one machine, one thread, shared desktop. Not a field-wide '
      'best-tool claim. CPPTRAJ was not evaluated (protocol §2); multithreaded and GPU paths were not evaluated.')
    w('- Workloads are reused development data (AdK, mdCATH 1ha8A00) and SYNTHETIC_tiled performance-only rows. '
      'This is engineering evidence, not independent scientific validation.')
    w('- The F32 tolerance was re-based on input resolution after adapter-development checks showed the first '
      '(a-priori msd-bound) proposal did not model float32 coefficient error. This was disclosed and frozen before '
      'any deciding measurement (protocol §4).')
    w('- An existing library being preferable for this task would not be a failure of the broader TopoFold mission; '
      'TopoFold\'s scientific questions are separate from which RMSD engine it uses.')
    w('')
    OUT.write_text('\n'.join(L) + '\n')
    print('wrote', OUT.name)


if __name__ == '__main__':
    main()
