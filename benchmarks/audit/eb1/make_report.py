#!/usr/bin/env python3
"""EB-1 report from results/eb1_measurements.json and results/eb1_analysis.json. Every number is read from those
files (AGENTS.md rule 3). Run: ../wpd_local_retrieval/build/venv/bin/python -B make_report.py
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
M = json.loads((HERE / 'results' / 'eb1_measurements.json').read_text())
A = json.loads((HERE / 'results' / 'eb1_analysis.json').read_text())
OUT = HERE / 'REPORT.md'


def f(x, d=3):
    return '—' if x is None else (f'{x:.{d}g}' if isinstance(x, float) else str(x))


def main():
    L = []
    w = L.append
    meta = M['_meta']
    w('# EB-1 engineering benchmark (development data only)')
    w('')
    w(f"Protocol `docs/EB1_ENGINEERING_BENCHMARK_PROTOCOL.md` (SHA-256 `{meta['protocol_sha256'][:16]}…`, frozen before "
      f"timing). {meta['label']}. Started {meta['started_utc']}; wall {meta.get('wall_s')} s; aggregate CPU "
      f"{f(meta.get('aggregate_cpu_h'))} h (budget 2 h). Cells: {meta.get('n_cells_ok')} of {meta.get('n_cells_planned')} "
      'completed. Engineering only: WP-D retrieval results are not re-estimated.')
    w('')
    d = A['decision']
    w(f"## Decision (frozen §7): **{d['outcome']}**" + (f" — order {d['order']}" if d.get('order') else '')
      + (f" — {d['reason']}" if d.get('reason') else ''))
    w('')
    w('## 1. Correctness at each processing boundary (mandatory)')
    w('')
    c = M['correctness']
    w(f"- C0 stored, D1: bit-identical {c['C0_D1']['bit_identical']} ({c['C0_D1']['reference']}; frames "
      f"{c['C0_D1']['frames']}). D2: bit-identical {c['C0_D2']['bit_identical']} ({c['C0_D2']['reference']}).")
    w(f"- C1 reconstructed: D1 max |Δ| {f(c['C1_D1']['max_abs_diff_A'])} Å vs {c['C1_D1']['reference']}; D2 max |Δ| "
      f"{f(c['C1_D2']['max_abs_diff_A'])} Å vs {c['C1_D2']['reference']}.")
    w('- C2 features (binding vs NumPy, identical float32-rounded inputs): ' + '; '.join(
        f"{k}: κ {f(v['max_kappa_diff'])}, τ {f(v['max_tau_circ_diff'])} rad" for k, v in c['C2'].items()))
    w('- C3 search (100 queries × window per cell; near ties by the frozen rule):')
    w('')
    w('| Cell | top-k agree | near ties | max |Δd| | tolerance | pass |')
    w('|---|---|---|---|---|---|')
    for k, v in c['C3'].items():
        w(f"| {k} | {v['topk_agree']}/{v['queries']} | {v['near_ties']} | {f(v['max_abs_dist_diff'])} | {v['dist_tol']} | {v['pass']} |")
    w('')
    w(f"Overall correctness pass: **{c['pass']}**.")
    w('')
    w('## 2. Measured stage costs (operational path; single process; BLAS threads 1)')
    w('')
    cells = [x for x in M['cells'] if x.get('status') == 'ok']
    w('| Block | Args | Work x | wall s (median) | CPU s | absolute peak RSS MB | MemAvailable at start MB |')
    w('|---|---|---|---|---|---|---|')
    groups = {}
    for x in cells:
        key = (x['block'], json.dumps(x['args'], sort_keys=True))
        groups.setdefault(key, []).append(x)
    for (b, a), xs in groups.items():
        w(f"| {b} | {a} | {f(xs[0].get('x'))} | {f(float(np.median([y['wall_s'] for y in xs])))} | "
          f"{f(float(np.median([y.get('cpu_s', y['wall_s']) for y in xs])))} | "
          f"{f(max(y['maxrss_kB'] for y in xs) / 1024)} | {f(min(y['mem_available_kB_at_start'] for y in xs) / 1024)} |")
    w('')
    s1 = [x for x in cells if x['block'] == 'S1']
    if s1:
        w('S1 first access vs repeats (wall s): ' + '; '.join(
            f"{x['args']} rep{x['rep']}: {f(x['wall_s'])}" for x in s1 if x['args'].get('n') in (2253, 4)))
        w('')
    w(f"Baseline import-only peak RSS: {f(A['_meta']['baseline_import_rss_kB'] / 1024 if A['_meta']['baseline_import_rss_kB'] else None)} MB "
      '(supplementary; decisions use absolute peak RSS).')
    w('')
    w('## 3. Fits and P2 projections (frozen §6; P2 = 1,000 residues, W = 8, 10,000 database frames, 1,000 batch queries)')
    w('')
    w('| Fit | supported | c0 | c1 | max rel. residual | P2 projection | extrapolation × |')
    w('|---|---|---|---|---|---|---|')
    for k, v in A['fits'].items():
        p = A['projections'].get(k)
        w(f"| {k} | {v['supported']} | {f(v.get('c0'))} | {f(v.get('c1'))} | {f(v.get('max_rel_residual'))} | "
          f"{f(p['value']) if p else '—'} | {f(p['extrapolation_factor']) if p else '—'} |")
    w('')
    w('## 4. Requirements (absolute)')
    w('')
    w('| Requirement | component | projected | limit | unit | ratio | violated |')
    w('|---|---|---|---|---|---|---|')
    for k, v in A['requirements'].items():
        w(f"| {k} | {v['component']} | {f(v['value'])} | {f(v['limit'])} | {v['unit']} | {f(v.get('ratio'))} | {v.get('violated')} |")
    w('')
    w('## 5. Research-only path (reported, not deciding)')
    w('')
    for k, v in A['research_only'].items():
        w(f"- {k}: wall {f(v['wall_s'])} s, CPU {f(v['cpu_s'])} s, peak RSS {f(v['maxrss_kB'] / 1024)} MB")
    w('')
    w(f"M8 parallel efficiency: {M.get('M8_parallel')}.")
    w('')
    w('## 6. Post-hoc observations (after results; NOT part of the frozen decision)')
    w('')
    f4 = A['fits']['R4_S4']
    s4 = [x for x in cells if x['block'] in ('S4single', 'S4batch')]
    w(f"- **Why R4_S4 is unsupported:** the S4 memory work variable `resident_bytes` was specified with the fixed batch "
      f"size 128 rather than each cell's actual query count, so it took only {len(set(f4['x']))} distinct values "
      f"(one per dataset) and failed the frozen ≥ 3-size rule. Measured absolute peak RSS nevertheless rose with the "
      f"batch query count (e.g. D1 W8 B2 batch: " + ', '.join(
          f"n_Q {x['args']['nq']} → {f(x['maxrss_kB'] / 1024)} MB" for x in s4
          if x['block'] == 'S4batch' and x['args'] == {**x['args'], 'd': 'D1', 'W': 8, 'm': 'B2'}) + '). '
      'This is a specification error in EB-1, not a measurement failure; the frozen outcome stands.')
    w(f"- Largest absolute peak RSS among all operational cells: {f(max(x['maxrss_kB'] for x in cells if x['block'] in ('S1','S2','S3','S4single','S4batch')) / 1024)} MB.")
    w('- Ratios of the supported projections to their limits (descriptive): ' + '; '.join(
        f"{k} {f(v['value'] / v['limit'])}" for k, v in A['requirements'].items() if v['value'] is not None) + '.')
    w('- Among supported projections, only R2_B2 (single-query, all-window latency for window Cα RMSD) exceeds its '
      'limit; no loader projection is close to its limit. Because one required projection is unsupported, this is not '
      'a decision.')
    w('')
    w('## 7. Limitations')
    w('')
    w('- Shared desktop under memory pressure (MemAvailable per cell above; swap full); timings are conditional on it and '
      'page caches could not be dropped (first access is labelled).')
    w('- C0 references use the same decoders (MDAnalysis, h5py): consistency checks, not independent decoders. C1 for D1 '
      'is a separate run of the same algorithm.')
    w('- P2 projections extrapolate linear fits far beyond the measured sizes (factors above); the linear models are '
      'analytically expected for these stages but are not measured at P2.')
    w('- Only existing implementations are measured; no saving is claimed for an implementation that does not exist.')
    w('- Two development datasets only; other formats and system sizes are untested.')
    w('')
    OUT.write_text('\n'.join(L) + '\n')
    print('wrote', OUT.name)


if __name__ == '__main__':
    main()
