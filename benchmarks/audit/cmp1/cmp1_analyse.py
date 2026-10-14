#!/usr/bin/env python3
"""CMP-1 analysis (protocol §7–§8). Applies the frozen criteria (cmp1_criteria) per workload and track, lists
competitor-specific numerical limitations separately, and summarizes kernel (L1) and workflow (L2) results.
Writes results/cmp1_analysis.json. Run: build/venv/bin/python -B cmp1_analyse.py
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cmp1_criteria as C  # noqa: E402

RES = HERE / 'results'
IN = RES / 'cmp1_measurements.json'
OUT = RES / 'cmp1_analysis.json'
WL = json.loads((RES / 'cmp1_workloads.json').read_text())
IMPLS = ('T-V2b', 'T-V0', 'E-MDT', 'E-MDA')


def med(x):
    return float(np.median(x)) if len(x) else None


def rng(x):
    return [float(min(x)), float(max(x))] if len(x) else None


def main():
    M = json.loads(IN.read_text())
    cells = {w['id']: {i: {'gate': M['gate'].get(w['id'], {}).get(i),
                           'reps': [r for r in M['reps'] if r['wid'] == w['id'] and r['impl'] == i]}
                       for i in IMPLS} for w in WL['workloads']}
    per = {}
    for wid, cs in cells.items():
        per[wid] = {'kind': next(w['kind'] for w in WL['workloads'] if w['id'] == wid),
                    'outcome_F32_primary': C.outcome(cs, 'F32'), 'outcome_F64_secondary': C.outcome(cs, 'F64'),
                    'impl': {}}
        for i, c in cs.items():
            reps = [r for r in c['reps'] if r.get('status') == 'ok']
            g = c['gate'] or {}
            per[wid]['impl'][i] = {
                'gate_status': g.get('status'), 'tracks': g.get('tracks'),
                'input_unchanged_gate': g.get('input_unchanged'),
                'reps_ok': len(reps), 'reps_total': len(c['reps']),
                'setup_import_s': med([r['setup_import_s'] for r in reps]),
                'load_canonical_s': med([r['load_canonical_s'] for r in reps]),
                'prepare_wall_s': med([r['prepare_wall_s'] for r in reps]),
                'prepare_wall_range': rng([r['prepare_wall_s'] for r in reps]),
                'first_query_span_s': med([r['first_query_span_s'] for r in reps]),
                'first_query_span_range': rng([r['first_query_span_s'] for r in reps]),
                'repeated_query_s': med([q for r in reps for q in r['query_wall_s'][1:]]),
                'repeated_query_range': rng([q for r in reps for q in r['query_wall_s'][1:]]),
                'kernel_L1_repeated_s': med([k for r in reps for k in r['kernel_s'][1:]]),
                'query_cpu_over_wall': med([sum(r['query_cpu_s']) / sum(r['query_wall_s']) for r in reps]),
                'peak_rss_kB': med([r['peak_rss_kB'] for r in reps]),
                'incremental_peak_rss_kB': med([r['peak_rss_kB'] - r['rss_after_setup_kB'] for r in reps]),
                'prepared_bytes': reps[0]['prepared_bytes'] if reps else None,
                'retained_output_bytes': reps[0]['retained_output_bytes'] if reps else None,
                'fallback_rows': [r['fallback_rows'] for r in reps],
                'os_threads': sorted({r['threads']['os_threads'] for r in reps}),
                'input_unchanged_timing': all(r['input_unchanged'] for r in reps) if reps else None,
                'outputs_equal_gate': (all(r['topk_sha256'] == g.get('topk_sha256') for r in reps) if reps else None),
                'statuses': [r.get('status') for r in c['reps']]}
        v = per[wid]['impl']['T-V2b']
        per[wid]['ratios_vs_V2b_descriptive'] = {
            i: {k: (per[wid]['impl'][i][k] / v[k] if per[wid]['impl'][i][k] and v[k] else None)
                for k in ('repeated_query_s', 'kernel_L1_repeated_s', 'first_query_span_s', 'prepare_wall_s',
                          'incremental_peak_rss_kB')}
            for i in IMPLS if i != 'T-V2b'}
    req = WL['required_for_conclusion']
    deg = {i: (M['degeneracy'].get(i) or {}) for i in IMPLS}
    v2b_deg = bool((deg['T-V2b'].get('F64') or {}).get('pass'))
    concl, outs = C.conclusion({w: per[w]['outcome_F32_primary'] for w in per}, req, v2b_deg)
    concl64, outs64 = C.conclusion({w: per[w]['outcome_F64_secondary'] for w in per}, req, v2b_deg)
    limits = []
    for wid in per:
        for i in IMPLS:
            t = per[wid]['impl'][i]['tracks'] or {}
            for T in ('F64', 'F32'):
                if t and not t.get(T, {}).get('pass'):
                    limits.append({'impl': i, 'track': T, 'where': wid,
                                   'max_abs_err_A': t[T]['max_abs_err_A'],
                                   'rank_nonequivalent': t[T]['rank_nonequivalent']})
    for i in IMPLS:
        for T in ('F64', 'F32'):
            d = deg[i].get(T) or {}
            fails = [k for k, x in d.items() if x is False and k != 'pass']
            if fails:
                limits.append({'impl': i, 'track': T, 'where': 'SYNTHETIC degeneracy suite', 'failed': fails,
                               'errors_A': {k: x for k, x in d.items() if k.endswith('_err_A')}})
    out = {'_meta': {'producer': 'cmp1_analyse.py', 'measurements_sha256': hashlib.sha256(IN.read_bytes()).hexdigest(),
                     'cpu_h_total': M['_meta'].get('cpu_h_total'), 'stopped': M.get('stopped'),
                     'primary_track': 'F32', 'criteria': {'adv_ratio': C.ADV_RATIO, 'paired_min': C.PAIRED_MIN,
                                                          'mem_ratio': C.MEM_RATIO}},
           'tests': M.get('tests'), 'references': M.get('ref'),
           'conclusion_primary_F32': concl, 'outcomes_primary_F32': outs,
           'conclusion_secondary_F64': concl64, 'outcomes_secondary_F64': outs64,
           'V2b_degeneracy_pass_F64': v2b_deg, 'degeneracy': deg,
           'competitor_and_candidate_numerical_limitations': limits, 'per_workload': per,
           'coverage': {'timed_cells_ok': sum(r.get('status') == 'ok' for r in M['reps']),
                        'timed_cells_total': len(M['reps'])}}
    OUT.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(concl, outs)
    print('F64:', concl64, outs64)
    return 0


if __name__ == '__main__':
    sys.exit(main())
