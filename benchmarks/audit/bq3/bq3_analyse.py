#!/usr/bin/env python3
"""BQ-3 analysis (protocol §3): applies the frozen acceptance rule (bq3_lib.decide) and summarizes every repetition.
Direct measurements only; the EB-1 R1 projection is shown separately and does not decide.
Writes results/bq3_SYNTHETIC_analysis.json. Run: ../wpd_local_retrieval/build/venv/bin/python -B bq3_analyse.py
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bq3_lib as Q  # noqa: E402

RES = HERE / 'results'
IN = RES / 'bq3_SYNTHETIC_measurements.json'
OUT = RES / 'bq3_SYNTHETIC_analysis.json'
A2 = json.loads((HERE.parent / 'eb1' / 'results' / 'eb1a2_analysis.json').read_text())


def summ(x):
    x = [float(v) for v in x]
    return {'n': len(x), 'median': float(np.median(x)), 'min': min(x), 'max': max(x), 'values': x} if x else None


def main():
    M = json.loads(IN.read_text())
    reps = M.get('reps', [])
    d = Q.decide(M.get('ref'), M.get('gate'), reps, M['_meta']['n_reps'])
    ok = [r for r in reps if r.get('status') == 'ok']
    direct = {
        'P1_initial_prep_wall_s': summ([r['P1_wall_s'] for r in ok]),
        'P1_initial_prep_cpu_s': summ([r['P1_cpu_s'] for r in ok]),
        'Q1_first_query_wall_s': summ([r['query_wall_s'][0] for r in ok]),
        'first_query_span_s': summ([r['first_query_span_s'] for r in ok]),
        'subsequent_query_wall_s': summ([q for r in ok for q in r['query_wall_s'][1:]]),
        'all_query_wall_s': summ([q for r in ok for q in r['query_wall_s']]),
        'repeated_prep_wall_s': summ([p for r in ok for p in r['P_repeat_wall_s']]),
        'repeated_prep_cpu_s': summ([p for r in ok for p in r['P_repeat_cpu_s']]),
        'query_cpu_over_wall': summ([sum(r['query_cpu_s']) / sum(r['query_wall_s']) for r in ok]),
        'peak_rss_kB': summ([r['peak_rss_kB'] for r in ok]),
        'build_rows_wall_s': summ([r['build_rows_wall_s'] for r in ok]),
        'fallback_rows_per_rep': [r['fallback_rows'] for r in ok],
        'os_threads': [r['os_threads'] for r in ok],
        'threadpool': ok[0]['threadpool'] if ok else None,
        'prepared_bytes': ok[0]['prepared_bytes'] if ok else None,
        'resident_bytes': ok[0]['resident_bytes'] if ok else None,
        'contention': [{'rep': r['rep'], 'start': r['start'], 'end': r['end']} for r in ok]}
    pc = direct['P1_initial_prep_cpu_s']
    mixed = {'R1_B2_EB1_A2_projection_cpu_h': A2['requirements']['R1_B2']['value'],
             'direct_P2_prep_cpu_s_median': pc['median'] if pc else None,
             'sum_cpu_h': (A2['requirements']['R1_B2']['value'] + pc['median'] / 3600) if pc else None,
             'note': 'one fitted EB-1 projection plus one direct BQ-3 measurement; descriptive, not deciding'}
    out = {'_meta': {'producer': 'bq3_analyse.py', 'measurements_sha256': hashlib.sha256(IN.read_bytes()).hexdigest(),
                     'label': 'SYNTHETIC_tiled performance-only', 'limits': {'R2_s': Q.R2_S, 'R4_kB': Q.R4_KB},
                     'cpu_h_total': M['_meta'].get('cpu_h_total'), 'stopped': M.get('stopped')},
           'decision': d, 'tests': M.get('tests'),
           'reference': {k: v for k, v in (M.get('ref') or {}).items() if k != 'stderr'},
           'gate': {k: v for k, v in (M.get('gate') or {}).items() if k != 'stderr'},
           'repetitions': [{k: v for k, v in r.items() if k not in ('threadpool',)} for r in reps],
           'direct_summary': direct, 'R1_mixed_descriptive': mixed}
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({k: d[k] for k in ('prepared_database', 'first_query', 'verification_ok', 'complete')}))
    return 0


if __name__ == '__main__':
    sys.exit(main())
