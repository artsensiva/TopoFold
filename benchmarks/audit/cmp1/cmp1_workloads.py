#!/usr/bin/env python3
"""CMP-1 workloads (read-only; writes results/cmp1_workloads.json). Real query frames are drawn (seed 20261020) from
query positions never used by EB-1, BQ-1, BQ-2, BQ-3 or the CMP-1 adapter-development checks.
Run: build/venv/bin/python -B cmp1_workloads.py
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
AUD = HERE.parent
B1 = json.loads((AUD / 'b2opt' / 'results' / 'b2opt_workloads.json').read_text())
B2 = json.loads((AUD / 'bq2' / 'results' / 'bq2_workloads.json').read_text())
DEV = {'D1': [557, 688, 720], 'D2': [263, 315, 371], 'tiled': [7, 8, 9]}
SEED = 20261020


def meta(d):
    src = AUD / 'eb1' / 'build' / ('cache_a1' if d == 'D1p' else 'cache')
    return json.loads((src / f'{d}_meta.json').read_text())


def main():
    used = {d: set(v) for d, v in B1['previously_used_queries'].items()}
    for w in B1['workloads'] + B2['workloads']:
        if w.get('kind', 'real') == 'real':
            used.setdefault(w['d'], set()).update(w['queries'])
    shared = set(used['D1']) | set(used['D1p']) | set(DEV['D1'])     # D1 and D1+ query positions are the same frames
    pool_d1 = sorted(set(meta('D1')['qi']) - shared)
    pool_d2 = sorted(set(meta('D2')['qi']) - used['D2'] - set(DEV['D2']))
    rng = np.random.default_rng(SEED)
    q_d1 = sorted(rng.choice(pool_d1, 6, replace=False).tolist())
    wl = [
        {'id': 'R1_D1', 'kind': 'real', 'd': 'D1', 'nd': 1015, 'nwin': 207, 'queries': q_d1[:3]},
        {'id': 'R2_D1plus', 'kind': 'real', 'd': 'D1p', 'nd': 2693, 'nwin': 207, 'queries': q_d1[3:]},
        {'id': 'R3_D2', 'kind': 'real', 'd': 'D2', 'nd': 405, 'nwin': 44,
         'queries': sorted(rng.choice(pool_d2, 3, replace=False).tolist())},
        {'id': 'S1_P2slice', 'kind': 'SYNTHETIC_tiled', 'd': 'D1p', 'nd': 10000, 'nwin': 100, 'queries': [20, 21, 22]},
        {'id': 'S2_P2D', 'kind': 'SYNTHETIC_tiled', 'd': 'D1p', 'nd': 10000, 'nwin': 993, 'queries': [23, 24, 25]}]
    out = {'seed': SEED, 'W': 8, 'k': 10, 'workloads': wl, 'required_for_conclusion': [w['id'] for w in wl],
           'pool_sizes': {'D1/D1+': len(pool_d1), 'D2': len(pool_d2)}, 'dev_queries_excluded': DEV,
           'note': ('real = reused development data (AdK Orädd 2021; mdCATH 1ha8A00 development domain), non-duplicated '
                    'MD frames; SYNTHETIC_tiled = performance-only repeated real frames (tiled query rows '
                    'QUERY_OFFSET + q); no held-out cohort domain is used')}
    (HERE / 'results').mkdir(exist_ok=True)
    (HERE / 'results' / 'cmp1_workloads.json').write_text(json.dumps(out, indent=1) + '\n')
    for w in wl:
        print(w)


if __name__ == '__main__':
    main()
