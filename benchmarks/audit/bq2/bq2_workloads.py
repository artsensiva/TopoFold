#!/usr/bin/env python3
"""BQ-2 workload definition (read-only; no timing, no RMSD computation). Writes results/bq2_workloads.json.

Real query frames are drawn (seed 20261018) from the BQ-1 eligible pools minus every query BQ-1 used, so no BQ-2 real
query frame was used by EB-1, its amendments or BQ-1. SYNTHETIC stress rows are tiled real frames (performance-only).
Run: ../wpd_local_retrieval/build/venv/bin/python -B bq2_workloads.py
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
B1 = json.loads((HERE.parent / 'b2opt' / 'results' / 'b2opt_workloads.json').read_text())
SEED = 20261018


def eligible(d):
    m = json.loads((HERE.parent / 'eb1' / 'build' / ('cache_a1' if d == 'D1p' else 'cache') / f'{d}_meta.json').read_text())
    used = set(B1['previously_used_queries'][d])
    used |= {q for w in B1['workloads'] if w['d'] in (d, 'D1' if d == 'D1p' else d) for q in w['queries']}
    if d == 'D1p':
        used |= set(B1['previously_used_queries']['D1'])
    return sorted(set(m['qi']) - used)


def main():
    rng = np.random.default_rng(SEED)
    pools = {d: eligible(d) for d in ('D1p', 'D2')}
    taken = {d: set() for d in pools}

    def pick(d):
        qs = sorted(rng.choice([q for q in pools[d] if q not in taken[d]], 3, replace=False).tolist())
        taken[d] |= set(qs)
        return qs
    qfit = pick('D1p')
    wl = []
    for nd in (100, 700, 2693):
        for nw in (5, 20, 60, 207):
            wl.append({'id': f'F_nD{nd}_nwin{nw}', 'kind': 'real', 'd': 'D1p', 'nd': nd, 'nwin': nw, 'role': 'fit',
                       'queries': qfit, 'baseline': False})
    for cid, d, nd, nw in (('RC1', 'D1p', 1500, 120), ('RC2', 'D1p', 2693, 207), ('RC3', 'D2', 405, 44)):
        wl.append({'id': cid, 'kind': 'real', 'd': d, 'nd': nd, 'nwin': nw, 'role': 'confirm_real', 'queries': pick(d),
                   'baseline': True})
    for cid, nd, nw, role in (('S1', 1000, 993, 'confirm_stress'), ('S2', 2500, 993, 'confirm_stress'),
                              ('S4_P2slice', 10000, 100, 'confirm_stress'),
                              ('S3', 4000, 993, 'confirm_stress_if_memory'), ('P2D', 10000, 993, 'direct_P2_if_memory')):
        wl.append({'id': cid, 'kind': 'SYNTHETIC_tiled', 'd': 'D1p', 'nd': nd, 'nwin': nw, 'role': role,
                   'queries': [0, 1, 2], 'baseline': role != 'direct_P2_if_memory'})
    out = {'W': 8, 'k': 10, 'seed': SEED, 'workloads': wl,
           'note': ('real query values are D1+/D2 query positions; database = first n_D database rows. SYNTHETIC_tiled '
                    'rows: bq2_lib.tiled_rows (database rows 0…n_D−1, query rows QUERY_OFFSET + q), L = 1000; '
                    'performance-only.'),
           'pools_after_exclusion': {d: len(v) for d, v in pools.items()}}
    (HERE / 'results' / 'bq2_workloads.json').write_text(json.dumps(out, indent=1) + '\n')
    for w in wl:
        print(w['id'], w['kind'], w['d'], w['nd'], w['nwin'], w['role'], w['queries'])


if __name__ == '__main__':
    main()
