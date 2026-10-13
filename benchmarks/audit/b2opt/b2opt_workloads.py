#!/usr/bin/env python3
"""B2 single-query optimization experiment: workload definition (read-only; no timing, no RMSD computation).

Reconstructs every query index previously used on D1, D1+ and D2 by EB-1, Amendment 1 and Amendment 2 (from their
code and seeds), marks the earlier single-query grid as REUSED comparisons, and draws fresh confirmation queries from
the remaining indices with a fixed seed. Writes results/b2opt_workloads.json.
Run: ../wpd_local_retrieval/build/venv/bin/python -B b2opt_workloads.py
"""
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EB1 = HERE.parent / 'eb1'
SEED = 20261012                 # EB-1 seed (eb1_run.SEED)
FRESH_SEED = 20261013           # this experiment
W = 8


def meta(d):
    src = EB1 / 'build' / ('cache_a1' if d == 'D1p' else 'cache')
    m = json.loads((src / f'{d}_meta.json').read_text())
    L = int(np.load(src / f'{d}_round.npy', mmap_mode='r').shape[1])
    return m, L


def main():
    used = {'D1': set(), 'D2': set(), 'D1p': set()}
    src = {}
    for d in ('D1', 'D2'):
        m, L = meta(d)
        qa = np.array(m['qi'])
        for Wx in (5, 8):                                       # EB-1 M4 single-query: rng(SEED + W), 5 queries
            s = np.random.default_rng(SEED + Wx).choice(qa, 5, replace=False).tolist()
            used[d] |= set(s)
        src[f'{d}_EB1_M4'] = 'rng(SEED+W).choice(qi, 5), W in {5, 8}'
    rng = np.random.default_rng(SEED)                           # EB-1 C3: windows then queries, D1 W5, D1 W8, D2 W5, D2 W8
    for d in ('D1', 'D2'):
        m, L = meta(d)
        qa = np.array(m['qi'])
        for Wx in (5, 8):
            rng.choice(L - Wx + 1, 5, replace=False)
            used[d] |= set(rng.choice(qa, 20, replace=False).tolist())
    m, L = meta('D1p')                                          # Amendment 1 C3': rng(SEED+1), 20 queries
    used['D1p'] |= set(np.random.default_rng(SEED + 1).choice(np.array(m['qi']), 20, replace=False).tolist())
    # batch prefixes (EB-1 / A1 / A2: largest n_Q per dataset), C4' (D2 first 256), single-query position 0
    prefix = {'D1': 512, 'D1p': 384, 'D2': 256}
    for d, n in prefix.items():
        mm, _ = meta(d)
        used[d] |= set(mm['qi'][:n])
    eligible = {}
    for d in used:
        mm, L = meta(d)
        eligible[d] = sorted(set(mm['qi']) - used[d])
    # D1 query frames are D1 positions; a frame used in D1 is the same structure as in D1+ (D1+ queries = D1 queries)
    shared = sorted(set(eligible['D1']) & set(eligible['D1p']))
    rng = np.random.default_rng(FRESH_SEED)
    reused = []
    for nd in (250, 500, 1015):
        for nw in (20, 207):
            reused.append({'id': f'G_nD{nd}_nwin{nw}', 'd': 'D1', 'nd': nd, 'nwin': nw, 'role': 'fit_reused'})
    reused += [{'id': 'CW3', 'd': 'D1p', 'nd': 2693, 'nwin': 207, 'role': 'comparison_reused'},
               {'id': 'CW5', 'd': 'D2', 'nd': 405, 'nwin': 44, 'role': 'comparison_reused'}]
    for c in reused:
        c['queries'] = [0, 1, 2]                                # query 0 was used before; 1, 2 are first-batch frames
    fresh_cfg = [('FQ1', 'D1', 750, 120), ('FQ2', 'D1p', 2000, 207), ('FQ3', 'D2', 405, 44), ('FQ4', 'D1p', 2693, 150)]
    fresh = []
    taken = set()
    for cid, d, nd, nw in fresh_cfg:
        pool = [q for q in (shared if d in ('D1', 'D1p') else eligible[d]) if q not in taken]
        qs = sorted(rng.choice(pool, 3, replace=False).tolist())
        taken |= set(qs)
        fresh.append({'id': cid, 'd': d, 'nd': nd, 'nwin': nw, 'role': 'confirmation_fresh', 'queries': qs})
    out = {'W': W, 'k': 10, 'previously_used_queries': {d: sorted(v) for d, v in used.items()},
           'previous_sources': {**src, 'EB1_C3': 'rng(SEED): D1 W5, D1 W8, D2 W5, D2 W8 (5 windows, 20 queries each)',
                                'A1_C3p': 'rng(SEED+1).choice(D1p qi, 20)',
                                'batch_and_C4p_prefixes': prefix, 'single_query_position': 0},
           'eligible_counts': {d: len(v) for d, v in eligible.items()}, 'workloads': reused + fresh,
           'note': 'Query values are positions in the dataset query list; database = first n_D database rows.'}
    (HERE / 'results' / 'b2opt_workloads.json').write_text(json.dumps(out, indent=1) + '\n')
    for c in reused + fresh:
        print(c['id'], c['d'], c['nd'], c['nwin'], c['role'], c['queries'])


if __name__ == '__main__':
    main()
