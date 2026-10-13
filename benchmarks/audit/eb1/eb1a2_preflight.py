#!/usr/bin/env python3
"""EB-1 Amendment 2 read-only workload preflight (NOT a performance measurement).

Derives frame counts from the actual cached arrays (before/after burn-in), validates every requested query, database
and window index of the Amendment 1 design, compares requested vs actual workload sizes, and classifies which
Amendment 1 memory results are affected by (i) retained backing allocations / load transients and (ii) repeated
outputs alive across repetitions. Writes results/eb1a2_preflight.json. Reads array headers and index metadata only.
Run: ../wpd_local_retrieval/build/venv/bin/python -B eb1a2_preflight.py
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import eb1a1_lib as A1  # noqa: E402
import eb1a1_run as R1  # noqa: E402

C0, C1 = HERE / 'build' / 'cache', HERE / 'build' / 'cache_a1'
OUT = HERE / 'results' / 'eb1a2_preflight.json'
COEF = json.loads((HERE / 'results' / 'eb1a1_temp_calibration.json').read_text())['coef']


def shape(p):
    with open(p, 'rb') as f:
        v = np.lib.format.read_magic(f)
        rd = np.lib.format.read_array_header_1_0 if v == (1, 0) else np.lib.format.read_array_header_2_0
        sh, fo, dt = rd(f)
    return list(sh), str(dt)


def main():
    inv = {'arrays': {}}
    for p in sorted(list(C0.glob('*.npy')) + list(C1.glob('*.npy'))):
        inv['arrays'][f'{p.parent.name}/{p.name}'] = dict(zip(('shape', 'dtype'), shape(p)))
    d2_stored = [inv['arrays'][f'cache/D2_stored_{i}.npy']['shape'][0] for i in range(4)]
    d2_post = [n - int(math.floor(0.10 * n)) for n in d2_stored]
    inv['frames'] = {
        'D1_stored_before_burn': inv['arrays']['cache/D1_stored.npy']['shape'][0], 'D1_burn_in': 0,
        'D1_round_all': inv['arrays']['cache/D1_round_all.npy']['shape'][0],
        'D2_stored_per_traj_before_burn': d2_stored, 'D2_stored_total_before_burn': sum(d2_stored),
        'D2_post_burn_per_traj': d2_post, 'D2_post_burn_total': sum(d2_post),
        'D2_round_all_actual': inv['arrays']['cache/D2_round_all.npy']['shape'][0],
        'D1p_round': inv['arrays']['cache_a1/D1p_round.npy']['shape'][0]}
    assert inv['frames']['D2_round_all_actual'] == inv['frames']['D2_post_burn_total']
    metas = {d: json.loads(((C1 if d == 'D1p' else C0) / f'{d}_meta.json').read_text()) for d in ('D1', 'D2', 'D1p')}
    inv['index_sets'] = {d: {'n_q': len(m['qi']), 'n_d': len(m['di']), 'q_range': [min(m['qi']), max(m['qi'])],
                             'd_range': [min(m['di']), max(m['di'])],
                             'q_contiguous': m['qi'] == list(range(m['qi'][0], m['qi'][0] + len(m['qi']))),
                             'd_contiguous': m['di'] == list(range(m['di'][0], m['di'][0] + len(m['di'])))}
                         for d, m in metas.items()}
    # validate the Amendment 1 design
    rows = []
    for c in R1.design():
        a = c['args']
        r = {'block': c['block'], 'args': a}
        if c['block'] in ('S2', 'S3'):
            avail = {('S2', 'D1'): inv['frames']['D1_stored_before_burn'], ('S3', 'D1'): inv['frames']['D1_round_all'],
                     ('S2', 'D2'): inv['frames']['D2_stored_total_before_burn'],
                     ('S3', 'D2'): inv['frames']['D2_round_all_actual']}[(c['block'], a['d'])]
            r.update({'requested_frames': a['n'], 'available_frames': avail, 'valid': a['n'] <= avail})
        else:
            m = metas[a['d']]
            L = inv['arrays'][f"{'cache_a1' if a['d'] == 'D1p' else 'cache'}/{a['d']}_round.npy"]['shape'][1]
            wins = R1.windows_for(L, a['nwin'])
            ok = a['nq'] <= len(m['qi']) and a['nd'] <= len(m['di']) and all(0 <= s <= L - 8 for s in wins)
            r.update({'requested': [a['nq'], a['nd'], a['nwin']], 'available': [len(m['qi']), len(m['di']), L - 8 + 1],
                      'actual_n_win': len(wins), 'n_win_matches': len(wins) == min(a['nwin'], L - 8 + 1), 'valid': ok})
        rows.append(r)
    inv['amendment1_design_validation'] = rows
    inv['amendment1_invalid'] = [r for r in rows if not r['valid'] or r.get('n_win_matches') is False]
    # classify Amendment 1 memory results
    m1 = json.loads((HERE / 'results' / 'eb1a1_measurements.json').read_text())
    cls = []
    for c in m1['cells']:
        a = c['args']
        if c['block'] == 'S2' and a['d'] == 'D2':
            why = 'retained backing: full concatenate of all D2 trajectories alive during the cell'
        elif c['block'] in ('S2', 'S3'):
            why = 'load transient: memory-mapped source pages coexisted with the copy (lifetime peak)'
        elif c['block'] == 'S4':
            res = c.get('resident_bytes', 0)
            transient = c.get('temp_bytes_pred', 0) + c.get('retained_bytes', 0)
            reasons = []
            if c.get('repeats', 1) > 1:
                reasons.append('repeated outputs: previous repetition outputs alive during the next')
            if res >= transient:
                reasons.append('load transient comparable to compute-phase transients')
            why = '; '.join(reasons) or None
        else:
            why = None
        cls.append({'block': c['block'], 'args': a, 'status': c.get('status'), 'memory_affected': why is not None,
                    'reason': why})
    inv['amendment1_memory_classification'] = cls
    inv['note'] = 'Read-only inventory: array headers and index metadata only; no timing.'
    OUT.write_text(json.dumps(inv, indent=1) + '\n')
    print('frames', inv['frames'])
    print('invalid in A1 design:', [(r['block'], r['args']) for r in inv['amendment1_invalid']])
    print('A1 memory cells affected:', sum(x['memory_affected'] for x in cls), 'of', len(cls))
    print('S4 unaffected:', [(x['args']['m'], x['args']['mode'], x['args']['role'], x['args']['nq'], x['args']['nd'],
                              x['args']['nwin']) for x in cls if x['block'] == 'S4' and not x['memory_affected']])


if __name__ == '__main__':
    main()
