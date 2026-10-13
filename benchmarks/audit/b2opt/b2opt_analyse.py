#!/usr/bin/env python3
"""BQ-1 analysis (protocol §6; Amendment 1). Written and hashed before any measurement.

Direct measurements (per workload) are kept separate from fitted P2 projections. Models and support rules reuse the
frozen EB-1 code (eb1_analyse.fit; eb1a1_analyse.s4_model / predict). Writes results/b2opt_analysis.json.
Run: ../wpd_local_retrieval/build/venv/bin/python -B b2opt_analyse.py
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EB1 = HERE.parent / 'eb1'
sys.path.insert(0, str(EB1))
sys.path.insert(0, str(HERE))
import eb1_analyse as AN  # noqa: E402
import eb1a1_analyse as AA  # noqa: E402

RES = HERE / 'results'
IN = RES / 'b2opt_measurements.json'
OUT = RES / 'b2opt_analysis.json'
WL = json.loads((RES / 'b2opt_workloads.json').read_text())
A2 = json.loads((EB1 / 'results' / 'eb1a2_analysis.json').read_text())
P2, LIM, CONF_MAX = AN.P2, AN.LIMITS, 0.25
K, W = 10, 8
ARMS = ('V1', 'V2', 'V2b')
FRESH = [w['id'] for w in WL['workloads'] if w['role'] == 'confirmation_fresh']
FIT = [w['id'] for w in WL['workloads'] if w['role'] == 'fit_reused']
REUSED_CMP = [w['id'] for w in WL['workloads'] if w['role'] == 'comparison_reused']


def med(x):
    return float(np.median(x))


def cell_stats(cells, arm, wid):
    cs = [c for c in cells if c['arm'] == arm and c['wid'] == wid and c.get('status') == 'ok']
    if not cs:
        return None
    w = next(x for x in WL['workloads'] if x['id'] == wid)
    q = [med(c['query_wall_s']) for c in cs]
    first = [c['prep_wall_s'] + c['query_wall_s'][0] for c in cs]
    return {'n_rounds': len(cs), 'prepared_query_wall_s': med(q), 'prepared_query_wall_range': [min(q), max(q)],
            'first_query_wall_s': med(first), 'first_query_range': [min(first), max(first)],
            'prep_wall_s': med([c['prep_wall_s'] for c in cs]), 'prep_cpu_s': med([c['prep_cpu_s'] for c in cs]),
            'query_cpu_s': med([med(c['query_cpu_s']) for c in cs]),
            'cpu_over_wall_query': med([sum(c['query_cpu_s']) / sum(c['query_wall_s']) for c in cs]),
            'maxrss_kB': max(c['rss_end_kB'] for c in cs), 'rss_after_load_kB': max(c['rss_after_load_kB'] for c in cs),
            'rss_after_prep_kB': max(c['rss_after_prep_kB'] for c in cs),
            'resident_bytes': cs[0]['resident_bytes'], 'prepared_bytes': cs[0]['prepared_bytes'],
            'retained_bytes': cs[0]['retained_bytes'], 'fallback_rows_per_run': cs[0]['fallback_rows'],
            'comparisons': w['nd'] * cs[0]['n_win'], 'n_win': cs[0]['n_win'], 'n_d': w['nd'],
            'topk_sha256': [c['topk_sha256'] for c in cs], 'seq': [c['seq'] for c in cs]}


def paired_speedups(cells, arm, wid):
    out = {'prepared': [], 'cold': []}
    for r in sorted({c['round'] for c in cells}):
        b = [c for c in cells if c['round'] == r and c['wid'] == wid and c['arm'] == 'V0' and c.get('status') == 'ok']
        a = [c for c in cells if c['round'] == r and c['wid'] == wid and c['arm'] == arm and c.get('status') == 'ok']
        if a and b:
            b0, a0 = b[0], a[0]
            out['prepared'].append(med(b0['query_wall_s']) / med(a0['query_wall_s']))
            out['cold'].append(b0['query_wall_s'][0] / (a0['prep_wall_s'] + a0['query_wall_s'][0]))
    return {k: ({'median': med(v), 'min': min(v), 'max': max(v), 'n': len(v)} if v else None) for k, v in out.items()}


def time_model(stats, key):
    """Query time c0 + c1·C + c2·n_win (eb1a1_analyse.s4_model, single_time) on fit cells; fresh confirmation."""
    fit = [{'comparisons': s['comparisons'], 'n_win': s['n_win'], 'wall_s': s[key]} for s in (stats[w] for w in FIT)]
    m = AA.s4_model(fit, 'single_time', 'wall_s')
    conf = []
    for wid in FRESH + REUSED_CMP:
        s = stats[wid]
        p = AA.predict(m, s)
        conf.append({'workload': wid, 'deciding': wid in FRESH, 'measured': s[key], 'predicted': p,
                     'abs_error': p - s[key], 'rel_error': abs(p - s[key]) / s[key]})
    m['confirmation'] = conf
    m['confirmation_supported'] = all(e['rel_error'] <= CONF_MAX for e in conf if e['deciding'])
    m['supported'] = bool(m['fit_supported'] and m['confirmation_supported'])
    feats = {'comparisons': P2['n_win'] * P2['n_D'], 'n_win': P2['n_win']}
    m['P2'] = {'value': AA.predict(m, feats), 'extrapolation_factor':
               feats['comparisons'] / max(stats[w]['comparisons'] for w in FIT + FRESH + REUSED_CMP)}
    return m


def prep_model(stats, key):
    xs = [stats[w]['comparisons'] for w in FIT]
    f = AN.fit(xs, [stats[w][key] for w in FIT])
    conf = []
    if 'c0' in f:
        for wid in FRESH + REUSED_CMP:
            s = stats[wid]
            p = f['c0'] + f['c1'] * s['comparisons']
            conf.append({'workload': wid, 'deciding': wid in FRESH, 'measured': s[key], 'predicted': p,
                         'abs_error': p - s[key], 'rel_error': abs(p - s[key]) / s[key]})
    f['confirmation'] = conf
    f['confirmation_supported'] = bool(conf) and all(e['rel_error'] <= CONF_MAX for e in conf if e['deciding'])
    f['supported_full'] = bool(f['supported'] and f['confirmation_supported'])
    x2 = P2['n_win'] * P2['n_D']
    f['P2'] = {'value': f['c0'] + f['c1'] * x2 if 'c0' in f else None,
               'extrapolation_factor': x2 / max(stats[w]['comparisons'] for w in FIT + FRESH + REUSED_CMP)}
    return f


def mem_model(stats, alpha, prepared):
    def P(s):
        return s['resident_bytes'] + s['prepared_bytes'] + s['retained_bytes'] + alpha * s['n_d']
    fit = [{'P': P(stats[w]), 'maxrss_kB': stats[w]['maxrss_kB']} for w in FIT]
    m = AA.s4_model(fit, 'mem', None)
    conf = []
    for wid in FRESH + REUSED_CMP:
        s = stats[wid]
        p = AA.predict(m, {'P': P(s)})
        meas = s['maxrss_kB'] * 1024.0
        conf.append({'workload': wid, 'deciding': wid in FRESH, 'measured': meas, 'predicted': p,
                     'abs_error': p - meas, 'rel_error': abs(p - meas) / meas})
    m['confirmation'] = conf
    m['confirmation_supported'] = all(e['rel_error'] <= CONF_MAX for e in conf if e['deciding'])
    m['supported'] = bool(m['fit_supported'] and m['confirmation_supported'])
    p2_prep = P2['n_D'] * P2['n_win'] * (W * 3 + 1) * 8 if prepared else 0
    p2_P = P2['frames_total'] * P2['L'] * 3 * 8 + p2_prep + P2['n_win'] * K * 8 + alpha * P2['n_D']
    m['P2'] = {'value_bytes': AA.predict(m, {'P': p2_P}), 'P2_predictor_bytes': p2_P, 'P2_prepared_bytes': p2_prep,
               'extrapolation_factor': p2_P / max(P(stats[w]) for w in FIT + FRESH + REUSED_CMP)}
    return m


def main():
    m = json.loads(IN.read_text())
    gate0 = m.get('gate_V0') or {}
    out = {'_meta': {'producer': 'b2opt_analyse.py', 'measurements_sha256': hashlib.sha256(IN.read_bytes()).hexdigest(),
                     'P2': P2, 'limits': LIM, 'conf_max_rel_error': CONF_MAX, 'cpu_h_total': m['_meta'].get('cpu_h_total'),
                     'stopped': m.get('stopped')},
           'tests': m.get('tests'), 'refcheck_pass': bool((m.get('refcheck') or {}).get('pass')),
           'gate_V0_pass': bool(gate0.get('pass')), 'arms': {}}
    r1_b2 = A2['requirements']['R1_B2']['value']
    out['amendment2_R2_B2_s'] = A2['requirements']['R2_B2']['value']
    for arm in ARMS:
        A = m['arms'].get(arm)
        res = {'status': 'untested' if A is None else A['status']}
        out['arms'][arm] = res
        if A is None or 'gate' not in A:
            continue
        g = A['gate']
        res['gate'] = {'pass': g.get('pass'), 'synthetic': g.get('synthetic'),
                       'workloads': {k: {kk: vv for kk, vv in v.items() if kk != 'topk_sha256'}
                                     for k, v in (g.get('workloads') or {}).items()}}
        if A['status'] != 'timed':
            continue
        cells = A['cells']
        stats = {w['id']: cell_stats(cells, arm, w['id']) for w in WL['workloads']}
        base = {w['id']: cell_stats(cells, 'V0', w['id']) for w in WL['workloads']}
        if any(v is None for v in list(stats.values()) + list(base.values())):
            res['status'] = 'timing_incomplete'
            continue
        # identity of timed outputs with the gate outputs (both the arm and the same-session V0)
        ident = all(h == g['workloads'][wid]['topk_sha256'] for wid, s in stats.items() for h in s['topk_sha256'])
        ident0 = all(h == gate0['workloads'][wid]['topk_sha256'] for wid, s in base.items() for h in s['topk_sha256'])
        res['timed_outputs_identical_to_gate'] = bool(ident and ident0)
        res['direct'] = {wid: {'arm': {k: v for k, v in stats[wid].items() if k != 'topk_sha256'},
                               'V0': {k: v for k, v in base[wid].items() if k != 'topk_sha256'},
                               'speedup_vs_V0': paired_speedups(cells, arm, wid)} for wid in stats}
        qm = time_model(stats, 'prepared_query_wall_s')
        q0 = time_model(base, 'prepared_query_wall_s')
        pw = prep_model(stats, 'prep_wall_s')
        pc = prep_model(stats, 'prep_cpu_s')
        mm = mem_model(stats, g['alpha']['alpha_bytes_per_row'], True)
        mm0 = mem_model(base, gate0['alpha']['alpha_bytes_per_row'], False)
        res['models'] = {'query_time': qm, 'prep_wall': pw, 'prep_cpu': pc, 'memory': mm,
                         'V0_same_session_query_time': q0, 'V0_same_session_memory': mm0}
        q2 = qm['P2']['value'] if qm['supported'] else None
        p2w = pw['P2']['value'] if pw['supported_full'] else None
        p2c = pc['P2']['value'] if pc['supported_full'] else None
        mem2 = mm['P2']['value_bytes'] if mm['supported'] else None
        r1 = None if p2c is None else r1_b2 + p2c / 3600
        crit = {'a_gate_and_identity': bool(g.get('pass') and res['timed_outputs_identical_to_gate']),
                'b_query_supported_and_le_10s': q2 is not None and q2 <= LIM['R2_wall_s'],
                'c_prep_models_supported': p2w is not None and p2c is not None,
                'd_memory_supported_and_le_8GB': mem2 is not None and mem2 / 1024 <= LIM['R4_rss_kB'],
                'e_R1_B2_plus_prep_le_96h': r1 is not None and r1 <= LIM['R1_cpu_h']}
        first2 = None if q2 is None or p2w is None else q2 + p2w
        res['projections_P2'] = {'prepared_query_s': q2, 'prep_wall_s': p2w, 'prep_cpu_s': p2c,
                                 'first_query_cold_s': first2, 'peak_rss_bytes': mem2,
                                 'R1_B2_plus_prep_cpu_h': r1, 'V0_same_session_query_s':
                                 q0['P2']['value'] if q0['supported'] else None,
                                 'amortized_s': ({str(n): p2w / n + q2 for n in (1, 10, 100, 1000)}
                                                 if first2 is not None else None),
                                 'break_even_n_for_10s': (p2w / (LIM['R2_wall_s'] - q2) if first2 is not None
                                                          and q2 < LIM['R2_wall_s'] else None)}
        res['criteria'] = crit
        res['R2_supported_prepared_database'] = all(crit.values())
        res['R2_supported_cold'] = bool(all(crit.values()) and first2 is not None and first2 <= LIM['R2_wall_s'])
    qual = [a for a, v in out['arms'].items() if v.get('R2_supported_prepared_database')]
    out['outcome'] = (f'R2_B2 supported (prepared database) by: {qual}' if qual else
                      'R2_B2 not supported by the tested variants on these workloads')
    out['qualifying_cold'] = [a for a, v in out['arms'].items() if v.get('R2_supported_cold')]
    OUT.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(out['outcome'], '| cold:', out['qualifying_cold'], '|', {a: v['status'] for a, v in out['arms'].items()})
    return 0


if __name__ == '__main__':
    sys.exit(main())
