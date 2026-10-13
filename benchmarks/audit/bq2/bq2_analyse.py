#!/usr/bin/env python3
"""BQ-2 analysis (protocol §8). Written and hashed before any BQ-2 measurement.

Direct measurements are kept separate from fitted projections. Models and support rules reuse frozen EB-1 code
(eb1a1_analyse.s4_model / predict). Writes results/bq2_analysis.json.
Run: ../wpd_local_retrieval/build/venv/bin/python -B bq2_analyse.py
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EB1 = HERE.parent / 'eb1'
B1DIR = HERE.parent / 'b2opt'
sys.path.insert(0, str(EB1))
import eb1_analyse as AN  # noqa: E402
import eb1a1_analyse as AA  # noqa: E402

RES = HERE / 'results'
IN = RES / 'bq2_measurements.json'
OUT = RES / 'bq2_analysis.json'
WL = {w['id']: w for w in json.loads((RES / 'bq2_workloads.json').read_text())['workloads']}
A2 = json.loads((EB1 / 'results' / 'eb1a2_analysis.json').read_text())
B1A = json.loads((B1DIR / 'results' / 'b2opt_analysis.json').read_text())
B1M = json.loads((B1DIR / 'results' / 'b2opt_measurements.json').read_text())
ALPHA = B1M['arms']['V2b']['gate']['alpha']['alpha_bytes_per_row']
P2, LIM, CONF_MAX = AN.P2, AN.LIMITS, 0.25
K, W = 10, 8
FIT = [k for k, w in WL.items() if w['role'] == 'fit']
REQUIRED = ['RC1', 'RC2', 'RC3', 'S1', 'S2', 'S4_P2slice']
OPTIONAL = ['S3']


def med(x):
    return float(np.median(x))


def stats(cells, arm, wid):
    cs = [c for c in cells if c['arm'] == arm and c['wid'] == wid and c.get('status') == 'ok']
    if not cs:
        return None
    q = [med(c['query_wall_s']) for c in cs]
    pw = [c['prep_wall_s'] for c in cs]
    allreps = [x for c in cs for x in c.get('prep_wall_reps_s', [])]
    return {'n_rounds': len(cs), 'prepared_query_wall_s': med(q), 'prepared_query_range': [min(q), max(q)],
            'query1_wall_s': med([c['query_wall_s'][0] for c in cs]),
            'prep_wall_s': med(pw), 'prep_cpu_s': med([c['prep_cpu_s'] for c in cs]),
            'prep_rep_quartiles_s': ([float(x) for x in np.percentile(allreps, [25, 50, 75])] if allreps else None),
            'prep_rep_rel_iqr': (float((np.percentile(allreps, 75) - np.percentile(allreps, 25)) / np.median(allreps))
                                 if allreps else None),
            'first_query_wall_s': med(pw) + med([c['query_wall_s'][0] for c in cs]),
            'cpu_over_wall_query': med([sum(c['query_cpu_s']) / sum(c['query_wall_s']) for c in cs]),
            'maxrss_kB': max(c['rss_end_kB'] for c in cs), 'rss_after_load_kB': max(c['rss_after_load_kB'] for c in cs),
            'rss_after_prep_kB': max(c['rss_after_prep_kB'] for c in cs), 'resident_bytes': cs[0]['resident_bytes'],
            'prepared_bytes': cs[0]['prepared_bytes'], 'retained_bytes': cs[0]['retained_bytes'],
            'fallback_rows_per_run': max(c['fallback_rows'] for c in cs), 'n_d': cs[0]['n_d'], 'n_win': cs[0]['n_win'],
            'comparisons': cs[0]['n_d'] * cs[0]['n_win'], 'topk_sha256': [c['topk_sha256'] for c in cs],
            'loadavg_1min': [float(c['loadavg_start'][0]) for c in cs], 'seq': [c['seq'] for c in cs]}


def paired(cells, wid):
    out = {'prepared': [], 'cold': []}
    for r in sorted({c['round'] for c in cells}):
        a = [c for c in cells if c['round'] == r and c['wid'] == wid and c['arm'] == 'V2b' and c.get('status') == 'ok']
        b = [c for c in cells if c['round'] == r and c['wid'] == wid and c['arm'] == 'V0' and c.get('status') == 'ok']
        if a and b:
            out['prepared'].append(med(b[0]['query_wall_s']) / med(a[0]['query_wall_s']))
            out['cold'].append(b[0]['query_wall_s'][0] / (a[0]['prep_wall_s'] + a[0]['query_wall_s'][0]))
    return {k: ({'median': med(v), 'min': min(v), 'max': max(v), 'n': len(v)} if v else None) for k, v in out.items()}


def mem_P(s):
    return s['resident_bytes'] + s['prepared_bytes'] + s['retained_bytes'] + ALPHA * s['n_d']


def model(S, key, kind, conf_ids, diag=None):
    if kind == 'mem':
        rows = [{'P': mem_P(S[w]), 'maxrss_kB': S[w]['maxrss_kB']} for w in FIT]
        m = AA.s4_model(rows, 'mem', None)
    else:
        rows = [{'comparisons': S[w]['comparisons'], 'n_win': S[w]['n_win'], 'y': S[w][key]} for w in FIT]
        m = AA.s4_model(rows, 'single_time', 'y')
    conf = []
    for wid in conf_ids:
        s = S[wid]
        if kind == 'mem':
            p, meas = AA.predict(m, {'P': mem_P(s)}), s['maxrss_kB'] * 1024.0
        else:
            p, meas = AA.predict(m, s), s[key]
        conf.append({'workload': wid, 'kind': WL[wid]['kind'], 'measured': meas, 'predicted': p,
                     'abs_error': p - meas, 'rel_error': abs(p - meas) / meas})
    m['confirmation'] = conf
    m['confirmation_supported'] = bool(conf) and all(e['rel_error'] <= CONF_MAX for e in conf)
    m['supported'] = bool(m['fit_supported'] and m['confirmation_supported'])
    m['diagnostic_BQ1'] = diag(m) if diag else None
    largest = max([S[w]['comparisons'] for w in FIT + conf_ids])
    if kind == 'mem':
        p2_P = (P2['frames_total'] * P2['L'] * 3 * 8 + P2['n_D'] * P2['n_win'] * (W * 3 + 1) * 8
                + P2['n_win'] * K * 8 + ALPHA * P2['n_D'])
        m['P2'] = {'value': AA.predict(m, {'P': p2_P}), 'P2_predictor_bytes': p2_P,
                   'extrapolation_factor': p2_P / max(mem_P(S[w]) for w in FIT + conf_ids)}
    else:
        f2 = {'comparisons': P2['n_D'] * P2['n_win'], 'n_win': P2['n_win']}
        m['P2'] = {'value': AA.predict(m, f2), 'extrapolation_factor': f2['comparisons'] / largest}
    return m


def bq1_diag(key):
    def f(m):
        rows = []
        for wid, d in B1A['arms']['V2b']['direct'].items():
            a = d['arm']
            feats = {'comparisons': a['comparisons'], 'n_win': a['n_win']}
            p = AA.predict(m, feats)
            meas = a[key]
            rows.append({'bq1_workload': wid, 'measured': meas, 'predicted': p, 'abs_error': p - meas,
                         'rel_error': abs(p - meas) / meas})
        return rows
    return f


def main():
    M = json.loads(IN.read_text())
    cells = M.get('cells', [])
    out = {'_meta': {'producer': 'bq2_analyse.py', 'measurements_sha256': hashlib.sha256(IN.read_bytes()).hexdigest(),
                     'P2': P2, 'limits': LIM, 'conf_max_rel_error': CONF_MAX, 'alpha_bytes_per_row_from_BQ1': ALPHA,
                     'cpu_h_total': M['_meta'].get('cpu_h_total'), 'stopped': M.get('stopped'),
                     'loadavg_start_end': [M['_meta'].get('loadavg_start'), M['_meta'].get('loadavg_end')]},
           'tests': M.get('tests')}
    gate_ok = {arm: {wid: bool((g or {}).get('status') == 'ok' and g.get('pass')) for wid, g in M['gate'][arm].items()}
               for arm in ('V2b', 'V0')}
    out['gates'] = {arm: {wid: {k: v for k, v in (g or {}).items() if k != 'topk_sha256'}
                          for wid, g in M['gate'][arm].items()} for arm in ('V2b', 'V0')}
    measured_correctness_failure = [f'{arm} {wid}' for arm in ('V2b', 'V0') for wid, g in M['gate'][arm].items()
                                    if (g or {}).get('status') == 'ok' and not g.get('pass')]
    S = {wid: stats(cells, 'V2b', wid) for wid in WL}
    S0 = {wid: stats(cells, 'V0', wid) for wid in WL if WL[wid]['baseline']}
    missing = [w for w in FIT + REQUIRED if S.get(w) is None or S[w]['n_rounds'] < 3]
    missing += [f'V0 {w}' for w in REQUIRED if S0.get(w) is None or S0[w]['n_rounds'] < 3]
    conf_ids = [w for w in REQUIRED + OPTIONAL if S.get(w) is not None and S[w]['n_rounds'] == 3]
    out['optional_measured'] = {w: S.get(w) is not None and S[w]['n_rounds'] == 3 for w in OPTIONAL}
    out['missing_required'] = missing
    # identity of timed outputs with gate outputs
    ident = True
    for arm, SS in (('V2b', S), ('V0', S0)):
        for wid, s in SS.items():
            if s is None or wid == 'P2D':
                continue
            g = M['gate'][arm].get(wid) or {}
            ident &= all(h == g.get('topk_sha256') for h in s['topk_sha256'])
    out['timed_outputs_identical_to_gate'] = bool(ident)
    out['direct'] = {wid: {'V2b': {k: v for k, v in S[wid].items() if k != 'topk_sha256'} if S.get(wid) else None,
                           'V0': ({k: v for k, v in S0[wid].items() if k != 'topk_sha256'} if S0.get(wid) else None),
                           'speedup_vs_V0': paired(cells, wid) if WL[wid]['baseline'] else None,
                           'kind': WL[wid]['kind'], 'role': WL[wid]['role']} for wid in WL}
    p2d = [c for c in cells if c['wid'] == 'P2D']
    out['P2D'] = ({k: v for k, v in p2d[0].items() if k != 'topk_sha256'} if p2d else {'status': 'not_attempted'})
    if p2d and p2d[0].get('status') == 'ok':
        g = M['gate']['V2b'].get('P2D') or {}
        out['P2D']['gate'] = {k: v for k, v in g.items() if k != 'topk_sha256'}
        out['P2D']['identical_to_gate'] = p2d[0]['topk_sha256'] == g.get('topk_sha256')
    if missing or any(S[w] is None for w in FIT):
        out['models'] = None
    else:
        out['models'] = {
            'query_time': model(S, 'prepared_query_wall_s', 'time', conf_ids, bq1_diag('prepared_query_wall_s')),
            'prep_wall': model(S, 'prep_wall_s', 'time', conf_ids, bq1_diag('prep_wall_s')),
            'prep_cpu': model(S, 'prep_cpu_s', 'time', conf_ids, bq1_diag('prep_cpu_s')),
            'memory': model(S, None, 'mem', conf_ids),
            'V0_query_time_on_confirmations': {wid: S0[wid]['prepared_query_wall_s'] for wid in S0 if S0[wid]}}
    # projections and conclusion
    r1_b2 = A2['requirements']['R1_B2']['value']
    proj, concl = {}, {}
    mods = out['models']
    if mods:
        sup = {k: mods[k]['supported'] for k in ('query_time', 'prep_wall', 'prep_cpu', 'memory')}
        q2 = mods['query_time']['P2']['value']
        pw2, pc2 = mods['prep_wall']['P2']['value'], mods['prep_cpu']['P2']['value']
        m2 = mods['memory']['P2']['value']
        sl = S['S4_P2slice']['prepared_query_wall_s'] * P2['n_win'] / S['S4_P2slice']['n_win']
        proj = {'prepared_query_s': q2, 'prep_wall_s': pw2, 'prep_cpu_s': pc2, 'first_query_s': q2 + pw2,
                'peak_rss_bytes': m2, 'R1_B2_plus_prep_cpu_h': r1_b2 + pc2 / 3600,
                'slice_scaled_prepared_query_s': sl,
                'amortized_s': {str(n): pw2 / n + q2 for n in (1, 10, 100, 1000)},
                'break_even_n_for_10s': (pw2 / (LIM['R2_wall_s'] - q2) if q2 < LIM['R2_wall_s'] else None),
                'models_supported': sup, 'values_valid_only_if_supported': True}
        direct_q = (out['P2D'].get('query_wall_s') or [None])[0] if out['P2D'].get('status') == 'ok' else None
        evidence_ok = (all(sup.values()) and not missing and out['timed_outputs_identical_to_gate']
                       and not measured_correctness_failure)
        within = (q2 <= LIM['R2_wall_s'] and m2 / 1024 <= LIM['R4_rss_kB']
                  and r1_b2 + pc2 / 3600 <= LIM['R1_cpu_h'] and (direct_q is None or direct_q <= LIM['R2_wall_s']))
        concl['A_prepared_database_supported'] = bool(evidence_ok and within)
        concl['B_first_query_warrants_optimization'] = bool(concl['A_prepared_database_supported']
                                                            and q2 + pw2 > LIM['R2_wall_s'])
        concl['C_insufficient_evidence'] = not evidence_ok
        concl['projected_violation_with_supported_models'] = bool(evidence_ok and not within)
    else:
        concl = {'A_prepared_database_supported': False, 'B_first_query_warrants_optimization': False,
                 'C_insufficient_evidence': True, 'projected_violation_with_supported_models': False}
    concl['measured_correctness_failure'] = measured_correctness_failure
    measured_violations = []
    for wid, d in out['direct'].items():
        s = d['V2b']
        if s and s['maxrss_kB'] > LIM['R4_rss_kB']:
            measured_violations.append(f'{wid} peak RSS')
    if out['P2D'].get('status') == 'ok' and out['P2D']['query_wall_s'][0] > LIM['R2_wall_s']:
        measured_violations.append('P2D prepared-query latency > 10 s (direct)')
    concl['measured_resource_limit_violation'] = measured_violations
    why = []
    if missing:
        why.append(f'missing required cells: {missing}')
    if mods:
        why += [f"{k} model unsupported (fit residual {mods[k]['max_fit_rel_residual']:.3g}; max confirmation error "
                f"{max((e['rel_error'] for e in mods[k]['confirmation']), default=float('nan')):.3g})"
                for k in ('query_time', 'prep_wall', 'prep_cpu', 'memory') if not mods[k]['supported']]
    if not out['timed_outputs_identical_to_gate']:
        why.append('timed outputs differ from gate outputs')
    concl['missing_evidence'] = why
    out['projections_P2'] = proj
    out['conclusion'] = concl
    OUT.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps({k: v for k, v in concl.items()}, default=str))
    return 0


if __name__ == '__main__':
    sys.exit(main())
