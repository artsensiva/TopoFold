#!/usr/bin/env python3
"""EB-1 Amendment 2 analysis (docs/EB1_AMENDMENT_2.md). Written and hashed before any Amendment 2 measurement.

- Time: Amendment 1 times are reused for every cell; the replaced S3 D2 size (1600) and FC1–FC4 use Amendment 2 times.
- Memory: re-measured Amendment 2 values replace every affected Amendment 1 memory result; the 10 unaffected S4 cells
  keep their Amendment 1 values; S1 is reused from the original run with the Amendment 1 predictor.
- A batch time uses the frozen replacement model t = c1·C + c2·call_db_rows + c3·calls (no intercept; negative terms
  dropped in the order c3, c2), decided by FC1–FC4; Amendment 1 CW1, CW2, CW4 errors are reported, not deciding.
- All other S4 models, support rules, thresholds, requirement composition and the decision (eb1_analyse.decide) are
  as in Amendment 1. Validity is reported separately for the loader, A and B2.
Writes results/eb1a2_analysis.json. Run: ../wpd_local_retrieval/build/venv/bin/python -B eb1a2_analyse.py
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import eb1_analyse as AN  # noqa: E402
import eb1a1_analyse as AA  # noqa: E402
import eb1a1_lib as A1  # noqa: E402

RES = HERE / 'results'
IN2 = RES / 'eb1a2_measurements.json'
IN1 = RES / 'eb1a1_measurements.json'
ORIG = RES / 'eb1_measurements.json'
OUT = RES / 'eb1a2_analysis.json'
COEF = json.loads((RES / 'eb1a1_temp_calibration.json').read_text())['coef']
P2, LIMITS, RESID_MAX, CONF_MAX = AN.P2, AN.LIMITS, 0.25, 0.25
W = 8
FC = ('FC1', 'FC2', 'FC3', 'FC4')
GROUPS = {'loader': ('S1', 'S2', 'R3', 'R4_S1', 'R4_S2'),
          'A': ('S3', 'S4_A', 'R1_A', 'R2_A', 'R4_S3', 'R4_S4_A'),
          'B2': ('S4_B2', 'R1_B2', 'R2_B2', 'R4_S4_B2')}


def group_of(name):
    for g, pre in GROUPS.items():
        if name.startswith(pre):
            return g
    raise ValueError(name)


def cell_group(c):
    return {'S2': 'loader', 'S3': 'A'}.get(c['block']) or c['args']['m']


def a_batch_time_model(fit_rows):
    """Frozen replacement: no intercept; drop c3 then c2 while negative. Same support rule as Amendment 1."""
    names = ['c1', 'c2', 'c3']
    keys = {'c1': 'comparisons', 'c2': 'call_db_rows', 'c3': 'calls'}
    y = np.array([c['cpu_s'] for c in fit_rows], float)
    n_distinct = len({(c['comparisons'], c['call_db_rows'], c['calls']) for c in fit_rows})
    dropped = []
    while True:
        X = np.array([[c[keys[n]] for n in names] for c in fit_rows], float)
        coef = np.linalg.lstsq(X, y, rcond=None)[0]
        neg = [n for n in ('c3', 'c2') if n in names and coef[names.index(n)] < 0]
        if not neg:
            break
        names.remove(neg[0])
        dropped.append(neg[0])
    rel = np.abs(y - X @ coef) / np.abs(y)
    cd = dict(zip(names, map(float, coef)))
    resid = [{'nq': c['args']['nq'], 'nd': c['args']['nd'], 'n_win': c['n_win'], 'measured': float(m),
              'predicted': float(p), 'abs_error': float(p - m), 'rel_error': float(r)}
             for c, m, p, r in zip(fit_rows, y, X @ coef, rel)]
    return {'kind': 'batch_time_no_intercept', 'coef': cd, 'dropped': dropped, 'n_fit': len(fit_rows),
            'n_distinct': n_distinct, 'max_fit_rel_residual': float(rel.max()), 'fit_residuals': resid,
            'fit_supported': bool(n_distinct >= 3 and cd['c1'] > 0 and rel.max() <= RESID_MAX)}


def predict_nointercept(model, f):
    keys = {'c1': 'comparisons', 'c2': 'call_db_rows', 'c3': 'calls'}
    return sum(v * f[keys[k]] for k, v in model['coef'].items())


def ckey(c):
    return json.dumps([c['block'], c['args'], c.get('rep', 0)], sort_keys=True)


def main():
    m2 = json.loads(IN2.read_text())
    m1 = json.loads(IN1.read_text())
    orig = json.loads(ORIG.read_text())
    missing = []                                                # (group or 'all', text)
    planned = m2['_meta'].get('n_cells_planned')
    if planned is None or len(m2.get('cells', [])) < planned:
        missing.append(('all', f"Amendment 2 cells never run ({m2.get('stopped')})"))
    missing += [(cell_group(c), f"A2 {c['block']} {c['args']}: {c.get('status')}") for c in m2.get('cells', [])
                if c.get('status') != 'ok']
    a2 = {ckey({**c, 'args': c['amends'] or c['args']}): c for c in m2.get('cells', []) if c.get('status') == 'ok'}
    # assemble the merged cell set: Amendment 1 cells (times) with corrected memory where affected
    merged, mem_source = [], {'A2_remeasured': 0, 'A1_unaffected_reused': 0}
    pre = json.loads((RES / 'eb1a2_preflight.json').read_text())
    affected = {json.dumps([x['block'], x['args']], sort_keys=True) for x in pre['amendment1_memory_classification']
                if x['memory_affected']}
    for c in m1['cells']:
        k = ckey(c)
        aff = json.dumps([c['block'], c['args']], sort_keys=True) in affected
        if not aff:
            if c.get('status') != 'ok':
                missing.append((cell_group(c), f"A1 {c['block']} {c['args']}: {c.get('status')}"))
                continue
            merged.append({**c, 'mem_source': 'A1'})
            mem_source['A1_unaffected_reused'] += 1
            continue
        n = a2.get(k)
        if n is None:
            continue                                            # already listed as missing above
        if c.get('status') == 'ok':                             # A1 time, A2 memory (and A2 P, identical definition)
            merged.append({**c, 'maxrss_kB': n['maxrss_kB'], 'P': n['P'], 'mem_source': 'A2',
                           'time_source': 'A1', 'a2_time_s': {'wall_s': n['wall_s'], 'cpu_s': n['cpu_s']}})
        else:                                                   # replaced invalid cell (S3 D2 1620 → 1600)
            merged.append({**n, 'mem_source': 'A2', 'time_source': 'A2'})
        mem_source['A2_remeasured'] += 1
    fresh = [c for c in m2.get('cells', []) if c.get('status') == 'ok' and c['args'].get('role') in FC]
    for c in fresh:
        merged.append({**c, 'mem_source': 'A2', 'time_source': 'A2'})
    lc = m2.get('loadcheck', {})
    correct = (bool(lc.get('pass')) and bool(m1.get('correctness', {}).get('pass'))
               and bool(orig.get('correctness', {}).get('pass')))
    fits, unsupported, proj = {}, [], {}

    def simple(name, rows, xkey, ykey, x_p2, scale=1.0):
        by = {}
        for r in rows:
            by.setdefault(r[xkey], []).append(r[ykey])
        xs = sorted(by)
        f = AN.fit(xs, [float(np.median(by[x])) * scale for x in xs])
        fits[name] = f
        p = AN.project(f, x_p2)
        proj[name] = p
        if p is None:
            unsupported.append(name)
        return None if p is None else p['value']

    o1 = [c for c in orig['cells'] if c.get('status') == 'ok' and c['block'] == 'S1']
    for c in o1:
        if c['args']['d'] == 'D1':
            c['P'] = c['args']['n'] * 214 * 12
        else:
            fr = [450, 450, 430, 450][:c['args']['n']]
            c['P'] = max(fr) * 51 * 12 * 8 + sum(fr) * 51 * 12 + 64 * 4 * 51 * 12
    s1x = simple('S1_xtc_wall', [c for c in o1 if c['args']['d'] == 'D1'], 'x', 'wall_s',
                 P2['frames_total'] * P2['atoms_xtc'])
    s1h = simple('S1_h5_wall', [c for c in o1 if c['args']['d'] == 'D2'], 'x', 'wall_s',
                 P2['frames_total'] * P2['atoms_h5_read'])
    s1xm = simple('R4_S1_xtc', [c for c in o1 if c['args']['d'] == 'D1'], 'P', 'maxrss_kB',
                  P2['frames_total'] * P2['L'] * 12, 1024.0)
    s1hm = simple('R4_S1_h5', [c for c in o1 if c['args']['d'] == 'D2'], 'P', 'maxrss_kB',
                  10000 * P2['L'] * 12 * 8 + P2['frames_total'] * P2['L'] * 12 + 64 * 4 * P2['L'] * 12, 1024.0)
    S2 = [c for c in merged if c['block'] == 'S2']
    S3 = [c for c in merged if c['block'] == 'S3']
    s2 = [simple(f'S2_{d}_wall', [c for c in S2 if c['args']['d'] == d], 'x', 'wall_s', P2['frames_total'] * P2['L'])
          for d in ('D1', 'D2')]
    s2m = simple('R4_S2', S2, 'P', 'maxrss_kB', P2['frames_total'] * P2['L'] * 72, 1024.0)
    s3cpu = [simple(f'S3_{d}_cpu', [c for c in S3 if c['args']['d'] == d], 'x', 'cpu_s', P2['frames_total'] * P2['L'])
             for d in ('D1', 'D2')]
    s3one = [simple(f'S3_{d}_wall_one_frame', [c for c in S3 if c['args']['d'] == d], 'x', 'wall_s', P2['L'])
             for d in ('D1', 'D2')]
    s3m = simple('R4_S3', S3, 'P', 'maxrss_kB',
                 P2['frames_total'] * P2['L'] * 36 + P2['frames_total'] * (2 * P2['L'] - 5) * 8, 1024.0)
    s4 = {}
    for meth in ('A', 'B2'):
        for mode in ('batch', 'single'):
            rows = [c for c in merged if c['block'] == 'S4' and c['args']['m'] == meth and c['args']['mode'] == mode]
            fit_rows = [c for c in rows if c['args']['role'] == 'fit']
            conf_rows = [c for c in rows if c['args']['role'] not in ('fit',)]
            nq = P2['n_Q'] if mode == 'batch' else 1
            b = A1.effective_batch(nq)
            wk = A1.work(nq, P2['n_D'], P2['n_win'])
            p2_P = (A1.resident_bytes(meth, P2['frames_total'], P2['L']) + A1.retained_bytes(nq, P2['n_win'])
                    + A1.temp_bytes(meth, b, P2['n_D'], W, COEF))
            p2_feats = {**wk, 'n_win': P2['n_win'], 'P': p2_P}
            for kind, tkey in ((f'{mode}_time', 'cpu_s' if mode == 'batch' else 'wall_s'), ('mem', None)):
                name = f"S4_{meth}_{mode}_{'mem' if kind == 'mem' else 'time'}"
                if not fit_rows:
                    unsupported.append(name)
                    continue
                replacement = meth == 'A' and kind == 'batch_time'
                if replacement:
                    mod = a_batch_time_model(fit_rows)
                    mod['amendment1_model_for_reference'] = AA.s4_model(fit_rows, kind, tkey)
                    pred = predict_nointercept
                else:
                    mod = AA.s4_model(fit_rows, kind, tkey)
                    pred = AA.predict
                errs, seen = [], []
                for c in conf_rows:
                    meas = c['maxrss_kB'] * 1024.0 if kind == 'mem' else c[tkey]
                    p = pred(mod, c)
                    e = {'workload': c['args']['role'], 'data': c['args']['d'], 'measured': meas, 'predicted': p,
                         'abs_error': p - meas, 'rel_error': abs(p - meas) / meas}
                    (seen if replacement and c['args']['role'] not in FC else errs).append(e)
                mod['confirmation'] = errs
                if replacement:
                    mod['amendment1_confirmations_not_deciding'] = seen
                mod['confirmation_supported'] = bool(errs) and all(e['rel_error'] <= CONF_MAX for e in errs)
                mod['supported'] = mod['fit_supported'] and mod['confirmation_supported']
                xs = [c['comparisons'] if kind != 'mem' else c['P'] for c in fit_rows + conf_rows]
                val = pred(mod, p2_feats)
                mod['P2'] = {'value': val, 'P2_predictor_bytes': p2_P, 'extrapolation_factor_vs_largest_measured':
                             (p2_feats['comparisons'] if kind != 'mem' else p2_P) / max(xs)}
                fits[name] = mod
                s4[name] = val if mod['supported'] else None
                if not mod['supported']:
                    unsupported.append(name)
    req = {}
    for meth in ('A', 'B2'):
        bt, st_ = s4.get(f'S4_{meth}_batch_time'), s4.get(f'S4_{meth}_single_time')
        if meth == 'A':
            r1 = None if None in s3cpu + [bt] else (max(s3cpu) + bt) / 3600
            r2 = None if None in s3one + [st_] else max(s3one) + st_
        else:
            r1 = None if bt is None else bt / 3600
            r2 = st_
        req[f'R1_{meth}'] = {'value': r1, 'limit': LIMITS['R1_cpu_h'], 'component': 'engine', 'unit': 'CPU-h'}
        req[f'R2_{meth}'] = {'value': r2, 'limit': LIMITS['R2_wall_s'], 'component': 'engine', 'unit': 's'}
        for mode in ('batch', 'single'):
            v = s4.get(f'S4_{meth}_{mode}_mem')
            req[f'R4_S4_{meth}_{mode}'] = {'value': None if v is None else v / 1024, 'limit': LIMITS['R4_rss_kB'],
                                          'component': 'engine', 'unit': 'kB'}
    for lay, v in (('xtc', s1x), ('h5', s1h)):
        req[f'R3_{lay}'] = {'value': None if v is None or None in s2 else v + max(s2), 'limit': LIMITS['R3_wall_s'],
                            'component': 'loader', 'unit': 's'}
    for k, v, comp in (('R4_S1_xtc', s1xm, 'loader'), ('R4_S1_h5', s1hm, 'loader'), ('R4_S2', s2m, 'loader'),
                       ('R4_S3_A', s3m, 'engine')):
        req[k] = {'value': None if v is None else v / 1024, 'limit': LIMITS['R4_rss_kB'], 'component': comp, 'unit': 'kB'}
    unsupported = sorted(set(unsupported))
    decision = AN.decide(json.loads(json.dumps(req)), [t for _, t in missing], correct, unsupported)
    ratios = {k: (v['value'] / v['limit'] if v['value'] is not None else None) for k, v in req.items()}
    # per-group validity (reported beside the global outcome; does not replace it)
    groups = {}
    for g in GROUPS:
        g_unsup = [u for u in unsupported if group_of(u) == g]
        g_req = {k: v for k, v in req.items() if group_of(k) == g}
        g_miss = [t for gg, t in missing if gg in (g, 'all')]
        valid = correct and not g_unsup and not g_miss and all(v['value'] is not None for v in g_req.values())
        groups[g] = {'valid': bool(valid), 'unsupported': g_unsup, 'missing': g_miss,
                     'requirements': {k: {'value': v['value'], 'limit': v['limit'], 'unit': v['unit'],
                                          'ratio': ratios[k]} for k, v in g_req.items()},
                     'violated': sorted(k for k, v in g_req.items() if ratios[k] is not None and ratios[k] > 1.0),
                     'status': ('insufficient_evidence' if not valid else
                                ('requirements_violated' if any(ratios[k] > 1.0 for k in g_req) else
                                 'within_requirements'))}
    # descriptive: Amendment 2 re-measured times vs reused Amendment 1 times
    time_check = []
    for c in merged:
        if c.get('time_source') == 'A1' and 'a2_time_s' in c:
            tk = 'cpu_s' if c['block'] == 'S3' or c['args'].get('mode') == 'batch' else 'wall_s'
            time_check.append({'block': c['block'], 'args': c['args'], 'rep': c.get('rep', 0), 'key': tk,
                               'A1': c[tk], 'A2': c['a2_time_s'][tk], 'ratio_A2_over_A1': c['a2_time_s'][tk] / c[tk]})
    mem_change = []
    a1idx = {ckey(c): c for c in m1['cells'] if c.get('status') == 'ok'}
    for c in merged:
        if c.get('mem_source') == 'A2' and c.get('time_source') == 'A1':
            o = a1idx[ckey(c)]
            mem_change.append({'block': c['block'], 'args': c['args'], 'rep': c.get('rep', 0),
                               'A1_maxrss_kB': o['maxrss_kB'], 'A2_maxrss_kB': c['maxrss_kB'],
                               'delta_kB': c['maxrss_kB'] - o['maxrss_kB']})
    decomposition = {}
    for c in merged:
        if c['block'] == 'S4' and 'decomposition_s' in c:
            d = c['decomposition_s']
            key = (f"{c['args']['m']}/{c['args']['mode']}/{c['args']['role']}/nQ{c['args']['nq']}/nD{c['args']['nd']}"
                   f"/nwin{c['n_win']}")
            decomposition[key] = {k: d[k] / d['total'] for k in ('slice_copy', 'kernel', 'topk', 'loop_other')}
    out = {'_meta': {'producer': 'eb1a2_analyse.py', 'measurements_sha256': hashlib.sha256(IN2.read_bytes()).hexdigest(),
                     'amendment1_measurements_sha256': hashlib.sha256(IN1.read_bytes()).hexdigest(),
                     'original_measurements_sha256': hashlib.sha256(ORIG.read_bytes()).hexdigest(), 'P2': P2,
                     'limits': LIMITS, 'conf_max_rel_error': CONF_MAX, 'fit_resid_max': RESID_MAX,
                     'memory_sources': mem_source},
           'correctness': {'loader_check_A2': bool(lc.get('pass')),
                           'amendment1': bool(m1.get('correctness', {}).get('pass')),
                           'original': bool(orig.get('correctness', {}).get('pass')), 'pass': correct},
           'missing_cells': [t for _, t in missing], 'unsupported': unsupported, 'fits': fits, 'projections_simple': proj,
           'requirements': req, 'ratios_descriptive': ratios, 'decision': decision, 'groups': groups,
           'time_reuse_check_descriptive': time_check, 'memory_change_A1_to_A2': mem_change,
           'S4_time_share_decomposition': decomposition}
    OUT.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps(decision))
    print({g: v['status'] for g, v in groups.items()})
    return 0


if __name__ == '__main__':
    sys.exit(main())
