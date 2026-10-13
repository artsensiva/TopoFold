#!/usr/bin/env python3
"""EB-1 Amendment 1 analysis (docs/EB1_AMENDMENT_1.md §6–§7). Written and hashed before any amended measurement.

Reuses the original decision function (eb1_analyse.decide) and limits unchanged. Writes results/eb1a1_analysis.json.
Run: ../wpd_local_retrieval/build/venv/bin/python -B eb1a1_analyse.py
"""
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import eb1_analyse as AN  # noqa: E402
import eb1a1_lib as A1  # noqa: E402

RES = HERE / 'results'
IN = RES / 'eb1a1_measurements.json'
ORIG = RES / 'eb1_measurements.json'
OUT = RES / 'eb1a1_analysis.json'
COEF = json.loads((RES / 'eb1a1_temp_calibration.json').read_text())['coef']
P2, LIMITS, RESID_MAX, CONF_MAX = AN.P2, AN.LIMITS, 0.25, 0.25
W = 8


def lin_fit(X, y):
    """Least squares y = X @ c (X includes an intercept column). Returns coefficients and fit diagnostics."""
    X, y = np.asarray(X, float), np.asarray(y, float)
    c = np.linalg.lstsq(X, y, rcond=None)[0]
    rel = np.abs(y - X @ c) / np.abs(y)
    return c, float(rel.max())


def s4_model(cells, kind, tkey):
    """kind: 'batch_time', 'single_time', 'mem'. Returns model dict with coefficients and support (fit part)."""
    if kind == 'mem':
        X = [[1.0, c['P']] for c in cells]
        y = [c['maxrss_kB'] * 1024.0 for c in cells]
        names = ['m0', 'm1']
    elif kind == 'batch_time':
        X = [[1.0, c['comparisons'], c['call_db_rows']] for c in cells]
        y = [c[tkey] for c in cells]
        names = ['c0', 'c1', 'c2']
    else:
        X = [[1.0, c['comparisons'], c['n_win']] for c in cells]
        y = [c[tkey] for c in cells]
        names = ['c0', 'c1', 'c2']
    n_distinct = len({tuple(r[1:]) for r in X})
    c, res = lin_fit(X, y)
    if kind != 'mem' and c[2] < 0:                          # pre-declared: drop the c2 term
        X = [r[:2] for r in X]
        c, res = lin_fit(X, y)
        names = names[:2]
    coef = dict(zip(names, map(float, c)))
    lead = coef['m1'] if kind == 'mem' else coef['c1']
    sup = n_distinct >= 3 and lead > 0 and res <= RESID_MAX
    return {'kind': kind, 'coef': coef, 'n_fit': len(cells), 'n_distinct': n_distinct, 'max_fit_rel_residual': res,
            'fit_supported': bool(sup)}


def predict(model, feats):
    c = model['coef']
    if model['kind'] == 'mem':
        return c['m0'] + c['m1'] * feats['P']
    second = feats['call_db_rows'] if model['kind'] == 'batch_time' else feats['n_win']
    return c['c0'] + c['c1'] * feats['comparisons'] + c.get('c2', 0.0) * second


def main():
    m = json.loads(IN.read_text())
    orig = json.loads(ORIG.read_text())
    cells = [c for c in m.get('cells', []) if c.get('status') == 'ok']
    planned = m['_meta'].get('n_cells_planned')
    missing = []
    if planned is None or len(cells) != planned:
        missing = [f"{c['block']} {c['args']}: {c.get('status')}" for c in m.get('cells', []) if c.get('status') != 'ok']
        if planned is not None and len(m.get('cells', [])) < planned:
            missing.append(f"{planned - len(m['cells'])} cells never run ({m.get('stopped')})")
    correct = bool(m.get('correctness', {}).get('pass')) and bool(orig.get('correctness', {}).get('pass'))
    fits, unsupported, proj, conf_report = {}, [], {}, {}

    def med_by(rows, xkey, ykey):
        by = {}
        for r in rows:
            by.setdefault(r[xkey], []).append(r[ykey])
        xs = sorted(by)
        return xs, [float(np.median(by[x])) for x in xs]

    def simple(name, rows, xkey, ykey, x_p2, scale=1.0):
        xs, ys = med_by(rows, xkey, ykey)
        f = AN.fit(xs, [y * scale for y in ys])
        fits[name] = f
        p = AN.project(f, x_p2)
        proj[name] = p
        if p is None:
            unsupported.append(name)
        return None if p is None else p['value']

    # S1 (reused, original measurements; corrected memory predictor for D2)
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
    # S2 and S3 (re-measured)
    s2 = [simple(f'S2_{d}_wall', [c for c in cells if c['block'] == 'S2' and c['args']['d'] == d], 'x', 'wall_s',
                 P2['frames_total'] * P2['L']) for d in ('D1', 'D2')]
    s2m = simple('R4_S2', [c for c in cells if c['block'] == 'S2'], 'P', 'maxrss_kB', P2['frames_total'] * P2['L'] * 72, 1024.0)
    s3cpu = [simple(f'S3_{d}_cpu', [c for c in cells if c['block'] == 'S3' and c['args']['d'] == d], 'x', 'cpu_s',
                    P2['frames_total'] * P2['L']) for d in ('D1', 'D2')]
    s3one = [simple(f'S3_{d}_wall_one_frame', [c for c in cells if c['block'] == 'S3' and c['args']['d'] == d], 'x',
                    'wall_s', P2['L']) for d in ('D1', 'D2')]
    s3m = simple('R4_S3', [c for c in cells if c['block'] == 'S3'], 'P', 'maxrss_kB',
                 P2['frames_total'] * P2['L'] * 36 + P2['frames_total'] * (2 * P2['L'] - 5) * 8, 1024.0)
    # S4 per method and mode, with confirmation
    s4 = {}
    for meth in ('A', 'B2'):
        for mode in ('batch', 'single'):
            rows = [c for c in cells if c['block'] == 'S4' and c['args']['m'] == meth and c['args']['mode'] == mode]
            fit_rows = [c for c in rows if c['args']['role'] == 'fit']
            conf_rows = [c for c in rows if c['args']['role'] != 'fit']
            b = A1.effective_batch(P2['n_Q'] if mode == 'batch' else 1)
            nq = P2['n_Q'] if mode == 'batch' else 1
            wk = A1.work(nq, P2['n_D'], P2['n_win'])
            p2_P = (A1.resident_bytes(meth, P2['frames_total'], P2['L']) + A1.retained_bytes(nq, P2['n_win'])
                    + A1.temp_bytes(meth, b, P2['n_D'], W, COEF))
            p2_feats = {**wk, 'n_win': P2['n_win'], 'P': p2_P}
            for kind, tkey in ((f'{mode}_time', 'cpu_s' if mode == 'batch' else 'wall_s'), ('mem', None)):
                name = f"S4_{meth}_{mode}_{'mem' if kind == 'mem' else 'time'}"
                if not fit_rows:
                    unsupported.append(name)
                    continue
                mod = s4_model(fit_rows, kind, tkey)
                errs = []
                for c in conf_rows:
                    meas = c['maxrss_kB'] * 1024.0 if kind == 'mem' else c[tkey]
                    pred = predict(mod, c)
                    errs.append({'workload': c['args']['role'], 'measured': meas, 'predicted': pred,
                                 'rel_error': abs(pred - meas) / meas})
                mod['confirmation'] = errs
                mod['confirmation_supported'] = bool(errs) and all(e['rel_error'] <= CONF_MAX for e in errs)
                mod['supported'] = mod['fit_supported'] and mod['confirmation_supported']
                max_fit_x = max(c['comparisons'] for c in fit_rows + conf_rows) if kind != 'mem' else max(
                    c['P'] for c in fit_rows + conf_rows)
                val = predict(mod, p2_feats)
                mod['P2'] = {'value': val, 'extrapolation_factor_vs_largest_measured':
                             (p2_feats['comparisons'] if kind != 'mem' else p2_P) / max_fit_x, 'P2_predictor_bytes': p2_P}
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
    decision = AN.decide(json.loads(json.dumps(req)), missing, correct, unsupported)
    ratios = {k: (v['value'] / v['limit'] if v['value'] is not None else None) for k, v in req.items()}
    decomposition = {}
    for c in cells:
        if c['block'] == 'S4':
            d = c['decomposition_s']
            key = f"{c['args']['m']}/{c['args']['mode']}/{c['args']['role']}/nQ{c['args']['nq']}/nD{c['args']['nd']}/nwin{c['n_win']}"
            decomposition[key] = {k: d[k] / d['total'] for k in ('slice_copy', 'kernel', 'topk', 'loop_other')}
    out = {'_meta': {'producer': 'eb1a1_analyse.py', 'measurements_sha256': hashlib.sha256(IN.read_bytes()).hexdigest(),
                     'original_measurements_sha256': hashlib.sha256(ORIG.read_bytes()).hexdigest(), 'P2': P2,
                     'limits': LIMITS, 'conf_max_rel_error': CONF_MAX},
           'correctness_pass': correct, 'missing_cells': missing, 'unsupported': unsupported, 'fits': fits,
           'projections_simple': proj, 'requirements': req, 'ratios_descriptive': ratios, 'decision': decision,
           'S4_time_share_decomposition': decomposition}
    OUT.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps(decision))
    return 0


if __name__ == '__main__':
    sys.exit(main())
