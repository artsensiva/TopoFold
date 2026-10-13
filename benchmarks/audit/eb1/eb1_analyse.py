#!/usr/bin/env python3
"""EB-1 analysis: fits, support rules, P2 projections, requirements R1–R4 and the frozen decision rule (§6–§7).

Written and hashed BEFORE any measurement. Interpretations of protocol details, fixed here in advance:
- repeated measurements of a cell are summarized by their median (time) and maximum (absolute peak RSS);
- time fits are per dataset (D1, D2) and per W/method where applicable; memory fits pool both datasets per stage
  (S1 per storage layout), since resident-array bytes differ strongly between datasets;
- R1/R2 are evaluated per retained method (A, B2) at W = 8; the engine is flagged if either method violates;
- the conservative P2 value is the larger of the D1- and D2-based projections where both exist.
Run: ../wpd_local_retrieval/build/venv/bin/python -B eb1_analyse.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
IN = RES / 'eb1_measurements.json'
OUT = RES / 'eb1_analysis.json'
RESID_MAX = 0.25
P2 = {'L': 1000, 'W': 8, 'n_win': 993, 'n_D': 10000, 'n_Q': 1000, 'frames_total': 11000,
      'atoms_xtc': 160000, 'atoms_h5_read': 4 * 1000}
LIMITS = {'R1_cpu_h': 96.0, 'R2_wall_s': 10.0, 'R3_wall_s': 7200.0, 'R4_rss_kB': 8 * 1024 ** 2}
BATCH, FEAT_BYTES_PER_RES = 128, 40


def fit(xs, ys):
    """Least squares y = c0 + c1 x. Supported iff ≥3 distinct x, c1 > 0, all |relative residual| ≤ 0.25."""
    x, y = np.asarray(xs, float), np.asarray(ys, float)
    out = {'n': int(len(x)), 'x': x.tolist(), 'y': y.tolist()}
    if len(np.unique(x)) < 3:
        return {**out, 'supported': False, 'reason': 'fewer than 3 distinct sizes'}
    A = np.vstack([np.ones_like(x), x]).T
    c0, c1 = np.linalg.lstsq(A, y, rcond=None)[0]
    pred = c0 + c1 * x
    rel = np.abs(y - pred) / np.maximum(np.abs(y), 1e-300)
    sup = bool(c1 > 0 and np.all(rel <= RESID_MAX))
    return {**out, 'c0': float(c0), 'c1': float(c1), 'max_rel_residual': float(rel.max()), 'supported': sup,
            'reason': None if sup else ('c1 <= 0' if c1 <= 0 else 'residual > 0.25')}


def project(f, x):
    if not f.get('supported'):
        return None
    return {'value': f['c0'] + f['c1'] * x, 'extrapolation_factor': x / max(f['x'])}


def decide(req, missing, correctness_pass, unsupported):
    """Frozen §7. req: {name: {'value', 'limit', 'component'}} (value None if unsupported)."""
    if not correctness_pass:
        return {'outcome': 'insufficient_evidence', 'reason': 'mandatory correctness failure'}
    if missing:
        return {'outcome': 'insufficient_evidence', 'reason': f'missing cells: {missing}'}
    if unsupported:
        return {'outcome': 'insufficient_evidence', 'reason': f'unsupported projections: {unsupported}'}
    flag = {'loader': 0.0, 'engine': 0.0}
    for r in req.values():
        ratio = r['value'] / r['limit']
        r['ratio'] = ratio
        r['violated'] = ratio > 1.0
        if r['violated']:
            flag[r['component']] = max(flag[r['component']], ratio)
    fl = [k for k, v in flag.items() if v > 0]
    if not fl:
        return {'outcome': 'neither'}
    if len(fl) == 1:
        return {'outcome': fl[0]}
    order = sorted(fl, key=lambda k: -flag[k])
    return {'outcome': 'both', 'order': order}


def med(vals):
    return float(np.median(vals))


def main():
    m = json.loads(IN.read_text())
    corr = m.get('correctness', {})
    cells = [c for c in m.get('cells', []) if c.get('status') == 'ok']
    planned = m['_meta'].get('n_cells_planned')
    missing = [] if planned is not None and len(cells) == planned else [
        f"{c['block']} {c['args']} rep{c['rep']}: {c.get('status')}" for c in m.get('cells', []) if c.get('status') != 'ok']
    if planned is not None and len(m.get('cells', [])) < planned:
        missing.append(f"{planned - len(m['cells'])} cells never run ({m.get('stopped')})")
    base = max(c['maxrss_kB'] for c in cells if c['block'] == 'M0') if any(c['block'] == 'M0' for c in cells) else None

    def group(block, **kw):
        return [c for c in cells if c['block'] == block and all(c['args'].get(k) == v for k, v in kw.items())]

    def tfit(block, xkey='x', tkey='wall_s', **kw):
        g = group(block, **kw)
        by = {}
        for c in g:
            by.setdefault(c[xkey], []).append(c[tkey])
        xs = sorted(by)
        return fit(xs, [med(by[x]) for x in xs])

    def mfit(blocks, bytes_key, **kw):
        g = [c for b in blocks for c in group(b, **kw)]
        by = {}
        for c in g:
            by.setdefault(c[bytes_key], []).append(c['maxrss_kB'])
        xs = sorted(by)
        return fit(xs, [max(by[x]) for x in xs])

    fits, proj, unsupported = {}, {}, []

    def need(name, f, x):
        fits[name] = f
        p = project(f, x)
        if p is None:
            unsupported.append(name)
        proj[name] = p
        return None if p is None else p['value']

    xtc = need('S1_xtc_wall', tfit('S1', d='D1'), P2['frames_total'] * P2['atoms_xtc'])
    h5 = need('S1_h5_wall', tfit('S1', d='D2'), P2['frames_total'] * P2['atoms_h5_read'])
    s2 = [need(f'S2_{d}_wall', tfit('S2', d=d), P2['frames_total'] * P2['L']) for d in ('D1', 'D2')]
    s3cpu = [need(f'S3_{d}_cpu', tfit('S3', tkey='cpu_s', d=d), P2['frames_total'] * P2['L']) for d in ('D1', 'D2')]
    s3one = [need(f'S3_{d}_wall_one_frame', tfit('S3', d=d), P2['L']) for d in ('D1', 'D2')]
    req = {}
    for meth in ('A', 'B2'):
        single = [need(f'S4single_{d}_W8_{meth}_wall', tfit('S4single', d=d, W=8, m=meth), P2['n_win'] * P2['n_D'])
                  for d in ('D1', 'D2')]
        batch = [need(f'S4batch_{d}_W8_{meth}_cpu', tfit('S4batch', tkey='cpu_s', d=d, W=8, m=meth),
                      P2['n_Q'] * P2['n_win'] * P2['n_D']) for d in ('D1', 'D2')]
        r1 = None if None in s3cpu + batch else (max(s3cpu) + max(batch)) / 3600
        r2 = None if None in s3one + single else max(s3one) + max(single)
        req[f'R1_{meth}'] = {'value': r1, 'limit': LIMITS['R1_cpu_h'], 'component': 'engine', 'unit': 'CPU-h'}
        req[f'R2_{meth}'] = {'value': r2, 'limit': LIMITS['R2_wall_s'], 'component': 'engine', 'unit': 's'}
    for lay, v in (('xtc', xtc), ('h5', h5)):
        req[f'R3_{lay}'] = {'value': None if v is None or None in s2 else v + max(s2), 'limit': LIMITS['R3_wall_s'],
                            'component': 'loader', 'unit': 's'}
    # R4 memory (absolute peak RSS)
    b2 = {'S1_xtc': P2['frames_total'] * P2['L'] * 12, 'S1_h5': P2['frames_total'] * P2['L'] * 12,
          'S2': P2['frames_total'] * P2['L'] * 60, 'S3': P2['frames_total'] * P2['L'] * FEAT_BYTES_PER_RES,
          'S4': P2['frames_total'] * P2['L'] * FEAT_BYTES_PER_RES + BATCH * P2['n_D'] * 25 * 8}
    mem = {'S1_xtc': need('R4_S1_xtc', mfit(['S1'], 'out_bytes', d='D1'), b2['S1_xtc']),
           'S1_h5': need('R4_S1_h5', mfit(['S1'], 'out_bytes', d='D2'), b2['S1_h5']),
           'S2': need('R4_S2', mfit(['S2'], 'resident_bytes'), b2['S2']),
           'S3': need('R4_S3', mfit(['S3'], 'resident_bytes'), b2['S3']),
           'S4': need('R4_S4', mfit(['S4single', 'S4batch'], 'resident_bytes'), b2['S4'])}
    for k, v in mem.items():
        req[f'R4_{k}'] = {'value': v, 'limit': LIMITS['R4_rss_kB'], 'component': 'loader' if k.startswith(('S1', 'S2'))
                          else 'engine', 'unit': 'kB'}
    decision = decide(req, missing, bool(corr.get('pass')), unsupported)
    # descriptive shares on the measured development workloads (operational path, single process)
    shares = {}
    out = {'_meta': {'producer': 'eb1_analyse.py', 'measurements_sha256': __import__('hashlib').sha256(IN.read_bytes()).hexdigest(),
                     'P2': P2, 'limits': LIMITS, 'baseline_import_rss_kB': base},
           'correctness_pass': bool(corr.get('pass')), 'missing_cells': missing, 'unsupported_projections': unsupported,
           'fits': fits, 'projections': proj, 'requirements': req, 'decision': decision, 'shares': shares,
           'research_only': {c['args']['part']: {'wall_s': c['wall_s'], 'cpu_s': c['cpu_s'], 'maxrss_kB': c['maxrss_kB']}
                             for c in group('R')}}
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(decision))
    return 0


if __name__ == '__main__':
    sys.exit(main())
