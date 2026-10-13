#!/usr/bin/env python3
"""Generate the closed SYNTHETIC fixture matrix of protocol §5 (with amendments §9).

Writes regenerable binaries to build/fixtures/ and the manifest fixtures/SYNTHETIC_fixtures.json
(parameters, seeds, RNG draw order, SHA-256). Descriptor-level constructions only; the traces are
geometric test objects, not proteins.
"""
import sys

import numpy as np

from build_traces import build_traces
from common import (FIXTURE_BIN, FIXTURES_JSON, FIXTURES_TSV, KAPPA, dump_json, wrap,
                    write_f64)

DEG = np.pi / 180.0
SEED_BASE = 20261004
ALPHAS_DEG = [0, 45, 90, 135, 170, 179, 180, -179]
F3_EPS = [1e-1, 1e-3, 1e-6, 1e-9, 0.0, -1e-9, -1e-6, -1e-3, -1e-1]
F3_NOISE = [('0', 0.0), ('1e-6rad', 1e-6), ('1deg', 1.0 * DEG)]
F4_P = [0.5, 0.49, 0.45]
F4_EPS = [0.0, 1e-6, -1e-6, 1e-3, -1e-3, 1e-1, -1e-1]
F4_NOISE = [('0', 0.0), ('1deg', 1.0 * DEG)]


def efmt(x):
    return '0' if x == 0 else f'{x:+.0e}'


class Matrix:
    def __init__(self):
        self.cases = []
        self.k = 0

    def rng(self):
        seed = SEED_BASE + self.k
        self.k += 1
        return seed, np.random.default_rng(seed)

    def add(self, cid, family, tors, W, params, seed, draws, coords=False, perturb=False,
            derived_from=None, alpha_deg=None, kind='angular', labels=None):
        self.cases.append(dict(id=cid, family=family, kind=kind, W=W, tors=tors, params=params,
                               seed=seed, rng_draw_order=draws, coords=coords, perturb=perturb,
                               derived_from=derived_from, alpha_deg=alpha_deg, labels=labels))


def f1(m):
    F = 1000
    specs = [('F1a', 5, [30, 90], [90, 30]),
             ('F1b', 8, [30, 90, 60, 60, 60], [90, 30, 60, 60, 60]),
             ('F1c', 5, [30, 30], [90, 90])]
    bases = {}
    for name, W, A, B in specs:
        for sig_name, sig in [('s0', 0.0), ('s2deg', 2.0 * DEG)]:
            seed, rng = m.rng()
            labels = np.arange(F) % 2                      # 0 = A, 1 = B, interleaved
            clean = np.where(labels[:, None] == 0, np.array(A) * DEG, np.array(B) * DEG)
            draws = []
            noise = np.zeros_like(clean)
            if sig > 0:
                noise = rng.normal(0.0, sig, size=clean.shape)
                draws = [f'normal(0, {sig!r}, size=({F}, {W - 3})) iid per torsion']
            tors = wrap(clean + noise)
            cid = f'{name}_{sig_name}'
            m.add(cid, 'F1', tors, W, dict(A_deg=A, B_deg=B, sigma_rad=sig, F=F,
                                           states='frame % 2 (0=A,1=B)'),
                  seed, draws, coords=True, labels=labels)
            bases[cid] = tors
    return bases


def f2(m, f1_bases):
    F = 1000
    bases = {'F2a': (f1_bases['F1a_s2deg'], None, ['reuses F1a_s2deg draws (paired)'], {})}
    for name, centre, sig in [('F2b', 40.0, 10.0), ('F2c', 100.0, 30.0)]:
        seed, rng = m.rng()
        state = rng.random(F) < 0.5
        x = np.where(state, centre * DEG, -centre * DEG) + rng.normal(0.0, sig * DEG, size=F)
        tors = wrap(np.stack([x, x], axis=1))
        bases[name] = (tors, seed, [f'random({F}) < 0.5 -> state (+{centre} deg if True)',
                                    f'normal(0, {sig} deg in rad, size={F}) per frame'],
                       dict(centres_deg=[-centre, centre], sigma_deg=sig, F=F))
    for name, (tors, seed, draws, params) in bases.items():
        for a in ALPHAS_DEG:
            rot = wrap(tors + a * DEG)
            m.add(f'{name}_a{a:+04d}', 'F2', rot, 5, dict(params, alpha_deg=a), seed, draws,
                  coords=True, derived_from=name, alpha_deg=a)


def f3(m):
    F, th0 = 200, 20.0 * DEG
    for eps in F3_EPS:
        for nname, sn in F3_NOISE:
            seed, rng = m.rng()
            eta = np.zeros(F)
            draws = []
            if sn > 0:
                eta = rng.normal(0.0, sn, size=F)
                draws = [f'normal(0, {sn!r}, size={F}) per frame, added to both torsions']
            t1 = th0 + eta
            t2 = th0 + eta + np.pi - eps
            tors = wrap(np.stack([t1, t2], axis=1))
            m.add(f'F3_e{efmt(eps)}_n{nname}', 'F3', tors, 5,
                  dict(theta0_rad=th0, eps_rad=eps, sigma_n_rad=sn, F=F), seed, draws,
                  coords=True, perturb=True)


def f4(m):
    F = 1000
    specs = [(p, eps, np.pi - eps, 'antipodal') for p in F4_P for eps in F4_EPS]
    specs += [(0.5, None, 120.0 * DEG, 'control120')]
    for p, eps, phi2, tag in specs:
        for nname, sig in F4_NOISE:
            seed, rng = m.rng()
            perm = rng.permutation(F)
            n2 = int(round(p * F))
            state = np.zeros(F, dtype=bool)
            state[perm[:n2]] = True
            draws = [f'permutation({F}); first {n2} indices -> phi2']
            nu = np.zeros(F)
            if sig > 0:
                nu = rng.normal(0.0, sig, size=F)
                draws.append(f'normal(0, {sig!r}, size={F}) per frame, added to both torsions')
            x = np.where(state, phi2, 0.0) + nu
            tors = wrap(np.stack([x, x], axis=1))
            base = (f'F4_p{p:.2f}_e{efmt(eps)}_n{nname}' if tag == 'antipodal'
                    else f'F4_control120_n{nname}')
            params = dict(p=p, eps_rad=eps, phi2_rad=phi2, sigma_rad=sig, F=F, n_phi2=n2)
            for a in ALPHAS_DEG:
                m.add(f'{base}_a{a:+04d}', 'F4', wrap(tors + a * DEG), 5,
                      dict(params, alpha_deg=a), seed, draws, coords=(a == 0),
                      derived_from=base, alpha_deg=a)


def f6(m):
    """Replay bc_check.py draw order on default_rng(1) without executing that script."""
    N, W = 2500, 5
    rng = np.random.default_rng(1)
    for mu, sd in [(-170, 15), (-175, 15), (180, 10), (-120, 15), (50, 10)]:
        rng.normal(mu, sd, size=(N, W))                   # section 1: consumed, discarded
    draws_common = ['default_rng(1)',
                    'section 1 replay: 5 x normal(mu, sd, size=(2500, 5)) discarded']
    rare = {}
    for p in [0.5, 0.1, 0.02, 0.005]:
        x = rng.normal(0, 0.05, N)
        x[rng.random(N) < p] += 1.0
        rare[p] = x
        m.add(f'F6a_rare_p{p}', 'F6a', x[:, None], None, dict(p=p, N=N), 1,
              draws_common + [f'section 2: normal(0, 0.05, {N}); random({N}) < {p} -> +1.0'],
              kind='scalar')
    for name, x in [('exponential', rng.exponential(1, N)),
                    ('lognormal_s0.8', rng.lognormal(0, 0.8, N)),
                    ('gamma_k1', rng.gamma(1, 1, N)),
                    ('chi2_df1', rng.chisquare(1, N)),
                    ('gaussian', rng.normal(size=N))]:
        m.add(f'F6a_{name}', 'F6a', x[:, None], None, dict(N=N), 1,
              draws_common + ['section 3 list order: exponential, lognormal, gamma, chisquare, '
                              'normal'], kind='scalar')
    for p, x in rare.items():
        assert np.all(np.abs(x) < np.pi), 'F6b assumption: values inside (−π, π]'
        m.add(f'F6b_rare_p{p}', 'F6b', np.stack([x, x], axis=1), 5, dict(p=p, F=N), 1,
              ['values of F6a_rare_p%s used as both torsions' % p], coords=True,
              derived_from=f'F6a_rare_p{p}')


def main():
    m = Matrix()
    b1 = f1(m)
    f2(m, b1)
    f3(m)
    f4(m)
    f6(m)

    records, tsv = [], []
    for c in m.cases:
        tors = np.asarray(c['tors'], dtype=np.float64)
        F, ncols = tors.shape
        tpath = f"{c['id']}.tors.f64"
        rec = {k: c[k] for k in ['id', 'family', 'kind', 'W', 'params', 'seed', 'rng_draw_order',
                                 'derived_from', 'alpha_deg', 'perturb']}
        rec.update(F=F, ncols=ncols, kappa=KAPPA if c['kind'] == 'angular' else None,
                   torsion_file=tpath, torsion_sha256=write_f64(FIXTURE_BIN / tpath, tors),
                   coords_file=None, coords_sha256=None, r2_p1=bool(c['coords']))
        if c['labels'] is not None:
            rec['labels_file'] = f"{c['id']}.labels.f64"
            rec['labels_sha256'] = write_f64(FIXTURE_BIN / rec['labels_file'], c['labels'])
        if c['coords']:
            W = c['W']
            X = build_traces(np.full((F, W - 2), KAPPA), tors)
            rec['coords_file'] = f"{c['id']}.coords.f64"
            rec['coords_sha256'] = write_f64(FIXTURE_BIN / rec['coords_file'], X)
        records.append(rec)
        tsv.append('\t'.join(str(v) for v in [
            rec['id'], rec['kind'], F, ncols, rec['W'] or 0, rec['kappa'] or 0.0,
            tpath, rec['coords_file'] or '-', int(rec['perturb'])]))

    ids = [r['id'] for r in records]
    assert len(ids) == len(set(ids)), 'duplicate case ids'
    FIXTURES_TSV.parent.mkdir(parents=True, exist_ok=True)
    FIXTURES_TSV.write_text('\n'.join(tsv) + '\n')
    dump_json({'_meta': {'type': 'SYNTHETIC', 'protocol': 'docs/ANGULAR_REPRESENTATION_'
                         'CHARACTERIZATION_PROTOCOL.md §5, §9',
                         'seed_rule': f'default_rng({SEED_BASE} + k), k = base-fixture index; '
                                      'rotations reuse base draws; F6 uses default_rng(1) replay',
                         'binary_format': 'little-endian float64, C order; torsions (F, ncols), '
                                          'coords (F, W, 3) in Å',
                         'kappa_rad': KAPPA, 'bond_A': 3.8,
                         'n_cases': len(records),
                         'n_r2_p1_cases': sum(r['r2_p1'] for r in records)},
               'cases': records}, FIXTURES_JSON)
    print(f'{len(records)} cases ({sum(r["r2_p1"] for r in records)} with coordinates) -> '
          f'{FIXTURES_JSON}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
