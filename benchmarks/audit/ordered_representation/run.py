#!/usr/bin/env python3
"""WP-C: ordered (cos τ, sin τ) representation on SYNTHETIC fixtures (protocol §4–§9, §11).

Default: writes N1/N2 controls and comparator outputs into build/, computes all checks and
retrieval scores, writes results/ and provenance. `--verify`: recomputes from inputs (no writes)
and exits nonzero if the stored results differ or a case is missing.

Run with the WP-A isolated interpreter (it holds the freshly built binding):
  benchmarks/audit/angular_characterization/build/venv/bin/python -B \
      benchmarks/audit/ordered_representation/run.py
"""
import argparse
import datetime
import hashlib
import json
import math
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np

from controls import N1_EPS, build_controls
from features import (bond_labels, ca_distance_vectors, d_chord, d_circ, d_euclid_rms, d_scalar,
                      feature_contracts, ordered_features, wrap)
from inputs import (REPO, WPA, InputHashError, arr_sha, f1_labels, f4_labels, load_case, manifest,
                    sha_bytes, wpa_build, wpa_common)
from retrieval import knn_purity

HERE = Path(__file__).resolve().parent
BUILD = HERE / 'build'
CTRL_DIR = BUILD / 'controls'
CMP_DIR = BUILD / 'comparator'
CTRL_JSON = HERE / 'fixtures' / 'SYNTHETIC_wpc_controls.json'
RES_JSON = HERE / 'results' / 'SYNTHETIC_ordered_representation_results.json'
PROV_JSON = HERE / 'results' / 'SYNTHETIC_ordered_representation_provenance.json'
KS = (1, 10)
MARGIN = 0.05                      # protocol §9 descriptive margin (not a production threshold)
KAPPA = wpa_common.KAPPA
REPS = ['D_ch', 'D_circ', 'scalar_comparator', 'ca_distance', 'D_kappa']

log = []


def run_cmd(cmd, cwd=REPO):
    p = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    log.append({'cmd': ' '.join(str(c) for c in cmd), 'exit': p.returncode})
    print(f'[exit {p.returncode}] {log[-1]["cmd"][:160]}')
    return p


def circ(a, b):
    return np.abs(wrap(np.asarray(a) - np.asarray(b)))


def case_lists():
    man = manifest()
    fam = lambda f, cond=lambda c: True: sorted(i for i, c in man.items()
                                                if c['family'] == f and cond(c))
    return {'F1': fam('F1'), 'F2': fam('F2'), 'F3': fam('F3'),
            'F4a0': fam('F4', lambda c: c['alpha_deg'] == 0), 'F6b': fam('F6b')}


# --- inputs ---------------------------------------------------------------------------------------

def controls(write):
    ctrl = build_controls()
    recs = {}
    for cid, c in ctrl.items():
        recs[cid] = {'W': c['W'], 'F': int(c['tors'].shape[0]), 'params': c['params'],
                     'draw_order': c['draw_order'], 'torsion_sha256': arr_sha(c['tors']),
                     'coords_sha256': arr_sha(c['coords']), 'labels_sha256': arr_sha(
                         c['labels'].astype(np.float64)),
                     'torsion_file': f'{cid}.tors.f64', 'coords_file': f'{cid}.coords.f64'}
    if write:
        CTRL_DIR.mkdir(parents=True, exist_ok=True)
        for cid, c in ctrl.items():
            np.ascontiguousarray(c['tors'], '<f8').tofile(CTRL_DIR / recs[cid]['torsion_file'])
            np.ascontiguousarray(c['coords'], '<f8').tofile(CTRL_DIR / recs[cid]['coords_file'])
        wpa_common.dump_json({'_meta': {'type': 'SYNTHETIC', 'protocol':
                                        'docs/ORDERED_ANGULAR_REPRESENTATION_PROTOCOL.md §4'},
                              'cases': recs}, CTRL_JSON)
    else:
        stored = wpa_common.load_json(CTRL_JSON)['cases']
        for cid, r in recs.items():
            if stored[cid]['torsion_sha256'] != r['torsion_sha256'] or \
                    stored[cid]['coords_sha256'] != r['coords_sha256']:
                raise InputHashError(f'control {cid} does not reproduce its recorded hash')
    return ctrl, recs


def gather_inputs():
    lists = case_lists()
    cases, sources = {}, {}
    for fam, ids in lists.items():
        for cid in ids:
            tors, coords, src = load_case(cid)
            cases[cid] = {'family': fam, 'tors': tors, 'coords': coords, 'W': manifest()[cid]['W']}
            sources[cid] = src
    for cid in lists['F1']:
        cases[cid]['labels'] = f1_labels(cid)
    for cid in lists['F4a0']:
        cases[cid]['labels'] = f4_labels(cid, cases[cid]['tors'])
    return lists, cases, sources


# --- historical comparator: unchanged WP-A Rust example (production functions) ---------------------

def comparator_cases(lists, ctrl):
    return lists['F1'] + lists['F3'] + lists['F4a0'] + sorted(ctrl)


def run_comparator(lists, cases, ctrl, recs):
    CMP_DIR.mkdir(parents=True, exist_ok=True)
    regen = BUILD / 'regenerated_wpa_inputs'
    rows = []
    man = manifest()
    for cid in comparator_cases(lists, ctrl):
        if cid in ctrl:
            t = (CTRL_DIR / recs[cid]['torsion_file']).relative_to(REPO)
            x = (CTRL_DIR / recs[cid]['coords_file']).relative_to(REPO)
            F, W = recs[cid]['F'], recs[cid]['W']
        else:
            c = man[cid]
            t = WPA / 'build' / 'fixtures' / c['torsion_file']
            x = WPA / 'build' / 'fixtures' / c['coords_file']
            if not t.exists():          # scratch copy of in-memory regenerated (hash-verified) input
                regen.mkdir(parents=True, exist_ok=True)
                t, x = regen / c['torsion_file'], regen / c['coords_file']
                np.ascontiguousarray(cases[cid]['tors'], '<f8').tofile(t)
                np.ascontiguousarray(cases[cid]['coords'], '<f8').tofile(x)
            t, x = t.relative_to(REPO), x.relative_to(REPO)
            F, W = c['F'], c['W']
        rows.append('\t'.join(map(str, [cid, 'angular', F, W - 3, W, KAPPA, t, x, 0])))
    tsv = CMP_DIR / 'comparator_manifest.tsv'
    tsv.write_text('\n'.join(rows) + '\n')
    p = run_cmd(['cargo', 'run', '--release', '-q', '-p', 'topofold-core', '--example',
                 'angular_characterization', '--', tsv.relative_to(REPO), '.',
                 (CMP_DIR / 'comparator_rust.json').relative_to(REPO),
                 (CMP_DIR / 'comparator_arrays.f64').relative_to(REPO)])
    if p.returncode != 0:
        print(p.stdout[-2000:], p.stderr[-2000:])
        raise SystemExit('comparator run failed')


def load_comparator(expected=None):
    jpath, apath = CMP_DIR / 'comparator_rust.json', CMP_DIR / 'comparator_arrays.f64'
    hashes = {'json': sha_bytes(jpath.read_bytes()), 'arrays': sha_bytes(apath.read_bytes())}
    if expected is not None and hashes != expected:
        raise InputHashError('comparator outputs do not match the hashes recorded in results')
    doc = json.loads(jpath.read_text())
    arrays = np.fromfile(apath, '<f8')
    out = {}
    for cid, v in doc['cases'].items():
        off, n = v['r1']['means']
        out[cid] = {'means': arrays[off:off + n], 'scan_bc_tau_r2': v['r2']['profile']['bc_tau'],
                    'scan_bc_tau_r1': v['r1']['moments']['bc'],
                    'r2_bits_equal': v['r2']['replication']['bc_tau_bits_equal']}
    return out, hashes


# --- binding identity ---------------------------------------------------------------------------

def binding_identity():
    import topofold
    import topofold.topofold as ext
    so = Path(ext.__file__).resolve()
    so_sha = sha_bytes(so.read_bytes())
    wpa_prov = wpa_common.load_json(WPA / 'results' /
                                    'SYNTHETIC_angular_characterization_provenance.json')
    src_now = {k: sha_bytes((REPO / k).read_bytes()) for k in wpa_prov['source_sha256']
               if k.startswith('crates/') and '/src/' in k or k == 'Cargo.lock'}
    src_same = all(wpa_prov['source_sha256'][k] == v for k, v in src_now.items())
    return topofold, {'extension_path': str(so), 'extension_sha256': so_sha,
                      'equals_wpa_freshly_built_extension':
                          so_sha == wpa_prov['extension']['loaded_extension_sha256'],
                      'production_sources_unchanged_since_wpa_build': src_same,
                      'compared_sources': sorted(src_now)}


# --- per-case computation -----------------------------------------------------------------------

def recover_paths(tf, tors0, coords):
    kap2, tau2 = wpa_build.recover(coords)
    F, W = coords.shape[:2]
    tau1 = np.empty_like(tors0)
    kap1 = np.empty((F, W - 2))
    for f in range(F):
        k, t, _ = tf.compute_invariants(coords[f].astype(np.float32))
        tau1[f], kap1[f] = np.asarray(t)[:W - 3], np.asarray(k)[:W - 2]
    return {'T2': (tau2, kap2), 'T1': (tau1, kap1)}


def feature_recovery_check(tors0, tau_t, delta):
    dx = np.abs(ordered_features(tau_t) - ordered_features(tors0)).max()
    return float(dx), bool(dx <= delta + 1e-15)


def retrieval_block(tors, kap, ca_vecs, means, labels):
    mats = {'D_ch': d_chord(tors, tors), 'D_circ': d_circ(tors, tors),
            'scalar_comparator': d_scalar(means, means),
            'ca_distance': d_euclid_rms(ca_vecs, ca_vecs),
            'D_kappa': d_euclid_rms(kap, kap)}
    scores = {r: {f'k{k}': knn_purity(mats[r], labels, k) for k in KS} for r in REPS}
    return mats, scores


def pair_chord(u, v):
    """Row-wise D_ch between matching frames of u and v (F, m)."""
    return np.sqrt(np.sum((ordered_features(u) - ordered_features(v)) ** 2, axis=-1) / u.shape[1])


def bounds_ok(Dch, Dci):
    return bool(np.all(Dch <= Dci + 1e-12) and np.all(Dch >= (2 / math.pi) * Dci - 1e-12))


def compute(tf, lists, cases, ctrl, cmp):
    R = {'checks': {}, 'cases': {}, 'retrieval': {}}
    chk = R['checks']

    # V2a feature contracts (T0 functions)
    chk['feature_contracts'] = [{'name': n, 'passed': ok, 'detail': d}
                                for n, ok, d in feature_contracts()]

    # V2b analytic distance checks independent of fixtures
    rng = np.random.default_rng(0)
    a = rng.uniform(-np.pi, np.pi, 2000)
    b = np.concatenate([a[:1000] + rng.uniform(-1e-9, 1e-9, 1000),
                        a[1000:] + np.pi - rng.uniform(0, 1e-9, 1000)])
    lhs = np.linalg.norm(ordered_features(a[:, None]) - ordered_features(b[:, None]), axis=1)
    rhs = 2 * np.abs(np.sin(circ(a, b) / 2))
    chk['AD1_chord_identity'] = {'max_abs_err': float(np.abs(lhs - rhs).max()),
                                 'passed': bool(np.abs(lhs - rhs).max() <= 1e-12)}
    br = float(d_chord(np.array([[np.pi - 1e-12]]), np.array([[-np.pi + 1e-12]]))[0, 0])
    chk['AD7_branch'] = {'D_ch': br, 'passed': br <= 3e-12}

    rec_fail = []
    # every case with coordinates: recovery-aware checks (T2, T1)
    paths = {}
    all_cases = {**cases, **{k: dict(v, family='N1' if k.startswith('N1') else 'N2')
                             for k, v in ctrl.items()}}
    for cid, c in all_cases.items():
        if c['coords'] is None:
            continue
        p = recover_paths(tf, c['tors'], c['coords'])
        paths[cid] = p
        entry = {'family': c['family']}
        for src, (tau, kap) in p.items():
            delta = float(circ(tau, c['tors']).max())
            dk = float(np.abs(kap - KAPPA).max())
            dx, ok = feature_recovery_check(c['tors'], tau, delta)
            entry[src] = {'max_tau_recovery_err': delta, 'max_kappa_recovery_err': dk,
                          'max_feature_dev': dx, 'feature_within_delta': ok}
            if not ok or (src == 'T2' and delta > wpa_common.R2_RECOVERY_TOL):
                rec_fail.append((cid, src))
        R['cases'][cid] = entry
    chk['recovery'] = {'failures': rec_fail, 'passed': not rec_fail,
                       'max_T2_err': max(v['T2']['max_tau_recovery_err'] for v in R['cases'].values()),
                       'max_T1_err': max(v['T1']['max_tau_recovery_err'] for v in R['cases'].values())}

    # AD3: F1 noiseless analytic distances, scalar comparator from production means
    exp = {'F1a_s0': (1.0, math.pi / 3, 0.0), 'F1b_s0': (math.sqrt(0.4), math.pi / 3 * math.sqrt(0.4), 0.0),
           'F1c_s0': (1.0, math.pi / 3, math.pi / 3)}
    ad3 = {}
    for cid, (ech, eci, esc) in exp.items():
        t = cases[cid]['tors']
        ch, ci = float(d_chord(t[:1], t[1:2])[0, 0]), float(d_circ(t[:1], t[1:2])[0, 0])
        sc = float(circ(cmp[cid]['means'][0], cmp[cid]['means'][1]))
        ad3[cid] = {'D_ch': ch, 'D_circ': ci, 'scalar': sc,
                    'passed': abs(ch - ech) <= 1e-12 and abs(ci - eci) <= 1e-12 and abs(sc - esc) <= 1e-12}
    chk['AD3_F1_analytic'] = {'cases': ad3, 'passed': all(v['passed'] for v in ad3.values())}

    # AD4: paired shift invariance on F2 consecutive frame pairs (T0)
    worst = 0.0
    for cid in lists['F2']:
        base = cases[f"{manifest()[cid]['derived_from']}_a+000"]['tors']
        t = cases[cid]['tors']
        for fn in (pair_chord, lambda u, v: np.sqrt(np.mean(wrap(u - v) ** 2, axis=1))):
            worst = max(worst, float(np.abs(fn(t[:-1], t[1:]) - fn(base[:-1], base[1:])).max()))
    chk['AD4_paired_shift'] = {'max_abs_diff': worst, 'passed': worst <= 1e-12}

    # AD5: continuity on F3 (Lipschitz) and ±ε noiseless pairs vs scalar comparator
    lip, pm = 0.0, {}
    lip_ok = True
    for cid in lists['F3']:
        t = cases[cid]['tors'][:1]
        for d in (1e-9, 1e-6):
            tp = t.copy()
            tp[0, 1] += d
            dev = float(np.linalg.norm(ordered_features(tp) - ordered_features(t)))
            lip = max(lip, dev / d)
            lip_ok &= dev <= d + 1e-15
    for e in [1e-1, 1e-3, 1e-6, 1e-9]:
        tp, tn = cases[f'F3_e{e:+.0e}_n0']['tors'], cases[f'F3_e{-e:+.0e}_n0']['tors']
        ch = float(d_chord(tp[:1], tn[:1])[0, 0])
        sc = float(circ(cmp[f'F3_e{e:+.0e}_n0']['means'][0], cmp[f'F3_e{-e:+.0e}_n0']['means'][0]))
        pm[f'{e:.0e}'] = {'D_ch_plus_vs_minus': ch, 'expected': math.sqrt(2) * abs(math.sin(e)),
                          'scalar_comparator_plus_vs_minus': sc,
                          'passed': abs(ch - math.sqrt(2) * abs(math.sin(e))) <= 1e-12}
    chk['AD5_continuity'] = {'max_feature_change_per_unit_delta': lip,
                             'lipschitz_passed': bool(lip_ok), 'plus_minus_eps': pm,
                             'passed': bool(lip_ok) and all(v['passed'] for v in pm.values())}

    # AD6: N2 mirror analytic; Cα-distance baseline cannot separate mirror images
    ta, tm = cases['F1a_s0']['tors'], ctrl['N2_mirror_F1a_s0']['tors']
    ch = pair_chord(ta, tm)
    expd = np.sqrt(np.mean(4 * np.sin(ta) ** 2, axis=1))
    cav, _ = ca_distance_vectors(cases['F1a_s0']['coords'])
    cmv, _ = ca_distance_vectors(ctrl['N2_mirror_F1a_s0']['coords'])
    chk['AD6_mirror'] = {'max_abs_err_D_ch': float(np.abs(ch - expd).max()),
                         'max_ca_distance_diff_A': float(np.abs(cav - cmv).max()),
                         'passed': bool(np.abs(ch - expd).max() <= 1e-12 and np.abs(cav - cmv).max() <= 1e-9)}

    # AD8: F6b equal torsions -> D_ch = 2|sin(Δx/2)| on consecutive frames
    worst = 0.0
    for cid in lists['F6b']:
        t = cases[cid]['tors']
        got = np.sqrt(np.sum((ordered_features(t[:-1]) - ordered_features(t[1:])) ** 2, -1) / 2)
        worst = max(worst, float(np.abs(got - 2 * np.abs(np.sin((t[:-1, 0] - t[1:, 0]) / 2))).max()))
    chk['AD8_F6b'] = {'max_abs_err': worst, 'passed': worst <= 1e-12}

    # comparator integrity: production R2 replication and REF agreement where conditioned
    ref_fail = []
    for cid, v in cmp.items():
        t = all_cases[cid]['tors']
        z = np.mean(np.exp(1j * t), axis=1)
        b = wpa_common.direction_bound(t.shape[1], np.abs(z))
        ok = b <= wpa_common.UNDEFINED_BOUND
        if not v['r2_bits_equal'] or np.any(circ(v['means'], np.angle(z))[ok] > 4 * b[ok] + 1e-12):
            ref_fail.append(cid)
    chk['comparator_integrity'] = {'failures': ref_fail, 'passed': not ref_fail}

    # retrieval (T0 primary, T1 secondary)
    retr_sets = {cid: cid for cid in lists['F1'] + lists['F4a0'] + [k for k in ctrl if k.startswith('N1')]}
    bounds_fail, d_rec_fail = [], []
    for name in list(retr_sets) + ['N2_pooled_F1a_s0_mirror']:
        if name == 'N2_pooled_F1a_s0_mirror':
            c1, c2 = cases['F1a_s0'], ctrl['N2_mirror_F1a_s0']
            tors = np.concatenate([c1['tors'], c2['tors']])
            coords = np.concatenate([c1['coords'], c2['coords']])
            labels = np.concatenate([c1['labels'], c2['labels'] + 2])
            means = np.concatenate([cmp['F1a_s0']['means'], cmp['N2_mirror_F1a_s0']['means']])
            p = {s: tuple(np.concatenate([paths['F1a_s0'][s][i], paths['N2_mirror_F1a_s0'][s][i]])
                          for i in range(2)) for s in ('T1', 'T2')}
        else:
            c = all_cases[name]
            tors, coords, labels, means, p = c['tors'], c['coords'], c['labels'], cmp[name]['means'], paths[name]
        F, W = coords.shape[:2]
        ca0, _ = ca_distance_vectors(coords)
        mats0, s0 = retrieval_block(tors, np.full((F, W - 2), KAPPA), ca0, means, labels)
        if not bounds_ok(mats0['D_ch'], mats0['D_circ']):
            bounds_fail.append(name)
        tau1, kap1 = p['T1']
        ca1, _ = ca_distance_vectors(coords.astype(np.float32).astype(np.float64))
        means1 = np.angle(np.mean(np.exp(1j * tau1), axis=1))     # REF formula on T1 torsions
        mats1, s1 = retrieval_block(tau1, kap1, ca1, means1, labels)
        for src, tau in (('T1', tau1), ('T2', p['T2'][0])):
            delta = float(circ(tau, tors).max())
            Dt = mats1['D_ch'] if src == 'T1' else d_chord(tau, tau)
            if np.abs(Dt - mats0['D_ch']).max() > 2 * delta + 1e-12:
                d_rec_fail.append((name, src))
        R['retrieval'][name] = {'n_frames': int(F), 'n_labels': int(len(np.unique(labels))),
                                'T0': s0, 'T1': s1,
                                'scan_bc_tau_r2': cmp[name]['scan_bc_tau_r2'] if name in cmp else None}
    chk['AD2_bounds_on_retrieval_matrices'] = {'failures': bounds_fail, 'passed': not bounds_fail}
    chk['distance_recovery_aware'] = {'failures': d_rec_fail, 'passed': not d_rec_fail}

    # tie-handling self-check (fully tied set -> chance purity)
    lab = np.array([0, 1, 0, 1, 0, 1])
    t = knn_purity(np.zeros((6, 6)), lab, 2)
    chk['tie_self_check'] = {'purity': t['purity'], 'chance': t['chance'],
                             'passed': abs(t['purity'] - t['chance']) <= 1e-15 and t['tie_frequency'] == 1.0}

    # scan outputs for N1 (paired ±ε): historical component, now with paired draws
    R['n1_scan_bc_tau'] = {k: cmp[k]['scan_bc_tau_r2'] for k in sorted(cmp) if k.startswith('N1')}
    return R


def verdicts(R, coverage_ok, inputs_ok, binding_ok):
    c = R['checks']
    v1 = {'coverage': coverage_ok, 'input_hashes': inputs_ok, 'binding_identity': binding_ok,
          'recovery': c['recovery']['passed'], 'comparator_integrity': c['comparator_integrity']['passed'],
          'tie_self_check': c['tie_self_check']['passed']}
    v2 = {'feature_contracts': all(x['passed'] for x in c['feature_contracts']),
          **{k: c[k]['passed'] for k in ['AD1_chord_identity', 'AD2_bounds_on_retrieval_matrices',
                                         'AD3_F1_analytic', 'AD4_paired_shift', 'AD5_continuity',
                                         'AD6_mirror', 'AD7_branch', 'AD8_F6b',
                                         'distance_recovery_aware']}}
    v3 = {}
    for name, r in R['retrieval'].items():
        ch = r['T0']['D_ch']['k10']['purity']
        sc = r['T0']['scalar_comparator']['k10']['purity']
        ca = r['T0']['ca_distance']['k10']['purity']
        v3[name] = {'D_ch': ch, 'scalar': sc, 'ca_distance': ca,
                    'deficit_vs_scalar': ch < sc - MARGIN, 'deficit_vs_ca': ch < ca - MARGIN,
                    'diff_vs_scalar': ch - sc, 'diff_vs_ca': ch - ca}
    V1, V2 = all(v1.values()), all(v2.values())
    scalar_deficits = sorted(n for n, v in v3.items() if v['deficit_vs_scalar'])
    if not V1:
        decision = 'INSUFFICIENT EVIDENCE'
    elif not V2:
        decision = 'REJECT'
    elif scalar_deficits:
        decision = 'RESTRICT'
    else:
        decision = 'RETAIN: eligible for a real-MD pilot'
    return {'V1_harness': {'passed': V1, 'components': v1},
            'V2_fidelity': {'passed': V2, 'components': v2},
            'V3_retrieval_k10_T0': {'margin': MARGIN, 'fixtures': v3,
                                    'deficits_vs_scalar': scalar_deficits,
                                    'deficits_vs_ca_distance': sorted(n for n, v in v3.items()
                                                                      if v['deficit_vs_ca'])},
            'decision': decision}


def coverage(R, lists, ctrl):
    expected_cases = set(sum(lists.values(), [])) | set(ctrl)
    expected_retr = set(lists['F1']) | set(lists['F4a0']) | {k for k in ctrl if k.startswith('N1')} \
        | {'N2_pooled_F1a_s0_mirror'}
    return set(R['cases']) == expected_cases and set(R['retrieval']) == expected_retr


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true', help='recompute and compare; write nothing')
    ap.add_argument('--results', default=str(RES_JSON))
    args = ap.parse_args(argv)
    prov = {'_meta': {'type': 'SYNTHETIC', 'producer': 'run.py',
                      'started_utc': datetime.datetime.now(datetime.UTC).isoformat()}}
    if not args.verify:
        prov['git'] = {'head': run_cmd(['git', 'rev-parse', 'HEAD']).stdout.strip(),
                       'status_porcelain': run_cmd(['git', 'status', '--porcelain']).stdout.splitlines()}
        prov['versions'] = {'python': sys.version.split()[0], 'numpy': np.__version__,
                            'platform': platform.platform(),
                            'rustc': run_cmd(['rustc', '--version']).stdout.strip()}
    inputs_ok = True
    try:
        lists, cases, sources = gather_inputs()
        ctrl, recs = controls(write=not args.verify)
    except InputHashError as e:
        print('INPUT HASH FAILURE:', e)
        return 1
    tf, ident = binding_identity()
    binding_ok = ident['equals_wpa_freshly_built_extension'] and \
        ident['production_sources_unchanged_since_wpa_build']
    if args.verify:
        stored = json.loads(Path(args.results).read_text())
        cmp, cmp_hashes = load_comparator(stored['inputs']['comparator_output_sha256'])
    else:
        run_comparator(lists, cases, ctrl, recs)
        cmp, cmp_hashes = load_comparator()
    R = compute(tf, lists, cases, ctrl, cmp)
    cov = coverage(R, lists, ctrl)
    R['verdicts'] = verdicts(R, cov, inputs_ok, binding_ok)
    R['inputs'] = {'wpa_manifest_sha256': sha_bytes((WPA / 'fixtures' / 'SYNTHETIC_fixtures.json').read_bytes()),
                   'wpa_input_source': sorted(set(sources.values())),
                   'controls_manifest_case_hashes': {k: v['torsion_sha256'] for k, v in recs.items()},
                   'comparator_output_sha256': cmp_hashes, 'binding': ident,
                   'case_families': {k: v for k, v in lists.items()},
                   'paths_per_case': {'T0/T2/T1 + retrieval': 'see cases and retrieval keys'}}
    R['_meta'] = {'type': 'SYNTHETIC', 'producer': 'run.py',
                  'protocol': 'docs/ORDERED_ANGULAR_REPRESENTATION_PROTOCOL.md (amended §11)'}
    fresh = json.loads(json.dumps(R, default=wpa_common._json_default))
    if args.verify:
        stored = json.loads(Path(args.results).read_text())
        missing = sorted(set(fresh['cases']) - set(stored.get('cases', {})))
        diff = sorted(k for k in fresh if stored.get(k) != fresh.get(k))
        print(f'verify: missing cases {missing}; differing sections {diff}')
        return 0 if not missing and not diff else 1
    wpa_common.dump_json(fresh, RES_JSON)
    prov['inputs_sha256'] = {'results': sha_bytes(RES_JSON.read_bytes()),
                             'controls_manifest': sha_bytes(CTRL_JSON.read_bytes())}
    prov['harness_sha256'] = {p.name: sha_bytes(p.read_bytes()) for p in sorted(HERE.glob('*.py'))}
    prov['commands'] = log
    prov['binding'] = ident
    prov['input_sources'] = sources
    wpa_common.dump_json(prov, PROV_JSON)
    v = fresh['verdicts']
    print(f"V1 harness {v['V1_harness']['passed']}; V2 fidelity {v['V2_fidelity']['passed']}; "
          f"deficits vs scalar {v['V3_retrieval_k10_T0']['deficits_vs_scalar']}; "
          f"deficits vs Cα {v['V3_retrieval_k10_T0']['deficits_vs_ca_distance']}")
    print('DECISION:', v['decision'])
    return 0 if v['V1_harness']['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
