#!/usr/bin/env python3
"""Compare R1/R2/P1/REF outputs (protocol §4, §6) and write the SYNTHETIC summary.

Harness-correctness checks fail the run (exit 1). Sensitivity and information-loss findings are
recorded, never failed. Conditioning rules are those of protocol §4, fixed before execution.
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np

from common import (FIXTURE_BIN, FIXTURES_JSON, GUARD_VAR, HERE, P1_RECOVERY_FLAG, PY_JSON,
                    R2_RECOVERY_TOL, REL_PERTURB_MAX, RUST_JSON, SUMMARY_JSON, U,
                    UNDEFINED_BOUND, circ_dist, direction_bound, dump_json, finite_or_none,
                    load_json, sha256_file)

DEG = math.pi / 180
WRAP_IN = 8 * U * math.pi            # input rounding from wrapping constructed angles
DOCUMENTED_F6 = {'F6a_rare_p0.005': 0.744, 'F6a_rare_p0.02': 0.839,
                 'F6a_exponential': 0.577, 'F6a_chi2_df1': 0.590}   # KNOWN_ISSUES M3


class Run:
    def __init__(self):
        self.checks = []

    def check(self, name, ok, detail=None):
        self.checks.append({'name': name, 'passed': bool(ok), 'detail': detail})
        return ok


def arr(store, ref):
    off, n = ref
    return store[off:off + n]


def nn(x):
    return None if x is None else finite_or_none(float(x))


def bc_comparison(means_a, R_ref, means_ref, devs_ref, var_a, var_ref, bc_a, bc_ref, ncols, F,
                  delta_in):
    """Protocol §4 conditioning diagnostic for one implementation vs REF."""
    b_frame = direction_bound(ncols, R_ref, delta_in)
    undefined_frames = int(np.sum(b_frame > UNDEFINED_BOUND))
    defined = b_frame <= UNDEFINED_BOUND
    d_frame = circ_dist(means_a, means_ref)
    frame_fail = int(np.sum(d_frame[defined] > 4 * b_frame[defined] + 1e-12))
    out = {'undefined_frames': undefined_frames, 'frame_mean_failures': frame_fail,
           'max_frame_circ_diff_defined': nn(d_frame[defined].max()) if defined.any() else None,
           'max_frame_circ_diff_undefined': nn(d_frame[~defined].max()) if (~defined).any()
           else None}
    _, R_across = np.angle(np.mean(np.exp(1j * means_ref))), np.abs(np.mean(np.exp(1j * means_ref)))
    max_b = float(b_frame.max()) if undefined_frames == 0 else math.inf
    b_across = float(direction_bound(F, R_across, max_b if math.isfinite(max_b) else 0.0)) \
        if undefined_frames == 0 else math.inf
    out['b_across'] = nn(b_across)
    out['across_R_ref'] = float(R_across)
    s = math.sqrt(var_ref) if var_ref > 0 else 0.0
    e_dev = (max_b + b_across) if math.isfinite(b_across) else math.inf
    near_pi = bool(np.any(np.abs(devs_ref) > math.pi - e_dev)) if math.isfinite(e_dev) else True
    r = e_dev / s if (s > 0 and math.isfinite(e_dev)) else math.inf
    guard_ambiguous = (var_a <= GUARD_VAR) != (var_ref <= GUARD_VAR)
    conditioned = (math.isfinite(b_across) and b_across <= UNDEFINED_BOUND and not near_pi
                   and r <= REL_PERTURB_MAX and not guard_ambiguous)
    out.update(e_dev=nn(e_dev), rel_perturbation=nn(r), near_pi_ambiguity=near_pi,
               guard_ambiguous=guard_ambiguous, conditioned=conditioned,
               bc_a=bc_a, bc_ref=bc_ref, bc_abs_diff=abs(bc_a - bc_ref))
    if conditioned:
        z = float(np.max(np.abs(devs_ref - devs_ref.mean())) / s)
        tol = 10 * z ** 3 * r + 1e-12
        out.update(bc_tol=tol, bc_within_tol=abs(bc_a - bc_ref) <= tol)
    else:
        out.update(bc_tol=None, bc_within_tol=None)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--rust', default=str(RUST_JSON))
    ap.add_argument('--python', default=str(PY_JSON))
    ap.add_argument('--out', default=str(SUMMARY_JSON))
    ap.add_argument('--fixtures-json', default=str(FIXTURES_JSON))
    ap.add_argument('--fixture-dir', default=str(FIXTURE_BIN))
    args = ap.parse_args(argv)
    fixture_dir = Path(args.fixture_dir)

    run = Run()
    man = {c['id']: c for c in load_json(args.fixtures_json)['cases']}
    rust_doc, py_doc = load_json(args.rust), load_json(args.python)
    rust, py = rust_doc['cases'], py_doc['cases']
    rust_arr = np.fromfile(HERE.parents[2] / rust_doc['_meta']['arrays_file'], dtype='<f8')
    ref_arr = np.fromfile(HERE / py_doc['_meta']['arrays_file'], dtype='<f8')
    run.check('ref arrays hash matches python results',
              sha256_file(HERE / py_doc['_meta']['arrays_file']) == py_doc['_meta']['arrays_sha256'])

    # --- coverage and fixture integrity -------------------------------------------------------
    ids = set(man)
    run.check('coverage: rust results cover all fixtures', set(rust) == ids,
              sorted(ids ^ set(rust))[:5])
    run.check('coverage: python results cover all fixtures', set(py) == ids,
              sorted(ids ^ set(py))[:5])
    r2_ids = {i for i, c in man.items() if c['r2_p1']}
    run.check('coverage: R2 present exactly for coordinate fixtures',
              {i for i in rust if 'r2' in rust[i]} == r2_ids)
    run.check('coverage: P1 present exactly for coordinate fixtures',
              {i for i in py if 'p1' in py[i]} == r2_ids)
    bad_hash = [i for i, c in man.items()
                if sha256_file(fixture_dir / c['torsion_file']) != c['torsion_sha256']
                or (c['coords_file'] and
                    sha256_file(fixture_dir / c['coords_file']) != c['coords_sha256'])]
    run.check('fixture files match manifest hashes', not bad_hash, bad_hash[:5])

    # --- construction recovery and replication ------------------------------------------------
    rec_r2 = {i: rust[i]['r2']['recovery'] for i in r2_ids}
    worst_r2 = max(max(v['max_tau_circ_err'], v['max_kappa_err']) for v in rec_r2.values())
    run.check(f'R2 construction recovery <= {R2_RECOVERY_TOL} rad', worst_r2 <= R2_RECOVERY_TOL,
              worst_r2)
    rec_p1 = {i: py[i]['p1']['recovery'] for i in r2_ids}
    worst_p1 = max(max(v['max_tau_circ_err'], v['max_kappa_err']) for v in rec_p1.values())
    p1_flags = sorted(i for i, v in rec_p1.items()
                      if max(v['max_tau_circ_err'], v['max_kappa_err']) > P1_RECOVERY_FLAG)
    bit_fail = [i for i in r2_ids if not (rust[i]['r2']['replication']['bc_tau_bits_equal'] and
                                          rust[i]['r2']['replication']['bc_kappa_bits_equal'])]
    run.check('R2 replication bitwise equal to compute_bimodality_profile', not bit_fail,
              bit_fail[:5])
    sar_fail = [i for i in ids if not rust[i]['r1']['moments']['sarles_fn_bits_equal']]
    run.check('R1 sarles_bimodality_coefficient equals StreamingMoments BC bitwise',
              not sar_fail, sar_fail[:5])

    # --- R1/R2 vs REF under the conditioning diagnostic ---------------------------------------
    per_case = {}
    cmp_fail = []
    for i, c in man.items():
        if c['kind'] == 'scalar':
            a, b = rust[i]['r1']['moments'], py[i]['ref']['moments']
            d = abs(a['bc'] - b['bc'])
            per_case[i] = {'r1_vs_ref': {'bc_abs_diff': d}}
            if d > 1e-12:
                cmp_fail.append((i, 'scalar', d))
            continue
        ref = py[i]['ref']
        means_ref, R_ref, devs_ref = (arr(ref_arr, ref['means']), arr(ref_arr, ref['R']),
                                      arr(ref_arr, ref['deviations']))
        entry = {}
        for path in ['r1', 'r2']:
            if path not in rust[i]:
                continue
            p = rust[i][path]
            mom = p['moments'] if path == 'r1' else p['tau_moments']
            delta_in = 0.0 if path == 'r1' else p['recovery']['max_tau_circ_err']
            res = bc_comparison(arr(rust_arr, p['means']), R_ref, means_ref, devs_ref,
                                mom['variance'] or 0.0, ref['moments']['variance'] or 0.0,
                                mom['bc'] or 0.0, ref['moments']['bc'] or 0.0, c['ncols'], c['F'],
                                delta_in)
            entry[f'{path}_vs_ref'] = res
            if res['frame_mean_failures'] or res['bc_within_tol'] is False:
                cmp_fail.append((i, path, res['frame_mean_failures'], res['bc_abs_diff']))
        per_case[i] = entry
    run.check('R1/R2 vs REF agree wherever conditioned (protocol §4)', not cmp_fail, cmp_fail[:5])
    n_cond = sum(1 for e in per_case.values() for k, v in e.items()
                 if k.endswith('_vs_ref') and v.get('conditioned'))
    n_uncond = sum(1 for e in per_case.values() for k, v in e.items()
                   if k.endswith('_vs_ref') and v.get('conditioned') is False)

    # --- analytical controls (noiseless constructions) ----------------------------------------
    an_fail = []

    def frame_means(i, path):
        src, store = (py[i]['ref'], ref_arr) if path == 'ref' else (rust[i]['r1'], rust_arr)
        return arr(store, src['means']), arr(store, src['R'])

    f1_expect = {'F1a_s0': ([60.0], [math.cos(30 * DEG)]),
                 'F1b_s0': ([60.0], [(2 * math.cos(30 * DEG) + 3) / 5]),
                 'F1c_s0': ([30.0, 90.0], [1.0])}
    for i, (mdeg, Rv) in f1_expect.items():
        for path in ['ref', 'r1']:
            m, R = frame_means(i, path)
            exp_m = np.array([mdeg[f % len(mdeg)] * DEG for f in range(len(m))])
            if circ_dist(m, exp_m).max() > 1e-12 or np.abs(R - Rv[0]).max() > 1e-12:
                an_fail.append((i, path, float(circ_dist(m, exp_m).max())))
        dist = py[i]['ref']['states']['embedding_state_mean_distance']
        if abs(dist - math.sqrt(2)) > 1e-12:
            an_fail.append((i, 'embedding', dist))
    for i, c in man.items():
        if c['family'] == 'F3' and c['params']['sigma_n_rad'] == 0:
            eps, th0 = c['params']['eps_rad'], c['params']['theta0_rad']
            for path in ['ref', 'r1']:
                m, R = frame_means(i, path)
                if np.abs(R - abs(math.sin(eps / 2))).max() > 1e-12:
                    an_fail.append((i, path, 'R'))
                if eps != 0:
                    exp_m = th0 + (math.pi - eps) / 2 + (math.pi if eps < 0 else 0.0)
                    b = float(direction_bound(2, abs(math.sin(eps / 2)), WRAP_IN))
                    if b <= UNDEFINED_BOUND and circ_dist(m, exp_m).max() > 4 * b + 1e-12:
                        an_fail.append((i, path, 'mean'))
        if c['family'] == 'F4' and c['params']['sigma_rad'] == 0 and c['alpha_deg'] == 0:
            p_eff = c['params']['n_phi2'] / c['F']
            exp_R = abs((1 - p_eff) + p_eff * np.exp(1j * c['params']['phi2_rad']))
            for path, val in [('ref', py[i]['ref']['across_R']), ('r1', rust[i]['r1']['across_R'])]:
                if abs(val - exp_R) > 1e-12:
                    an_fail.append((i, path, val, exp_R))
    run.check('analytical controls (F1, F3, F4 noiseless) reproduced by REF and R1', not an_fail,
              an_fail[:5])

    f6_obs = {i: round(py[i]['ref']['moments']['bc'], 3) for i in DOCUMENTED_F6}
    f6_ok = all(f6_obs[i] == v for i, v in DOCUMENTED_F6.items())
    run.check('F6 replay consistent with values documented in KNOWN_ISSUES M3', f6_ok,
              {'observed_3dp': f6_obs, 'documented': DOCUMENTED_F6})

    # F6b: coordinate path vs scalar path on the same values (shift invariance, no wrap)
    f6b = {}
    for i, c in man.items():
        if c['family'] != 'F6b':
            continue
        src = c['derived_from']
        devs = arr(ref_arr, py[i]['ref']['deviations'])
        cmpd = per_case[i]['r2_vs_ref']
        f6b[i] = {'max_abs_deviation': float(np.abs(devs).max()),
                  'bc_r2': rust[i]['r2']['profile']['bc_tau'],
                  'bc_scalar_ref': py[src]['ref']['moments']['bc'],
                  'bc_p1': py[i]['p1']['profile']['bc_tau'], 'r2_vs_ref_tol': cmpd['bc_tol']}
        f6b[i]['r2_minus_scalar'] = f6b[i]['bc_r2'] - f6b[i]['bc_scalar_ref']
    f6b_fail = [i for i, v in f6b.items()
                if v['max_abs_deviation'] >= math.pi or v['r2_vs_ref_tol'] is None
                or abs(v['r2_minus_scalar']) > v['r2_vs_ref_tol']]
    run.check('F6b coordinate path equals scalar path within conditioned tolerance (no wraps)',
              not f6b_fail, f6b_fail)

    # --- findings (recorded, not failed) --------------------------------------------------------
    def bcs(i):
        out = {'r1': rust[i]['r1']['moments']['bc'], 'ref': py[i]['ref']['moments']['bc'],
               'r1_variance': rust[i]['r1']['moments']['variance']}
        if 'r2' in rust[i]:
            out.update(r2=rust[i]['r2']['profile']['bc_tau'], p1=py[i]['p1']['profile']['bc_tau'],
                       r2_variance=rust[i]['r2']['tau_moments']['variance'],
                       r2_guard_triggered=(rust[i]['r2']['tau_moments']['variance'] or 0)
                       <= GUARD_VAR)
        return out

    findings = {}
    findings['F1'] = {i: dict(bcs(i), **py[i]['ref']['states'])
                      for i, c in man.items() if c['family'] == 'F1'}

    f2 = {}
    for i, c in man.items():
        if c['family'] != 'F2':
            continue
        base = f"{c['derived_from']}_a+000"
        d0 = arr(rust_arr, rust[base]['r1']['deviations'])
        da = arr(rust_arr, rust[i]['r1']['deviations'])
        diff = np.abs(da - d0)
        tol0 = per_case[base]['r1_vs_ref']['bc_tol']
        dbc = rust[i]['r1']['moments']['bc'] - rust[base]['r1']['moments']['bc']
        f2[i] = dict(bcs(i), alpha_deg=c['alpha_deg'],
                     r1_bc_minus_alpha0=dbc,
                     r2_bc_minus_alpha0=rust[i]['r2']['profile']['bc_tau']
                     - rust[base]['r2']['profile']['bc_tau'],
                     p1_bc_minus_alpha0=py[i]['p1']['profile']['bc_tau']
                     - py[base]['p1']['profile']['bc_tau'],
                     r1_max_abs_dev_change_nonwrap=float(diff[diff <= math.pi].max()),
                     r1_wrap_events=int(np.sum(diff > math.pi)),
                     r1_devs_near_pi=int(np.sum(np.abs(da) > math.pi - 1e-6)),
                     alpha0_conditioned_tol=tol0,
                     r1_within_alpha0_tol=(abs(dbc) <= tol0) if tol0 is not None else None)
    findings['F2'] = f2

    f3 = {}
    for i, c in man.items():
        if c['family'] != 'F3':
            continue
        ref, r1 = py[i]['ref'], rust[i]['r1']
        R = arr(ref_arr, ref['R'])
        b = direction_bound(2, R)
        m0 = r1['perturbation_frame0']
        base0 = arr(rust_arr, r1['means'])[0]
        f3[i] = dict(bcs(i), eps=c['params']['eps_rad'], sigma_n=c['params']['sigma_n_rad'],
                     frame_R_min=float(R.min()), frame_R_max=float(R.max()),
                     frames_undefined=int(np.sum(b > UNDEFINED_BOUND)),
                     across_R_r1=r1['across_R'],
                     perturbation_frame0_shift={f"{p['delta']:+.0e}": float(
                         circ_dist(p['frame0_mean'], base0)) for p in m0},
                     r1_vs_ref=per_case[i]['r1_vs_ref'])
    findings['F3'] = f3

    f4 = {}
    for i, c in man.items():
        if c['family'] != 'F4':
            continue
        r1 = rust[i]['r1']
        d = arr(rust_arr, r1['deviations'])
        vals, counts = np.unique(d, return_counts=True)
        top = np.argsort(-counts)[:4]
        f4[i] = dict(bcs(i), base=c['derived_from'], alpha_deg=c['alpha_deg'], p=c['params']['p'],
                     eps=c['params']['eps_rad'], sigma=c['params']['sigma_rad'],
                     across_R_r1=r1['across_R'], across_mean_r1=r1['across_mean'],
                     b_across=nn(direction_bound(c['F'], r1['across_R'])),
                     distinct_deviations=int(vals.size),
                     top_deviation_values=[[float(vals[k]), int(counts[k])] for k in top]
                     if c['params']['sigma_rad'] == 0 else None)
    findings['F4'] = f4
    by_base = {}
    for i, v in f4.items():
        by_base.setdefault(v['base'], []).append(v['r1'])
    findings['F4_rotation_bc_range'] = {k: {'min': min(v), 'max': max(v), 'range': max(v) - min(v)}
                                        for k, v in by_base.items()}

    findings['F6'] = {i: {'ref': py[i]['ref']['moments'], 'r1_bc': rust[i]['r1']['moments']['bc']}
                      for i, c in man.items() if c['family'] == 'F6a'}
    findings['F6b'] = f6b

    kap = {}
    for i in r2_ids:
        kap[i] = {'r2_kappa_variance': rust[i]['r2']['kappa_variance'],
                  'r2_bc_kappa': rust[i]['r2']['profile']['bc_kappa'],
                  'p1_bc_kappa': py[i]['p1']['profile']['bc_kappa'],
                  'p1_kappa_mean_variance_ref': py[i]['p1']['kappa_mean_variance_ref'],
                  'p1_score': py[i]['p1']['profile']['score'],
                  'p1_bc_tau': py[i]['p1']['profile']['bc_tau']}
    findings['kappa_channel'] = {
        'cases': kap,
        'r2_cases_bc_kappa_nonzero': sum(v['r2_bc_kappa'] > 0 for v in kap.values()),
        'p1_cases_bc_kappa_nonzero': sum(v['p1_bc_kappa'] > 0 for v in kap.values()),
        'p1_cases_score_from_kappa': sum(v['p1_bc_kappa'] > v['p1_bc_tau'] for v in kap.values()),
        'n_cases': len(kap)}
    p1d = np.array([py[i]['p1']['profile']['bc_tau'] - rust[i]['r2']['profile']['bc_tau']
                    for i in sorted(r2_ids)])
    findings['p1_minus_r2_bc_tau'] = {'max_abs': float(np.abs(p1d).max()),
                                      'median_abs': float(np.median(np.abs(p1d))),
                                      'n_exactly_equal': int(np.sum(p1d == 0)), 'n': int(p1d.size)}

    harness_pass = all(c['passed'] for c in run.checks)
    dump_json({'_meta': {'type': 'SYNTHETIC', 'producer': 'compare.py',
                         'inputs': {'fixtures': sha256_file(Path(args.fixtures_json)),
                                    'rust_results': sha256_file(Path(args.rust)),
                                    'python_results': sha256_file(Path(args.python))},
                         'conditioning': {'unit_roundoff': U, 'undefined_bound_rad': UNDEFINED_BOUND,
                                          'rel_perturbation_max': REL_PERTURB_MAX,
                                          'r2_recovery_tol': R2_RECOVERY_TOL,
                                          'p1_recovery_flag': P1_RECOVERY_FLAG}},
               'harness': {'passed': harness_pass, 'checks': run.checks,
                           'r2_worst_recovery_rad': worst_r2, 'p1_worst_recovery_rad': worst_p1,
                           'p1_recovery_flagged_cases': p1_flags,
                           'comparisons_conditioned': n_cond,
                           'comparisons_not_conditioned': n_uncond},
               'per_case_comparison': per_case,
               'findings': findings}, Path(args.out))
    for c in run.checks:
        print(f"{'PASS' if c['passed'] else 'FAIL'}  {c['name']}"
              + ('' if c['passed'] else f"  {c['detail']}"))
    print(f"conditioned comparisons: {n_cond}; not conditioned (recorded only): {n_uncond}")
    print(f"HARNESS {'PASS' if harness_pass else 'FAIL'} -> {args.out}")
    return 0 if harness_pass else 1


if __name__ == '__main__':
    sys.exit(main())
