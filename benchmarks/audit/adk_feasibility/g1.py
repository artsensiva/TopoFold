#!/usr/bin/env python3
"""WP-B G1: independent reference-based assessment, without TopoFold (Amendments 1–2 govern).

Third-party MD data (Zenodo 10.5281/zenodo.5583119, CC BY 4.0): not SYNTHETIC. Trajectories are
kept separate. Proximity classes are descriptive and provisional. No populations, rates or free
energies. Run after g0.py passes: build/venv/bin/python -B g1.py
"""
import csv
import datetime
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from adk_lib import (OUTCOME_TEXT, classify_frames, derive_core, final_outcome, fit_rmsd,
                     make_whole, map_reference, max_consecutive, parse_pdb_ca, runs_with_time,
                     starts_near, summarize_runs, tau_int, variant_category)
from g0 import AMD, RAW, RES, RESIDS, SPEC, SPEC_FILES, sha256

HERE = Path(__file__).resolve().parent
DERIVED = HERE / 'data' / 'derived'
OUT = RES / 'adk_g1_summary.json'
START_REF = {'open_start': 'open', 'closed_start': 'closed'}


def fail(summary, msg):
    summary['stopped'] = msg
    OUT.write_text(json.dumps(summary, indent=1) + '\n')
    print('STOP:', msg)
    return 1


def main():
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms
    summary = {'_meta': {'producer': 'g1.py', 'data': 'third-party MD (not SYNTHETIC)',
                         'source_doi': SPEC['source_doi'],
                         'governing_spec': 'SPEC.md + AMENDMENTS.md (Amendments 1–2)',
                         'started_utc': datetime.datetime.now(datetime.UTC).isoformat()}}

    # --- A6: bind to G0 -------------------------------------------------------------------------
    g0p = RES / 'adk_g0_provenance.json'
    g0 = json.loads(g0p.read_text())
    freeze = json.loads((RES / 'adk_spec_amendment2_freeze.json').read_text())
    now_spec = {f: sha256(HERE / f) for f in SPEC_FILES}
    checks = {'g0_passed': g0.get('g0_passed') is True,
              'g0_spec_hashes_equal_current': g0.get('spec_sha256') == now_spec,
              'current_spec_equals_freeze': now_spec == freeze['files'],
              'g0_after_amendment_freeze': g0['_meta'].get('finished_utc', '') > freeze['frozen_utc'],
              'raw_files_equal_g0': all(sha256(RAW / n) == g0['files'][n]['sha256'] for n in SPEC['files']),
              'references_equal_g0': all(sha256(RAW / f'{r}.pdb') == g0['references'][r]['sha256']
                                         for r in ('4AKE', '1AKE'))}
    summary['g0_binding'] = dict(checks, g0_provenance_sha256=sha256(g0p))
    if not all(checks.values()):
        return fail(summary, f'G0 binding failed: {checks}')

    # --- A1: references -------------------------------------------------------------------------
    refs = {}
    for rid in ('4AKE', '1AKE'):
        chains, _, _ = parse_pdb_ca(RAW / f'{rid}.pdb')
        for ch in (AMD['primary_chains'][rid], AMD['sensitivity_chains'][rid]):
            xyz, rep = map_reference(chains.get(ch, {}), AMD['declared_sequence'], RESIDS)
            refs[(rid, ch)] = xyz
    if refs[('4AKE', 'A')] is None or refs[('1AKE', 'A')] is None:
        return fail(summary, 'primary reference mapping invalid')
    chain_b_ok = refs[('4AKE', 'B')] is not None and refs[('1AKE', 'B')] is not None

    # --- A2: cores --------------------------------------------------------------------------------
    cdef = SPEC['core']
    alg = AMD['core_algorithm']
    core_specs = {'A_1.0': ('A', cdef['cutoff_A']), 'A_0.5': ('A', 0.5), 'A_2.0': ('A', 2.0)}
    if chain_b_ok:
        core_specs['B_1.0'] = ('B', cdef['cutoff_A'])
    cores, core_meta = {}, {}
    for name, (ch, cut) in core_specs.items():
        idx, meta = derive_core(refs[('4AKE', ch)], refs[('1AKE', ch)], cut, alg['phase2_max_iter'],
                                cdef['min_size'], alg['phase1_trim_fraction'], alg['phase1_max_steps'])
        meta['resids'] = [RESIDS[i] for i in idx]
        core_meta[name] = meta
        if meta['valid']:
            cores[name] = idx
    if not chain_b_ok:
        core_meta['B_1.0'] = {'valid': False, 'invalid_reason': 'chain B reference mapping invalid'}
    summary['cores'] = core_meta
    if 'A_1.0' not in cores:
        return fail(summary, f"primary core invalid: {core_meta['A_1.0']['invalid_reason']}")
    primary = cores['A_1.0']
    O, C = refs[('4AKE', 'A')], refs[('1AKE', 'A')]
    all_idx = np.arange(len(RESIDS))
    non = sorted(set(all_idx) - set(primary))
    segs, cur = [], []
    for i in non:
        if cur and i != cur[-1] + 1:
            segs.append(cur)
            cur = []
        cur.append(i)
    segs = [s for s in segs + ([cur] if cur else []) if len(s) >= SPEC['segments_min_len']]
    summary['segments_resids'] = [[RESIDS[s[0]], RESIDS[s[-1]]] for s in segs]

    # variants: (open ref, closed ref, core idx, fraction)
    variants = {}
    for frac in AMD['variants']['fractions_primary_core']:
        variants[f'A_core1.0_f{frac}'] = ('A', 'A_1.0', frac)
    for cut in AMD['variants']['core_cutoff_sensitivity_A']:
        key = f'A_{cut}'
        if key in cores:
            variants[f'A_core{cut}_f0.25'] = ('A', key, 0.25)
    if 'B_1.0' in cores:
        variants['B_core1.0_f0.25'] = ('B', 'B_1.0', 0.25)
    summary['variants_used'] = sorted(variants)
    summary['variants_excluded'] = {k: v['invalid_reason'] for k, v in core_meta.items() if not v.get('valid')}
    d_ref = {}
    for name, (ch, ck, _) in variants.items():
        d_ref[(ch, ck)] = fit_rmsd(refs[('1AKE', ch)], refs[('4AKE', ch)], cores[ck], [all_idx])[0]
    summary['d_ref_A'] = {f'{ch}/{ck}': v for (ch, ck), v in d_ref.items()}

    # --- A4: Kabsch cross-check, crystal pair -------------------------------------------------------
    tol = AMD['kabsch_crosscheck']['tolerance_A']
    xc = []
    ours = fit_rmsd(C, O, primary, [primary])[0]
    theirs = float(rms.rmsd(C[primary], O[primary], center=True, superposition=True))
    xc.append({'case': 'crystal 1AKE:A vs 4AKE:A core', 'ours': ours, 'mdanalysis': theirs,
               'abs_diff': abs(ours - theirs)})

    DERIVED.mkdir(parents=True, exist_ok=True)
    trajs, rmsd_cache = {}, {}
    nominal = {lab: g0['trajectories'][lab]['nominal_dt_ps'] for lab in SPEC['pairs']}
    for lab, (gro, xtc) in SPEC['pairs'].items():
        u = mda.Universe(str(RAW / gro), str(RAW / xtc))
        ca = u.select_atoms(SPEC['selection'])
        n = len(u.trajectory)
        T, X = np.empty(n), np.empty((n, len(RESIDS), 3))
        broken, max_after = 0, 0.0
        for k, ts in enumerate(u.trajectory):
            raw = ca.positions.astype(np.float64)
            broken += max_consecutive(raw) > SPEC['pbc']['max_consecutive_ca_A']
            X[k] = make_whole(raw, ts.dimensions.astype(np.float64))
            max_after = max(max_after, max_consecutive(X[k]))
            T[k] = ts.time
        if max_after > SPEC['pbc']['max_consecutive_ca_A']:
            return fail(summary, f'{lab}: chain not whole after reconstruction ({max_after:.3f} Å)')
        np.save(DERIVED / f'{lab}_ca_whole.npy', X)
        for k in range(0, n, AMD['kabsch_crosscheck']['every_nth_frame']):
            for rname, R in (('4AKE:A', O), ('1AKE:A', C)):
                a = fit_rmsd(X[k], R, primary, [primary])[0]
                b = float(rms.rmsd(X[k][primary], R[primary], center=True, superposition=True))
                xc.append({'case': f'{lab} frame {k} vs {rname} core', 'ours': a, 'mdanalysis': b,
                           'abs_diff': abs(a - b)})
        for ch in {v[0] for v in variants.values()}:
            for ck in {v[1] for v in variants.values() if v[0] == ch}:
                ro = np.array([fit_rmsd(x, refs[('4AKE', ch)], cores[ck], [all_idx])[0] for x in X])
                rc = np.array([fit_rmsd(x, refs[('1AKE', ch)], cores[ck], [all_idx])[0] for x in X])
                rmsd_cache[(lab, ch, ck)] = (ro, rc)
        prim = {'time_ps': T}
        for name, R in (('open', O), ('closed', C)):
            vals = np.array([fit_rmsd(x, R, primary, [all_idx, primary] + [np.array(s) for s in segs]) for x in X])
            prim[f'rmsd_{name}'], prim[f'core_rmsd_{name}'] = vals[:, 0], vals[:, 1]
            for j, s in enumerate(segs):
                prim[f'seg{RESIDS[s[0]]}_{RESIDS[s[-1]]}_rmsd_{name}'] = vals[:, 2 + j]
        prim['delta'] = prim['rmsd_closed'] - prim['rmsd_open']
        with open(RES / f'adk_g1_frames_{lab.split("_")[0]}.csv', 'w', newline='') as fh:
            w = csv.writer(fh)
            w.writerow(list(prim))
            for row in zip(*prim.values()):
                w.writerow([f'{v:.6f}' for v in row])
        stats = {}
        for key in ('rmsd_open', 'rmsd_closed', 'delta', 'core_rmsd_open', 'core_rmsd_closed'):
            v = prim[key]
            m = max(1, n // 10)
            ti = tau_int(v)
            stats[key] = {'min': float(v.min()), 'q05': float(np.quantile(v, .05)),
                          'median': float(np.median(v)), 'q95': float(np.quantile(v, .95)),
                          'max': float(v.max()), 'first10pct_mean': float(v[:m].mean()),
                          'last10pct_mean': float(v[-m:].mean()), 'tau_int_frames_est': ti,
                          'ess_est': (n / ti) if ti else None}
        trajs[lab] = {'n_frames': n, 'times': T, 'frames_chain_split_in_raw': int(broken),
                      'max_consecutive_ca_after_whole_A': max_after, 'observables': stats,
                      'derived_ca_sha256': sha256(DERIVED / f'{lab}_ca_whole.npy'),
                      'frames_csv_sha256': sha256(RES / f'adk_g1_frames_{lab.split("_")[0]}.csv')}
        print(f'{lab}: {n} frames processed')

    summary['kabsch_crosscheck'] = {'tolerance_A': tol, 'n_comparisons': len(xc),
                                    'max_abs_diff': max(x['abs_diff'] for x in xc),
                                    'passed': all(x['abs_diff'] <= tol for x in xc),
                                    'worst': max(xc, key=lambda x: x['abs_diff'])}
    if not summary['kabsch_crosscheck']['passed']:
        return fail(summary, 'Kabsch cross-check exceeded the predeclared tolerance')

    # --- A3/A5: classification per variant -------------------------------------------------------
    sustained_ps = AMD['sustained_span_ns'] * 1000.0
    window_ps = AMD['start_window_ns'] * 1000.0
    per_variant, categories = {}, {}
    for vname, (ch, ck, frac) in variants.items():
        radius = frac * d_ref[(ch, ck)]
        vt = {}
        for lab in SPEC['pairs']:
            ro, rc = rmsd_cache[(lab, ch, ck)]
            labels = classify_frames(ro, rc, radius)
            runs = runs_with_time(labels, trajs[lab]['times'], nominal[lab], AMD['cadence_tolerance_ps'])
            vt[lab] = {'summary': summarize_runs(runs, sustained_ps), 'n_frames': trajs[lab]['n_frames'],
                       'start_ref': START_REF[lab],
                       'starts_near_own_reference': starts_near(labels, trajs[lab]['times'],
                                                                START_REF[lab], window_ps),
                       'frame0_class': str(labels[0])}
        categories[vname] = variant_category(vt, AMD['ambiguity_dominance_fraction'])
        per_variant[vname] = {'radius_A': radius, 'category': categories[vname], 'trajectories': vt}
    out, agree = final_outcome(categories)
    other = {'open': 'closed', 'closed': 'open'}
    evidence = {}
    for lab in SPEC['pairs']:
        obs = {v: per_variant[v]['trajectories'][lab]['summary'][other[START_REF[lab]]]['sustained_runs'] > 0
               for v in variants}
        claimed = all(per_variant[v]['trajectories'][lab]['starts_near_own_reference'] and obs[v]
                      for v in variants)
        evidence[lab] = {'observed_cross_visit_by_variant': obs, 'transition_evidence_claimed': claimed}
    summary['classification'] = per_variant
    summary['cross_visits_and_transition_evidence'] = evidence
    summary['outcome'] = {'categories_by_variant': categories, 'variants_agree': agree,
                          'final': OUTCOME_TEXT[out]}
    for lab in trajs:
        trajs[lab].pop('times')
    summary['trajectories'] = trajs
    summary['_meta']['finished_utc'] = datetime.datetime.now(datetime.UTC).isoformat()
    OUT.write_text(json.dumps(summary, indent=1) + '\n')
    print('OUTCOME:', summary['outcome'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
