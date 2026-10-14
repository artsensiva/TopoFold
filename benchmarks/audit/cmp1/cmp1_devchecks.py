#!/usr/bin/env python3
"""CMP-1 adapter-development checks (NOT deciding; separate queries, sizes and seeds from the frozen workloads).

Verifies documented semantics of the pinned libraries, input immutability, unit-conversion/representation error,
accuracy-track behaviour and per-pair cost for budget planning. Writes results/cmp1_adapter_dev_checks.json.
Run: build/venv/bin/python -B cmp1_devchecks.py
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

import cmp1_accuracy as ACC  # noqa: E402
import cmp1_adapters as AD  # noqa: E402
import cmp1_data as DATA  # noqa: E402

HERE = Path(__file__).resolve().parent
DEV = [{'id': 'dev_D1', 'kind': 'real', 'd': 'D1', 'nd': 200, 'nwin': 20, 'queries': [557, 688, 720]},
       {'id': 'dev_D2', 'kind': 'real', 'd': 'D2', 'nd': 100, 'nwin': 10, 'queries': [263, 315, 371]},
       {'id': 'dev_tiled', 'kind': 'SYNTHETIC_tiled', 'd': 'D1p', 'nd': 300, 'nwin': 50, 'queries': [7, 8, 9]}]


def h(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def main():
    out = {'note': 'adapter-development checks; not deciding; separate from frozen workloads', 'cases': {}}
    for w in DEV:
        X, di, wins, info = DATA.load(w)
        h0 = h(X)
        c = {'canonical_input': info}
        Dr = ACC.ref_distances(X, di, wins, 0)
        Gq = [ACC.window_G(X, [0], s)[0] for s in wins]
        Gd = [ACC.window_G(X, di, s) for s in wins]
        for name, cls in AD.IMPLS.items():
            impl = cls()
            t0 = time.perf_counter()
            hd = impl.prepare(X, di, wins)
            tp = time.perf_counter() - t0
            t0 = time.perf_counter()
            impl.query(hd, X, 0)
            tq = time.perf_counter() - t0
            D = impl.distances(hd, X, 0)
            TI = np.stack([AD.B.topk_lex(D[i]) for i in range(len(wins))])
            c[name] = {t: ACC.evaluate(t, D, Dr, TI, Gq, Gd) for t in ('F64', 'F32')}
            c[name]['dev_timing_s'] = {'prepare': tp, 'query': tq,
                                       'query_per_pair_us': 1e6 * tq / (len(di) * len(wins))}
            c[name]['input_unchanged'] = h(X) == h0
        # MDTraj representation error: reference on the converted representation vs canonical reference
        mdt = AD.MDT()
        Xr = mdt.to_nm32(X).astype(np.float64) * 10.0
        c['E-MDT_representation_error_A'] = float(np.max(np.abs(ACC.ref_distances(Xr, di, wins, 0) - Dr)))
        out['cases'][w['id']] = c
    # MDTraj semantics (pinned 1.11.2): precentered ignored with atom_indices; pre-sliced precentered valid
    import mdtraj as md
    X, di, wins, _ = DATA.load(DEV[0])
    mdt = AD.MDT()
    hd = mdt.prepare(X, di, wins)
    top = md.Topology()
    ch = top.add_chain()
    for _ in range(X.shape[1]):
        top.add_atom('CA', md.element.carbon, top.add_residue('ALA', ch))
    tr = md.Trajectory(mdt.to_nm32(X[np.asarray(di)]), top)
    ref = md.Trajectory(mdt.to_nm32(X[[0]]), top)
    s = wins[3]
    idx = np.arange(s, s + 8)
    tr.center_coordinates()
    ref.center_coordinates()
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter('always')
        d_sub_pre = md.rmsd(tr, ref, 0, atom_indices=idx, precentered=True, parallel=False)
    tr2 = md.Trajectory(mdt.to_nm32(X[np.asarray(di)]), top)
    xyz_before = tr2.xyz.copy()
    d_sub = md.rmsd(tr2, ref, 0, atom_indices=idx, parallel=False)
    d_slice = mdt._dist(hd, X, 0, 3)[0] / 10.0
    tr3 = md.Trajectory(mdt.to_nm32(X[np.asarray(di)]), top)
    b3 = tr3.xyz.copy()
    md.rmsd(tr3, ref, 0, parallel=False)                         # atom_indices=None, not precentered
    out['mdtraj_semantics'] = {
        'precentered_with_atom_indices_emits_warning': any('precentered is ignored' in str(x.message) for x in wl),
        'precentered_with_atom_indices_vs_unprecentered_max_abs_nm': float(np.max(np.abs(d_sub_pre - d_sub))),
        'atom_indices_path_leaves_trajectory_unchanged': bool(np.array_equal(xyz_before, tr2.xyz)),
        'atom_indices_None_path_centres_input_in_place': bool(not np.array_equal(b3, tr3.xyz)),
        'pre_sliced_precentered_vs_subset_max_abs_nm': float(np.max(np.abs(d_slice - d_sub))),
        'distances_dtype': str(d_sub.dtype)}
    out['degeneracy_suite_dev_seed'] = {n: ACC.degeneracy_suite(c(), seed=777) for n, c in AD.IMPLS.items()}
    (HERE / 'results').mkdir(exist_ok=True)
    (HERE / 'results' / 'cmp1_adapter_dev_checks.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    for cid, c in out['cases'].items():
        print(cid, {n: (c[n]['F64']['pass'], c[n]['F32']['pass'], f"{c[n]['F64']['max_abs_err_A']:.2e}", round(c[n]['F64']['native_msd_ratio_descriptive'], 1), round(c[n]['dev_timing_s']['query_per_pair_us'], 3),
                        c[n]['input_unchanged']) for n in AD.IMPLS}, 'repr_err', c['E-MDT_representation_error_A'])
    print(out['mdtraj_semantics'])
    print({n: (v['F64']['pass'], v['F32']['pass'], {k: x for k, x in v['F32'].items() if x is False},
               {k: x for k, x in v['F32'].items() if k.endswith('_err_A')})
           for n, v in out['degeneracy_suite_dev_seed'].items()})


if __name__ == '__main__':
    main()
