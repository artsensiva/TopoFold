#!/usr/bin/env python3
"""P1 (Python binding, float32 coordinates) and REF (independent NumPy) paths. SYNTHETIC fixtures.

Run with the interpreter of the isolated environment into which the freshly built wheel was
installed (see run_wpa.py); the loaded extension path and hash are recorded.
"""
import hashlib
import platform
import sys
from pathlib import Path

import numpy as np

from common import (ARRAY_BIN, FIXTURE_BIN, FIXTURES_JSON, KAPPA, PY_JSON, circ_dist, dump_json,
                    finite_or_none, load_json, sha256_file)
from reference import angular_pipeline, bc_from_values, mean_and_resultant

PERTURBATIONS = [1e-9, -1e-9, 1e-6, -1e-6]


class Arrays:
    def __init__(self):
        self.chunks, self.n = [], 0

    def push(self, v):
        v = np.asarray(v, dtype=np.float64).ravel()
        ref = [self.n, int(v.size)]
        self.chunks.append(v)
        self.n += v.size
        return ref


def f(x):
    return finite_or_none(float(x))


def bc_json(b):
    return {k: (f(v) if isinstance(v, (float, np.floating)) else v) for k, v in b.items()}


def ref_case(case, tors, arrays):
    if case['kind'] == 'scalar':
        return {'moments': bc_json(bc_from_values(tors[:, 0]))}
    r = angular_pipeline(tors)
    out = {'across_mean': f(r['across_mean']), 'across_R': f(r['across_R']),
           'moments': bc_json(r['bc']), 'means': arrays.push(r['means']),
           'R': arrays.push(r['R']), 'deviations': arrays.push(r['deviations'])}
    if case['perturb']:
        pert = []
        for d in PERTURBATIONS:
            row = tors[0].copy()
            row[1] += d
            pert.append({'delta': d, 'frame0_mean': f(mean_and_resultant(row)[0])})
        out['perturbation_frame0'] = pert
    if case['family'] == 'F1':
        labels = np.fromfile(FIXTURE_BIN / case['labels_file'], dtype='<f8').astype(int)
        emb = np.concatenate([np.cos(tors), np.sin(tors)], axis=1)      # (cos τ_j…, sin τ_j…)
        ea, eb = emb[labels == 0].mean(0), emb[labels == 1].mean(0)
        sa = mean_and_resultant(r['means'][labels == 0])[0]
        sb = mean_and_resultant(r['means'][labels == 1])[0]
        out['states'] = {'embedding_state_mean_distance': f(np.linalg.norm(ea - eb)),
                         'scalar_state_mean_A': f(sa), 'scalar_state_mean_B': f(sb),
                         'scalar_state_mean_circ_dist': f(circ_dist(sa, sb))}
    return out


def p1_case(case, tors, tf):
    W, F = case['W'], case['F']
    X = np.fromfile(FIXTURE_BIN / case['coords_file'], dtype='<f8').reshape(F, W, 3)
    X32 = X.astype(np.float32)
    scores, bc_tau, bc_kappa = tf.compute_bimodality_profile(X32, window_size=W)
    assert len(bc_tau) == 1
    tau_err, kap_err, kap_means = 0.0, 0.0, np.empty(F)
    for i in range(F):
        kap, tau, _ = tf.compute_invariants(X32[i])
        tau_err = max(tau_err, float(circ_dist(np.asarray(tau)[:W - 3], tors[i]).max()))
        kap = np.asarray(kap, dtype=np.float64)
        kap_err = max(kap_err, float(np.abs(kap - KAPPA).max()))
        kap_means[i] = kap[:W - 2].mean()
    return {'profile': {'bc_tau': f(bc_tau[0]), 'bc_kappa': f(bc_kappa[0]), 'score': f(scores[0])},
            'recovery': {'max_tau_circ_err': tau_err, 'max_kappa_err': kap_err},
            'kappa_mean_variance_ref': f(np.var(kap_means))}


def module_identity(tf):
    import topofold.topofold as ext
    so = Path(ext.__file__)
    return {'package_file': str(Path(tf.__file__).resolve()), 'extension_file': str(so.resolve()),
            'extension_sha256': sha256_file(so)}


def main():
    import topofold as tf
    manifest = load_json(FIXTURES_JSON)
    arrays = Arrays()
    cases = {}
    for case in manifest['cases']:
        tpath = FIXTURE_BIN / case['torsion_file']
        assert sha256_file(tpath) == case['torsion_sha256'], f"hash mismatch {case['id']}"
        tors = np.fromfile(tpath, dtype='<f8').reshape(case['F'], case['ncols'])
        entry = {'ref': ref_case(case, tors, arrays)}
        if case['r2_p1']:
            assert sha256_file(FIXTURE_BIN / case['coords_file']) == case['coords_sha256']
            entry['p1'] = p1_case(case, tors, tf)
        cases[case['id']] = entry
    arr_path = ARRAY_BIN / 'SYNTHETIC_ref_arrays.f64'
    arr_path.parent.mkdir(parents=True, exist_ok=True)
    np.concatenate(arrays.chunks).astype('<f8').tofile(arr_path)
    dump_json({'_meta': {'type': 'SYNTHETIC', 'producer': 'run_python_path.py',
                         'python': sys.version.split()[0], 'numpy': np.__version__,
                         'platform': platform.platform(), 'topofold': module_identity(tf),
                         'arrays_file': str(arr_path.relative_to(Path(__file__).parent)),
                         'arrays_sha256': sha256_file(arr_path),
                         'n_cases': len(cases)},
               'cases': cases}, PY_JSON)
    print(f'{len(cases)} cases -> {PY_JSON}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
