#!/usr/bin/env python3
"""LE-1 computation measurements (protocol §9). Separate from retrieval quality.

Subset: open_start, W = 8, windows i = 1…20, all query × database frames. Each implementation runs in its own
subprocess (single thread) so peak RSS is attributable; the median of 3 repeats is reported. Setup (index fit,
feature extraction) is timed separately from steady-state distance evaluation and exhaustive top-10 retrieval.
NumPy vs Rust differences are implementation differences, not representation effects.
Run: ../wpd_local_retrieval/build/venv/bin/python -B timing_dev.py
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import json  # noqa: E402
import platform  # noqa: E402
import resource  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
BUILD = HERE / 'build'
SUBSET = BUILD / 'le_timing_subset.npz'
OUT = HERE / 'results' / 'le_compute_timing.json'
W, WINDOWS, REPEATS, K = 8, range(1, 21), 3, 10
IMPLS = ('baseline', 'np_features', 'M1', 'M2', 'K', 'A', 'F_lock', 'F_numpy', 'B1', 'B2', 'F_rust_fit',
         'F_rust_query', 'GF_rust_query', 'cascade_query')


def prepare():
    sys.path.insert(0, str(HERE.parent / 'wpd_local_retrieval'))
    import run_adk as WD
    CA, BB, _, _ = WD.load_traj('open_start')
    q, _, d = WD.split_ranges(CA.shape[0])
    BUILD.mkdir(exist_ok=True)
    np.savez(SUBSET, CA=CA, q=np.array(q), d=np.array(d))


def topk(R):
    return np.argpartition(R, K - 1, axis=1)[:, :K]


def run(impl):
    import topofold as tf
    from legacy_lib import (ca_pair_vectors, d_aligned, d_chord, d_circ, d_kappa, d_rms, frechet_legacy,
                            frechet_lockstep, horn_rmsd)
    Z = np.load(SUBSET)
    CA, (q0, q1), (d0, d1) = Z['CA'], Z['q'], Z['d']
    CA32 = CA.astype(np.float32)
    out = {'impl': impl}
    # feature extraction (binding, all frames used) — shared setup for the angular methods
    t = time.perf_counter()
    kap = np.empty((len(CA), 212))
    tau = np.empty((len(CA), 211))
    for f in range(len(CA)):
        k_, t_, _ = tf.compute_invariants(CA32[f])
        kap[f], tau[f] = k_, t_
    out['feature_extraction_binding_s'] = time.perf_counter() - t
    Qi, Di = np.arange(q0, q1), np.arange(d0, d1)
    times = []
    if impl == 'baseline':
        times = [0.0]
    elif impl in ('M1', 'M2', 'K', 'A', 'F_lock', 'F_numpy', 'B1', 'B2', 'np_features'):
        for _ in range(REPEATS):
            t = time.perf_counter()
            for i in WINDOWS:
                s = i - 1
                tw, kw, X = tau[:, s:s + W - 3], kap[:, s:s + W - 2], CA[:, s:s + W]
                if impl == 'np_features':
                    ca_pair_vectors(X)
                    continue
                for b in range(0, len(Qi), 128):
                    qb = Qi[b:b + 128]
                    if impl == 'M1':
                        R = d_chord(tw[qb], tw[Di])
                    elif impl == 'M2':
                        R = d_circ(tw[qb], tw[Di])
                    elif impl == 'K':
                        R = d_kappa(kw[qb], kw[Di])
                    elif impl == 'A':
                        R = d_aligned(kw[qb], tw[qb], kw[Di], tw[Di])
                    elif impl == 'F_lock':
                        R = frechet_lockstep(kw[qb], tw[qb], kw[Di], tw[Di])
                    elif impl == 'F_numpy':
                        R = frechet_legacy(kw[qb], tw[qb], kw[Di], tw[Di])
                    elif impl == 'B1':
                        P = ca_pair_vectors(X)
                        R = d_rms(P[qb], P[Di])
                    elif impl == 'B2':
                        R = horn_rmsd(X[qb], X[Di])
                    topk(R)
            times.append(time.perf_counter() - t)
    elif impl == 'F_rust_fit':
        for _ in range(REPEATS):
            t = time.perf_counter()
            idx = tf.ConformationalIndex(window_size=8, coarse_radius=0.5)
            idx.fit(np.ascontiguousarray(CA32[d0:d1]))
            times.append(time.perf_counter() - t)
    else:
        idx = tf.ConformationalIndex(window_size=8, coarse_radius=0.5)
        idx.fit(np.ascontiguousarray(CA32[d0:d1]))
        for _ in range(REPEATS):
            t = time.perf_counter()
            if impl == 'F_rust_query':
                for i in WINDOWS:
                    for qf in Qi:
                        idx.query_subcurve(i - 1, i + W - 2, CA32[qf], k=K)
            elif impl == 'GF_rust_query':
                for qf in Qi[:50]:
                    idx.query_subcurve(0, 213, CA32[qf], k=K)
            elif impl == 'cascade_query':
                for qf in Qi[:50]:
                    idx.query(CA32[qf], k=K)
            times.append(time.perf_counter() - t)
    out['steady_state_s_median'] = float(np.median(times))
    out['repeats_s'] = times
    out['units'] = ('20 windows × all queries × all candidates' if impl not in ('GF_rust_query', 'cascade_query',
                                                                            'F_rust_fit', 'baseline')
                    else ('50 full-chain queries' if impl != 'F_rust_fit' else 'one index fit over the database'))
    out['n_queries'], out['n_candidates'] = int(len(Qi)), int(len(Di))
    out['maxrss_kB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(json.dumps(out))


def main():
    if len(sys.argv) > 2 and sys.argv[1] == 'run':
        return run(sys.argv[2])
    if not SUBSET.exists():
        prepare()
    res = {}
    for impl in IMPLS:
        p = subprocess.run([sys.executable, '-B', __file__, 'run', impl], capture_output=True, text=True)
        line = [x for x in p.stdout.splitlines() if x.startswith('{')]
        res[impl] = json.loads(line[-1]) if line else {'error': p.stderr[-500:]}
        print(impl, res[impl].get('steady_state_s_median'), res[impl].get('maxrss_kB'), flush=True)
    base = res['baseline']['maxrss_kB']
    for v in res.values():
        if 'maxrss_kB' in v:
            v['maxrss_above_baseline_kB'] = v['maxrss_kB'] - base
    cpu = next((l.split(':', 1)[1].strip() for l in open('/proc/cpuinfo') if l.startswith('model name')), None)
    OUT.write_text(json.dumps({'_meta': {'producer': 'timing_dev.py', 'subset': 'open_start, W=8, windows 1-20',
                                         'cpu': cpu, 'python': platform.python_version(), 'numpy': np.__version__,
                                         'threads': 1, 'repeats': REPEATS,
                                         'note': 'NumPy vs Rust differences are implementation differences'},
                               'results': res}, indent=1) + '\n')


if __name__ == '__main__':
    main()
