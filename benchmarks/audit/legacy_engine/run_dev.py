#!/usr/bin/env python3
"""LE-1 development run on AdK (docs/LEGACY_ENGINE_EVALUATION_PROTOCOL.md). Third-party MD derivative.

Reuses the WP-D AdK loader, geometry, splits and exclusion rules unchanged (wpd_local_retrieval/run_adk.py).
Legacy F comes from the actual Rust path (ConformationalIndex.query_subcurve); the NumPy Fréchet is computed
on every (window, query, candidate) as a cross-check. Per-(window, query) NR@10 arrays go to build/ (git-ignored);
summaries and provenance to results/.
Run: ../wpd_local_retrieval/build/venv/bin/python -B run_dev.py --workers 4
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    os.environ.setdefault(_v, '1')

import argparse  # noqa: E402
import datetime  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import platform  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
WPD = HERE.parent / 'wpd_local_retrieval'
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(WPD))
import run_adk as WD  # noqa: E402
from legacy_lib import (ca_pair_vectors, d_aligned, d_chord, d_circ, d_kappa, d_rms,  # noqa: E402
                        frechet_legacy, frechet_lockstep, horn_rmsd, peptide_plane_atoms,
                        query_metrics)

RES = HERE / 'results'
BUILD = HERE / 'build'
PROTOCOL = REPO / 'docs' / 'LEGACY_ENGINE_EVALUATION_PROTOCOL.md'
FREEZE = RES / 'le_protocol_freeze.json'
OUT = RES / 'le_adk_THIRD_PARTY_MD_DERIVATIVE_summary.json'
OUT_PROV = RES / 'le_adk_THIRD_PARTY_MD_DERIVATIVE_provenance.json'
LOCAL = ('M1', 'M2', 'K', 'A', 'F', 'F_lock', 'B1', 'B2')
GLOBAL = ('G_RMSD', 'G_F')
METHODS = LOCAL + GLOBAL
TARGETS = ('phipsi', 'pp')
KS = (1, 10)
BATCH = 128
FRECHET_TOL = 1e-9
INDEX_KW = {'window_size': 8, 'coarse_radius': 0.5}

G = {}


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 22), b''):
            h.update(c)
    return h.hexdigest()


def load_oxygen(label):
    """Backbone carbonyl O of residues 1–213, made whole relative to its own C by minimum image."""
    import MDAnalysis as mda
    gro, xtc = WD.PAIRS[label]
    u = mda.Universe(str(WD.RAW / gro), str(WD.RAW / xtc))
    prot = u.select_atoms('protein')
    o = prot.select_atoms('name O')
    o = o[np.argsort(o.resids, kind='stable')]
    c = prot.select_atoms('name C')
    c = c[np.argsort(c.resids, kind='stable')]
    if list(o.resids) != list(range(1, 214)):
        raise WD.Stop(f'{label}: O atoms are not exactly residues 1–213')
    F = len(u.trajectory)
    O = np.full((F, WD.N_RES, 3), np.nan)
    changed = 0
    for f, ts in enumerate(u.trajectory):
        L = np.asarray(ts.dimensions[:3], float)
        oc = o.positions.astype(np.float64)
        cc = c.positions.astype(np.float64)[:213]
        ow = oc - L * np.round((oc - cc) / L)
        changed += bool(np.any(ow != oc))
        O[f, :213] = ow
    return O, changed


def legacy_index(label):
    """Per-process Rust index over the trajectory's database frames (fitted lazily in the worker)."""
    key = ('index', label)
    if key not in G:
        import topofold as tf
        d0, d1 = G[label]['d']
        idx = tf.ConformationalIndex(**INDEX_KW)
        idx.fit(np.ascontiguousarray(G[label]['CA'][d0:d1].astype(np.float32)))
        G[key] = idx
    return G[key]


def rust_rows(label, s, e, qframes, n_db):
    idx = legacy_index(label)
    CA32 = G[label]['CA32']
    out = np.empty((len(qframes), n_db))
    for r, q in enumerate(qframes):
        hits = idx.query_subcurve(s, e, CA32[q], k=n_db)
        if len(hits) != n_db:
            raise WD.Stop(f'{label}: query_subcurve returned {len(hits)} of {n_db}')
        for fid, d in hits:
            out[r, fid] = d
    return out


def run_global(task):
    """Exhaustive legacy full-chain Fréchet rows and writhe-cascade top-10 for a chunk of query frames."""
    label, q_lo, q_hi = task
    d0, d1 = G[label]['d']
    n_db = d1 - d0
    qf = np.arange(q_lo, q_hi)
    rows = rust_rows(label, 0, WD.N_RES - 1, qf, n_db)
    idx = legacy_index(label)
    casc = []
    for q in qf:
        hits = idx.query(G[label]['CA32'][q], k=10)
        casc.append([int(f) for f, _ in hits])
    return label, q_lo, rows, casc


def run_window(task):
    label, W, i = task
    g = G[label]
    geo, CA, (q0, q1), (d0, d1) = g['geo'], g['CA'], g['q'], g['d']
    s = i - 1
    valid, why = WD.window_validity(geo, s, W)
    qidx = np.arange(q0, q1)[valid[q0:q1]]
    didx = np.arange(d0, d1)[valid[d0:d1]]
    nq_all = q1 - q0
    rec = {'trajectory': label, 'W': W, 'i': i, 'n_query_valid': int(len(qidx)), 'n_db_valid': int(len(didx)),
           'excluded_query': {k: int(v[q0:q1].sum()) for k, v in why.items()},
           'excluded_db': {k: int(v[d0:d1].sum()) for k, v in why.items()}}
    nr10 = {(m, t): np.full(nq_all, np.nan, np.float32) for m in METHODS for t in TARGETS}
    if len(qidx) == 0 or len(didx) < max(KS):
        rec['undefined'] = 'too few valid frames'
        return rec, nr10
    tau = geo['tau'][:, s:s + W - 3]
    kap = geo['kap'][:, s:s + W - 2]
    X = CA[:, s:s + W]
    v = np.arange(s + 1, s + W - 1)
    tgt = np.concatenate([geo['phi'][:, v], geo['psi'][:, v]], axis=1)
    PP = peptide_plane_atoms(g['BB'], g['O'], s, W)
    pairs = ca_pair_vectors(X)
    acc = {(m, t, k): {'overlap': [], 'regret': [], 'norm': [], 'ndef': [], 'valid': []}
           for m in METHODS for t in TARGETS for k in KS}
    cnt = {(m, t): {'rounded': 0, 'failures': 0, 'nonfinite_distances': 0} for m in METHODS for t in TARGETS}
    fr_maxdiff = 0.0
    for b0 in range(0, len(qidx), BATCH):
        qb = qidx[b0:b0 + BATCH]
        T = {'phipsi': d_circ(tgt[qb], tgt[didx]), 'pp': horn_rmsd(PP[qb], PP[didx])}
        f_rust_all = rust_rows(label, s, s + W - 1, qb, d1 - d0)
        f_rust = f_rust_all[:, didx - d0]
        f_np = frechet_legacy(kap[qb], tau[qb], kap[didx], tau[didx])
        fr_maxdiff = max(fr_maxdiff, float(np.nanmax(np.abs(f_np - f_rust))))
        Rm = {'M1': d_chord(tau[qb], tau[didx]), 'M2': d_circ(tau[qb], tau[didx]),
              'K': d_kappa(kap[qb], kap[didx]), 'A': d_aligned(kap[qb], tau[qb], kap[didx], tau[didx]),
              'F': f_rust, 'F_lock': frechet_lockstep(kap[qb], tau[qb], kap[didx], tau[didx]),
              'B1': d_rms(pairs[qb], pairs[didx]), 'B2': horn_rmsd(X[qb], X[didx]),
              'G_RMSD': g['G_RMSD'][np.ix_(qb - q0, didx - d0)], 'G_F': g['G_F'][np.ix_(qb - q0, didx - d0)]}
        for m, R in Rm.items():
            nonfin = ~np.isfinite(R)
            for t in TARGETS:
                cnt[(m, t)]['nonfinite_distances'] += int(nonfin.sum())
            R = np.where(nonfin, np.inf, R)
            for t in TARGETS:
                for k in KS:
                    r = query_metrics(R, T[t], k)
                    a = acc[(m, t, k)]
                    a['overlap'].append(r['overlap'])
                    a['regret'].append(r['regret'])
                    a['norm'].append(r['norm_regret'])
                    a['ndef'].append(r['norm_defined'])
                    a['valid'].append(r['valid'])
                    if k == 10:
                        nr = np.where(r['norm_defined'], r['norm_regret'], np.nan)
                        nr10[(m, t)][qb - q0] = nr.astype(np.float32)
                        cnt[(m, t)]['rounded'] += r['n_regret_rounded']
                        cnt[(m, t)]['failures'] += r['n_regret_failures']
    if fr_maxdiff > FRECHET_TOL:
        rec['frechet_rust_numpy_violation'] = fr_maxdiff
    rec['frechet_rust_numpy_max_abs_diff'] = fr_maxdiff
    rec['random_overlap'] = {f'k{k}': k / len(didx) for k in KS}
    rec['methods'] = {}
    for m in METHODS:
        rec['methods'][m] = {}
        for t in TARGETS:
            e = {}
            for k in KS:
                a = {key: np.concatenate(v) for key, v in acc[(m, t, k)].items()}
                ok, nd = a['valid'], a['valid'] & a['ndef']
                e[f'k{k}'] = {'mean_overlap': float(a['overlap'][ok].mean()),
                              'mean_regret': float(a['regret'][ok].mean()),
                              'mean_norm_regret': float(a['norm'][nd].mean()) if nd.any() else None,
                              'n_norm_undefined': int((ok & ~a['ndef']).sum()), 'n_invalid': int((~ok).sum())}
            e.update(cnt[(m, t)])
            rec['methods'][m][t] = e
    return rec, nr10


def init_worker(g):
    G.update(g)


def git_state():
    head = subprocess.run(['git', 'rev-parse', 'HEAD'], capture_output=True, text=True, cwd=REPO).stdout.strip()
    st = subprocess.run(['git', 'status', '--porcelain'], capture_output=True, text=True, cwd=REPO).stdout
    return {'head': head, 'status_porcelain': st.rstrip('\n').splitlines()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=4)
    args = ap.parse_args()
    t0 = time.time()
    started = datetime.datetime.now(datetime.UTC).isoformat()
    fr = json.loads(FREEZE.read_text())
    if sha256(PROTOCOL) != fr['protocol_sha256']:
        raise WD.Stop('LE-1 protocol changed after freeze')
    for p, h in fr['legacy_source_sha256'].items():
        if sha256(REPO / p) != h:
            raise WD.Stop(f'legacy source {p} changed after freeze')
    tests = subprocess.run([sys.executable, '-B', '-m', 'pytest', '-p', 'no:cacheprovider', '-q',
                            str(HERE / 'test_legacy.py')], capture_output=True, text=True, cwd=HERE)
    test_summary = tests.stdout.strip().splitlines()[-1] if tests.stdout else ''
    if tests.returncode != 0:
        raise WD.Stop(f'synthetic verification failed: {test_summary}')
    inputs = WD.verify_inputs()
    g, traj = {}, {}
    for label in WD.PAIRS:
        CA, BB, times, wstats = WD.load_traj(label)
        geo, diag = WD.frame_geometry(label, CA, BB)
        O, o_changed = load_oxygen(label)
        q, gap, d = WD.split_ranges(CA.shape[0])
        Xq, Xd = CA[q[0]:q[1]], CA[d[0]:d[1]]
        g_rmsd = np.vstack([horn_rmsd(Xq[b:b + 64], Xd) for b in range(0, len(Xq), 64)])
        g[label] = {'geo': geo, 'CA': CA, 'CA32': CA.astype(np.float32), 'BB': BB, 'O': O, 'q': q, 'd': d,
                    'G_RMSD': g_rmsd}
        traj[label] = {'n_frames': int(CA.shape[0]), 'splits_half_open': {'query': q, 'gap': gap, 'database': d},
                       'whole_reconstruction': wstats, 'oxygen_whole_changed_frames': int(o_changed),
                       'reference_check': diag['reference_check'],
                       'binding_failure_frames': diag['binding_failure_frames']}
        print(label, 'geometry + global RMSD done', f'{time.time() - t0:.0f}s', flush=True)
    ctx = mp.get_context('fork')
    gtasks = [(label, lo, min(lo + 64, g[label]['q'][1])) for label in WD.PAIRS
              for lo in range(g[label]['q'][0], g[label]['q'][1], 64)]
    with ctx.Pool(args.workers, initializer=init_worker, initargs=(g,)) as pool:
        for label in WD.PAIRS:
            q0, q1 = g[label]['q']
            g[label]['G_F'] = np.empty((q1 - q0, g[label]['d'][1] - g[label]['d'][0]))
            g[label]['cascade'] = [None] * (q1 - q0)
        tg = time.time()
        for label, lo, rows, casc in pool.imap_unordered(run_global, gtasks):
            q0 = g[label]['q'][0]
            g[label]['G_F'][lo - q0:lo - q0 + len(rows)] = rows
            g[label]['cascade'][lo - q0:lo - q0 + len(casc)] = casc
        global_seconds = time.time() - tg
        print('global legacy Fréchet + cascade done', f'{time.time() - t0:.0f}s', flush=True)
    # cascade recall@10 vs exhaustive full-chain Fréchet (tie-aware membership)
    from wpd_lib import membership
    recall = {}
    for label in WD.PAIRS:
        Pm, _, _ = membership(g[label]['G_F'], 10)
        rec_q = np.array([Pm[r, c].sum() / 10 if len(c) else 0.0 for r, c in enumerate(g[label]['cascade'])])
        returned = np.array([len(c) for c in g[label]['cascade']])
        recall[label] = {'mean_recall_at_10': float(rec_q.mean()), 'min_recall_at_10': float(rec_q.min()),
                         'fraction_queries_recall_1': float(np.mean(rec_q >= 1 - 1e-12)),
                         'queries_with_fewer_than_10_hits': int(np.sum(returned < 10))}
    tasks = [(label, W, i) for label in WD.PAIRS for W in WD.WINDOWS for i in range(1, WD.N_RES + 2 - W)]
    windows = []
    store = {(label, W): {} for label in WD.PAIRS for W in WD.WINDOWS}
    with ctx.Pool(args.workers, initializer=init_worker, initargs=(g,)) as pool:
        for n, (rec, nr10) in enumerate(pool.imap_unordered(run_window, tasks, chunksize=1), 1):
            if 'frechet_rust_numpy_violation' in rec:
                raise WD.Stop(f"Rust/NumPy Fréchet mismatch {rec['frechet_rust_numpy_violation']} "
                              f"at {rec['trajectory']} W={rec['W']} i={rec['i']}")
            windows.append(rec)
            store[(rec['trajectory'], rec['W'])][rec['i']] = nr10
            if n % 50 == 0 or n == len(tasks):
                print(f'{n}/{len(tasks)} windows {time.time() - t0:.0f}s', flush=True)
    windows.sort(key=lambda w: (list(WD.PAIRS).index(w['trajectory']), w['W'], w['i']))
    BUILD.mkdir(exist_ok=True)
    arrays = {}
    for (label, W), d in store.items():
        ii = sorted(d)
        for m in METHODS:
            for t in TARGETS:
                arrays[f'{label}__W{W}__{m}__{t}'] = np.stack([d[i][(m, t)] for i in ii])
        arrays[f'{label}__W{W}__window_i'] = np.array(ii)
    npz = BUILD / 'le_adk_nr10_per_window_query.npz'
    np.savez_compressed(npz, **arrays)
    out = {'_meta': {'producer': 'run_dev.py', 'label': 'THIRD-PARTY MD DERIVATIVE (not SYNTHETIC)',
                     'data': 'Orädd et al. 2021 AdK MD, Zenodo 10.5281/zenodo.5583119, CC BY 4.0',
                     'tier': 'development (AdK); geometric retrieval benchmark, not biological validation',
                     'protocol': str(PROTOCOL.relative_to(REPO)), 'protocol_sha256': fr['protocol_sha256'],
                     'started_utc': started, 'finished_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                     'wall_seconds': round(time.time() - t0, 1), 'global_stage_seconds': round(global_seconds, 1),
                     'methods': METHODS, 'targets': TARGETS, 'index_parameters': INDEX_KW,
                     'per_window_query_nr10_npz': str(npz.relative_to(REPO)), 'npz_sha256': sha256(npz)},
           'trajectories': traj, 'cascade_recall_vs_exhaustive_G_F': recall, 'windows': windows}
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    import MDAnalysis
    import scipy
    prov = {'command': ' '.join([Path(sys.executable).name, '-B'] + sys.argv), 'tests': {
        'exit_code': tests.returncode, 'summary': test_summary}, 'wpd_inputs': inputs, 'git': git_state(),
        'harness_sha256': {n: sha256(HERE / n) for n in ('legacy_lib.py', 'run_dev.py', 'test_legacy.py')},
        'reused_sha256': {n: sha256(WPD / n) for n in ('run_adk.py', 'wpd_lib.py')},
        'versions': {'python': platform.python_version(), 'numpy': np.__version__,
                     'mdanalysis': MDAnalysis.__version__, 'scipy': scipy.__version__, 'workers': args.workers},
        'result_file': {'path': str(OUT.relative_to(REPO)), 'sha256': sha256(OUT)}}
    OUT_PROV.write_text(json.dumps(prov, indent=1) + '\n')
    print('wrote', OUT.name, f'{time.time() - t0:.0f}s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
