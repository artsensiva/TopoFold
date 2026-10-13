#!/usr/bin/env python3
"""EB-1 runner (frozen protocol docs/EB1_ENGINEERING_BENCHMARK_PROTOCOL.md). Development data only.

Order: verify freeze → prepare stage-input caches (build/, git-ignored; CPU counted) → M7 correctness (mandatory;
stop on failure) → M0–M6 cells, each in a fresh subprocess (BLAS threads 1, RLIMIT_AS 6 GB) → M8 if memory allows.
Writes results/eb1_measurements.json. Run: ../wpd_local_retrieval/build/venv/bin/python -B eb1_run.py
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    os.environ[_v] = '1'

import datetime  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import resource  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
RES = HERE / 'results'
CACHE = HERE / 'build' / 'cache'
FREEZE = RES / 'eb1_protocol_freeze.json'
PROTOCOL = REPO / 'docs' / 'EB1_ENGINEERING_BENCHMARK_PROTOCOL.md'
OUT = RES / 'eb1_measurements.json'
AS_LIMIT = 6 * 1024 ** 3
CPU_STOP_S = 1.8 * 3600
MEM_MIN_KB, MEM_WAIT_S = 500 * 1024, 600
SEED = 20261012


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 22), b''):
            h.update(c)
    return h.hexdigest()


def mem_available_kb():
    for line in open('/proc/meminfo'):
        if line.startswith('MemAvailable'):
            return int(line.split()[1])


def io_counters():
    d = {}
    for line in open('/proc/self/io'):
        k, v = line.split(':')
        d[k.strip()] = int(v)
    return d


# --- cell implementations (run inside a child process) ------------------------------------------------------------

def _load_feats(d):
    X = np.load(CACHE / f'{d}_round.npy')
    kap = np.load(CACHE / f'{d}_kap.npy')
    tau = np.load(CACHE / f'{d}_tau.npy')
    meta = json.loads((CACHE / f'{d}_meta.json').read_text())
    return X, kap, tau, meta


def run_cell(spec):
    import eb1_lib as E
    blk, a = spec['block'], spec['args']
    res = {}
    if blk == 'M0':
        import h5py  # noqa: F401
        import MDAnalysis  # noqa: F401
        import topofold  # noqa: F401
        t0 = time.perf_counter()
        res['wall_s'] = time.perf_counter() - t0
        return res
    if blk == 'S1':
        io0 = io_counters()
        t0, c0 = time.perf_counter(), time.process_time()
        if a['d'] == 'D1':
            ca, _, info = E.ingest_d1(a['n'])
            res['x'] = a['n'] * info['atoms_decoded_per_frame']
            res['out_bytes'] = ca.nbytes
        else:
            trajs = E.ingest_d2(a['n'])
            nfr = sum(t[0].shape[0] for t in trajs)
            res['x'] = nfr * 4 * trajs[0][0].shape[1]
            res['out_bytes'] = sum(t[0].nbytes for t in trajs)
        res['wall_s'], res['cpu_s'] = time.perf_counter() - t0, time.process_time() - c0
        io1 = io_counters()
        res['read_bytes'], res['rchar'] = io1['read_bytes'] - io0['read_bytes'], io1['rchar'] - io0['rchar']
        return res
    if blk == 'S2':
        if a['d'] == 'D1':
            ca = np.load(CACHE / 'D1_stored.npy')[:a['n']]
            box = np.load(CACHE / 'D1_box.npy')[:a['n']]
            t0, c0 = time.perf_counter(), time.process_time()
            E.reconstruct(ca, box)
            fr = ca.shape[0]
        else:
            fr = 0
            items = [(np.load(CACHE / f'D2_stored_{i}.npy'), np.load(CACHE / f'D2_box_{i}.npy')) for i in range(a['n'])]
            t0, c0 = time.perf_counter(), time.process_time()
            for ca, box in items:
                E.reconstruct(ca, box, 0.10)
                fr += ca.shape[0]
        res['wall_s'], res['cpu_s'] = time.perf_counter() - t0, time.process_time() - c0
        L = 214 if a['d'] == 'D1' else 51
        res['x'] = fr * L
        res['resident_bytes'] = fr * L * (12 + 24 + 24)
        return res
    if blk == 'S3':
        X = np.load(CACHE / f"{a['d']}_round_all.npy")[:a['n']]
        E.features(X[:5])                                   # warm-up, excluded
        t0, c0 = time.perf_counter(), time.process_time()
        E.features(X)
        res['wall_s'], res['cpu_s'] = time.perf_counter() - t0, time.process_time() - c0
        t1 = time.perf_counter()
        E.features_numpy(X)
        res['numpy_reference_wall_s'] = time.perf_counter() - t1
        res['x'] = X.shape[0] * X.shape[1]
        res['resident_bytes'] = X.nbytes + X.shape[0] * (2 * X.shape[1]) * 8
        return res
    if blk in ('S4single', 'S4batch'):
        X, kap, tau, meta = _load_feats(a['d'])
        q_all, d_all = np.array(meta['qi']), np.array(meta['di'])
        L, W = X.shape[1], a['W']
        if blk == 'S4single':
            nd = int(np.floor(a['frac'] * len(d_all)))
            di = d_all[:nd]
            windows = list(range(L - W + 1))
            rng = np.random.default_rng(SEED + W)
            qs = rng.choice(q_all, 5, replace=False)
            E.search(a['m'], W, windows[:2], X, kap, tau, qs[:1], di, batch=1)      # warm-up, excluded
            walls, cpus = [], []
            for q in qs:
                t0, c0 = time.perf_counter(), time.process_time()
                E.search(a['m'], W, windows, X, kap, tau, np.array([q]), di, batch=1)
                walls.append(time.perf_counter() - t0)
                cpus.append(time.process_time() - c0)
            res.update({'wall_s': float(np.median(walls)), 'cpu_s': float(np.median(cpus)), 'repeats_wall_s': walls,
                        'x': len(windows) * nd, 'n_D': nd, 'n_win': len(windows)})
        else:
            windows = (list(range(0, L - W + 1, 10)) if a['d'] == 'D1' else list(range(L - W + 1)))
            qi = q_all[:a['nq']]
            E.search(a['m'], W, windows[:1], X, kap, tau, qi[:8], d_all)                   # warm-up, excluded
            t0, c0 = time.perf_counter(), time.process_time()
            E.search(a['m'], W, windows, X, kap, tau, qi, d_all)
            res.update({'wall_s': time.perf_counter() - t0, 'cpu_s': time.process_time() - c0,
                        'x': len(qi) * len(windows) * len(d_all), 'n_Q': len(qi), 'n_win': len(windows),
                        'n_D': len(d_all)})
        res['resident_bytes'] = X.nbytes + kap.nbytes + tau.nbytes + E.BATCH * len(d_all) * 25 * 8
        return res
    if blk == 'R':
        sys.path.insert(0, str(E.WPD))
        if a['part'] == 'bootstrap':
            import heldout_decision_a2 as HD2
            r = np.random.default_rng(SEED)
            nr = {m: r.uniform(0.4, 0.9, 30) for m in ('B2', 'B1', 'B3', 'M1', 'M2', 'A', 'C1')}
            t0, c0 = time.perf_counter(), time.process_time()
            HD2.analyse_cell(nr, {m: v.copy() for m, v in nr.items()}, n_boot=2000)
            res.update({'wall_s': time.perf_counter() - t0, 'cpu_s': time.process_time() - c0})
            return res
        from wpd_lib import d_circ, query_metrics
        X, kap, tau, meta = _load_feats('D1')
        BB = np.load(CACHE / 'D1_backbone_round_targets.npz')
        phi, psi = BB['phi'], BB['psi']
        qi, di = np.array(meta['qi'])[:128], np.array(meta['di'])
        windows = list(range(0, 214 - 8 + 1, 10))
        t0, c0 = time.perf_counter(), time.process_time()
        for s in windows:
            v = np.arange(s + 1, s + 7)
            tq = np.concatenate([phi[qi][:, v], psi[qi][:, v]], 1)
            td = np.concatenate([phi[di][:, v], psi[di][:, v]], 1)
            T = d_circ(tq, td)
            if a['part'] == 'metrics':
                R = d_aligned_cache = E.window_dist('A', 8, s, X, kap, tau, qi, di)  # noqa: F841
                for k in (1, 10):
                    query_metrics(R, T, k)
        res.update({'wall_s': time.perf_counter() - t0, 'cpu_s': time.process_time() - c0, 'n_win': len(windows)})
        return res
    raise ValueError(blk)


def child_main():
    spec = json.loads(sys.argv[2])
    out = run_cell(spec)
    out['maxrss_kB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print('EB1RESULT ' + json.dumps(out))


# --- preparation and correctness (also run in children) -----------------------------------------------------------

def prep():
    import eb1_lib as E
    CACHE.mkdir(parents=True, exist_ok=True)
    ca, box, info = E.ingest_d1(2253)
    np.save(CACHE / 'D1_stored.npy', ca)
    np.save(CACHE / 'D1_box.npy', box)
    whole, rnd = E.reconstruct(ca, box)
    np.save(CACHE / 'D1_round_all.npy', rnd)
    qi, di = list(range(*E.D1_Q)), list(range(*E.D1_DB))
    np.save(CACHE / 'D1_round.npy', rnd)
    kap, tau = E.features(rnd)
    np.save(CACHE / 'D1_kap.npy', kap)
    np.save(CACHE / 'D1_tau.npy', tau)
    (CACHE / 'D1_meta.json').write_text(json.dumps({'qi': qi, 'di': di, 'atoms_decoded_per_frame':
                                                    info['atoms_decoded_per_frame']}))
    # research targets (φ/ψ from the backbone, as WP-D): reuse WP-D loader for N/CA/C
    sys.path.insert(0, str(E.WPD))
    import run_adk as WD
    CA, BBk, _, _ = WD.load_traj('open_start')
    N, A, C = BBk[:, :, 0], BBk[:, :, 1], BBk[:, :, 2]
    phi = np.full((CA.shape[0], 214), np.nan)
    psi = np.full((CA.shape[0], 214), np.nan)
    phi[:, 1:] = E.dihedral(C[:, :-1], N[:, 1:], A[:, 1:], C[:, 1:])
    psi[:, :-1] = E.dihedral(N[:, :-1], A[:, :-1], C[:, :-1], N[:, 1:])
    np.savez(CACHE / 'D1_backbone_round_targets.npz', phi=phi, psi=psi)
    trajs = E.ingest_d2(4)
    rounded = []
    for i, (c, b) in enumerate(trajs):
        np.save(CACHE / f'D2_stored_{i}.npy', c)
        np.save(CACHE / f'D2_box_{i}.npy', b)
        rounded.append(E.reconstruct(c, b, 0.10)[1])
    np.save(CACHE / 'D2_round_all.npy', np.concatenate(rounded))
    q, d = rounded[0], rounded[1]
    X = np.concatenate([q, d])
    np.save(CACHE / 'D2_round.npy', X)
    kap, tau = E.features(X)
    np.save(CACHE / 'D2_kap.npy', kap)
    np.save(CACHE / 'D2_tau.npy', tau)
    (CACHE / 'D2_meta.json').write_text(json.dumps({'qi': list(range(len(q))), 'di': list(range(len(q), len(X)))}))
    print('EB1RESULT ' + json.dumps({'prepared': sorted(p.name for p in CACHE.iterdir())}))


def correctness():
    import eb1_lib as E
    out = {}
    # C0 stored: D1 selection read vs full-system positions indexed at Cα indices (same decoder: consistency only)
    import MDAnalysis as mda
    u = mda.Universe(str(E.D1_GRO), str(E.D1_XTC))
    ca = u.select_atoms('protein and name CA')
    ca = ca[np.argsort(ca.resids, kind='stable')]
    stored = np.load(CACHE / 'D1_stored.npy')
    frames = [0, 1, 500, 1237, 1238, 2252]
    eq = []
    for f in frames:
        u.trajectory[f]
        eq.append(bool(np.array_equal(u.atoms.positions[ca.indices], stored[f])))
    out['C0_D1'] = {'frames': frames, 'bit_identical': all(eq), 'atoms': 'CA', 'order': 'resid ascending',
                    'dtype': 'float32', 'units': 'Angstrom', 'pbc': 'as stored (none applied)',
                    'reference': 'MDAnalysis full-system positions indexed at Cα indices (same decoder)'}
    import h5py
    from mdcath_reader_r1 import MdcathDomainR1
    d = MdcathDomainR1(E.D2_H5)
    ok = []
    for i, (t, r) in enumerate(E.D2_TRAJ):
        full = d.g[t][r]['coords'][:]
        ok.append(bool(np.array_equal(full[:, d.mapping['index']['CA']], np.load(CACHE / f'D2_stored_{i}.npy'))))
    d.close()
    out['C0_D2'] = {'trajectories': [f'{t}/{r}' for t, r in E.D2_TRAJ], 'bit_identical': all(ok),
                    'reference': 'h5py whole coords dataset read then indexed at Cα indices', 'h5py': h5py.__version__}
    # C1 reconstructed (float64 before rounding)
    whole, _ = E.reconstruct(stored, np.load(CACHE / 'D1_box.npy'))
    ref = np.load(E.D1_DERIVED)
    out['C1_D1'] = {'max_abs_diff_A': float(np.abs(whole - ref).max()), 'reference': 'WP-B derived open_start_ca_whole.npy '
                    '(separate run of the same reader and make_whole algorithm)', 'frames': int(ref.shape[0])}
    from mdcath_reader_r1 import make_whole_r1
    d = MdcathDomainR1(E.D2_H5)
    mx = 0.0
    for i, (t, r) in enumerate(E.D2_TRAJ):
        bb, O, _ = d.backbone_r1(t, r)
        wr, _ = make_whole_r1(bb, O, d.box_A(t, r))
        op, _ = E.reconstruct(np.load(CACHE / f'D2_stored_{i}.npy'), np.load(CACHE / f'D2_box_{i}.npy'), 0.0)
        mx = max(mx, float(np.abs(op - wr[:, :, 1]).max()))
    d.close()
    out['C1_D2'] = {'max_abs_diff_A': mx, 'reference': 'frozen harness make_whole_r1 along N–CA–C, Cα column'}
    # C2 features
    c2 = {}
    for dset in ('D1', 'D2'):
        X, kap, tau, _ = _load_feats(dset)
        kn, tn = E.features_numpy(X)
        c2[dset] = {'max_kappa_diff': float(np.abs(kap - kn).max()),
                    'max_tau_circ_diff': float(np.abs(E.wrap(tau - tn)).max())}
    out['C2'] = c2
    # C3 search
    c3 = {}
    rng = np.random.default_rng(SEED)
    for dset in ('D1', 'D2'):
        X, kap, tau, meta = _load_feats(dset)
        qa, da = np.array(meta['qi']), np.array(meta['di'])
        for W in (5, 8):
            wins = rng.choice(X.shape[1] - W + 1, 5, replace=False)
            qs = rng.choice(qa, 20, replace=False)
            for m in ('A', 'B2'):
                n_ok, n_tie, mxd = 0, 0, 0.0
                for s in wins:
                    D = E.window_dist(m, W, int(s), X, kap, tau, qs, da)
                    for j, q in enumerate(qs):
                        ref = E.reference_distances(m, W, int(s), X, kap, tau, int(q), da)
                        ok, tie, diff = E.compare_topk(D[j], ref)
                        n_ok += ok
                        n_tie += tie
                        mxd = max(mxd, diff)
                tol = 1e-10 if m == 'A' else 1e-6
                c3[f'{dset}/W{W}/{m}'] = {'queries': 100, 'topk_agree': int(n_ok), 'near_ties': int(n_tie),
                                         'max_abs_dist_diff': mxd, 'dist_tol': tol,
                                         'pass': n_ok == 100 and mxd <= tol}
    out['C3'] = c3
    out['pass'] = (out['C0_D1']['bit_identical'] and out['C0_D2']['bit_identical']
                   and out['C1_D1']['max_abs_diff_A'] <= 1e-6 and out['C1_D2']['max_abs_diff_A'] <= 1e-6
                   and all(v['max_kappa_diff'] <= 1e-9 and v['max_tau_circ_diff'] <= 1e-9 for v in c2.values())
                   and all(v['pass'] for v in c3.values()))
    print('EB1RESULT ' + json.dumps(out))


# --- orchestrator -------------------------------------------------------------------------------------------------

def matrix():
    cells = [{'block': 'M0', 'args': {}, 'rep': r} for r in range(3)]
    for blk in ('S1', 'S2', 'S3'):
        for n in (563, 1126, 2253):
            cells += [{'block': blk, 'args': {'d': 'D1', 'n': n}, 'rep': r} for r in range(3)]
        for n in (1, 2, 4):
            if blk == 'S3':
                cells += [{'block': blk, 'args': {'d': 'D2', 'n': n * 405}, 'rep': r} for r in range(3)]
            else:
                cells += [{'block': blk, 'args': {'d': 'D2', 'n': n}, 'rep': r} for r in range(3)]
    for d in ('D1', 'D2'):
        for W in (5, 8):
            for m in ('A', 'B2'):
                cells += [{'block': 'S4single', 'args': {'d': d, 'W': W, 'm': m, 'frac': fr}, 'rep': 0}
                          for fr in (0.25, 0.5, 1.0)]
                cells += [{'block': 'S4batch', 'args': {'d': d, 'W': W, 'm': m, 'nq': nq}, 'rep': 0}
                          for nq in (32, 64, 128)]
    cells += [{'block': 'R', 'args': {'part': p}, 'rep': 0} for p in ('targets', 'metrics', 'bootstrap')]
    return cells


def spawn(mode, spec=None, timeout=3600):
    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (AS_LIMIT, AS_LIMIT))
    waited = 0
    while mem_available_kb() < MEM_MIN_KB and waited < MEM_WAIT_S:
        time.sleep(15)
        waited += 15
    mem0 = mem_available_kb()
    if mem0 < MEM_MIN_KB:
        return {'status': 'skipped_low_memory', 'mem_available_kB': mem0}
    args = [sys.executable, '-B', __file__, mode] + ([json.dumps(spec)] if spec is not None else [])
    t0 = time.perf_counter()
    p = subprocess.run(args, capture_output=True, text=True, preexec_fn=limits, timeout=timeout, cwd=HERE)
    line = [x for x in p.stdout.splitlines() if x.startswith('EB1RESULT ')]
    if p.returncode != 0 or not line:
        return {'status': 'failed', 'returncode': p.returncode, 'stderr': p.stderr[-1500:], 'mem_available_kB': mem0}
    out = json.loads(line[-1][10:])
    out.update({'status': 'ok', 'mem_available_kB_at_start': mem0, 'process_wall_s': time.perf_counter() - t0})
    return out


def children_cpu():
    r = resource.getrusage(resource.RUSAGE_CHILDREN)
    s = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime + s.ru_utime + s.ru_stime


def main():
    if len(sys.argv) >= 2 and sys.argv[1] in ('cell', 'prep', 'correctness'):
        return {'cell': child_main, 'prep': prep, 'correctness': correctness}[sys.argv[1]]()
    fr = json.loads(FREEZE.read_text())
    if sha256(PROTOCOL) != fr['protocol_sha256']:
        raise SystemExit('STOP: EB-1 protocol changed after freeze')
    import eb1_lib as E
    for key, p in (('adk_open_xtc', E.D1_XTC), ('adk_open_gro', E.D1_GRO), ('adk_derived_open_ca_whole', E.D1_DERIVED),
                   ('mdcath_1ha8A00', E.D2_H5)):
        if sha256(p) != fr['inputs_sha256'][key]:
            raise SystemExit(f'STOP: input {key} changed')
    started = datetime.datetime.now(datetime.UTC).isoformat()
    t0 = time.time()
    rec = {'_meta': {'producer': 'eb1_run.py', 'protocol_sha256': fr['protocol_sha256'], 'started_utc': started,
                     'label': 'development data only (AdK open_start; mdCATH 1ha8A00 schema sample)',
                     'harness_sha256': {n: sha256(HERE / n) for n in ('eb1_lib.py', 'eb1_run.py')}}}
    print('prep', flush=True)
    rec['prep'] = spawn('prep', timeout=3600)
    if rec['prep'].get('status') != 'ok':
        rec['outcome_note'] = 'preparation failed'
        OUT.write_text(json.dumps(rec, indent=1) + '\n')
        return 2
    rec['prep_cache_sha256'] = {p.name: sha256(p) for p in sorted(CACHE.iterdir())}
    print('correctness', flush=True)
    rec['correctness'] = spawn('correctness', timeout=3600)
    OUT.write_text(json.dumps(rec, indent=1) + '\n')
    if rec['correctness'].get('status') != 'ok' or not rec['correctness'].get('pass'):
        rec['stopped'] = 'mandatory correctness failure'
        OUT.write_text(json.dumps(rec, indent=1) + '\n')
        print('STOP correctness', flush=True)
        return 2
    cells = matrix()
    rec['cells'] = []
    for i, c in enumerate(cells):
        if children_cpu() > CPU_STOP_S:
            rec['stopped'] = f'CPU budget reached before cell {i}'
            break
        r = spawn('cell', c)
        rec['cells'].append({**c, **r})
        if i % 10 == 0:
            print(f'{i}/{len(cells)} cells, cpu {children_cpu() / 3600:.3f} h', flush=True)
        OUT.write_text(json.dumps(rec, indent=1) + '\n')
    mem = mem_available_kb()
    rec['M8_parallel'] = ({'status': 'skipped', 'reason': f'MemAvailable {mem} kB < 2 GB (protocol §5)'}
                          if mem < 2 * 1024 ** 2 else {'status': 'not_implemented_in_this_run'})
    rec['_meta'].update({'finished_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                         'wall_s': round(time.time() - t0, 1), 'aggregate_cpu_h': children_cpu() / 3600,
                         'n_cells_planned': len(cells), 'n_cells_ok': sum(c.get('status') == 'ok' for c in rec['cells'])})
    OUT.write_text(json.dumps(rec, indent=1) + '\n')
    print('done', rec['_meta']['n_cells_ok'], '/', len(cells), f"cpu {rec['_meta']['aggregate_cpu_h']:.3f} h")
    return 0


if __name__ == '__main__':
    sys.exit(main())
