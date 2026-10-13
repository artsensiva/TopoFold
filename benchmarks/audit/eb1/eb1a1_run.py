#!/usr/bin/env python3
"""EB-1 Amendment 1 runner (docs/EB1_AMENDMENT_1.md). Development data only; outputs separate from the original run.

Order: verify freezes and reused caches → prep D1+ (closed_start via the operational S1/S2/S3 path) → correctness
C1′–C4′ (separate subprocess; stop on failure) → S2, S3 (A only), S4 fit and confirmation cells, each in a fresh
subprocess. Writes results/eb1a1_measurements.json. Run: ../wpd_local_retrieval/build/venv/bin/python -B eb1a1_run.py
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    os.environ[_v] = '1'

import datetime  # noqa: E402
import json  # noqa: E402
import resource  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import eb1_run as R0  # noqa: E402  (reuse sha256, mem/io helpers, constants)

RES = HERE / 'results'
C0 = HERE / 'build' / 'cache'
C1 = HERE / 'build' / 'cache_a1'
OUT = RES / 'eb1a1_measurements.json'
FREEZE = RES / 'eb1a1_freeze.json'
COEF = json.loads((RES / 'eb1a1_temp_calibration.json').read_text())['coef']
W = 8
ORIGINAL_CPU_H = json.loads((RES / 'eb1_measurements.json').read_text())['_meta']['aggregate_cpu_h']
CPU_STOP_H = 1.8
CLOSED_GRO = R0.REPO / 'benchmarks/audit/adk_feasibility/data/raw/ORADD_closed_state.gro'
CLOSED_XTC = R0.REPO / 'benchmarks/audit/adk_feasibility/data/raw/ORADD_closed_state_traj.xtc'
CLOSED_DERIVED = R0.REPO / 'benchmarks/audit/adk_feasibility/data/derived/closed_start_ca_whole.npy'


def windows_for(L, n_win):
    nmax = L - W + 1
    return list(range(nmax)) if n_win >= nmax else sorted(set(np.linspace(0, nmax - 1, n_win).astype(int).tolist()))


def design():
    cells = []
    for n in (563, 1126, 2253):
        cells += [{'block': 'S2', 'args': {'d': 'D1', 'n': n}, 'rep': r} for r in range(3)]
        cells += [{'block': 'S3', 'args': {'d': 'D1', 'n': n}, 'rep': r} for r in range(3)]
    for n in (405, 810, 1620):
        cells += [{'block': 'S2', 'args': {'d': 'D2', 'n': n}, 'rep': r} for r in range(3)]
        cells += [{'block': 'S3', 'args': {'d': 'D2', 'n': n}, 'rep': r} for r in range(3)]
    for m in ('A', 'B2'):
        fit = [('batch', nq, 1015, 20) for nq in (32, 64, 128, 192, 256)]
        fit += [('batch', nq, nd, 20) for nq in (128, 256) for nd in (250, 500)]
        fit += [('batch', 64, 1015, nw) for nw in (5, 80)]
        fit += [('single', 1, nd, nw) for nd in (250, 500, 1015) for nw in (20, 207)]
        for mode, nq, nd, nw in fit:
            cells.append({'block': 'S4', 'args': {'m': m, 'mode': mode, 'd': 'D1', 'nq': nq, 'nd': nd, 'nwin': nw,
                                                  'role': 'fit'}, 'rep': 0})
        conf = [('CW1', 'batch', 'D1', 512, 1015, 207), ('CW2', 'batch', 'D1p', 256, 2693, 40),
                ('CW3', 'single', 'D1p', 1, 2693, 207), ('CW4', 'batch', 'D2', 256, 405, 44),
                ('CW5', 'single', 'D2', 1, 405, 44)]
        for cid, mode, d, nq, nd, nw in conf:
            cells.append({'block': 'S4', 'args': {'m': m, 'mode': mode, 'd': d, 'nq': nq, 'nd': nd, 'nwin': nw,
                                                  'role': cid}, 'rep': 0})
    return cells


# --- child functions --------------------------------------------------------------------------------------------

def _dataset(d):
    src = {'D1': C0, 'D2': C0, 'D1p': C1}[d]
    meta = json.loads((src / f'{d}_meta.json').read_text())
    return src, meta


def run_cell(spec):
    import eb1_lib as E
    import eb1a1_lib as A1
    blk, a = spec['block'], spec['args']
    res = {}
    if blk == 'S2':
        if a['d'] == 'D1':
            st = A1.owned_subset(np.load(C0 / 'D1_stored.npy', mmap_mode='r'), np.arange(a['n']))
            bx = A1.owned_subset(np.load(C0 / 'D1_box.npy', mmap_mode='r'), np.arange(a['n']))
        else:
            parts = [np.load(C0 / f'D2_stored_{i}.npy') for i in range(4)]
            boxes = [np.repeat(np.load(C0 / f'D2_box_{i}.npy')[None], len(p), 0) for i, p in enumerate(parts)]
            st = np.ascontiguousarray(np.concatenate(parts)[:a['n']])
            bx = np.ascontiguousarray(np.concatenate(boxes)[:a['n']])
            del parts, boxes
        t0, c0 = time.perf_counter(), time.process_time()
        E.reconstruct(st, bx)
        res.update({'wall_s': time.perf_counter() - t0, 'cpu_s': time.process_time() - c0})
        L = st.shape[1]
        res.update({'x': a['n'] * L, 'P': a['n'] * L * (12 + 24 + 12 + 24)})
        return res
    if blk == 'S3':
        X = A1.owned_subset(np.load(C0 / f"{a['d']}_round_all.npy", mmap_mode='r'), np.arange(a['n']))
        E.features(X[:5])
        t0, c0 = time.perf_counter(), time.process_time()
        E.features(X)                                           # operational only: no NumPy reference here
        res.update({'wall_s': time.perf_counter() - t0, 'cpu_s': time.process_time() - c0})
        L = X.shape[1]
        res.update({'x': a['n'] * L, 'P': a['n'] * L * (24 + 12) + a['n'] * (2 * L - 5) * 8})
        return res
    if blk == 'S4':
        src, meta = _dataset(a['d'])
        qa, da = np.array(meta['qi']), np.array(meta['di'])
        rows = np.concatenate([qa[:a['nq']], da[:a['nd']]])
        if a['m'] == 'A':
            feat = {'kap': A1.owned_subset(np.load(src / f"{a['d']}_kap.npy", mmap_mode='r'), rows),
                    'tau': A1.owned_subset(np.load(src / f"{a['d']}_tau.npy", mmap_mode='r'), rows)}
            L = feat['kap'].shape[1] + 2
        else:
            feat = {'X': A1.owned_subset(np.load(src / f"{a['d']}_round.npy", mmap_mode='r'), rows)}
            L = feat['X'].shape[1]
        qi, di = np.arange(a['nq']), np.arange(a['nq'], a['nq'] + a['nd'])
        wins = windows_for(L, a['nwin'])
        batch = 1 if a['mode'] == 'single' else A1.BATCH_LIMIT
        A1.search_instrumented(a['m'], W, wins[:1], feat, qi[:min(8, len(qi))], di, batch)      # warm-up, excluded
        reps, parts = [], []
        for r in range(3):
            q = qi if a['mode'] == 'batch' else qi[:1]
            t0, c0 = time.perf_counter(), time.process_time()
            out, tp = A1.search_instrumented(a['m'], W, wins, feat, q, di, batch)
            reps.append((time.perf_counter() - t0, time.process_time() - c0))
            parts.append(tp)
            if r == 0 and reps[0][0] >= 5.0:
                break
        wall = float(np.median([x[0] for x in reps]))
        cpu = float(np.median([x[1] for x in reps]))
        wk = A1.work(a['nq'], a['nd'], len(wins))
        resident = sum(v.nbytes for v in feat.values())
        retained = sum(v.nbytes for v in out.values())
        res.update({'wall_s': wall, 'cpu_s': cpu, 'repeats': len(reps), 'repeat_wall_s': [x[0] for x in reps],
                    'decomposition_s': {k: float(np.median([p[k] for p in parts])) for k in parts[0]},
                    'n_win': len(wins), 'L': L, **wk, 'resident_bytes': resident, 'retained_bytes': retained,
                    'temp_bytes_pred': A1.temp_bytes(a['m'], wk['b'], a['nd'], W, COEF)})
        res['P'] = resident + retained + res['temp_bytes_pred']
        return res
    raise ValueError(blk)


def child_cell():
    out = run_cell(json.loads(sys.argv[2]))
    out['maxrss_kB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print('EB1RESULT ' + json.dumps(out))


def prep():
    """D1+ = D1 queries (open frames 0–1012) + database (open 1238–2252, then all closed frames), operational path."""
    import MDAnalysis as mda
    import eb1_lib as E
    C1.mkdir(parents=True, exist_ok=True)
    u = mda.Universe(str(CLOSED_GRO), str(CLOSED_XTC))
    ca = u.select_atoms('protein and name CA')
    ca = ca[np.argsort(ca.resids, kind='stable')]
    n = len(u.trajectory)
    st = np.empty((n, ca.n_atoms, 3), np.float32)
    box = np.empty((n, 3))
    for k, ts in enumerate(u.trajectory):
        st[k] = ca.positions
        box[k] = ts.dimensions[:3]
    whole, rnd = E.reconstruct(st, box)
    np.save(C1 / 'closed_whole_f64.npy', whole)
    kap, tau = E.features(rnd)
    X0, k0, t0 = np.load(C0 / 'D1_round.npy'), np.load(C0 / 'D1_kap.npy'), np.load(C0 / 'D1_tau.npy')
    m0 = json.loads((C0 / 'D1_meta.json').read_text())
    q, d = np.array(m0['qi']), np.array(m0['di'])
    np.save(C1 / 'D1p_round.npy', np.concatenate([X0[q], X0[d], rnd]))
    np.save(C1 / 'D1p_kap.npy', np.concatenate([k0[q], k0[d], kap]))
    np.save(C1 / 'D1p_tau.npy', np.concatenate([t0[q], t0[d], tau]))
    nq, nd = len(q), len(d) + n
    (C1 / 'D1p_meta.json').write_text(json.dumps({'qi': list(range(nq)), 'di': list(range(nq, nq + nd)),
                                                  'closed_frames': n, 'label': 'real AdK frames (open + closed)'}))
    np.save(C1 / 'closed_round.npy', rnd)
    np.save(C1 / 'closed_kap.npy', kap)
    np.save(C1 / 'closed_tau.npy', tau)
    print('EB1RESULT ' + json.dumps({'closed_frames': n, 'n_db': nd}))


def correctness():
    import eb1_lib as E
    import eb1a1_lib as A1
    out = {}
    whole = np.load(C1 / 'closed_whole_f64.npy')
    out['C1p'] = {'max_abs_diff_A': float(np.abs(whole - np.load(CLOSED_DERIVED)).max())}
    rnd, kap, tau = np.load(C1 / 'closed_round.npy'), np.load(C1 / 'closed_kap.npy'), np.load(C1 / 'closed_tau.npy')
    kn, tn = E.features_numpy(rnd)
    out['C2p'] = {'max_kappa_diff': float(np.abs(kap - kn).max()), 'max_tau_circ_diff': float(np.abs(E.wrap(tau - tn)).max())}
    X, K, T = np.load(C1 / 'D1p_round.npy'), np.load(C1 / 'D1p_kap.npy'), np.load(C1 / 'D1p_tau.npy')
    meta = json.loads((C1 / 'D1p_meta.json').read_text())
    qa, da = np.array(meta['qi']), np.array(meta['di'])
    rng = np.random.default_rng(R0.SEED + 1)
    qs, wins = rng.choice(qa, 20, replace=False), rng.choice(X.shape[1] - W + 1, 3, replace=False)
    c3 = {}
    for m in ('A', 'B2'):
        feat = {'kap': K, 'tau': T} if m == 'A' else {'X': X}
        n_ok = n_tie = 0
        mx = 0.0
        for s in wins:
            D = E.window_dist(m, W, int(s), X, K, T, qs, da)
            for j, q in enumerate(qs):
                ref = E.reference_distances(m, W, int(s), X, K, T, int(q), da)
                ok, tie, diff = E.compare_topk(D[j], ref)
                n_ok += ok
                n_tie += tie
                mx = max(mx, diff)
            got, _ = A1.search_instrumented(m, W, [int(s)], feat, qs, da)
            ref_s = E.search(m, W, [int(s)], X, K, T, qs, da)
            assert np.array_equal(got[int(s)], ref_s[int(s)])
        tol = 1e-10 if m == 'A' else 1e-6
        c3[m] = {'queries': 60, 'topk_agree': int(n_ok), 'near_ties': int(n_tie), 'max_abs_dist_diff': mx, 'tol': tol,
                 'pass': n_ok == 60 and mx <= tol}
    out['C3p'] = c3
    src, m2 = _dataset('D2')
    X2, K2, T2 = np.load(src / 'D2_round.npy'), np.load(src / 'D2_kap.npy'), np.load(src / 'D2_tau.npy')
    qi, di = np.array(m2['qi'])[:256], np.array(m2['di'])
    same = True
    for m in ('A', 'B2'):
        feat = {'kap': K2, 'tau': T2} if m == 'A' else {'X': X2}
        wins = windows_for(X2.shape[1], 44)
        g, _ = A1.search_instrumented(m, W, wins, feat, qi, di)
        r = E.search(m, W, wins, X2, K2, T2, qi, di)
        same &= all(np.array_equal(g[s], r[s]) for s in wins)
    out['C4p'] = {'instrumented_identical_to_harness_CW4': bool(same)}
    out['pass'] = (out['C1p']['max_abs_diff_A'] <= 1e-6 and out['C2p']['max_kappa_diff'] <= 1e-9
                   and out['C2p']['max_tau_circ_diff'] <= 1e-9 and all(v['pass'] for v in c3.values()) and same)
    print('EB1RESULT ' + json.dumps(out))


def main():
    if len(sys.argv) >= 2 and sys.argv[1] in ('cell', 'prep', 'correctness'):
        return {'cell': child_cell, 'prep': prep, 'correctness': correctness}[sys.argv[1]]()
    fr = json.loads(FREEZE.read_text())
    for rel, h in fr['sha256'].items():
        if R0.sha256(R0.REPO / rel) != h:
            raise SystemExit(f'STOP: {rel} changed after Amendment 1 freeze')
    orig = json.loads((RES / 'eb1_measurements.json').read_text())
    for name, h in orig['prep_cache_sha256'].items():
        if R0.sha256(C0 / name) != h:
            raise SystemExit(f'STOP: reused cache {name} differs from the original run')
    t0 = time.time()
    rec = {'_meta': {'producer': 'eb1a1_run.py', 'amendment_sha256': fr['sha256']['docs/EB1_AMENDMENT_1.md'],
                     'started_utc': datetime.datetime.now(datetime.UTC).isoformat(), 'original_cpu_h': ORIGINAL_CPU_H,
                     'reused': 'original EB-1 caches (hash-verified) and S1 measurements (analysis only)'}}

    def cpu_h():
        return ORIGINAL_CPU_H + R0.children_cpu() / 3600
    rec['prep'] = _spawn('prep')
    if rec['prep'].get('status') != 'ok':
        OUT.write_text(json.dumps(rec, indent=1) + '\n')
        return 2
    rec['cache_a1_sha256'] = {p.name: R0.sha256(p) for p in sorted(C1.iterdir())}
    rec['correctness'] = _spawn('correctness')
    OUT.write_text(json.dumps(rec, indent=1) + '\n')
    if rec['correctness'].get('status') != 'ok' or not rec['correctness'].get('pass'):
        rec['stopped'] = 'mandatory correctness failure'
        OUT.write_text(json.dumps(rec, indent=1) + '\n')
        return 2
    cells = design()
    rec['cells'] = []
    for i, c in enumerate(cells):
        if cpu_h() > CPU_STOP_H:
            rec['stopped'] = f'CPU budget reached before cell {i}'
            break
        rec['cells'].append({**c, **_spawn('cell', c)})
        if i % 10 == 0:
            print(f'{i}/{len(cells)} cpu_total {cpu_h():.3f} h', flush=True)
        OUT.write_text(json.dumps(rec, indent=1) + '\n')
    rec['_meta'].update({'finished_utc': datetime.datetime.now(datetime.UTC).isoformat(), 'wall_s': round(time.time() - t0, 1),
                         'amendment_cpu_h': R0.children_cpu() / 3600, 'aggregate_cpu_h_including_original': cpu_h(),
                         'n_cells_planned': len(cells), 'n_cells_ok': sum(c.get('status') == 'ok' for c in rec['cells'])})
    OUT.write_text(json.dumps(rec, indent=1) + '\n')
    print('done', rec['_meta']['n_cells_ok'], '/', len(cells), f"cpu_total {cpu_h():.3f} h")
    return 0


def _spawn(mode, spec=None):
    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (R0.AS_LIMIT, R0.AS_LIMIT))
    waited = 0
    while R0.mem_available_kb() < R0.MEM_MIN_KB and waited < R0.MEM_WAIT_S:
        time.sleep(15)
        waited += 15
    mem0 = R0.mem_available_kb()
    if mem0 < R0.MEM_MIN_KB:
        return {'status': 'skipped_low_memory', 'mem_available_kB': mem0}
    args = [sys.executable, '-B', __file__, mode] + ([json.dumps(spec)] if spec is not None else [])
    p = subprocess.run(args, capture_output=True, text=True, preexec_fn=limits, timeout=3600, cwd=HERE)
    line = [x for x in p.stdout.splitlines() if x.startswith('EB1RESULT ')]
    if p.returncode != 0 or not line:
        return {'status': 'failed', 'returncode': p.returncode, 'stderr': p.stderr[-1500:], 'mem_available_kB': mem0}
    out = json.loads(line[-1][10:])
    out.update({'status': 'ok', 'mem_available_kB_at_start': mem0})
    return out


if __name__ == '__main__':
    sys.exit(main())
