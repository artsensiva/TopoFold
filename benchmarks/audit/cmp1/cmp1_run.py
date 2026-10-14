#!/usr/bin/env python3
"""CMP-1 runner (docs/CMP1_COMPETITIVE_EVALUATION_PROTOCOL.md). Isolated env: build/venv (MDTraj 1.11.2, MDAnalysis
2.10.0, NumPy 2.5.3). Development data + SYNTHETIC_tiled performance-only rows.

Order: freeze checks → tests → references → gates (both tracks) → degeneracy suites → 5 seeded timing rounds (fresh
process per cell) . Memory start rule as BQ-2/BQ-3; stops before any step once CMP-1 CPU exceeds 1.8 h.
Writes results/cmp1_measurements.json. Run: build/venv/bin/python -B cmp1_run.py
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    os.environ[_v] = '1'

import datetime  # noqa: E402
import gc  # noqa: E402
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
RES = HERE / 'results'
REF = HERE / 'build' / 'ref'
OUT = RES / 'cmp1_measurements.json'
WL = json.loads((RES / 'cmp1_workloads.json').read_text())
IMPLS = ('T-V2b', 'T-V0', 'E-MDT', 'E-MDA')
ROUNDS, ORDER_SEED, DEGEN_SEED = 5, 20261022, 20261021
CPU_STOP_H = 1.8
AS_LIMIT = 6 * 1024 ** 3
HEADROOM = 600 * 1024 ** 2
EXT_IMPORT = 300 * 1024 ** 2


def wl(wid):
    return next(w for w in WL['workloads'] if w['id'] == wid)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def arr_hash(*arrs):
    h = hashlib.sha256()
    for a in arrs:
        h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


def maxrss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


def state():
    mi = {l.split(':')[0]: int(l.split()[1]) for l in open('/proc/meminfo') if ':' in l}
    return {'loadavg': open('/proc/loadavg').read().split()[:3], 'MemAvailable_kB': mi['MemAvailable'],
            'SwapFree_kB': mi['SwapFree']}


def threads():
    try:
        import threadpoolctl
        tp = [{'api': x.get('internal_api'), 'num_threads': x.get('num_threads')} for x in threadpoolctl.threadpool_info()]
    except Exception as e:
        tp = [{'error': repr(e)}]
    n = next(int(l.split()[1]) for l in open('/proc/self/status') if l.startswith('Threads:'))
    return {'threadpool': tp, 'os_threads': n}


def dims(w):
    L = 1000 if w['kind'] != 'real' else (51 if w['d'] == 'D2' else 214)
    return L, min(w['nwin'], L - 7)


def predicted(w, impl, extra=0):
    L, nw = dims(w)
    b = 124 * 1024 ** 2 + (3 + w['nd']) * L * 24 + (3706 * 214 * 24 if w['kind'] != 'real' else 0)
    prep = {'T-V2b': w['nd'] * nw * 200, 'T-V0': 0, 'E-MDA': w['nd'] * nw * 192,
            'E-MDT': w['nd'] * L * 12 + w['nd'] * nw * 100}[impl]
    return b + prep + (EXT_IMPORT if impl.startswith('E') else 0) + extra


# --- children ----------------------------------------------------------------------------------------------------------

def child_tests():
    p = subprocess.run([sys.executable, '-B', '-m', 'pytest', '-p', 'no:cacheprovider', '-q', 'test_cmp1.py'],
                       capture_output=True, text=True, cwd=HERE)
    tail = p.stdout.strip().splitlines()[-1] if p.stdout.strip() else ''
    return {'summary': tail, 'pass': p.returncode == 0}


def child_ref(wid):
    import cmp1_accuracy as ACC
    import cmp1_data as DATA
    w = wl(wid)
    X, di, wins, info = DATA.load(w)
    REF.mkdir(parents=True, exist_ok=True)
    D = np.lib.format.open_memmap(REF / f'{wid}.npy', mode='w+', dtype=np.float64, shape=(3, len(wins), len(di)))
    for q in range(3):
        D[q] = ACC.ref_distances(X, di, wins, q)
    fin = bool(all(np.isfinite(D[q]).all() for q in range(3)))
    D.flush()
    db = X[3:]
    dup = int(len(db) - len(np.unique(db.reshape(len(db), -1), axis=0)))
    return {'finite': fin, 'input': info, 'duplicate_database_rows': dup,
            'queries_in_database': int(sum(any(np.array_equal(X[q], r) for r in db) for q in range(3)))}


def child_gate(impl_name, wid):
    import cmp1_accuracy as ACC
    import cmp1_adapters as AD
    import cmp1_data as DATA
    w = wl(wid)
    X, di, wins, _ = DATA.load(w)
    h0 = arr_hash(X)
    impl = AD.IMPLS[impl_name]()
    Dr_all = np.load(REF / f'{wid}.npy', mmap_mode='r')
    hd = impl.prepare(X, di, wins)
    tracks = {t: {'pass': True, 'max_abs_err_A': 0.0, 'native_msd_ratio_descriptive': 0.0, 'rank_nonequivalent': 0,
                  'strict_identical': 0, 'near_tie': 0, 'windows': 0} for t in ('F64', 'F32')}
    hashes, consistent, finite = [], True, True
    for q in range(3):
        ti, td, _, _ = impl.query(hd, X, q)
        D = impl.distances(hd, X, q)
        finite &= bool(np.isfinite(D).all())
        TI = np.stack([AD.B.topk_lex(D[i]) for i in range(len(wins))])
        consistent &= bool(np.array_equal(TI, ti) and np.array_equal(np.take_along_axis(D, TI, 1), td))
        Dr = np.asarray(Dr_all[q])
        Gq = [ACC.window_G(X, [q], s)[0] for s in wins]
        Gd = [ACC.window_G(X, di, s) for s in wins]
        for t in tracks:
            e = ACC.evaluate(t, D, Dr, ti, Gq, Gd)
            r = tracks[t]
            r['pass'] &= e['pass']
            for k in ('max_abs_err_A', 'native_msd_ratio_descriptive'):
                r[k] = max(r[k], e[k])
            for k in ('rank_nonequivalent', 'strict_identical', 'near_tie', 'windows'):
                r[k] += e[k]
        hashes.append(arr_hash(ti, td))
        del D, Dr
    for t in tracks:
        tracks[t]['pass'] = bool(tracks[t]['pass'] and finite and consistent)
    return {'tracks': tracks, 'topk_sha256': hashes, 'query_path_consistent_with_distances': consistent,
            'finite': finite, 'input_unchanged': arr_hash(X) == h0, 'peak_rss_kB': maxrss()}


def child_degen(impl_name):
    import cmp1_accuracy as ACC
    import cmp1_adapters as AD
    return ACC.degeneracy_suite(AD.IMPLS[impl_name](), seed=DEGEN_SEED)


def child_time(impl_name, wid, rnd):
    t_imp = time.perf_counter()
    import cmp1_adapters as AD
    import cmp1_data as DATA
    impl = AD.IMPLS[impl_name]()                                  # library import happens here
    setup_s = time.perf_counter() - t_imp
    res = {'round': rnd, 'setup_import_s': setup_s, 'start': state(), 'threads': threads()}
    rng = np.random.default_rng(1)                                # identical SYNTHETIC warm-up for every implementation
    v = rng.normal(size=(24, 12, 3))
    Xw = np.cumsum(3.8 * v / np.linalg.norm(v, axis=-1, keepdims=True), axis=1).astype(np.float32).astype(np.float64)
    Xw.setflags(write=False)
    hw = impl.prepare(Xw, np.arange(3, 24), list(range(5)))
    impl.query(hw, Xw, 0)
    del hw, Xw
    gc.collect()
    res['rss_after_setup_kB'] = maxrss()
    t0 = time.perf_counter()
    X, di, wins, _ = DATA.load(wl(wid))
    res.update({'load_canonical_s': time.perf_counter() - t0, 'rss_after_load_kB': maxrss()})
    h0 = arr_hash(X)
    s0, c0 = time.perf_counter(), time.process_time()
    hd = impl.prepare(X, di, wins)
    p_w, p_c = time.perf_counter() - s0, time.process_time() - c0
    qw, qc, ks, hashes, nfb, out = [], [], [], [], 0, None
    for q in range(3):
        if out is not None:
            del out
            out = None
            gc.collect()
        a, c = time.perf_counter(), time.process_time()
        out = impl.query(hd, X, q)
        qw.append(time.perf_counter() - a)
        qc.append(time.process_time() - c)
        if q == 0:
            res['first_query_span_s'] = time.perf_counter() - s0
            res['rss_after_prepare_Q1_kB'] = maxrss()
        ks.append(out[2])
        nfb += out[3]
        hashes.append(arr_hash(out[0], out[1]))
    res.update({'prepare_wall_s': p_w, 'prepare_cpu_s': p_c, 'query_wall_s': qw, 'query_cpu_s': qc,
                'kernel_s': ks, 'fallback_rows': nfb, 'topk_sha256': hashes,
                'retained_output_bytes': int(out[0].nbytes + out[1].nbytes),
                'prepared_bytes': impl.prepared_bytes(hd), 'peak_rss_kB': maxrss(),
                'input_unchanged': arr_hash(X) == h0, 'end': state()})
    return res


# --- parent ------------------------------------------------------------------------------------------------------------

def spawn(args, pred):
    def lim():
        resource.setrlimit(resource.RLIMIT_AS, (AS_LIMIT, AS_LIMIT))

    def avail():
        return int(open('/proc/meminfo').read().split('MemAvailable:')[1].split()[0]) * 1024
    waited = 0
    while avail() - pred < HEADROOM and waited < 600:
        time.sleep(15)
        waited += 15
    if avail() - pred < HEADROOM:
        return {'status': 'skipped_memory_rule', 'predicted_peak_bytes': pred, 'mem_available_bytes': avail()}
    t = datetime.datetime.now(datetime.UTC).isoformat()
    p = subprocess.run([sys.executable, '-B', __file__] + args, capture_output=True, text=True, preexec_fn=lim,
                       timeout=5400, cwd=HERE)
    line = [x for x in p.stdout.splitlines() if x.startswith('CMP1RESULT ')]
    if p.returncode != 0 or not line:
        return {'status': 'failed', 'returncode': p.returncode, 'stderr': p.stderr[-1500:], 'started_utc': t}
    r = json.loads(line[-1][11:])
    r.update({'status': 'ok', 'started_utc': t, 'predicted_peak_bytes': pred, 'waited_s': waited})
    return r


def cpu_h():
    r, s = resource.getrusage(resource.RUSAGE_CHILDREN), resource.getrusage(resource.RUSAGE_SELF)
    return (r.ru_utime + r.ru_stime + s.ru_utime + s.ru_stime) / 3600


def main():
    if len(sys.argv) >= 2 and sys.argv[1] in ('tests', 'ref', 'gate', 'degen', 'time'):
        fn = {'tests': child_tests, 'ref': child_ref, 'gate': child_gate, 'degen': child_degen, 'time': child_time}
        a = sys.argv[2:]
        r = fn[sys.argv[1]](*([a[0], a[1], int(a[2])] if sys.argv[1] == 'time' else a))
        print('CMP1RESULT ' + json.dumps(r, default=float))
        return 0
    for rec in ('cmp1_protocol_freeze.json', 'cmp1_pre_measurement_record.json'):
        for rel, h in json.loads((RES / rec).read_text())['sha256'].items():
            if sha(REPO / rel) != h:
                raise SystemExit(f'STOP: {rel} changed ({rec})')
    M = {'_meta': {'producer': 'cmp1_run.py', 'started_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                   'cpu_cap_h': 2.0, 'cpu_stop_h': CPU_STOP_H, 'rounds': ROUNDS, 'order_seed': ORDER_SEED,
                   'start_state': state()}, 'steps': [], 'ref': {}, 'gate': {}, 'degeneracy': {}, 'reps': []}

    def save():
        M['_meta']['cpu_h_so_far'] = cpu_h()
        OUT.write_text(json.dumps(M, indent=1, default=float) + '\n')

    def run(name, args, pred):
        if cpu_h() > CPU_STOP_H:
            M['stopped'] = f'CPU budget reached before {name}'
            return None
        r = spawn(args, pred)
        M['steps'].append({'step': name, 'status': r.get('status'), 'cpu_h_after': cpu_h(), 'started_utc': r.get('started_utc')})
        save()
        return r
    M['tests'] = run('tests', ['tests'], 300 * 1024 ** 2)
    if not M['tests'] or not M['tests'].get('pass'):
        M.setdefault('stopped', 'tests failed')
        save()
        return 2
    for w in WL['workloads']:
        L, nw = dims(w)
        M['ref'][w['id']] = run(f"ref {w['id']}", ['ref', w['id']], predicted(w, 'T-V0', 3 * nw * w['nd'] * 8))
    for w in WL['workloads']:
        L, nw = dims(w)
        M['gate'][w['id']] = {}
        for i in IMPLS:
            if (M['ref'][w['id']] or {}).get('status') != 'ok':
                M['gate'][w['id']][i] = {'status': 'no_reference'}
                continue
            M['gate'][w['id']][i] = run(f"gate {i} {w['id']}", ['gate', i, w['id']], predicted(w, i, 2 * nw * w['nd'] * 8))
    for i in IMPLS:
        M['degeneracy'][i] = run(f'degeneracy {i}', ['degen', i], 600 * 1024 ** 2)
    rng = np.random.default_rng(ORDER_SEED)
    cells = [(w['id'], i) for w in WL['workloads'] for i in IMPLS]
    seq = 0
    for r in range(ROUNDS):
        for k in rng.permutation(len(cells)).tolist():
            wid, i = cells[k]
            g = M['gate'][wid].get(i) or {}
            if g.get('status') != 'ok' or not g['tracks']['F32']['pass']:
                M['reps'].append({'round': r, 'wid': wid, 'impl': i, 'status': 'not_timed_gate_F32_not_passed'})
                continue
            if cpu_h() > CPU_STOP_H:
                M['stopped'] = f'CPU budget reached in round {r}'
                M['reps'].append({'round': r, 'wid': wid, 'impl': i, 'status': 'not_run_budget'})
                continue
            seq += 1
            x = spawn(['time', i, wid, str(r)], predicted(wl(wid), i))
            M['reps'].append({'seq': seq, 'wid': wid, 'impl': i, **x, 'round': r})
            save()
    M['_meta'].update({'finished_utc': datetime.datetime.now(datetime.UTC).isoformat(), 'cpu_h_total': cpu_h(),
                       'end_state': state()})
    save()
    print('done', sum(x.get('status') == 'ok' for x in M['reps']), 'timed cells ok;', f'cpu {cpu_h():.3f} h')
    return 0


if __name__ == '__main__':
    sys.exit(main())
