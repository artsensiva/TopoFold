#!/usr/bin/env python3
"""BQ-2 runner (docs/B2_SCALING_CONFIRMATION_PROTOCOL.md). Development data + SYNTHETIC_tiled performance-only rows.

Order: freeze checks → tests → references (per workload) → V2b gate → V0 gate → timing rounds → P2D (if memory
permits) . Every subprocess obeys the frozen memory start rule (predicted peak + 600 MB ≤ MemAvailable; wait ≤ 10
min, else skipped). Stops before any step once BQ-2 CPU (parent + children) exceeds 0.9 h.
Writes results/bq2_measurements.json. Run: ../wpd_local_retrieval/build/venv/bin/python -B bq2_run.py
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
EB1 = HERE.parent / 'eb1'
B1DIR = HERE.parent / 'b2opt'
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(EB1))
import eb1_run as R0  # noqa: E402

RES = HERE / 'results'
REF = HERE / 'build' / 'ref'
OUT = RES / 'bq2_measurements.json'
WL = json.loads((RES / 'bq2_workloads.json').read_text())
CPU_STOP_H = 0.9
ROUNDS = 3
ORDER_SEED = 20261019


def workload(wid):
    return next(w for w in WL['workloads'] if w['id'] == wid)


def b1_prep_pred(w):
    a = json.loads((B1DIR / 'results' / 'b2opt_analysis.json').read_text())['arms']['V2b']['models']['prep_wall']
    return a['c0'] + a['c1'] * w['nd'] * dims(w)[1]


def dims(w):
    L = 1000 if w['kind'] != 'real' else (51 if w['d'] == 'D2' else 214)
    return L, min(w['nwin'], L - 7)


def load(w):
    """Owned array: rows 0–2 queries, rows 3… the first n_D database rows. Returns (X, di, windows)."""
    import eb1a1_run as R1
    import eb1a2_lib as L2
    src = EB1 / 'build' / ('cache_a1' if w['d'] == 'D1p' else 'cache')
    meta = json.loads((src / f"{w['d']}_meta.json").read_text())
    if w['kind'] == 'real':
        L2.validate_size(w['nd'], len(meta['di']), 'n_D')
        rows = [meta['qi'][q] for q in w['queries']] + meta['di'][:w['nd']]
        runs = []
        for i in rows:
            if runs and runs[-1][0] + runs[-1][1] == i:
                runs[-1][1] += 1
            else:
                runs.append([i, 1])
        X = L2.read_rows(src / f"{w['d']}_round.npy", [tuple(r) for r in runs])
    else:
        import bq2_lib as Q
        full = L2.read_rows(src / 'D1p_round.npy', [(0, L2.npy_shape(src / 'D1p_round.npy')[0][0])])
        X = np.empty((3 + w['nd'], Q.L_STRESS, 3))
        X[:3] = Q.tiled_rows(full, [Q.QUERY_OFFSET + q for q in w['queries']])
        for s in range(0, w['nd'], 500):
            r = list(range(s, min(s + 500, w['nd'])))
            X[3 + s:3 + s + len(r)] = Q.tiled_rows(full, r)
        del full
        gc.collect()
    wins = R1.windows_for(X.shape[1], w['nwin'])
    assert len(wins) == dims(w)[1]
    return X, np.arange(3, 3 + w['nd']), wins


def topk_hash(out):
    return hashlib.sha256(np.ascontiguousarray(out, np.int64).tobytes()).hexdigest()


def maxrss_kb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


def loadavg():
    return open('/proc/loadavg').read().split()[:3]


def threads_info():
    try:
        import threadpoolctl
        return [{'api': x.get('internal_api'), 'num_threads': x.get('num_threads')}
                for x in threadpoolctl.threadpool_info()]
    except Exception as e:
        return [{'error': repr(e)}]


# --- children --------------------------------------------------------------------------------------------------------

def child_tests():
    p = subprocess.run([sys.executable, '-B', '-m', 'pytest', '-p', 'no:cacheprovider', '-q', 'test_bq2.py'],
                       capture_output=True, text=True, cwd=HERE)
    tail = p.stdout.strip().splitlines()[-1] if p.stdout.strip() else ''
    print('EB1RESULT ' + json.dumps({'returncode': p.returncode, 'summary': tail, 'pass': p.returncode == 0}))


def child_ref():
    import bq2_lib as Q
    w = workload(sys.argv[2])
    REF.mkdir(parents=True, exist_ok=True)
    X, di, wins = load(w)
    nq = 1 if w['role'] == 'direct_P2_if_memory' else 3
    D = np.empty((nq, len(wins), len(di)))
    for qi in range(nq):
        for i, s in enumerate(wins):
            D[qi, i] = Q.B.ref_svd(X[qi, s:s + 8], X[di, s:s + 8])
    np.save(REF / f"{w['id']}.npy", D)
    print('EB1RESULT ' + json.dumps({'finite': bool(np.isfinite(D).all()), 'pairs': int(D.size),
                                     'sha256': hashlib.sha256(D.tobytes()).hexdigest()}))


def child_gate():
    import bq2_lib as Q
    arm, w = sys.argv[2], workload(sys.argv[3])
    X, di, wins = load(w)
    D_ref = np.load(REF / f"{w['id']}.npy", mmap_mode='r')        # one query slice resident at a time
    prep = None if arm == 'V0' else Q.B.prepare(X, di, wins)
    err, n_bad, n_strict, n_near, nfb, hashes = 0.0, 0, 0, 0, 0, []
    for qi in range(D_ref.shape[0]):
        out, D, f = Q.B.search(arm, X, qi, di, wins, prep, keep_dist=True)
        nfb += f
        Dr = np.asarray(D_ref[qi])
        err = max(err, float(np.max(np.abs(D - Dr))) if np.isfinite(D).all() else float('inf'))
        for i in range(len(wins)):
            eq, st, nr = Q.B.rank_check(out[i], Dr[i])
            n_bad += not eq
            n_strict += st
            n_near += nr
        hashes.append(topk_hash(out))
        del out, D, Dr
    qw = D_ref.shape[0] * len(wins)
    print('EB1RESULT ' + json.dumps({'arm': arm, 'max_abs_dist_err_A': err, 'query_windows': qw,
                                     'rank_nonequivalent': n_bad, 'rank_strict_identical': n_strict,
                                     'near_tie_query_windows': n_near, 'fallback_rows': nfb,
                                     'fallback_fraction': nfb / (qw * len(di)), 'topk_sha256': hashes,
                                     'pass': bool(err <= Q.B.TOL_D and n_bad == 0)}))


def child_time():
    import bq2_lib as Q
    spec = json.loads(sys.argv[2])
    arm, w = spec['arm'], workload(spec['wid'])
    la0, ma0 = loadavg(), R0.mem_available_kb()
    rng = np.random.default_rng(1)
    Xs = Q.B._chains(rng, 24, L=12)
    ws = list(range(5))
    Q.B.search(arm, Xs, 0, np.arange(1, 24), ws, None if arm == 'V0' else Q.B.prepare(Xs, np.arange(1, 24), ws))
    del Xs
    gc.collect()
    rss0 = maxrss_kb()
    X, di, wins = load(w)
    res = {'threads': threads_info(), 'rss_after_warmup_kB': rss0, 'rss_after_load_kB': maxrss_kb(),
           'resident_bytes': int(X.nbytes), 'n_win': len(wins), 'n_d': len(di), 'loadavg_start': la0,
           'mem_available_kB_start': ma0}
    prep = None
    if arm != 'V0':
        reps = Q.prep_reps_for(spec['pred_prep_s'])
        prep, pw, pc = Q.prepare_reps(X, di, wins, reps)
        res.update({'prep_reps': reps, 'prep_wall_reps_s': pw, 'prep_cpu_reps_s': pc,
                    'prep_wall_s': float(np.median(pw)), 'prep_cpu_s': float(np.median(pc)),
                    'prepared_bytes': int(prep[0].nbytes + prep[1].nbytes)})
    else:
        res.update({'prep_reps': 0, 'prep_wall_s': 0.0, 'prep_cpu_s': 0.0, 'prepared_bytes': 0})
    res['rss_after_prep_kB'] = maxrss_kb()
    walls, cpus, hashes, nfb, out = [], [], [], 0, None
    for qi in range(spec.get('n_queries', 3)):
        if out is not None:
            del out
            out = None
            gc.collect()
        t0, c0 = time.perf_counter(), time.process_time()
        out, _, f = Q.B.search(arm, X, qi, di, wins, prep)
        walls.append(time.perf_counter() - t0)
        cpus.append(time.process_time() - c0)
        nfb += f
        hashes.append(topk_hash(out))
    res.update({'query_wall_s': walls, 'query_cpu_s': cpus, 'topk_sha256': hashes, 'fallback_rows': nfb,
                'retained_bytes': int(out.nbytes), 'rss_end_kB': maxrss_kb(), 'loadavg_end': loadavg(),
                'mem_available_kB_end': R0.mem_available_kb()})
    print('EB1RESULT ' + json.dumps(res))


# --- parent ----------------------------------------------------------------------------------------------------------

def predicted(w, arm, extra=0):
    import bq2_lib as Q
    L, nwin = dims(w)
    nrows = 3 + w['nd'] + (3706 if w['kind'] != 'real' else 0)
    nd_prep = w['nd'] if arm != 'V0' else 0
    return Q.predicted_peak_bytes(nrows, L, nd_prep, nwin) + extra


def _spawn(mode, args, pred_bytes):
    import bq2_lib as Q

    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (R0.AS_LIMIT, R0.AS_LIMIT))
    waited = 0
    while not Q.memory_ok(pred_bytes, R0.mem_available_kb()) and waited < R0.MEM_WAIT_S:
        time.sleep(15)
        waited += 15
    mem0 = R0.mem_available_kb()
    if not Q.memory_ok(pred_bytes, mem0):
        return {'status': 'skipped_memory_rule', 'mem_available_kB': mem0, 'predicted_peak_bytes': pred_bytes,
                'waited_s': waited}
    t_start = datetime.datetime.now(datetime.UTC).isoformat()
    p = subprocess.run([sys.executable, '-B', __file__, mode] + args, capture_output=True, text=True,
                       preexec_fn=limits, timeout=3600, cwd=HERE)
    line = [x for x in p.stdout.splitlines() if x.startswith('EB1RESULT ')]
    if p.returncode != 0 or not line:
        return {'status': 'failed', 'returncode': p.returncode, 'stderr': p.stderr[-1500:], 'started_utc': t_start}
    out = json.loads(line[-1][10:])
    out.update({'status': 'ok', 'mem_available_kB_at_start': mem0, 'predicted_peak_bytes': pred_bytes,
                'waited_s': waited, 'started_utc': t_start})
    return out


def main():
    if len(sys.argv) >= 2 and sys.argv[1] in ('tests', 'ref', 'gate', 'time'):
        return {'tests': child_tests, 'ref': child_ref, 'gate': child_gate, 'time': child_time}[sys.argv[1]]()
    for rec_name in ('bq2_protocol_freeze.json', 'bq2_pre_measurement_record.json'):
        fr = json.loads((RES / rec_name).read_text())
        for key in ('sha256', 'bq1_preserved'):
            for rel, h in fr.get(key, {}).items():
                if R0.sha256(R0.REPO / rel) != h:
                    raise SystemExit(f'STOP: {rel} changed ({rec_name})')
    orig = json.loads((EB1 / 'results' / 'eb1_measurements.json').read_text())
    for name, h in orig['prep_cache_sha256'].items():
        if R0.sha256(EB1 / 'build' / 'cache' / name) != h:
            raise SystemExit(f'STOP: cache {name} differs from EB-1')
    a1 = json.loads((EB1 / 'results' / 'eb1a1_measurements.json').read_text())
    for name, h in a1['cache_a1_sha256'].items():
        if R0.sha256(EB1 / 'build' / 'cache_a1' / name) != h:
            raise SystemExit(f'STOP: cache_a1 {name} differs from Amendment 1')
    t_all = time.time()
    rec = {'_meta': {'producer': 'bq2_run.py', 'started_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                     'cpu_cap_h': 1.0, 'cpu_stop_h': CPU_STOP_H, 'rounds': ROUNDS, 'order_seed': ORDER_SEED,
                     'loadavg_start': loadavg(), 'mem_available_kB_start': R0.mem_available_kb()},
           'steps': [], 'ref': {}, 'gate': {'V2b': {}, 'V0': {}}, 'cells': []}

    def cpu_h():
        return R0.children_cpu() / 3600

    def save():
        rec['_meta']['cpu_h_so_far'] = cpu_h()
        OUT.write_text(json.dumps(rec, indent=1) + '\n')

    def run(name, mode, args, pred):
        if cpu_h() > CPU_STOP_H:
            rec['stopped'] = f'CPU budget reached before {name}'
            return None
        r = _spawn(mode, args, pred)
        rec['steps'].append({'step': name, 'status': r.get('status'), 'cpu_h_after': cpu_h()})
        save()
        return r
    rec['tests'] = run('tests', 'tests', [], 300 * 1024 ** 2)
    save()
    if not rec['tests'] or not rec['tests'].get('pass'):
        rec.setdefault('stopped', 'tests failed')
        save()
        return 2
    gated = [w for w in WL['workloads'] if w['role'] != 'direct_P2_if_memory']
    for w in gated:
        L, nwin = dims(w)
        rec['ref'][w['id']] = run(f"ref {w['id']}", 'ref', [w['id']], predicted(w, 'V0', 2 * 3 * nwin * w['nd'] * 8))
    for arm in ('V2b', 'V0'):
        for w in gated:
            if arm == 'V0' and not w['baseline']:
                continue
            if (rec['ref'].get(w['id']) or {}).get('status') != 'ok':
                rec['gate'][arm][w['id']] = {'status': 'no_reference'}
                continue
            L, nwin = dims(w)
            rec['gate'][arm][w['id']] = run(f"gate {arm} {w['id']}", 'gate', [arm, w['id']],
                                            predicted(w, arm, 2 * nwin * w['nd'] * 8))
        if arm == 'V2b' and any(g.get('status') == 'ok' and not g.get('pass') for g in rec['gate']['V2b'].values()):
            rec['stopped'] = 'V2b gate failure (measured correctness failure)'
            save()
            return 2
    save()
    rng = np.random.default_rng(ORDER_SEED)
    seq = 0
    stopped = False
    for r in range(ROUNDS):
        for wid in rng.permutation([w['id'] for w in gated]).tolist():
            w = workload(wid)
            arms = ['V2b', 'V0'] if w['baseline'] else ['V2b']
            for arm in rng.permutation(arms).tolist():
                g = rec['gate'][arm].get(wid) or {}
                if g.get('status') != 'ok' or not g.get('pass'):
                    rec['cells'].append({'seq': None, 'round': r, 'wid': wid, 'arm': arm, 'status': 'not_gated'})
                    continue
                if cpu_h() > CPU_STOP_H:
                    rec['stopped'] = 'CPU budget reached during timing'
                    stopped = True
                    break
                seq += 1
                spec = {'arm': arm, 'wid': wid, 'pred_prep_s': b1_prep_pred(w)}
                c = _spawn('time', [json.dumps(spec)], predicted(w, arm))
                rec['cells'].append({'seq': seq, 'round': r, 'wid': wid, 'arm': arm, **c})
                save()
            if stopped:
                break
        if stopped:
            break
    if not stopped:                                             # P2D: direct P2 check, only if the memory rule allows
        w = workload('P2D')
        L, nwin = dims(w)
        pred = predicted(w, 'V2b', nwin * w['nd'] * 8 * 2)
        spec = {'arm': 'V2b', 'wid': 'P2D', 'pred_prep_s': b1_prep_pred(w), 'n_queries': 1}
        if cpu_h() <= CPU_STOP_H:
            seq += 1
            c = _spawn('time', [json.dumps(spec)], pred)
            rec['cells'].append({'seq': seq, 'round': 0, 'wid': 'P2D', 'arm': 'V2b', **c})
            if c.get('status') == 'ok':
                rec['ref']['P2D'] = run('ref P2D', 'ref', ['P2D'], pred)
                if (rec['ref']['P2D'] or {}).get('status') == 'ok':
                    rec['gate']['V2b']['P2D'] = run('gate V2b P2D', 'gate', ['V2b', 'P2D'], pred)
    rec['_meta'].update({'finished_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                         'wall_s': round(time.time() - t_all, 1), 'cpu_h_total': cpu_h(), 'loadavg_end': loadavg(),
                         'mem_available_kB_end': R0.mem_available_kb()})
    save()
    print('done', sum(c.get('status') == 'ok' for c in rec['cells']), 'cells ok;',
          {c['wid']: c['status'] for c in rec['cells'] if c.get('status') != 'ok'}, f'cpu {cpu_h():.3f} h')
    return 0


if __name__ == '__main__':
    sys.exit(main())
