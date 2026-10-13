#!/usr/bin/env python3
"""BQ-1 runner (docs/B2_SINGLE_QUERY_OPTIMIZATION_PROTOCOL.md + AMENDMENT_1). Development data only.

Order: freeze checks → tests → reference check → V0 gate → [V1 gate → V1 timing] → [V2 gate → V2 timing] →
[V2b gate → V2b timing]. Timing rounds interleave V0 and the arm in fresh subprocesses (seeded order, recorded).
Stops before any step once the experiment's cumulative CPU (parent + children) exceeds 0.9 h.
Writes results/b2opt_measurements.json. Run: ../wpd_local_retrieval/build/venv/bin/python -B b2opt_run.py
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
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(EB1))
import eb1_run as R0  # noqa: E402  (sha256, memory helpers, limits, CPU accounting)

RES = HERE / 'results'
REF = HERE / 'build' / 'ref'
OUT = RES / 'b2opt_measurements.json'
WL = json.loads((RES / 'b2opt_workloads.json').read_text())
CPU_STOP_H = 0.9
ROUNDS = 3
ORDER_SEED = 20261017
SAMPLE_SEED = 20261014
ARMS = ('V1', 'V2', 'V2b')


# --- data ------------------------------------------------------------------------------------------------------------

def workload(wid):
    return next(w for w in WL['workloads'] if w['id'] == wid)


def load(w):
    """Owned array: rows 0–2 = the three query frames, rows 3… = first n_D database frames. Returns (X, di, windows)."""
    import eb1a1_run as R1
    import eb1a2_lib as L2
    src = EB1 / 'build' / ('cache_a1' if w['d'] == 'D1p' else 'cache')
    meta = json.loads((src / f"{w['d']}_meta.json").read_text())
    L2.validate_size(w['nd'], len(meta['di']), 'n_D')
    rows = [meta['qi'][q] for q in w['queries']] + meta['di'][:w['nd']]
    runs = []
    for i in rows:
        if runs and runs[-1][0] + runs[-1][1] == i:
            runs[-1][1] += 1
        else:
            runs.append([i, 1])
    X = L2.read_rows(src / f"{w['d']}_round.npy", [tuple(r) for r in runs])
    wins = R1.windows_for(X.shape[1], w['nwin'])
    assert len(wins) == min(w['nwin'], X.shape[1] - 7)
    return X, np.arange(3, 3 + w['nd']), wins


def topk_hash(out):
    return hashlib.sha256(np.ascontiguousarray(out, np.int64).tobytes()).hexdigest()


def maxrss_kb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


def threads_info():
    try:
        import threadpoolctl
        return [{'api': x.get('internal_api'), 'num_threads': x.get('num_threads')}
                for x in threadpoolctl.threadpool_info()]
    except Exception as e:                                      # recorded, never silently ignored
        return [{'error': repr(e)}]


# --- children --------------------------------------------------------------------------------------------------------

def child_tests():
    p = subprocess.run([sys.executable, '-B', '-m', 'pytest', '-p', 'no:cacheprovider', '-q', 'test_b2opt.py'],
                       capture_output=True, text=True, cwd=HERE)
    tail = p.stdout.strip().splitlines()[-1] if p.stdout.strip() else ''
    print('EB1RESULT ' + json.dumps({'returncode': p.returncode, 'summary': tail, 'pass': p.returncode == 0}))


def child_refcheck():
    import b2opt_lib as B
    REF.mkdir(parents=True, exist_ok=True)
    rows = {}
    ok = True
    for w in WL['workloads']:
        X, di, wins = load(w)
        D = np.empty((3, len(wins), len(di)))
        for qi in range(3):
            for i, s in enumerate(wins):
                D[qi, i] = B.ref_svd(X[qi, s:s + B.W], X[di, s:s + B.W])
        np.save(REF / f"{w['id']}.npy", D)
        rng = np.random.default_rng(SAMPLE_SEED)
        e1 = e2 = 0.0
        for _ in range(2000):
            qi, i, j = rng.integers(3), rng.integers(len(wins)), rng.integers(len(di))
            s = wins[i]
            a, b = X[qi, s:s + B.W], X[di[j], s:s + B.W]
            e1 = max(e1, abs(B.ref_kabsch_pair(a, b) - D[qi, i, j]))
            e2 = max(e2, abs(B.ref_qcprot_pair(a, b) - D[qi, i, j]))
        finite = bool(np.isfinite(D).all())
        rows[w['id']] = {'pairs': int(D.size), 'max_vs_kabsch_explicit': e1, 'max_vs_mdanalysis_qcprot': e2,
                         'finite': finite, 'sha256': hashlib.sha256(D.tobytes()).hexdigest()}
        ok &= finite and e1 <= B.TOL_D and e2 <= B.TOL_D
    print('EB1RESULT ' + json.dumps({'workloads': rows, 'pass': bool(ok)}))


def child_gate():
    import b2opt_lib as B
    arm = sys.argv[2]
    res = {'arm': arm, 'workloads': {}}
    ok = True
    for w in WL['workloads']:
        X, di, wins = load(w)
        D_ref = np.load(REF / f"{w['id']}.npy")
        prep = None if arm == 'V0' else B.prepare(X, di, wins)
        err, n_eq_fail, n_strict, n_near, nfb, hashes = 0.0, 0, 0, 0, 0, []
        for qi in range(3):
            out, D, f = B.search(arm, X, qi, di, wins, prep, keep_dist=True)
            nfb += f
            err = max(err, float(np.max(np.abs(D - D_ref[qi]))) if np.isfinite(D).all() else float('inf'))
            for i in range(len(wins)):
                eq, st, nr = B.rank_check(out[i], D_ref[qi, i])
                n_eq_fail += not eq
                n_strict += st
                n_near += nr
            hashes.append(topk_hash(out))
            del out, D
        qw = 3 * len(wins)
        res['workloads'][w['id']] = {'max_abs_dist_err_A': err, 'query_windows': qw, 'rank_nonequivalent': n_eq_fail,
                                     'rank_strict_identical': n_strict, 'near_tie_query_windows': n_near,
                                     'fallback_rows': nfb, 'fallback_fraction': nfb / (qw * len(di)),
                                     'topk_sha256': hashes}
        ok &= err <= B.TOL_D and n_eq_fail == 0
        del prep
        gc.collect()
    res['synthetic'] = B.synthetic_suite(arm)
    res['alpha'] = B.calibrate_alpha(arm)
    res['pass'] = bool(ok and res['synthetic']['pass'])
    print('EB1RESULT ' + json.dumps(res))


def child_time():
    import b2opt_lib as B
    spec = json.loads(sys.argv[2])
    arm, w = spec['arm'], workload(spec['wid'])
    rng = np.random.default_rng(1)
    Xs = B._chains(rng, 24, L=12)                               # identical SYNTHETIC warm-up for every arm
    ws = list(range(5))
    B.search(arm, Xs, 0, np.arange(1, 24), ws, None if arm == 'V0' else B.prepare(Xs, np.arange(1, 24), ws))
    del Xs
    gc.collect()
    rss0 = maxrss_kb()
    X, di, wins = load(w)
    res = {'threads': threads_info(), 'rss_after_warmup_kB': rss0, 'rss_after_load_kB': maxrss_kb(),
           'resident_bytes': int(X.nbytes), 'n_win': len(wins), 'n_d': len(di)}
    prep = None
    if arm != 'V0':
        t0, c0 = time.perf_counter(), time.process_time()
        prep = B.prepare(X, di, wins)
        res.update({'prep_wall_s': time.perf_counter() - t0, 'prep_cpu_s': time.process_time() - c0,
                    'prepared_bytes': int(prep[0].nbytes + prep[1].nbytes)})
    else:
        res.update({'prep_wall_s': 0.0, 'prep_cpu_s': 0.0, 'prepared_bytes': 0})
    res['rss_after_prep_kB'] = maxrss_kb()
    walls, cpus, hashes, nfb = [], [], [], 0
    out = None
    for qi in range(3):
        if out is not None:
            del out
            out = None
            gc.collect()
        t0, c0 = time.perf_counter(), time.process_time()
        out, _, f = B.search(arm, X, qi, di, wins, prep)
        walls.append(time.perf_counter() - t0)
        cpus.append(time.process_time() - c0)
        nfb += f
        hashes.append(topk_hash(out))
    res.update({'query_wall_s': walls, 'query_cpu_s': cpus, 'topk_sha256': hashes, 'fallback_rows': nfb,
                'retained_bytes': int(out.nbytes), 'rss_end_kB': maxrss_kb()})
    print('EB1RESULT ' + json.dumps(res))


# --- parent ----------------------------------------------------------------------------------------------------------

def _spawn(mode, arg=None, timeout=3600):
    def limits():
        resource.setrlimit(resource.RLIMIT_AS, (R0.AS_LIMIT, R0.AS_LIMIT))
    waited = 0
    while R0.mem_available_kb() < R0.MEM_MIN_KB and waited < R0.MEM_WAIT_S:
        time.sleep(15)
        waited += 15
    mem0 = R0.mem_available_kb()
    if mem0 < R0.MEM_MIN_KB:
        return {'status': 'skipped_low_memory', 'mem_available_kB': mem0}
    args = [sys.executable, '-B', __file__, mode] + ([arg] if arg is not None else [])
    t_start = datetime.datetime.now(datetime.UTC).isoformat()
    p = subprocess.run(args, capture_output=True, text=True, preexec_fn=limits, timeout=timeout, cwd=HERE)
    line = [x for x in p.stdout.splitlines() if x.startswith('EB1RESULT ')]
    if p.returncode != 0 or not line:
        return {'status': 'failed', 'returncode': p.returncode, 'stderr': p.stderr[-1500:], 'started_utc': t_start}
    out = json.loads(line[-1][10:])
    out.update({'status': 'ok', 'mem_available_kB_at_start': mem0, 'started_utc': t_start})
    return out


def main():
    if len(sys.argv) >= 2 and sys.argv[1] in ('tests', 'refcheck', 'gate', 'time'):
        return {'tests': child_tests, 'refcheck': child_refcheck, 'gate': child_gate, 'time': child_time}[sys.argv[1]]()
    for rec_name in ('b2opt_protocol_freeze.json', 'b2opt_pre_measurement_record.json'):
        fr = json.loads((RES / rec_name).read_text())
        for rel, h in fr['sha256'].items():
            if R0.sha256(R0.REPO / rel) != h:
                raise SystemExit(f'STOP: {rel} changed after {rec_name}')
    orig = json.loads((EB1 / 'results' / 'eb1_measurements.json').read_text())
    for name, h in orig['prep_cache_sha256'].items():
        if R0.sha256(EB1 / 'build' / 'cache' / name) != h:
            raise SystemExit(f'STOP: cache {name} differs from EB-1')
    a1 = json.loads((EB1 / 'results' / 'eb1a1_measurements.json').read_text())
    for name, h in a1['cache_a1_sha256'].items():
        if R0.sha256(EB1 / 'build' / 'cache_a1' / name) != h:
            raise SystemExit(f'STOP: cache_a1 {name} differs from Amendment 1')
    t_all = time.time()
    rec = {'_meta': {'producer': 'b2opt_run.py', 'started_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                     'cpu_cap_h': 1.0, 'cpu_stop_h': CPU_STOP_H, 'rounds': ROUNDS, 'order_seed': ORDER_SEED},
           'steps': [], 'arms': {}}

    def cpu_h():
        return R0.children_cpu() / 3600

    def save():
        rec['_meta']['cpu_h_so_far'] = cpu_h()
        OUT.write_text(json.dumps(rec, indent=1) + '\n')

    def step(name, mode, arg=None):
        if cpu_h() > CPU_STOP_H:
            rec['stopped'] = f'CPU budget reached before {name}'
            save()
            return None
        r = _spawn(mode, arg)
        rec['steps'].append({'step': name, 'status': r.get('status'), 'cpu_h_after': cpu_h()})
        save()
        return r
    rec['tests'] = step('tests', 'tests')
    if not rec['tests'] or not rec['tests'].get('pass'):
        rec.setdefault('stopped', 'tests failed')
        save()
        return 2
    rec['refcheck'] = step('refcheck', 'refcheck')
    if not rec['refcheck'] or not rec['refcheck'].get('pass'):
        rec.setdefault('stopped', 'reference check failed')
        save()
        return 2
    rec['gate_V0'] = step('gate V0', 'gate', 'V0')
    if not rec['gate_V0'] or not rec['gate_V0'].get('pass'):
        rec.setdefault('stopped', 'V0 gate failed (baseline and reference disagree)')
        save()
        return 2
    rng = np.random.default_rng(ORDER_SEED)
    seq = 0
    for arm in ARMS:
        A = rec['arms'][arm] = {'status': 'untested'}
        g = step(f'gate {arm}', 'gate', arm)
        if g is None:
            break
        A['gate'] = g
        if g.get('status') != 'ok' or not g.get('pass'):
            A['status'] = 'rejected_at_gate'
            save()
            continue
        A['status'] = 'timing_incomplete'
        A['cells'] = []
        wids = [w['id'] for w in WL['workloads']]
        stopped = False
        for r in range(ROUNDS):
            for wid in rng.permutation(wids).tolist():
                for a in rng.permutation(['V0', arm]).tolist():
                    if cpu_h() > CPU_STOP_H:
                        rec['stopped'] = f'CPU budget reached during {arm} timing'
                        stopped = True
                        break
                    seq += 1
                    c = _spawn('time', json.dumps({'arm': a, 'wid': wid}))
                    A['cells'].append({'seq': seq, 'round': r, 'wid': wid, 'arm': a, **c})
                    save()
                if stopped:
                    break
            if stopped:
                break
        if stopped:
            save()
            break
        A['status'] = 'timed'
        save()
    rec['_meta'].update({'finished_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                         'wall_s': round(time.time() - t_all, 1), 'cpu_h_total': cpu_h()})
    save()
    print('done', {a: v['status'] for a, v in rec['arms'].items()}, f'cpu {cpu_h():.3f} h')
    return 0


if __name__ == '__main__':
    sys.exit(main())
