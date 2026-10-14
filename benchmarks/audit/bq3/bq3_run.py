#!/usr/bin/env python3
"""BQ-3 runner (docs/B2_DIRECT_P2_QUALIFICATION_PROTOCOL.md). SYNTHETIC_tiled, performance-only P2 workload.

Order: freeze checks → tests → reference (untimed) → gate (untimed) → 5 timing repetitions, each a fresh process.
Every process obeys the BQ-2 memory start rule; stops before any step once BQ-3 CPU exceeds 0.9 h.
Writes results/bq3_SYNTHETIC_measurements.json. Run: ../wpd_local_retrieval/build/venv/bin/python -B bq3_run.py
"""
import os

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'RAYON_NUM_THREADS'):
    os.environ[_v] = '1'

import datetime  # noqa: E402
import gc  # noqa: E402
import json  # noqa: E402
import resource  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
sys.path.insert(0, str(AUDIT / 'bq2'))
sys.path.insert(0, str(AUDIT / 'eb1'))
import bq2_run as B2R  # noqa: E402  (P2D workload definition, loader, hashing helpers; unchanged)
import eb1_run as R0  # noqa: E402

RES = HERE / 'results'
REF = HERE / 'build' / 'ref_P2D_SYNTHETIC.npy'
OUT = RES / 'bq3_SYNTHETIC_measurements.json'
N_REPS = 5
CPU_STOP_H = 0.9
W = B2R.workload('P2D')


def maxrss_kb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss


def os_threads():
    for line in open('/proc/self/status'):
        if line.startswith('Threads:'):
            return int(line.split()[1])


def mem_state():
    mi = {l.split(':')[0]: int(l.split()[1]) for l in open('/proc/meminfo') if ':' in l}
    return {'loadavg': open('/proc/loadavg').read().split()[:3], 'MemAvailable_kB': mi['MemAvailable'],
            'SwapFree_kB': mi['SwapFree']}


# --- children --------------------------------------------------------------------------------------------------------

def child_tests():
    p = subprocess.run([sys.executable, '-B', '-m', 'pytest', '-p', 'no:cacheprovider', '-q', 'test_bq3.py'],
                       capture_output=True, text=True, cwd=HERE)
    tail = p.stdout.strip().splitlines()[-1] if p.stdout.strip() else ''
    print('EB1RESULT ' + json.dumps({'returncode': p.returncode, 'summary': tail, 'pass': p.returncode == 0}))


def child_ref():
    import bq2_lib as Q
    REF.parent.mkdir(parents=True, exist_ok=True)
    X, di, wins = B2R.load(W)
    D = np.lib.format.open_memmap(REF, mode='w+', dtype=np.float64, shape=(3, len(wins), len(di)))
    for qi in range(3):
        for i, s in enumerate(wins):
            D[qi, i] = Q.B.ref_svd(X[qi, s:s + 8], X[di, s:s + 8])
    finite = bool(all(np.isfinite(D[qi]).all() for qi in range(3)))
    D.flush()
    del D
    print('EB1RESULT ' + json.dumps({'finite': finite, 'shape': [3, len(wins), len(di)], 'sha256': R0.sha256(REF)}))


def child_gate():
    import bq2_lib as Q
    X, di, wins = B2R.load(W)
    D_ref = np.load(REF, mmap_mode='r')
    prep = Q.B.prepare(X, di, wins)
    err, fin, n_bad, n_strict, n_near, nfb, hashes = 0.0, True, 0, 0, 0, 0, []
    for qi in range(3):
        out, D, f = Q.B.search('V2b', X, qi, di, wins, prep, keep_dist=True)
        Dr = np.asarray(D_ref[qi])
        nfb += f
        fin &= bool(np.isfinite(D).all())
        err = max(err, float(np.nanmax(np.abs(D - Dr))))
        for i in range(len(wins)):
            eq, st, nr = Q.B.rank_check(out[i], Dr[i])
            n_bad += not eq
            n_strict += st
            n_near += nr
        hashes.append(B2R.topk_hash(out))
        del out, D, Dr
    print('EB1RESULT ' + json.dumps({'candidate_finite': fin, 'max_abs_dist_err_A': err, 'query_windows': 3 * len(wins),
                                     'rank_nonequivalent': n_bad, 'rank_strict_identical': n_strict,
                                     'near_tie_query_windows': n_near, 'fallback_rows': nfb, 'topk_sha256': hashes,
                                     'peak_rss_kB': maxrss_kb()}))


def child_rep():
    import bq2_lib as Q
    B = Q.B
    res = {'start': mem_state(), 'threads_env': {v: os.environ.get(v) for v in ('OMP_NUM_THREADS',
                                                                                'OPENBLAS_NUM_THREADS')},
           'threadpool': B2R.threads_info()}
    rng = np.random.default_rng(1)
    Xs = B._chains(rng, 24, L=12)                                 # identical SYNTHETIC warm-up (BQ-1/BQ-2)
    ws = list(range(5))
    B.search('V2b', Xs, 0, np.arange(1, 24), ws, B.prepare(Xs, np.arange(1, 24), ws))
    del Xs
    gc.collect()
    t0 = time.perf_counter()
    X, di, wins = B2R.load(W)
    res.update({'build_rows_wall_s': time.perf_counter() - t0, 'rss_after_load_kB': maxrss_kb(),
                'resident_bytes': int(X.nbytes)})

    def timed(fn):
        a, c = time.perf_counter(), time.process_time()
        v = fn()
        return v, time.perf_counter() - a, time.process_time() - c
    s0 = time.perf_counter()                                     # actual first-query sequence: P1 then Q1
    prep, p1w, p1c = timed(lambda: B.prepare(X, di, wins))
    (out, _, f1), q1w, q1c = timed(lambda: B.search('V2b', X, 0, di, wins, prep))
    span = time.perf_counter() - s0
    res.update({'P1_wall_s': p1w, 'P1_cpu_s': p1c, 'rss_after_P1_Q1_kB': maxrss_kb(),
                'prepared_bytes': int(prep[0].nbytes + prep[1].nbytes), 'first_query_span_s': span})
    qw, qc, hashes, nfb = [q1w], [q1c], [B2R.topk_hash(out)], f1
    for qi in (1, 2):                                            # subsequent queries on the same prepared database
        del out
        gc.collect()
        (out, _, f), w_, c_ = timed(lambda: B.search('V2b', X, qi, di, wins, prep))
        qw.append(w_)
        qc.append(c_)
        hashes.append(B2R.topk_hash(out))
        nfb += f
    res.update({'query_wall_s': qw, 'query_cpu_s': qc, 'topk_sha256': hashes, 'fallback_rows': nfb,
                'retained_bytes': int(out.nbytes), 'rss_after_queries_kB': maxrss_kb(), 'os_threads': os_threads()})
    del out, prep
    gc.collect()
    pw, pc = [], []
    for _ in range(2):                                           # repeated preparations, each released
        p, w_, c_ = timed(lambda: B.prepare(X, di, wins))
        pw.append(w_)
        pc.append(c_)
        del p
        gc.collect()
    res.update({'P_repeat_wall_s': pw, 'P_repeat_cpu_s': pc, 'peak_rss_kB': maxrss_kb(), 'end': mem_state()})
    print('EB1RESULT ' + json.dumps(res))


# --- parent ----------------------------------------------------------------------------------------------------------

def _spawn(mode, args, pred_bytes):
    """Fresh subprocess of this script under RLIMIT_AS, started only if the BQ-2 memory rule holds (wait ≤ 10 min)."""
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
    if len(sys.argv) >= 2 and sys.argv[1] in ('tests', 'ref', 'gate', 'rep'):
        return {'tests': child_tests, 'ref': child_ref, 'gate': child_gate, 'rep': child_rep}[sys.argv[1]]()
    for rec_name in ('bq3_protocol_freeze.json', 'bq3_pre_measurement_record.json'):
        for rel, h in json.loads((RES / rec_name).read_text())['sha256'].items():
            if R0.sha256(R0.REPO / rel) != h:
                raise SystemExit(f'STOP: {rel} changed ({rec_name})')
    for rec in (AUDIT / 'b2opt' / 'results' / 'b2opt_pre_measurement_record.json',
                AUDIT / 'bq2' / 'results' / 'bq2_pre_measurement_record.json'):
        for rel, h in json.loads(rec.read_text())['sha256'].items():
            if rel.endswith(('b2opt_lib.py', 'bq2_lib.py', 'bq2_run.py')) and R0.sha256(R0.REPO / rel) != h:
                raise SystemExit(f'STOP: {rel} differs from {rec.name}')
    import bq2_lib as Q
    L, nwin = B2R.dims(W)
    base = Q.predicted_peak_bytes(3 + W['nd'] + 3706, L, W['nd'], nwin)
    rec = {'_meta': {'producer': 'bq3_run.py', 'label': 'SYNTHETIC_tiled performance-only',
                     'started_utc': datetime.datetime.now(datetime.UTC).isoformat(), 'n_reps': N_REPS,
                     'cpu_cap_h': 1.0, 'start_state': mem_state()}, 'steps': [], 'reps': []}

    def cpu_h():
        return R0.children_cpu() / 3600

    def save():
        rec['_meta']['cpu_h_so_far'] = cpu_h()
        OUT.write_text(json.dumps(rec, indent=1) + '\n')

    def run(name, mode, pred, extra=()):
        if cpu_h() > CPU_STOP_H:
            rec['stopped'] = f'CPU budget reached before {name}'
            return None
        r = _spawn(mode, list(extra), pred)
        rec['steps'].append({'step': name, 'status': r.get('status'), 'cpu_h_after': cpu_h(),
                             'started_utc': r.get('started_utc')})
        save()
        return r
    rec['tests'] = run('tests', 'tests', 300 * 1024 ** 2)
    if not rec['tests'] or not rec['tests'].get('pass'):
        rec.setdefault('stopped', 'failure-path tests failed')
        save()
        return 2
    ref_pred = Q.predicted_peak_bytes(3 + W['nd'] + 3706, L, 0, nwin) + nwin * W['nd'] * 8
    rec['ref'] = run('reference', 'ref', ref_pred)
    rec['gate'] = run('gate', 'gate', base + 2 * nwin * W['nd'] * 8) if (rec['ref'] or {}).get('status') == 'ok' else None
    if not rec['gate'] or rec['gate'].get('status') != 'ok':
        rec.setdefault('stopped', 'gate missing or failed; repetitions not run')
        save()
        return 2
    for i in range(N_REPS):
        r = run(f'rep {i}', 'rep', base)
        if r is None:
            break
        rec['reps'].append({'rep': i, **r})
        save()
    rec['_meta'].update({'finished_utc': datetime.datetime.now(datetime.UTC).isoformat(), 'cpu_h_total': cpu_h(),
                         'end_state': mem_state()})
    save()
    print('done', [r.get('status') for r in rec['reps']], f'cpu {cpu_h():.3f} h')
    return 0


if __name__ == '__main__':
    sys.exit(main())
