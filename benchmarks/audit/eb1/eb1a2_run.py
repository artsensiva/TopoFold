#!/usr/bin/env python3
"""EB-1 Amendment 2 runner (docs/EB1_AMENDMENT_2.md). Development data only; outputs separate from EB-1 and Amendment 1.

Order: verify the Amendment 2 freeze, the Amendment 1 freeze and the reused caches → loader correctness (separate
subprocess; new direct-read loaders vs the earlier load paths; stop on failure) → affected cells only, each in a fresh
subprocess: S2 and S3 at the Amendment 1 sizes with S3 D2 1620 → 1600, affected S4 cells, fresh A confirmation
workloads FC1–FC4. Inputs are read directly into owned arrays (eb1a2_lib.read_rows / read_d2_stored); repeated S4
measurements release each output first (eb1a2_lib.measure_reps). Writes results/eb1a2_measurements.json.
Run: ../wpd_local_retrieval/build/venv/bin/python -B eb1a2_run.py
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
import eb1_run as R0  # noqa: E402  (sha256, mem helpers, limits)
import eb1a1_run as R1  # noqa: E402  (Amendment 1 design and windows_for, unchanged)

RES = HERE / 'results'
C0 = HERE / 'build' / 'cache'
C1 = HERE / 'build' / 'cache_a1'
OUT = RES / 'eb1a2_measurements.json'
FREEZE = RES / 'eb1a2_freeze.json'
PREFLIGHT = RES / 'eb1a2_preflight.json'
COEF = json.loads((RES / 'eb1a1_temp_calibration.json').read_text())['coef']
W = 8
CPU_STOP_H = 1.8
S3_D2_REPLACEMENT = {1620: 1600}
FRESH = [('FC1', 'D1', 48, 1015, 10), ('FC2', 'D1', 160, 750, 30), ('FC3', 'D1p', 384, 2693, 60),
         ('FC4', 'D2', 200, 405, 44)]


def prior_cpu_h():
    a1 = json.loads((RES / 'eb1a1_measurements.json').read_text())['_meta']
    return a1['aggregate_cpu_h_including_original']


def design():
    """Affected Amendment 1 cells (from the frozen preflight classification) with the S3 D2 size replaced, plus FC1–FC4."""
    pre = json.loads(PREFLIGHT.read_text())
    affected = [x for x in pre['amendment1_memory_classification'] if x['memory_affected']]
    a1 = R1.design()
    cells = []
    for c in a1:
        hit = [x for x in affected if x['block'] == c['block'] and x['args'] == c['args']]
        if not hit:
            continue
        a = dict(c['args'])
        if c['block'] == 'S3' and a['d'] == 'D2' and a['n'] in S3_D2_REPLACEMENT:
            a['n'] = S3_D2_REPLACEMENT[a['n']]
        cells.append({'block': c['block'], 'args': a, 'rep': c['rep'], 'amends': c['args']})
    for cid, d, nq, nd, nw in FRESH:
        cells.append({'block': 'S4', 'args': {'m': 'A', 'mode': 'batch', 'd': d, 'nq': nq, 'nd': nd, 'nwin': nw,
                                              'role': cid}, 'rep': 0, 'amends': None})
    return cells


# --- loaders (direct reads into owned arrays) -----------------------------------------------------------------------

def runs(idx):
    """Contiguous (start, n) blocks of an index list, in order."""
    out = []
    for i in idx:
        if out and out[-1][0] + out[-1][1] == i:
            out[-1][1] += 1
        else:
            out.append([int(i), 1])
    return [tuple(x) for x in out]


def _dataset(d):
    src = {'D1': C0, 'D2': C0, 'D1p': C1}[d]
    return src, json.loads((src / f'{d}_meta.json').read_text())


def load_inputs(blk, a):
    import eb1a2_lib as L2
    if blk == 'S2':
        if a['d'] == 'D1':
            return {'st': L2.read_rows(C0 / 'D1_stored.npy', [(0, a['n'])]),
                    'bx': L2.read_rows(C0 / 'D1_box.npy', [(0, a['n'])])}
        st, bx = L2.read_d2_stored(C0, a['n'])
        return {'st': st, 'bx': bx}
    if blk == 'S3':
        return {'X': L2.read_rows(C0 / f"{a['d']}_round_all.npy", [(0, a['n'])])}
    src, meta = _dataset(a['d'])
    L2.validate_size(a['nq'], len(meta['qi']), 'n_Q')
    L2.validate_size(a['nd'], len(meta['di']), 'n_D')
    rr = runs(meta['qi'][:a['nq']] + meta['di'][:a['nd']])
    if a['m'] == 'A':
        return {'kap': L2.read_rows(src / f"{a['d']}_kap.npy", rr), 'tau': L2.read_rows(src / f"{a['d']}_tau.npy", rr)}
    return {'X': L2.read_rows(src / f"{a['d']}_round.npy", rr)}


def load_old(blk, a):
    """The Amendment 1 load paths, used only in the correctness subprocess."""
    import eb1a1_lib as A1
    if blk == 'S2':
        if a['d'] == 'D1':
            return {'st': A1.owned_subset(np.load(C0 / 'D1_stored.npy', mmap_mode='r'), np.arange(a['n'])),
                    'bx': A1.owned_subset(np.load(C0 / 'D1_box.npy', mmap_mode='r'), np.arange(a['n']))}
        parts = [np.load(C0 / f'D2_stored_{i}.npy') for i in range(4)]
        boxes = [np.repeat(np.load(C0 / f'D2_box_{i}.npy')[None], len(p), 0) for i, p in enumerate(parts)]
        return {'st': np.ascontiguousarray(np.concatenate(parts)[:a['n']]),
                'bx': np.ascontiguousarray(np.concatenate(boxes)[:a['n']])}
    if blk == 'S3':
        return {'X': A1.owned_subset(np.load(C0 / f"{a['d']}_round_all.npy", mmap_mode='r'), np.arange(a['n']))}
    src, meta = _dataset(a['d'])
    rows = np.concatenate([np.array(meta['qi'])[:a['nq']], np.array(meta['di'])[:a['nd']]])
    names = ('kap', 'tau') if a['m'] == 'A' else ('round',)
    return {('X' if k == 'round' else k): A1.owned_subset(np.load(src / f"{a['d']}_{k}.npy", mmap_mode='r'), rows)
            for k in names}


# --- child functions ------------------------------------------------------------------------------------------------

def run_cell(spec):
    import eb1_lib as E
    import eb1a1_lib as A1
    import eb1a2_lib as L2
    blk, a = spec['block'], spec['args']
    inp = load_inputs(blk, a)
    res = {'inputs_owned': all(v.flags['OWNDATA'] and v.base is None for v in inp.values()),
           'input_bytes': int(sum(v.nbytes for v in inp.values()))}
    if blk == 'S2':
        L = inp['st'].shape[1]
        t0, c0 = time.perf_counter(), time.process_time()
        E.reconstruct(inp['st'], inp['bx'])                  # outputs discarded at once (as in Amendment 1)
        res.update({'wall_s': time.perf_counter() - t0, 'cpu_s': time.process_time() - c0,
                    'x': a['n'] * L, 'P': a['n'] * L * (12 + 24 + 12 + 24)})
        return res
    if blk == 'S3':
        X = inp['X']
        L = X.shape[1]
        E.features(X[:5])
        t0, c0 = time.perf_counter(), time.process_time()
        E.features(X)
        res.update({'wall_s': time.perf_counter() - t0, 'cpu_s': time.process_time() - c0,
                    'x': a['n'] * L, 'P': a['n'] * L * (24 + 12) + a['n'] * (2 * L - 5) * 8})
        return res
    L = inp['kap'].shape[1] + 2 if a['m'] == 'A' else inp['X'].shape[1]
    qi, di = np.arange(a['nq']), np.arange(a['nq'], a['nq'] + a['nd'])
    wins = R1.windows_for(L, a['nwin'])
    batch = 1 if a['mode'] == 'single' else A1.BATCH_LIMIT
    A1.search_instrumented(a['m'], W, wins[:1], inp, qi[:min(8, len(qi))], di, batch)      # warm-up, excluded
    q = qi if a['mode'] == 'batch' else qi[:1]
    parts = []

    def one():
        o, tp = A1.search_instrumented(a['m'], W, wins, inp, q, di, batch)
        parts.append(tp)
        return o
    out, walls, cpus = L2.measure_reps(one, max_reps=3, long_s=5.0)
    wk = A1.work(a['nq'], a['nd'], len(wins))
    resident = sum(v.nbytes for v in inp.values())
    retained = sum(v.nbytes for v in out.values())
    res.update({'wall_s': float(np.median(walls)), 'cpu_s': float(np.median(cpus)), 'repeats': len(walls),
                'repeat_wall_s': walls, 'repeat_cpu_s': cpus,
                'decomposition_s': {k: float(np.median([p[k] for p in parts])) for k in parts[0]},
                'n_win': len(wins), 'L': L, **wk, 'resident_bytes': resident, 'retained_bytes': retained,
                'temp_bytes_pred': A1.temp_bytes(a['m'], wk['b'], a['nd'], W, COEF)})
    res['P'] = resident + retained + res['temp_bytes_pred']
    return res


def child_cell():
    out = run_cell(json.loads(sys.argv[2]))
    out['maxrss_kB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print('EB1RESULT ' + json.dumps(out))


def loadcheck():
    """Mandatory: every new-loader input of the design equals the earlier load path (identical dtype, shape, values),
    is owned, and an FC search through measure_reps equals the EB-1 harness search."""
    import eb1_lib as E
    import eb1a1_lib as A1
    import eb1a2_lib as L2
    seen, rows = set(), []
    for c in design():
        key = json.dumps([c['block'], {k: v for k, v in c['args'].items() if k not in ('role', 'mode', 'nwin')}],
                         sort_keys=True)
        if key in seen:
            continue
        seen.add(key)
        new, old = load_inputs(c['block'], c['args']), load_old(c['block'], c['args'])
        same = (new.keys() == old.keys() and all(new[k].dtype == old[k].dtype and new[k].shape == old[k].shape
                                                 and np.array_equal(new[k], old[k]) for k in new))
        owned = all(v.flags['OWNDATA'] and v.base is None for v in new.values())
        rows.append({'block': c['block'], 'args': c['args'], 'identical': bool(same), 'owned': bool(owned)})
        del new, old
    a = {'m': 'A', 'mode': 'batch', 'd': 'D1', 'nq': 48, 'nd': 1015, 'nwin': 10}
    inp = load_inputs('S4', a)
    src, meta = _dataset('D1')
    X, K, T = np.load(src / 'D1_round.npy'), np.load(src / 'D1_kap.npy'), np.load(src / 'D1_tau.npy')
    wins = R1.windows_for(X.shape[1], a['nwin'])
    got, _, _ = L2.measure_reps(lambda: A1.search_instrumented('A', W, wins, inp, np.arange(48),
                                                                np.arange(48, 48 + 1015))[0], max_reps=2)
    ref = E.search('A', W, wins, X, K, T, np.array(meta['qi'])[:48], np.array(meta['di'])[:1015])
    # harness indices are database positions; both paths index the database 0..n_D−1
    fc_same = all(np.array_equal(got[s], ref[s]) for s in wins)
    out = {'loads': rows, 'FC1_search_identical_to_harness': bool(fc_same),
           'pass': all(r['identical'] and r['owned'] for r in rows) and bool(fc_same)}
    print('EB1RESULT ' + json.dumps(out))


def main():
    if len(sys.argv) >= 2 and sys.argv[1] in ('cell', 'loadcheck'):
        return {'cell': child_cell, 'loadcheck': loadcheck}[sys.argv[1]]()
    for fpath, what in ((FREEZE, 'Amendment 2'), (RES / 'eb1a1_freeze.json', 'Amendment 1')):
        fr = json.loads(fpath.read_text())
        for rel, h in fr['sha256'].items():
            if R0.sha256(R0.REPO / rel) != h:
                raise SystemExit(f'STOP: {rel} changed after the {what} freeze')
    for rel, h in json.loads(FREEZE.read_text())['preserved'].items():
        if R0.sha256(R0.REPO / rel) != h:
            raise SystemExit(f'STOP: preserved record {rel} changed')
    orig = json.loads((RES / 'eb1_measurements.json').read_text())
    for name, h in orig['prep_cache_sha256'].items():
        if R0.sha256(C0 / name) != h:
            raise SystemExit(f'STOP: reused cache {name} differs from the original run')
    a1m = json.loads((RES / 'eb1a1_measurements.json').read_text())
    for name, h in a1m['cache_a1_sha256'].items():
        if R0.sha256(C1 / name) != h:
            raise SystemExit(f'STOP: reused cache_a1 {name} differs from Amendment 1')
    prior = prior_cpu_h()
    t0 = time.time()
    rec = {'_meta': {'producer': 'eb1a2_run.py',
                     'amendment_sha256': json.loads(FREEZE.read_text())['sha256']['docs/EB1_AMENDMENT_2.md'],
                     'started_utc': datetime.datetime.now(datetime.UTC).isoformat(), 'prior_cpu_h': prior,
                     'reused': 'EB-1 and Amendment 1 caches (hash-verified); Amendment 1 times; unaffected memory cells'}}

    def cpu_h():
        return prior + R0.children_cpu() / 3600
    rec['loadcheck'] = _spawn('loadcheck')
    OUT.write_text(json.dumps(rec, indent=1) + '\n')
    if rec['loadcheck'].get('status') != 'ok' or not rec['loadcheck'].get('pass'):
        rec['stopped'] = 'mandatory loader correctness failure'
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
    rec['_meta'].update({'finished_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                         'wall_s': round(time.time() - t0, 1), 'amendment_cpu_h': R0.children_cpu() / 3600,
                         'aggregate_cpu_h_including_prior': cpu_h(), 'n_cells_planned': len(cells),
                         'n_cells_ok': sum(c.get('status') == 'ok' for c in rec['cells'])})
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
