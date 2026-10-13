"""SYNTHETIC regression tests for EB-1 Amendment 2 (impossible sizes, silent truncation, retained backing allocations,
repeated-output lifetime). Uses temporary files; no development data.

Run: ../wpd_local_retrieval/build/venv/bin/python -B -m pytest -p no:cacheprovider test_eb1a2.py
"""
import sys
import tracemalloc
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import eb1a2_lib as L2  # noqa: E402


def test_impossible_size_raises_instead_of_silent_truncation(tmp_path):
    a = np.arange(1602 * 6, dtype=np.float64).reshape(1602, 2, 3)
    p = tmp_path / 'x.npy'
    np.save(p, a)
    assert np.load(p)[:1620].shape[0] == 1602                  # numpy slicing truncates silently
    with pytest.raises(L2.SizeError):
        L2.read_rows(p, [(0, 1620)])
    with pytest.raises(L2.SizeError):
        L2.read_rows(p, [(1600, 3)])
    with pytest.raises(L2.SizeError):
        L2.validate_size(1620, 1602, 'S3 D2')
    assert L2.validate_size(1600, 1602, 'S3 D2') == 1600


def test_read_rows_owned_exact_and_identical(tmp_path):
    a = np.random.default_rng(0).random((500, 7, 3)).astype(np.float32)
    p = tmp_path / 'y.npy'
    np.save(p, a)
    got = L2.read_rows(p, [(10, 40), (300, 100)])
    assert got.base is None and got.flags['OWNDATA'] and got.nbytes == 140 * 7 * 3 * 4
    assert np.array_equal(got, np.concatenate([a[10:50], a[300:400]]))


def test_old_concatenate_pattern_retains_backing_but_direct_reader_does_not(tmp_path):
    parts = [np.random.default_rng(i).random((450, 51, 3)).astype(np.float32) for i in range(4)]
    for i, x in enumerate(parts):
        np.save(tmp_path / f'D2_stored_{i}.npy', x)
        np.save(tmp_path / f'D2_box_{i}.npy', np.array([40.0, 40.0, 40.0]))
    old = np.concatenate(parts)[:500]                          # view: keeps the full concatenation alive
    assert old.base is not None and old.base.nbytes > old.nbytes
    full = sum(x.nbytes for x in parts)
    for n in (500, 1700):
        L2.read_d2_stored(tmp_path, n)
        tracemalloc.start()
        base = tracemalloc.get_traced_memory()[0]
        st, bx = L2.read_d2_stored(tmp_path, n)
        extra = tracemalloc.get_traced_memory()[1] - base - st.nbytes - bx.nbytes
        tracemalloc.stop()
        assert st.flags['OWNDATA'] and st.shape == (n, 51, 3)
        assert np.array_equal(st, np.concatenate(parts)[:n])
        # beyond the output only a constant .npy-header parsing allowance, never data-proportional backing
        assert extra < 256 * 1024 < full
    with pytest.raises(L2.SizeError):
        L2.read_d2_stored(tmp_path, 1801)


def test_repeated_output_lifetime_released():
    size = 2_000_000

    def fn():
        return np.ones(size)                                    # 16 MB output per repetition
    tracemalloc.start()
    out, walls, _ = L2.measure_reps(fn, max_reps=3, long_s=1e9)
    peak_new = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert len(walls) == 3 and peak_new < 1.3 * size * 8
    tracemalloc.start()
    prev = None
    for _ in range(3):                                         # Amendment 1 pattern: previous output alive
        cur = fn()
        prev, cur = cur, None
        cur = fn() if prev is not None else None
        prev = cur
    peak_old = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert peak_old > 1.8 * size * 8


def test_amendment2_design_sizes_valid_against_preflight():
    import json
    import eb1a2_run as R2
    pre = json.loads((Path(__file__).resolve().parent / 'results' / 'eb1a2_preflight.json').read_text())
    fr = pre['frames']
    idx = pre['index_sets']
    for c in R2.design():
        a = c['args']
        if c['block'] == 'S3':
            assert a['n'] <= {'D1': fr['D1_round_all'], 'D2': fr['D2_round_all_actual']}[a['d']]
        elif c['block'] == 'S2':
            assert a['n'] <= {'D1': fr['D1_stored_before_burn'], 'D2': fr['D2_stored_total_before_burn']}[a['d']]
        else:
            assert a['nq'] <= idx[a['d']]['n_q'] and a['nd'] <= idx[a['d']]['n_d']
    assert not any(c['block'] == 'S3' and c['args'] == {'d': 'D2', 'n': 1620} for c in R2.design())


def test_runs_blocks_and_a_batch_model_drop_order():
    import eb1a2_analyse as AZ
    import eb1a2_run as R2
    assert R2.runs([0, 1, 2, 5, 6, 9]) == [(0, 3), (5, 2), (9, 1)]
    rng = np.random.default_rng(1)
    rows = []
    for _ in range(8):                                          # SYNTHETIC: t = 2e-7 C + 1e-6 rows, no calls term
        C, r, k = rng.integers(10**6, 10**8), rng.integers(10**3, 10**5), rng.integers(5, 300)
        rows.append({'comparisons': int(C), 'call_db_rows': int(r), 'calls': int(k), 'n_win': 1,
                     'args': {'nq': 1, 'nd': 1}, 'cpu_s': 2e-7 * C + 1e-6 * r})
    m = AZ.a_batch_time_model(rows)
    assert 'c0' not in m['coef'] and abs(m['coef']['c1'] - 2e-7) < 1e-9
    for r in rows:                                              # force a negative calls coefficient
        r['cpu_s'] = 2e-7 * r['comparisons'] + 1e-6 * r['call_db_rows'] - 1e-5 * r['calls'] + 1e-3
    m = AZ.a_batch_time_model(rows)
    assert m['dropped'][:1] == ['c3'] and 'c3' not in m['coef']
