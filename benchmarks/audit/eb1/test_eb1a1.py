"""SYNTHETIC regression tests for EB-1 Amendment 1 accounting (no development data).

Run: ../wpd_local_retrieval/build/venv/bin/python -B -m pytest -p no:cacheprovider test_eb1a1.py
"""
import json
import sys
import tracemalloc
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import eb1_lib as E  # noqa: E402
import eb1a1_lib as A1  # noqa: E402

COEF = json.loads((HERE / 'results' / 'eb1a1_temp_calibration.json').read_text())['coef']
R = np.random.default_rng(1)


def feats(method, n, L=20):
    if method == 'A':
        return {'kap': R.uniform(0.5, 2, (n, L - 2)), 'tau': R.uniform(-3, 3, (n, L - 3))}
    return {'X': R.normal(size=(n, L, 3)) * 5}


def test_work_counts_and_batch_saturation():
    assert A1.effective_batch(1) == 1 and A1.effective_batch(64) == 64 and A1.effective_batch(500) == 128
    w = A1.work(256, 1000, 20)
    assert w == {'comparisons': 256 * 1000 * 20, 'calls': 40, 'call_db_rows': 40 * 1000, 'b': 128}
    assert A1.work(129, 10, 1)['calls'] == 2 and A1.work(128, 10, 1)['calls'] == 1


def test_retained_outputs_match_actual_arrays():
    f = feats('A', 300)
    qi, di = np.arange(37), np.arange(37, 300)
    out, _ = A1.search_instrumented('A', 8, list(range(5)), f, qi, di)
    assert sum(v.nbytes for v in out.values()) == A1.retained_bytes(37, 5)


@pytest.mark.parametrize('method', ['A', 'B2'])
def test_transient_bytes_predicted_by_calibration_at_unseen_sizes(method):
    for b, nd in ((40, 1200), (128, 2500)):
        f = feats(method, b + nd)
        qi, di = np.arange(b), np.arange(b, b + nd)
        A1.search_instrumented(method, 8, [4], f, qi, di)
        tracemalloc.start()
        base = tracemalloc.get_traced_memory()[0]
        out, _ = A1.search_instrumented(method, 8, [4], f, qi, di)
        peak = tracemalloc.get_traced_memory()[1] - base - out[4].nbytes
        tracemalloc.stop()
        pred = A1.temp_bytes(method, b, nd, 8, COEF)
        assert abs(pred - peak) / peak < 0.15


@pytest.mark.parametrize('method', ['A', 'B2'])
def test_instrumented_search_identical_to_harness(method):
    f = feats(method, 400)
    qi, di = np.arange(150), np.arange(150, 400)
    X = f.get('X', np.zeros((400, 20, 3)))
    kap, tau = f.get('kap', np.zeros((400, 18))), f.get('tau', np.zeros((400, 17)))
    ref = E.search(method, 8, [0, 5, 9], X, kap, tau, qi, di)
    got, t = A1.search_instrumented(method, 8, [0, 5, 9], f, qi, di)
    assert all(np.array_equal(ref[s], got[s]) for s in ref)
    assert t['total'] >= t['kernel'] > 0 and t['loop_other'] > -1e-6


def test_method_specific_stages_and_resident_bytes():
    assert 'S3' not in A1.STAGES['B2'] and 'S3' in A1.STAGES['A']
    assert A1.resident_bytes('A', 10, 214) == 10 * (2 * 214 - 5) * 8
    assert A1.resident_bytes('B2', 10, 214) == 10 * 214 * 3 * 8


def test_owned_subset_has_no_retained_backing():
    big = np.zeros((1000, 50))
    sub = A1.owned_subset(big, np.arange(10))
    assert sub.base is None or sub.base is not big
    assert sub.nbytes == 10 * 50 * 8
    view = big[:10]
    assert view.base is big                              # the original EB-1 pattern kept `big` alive


def test_s4_model_fit_drop_rule_and_prediction():
    import eb1a1_analyse as AA
    rows = [{'comparisons': c, 'call_db_rows': r, 'n_win': 1, 'cpu_s': 1e-6 * c + 1e-7 * r + 0.01}
            for c, r in ((1e5, 1e3), (2e5, 1e3), (4e5, 4e3), (8e5, 2e3))]
    m = AA.s4_model(rows, 'batch_time', 'cpu_s')
    assert m['fit_supported'] and abs(m['coef']['c1'] - 1e-6) < 1e-9
    assert abs(AA.predict(m, {'comparisons': 1e6, 'call_db_rows': 5e3, 'n_win': 1}) - (1 + 5e-4 + 0.01)) < 1e-6
    neg = [{'comparisons': c, 'call_db_rows': r, 'n_win': 1, 'cpu_s': 1e-6 * c - 1e-6 * r + 0.5}
           for c, r in ((1e5, 1e3), (2e5, 1e3), (4e5, 4e3), (8e5, 2e3))]
    assert 'c2' not in AA.s4_model(neg, 'batch_time', 'cpu_s')['coef']           # pre-declared drop of c2 < 0
    mem = [{'P': p, 'maxrss_kB': (5e7 + 1.1 * p) / 1024} for p in (1e7, 3e7, 9e7)]
    mm = AA.s4_model(mem, 'mem', None)
    assert mm['fit_supported'] and abs(mm['coef']['m1'] - 1.1) < 1e-9
