"""SYNTHETIC tests for EB-1 analysis rules and correctness helpers (no measurements, no AdK/mdCATH data).

Run: ../wpd_local_retrieval/build/venv/bin/python -B -m pytest -p no:cacheprovider test_eb1.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import eb1_analyse as AN  # noqa: E402
import eb1_lib as E  # noqa: E402


def test_fit_support_rules():
    good = AN.fit([1, 2, 4], [2.0, 3.0, 5.0])
    assert good['supported'] and abs(good['c1'] - 1.0) < 1e-12
    assert not AN.fit([1, 2], [1.0, 2.0])['supported']                       # < 3 sizes
    assert not AN.fit([1, 2, 4], [5.0, 3.0, 1.0])['supported']               # c1 <= 0
    assert not AN.fit([1, 2, 4], [1.0, 10.0, 4.0])['supported']              # large residual
    assert AN.project(AN.fit([1, 2], [1.0, 2.0]), 10) is None


def _req(r1, r2, r3, r4l, r4e):
    return {'R1': {'value': r1, 'limit': 1.0, 'component': 'engine'}, 'R2': {'value': r2, 'limit': 1.0, 'component': 'engine'},
            'R3': {'value': r3, 'limit': 1.0, 'component': 'loader'}, 'R4a': {'value': r4l, 'limit': 1.0, 'component': 'loader'},
            'R4b': {'value': r4e, 'limit': 1.0, 'component': 'engine'}}


def test_decision_rule_order_and_outcomes():
    assert AN.decide(_req(.5, .5, .5, .5, .5), [], True, [])['outcome'] == 'neither'
    assert AN.decide(_req(.5, .5, 2, .5, .5), [], True, [])['outcome'] == 'loader'
    assert AN.decide(_req(3, .5, .5, .5, .5), [], True, [])['outcome'] == 'engine'
    assert AN.decide(_req(.5, .5, .5, .5, 1.5), [], True, [])['outcome'] == 'engine'       # memory by S3/S4
    both = AN.decide(_req(3, .5, 2, .5, .5), [], True, [])
    assert both['outcome'] == 'both' and both['order'] == ['engine', 'loader']
    assert AN.decide(_req(.5, .5, .5, .5, .5), ['x'], True, [])['outcome'] == 'insufficient_evidence'
    assert AN.decide(_req(.5, .5, .5, .5, .5), [], False, [])['outcome'] == 'insufficient_evidence'
    assert AN.decide(_req(.5, .5, .5, .5, .5), [], True, ['S1'])['outcome'] == 'insufficient_evidence'


def test_topk_and_near_tie_comparison():
    r = np.random.default_rng(0)
    d = r.random(200)
    ok, tie, diff = E.compare_topk(d + 1e-13, d)
    assert ok and not tie and diff < 1e-12
    dt = d.copy()
    srt = np.argsort(dt)
    dt[srt[10]] = dt[srt[9]]                                              # exact boundary tie
    ok, tie, _ = E.compare_topk(dt.copy(), dt)
    assert tie and ok
    wrong = d.copy()
    wrong[srt[0]] = 10.0                                                  # operational result drops a true neighbour
    assert not E.compare_topk(wrong, d)[0]
    D = r.random((3, 50))
    idx = E.topk(D)
    assert all(set(idx[i]) == set(np.argsort(D[i])[:10]) for i in range(3))


def test_comparison_count_model():
    assert AN.P2['n_win'] == AN.P2['L'] - AN.P2['W'] + 1
