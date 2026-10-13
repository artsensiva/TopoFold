"""SYNTHETIC regression tests for BQ-2 (tiled stress rows, preparation repetitions, memory rule, model form, V2b
unchanged). Run: ../wpd_local_retrieval/build/venv/bin/python -B -m pytest -p no:cacheprovider test_bq2.py
"""
import gc
import hashlib
import json
import sys
import tracemalloc
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'eb1'))
import bq2_lib as Q  # noqa: E402

RNG = np.random.default_rng(11)


def fake_src(F=37):
    v = RNG.normal(size=(F, Q.SEG, 3))
    v *= 3.8 / np.linalg.norm(v, axis=-1, keepdims=True)
    return np.cumsum(v, axis=1)


def test_v2b_module_is_the_bq1_frozen_file():
    rec = json.loads((HERE.parent / 'b2opt' / 'results' / 'b2opt_pre_measurement_record.json').read_text())
    got = hashlib.sha256(Path(Q.B.__file__).read_bytes()).hexdigest()
    assert got == rec['sha256']['benchmarks/audit/b2opt/b2opt_lib.py']


def test_tiled_rows_shape_determinism_junctions_and_real_windows():
    src = fake_src()
    X = Q.tiled_rows(src, [0, 1, 40])
    assert X.shape == (3, 1000, 3) and X.flags['C_CONTIGUOUS'] and X.dtype == np.float64
    assert np.array_equal(X, Q.tiled_rows(src, [0, 1, 40]))
    assert np.allclose(np.linalg.norm(X[0, Q.SEG] - X[0, Q.SEG - 1]), Q.JUNCTION_A)
    # a window inside a segment is a rigid translation of the real frame window: RMSD 0 to it
    seg = src[(5 * 40 + 2) % len(src)]
    a, b = X[2, 2 * Q.SEG + 10:2 * Q.SEG + 18], seg[10:18]
    assert Q.B.ref_svd(a, b[None])[0] <= 1e-6
    assert np.array_equal(Q.tiled_rows(src, [0])[0], Q.tiled_rows(src, [len(src)])[0])   # repeats after F rows


def test_prepare_reps_releases_previous_prepared_arrays():
    X = fake_src(301)[:, :60]
    di = np.arange(1, 301)
    wins = list(range(40))
    one = len(di) * len(wins) * 25 * 8
    tracemalloc.start()
    prep, walls, _ = Q.prepare_reps(X, di, wins, 5)
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    assert len(walls) == 5 and prep[0].shape == (40, 300, 8, 3)
    assert peak < 1.3 * one + 2 * 1024 ** 2                    # only one prepared array alive at a time
    del prep
    gc.collect()


def test_memory_rule_and_repetition_rule():
    pred = Q.predicted_peak_bytes(4003, 1000, 4000, 993)
    assert pred == Q.BASE_RSS + Q.TRANSIENT + 4003 * 1000 * 24 + 4000 * 993 * 25 * 8
    assert Q.memory_ok(pred, (pred + Q.HEADROOM) // 1024 + 1)
    assert not Q.memory_ok(pred, (pred + Q.HEADROOM) // 1024 - 1)
    assert [Q.prep_reps_for(x) for x in (0.01, 0.5, 2.0)] == [9, 5, 3]


def test_prep_model_form_recovers_per_window_term():
    import eb1a1_analyse as AA
    rows = []
    for nd in (100, 700, 2693):
        for nw in (5, 20, 60, 207):                             # SYNTHETIC: 4e-7 s per row-window + 4e-5 s per window
            rows.append({'comparisons': nd * nw, 'n_win': nw, 'y': 4e-7 * nd * nw + 4e-5 * nw + 1e-4})
    m = AA.s4_model(rows, 'single_time', 'y')
    assert abs(m['coef']['c1'] - 4e-7) < 1e-12 and abs(m['coef']['c2'] - 4e-5) < 1e-10
    assert m['max_fit_rel_residual'] < 1e-6
