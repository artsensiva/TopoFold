"""WP-C harness tests: feature contracts, tie handling, input integrity, corruption detection.

Each corruption violates a stated contract (protocol §7) and must be detected. Repository
artifacts are only read; corruption tests use substituted functions or tmp copies.
Run with the WP-A isolated interpreter:
  ../angular_characterization/build/venv/bin/python -B -m pytest -p no:cacheprovider test_ordered.py
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from features import (bond_labels, d_chord, feature_contracts, ordered_features,  # noqa: E402
                      window_torsions)
from inputs import InputHashError, WPA_FIXTURE_DIR, load_case, manifest, read_verified  # noqa: E402
from retrieval import knn_purity  # noqa: E402

RESULTS = HERE / 'results' / 'SYNTHETIC_ordered_representation_results.json'


def failed(checks):
    return {name for name, ok, _ in checks if not ok}


# --- feature contracts (§5a) and the corruption set (§7, items 1–4) ------------------------------

def test_feature_contracts_pass():
    assert not failed(feature_contracts())


def test_distance_checks_cannot_see_consistent_permutation():
    """Why amendment A exists: a consistent permutation leaves every distance unchanged."""
    rng = np.random.default_rng(1)
    t = rng.uniform(-np.pi, np.pi, (50, 5))
    # equal up to summation-order round-off (not bitwise: the sum order differs)
    assert np.abs(d_chord(t, t) - d_chord(t[:, ::-1], t[:, ::-1])).max() <= 1e-15


def test_corruption_permuted_positions_detected():
    bad = lambda t: ordered_features(np.asarray(t)[..., ::-1])
    assert any(n.startswith('FC1') for n in failed(feature_contracts(feature_fn=bad)))


def test_corruption_cos_sin_swap_detected():
    def bad(t):
        x = ordered_features(t)
        y = x.copy()
        y[..., 0::2], y[..., 1::2] = x[..., 1::2], x[..., 0::2]
        return y
    assert any(n.startswith('FC1') for n in failed(feature_contracts(feature_fn=bad)))


def test_corruption_sine_sign_detected():
    def bad(t):
        x = ordered_features(t)
        x[..., 1::2] *= -1
        return x
    f = failed(feature_contracts(feature_fn=bad))
    assert any(n.startswith('FC1') for n in f) and any(n.startswith('FC4') for n in f)


def test_corruption_window_offset_detected():
    f = failed(feature_contracts(window_fn=lambda c, i, W: window_torsions(c, i + 1, W)))
    assert any(n.startswith('FC2') for n in f)
    f = failed(feature_contracts(bond_fn=lambda i, W: bond_labels(i + 1, W)))
    assert any(n.startswith('FC3') for n in f)


# --- tie handling (§6, corruption item 5) ---------------------------------------------------------

def test_partial_boundary_tie_expected_purity():
    """Query 0: neighbour 1 (same label) below the boundary; boundary set {2, 3, 4} at d_k = 2
    holds one same-label frame; k = 2 -> purity = (1 + 1·(1/3)) / 2 = 2/3."""
    D = np.full((5, 5), 9.0)
    D[0, 1:] = [1.0, 2.0, 2.0, 2.0]
    D[1:, 0] = D[0, 1:]
    labels = np.array([0, 0, 1, 1, 0])
    res = knn_purity(D, labels, 2, per_query=True)
    assert abs(res['per_query_purity'][0] - 2 / 3) < 1e-15


def test_fully_tied_set_gives_chance_not_index_order():
    """Contract: ties are label- and order-independent. A fully tied set must score chance under
    any frame ordering; the corruption (breaking ties by frame index) depends on the ordering."""
    def naive(labels, k=2):
        n = len(labels)
        return np.mean([np.mean(labels[[j for j in range(n) if j != q][:k]] == labels[q])
                        for q in range(n)])
    blocks = np.array([0, 0, 0, 1, 1, 1])
    alternating = np.array([0, 1, 0, 1, 0, 1])          # same label multiset, different order
    for labels in (blocks, alternating):
        ours = knn_purity(np.zeros((6, 6)), labels, 2)
        assert abs(ours['purity'] - ours['chance']) < 1e-15 and ours['tie_frequency'] == 1.0
    assert abs(naive(blocks) - naive(alternating)) > 0.05


def test_purity_invariant_to_frame_order():
    rng = np.random.default_rng(4)
    X = rng.normal(size=(60, 3))
    X[:30] += 1.0
    labels = np.r_[np.zeros(30), np.ones(30)]
    D = np.linalg.norm(X[:, None] - X[None], axis=-1)
    perm = rng.permutation(60)
    a = knn_purity(D, labels, 10)['purity']
    b = knn_purity(D[np.ix_(perm, perm)], labels[perm], 10)['purity']
    assert abs(a - b) < 1e-12


# --- input integrity (corruption item 6) and regeneration ----------------------------------------

def test_altered_fixture_byte_detected(tmp_path):
    c = manifest()['F1a_s0']
    src = WPA_FIXTURE_DIR / c['torsion_file']
    dst = tmp_path / c['torsion_file']
    shutil.copy2(src, dst)
    b = bytearray(dst.read_bytes())
    b[7] ^= 1
    dst.write_bytes(bytes(b))
    with pytest.raises(InputHashError):
        read_verified(dst, c['torsion_sha256'], (c['F'], c['ncols']))


@pytest.mark.parametrize('cid', ['F1a_s0', 'F3_e-1e-03_n1deg', 'F4_p0.49_e+1e-03_n1deg_a+000'])
def test_in_memory_regeneration_matches_frozen_manifest(tmp_path, cid):
    t_file, x_file, _ = load_case(cid)
    t_mem, x_mem, src = load_case(cid, fixture_dir=tmp_path)    # empty dir -> regenerate in memory
    assert src == 'in_memory_regeneration'
    assert np.array_equal(t_file, t_mem) and np.array_equal(x_file, x_mem)


# --- stored results: reproducibility and coverage (corruption items 7–8) --------------------------

def run_verify(path):
    return subprocess.run([sys.executable, '-B', str(HERE / 'run.py'), '--verify', '--results',
                           str(path)], capture_output=True, text=True, cwd=HERE)


def test_verify_passes_on_stored_results():
    p = run_verify(RESULTS)
    assert p.returncode == 0, p.stdout[-1500:]


def test_verify_detects_altered_value_and_missing_case(tmp_path):
    doc = json.loads(RESULTS.read_text())
    doc['retrieval']['F1a_s0']['T0']['D_ch']['k10']['purity'] -= 0.01
    bad = tmp_path / 'altered.json'
    bad.write_text(json.dumps(doc))
    assert run_verify(bad).returncode == 1
    doc = json.loads(RESULTS.read_text())
    doc['cases'].pop('F2a_a+045')
    bad2 = tmp_path / 'missing.json'
    bad2.write_text(json.dumps(doc))
    assert run_verify(bad2).returncode == 1
