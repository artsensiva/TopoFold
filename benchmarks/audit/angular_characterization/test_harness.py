"""Tests of the WP-A harness itself (SYNTHETIC fixtures; protocol §7).

Corruption tests operate on copies in pytest tmp directories; repository artifacts are only read.
Run with the isolated environment:
  build/venv/bin/python -B -m pytest -p no:cacheprovider test_harness.py
"""
import hashlib
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import generate_fixtures as gen  # noqa: E402
from build_traces import build_traces, recover  # noqa: E402
from common import (FIXTURE_BIN, FIXTURES_JSON, KAPPA, PY_JSON, RUST_JSON, circ_dist,  # noqa: E402
                    direction_bound, load_json, wrap)
from reference import angular_pipeline, bc_from_values, mean_and_resultant  # noqa: E402

MANIFEST = load_json(FIXTURES_JSON)
CASES = {c['id']: c for c in MANIFEST['cases']}


def run_compare(tmp_path, **paths):
    args = [sys.executable, '-B', str(HERE / 'compare.py'), '--out', str(tmp_path / 'summary.json')]
    for k, v in paths.items():
        args += [f"--{k.replace('_', '-')}", str(v)]
    return subprocess.run(args, capture_output=True, text=True, cwd=HERE)


# --- fixture coverage and determinism ---------------------------------------------------------

def test_fixture_matrix_is_the_closed_protocol_matrix():
    fam = {}
    for c in MANIFEST['cases']:
        fam[c['family']] = fam.get(c['family'], 0) + 1
    n_alpha = len(gen.ALPHAS_DEG)
    assert fam == {'F1': 6, 'F2': 3 * n_alpha,
                   'F3': len(gen.F3_EPS) * len(gen.F3_NOISE),
                   'F4': (len(gen.F4_P) * len(gen.F4_EPS) + 1) * len(gen.F4_NOISE) * n_alpha,
                   'F6a': 9, 'F6b': 4}
    eps4 = {c['params']['eps_rad'] for c in MANIFEST['cases'] if c['family'] == 'F4'}
    assert {-1e-6, -1e-3, -1e-1} <= eps4, 'amendment C: negative eps present in F4'
    assert sum(c['r2_p1'] for c in MANIFEST['cases']) == 105


def test_regeneration_reproduces_every_fixture_hash():
    m = gen.Matrix()
    b1 = gen.f1(m)
    gen.f2(m, b1)
    gen.f3(m)
    gen.f4(m)
    gen.f6(m)
    assert len(m.cases) == len(CASES)
    for c in m.cases:
        digest = hashlib.sha256(np.ascontiguousarray(c['tors'], '<f8').tobytes()).hexdigest()
        assert digest == CASES[c['id']]['torsion_sha256'], c['id']


def test_rotated_cases_share_base_draws():
    for c in MANIFEST['cases']:
        if c['family'] not in ('F2', 'F4') or c['alpha_deg'] == 0:
            continue
        base = CASES[f"{c['derived_from']}_a+000"]
        t0 = np.fromfile(FIXTURE_BIN / base['torsion_file'], '<f8')
        ta = np.fromfile(FIXTURE_BIN / c['torsion_file'], '<f8')
        assert circ_dist(ta, t0 + math.radians(c['alpha_deg'])).max() < 1e-14, c['id']


# --- construction recovery ---------------------------------------------------------------------

def test_construction_recovers_prescribed_angles_including_near_pi():
    rng = np.random.default_rng(7)
    tau = rng.uniform(-np.pi, np.pi, (200, 5))
    tau[:, 0] = np.pi - 1e-12
    tau[:, 1] = -np.pi + 1e-12
    X = build_traces(np.full((200, 6), KAPPA), tau)
    k, t = recover(X)
    assert np.abs(k - KAPPA).max() < 1e-12
    assert circ_dist(t, tau).max() < 1e-12
    assert np.allclose(np.linalg.norm(np.diff(X, axis=1), axis=2), 3.8, atol=1e-12)


def test_independent_dihedral_detects_sign_error():
    """A mirrored construction must be detected: recovery is a real check, not a tautology."""
    tau = np.full((1, 5), 0.7)
    X = build_traces(np.full((1, 6), KAPPA), tau)
    _, t = recover(X * np.array([1.0, 1.0, -1.0]))
    assert circ_dist(t, tau).max() > 1.0


def test_binding_recovery_float32_if_available():
    tf = pytest.importorskip('topofold')
    tau = np.random.default_rng(3).uniform(-np.pi, np.pi, (1, 5))
    X = build_traces(np.full((1, 6), KAPPA), tau)[0].astype(np.float32)
    kap, t, _ = tf.compute_invariants(X)
    assert circ_dist(np.asarray(t), tau[0]).max() < 1e-5
    assert np.abs(np.asarray(kap) - KAPPA).max() < 1e-5


# --- REF analytical controls -------------------------------------------------------------------

def test_ref_f1_mean_and_resultant():
    m, R = mean_and_resultant(np.radians([[30.0, 90.0], [90.0, 30.0]]), axis=1)
    assert circ_dist(m, math.radians(60)).max() < 1e-15
    assert np.abs(R - math.cos(math.radians(30))).max() < 1e-15


@pytest.mark.parametrize('eps', [0.1, 1e-3, -1e-3, -0.1])
def test_ref_f3_branch(eps):
    th0 = math.radians(20)
    m, R = mean_and_resultant(wrap(np.array([th0, th0 + math.pi - eps])))
    expected = th0 + (math.pi - eps) / 2 + (math.pi if eps < 0 else 0.0)
    assert abs(R - abs(math.sin(eps / 2))) < 1e-15
    assert circ_dist(m, expected) < 1e-12


def test_ref_f4_resultant():
    x = np.concatenate([np.zeros(510), np.full(490, math.pi)])
    _, R = mean_and_resultant(x)
    assert abs(R - 0.02) < 1e-13


def test_ref_bc_guard_and_two_point():
    assert bc_from_values(np.zeros(100))['guard'] == 'variance<=1e-14'
    two = bc_from_values(np.r_[np.zeros(500), np.ones(500)])
    assert 0.99 < two['bc'] <= 1.0
    r = angular_pipeline(np.radians([[30.0, 90.0]] * 10))
    assert r['bc']['guard'] == 'variance<=1e-14'


def test_direction_bound_classifies_cancellation():
    assert direction_bound(2, 1.0) < 1e-15
    assert direction_bound(2, 1e-16) > 0.1
    assert math.isinf(float(direction_bound(2, 0.0)))


# --- executed results: evidence and corruption detection --------------------------------------

def test_executed_r2_replication_and_recovery():
    rust = load_json(RUST_JSON)['cases']
    r2 = [v['r2'] for v in rust.values() if 'r2' in v]
    assert len(r2) == 105
    assert all(v['replication']['bc_tau_bits_equal'] and v['replication']['bc_kappa_bits_equal']
               for v in r2)
    assert max(v['recovery']['max_tau_circ_err'] for v in r2) <= 1e-9


def test_compare_passes_on_executed_results(tmp_path):
    assert run_compare(tmp_path).returncode == 0


def test_compare_fails_on_broken_replication_flag(tmp_path):
    doc = load_json(RUST_JSON)
    cid = next(k for k, v in doc['cases'].items() if 'r2' in v)
    doc['cases'][cid]['r2']['replication']['bc_tau_bits_equal'] = False
    bad = tmp_path / 'rust.json'
    bad.write_text(json.dumps(doc))
    assert run_compare(tmp_path, rust=bad).returncode == 1


def test_compare_fails_on_altered_conditioned_bc(tmp_path):
    doc = load_json(PY_JSON)
    doc['cases']['F2b_a+000']['ref']['moments']['bc'] += 1e-6
    bad = tmp_path / 'py.json'
    bad.write_text(json.dumps(doc))
    proc = run_compare(tmp_path, python=bad)
    assert proc.returncode == 1 and 'R1/R2 vs REF' in proc.stdout


def test_compare_fails_on_missing_case(tmp_path):
    doc = load_json(RUST_JSON)
    doc['cases'].pop('F1a_s0')
    bad = tmp_path / 'rust.json'
    bad.write_text(json.dumps(doc))
    assert run_compare(tmp_path, rust=bad).returncode == 1


def test_compare_fails_on_corrupted_fixture(tmp_path):
    d = tmp_path / 'fixtures'
    shutil.copytree(FIXTURE_BIN, d)
    f = d / CASES['F3_e+1e-01_n0']['torsion_file']
    b = bytearray(f.read_bytes())
    b[0] ^= 1
    f.write_bytes(bytes(b))
    proc = run_compare(tmp_path, fixture_dir=d)
    assert proc.returncode == 1 and 'fixture files match manifest hashes' in proc.stdout
