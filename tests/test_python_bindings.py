"""
Integration and regression tests for TopoFold Python bindings.
"""

import os
import numpy as np
import pytest
import topofold as tf


def test_read_pdb_and_compute_invariants():
    pdb_path = os.path.join(os.path.dirname(__file__), "..", "benchmarks", "data", "5PTI.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip("5PTI.pdb not present")

    coords = tf.read_pdb(pdb_path, chain="A")
    assert coords.shape == (58, 3)
    assert coords.dtype == np.float32

    kappa, tau, writhe = tf.compute_invariants(coords)
    assert len(kappa) == 56
    assert len(tau) == 55
    assert len(writhe) > 0

    # Curvature values in [0, pi]
    assert np.all(kappa >= 0.0)
    assert np.all(kappa <= np.pi)

    # Torsion values in (-pi, pi]
    assert np.all(tau >= -np.pi)
    assert np.all(tau <= np.pi)


def test_dcd_indexing_and_subcurve_query():
    dcd_path = os.path.join(os.path.dirname(__file__), "..", "benchmarks", "data", "bpti_equilibrium.dcd")
    pdb_path = os.path.join(os.path.dirname(__file__), "..", "benchmarks", "data", "5PTI.pdb")
    if not os.path.exists(dcd_path) or not os.path.exists(pdb_path):
        pytest.skip("BPTI trajectory files not present")

    coords = tf.read_pdb(pdb_path, chain="A")
    index = tf.ConformationalIndex.from_dcd(
        dcd_path,
        ca_indices=list(range(58)),
        window_size=8,
        coarse_radius=0.25,
    )
    assert len(index) == 2500

    # Subcurve query on active site loop (residues 10..18, 0-indexed 9..17)
    hits = index.query_subcurve(9, 17, coords, k=10)
    assert len(hits) == 10
    # Nearest neighbor should have very low Frechet distance
    assert hits[0][1] >= 0.0
    assert hits[0][1] < 1.0


def test_read_dcd():
    dcd_path = os.path.join(os.path.dirname(__file__), "..", "benchmarks", "data", "bpti_equilibrium.dcd")
    if not os.path.exists(dcd_path):
        pytest.skip("bpti_equilibrium.dcd not present")

    traj = tf.read_dcd(dcd_path)
    assert traj.shape == (2500, 58, 3)
    assert traj.dtype == np.float32


def test_scan_cryptic_pockets():
    dcd_path = os.path.join(os.path.dirname(__file__), "..", "benchmarks", "data", "bpti_equilibrium.dcd")
    if not os.path.exists(dcd_path):
        pytest.skip("bpti_equilibrium.dcd not present")

    traj = tf.read_dcd(dcd_path)
    candidates = tf.scan_cryptic_pockets(traj, window_size=8, bc_threshold=0.6)
    assert len(candidates) > 0

    # Top candidate must have high bimodality (> 0.8) and overlap with known loop (residues 10..18, 0-indexed 9..17)
    cand_start, cand_end, score = candidates[0]
    assert score >= 0.80
    overlaps = not (cand_end < 9 or cand_start > 17)
    assert overlaps

    scores, bc_tau, bc_kappa = tf.compute_bimodality_profile(traj, window_size=8)
    assert len(scores) == 58 - 8 + 1
    assert len(bc_tau) == len(scores)
    assert len(bc_kappa) == len(scores)
    assert np.all(scores >= 0.0)
    assert np.all(scores <= 1.0)


def test_read_pdb_cbeta_and_theta_beta():
    pdb_path = os.path.join(os.path.dirname(__file__), "..", "benchmarks", "data", "5PTI.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip("5PTI.pdb not present")

    ca_coords, cb_coords = tf.read_pdb(pdb_path, chain="A", extract_cbeta=True)
    assert ca_coords.shape == (58, 3)
    assert cb_coords.shape == (58, 3)
    assert ca_coords.dtype == np.float32
    assert cb_coords.dtype == np.float32

    # Compute theta_beta directly
    theta_beta = tf.compute_theta_beta(ca_coords, cb_coords)
    assert theta_beta.shape == (56,)
    assert theta_beta.dtype == np.float64
    assert np.all(theta_beta > -np.pi)
    assert np.all(theta_beta <= np.pi)

    # Compute invariants with cb_coords -> returns 4-tuple
    kappa, tau, writhe, theta_inv = tf.compute_invariants(ca_coords, cb_coords)
    assert len(kappa) == 56
    assert len(tau) == 55
    assert len(theta_inv) == 56
    assert np.allclose(theta_beta, theta_inv, atol=1e-6)


def test_se3_strict_invariance_cbeta():
    """Verify that theta_beta is strictly SE(3)-invariant under arbitrary 3D rotation and translation."""
    pdb_path = os.path.join(os.path.dirname(__file__), "..", "benchmarks", "data", "5PTI.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip("5PTI.pdb not present")

    ca_coords, cb_coords = tf.read_pdb(pdb_path, chain="A", extract_cbeta=True)
    theta_ref = tf.compute_theta_beta(ca_coords, cb_coords)

    # Random rotation matrix (Euler angles)
    rng = np.random.default_rng(42)
    angles = rng.uniform(0, 2 * np.pi, size=3)
    c, s = np.cos(angles), np.sin(angles)
    Rx = np.array([[1, 0, 0], [0, c[0], -s[0]], [0, s[0], c[0]]], dtype=np.float32)
    Ry = np.array([[c[1], 0, s[1]], [0, 1, 0], [-s[1], 0, c[1]]], dtype=np.float32)
    Rz = np.array([[c[2], -s[2], 0], [s[2], c[2], 0], [0, 0, 1]], dtype=np.float32)
    R = (Rz @ Ry @ Rx).astype(np.float32)
    t = rng.uniform(-50, 50, size=(1, 3)).astype(np.float32)

    ca_rot = (ca_coords @ R.T) + t
    cb_rot = (cb_coords @ R.T) + t

    theta_rot = tf.compute_theta_beta(ca_rot, cb_rot)
    diff = np.abs(theta_ref - theta_rot)
    # Account for 2*pi periodicity at boundary
    diff = np.minimum(diff, 2 * np.pi - diff)
    max_err = np.max(diff)
    assert max_err < 1e-4, f"SE(3) invariance violated: max error {max_err}"


def test_rotamer_gating_detection():
    """Verify that C-beta ribbon orientation detects rotamer flips where C-alpha is completely rigid."""
    # 20 residues ideal extended/helical backbone
    n_res = 20
    n_frames = 400
    ca_base = np.zeros((n_res, 3), dtype=np.float32)
    for i in range(n_res):
        ca_base[i] = [i * 3.8, np.sin(i * 0.5) * 1.5, np.cos(i * 0.5) * 1.5]

    # Rigid CA coordinates across all frames (thermal noise only ~0.001 A)
    rng = np.random.default_rng(123)
    ca_traj = np.repeat(ca_base[np.newaxis, :, :], n_frames, axis=0)
    ca_traj += rng.normal(0, 0.001, size=ca_traj.shape).astype(np.float32)

    # Base CB coordinates
    cb_traj = ca_traj.copy()
    cb_traj[:, :, 2] += 1.52

    # Residue 10 undergoes a bistable rotamer flip:
    # State A (frames 0..199): CB points along +Y
    # State B (frames 200..399): CB points along -Y
    cb_traj[:200, 10, 1] = ca_traj[:200, 10, 1] + 1.52
    cb_traj[200:, 10, 1] = ca_traj[200:, 10, 1] - 1.52

    # 1. Without C-beta: pure CA bimodality profile should be unimodal everywhere (BC < 0.555)
    ca_scores, _, _ = tf.compute_bimodality_profile(ca_traj, window_size=8)
    assert np.all(ca_scores < 0.555), f"Expected unimodal CA scores, got max {np.max(ca_scores)}"

    # 2. With C-beta: bimodality profile detects the rotamer flip at residue 10 with BC > 0.80
    scores, bc_tau, bc_kappa, bc_theta = tf.compute_bimodality_profile(ca_traj, window_size=8, cb_coords=cb_traj)
    # The windows covering residue 10 (windows 3..10) must show high bc_theta
    assert np.max(bc_theta) > 0.80, f"Expected high bc_theta, got {np.max(bc_theta)}"
    assert np.max(scores) > 0.80

    candidates = tf.scan_cryptic_pockets(ca_traj, window_size=8, bc_threshold=0.6, cb_coords=cb_traj)
    assert len(candidates) > 0
    c_start, c_end, score = candidates[0]
    assert score > 0.80
    assert c_start <= 10 <= c_end, f"Candidate [{c_start}, {c_end}] must cover rotamer gate at residue 10"

