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


def test_idp_spectral_topological_density():
    # Synthetic 50-residue polymer ensemble across 100 frames
    n_frames = 100
    n_res = 50
    rng = np.random.default_rng(42)

    traj = np.zeros((n_frames, n_res, 3), dtype=np.float32)
    for f in range(n_frames):
        # Straight chain
        for i in range(n_res):
            traj[f, i] = [i * 3.8, 0.0, 0.0]

    # For collinear chains, local writhe and curvature are identically zero
    s_topo = tf.compute_idp_topological_density(traj, window_radius=4)
    assert s_topo.shape == (n_res,)
    assert np.all(s_topo < 1e-6)

    mean_d, var_d, std_d, z_d = tf.compute_idp_density_profile(traj, window_radius=4)
    assert len(mean_d) == n_res
    assert len(var_d) == n_res
    assert len(std_d) == n_res
    assert len(z_d) == n_res


def test_idp_transient_motif_detection():
    # Build a 60-residue chain with a chiral helical turn in the center (residues 25..35)
    n_frames = 100
    n_res = 60
    traj = np.zeros((n_frames, n_res, 3), dtype=np.float32)

    for f in range(n_frames):
        # Base straight chain
        for i in range(n_res):
            traj[f, i] = [i * 3.8, 0.0, 0.0]

        # 30% of frames have a chiral helical loop at residues 25..35
        if f % 3 == 0:
            for i in range(25, 36):
                t = (i - 25)
                angle = t * 1.5
                traj[f, i] = [25 * 3.8 + t * 1.5, 3.0 * np.cos(angle), 3.0 * np.sin(angle)]

    s_topo = tf.compute_idp_topological_density(traj, window_radius=4)
    assert s_topo.shape == (n_res,)
    # Tail residues must have zero/low compactness
    assert s_topo[5] < 1e-4
    assert s_topo[55] < 1e-4
    # Central loop must have elevated compactness
    assert np.max(s_topo[25:36]) > 0.01

    motifs = tf.detect_transient_motifs(traj, window_size=8, z_threshold=1.0)
    assert len(motifs) > 0
    top = motifs[0]
    # top is (start, end, peak, mean, peak_val, z_score)
    assert top[0] <= 32 and top[1] >= 28, f"Top motif {top} should cover residues 25..35"
    assert top[5] > 1.5, f"Expected elevated Z-score, got {top[5]}"


def test_intermolecular_allosteric_network_and_cooperativity():
    # 2 chains: Chain A (12 residues) and Chain B (16 residues) across 20 frames
    n_frames = 20
    n_a = 12
    n_b = 16

    traj_a = np.zeros((n_frames, n_a, 3), dtype=np.float32)
    traj_b = np.zeros((n_frames, n_b, 3), dtype=np.float32)

    for f in range(n_frames):
        # Coupled conformational breathing in curvature and torsion
        rad_a = 2.0 + f * 0.15
        pitch_a = 3.8 + f * 0.05
        for i in range(n_a):
            ang = i * 0.6
            traj_a[f, i] = [i * pitch_a, rad_a * np.cos(ang), rad_a * np.sin(ang)]

        rad_b = 2.5 + f * 0.20
        pitch_b = 3.8 + f * 0.05
        for j in range(n_b):
            ang = j * 0.4
            traj_b[f, j] = [j * pitch_b, rad_b * np.cos(ang), rad_b * np.sin(ang)]

    inter_mat = tf.compute_intermolecular_allosteric_network(traj_a, traj_b)
    assert inter_mat.shape == (n_a, n_b)
    assert np.all(inter_mat >= 0.0)
    assert np.all(inter_mat <= 1.0)

    # All interface pairs between coupled breathing chains should have high cooperativity
    interface_pairs = [(2, 3), (3, 4), (4, 5)]
    coop = tf.compute_ternary_cooperativity_index(inter_mat, interface_pairs)
    assert coop > 0.50
    assert coop <= 1.0




