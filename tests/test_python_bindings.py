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
