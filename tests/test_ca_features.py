"""Tests for topofold_inputs.features (κ/τ from the TopoFold core with explicit validity). SYNTHETIC fixtures only.

Runs where the ``topofold`` extension is importable (e.g. .venv). The end-to-end pipeline test additionally needs an
interpreter with MDAnalysis, taken from $TOPOFOLD_ADAPTER_PYTHON or the WP-D audit venv if present.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

tf = pytest.importorskip("topofold")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from topofold_inputs.features import (  # noqa: E402
    COLLINEAR, CORE_ERROR, DEGENERATE, MISSING, NONFINITE, SHORT, TERMINAL, VALID, compute_ca_features,
)


# --- independent NumPy reference (textbook formulas, not the core) -------------------------------------------------

def ref_kappa(p0, p1, p2):
    a, b = p1 - p0, p2 - p1
    return float(np.arccos(np.clip(a @ b / (np.linalg.norm(a) * np.linalg.norm(b)), -1.0, 1.0)))


def ref_dihedral(p0, p1, p2, p3):
    """IUPAC dihedral: atan2(|b2| b1·(b2×b3), (b1×b2)·(b2×b3))."""
    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    return float(np.arctan2(np.linalg.norm(b2) * (b1 @ n2), n1 @ n2))


def circ_diff(a, b):
    return np.abs(np.angle(np.exp(1j * (np.asarray(a) - np.asarray(b)))))


def helix(n, noise=0.0, seed=0):
    """SYNTHETIC Cα-like helix (radius 2.3 Å, 100° per residue, 1.5 Å rise) with optional deterministic noise."""
    t = np.arange(n) * np.deg2rad(100.0)
    X = np.stack([2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)], axis=1)
    if noise:
        X = X + np.random.default_rng(seed).normal(scale=noise, size=X.shape)
    return X


def run(X, segments=None, present=None, **kw):
    X = np.asarray(X, dtype=np.float64)
    if X.ndim == 2:
        X = X[None]
    n = X.shape[1]
    present = np.ones(n, bool) if present is None else np.asarray(present)
    segments = [[0, n]] if segments is None else segments
    return compute_ca_features(X, present, segments, **kw)


def check_against_reference(feat, X, f=0, tol=1e-10):
    P = np.asarray(X, dtype=np.float32).astype(np.float64)  # the core receives float32
    for i in np.flatnonzero(feat.kappa_valid[f]):
        assert abs(feat.kappa[f, i] - ref_kappa(P[i - 1], P[i], P[i + 1])) < tol
    for h in np.flatnonzero(feat.tau_valid[f]):
        assert circ_diff(feat.tau[f, h], ref_dihedral(P[h - 1], P[h], P[h + 1], P[h + 2])) < tol


def invariants(feat):
    for val, reason in ((feat.kappa, feat.kappa_reason), (feat.tau, feat.tau_reason)):
        assert np.isfinite(val[reason == VALID]).all()
        assert np.isnan(val[reason != VALID]).all()


# --- tests ---------------------------------------------------------------------------------------------------------

def test_matches_independent_reference_and_attribution():
    X = helix(20, noise=0.3, seed=1)
    feat = run(X)
    invariants(feat)
    assert feat.kappa_valid[0, 1:-1].all() and not feat.kappa_valid[0, [0, -1]].any()
    assert feat.tau_valid[0, 1:-1].all() and not feat.tau_valid[0, [0, -1]].any()
    check_against_reference(feat, X)
    assert feat.kappa_atoms[5].tolist() == [4, 5, 6] and feat.tau_atoms[5].tolist() == [4, 5, 6, 7]
    assert feat.kappa_atoms[0].tolist() == [-1, 0, 1] and feat.tau_atoms[-1].tolist() == [17, 18, 19, -1]


def test_kink_is_attributed_to_its_residue():
    X = np.array([[3.8 * i, 0.3 * (i % 2), 0.0] for i in range(10)])
    X[6:, 1] += np.arange(1, 5) * 3.0              # direction changes at residue 6
    feat = run(X)
    k = np.where(feat.kappa_valid[0], feat.kappa[0], -1)
    assert int(np.argmax(k)) == 5 or int(np.argmax(k)) == 6
    check_against_reference(feat, X)


def test_rigid_motion_invariance():
    X = helix(16, noise=0.2, seed=2)
    rng = np.random.default_rng(3)
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    R = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                  [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                  [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
    Y = X @ R.T + np.array([40.0, -25.0, 13.0])
    a, b = run(X), run(Y)
    assert (a.kappa_reason == b.kappa_reason).all() and (a.tau_reason == b.tau_reason).all()
    # bound set by float32 rounding of the core input, not by the geometry
    assert np.nanmax(np.abs(a.kappa - b.kappa)) < 1e-4
    assert np.nanmax(circ_diff(a.tau[a.tau_valid], b.tau[b.tau_valid])) < 1e-4


def test_angular_wrapping_near_pi():
    n = 12
    X = np.array([[3.8 * i, 1.0 * (-1) ** i, 0.02 * ((-1) ** (i // 2))] for i in range(n)])  # near-trans
    feat = run(X)
    t = feat.tau[0, feat.tau_valid[0]]
    assert (t > 3.0).any() and (t < -3.0).any()     # both sides of the branch cut occur
    check_against_reference(feat, X)                # compared with circular differences


def test_breaks_are_never_crossed():
    X = helix(10, noise=0.2, seed=4)
    whole = run(X)
    split = run(X, segments=[[0, 5], [5, 10]])
    assert split.kappa_reason[0, 4] == TERMINAL and split.kappa_reason[0, 5] == TERMINAL
    assert (split.tau_reason[0, [3, 4, 5]] == TERMINAL).all()
    for i in (1, 2, 3, 6, 7, 8):
        assert split.kappa[0, i] == whole.kappa[0, i]
    # second segment computed on its own atoms only
    alone = run(X[5:])
    np.testing.assert_array_equal(split.kappa[0, 6:9], alone.kappa[0, 1:4])
    invariants(split)


def test_missing_residue_propagates():
    X = helix(10, noise=0.2, seed=5)
    X[5] = np.nan
    present = np.ones(10, bool)
    present[5] = False
    feat = run(X, segments=[[0, 5], [6, 10]], present=present)
    assert feat.kappa_reason[0, 5] == MISSING
    assert (feat.tau_reason[0, [3, 4, 5, 6]] == MISSING).all()
    assert feat.kappa_reason[0, 4] == TERMINAL and feat.kappa_reason[0, 6] == TERMINAL
    invariants(feat)
    check_against_reference(feat, X)


def test_degenerate_segment_and_collinear_vertex():
    X = helix(14, noise=0.2, seed=6)
    X[7] = X[6]                                    # zero-length segment 6 -> 7
    feat = run(X)
    assert feat.kappa_reason[0, 6] == DEGENERATE and feat.kappa_reason[0, 7] == DEGENERATE
    assert (feat.tau_reason[0, [5, 6, 7]] == DEGENERATE).all()
    assert feat.kappa_valid[0, 2] and feat.kappa_valid[0, 10]
    check_against_reference(feat, X)
    invariants(feat)

    Y = helix(14, noise=0.2, seed=7)
    Y[8] = Y[7] + (Y[7] - Y[6])                    # residue 7 is an exactly collinear vertex
    feat = run(Y)
    assert feat.kappa_reason[0, 7] == COLLINEAR
    assert (feat.tau_reason[0, [6, 7]] == COLLINEAR).all()
    assert feat.tau_reason[0, 5] == VALID and feat.tau_reason[0, 8] == VALID   # neighbours of the collinear vertex
    assert feat.kappa_reason[0, 6] == VALID and feat.kappa_reason[0, 8] == VALID
    check_against_reference(feat, Y)
    assert feat.diagnostics["n_core_errors"] == 0


def test_short_segment_and_input_validation():
    X = helix(3)
    feat = run(X)
    assert feat.kappa_reason[0, 1] == SHORT            # core needs at least 4 atoms
    with pytest.raises(ValueError):
        run(helix(6), segments=[[0, 6]], present=[True] * 5 + [False])
    bad = helix(6)
    bad[2, 0] = np.inf
    with pytest.raises(ValueError):
        run(bad)
    with pytest.raises(ValueError):
        run(helix(6), segments=[[0, 4], [3, 6]])


class _RaisingCore:
    @staticmethod
    def compute_invariants(_):
        raise ValueError("SYNTHETIC core failure")


class _NaNCore:
    @staticmethod
    def compute_invariants(c):
        n = len(c)
        return np.full(n - 2, np.nan), np.full(n - 3, np.nan), np.zeros(1)


def test_core_failures_are_recorded_not_substituted():
    X = helix(8, noise=0.1, seed=8)
    feat = run(X, core=_RaisingCore)
    assert (feat.kappa_reason[0, 1:-1] == CORE_ERROR).all() and feat.diagnostics["n_core_errors"] == 1
    assert np.isnan(feat.kappa).all() and np.isnan(feat.tau).all()
    feat = run(X, core=_NaNCore)
    assert (feat.kappa_reason[0, 1:-1] == NONFINITE).all()
    invariants(feat)


def test_multiple_frames_independent():
    X = np.stack([helix(10, noise=0.2, seed=s) for s in range(3)])
    feat = run(X)
    for f in range(3):
        check_against_reference(feat, X[f], f=f)


def test_core_provenance_distinguishes_wrapper_and_native_extension():
    import importlib.machinery
    import hashlib
    from topofold_inputs.features import core_provenance
    prov = core_provenance(tf)
    w, n = prov["wrapper"], prov["native_extension"]
    assert n["status"] == "resolved" and n["module"] == tf.compute_invariants.__module__
    assert any(n["path"].endswith(sfx) for sfx in importlib.machinery.EXTENSION_SUFFIXES)
    assert Path(n["path"]).resolve() != Path(w["path"]).resolve()
    assert n["sha256"] == hashlib.sha256(Path(n["path"]).read_bytes()).hexdigest()
    assert w["sha256"] == hashlib.sha256(Path(w["path"]).read_bytes()).hexdigest()
    assert n["sha256"] != w["sha256"]


def test_core_provenance_unresolved_extension_is_explicit():
    import types
    from topofold_inputs.features import core_provenance
    fake = types.ModuleType("SYNTHETIC_core")
    fake.__file__ = __file__

    def compute_invariants(_):  # a pure-Python function: defining module is not an extension
        raise NotImplementedError
    fake.compute_invariants = compute_invariants
    n = core_provenance(fake)["native_extension"]
    assert n["status"] == "unresolved" and n["sha256"] is None and n["path"] is None
    compute_invariants.__module__ = "SYNTHETIC_not_loaded"
    n = core_provenance(fake)["native_extension"]
    assert n["status"] == "unresolved" and "not loaded" in n["reason"]
    del fake.compute_invariants
    assert core_provenance(fake)["native_extension"]["status"] == "unresolved"


def _adapter_python():
    cand = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(
        ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
    if not Path(cand).is_file():
        return None
    ok = subprocess.run([cand, "-c", "import MDAnalysis"], capture_output=True).returncode == 0
    return cand if ok else None


def test_end_to_end_pipeline(tmp_path):
    py = _adapter_python()
    if py is None:
        pytest.skip("no interpreter with MDAnalysis available for the adapter stage")
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"))
    subprocess.run([py, str(ROOT / "python/examples/ca_adapter_example.py"), "--out", str(tmp_path)],
                   env=env, check=True, capture_output=True)
    rc = subprocess.run([sys.executable, "-m", "topofold_inputs.pipeline",
                         str(tmp_path / "SYNTHETIC_ca_adapter_fixture.pdb"), "--run-id", "SYNTHETIC_e2e",
                         "--out", str(tmp_path), "--adapter-python", py,
                         "--correspondence", str(tmp_path / "SYNTHETIC_correspondence.json"),
                         "--altloc", "highest_occupancy"], env=env, capture_output=True, text=True)
    assert rc.returncode == 0, rc.stderr
    doc = json.loads((tmp_path / "SYNTHETIC_e2e_features.json").read_text())
    data = np.load(tmp_path / "SYNTHETIC_e2e_features.npz")
    assert doc["counts"]["kappa"]["valid"] == int(data["kappa_valid"].sum())
    assert doc["non_finite_valid_values"] == 0 and doc["units"]["kappa"].startswith("radian")
    core = doc["provenance"]["core"]
    assert core["native_extension"]["status"] == "resolved"
    assert core["native_extension"]["sha256"] != core["wrapper"]["sha256"] and "binary_sha256" not in core
    assert data["label"][0] == "SYNTH_A:10" and data["kappa_rad"].shape == (3, 12)
    assert (tmp_path / "SYNTHETIC_e2e_report.md").read_text().startswith("# Cα feature inspection report")
    # tampering with the adapter export is detected
    from topofold_inputs.features import features_from_export
    a_npz = tmp_path / "SYNTHETIC_e2e.npz"
    a_npz.write_bytes(a_npz.read_bytes() + b"\0")
    with pytest.raises(ValueError, match="SHA-256"):
        features_from_export(a_npz, tmp_path / "SYNTHETIC_e2e.json")
