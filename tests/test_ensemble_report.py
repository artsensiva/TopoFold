"""Tests for topofold_inputs.ensemble (descriptive group comparison). SYNTHETIC inputs only."""

import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs import ensemble as en  # noqa: E402


def _member(run_id, n_res=8, n_frames=120, tau_center=math.pi, phi=0.0, seed=0, invalid=None):
    rng = np.random.default_rng(seed)
    z = np.empty((n_frames, n_res))
    z[0] = rng.standard_normal(n_res)
    for t in range(1, n_frames):
        z[t] = phi * z[t - 1] + math.sqrt(1 - phi * phi) * rng.standard_normal(n_res)
    coords = np.zeros((n_frames, n_res, 3))
    coords[:, :, 0] = 3.8 * np.arange(n_res)
    coords[:, :, 1] = 0.5 * (-1) ** np.arange(n_res) + 0.05 * z
    kappa = np.full((n_frames, n_res), np.nan)
    kappa[:, 1:-1] = 1.5 + 0.1 * z[:, 1:-1]
    tau = np.full((n_frames, n_res - 1), np.nan)
    tau[:, 1:-1] = np.angle(np.exp(1j * (tau_center + 0.3 * z[:, 1:-2])))
    kv, tv = np.isfinite(kappa), np.isfinite(tau)
    if invalid:
        kv[:, invalid] = False
    residues = [{"chain": "A", "resid": i + 1, "icode": "", "resname": "ALA"} for i in range(n_res)]
    return {"run_id": run_id, "coords": coords, "present": np.ones(n_res, bool), "frame_index": np.arange(n_frames),
            "kappa": kappa, "kappa_valid": kv, "tau": tau, "tau_valid": tv, "residues": residues,
            "provenance": {"source_files": [], "core": {"native_extension": {"status": "SYNTHETIC", "sha256": None}},
                           "frame_identity": "SYNTHETIC", "units": {}}}


OPTS = {"bins": 36, "max_lag": 50, "pbc": "none"}


def test_mapping_labels_rows_and_single_member_handling():
    groups = {"A": [_member("A_0", seed=1), _member("A_1", seed=2)], "B": [_member("B_0", seed=3)]}
    res = en.analyze(groups, OPTS)
    rows = res["rows"]
    k = [r for r in rows["distribution"] if r["channel"] == "kappa"]
    assert {r["residues"] for r in k} == {f"A:{i}-{i + 2}" for i in range(1, 7)}          # 3-residue windows
    t = [r for r in rows["distribution"] if r["channel"] == "tau"]
    assert {r["residues"] for r in t} == {f"A:{i}-{i + 3}" for i in range(1, 6)}          # 4-residue windows
    rb = [r for r in rows["replicas"] if r["group"] == "B"]
    assert all(r["member_center_spread"] is None and r["within_member_js_median"] is None for r in rb)
    ra = [r for r in rows["replicas"] if r["group"] == "A"]
    assert all(r["within_member_js_median"] is not None for r in ra)
    assert len(rows["between"]) == len(k) // 2 + len(t) // 2 + 6 + 5                      # 4 channels, one group pair


def test_exclusions_reported_not_silent():
    groups = {"A": [_member("A_0", invalid=[3])], "B": [_member("B_0", seed=5)]}
    res = en.analyze(groups, OPTS)
    ex = res["rows"]["exclusions"]
    assert {(e["channel"], e["window"]) for e in ex} == {("kappa", 3), ("d_i_i2", 3)}
    assert json.loads(ex[0]["invalid_frames_by_member"]) == {"A_0": 120}
    assert not any(r["channel"] == "kappa" and r["window"] == 3 for r in res["rows"]["distribution"])


def test_identity_mismatch_rejected():
    b = _member("B_0")
    b["residues"][2]["resid"] = 99
    with pytest.raises(ValueError, match="residue identity"):
        en.analyze({"A": [_member("A_0")], "B": [b]}, OPTS)


def test_circular_center_near_pi_and_wrapped_difference():
    groups = {"A": [_member("A_0", tau_center=math.pi - 0.05, seed=7)],
              "B": [_member("B_0", tau_center=-math.pi + 0.05, seed=8)]}
    res = en.analyze(groups, OPTS)
    d = [r for r in res["rows"]["distribution"] if r["channel"] == "tau" and r["group"] == "A"]
    assert all(abs(abs(r["center"]) - (math.pi - 0.05)) < 0.1 for r in d)
    bt = [r for r in res["rows"]["between"] if r["channel"] == "tau"]
    assert all(abs(r["center_difference"] + 0.1) < 0.15 for r in bt)                     # wrapped, not ~2π


def test_temporal_diagnostics_track_correlation():
    fast = en.temporal_diagnostics(np.random.default_rng(0).standard_normal(5000), False, 50)
    assert abs(fast["lag1"]) < 0.05 and fast["decorrelation_lag"] == 1
    rng = np.random.default_rng(1)
    x = np.empty(20000)
    x[0] = rng.standard_normal()
    for t in range(1, len(x)):
        x[t] = 0.9 * x[t - 1] + math.sqrt(1 - 0.81) * rng.standard_normal()
    slow = en.temporal_diagnostics(x, False, 50)
    assert abs(slow["lag1"] - 0.9) < 0.02
    assert slow["decorrelation_lag"] == math.ceil(-1 / math.log(0.9)) or abs(slow["decorrelation_lag"] - 10) <= 1
    circ = en.temporal_diagnostics(np.angle(np.exp(1j * (math.pi + 0.3 * x))), True, 50)
    assert circ["lag1"] > 0.8


def test_outputs_are_deterministic(tmp_path):
    groups = {"A": [_member("A_0", seed=1), _member("A_1", seed=2)], "B": [_member("B_0", seed=3)]}
    for d in ("r1", "r2"):
        en.write_outputs(en.analyze(groups, OPTS), groups, {"name": "SYNTHETIC"}, OPTS, tmp_path / d, None)
    for f in ["tables/distribution.csv", "tables/replicas.csv", "tables/between.csv", "tables/temporal.csv",
              "tables/exclusions.csv", "distributions.npz", "report.md"]:
        assert (tmp_path / "r1" / f).read_bytes() == (tmp_path / "r2" / f).read_bytes(), f
    m = json.loads((tmp_path / "r1" / "manifest.json").read_text())
    assert "no p-values" in m["statement"]
    rep = (tmp_path / "r1" / "report.md").read_text().lower()
    assert "p-value" in rep and "significant" not in rep


def _adapter_python():
    cand = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
    if not Path(cand).is_file():
        return None
    return cand if subprocess.run([cand, "-c", "import MDAnalysis"], capture_output=True).returncode == 0 else None


def test_end_to_end_with_adapter_and_reproducibility(tmp_path):
    py = _adapter_python()
    if py is None:
        pytest.skip("no interpreter with MDAnalysis available")
    sys.path.insert(0, str(ROOT / "python" / "examples"))
    from ca_adapter_example import write_synthetic_pdb
    pdb = tmp_path / "SYNTHETIC_group.pdb"
    write_synthetic_pdb(pdb, n_models=3)   # more models push Cα–Cα above 4.2 Å (geometric breaks)
    spec = {"name": "SYNTHETIC end-to-end", "groups": {
        "G1": [{"type": "mdanalysis", "topology": str(pdb), "altloc": "highest_occupancy"},
               {"type": "mdanalysis", "topology": str(pdb), "altloc": "highest_occupancy"}],
        "G2": [{"type": "mdanalysis", "topology": str(pdb), "altloc": "highest_occupancy", "start": 1}]}}
    sp = tmp_path / "spec.json"
    sp.write_text(json.dumps(spec))
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"))
    for d in ("o1", "o2"):
        r = subprocess.run([sys.executable, "-m", "topofold_inputs.ensemble", str(sp), "--out", str(tmp_path / d),
                            "--adapter-python", py], env=env, capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
    for f in ["tables/distribution.csv", "tables/between.csv", "tables/exclusions.csv", "distributions.npz", "report.md"]:
        assert (tmp_path / "o1" / f).read_bytes() == (tmp_path / "o2" / f).read_bytes(), f
    ex = list(csv.DictReader(open(tmp_path / "o1" / "tables/exclusions.csv")))
    assert ex, "the SYNTHETIC fixture has gaps and a missing Cα; windows touching them must be excluded"
    dist = list(csv.DictReader(open(tmp_path / "o1" / "tables/distribution.csv")))
    assert {r["group"] for r in dist} == {"G1", "G2"} and all(r["residues"].startswith(("A:", "B:")) for r in dist)
    m = json.loads((tmp_path / "o1" / "manifest.json").read_text())
    assert m["groups"]["G2"][0]["frame_index_first_last"] == [1, 2]
