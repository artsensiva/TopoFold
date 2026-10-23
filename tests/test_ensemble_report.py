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


def _member(run_id, n_res=8, n_frames=120, tau_center=math.pi, phi=0.0, seed=0, invalid=None, invalid_code=2,
            breaks=None, statuses=None, frame_offset=0):
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
    kr = np.where(kv, 0, 2).astype(np.uint8)
    tr = np.where(tv, 0, 2).astype(np.uint8)
    if invalid:
        if isinstance(invalid, dict):                      # {window: frame mask} for partial invalidity
            for w_, fm in invalid.items():
                kv[fm, w_] = False
                kr[fm, w_] = invalid_code
        else:
            kv[:, invalid] = False
            kr[:, invalid] = invalid_code
    residues = [{"chain": "A", "resid": i + 1, "icode": "", "resname": "ALA",
                 "status": (statuses or {}).get(i, "present")} for i in range(n_res)]
    return {"run_id": run_id, "coords": coords, "present": np.ones(n_res, bool),
            "frame_index": np.arange(n_frames) + frame_offset, "kappa": kappa, "kappa_valid": kv, "tau": tau,
            "tau_valid": tv, "kappa_reason": kr, "tau_reason": tr, "breaks": breaks or [], "residues": residues,
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
    groups = {"A": [_member("A_0", invalid=[3], breaks=[{"after_index": 3, "reasons": ["numbering_discontinuity"]}])],
              "B": [_member("B_0", seed=5)]}
    res = en.analyze(groups, OPTS)
    ex = res["rows"]["exclusions"]
    assert {(e["channel"], e["window"]) for e in ex} == {("kappa", 3), ("d_i_i2", 3)}
    r = [e for e in ex if e["channel"] == "kappa"][0]
    assert r["member"] == "A_0" and r["cause"] == "numbering_discontinuity" and r["location"] == "A:4→A:5"
    assert r["n_frames"] == 120 and r["source_frames"].startswith("0;1;2") and r["source_frames_truncated"]
    assert not any(r["channel"] == "kappa" and r["window"] == 3 for r in res["rows"]["distribution"])


def test_exclusion_causes_numbering_missing_and_degenerate():
    br = [{"after_index": 2, "reasons": ["numbering_discontinuity", "ca_ca_distance"], "n_frames": 4,
           "max_distance_angstrom": 7.5, "first_frames": [100, 101]}]
    a = _member("A_0", invalid=[2], breaks=br, frame_offset=100)
    deg = np.zeros(120, bool)
    deg[[5, 6, 7]] = True
    b = _member("B_0", seed=4, invalid={5: deg}, invalid_code=3)       # degenerate geometry in 3 frames only
    c = _member("C_0", seed=6, invalid=[4], invalid_code=1, statuses={4: "no_CA_atom"})
    res = en.analyze({"A": [a], "B": [b], "C": [c]}, OPTS)
    ex = res["rows"]["exclusions"]
    ka = [e for e in ex if e["channel"] == "kappa" and e["window"] == 2]
    assert {e["cause"] for e in ka} == {"numbering_discontinuity", "ca_ca_distance"}   # multiple causes kept
    assert all(e["location"] == "A:3→A:4" and e["member"] == "A_0" for e in ka)
    assert "max 7.50 Å" in [e for e in ka if e["cause"] == "ca_ca_distance"][0]["detail"]
    assert ka[0]["source_frames"].startswith("100;101")                                # source, not local, frames
    kb = [e for e in ex if e["channel"] == "kappa" and e["window"] == 5]
    assert len(kb) == 1 and kb[0]["cause"] == "degenerate_segment" and kb[0]["n_frames"] == 3
    assert kb[0]["source_frames"] == "5;6;7" and not kb[0]["source_frames_truncated"]
    kc = [e for e in ex if e["channel"] == "kappa" and e["window"] == 4]
    assert kc[0]["cause"] == "missing_coordinates" and "no_CA_atom" in kc[0]["detail"] and kc[0]["location"] == "A:5"
    # a break code with no adapter record is reported as insufficient provenance, not guessed
    d = _member("D_0", invalid=[6])
    ed = en.analyze({"D": [d], "E": [_member("E_0", seed=9)]}, OPTS)["rows"]["exclusions"]
    assert ed[0]["cause"] == "terminal_or_break" and "provenance insufficient" in ed[0]["detail"]


def test_exclusion_detail_does_not_change_included_windows_or_js():
    base = {"A": [_member("A_0", seed=1), _member("A_1", seed=2)], "B": [_member("B_0", seed=3)]}
    r1 = en.analyze(base, OPTS)
    flagged = {"A": [_member("A_0", seed=1, invalid=[3]), _member("A_1", seed=2)], "B": [_member("B_0", seed=3)]}
    r2 = en.analyze(flagged, OPTS)
    keep = lambda rows: {(r["channel"], r["window"], r["group_a"], r["group_b"]): r["pooled_js"] for r in rows  # noqa: E731
                         if not (r["window"] == 3 and r["channel"] in ("kappa", "d_i_i2"))}
    assert keep(r1["rows"]["between"]) == keep(r2["rows"]["between"])


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
    causes = {(r["cause"], r["location"]) for r in ex}
    assert ("numbering_discontinuity", "A:13→A:16") in causes
    assert ("missing_coordinates", "A:18") in causes
    assert ("chain_boundary", "A:18→B:10") in causes and ("missing_CA", "A:18→B:10") in causes   # both kept
    assert all(r["member"] and r["source_frames"] != "" for r in ex)
    rep = (tmp_path / "o1" / "report.md").read_text()
    assert "numbering_discontinuity" in rep and "missing_coordinates" in rep
    dist = list(csv.DictReader(open(tmp_path / "o1" / "tables/distribution.csv")))
    assert {r["group"] for r in dist} == {"G1", "G2"} and all(r["residues"].startswith(("A:", "B:")) for r in dist)
    m = json.loads((tmp_path / "o1" / "manifest.json").read_text())
    assert m["groups"]["G2"][0]["frame_index_first_last"] == [1, 2]


# --- corrections (2026-10-23): origin-invariant circular diagnostic, conditioning, degenerate series ----------------

def _ar1(n, phi, seed):
    rng = np.random.default_rng(seed)
    x = np.empty(n)
    x[0] = rng.standard_normal()
    for t in range(1, n):
        x[t] = phi * x[t - 1] + math.sqrt(1 - phi * phi) * rng.standard_normal()
    return x


def test_circular_autocorrelation_is_origin_invariant():
    theta = np.angle(np.exp(1j * (0.4 + 0.8 * _ar1(400, 0.85, 3))))
    base, _ = en.autocorrelation(theta, True, 30)
    for shift in (0.3, 1.7, -2.9, math.pi, 5.0):
        rot = np.angle(np.exp(1j * (theta + shift)))
        c, status = en.autocorrelation(rot, True, 30)
        assert status == "ok" and np.max(np.abs(c - base)) < 1e-12
    # the earlier average of separately normalized cos/sin autocorrelations depends on the origin
    def old(t):
        a = [en.autocorrelation(np.cos(t), False, 5)[0], en.autocorrelation(np.sin(t), False, 5)[0]]
        return 0.5 * (a[0] + a[1])
    assert abs(old(theta)[1] - old(np.angle(np.exp(1j * (theta + 1.0))))[1]) > 1e-3


def test_linear_autocorrelation_matches_direct_formula():
    x = _ar1(500, 0.6, 4)
    c, _ = en.autocorrelation(x, False, 5)
    u = x - x.mean()
    np.testing.assert_allclose(c, [np.dot(u[: len(u) - k], u[k:]) / np.dot(u, u) for k in range(6)], atol=1e-12)


def test_circular_summary_conditioning_and_rounding():
    balanced = np.array([0.0, 2 * math.pi / 3, -2 * math.pi / 3])
    s = en.circular_summary(balanced)
    assert s["center"] is None and s["spread"] is None and s["resultant_length"] < en.R_MIN
    assert "resultant length below" in s["status"]
    same = en.circular_summary(np.full(1000, 2.5))                 # R may round to slightly above 1
    assert same["status"] == "ok" and same["spread"] == 0.0 and same["resultant_length"] <= 1.0
    ok = en.circular_summary(np.array([3.1, -3.1, 3.05]))
    assert ok["status"] == "ok" and abs(abs(ok["center"]) - math.pi) < 0.1 and 0.99 < ok["resultant_length"] <= 1


def test_degenerate_series_have_explicit_reasons():
    one = en.temporal_diagnostics(np.array([1.0]), False, 50)
    assert one["lag1"] is None and one["status"] == one["decorrelation_status"] == "unavailable: fewer than 2 frames"
    const = en.temporal_diagnostics(np.full(100, 2.0), False, 50)
    assert const["status"] == "unavailable: constant series" and const["decorrelation_lag"] is None
    cconst = en.temporal_diagnostics(np.full(100, 3.0), True, 50)
    assert cconst["status"] == "unavailable: constant series"
    short = en.temporal_diagnostics(np.array([0.0, 0.1, 0.2, 0.3, 0.4]), False, 50)
    assert short["lags_available"] == 4 and short["status"] == "ok"
    slow = en.temporal_diagnostics(_ar1(60, 0.99, 5), False, 10)
    assert slow["decorrelation_status"] == "not reached within 10 lags" and slow["decorrelation_lag"] is None


def test_unavailable_values_reach_tables_and_report(tmp_path):
    a = _member("A_0", n_frames=3, seed=1)
    b = _member("B_0", n_frames=3, seed=2)
    for m in (a, b):                                                # τ window 2: exactly balanced angles, R = 0
        m["tau"][:, 2] = [0.0, 2 * math.pi / 3, -2 * math.pi / 3]
    single = _member("C_0", n_frames=1, seed=3)
    groups = {"A": [a], "B": [b], "C": [single]}
    res = en.analyze(groups, OPTS)
    d = [r for r in res["rows"]["distribution"] if r["channel"] == "tau" and r["window"] == 2 and r["group"] == "A"][0]
    assert d["center"] is None and "resultant length below" in d["status"]
    bt = [r for r in res["rows"]["between"] if r["channel"] == "tau" and r["window"] == 2
          and r["group_a"] == "A" and r["group_b"] == "B"][0]
    assert bt["center_difference"] is None and "center of A" in bt["center_difference_status"]
    assert bt["pooled_js"] == 0.0                                   # JS unaffected by the conditioning rule
    c = [r for r in res["rows"]["distribution"] if r["group"] == "C" and r["channel"] == "kappa"][0]
    assert c["spread"] is None and "fewer than 2 frames" in c["status"]
    t = [r for r in res["rows"]["temporal"] if r["group"] == "C"]
    assert all(r["status"] == "unavailable: fewer than 2 frames" for r in t)
    en.write_outputs(res, groups, {"name": "SYNTHETIC degenerate"}, OPTS, tmp_path / "o", None)
    rep = (tmp_path / "o" / "report.md").read_text()
    assert "unavailable: fewer than 2 frames" in rep and "beyond max lag" not in rep
    rows = list(csv.DictReader(open(tmp_path / "o" / "tables/between.csv")))
    assert any(r["center_difference"] == "" and r["center_difference_status"].startswith("unavailable") for r in rows)
