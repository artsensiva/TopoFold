"""Tests of the ABL-MYR-1b protocol writer: measured numbers come from the a2 results JSON, wording restrictions hold, and the
lineage files are the preserved ones. Pure Python; no trajectory and no study archive is involved."""

import copy
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_myr_1b"))
import write_protocol_1b as W  # noqa: E402

RES = json.loads(W.A2_RESULTS.read_text())


def test_measured_numbers_equal_the_a2_results_json():
    m = W.measured(RES)
    r2 = RES["R2_production_start_geometry"]
    assert m["sets"]["Set1"]["between_Myr_noMyr"]["ca_rmsd_A"]["median"] == r2["Set1"]["between_Myr_noMyr"]["ca_rmsd_A"]["median"]
    assert m["reader_agreement"]["runs_ok"] == RES["a2"]["agreement_ok"] == 40 and m["runs_validated"] == 40
    for s, d in m["sets"].items():
        w = max(d["within_Myr"]["ca_rmsd_A"]["median"], d["within_noMyr"]["ca_rmsd_A"]["median"])
        assert d["start_contrast_E_A"] == pytest.approx(d["between_Myr_noMyr"]["ca_rmsd_A"]["median"] - w)
    assert all(m["criteria"].values()) and m["outcome"] == "comparable_for_1b"


def test_protocol_is_built_from_the_results_and_rendered_markdown_matches():
    p = W.build(RES)
    assert p["post_hoc"] and "No trajectory has been opened or scored" in p["status"]
    assert W.render_md(p) == W.PROTOCOL_MD.read_text() and json.loads(W.PROTOCOL.read_text()) == p
    assert f"{p['measured_comparability_evidence']['sets']['Set1']['between_Myr_noMyr']['ca_rmsd_A']['median']:.3f}" in W.PROTOCOL_MD.read_text()


def test_restricted_wording_and_unchanged_thresholds():
    p = W.build(RES)
    labels = " ".join(o["label"] for o in p["outcomes_ordered_with_restricted_labels"])
    assert "reproduced" not in labels and "independently" not in labels
    assert all("post-hoc" in o["label"] for o in p["outcomes_ordered_with_restricted_labels"])
    orig = json.loads((ROOT / "benchmarks/audit/abl_myr_1/run_spec.json").read_text())
    steps = {o["step"]: o["rule"] for o in orig["outcomes_ordered"]}
    assert steps[3] == "|R_T \\ N2(R_B)| >= 3" and p["outcomes_ordered_with_restricted_labels"][2]["rule"].startswith(steps[3])
    assert p["unchanged_from_original"]["budget"] == {"download_bytes_max": orig["budget"]["download_bytes_max"],
                                                      "storage_bytes_max": orig["budget"]["storage_bytes_max"],
                                                      "cpu_hours_max": orig["budget"]["cpu_hours_max"]}
    for k in ("450", "5 %", "8 runs", "1000 ns"):
        assert any(k in (v if isinstance(v, str) else json.dumps(v)) for v in orig["gates"].values()), k
    assert "independent replication" in " ".join(p["interpretation_restrictions"])
    assert "No trajectory" in p["review_gate"] and p["implementation_plan"]["required_before_trajectories"]


def test_protocol_states_what_is_not_supported():
    p = W.build(RES)
    txt = " ".join(p["what_the_evidence_does_not_support"] + p["interpretation_restrictions"]).lower()
    for needle in ("identical", "no confound", "independent replicates", "causal", "alphai"):
        assert needle in txt, needle
    assert "Ayaz" in p["separate_from"] and "HDX" in p["separate_from"]


def test_writer_refuses_when_a2_did_not_pass(tmp_path, monkeypatch):
    bad = copy.deepcopy(RES)
    bad["outcome"] = "not_comparable"
    f = tmp_path / "r.json"
    f.write_text(json.dumps(bad))
    monkeypatch.setattr(W, "A2_RESULTS", f)
    with pytest.raises(SystemExit, match="no ABL-MYR-1b protocol"):
        W.cmd_write()


def test_lineage_files_exist_and_freeze_matches_when_present():
    for v in W.LINEAGE.values():
        assert (ROOT / v).is_file(), v
    if not W.FREEZE.exists():
        pytest.skip("protocol freeze not yet written")
    f = json.loads(W.FREEZE.read_text())
    assert f["protocol_sha256"] == W.sha256(W.PROTOCOL) and f["protocol_md_sha256"] == W.sha256(W.PROTOCOL_MD)
    assert f["lineage_sha256"] == {k: W.sha256(ROOT / v) for k, v in W.LINEAGE.items()}
    assert f["trajectories_opened"] is False and f["xtc_extracted"] is False
