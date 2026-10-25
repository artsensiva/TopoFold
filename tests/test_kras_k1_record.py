"""Consistency of the committed KRAS K1 endpoint record (S1) with the frozen rules; uses only committed files."""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
K1 = ROOT / "benchmarks/audit/kras_k1"
sys.path.insert(0, str(K1))
import k1_endpoint as K  # noqa: E402

REC = json.loads((K1 / "results/KRAS-K1_endpoint.json").read_text())
P = json.loads((K1 / "protocol_params.json").read_text())
EP = P["endpoint"]


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def test_freeze_hashes_still_match_the_frozen_files():
    f = K.verify_freeze()
    assert REC["freeze_record_sha256"] == _sha(K.FREEZE)
    assert f["state"].startswith("frozen before any crystal endpoint was computed")
    assert all(v["matches_spec"] for v in REC["inputs"].values())


def test_recorded_positives_follow_the_frozen_rule_from_the_recorded_arrays():
    pr = REC["per_residue"]
    thr, pos, _ = K.compute_positive(np.array(pr["b_mean"]), np.array(pr["w_max"]), EP["min_displacement_A"], EP["within_state_factor"])
    assert thr.tolist() == pytest.approx(pr["threshold"], abs=1e-12)
    assert pos.tolist() == pr["positive"]
    assert [r for r, p in zip(pr["resid"], pr["positive"]) if p] == REC["E"] and REC["n_E"] == len(REC["E"])
    assert pr["resid"] == REC["scored_residues"] and REC["n_scored"] == 153


def test_recorded_gate_and_decision_follow_from_the_recorded_endpoint():
    g = K.gate(REC["E"], REC["n_scored"], EP["gate"], EP["segments"]["merge_gap"])
    assert g == REC["gate"]
    assert REC["decision"] == ("K1_ENDPOINT_ADEQUATE" if g["passes"] else "K1_ENDPOINT_INADEQUATE")
    assert REC["md_opened"] is False
    assert REC["evidence_class"] == EP["evidence_class"]


def test_sensitivity_sets_are_reported_only():
    s = REC["sensitivity"]
    assert set(s) == {"primary_pair_only", "median_over_pairs"}
    assert all({"E", "n", "segments", "jaccard_with_E"} <= set(v) for v in s.values())
    assert all(v["n"] == len(v["E"]) for v in s.values())


def test_figure_and_report_are_generated_from_this_record():
    prov = json.loads((K1 / "figures/KRAS-K1_endpoint_residues_2026-10-25_provenance.json").read_text())
    assert prov["inputs_sha256"]["results/KRAS-K1_endpoint.json"] == _sha(K1 / "results/KRAS-K1_endpoint.json")
    assert prov["output"]["sha256"] == _sha(K1 / "figures/KRAS-K1_endpoint_residues_2026-10-25.png")
    assert prov["generator"]["sha256"] == _sha(K1 / "plot_k1_endpoint.py")
    assert prov["output"]["bytes"] < 5 * 1024 * 1024
    rep = (K1 / "results/KRAS-K1_endpoint_report.md").read_text()
    assert REC["decision"] in rep and REC["frozen_commit"] in rep and _sha(K1 / "results/KRAS-K1_endpoint.json") in rep


def test_independent_verification_agrees():
    v = json.loads((K1 / "results/KRAS-K1_endpoint_verification.json").read_text())
    assert v["agrees"] and v["E_equal"] and v["E_independent"] == REC["E"] and v["md_opened"] is False
    assert max(v["max_abs_difference_A"].values()) < v["tolerance_A"]
