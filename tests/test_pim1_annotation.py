"""Consistency of the recorded PIM1 structural-role annotation (benchmarks/audit/pim1_structural_annotation_2026-10-26): checks the record, not the biology."""

import hashlib
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "benchmarks/audit/pim1_structural_annotation_2026-10-26"
REC = D / "results/pim1_annotation.json"
pytestmark = pytest.mark.skipif(not REC.is_file(), reason="annotation record not present")


def _d():
    return json.loads(REC.read_text())


def test_reference_anchors_and_helper_checks_are_recorded_and_the_rules_were_met():
    d = _d()
    assert d["reference"]["accession"] == "P11309" and d["reference"]["length"] == 313 and d["reference"]["matches_selection_record"] is True
    assert d["alphaC_glu_rule"]["qualifying"] == [89] and d["alphaC_glu_rule"]["established"] is True and d["dfg_window"]["hits"] == [186]
    h = d["kinase_roles_helper"]
    assert h["ran"] and h["all_accepted"] and h["allow_role_mutations"] is False and h["anchors_declared"] == {"beta3_Lys": 67, "alphaC_Glu": 89, "HRD_H": 165, "DFG_Asp": [168, 290]}
    for e, v in h["per_structure"].items():
        assert v["accepted"] and v["problems"] == [] and v["all_unique"] and v["anchor_order"]["passed"], e
    sel = json.loads((ROOT / "benchmarks/audit/kin_pair_1/results/KIN-PAIR-1_selection.json").read_text())
    shas = {p[k]: p[f"sha256_{k}"] for p in sel["panel"] for k in "AB"}
    assert all(d["structures"][e]["sha256"] == shas[e] for e in ("8R0Q", "6NO9", "6VRV", "6L17"))
    assert all(s["author_minus_reference_offsets"] == [0] for s in d["structures"].values())


def test_every_new_window_contains_residue_261_and_is_annotated_by_complete_identifiers():
    d = _d()
    for pid, p in d["pairs"].items():
        assert p["n_new_windows"] == {"3": 3, "4": 4} and len(p["new_windows"]) == 7
        for w in p["new_windows"]:
            assert len(w["residues"]) == int(w["window_residues"]) and 261 in w["reference_indices"] and w["same_reference_index_in_both"]
            assert any(r["A_modification"] or r["B_modification"] for r in w["residues"]) and w["name_difference_flag"] and w["contains_modified_residue"]
            assert w["inside_uniprot_kinase_domain"] and not w["in_a_top10_of_any_metric"]
            for m, v in w["values"].items():
                rk = w["ranks_in_new_run_of_windows_with_values"][m]
                assert 1 <= rk["rank"] <= rk["of"] and v >= 0


def test_previous_descriptive_conclusions_are_reproduced_and_unchanged_and_the_record_does_not_overclaim():
    d = _d()
    for p in d["pairs"].values():
        ds = p["descriptive"]
        assert ds["old_reapplication_reproduces_recorded"] is True and ds["classification_changed_by_the_added_windows"] is False
        assert ds["top10_old_still_in_top10_and_no_new_entrant"] is True and "displacement" in ds["not_recomputed"]
    assert "not evidence of a functional state" in d["uniprot_annotation"]["note"] and "do not isolate" in d["scope"]
    assert any("conventional names" in u for u in d["unverified"]) and any("activation loop" in u for u in d["unverified"])
    raw_hash = d["uniprot_annotation"]["raw_sha256"]
    assert len(raw_hash) == 64 and hashlib.sha256(b"").hexdigest() != raw_hash
