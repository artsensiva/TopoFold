"""Consistency of the DMS candidate screen records (committed results only; the git-ignored raw files are not needed)."""
import json
from pathlib import Path

D = Path(__file__).resolve().parent.parent / "benchmarks/audit/dms_candidate_screen_2026-10-27/results"
met = json.loads((D / "met_screen.json").read_text())
scr = json.loads((D / "dms_candidate_screen.json").read_text())


def test_met_numbering_and_counts_are_consistent():
    m = met["mapping"]
    assert m["wild_type_identity_all_positions"] and m["dms_positions"] == m["uniprot_residue_range"][1] - m["uniprot_residue_range"][0] + 1
    assert len(met["conditions"]) == 12 and met["rows_with_non_numeric_effect"] == 0
    for k, v in met["per_condition"].items():
        assert v["paired_and_not_DMSO_LOF"] + v["DMSO_LOF_among_paired"] == v["paired_with_DMSO"]
        assert v["missense_variants_with_effect"] <= met["possible_missense_per_condition"] and v["positions"] <= m["dms_positions"]


def test_engineered_template_residues_are_the_only_identity_mismatches():
    mism = {p: v["residues_differing_from_P08581_in_span"] for p, v in met["structures"].items() if v["residues_differing_from_P08581_in_span"]}
    assert set(mism) == {"3QTI", "3DKC"}
    assert {n for c in mism["3QTI"].values() for n, _, _ in c} == {1194, 1234} and {n for c in mism["3DKC"].values() for n, _, _ in c} == {1194, 1234, 1235}


def test_download_ceiling_and_decision_record():
    assert met["download_ledger"]["total_bytes_upper_bound"] < scr["ceiling_bytes"] == 20_000_000
    assert len(scr["candidates"]) == 3 and scr["decision"]["value"].startswith("proceed")
    assert scr["candidates"][1]["verdict"].startswith("does not qualify") and scr["candidates"][2]["verdict"].startswith("not evaluated")


def test_hauser_supplement_and_github_copies_stay_unresolved():
    st = scr["hauser_permissions_closure"]["recorded_status"]
    assert "unresolved" in st["official_supplement_PDF"] and "unresolved" in st["GitHub_copies_and_prepared_models"] and "unresolved" in st["third_party_source_values"]
    assert "No legal conclusion" in scr["hauser_permissions_closure"]["note"]
