"""Consistency of the 2026-10-27 correction records for the Abl response intake (committed results only; no raw files needed)."""
import json
from pathlib import Path

D = Path(__file__).resolve().parent.parent / "benchmarks/audit/abl_response_intake_2026-10-27"
coh = json.loads((D / "results/endpoint_cohort_2026-10-27.json").read_text())
cor = json.loads((D / "results/abl_response_corrections_2026-10-27.json").read_text())


def test_each_record_belongs_to_exactly_one_cohort_and_assays_are_not_pooled():
    recs = coh["records"]
    assert len({r["id"] for r in recs}) == len(recs) == 144
    assert sum(v["records"] for v in coh["cohorts"].values()) == 144
    by = {}
    for r in recs:
        by.setdefault(r["cohort"], set()).add(r["inhibitor"])
    assert by["primary_candidate_bafs_viability"] == {"bosutinib", "dasatinib", "imatinib", "nilotinib", "ponatinib"}
    assert by["secondary_axitinib_cell_ic50"] == {"axitinib"} and by["secondary_binding_davis_kd"] == {"erlotinib", "gefitinib"}


def test_primary_candidate_counts_follow_from_the_records_and_the_exclusion():
    prim = [r for r in coh["records"] if r["cohort"] == "primary_candidate_bafs_viability"]
    kept = [r for r in prim if r["exclusion"] is None]
    assert [r["id"] for r in prim if r["exclusion"]] == ["ponatinib:G250E"]
    P = coh["primary_candidate_after_exclusion"]
    assert P["records"] == len(kept) == sum(P["records_per_position"].values())
    assert P["positions"] == len({r["position"] for r in kept}) and P["substitutions"] == len({r["mutation"] for r in kept})
    flagged = set(coh["primary_flagged_censored_or_bound_records"])
    assert coh["primary_sensitivity_subset_without_flagged"]["records"] == P["records"] - len(flagged)


def test_sheet_and_official_si_agree_and_wild_type_denominators_reproduce_the_sheet():
    sv = coh["si_vs_github_sheet"]
    assert sv["six_inhibitors_ic50_equal"] == sv["six_inhibitors_records"] == 131 and sv["table1_ddg_equal_within_0.006"] == sv["all_records"]
    for r in coh["records"]:
        if r["id"] != "ponatinib:G250E":
            assert abs(r["ddg_recomputed_kcal_mol"] - r["ddg_sheet_kcal_mol"]) <= 0.05


def test_permissions_are_unresolved_and_the_implementation_is_not_run():
    assert set(cor["permissions_status"].values()) <= {"unresolved", "unresolved for the data files (the article text is CC BY 4.0)"}
    assert cor["implementation_decision"]["data_dependent_feature_export"] == "not run"
    assert cor["downloads"]["new_bytes"] < cor["downloads"]["correction_ceiling_bytes"]
    assert any("terms" in m for m in cor["missing_evidence"])
