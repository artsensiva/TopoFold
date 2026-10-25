"""The Abl demonstration record (1OPJ chains A, B vs 1OPK chain A) and the correction record: consistency of committed files only."""

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "benchmarks/audit/structpair_abl_demo_2026-10-25/results"
REC = json.loads((DEMO / "abl_demo_record.json").read_text())
RUN = DEMO / "run"


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def test_inputs_are_the_recorded_files_and_the_header_facts_hold():
    kin = json.loads((ROOT / "benchmarks/audit/kinase_descriptors/run_spec.json").read_text())["structures"]
    for k, pid in (("A", "1OPJ"), ("B", "1OPK")):
        i = REC["inputs"][k]
        assert i["pdb"] == pid and i["sha256"] == kin[pid]["sha256"] and i["sha256_equals_kin_desc_1_record"]
        assert i["header"]["organism"] == "MUS MUSCULUS" and i["header"]["taxid"] == "10090"
        assert {d["accession"] for d in i["header"]["dbref"]} == {"P00520"} and i["rcsb"]["uniprot"] == ["P00520"]
    a, b = REC["inputs"]["A"]["header"], REC["inputs"]["B"]["header"]
    assert (a["space_group"], b["space_group"]) == ("P 1", "C 2 2 21") and (a["resolution_A"], b["resolution_A"]) == (1.75, 1.8)
    assert {h["name"] for h in a["het"]} == {"STI", "MYR", "CL"} and {h["name"] for h in b["het"]} == {"P16", "MYR", "GOL"}
    assert a["hetsyn"]["STI"] == "STI-571; IMATINIB" and b["hetsyn"]["P16"] == "PD166326"
    assert [x["resid"] for x in b["seqadv"] if "ENGINEERED" in x["comment"]] == [382] and not [x for x in a["seqadv"] if "ENGINEERED" in x["comment"]]
    assert "nucleotide" not in json.dumps([a["hetnam"], b["hetnam"]]).lower()


def test_every_input_check_passed_and_the_only_name_difference_is_the_documented_mutation():
    c = REC["input_checks"]
    assert c["remark_465_consistent_with_coordinates"] and c["no_altloc_or_icode_at_mapped_calpha"] and c["same_species_and_gene"]
    assert c["differences_all_documented_as_engineered_in_SEQADV"] and c["same_numbering_offset_to_database"]
    assert list(c["residue_name_differences"]) == ["382"] and c["residue_name_differences"]["382"]["1OPK:A"] == "ASN"
    assert all(all(v.values()) for v in c["rcsb_agrees_with_local_records"].values())
    assert "P00519" in REC["kin_desc_1_note"] and "not reopened" in REC["kin_desc_1_note"]


def test_partial_status_with_missing_not_zero_for_the_single_chain_entry():
    s = REC["result_summary"]
    assert s["within_entry_status"] == "partial" and s["n_pairs"] == {"cross": 2, "withinA": 1, "withinB": 0}
    assert s["positions_mapped"] == 282 and s["mapped_ranges"] == "248-529"
    rows = list(csv.DictReader(open(RUN / "chain_variation.csv")))
    assert len(rows) == 282 * 7
    assert all(r["withinB_n"] == "0" and r["withinB_min"] == "" and r["withinB_max"] == "" for r in rows)          # empty, never 0
    assert all(r["within_scope"] == "A_only" for r in rows if r["between_n"] != "0")
    assert all(r["withinA_n"] == "1" and r["withinA_max"] != "" for r in rows if r["channel"] == "disp")
    assert all(r["advisory_flag"] == "" for r in rows)                                    # no rule configured: nothing assessed, nothing true
    assert all("minimum difference" in r["advisory_not_assessed_reason"] for r in rows if r["between_n"] != "0")
    assert s["entry_B_within_cells_empty_not_zero"] and s["advisory_configured_channels"] == [] and s["advisory_flags_assessed"] == 0


def test_run_outputs_match_their_hashes_and_the_manifest():
    for name, h in REC["outputs_sha256"].items():
        assert _sha(RUN / name) == h, name
    man = json.loads((RUN / "manifest.json").read_text())
    assert man["schema"] == "topofold_inputs.structpair/3" and man["software"]["structpair_version"] == "0.3.0"
    assert man["within_entry_variation"]["advisory_semantics"] == "advisory/2" and man["inputs"]["B"]["chains"] == ["A"]
    for f in ("chain_variation.csv", "chains.csv", "windows.csv", "residues.csv", "exclusions.csv", "summary.json"):
        assert man["outputs_sha256"][f] == _sha(RUN / f), f
    assert (RUN / "chain_variation.png").stat().st_size < 5 * 1024 * 1024
    assert len(list(csv.DictReader(open(RUN / "exclusions.csv")))) == 0                   # no window excluded inside the mapping


def test_caption_states_the_limits_in_plain_words():
    cap = (RUN / "chain_variation_caption.md").read_text()
    for needle in ("missing, not zero", "nothing here can be attributed to imatinib, nucleotide state, mutation or resistance", "not human",
                   "D363N", "Neither entry\ncontains a nucleotide", "not independent experimental replicates", "not a noise floor",
                   "REMARK 465", "1OPJ chain A: 530-534", "P 1", "C 2 2 21"):
        assert needle in cap, needle
    assert "2F4J" not in cap and "separated" not in cap.lower()


def test_correction_record_shows_the_failing_before_and_passing_after_cases():
    r = json.loads((ROOT / "benchmarks/audit/structpair_advisory_correction_2026-10-25/results/regression_before_after.json").read_text())
    assert r["before"]["code"]["structpair_version"] == "0.2.0" and r["after"]["code"]["structpair_version"] == "0.3.0"
    assert r["counts"] == {"before": {"passed": 6, "failed": 13}, "after": {"passed": 19, "failed": 0}} and not r["still_failing_after"]
    assert len(r["failed_before_and_passing_after"]) == 13
    assert "test_identical_chains_with_all_values_zero_are_never_flagged_by_default" in r["failed_before_and_passing_after"]
    assert any(t.startswith("test_nan_and_infinity_are_rejected_at_parsing") for t in r["failed_before_and_passing_after"])
