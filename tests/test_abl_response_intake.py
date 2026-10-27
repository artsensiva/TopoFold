"""Consistency of the Hauser/Aldeghi Abl response intake record (reads committed results only; the git-ignored raw files are not needed)."""
import json
import re
import zipfile
from pathlib import Path

import pytest

D = Path(__file__).resolve().parent.parent / "benchmarks/audit/abl_response_intake_2026-10-27"


def _load(name):
    return json.loads((D / "results" / name).read_text())


def test_endpoint_counts_are_internally_consistent():
    t = _load("inputs_inspection.json")["endpoint_table"]
    assert t["rows_with_full_values"] + t["rows_without_values"] == t["rows_in_sheet"]
    assert sum(t["inhibitors"].values()) == t["rows_with_full_values"]
    assert sum(t["mutations_per_site"].values()) == t["mutations_with_values"]
    assert len(t["mutations_per_site"]) == t["sites_with_values"]
    assert t["sites_with_any_resistant_record"] == len(t["sites_with_any_resistant_record_list"])
    assert t["wildtype_residue_matches_P00519_1a_numbering"]["matches"] == t["mutations_with_values"]


def test_every_ddg_but_the_flagged_row_follows_ic50_and_the_implied_wild_type():
    t = _load("inputs_inspection.json")["endpoint_table"]
    flagged = t["rows_where_ddg_exp_disagrees_with_ic50_column"]
    assert [(r["inhibitor"], r["mutation"]) for r in flagged] == [("ponatinib", "G250E")]
    for k, v in t["implied_wt_ic50_consistency_per_inhibitor"].items():
        assert v["max_abs_residual_kcal"] < 0.02 or k == "ponatinib"


def test_deposit_has_wild_type_models_only_and_cas_differ_by_template():
    ins = _load("inputs_inspection.json")
    assert ins["mutant_structures_in_deposit"] == {"files_matching_mutant_or_mutation_in_names": [], "trajectory_files": []}
    assert len(ins["structures"]) == 8
    assert all(v["mismatches_to_P00519_1a"] == [] for v in ins["structures"].values())
    assert all(v["mutation_wt_identity_present"] == v["of"] for v in ins["encoded_wt_identity"].values())
    assert ins["ca_trace_comparison"]["max_A"] > 100 * ins["ca_trace_comparison"]["coordinate_rounding_floor_A"]


def test_record_keeps_the_download_ceiling_and_the_unresolved_terms():
    r = _load("abl_response_intake_record.json")
    assert r["download_ledger"]["scientific_download_total_bytes_upper_bound"] < 50_000_000
    assert r["permissions"]["redistribution_of_derived_files"]["status"] == "unknown"
    assert r["decision"]["value"] == "conditional-go" and r["sent_externally"] == "nothing"


def test_minimal_xlsx_reader_reads_a_generated_workbook(tmp_path):
    import sys

    sys.path.insert(0, str(D))
    from xlsx_min import read_xlsx

    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    p = tmp_path / "t.xlsx"
    with zipfile.ZipFile(p, "w") as z:
        z.writestr("xl/workbook.xml", f'<workbook {ns}><sheets><sheet name="S" sheetId="1" r:id="rId1"/></sheets></workbook>')
        z.writestr("xl/_rels/workbook.xml.rels", '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr("xl/sharedStrings.xml", f"<sst {ns}><si><t>TKI</t></si></sst>")
        z.writestr("xl/worksheets/sheet1.xml", f'<worksheet {ns}><sheetData><row r="1"><c r="A1" t="s"><v>0</v></c><c r="C1"><v>2.5</v></c></row></sheetData></worksheet>')
    assert read_xlsx(p) == {"S": [["TKI", None, 2.5]]}
