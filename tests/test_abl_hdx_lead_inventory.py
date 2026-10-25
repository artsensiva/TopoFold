"""ABL HDX lead inventory: helpers on SYNTHETIC input; the recorded inventory against its own record. No network."""

import io
import json
import random
import sys
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "benchmarks/audit/abl_hdx_lead_2026-10-25"
sys.path.insert(0, str(D))
import inventory_hdx_lead as H  # noqa: E402

REC = json.loads((D / "results/hdx_inventory.json").read_text())
DOC = (ROOT / "docs/ABL_HDX_PUBLIC_DATA_INVENTORY_2026-10-25.md").read_text()


def test_mhp_of_a_known_peptide():
    # monoisotopic [M+H]+ of GG: 2 x 57.02146 + 18.010565 + 1.007276
    assert H.mhp("GG") == pytest.approx(133.060761, abs=1e-6)


def test_parse_export_keeps_the_preamble_and_skips_empty_rows():
    p = H.parse_export(",,,,,\r\nSequence,Start,End,MHP,MaxUptake,RT\r\nGG,1,2,133.0608,1,3.5\r\n,,,,,\r\n")
    assert p["preamble"] == [["", "", "", "", "", ""]] and p["header"][0] == "Sequence" and p["rows"] == [["GG", "1", "2", "133.0608", "1", "3.5"]]


def test_describe_export_on_a_synthetic_reference_with_offset_and_variant_peptides():
    rng = random.Random(5)
    base = list("".join(rng.choice("ACDEFGHIKLMNQRSTVWY") for _ in range(500)))                 # unique 8-mers
    base[343], base[471] = "A", "M"                              # residues 344 and 472 as in the real reference
    ref = "".join(base)
    a344v = ref[:343] + "V" + ref[344:]
    wt = ref[100:108]                                            # file numbering = reference index - 90
    var = a344v[340:349]                                         # contains the substituted residue
    rows = [[wt, "11", "18", f"{H.mhp(wt):.4f}", str(len(wt) - 1 - wt.count("P")), "4.0"],
            [var, "251", "259", f"{H.mhp(var):.4f}", str(len(var) - 1 - var.count("P")), "5.0"]]
    d = H.describe_export({"preamble": [], "header": ["Sequence", "Start", "End", "MHP", "MaxUptake", "RT"], "rows": rows}, ref)
    assert d["n_peptides"] == 2 and d["length_equals_end_minus_start_plus_1"] == 2
    rc = d["reference_check"]
    assert rc["matched_in_wild_type"] == 1 and rc["matched_only_in_a_variant"] == {"A344V": 1, "M472I": 0} and rc["unmatched"] == 0
    assert d["mhp_is_monoisotopic_MH_plus"]["rows_within_0.001"] == 2 and d["maxuptake_equals_len_minus_1_minus_prolines"] == 2
    assert rc["offsets_seen"] == [90]


def test_locate_reports_the_offset_the_variant_and_no_match():
    assert H.locate("DEF", "ABCDEFGH", 1, {}) == {"matched_in": "wild type", "offset": 3}
    assert H.locate("DVF", "ABCDEFGH", 1, {"X1V": "ABCDVFGH"}) == {"matched_in": "X1V", "offset": 3}
    assert H.locate("QQQ", "ABCDEFGH", 1, {}) == {"matched_in": None, "offset": None}


def test_file_shape_detects_line_endings_and_rows():
    s = H.file_shape(b",,,,,\r\nSequence,Start,End,MHP,MaxUptake,RT\r\nGG,1,2,1,1,2\r\nAA,3,4,1,1,2\r\n")
    assert s["n_lines"] == 4 and s["line_ending"] == "CRLF" and s["trailing_newline"] and s["last_data_row"] == "AA,3,4,1,1,2"
    assert H.file_shape(b"a\nb\nc")["line_ending"] == "mixed or LF" and H.file_shape(b"a\nb\nc")["trailing_newline"] is False


def _xlsx(sheet_names):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        sheets = "".join(f'<sheet name="{n}" sheetId="{i + 1}" r:id="rId{i + 1}"/>' for i, n in enumerate(sheet_names))
        z.writestr("xl/workbook.xml", f"<workbook><sheets>{sheets}</sheets></workbook>")
        z.writestr("xl/_rels/workbook.xml.rels", "<Relationships>" + "".join(
            f'<Relationship Id="rId{i + 1}" Type="x" Target="worksheets/sheet{i + 1}.xml"/>' for i in range(len(sheet_names))) + "</Relationships>")
        z.writestr("xl/sharedStrings.xml", "<sst><si><t>peptide uptake</t></si><si><t>Mean</t></si></sst>")
        for i in range(len(sheet_names)):
            z.writestr(f"xl/worksheets/sheet{i + 1}.xml", '<worksheet><dimension ref="A1:B2"/><sheetData><row r="1"><c r="A1" t="s"><v>0</v></c></row>'
                       '<row r="2"><c r="A2" t="s"><v>1</v></c><c r="B2"><v>3.5</v></c></row></sheetData></worksheet>')
    return buf.getvalue()


def test_xlsx_inventory_lists_sheets_headers_and_flags_hdx_terms():
    inv = H.xlsx_inventory(_xlsx(["Figure4F-Tm", "FigureS4-uptake"]))
    assert inv["n_sheets"] == 2 and [s["name"] for s in inv["sheets"]] == ["Figure4F-Tm", "FigureS4-uptake"]
    assert inv["sheets"][0]["dimension"] == "A1:B2" and "peptide uptake" in inv["sheets"][0]["header_strings"]
    assert H.hdx_terms(inv) == ["Figure4F-Tm", "FigureS4-uptake"]      # the header string contains 'peptide'/'uptake'
    plain = H.xlsx_inventory(_xlsx(["Figure4F-Tm"]))
    plain["sheets"][0]["header_strings"] = ["Imatinib", "WT"]
    assert H.hdx_terms(plain) == []


def test_recorded_inventory_matches_the_provider_record_and_the_stated_counts():
    pr = REC["provider_record"]
    assert pr["size_matches"] and pr["sha1_matches"] and pr["file_record"]["fileSizeBytes"] == 3386 == REC["csv"]["bytes"]
    assert pr["file_record"]["checksum_sha1"] == REC["csv"]["sha1"] and pr["licence"].endswith("(CC0)") and pr["submission_type"] == "PARTIAL"
    i = REC["csv"]["inspection"]
    assert i["header"] == ["Sequence", "Start", "End", "MHP", "MaxUptake", "RT"] and i["n_peptides"] == 80
    rc = i["reference_check"]
    assert (rc["matched_in_wild_type"], rc["matched_only_in_a_variant"], rc["unmatched"], rc["offsets_seen"]) == (76, {"A344V": 2, "M472I": 2}, 0, [228])
    assert i["mhp_is_monoisotopic_MH_plus"]["rows_within_0.001"] == 80 and i["maxuptake_equals_len_minus_1_minus_prolines"] == 80
    assert i["file_shape"]["n_lines"] == 82 and i["file_shape"]["line_ending"] == "CRLF"
    assert REC["csv"]["fields_present"] == {k: True for k in ("Sequence", "Start", "End", "MHP", "MaxUptake", "RT")}
    assert any("exchange time" in a for a in REC["csv"]["fields_absent"]) and any("uptake" in a for a in REC["csv"]["fields_absent"])
    assert REC["source_data_workbook"]["n_sheets"] == 8 and REC["source_data_workbook"]["sheets_with_hdx_terms"] == []
    pdf = REC["supplementary_information_pdf"]
    assert pdf["supp_fig_4B_caption_found"] and pdf["text_layer_lines_that_look_like_peptide_rows"] == 0 and pdf["raster_images_in_pdf"] == 10


def test_inventory_document_quotes_the_recorded_values():
    for needle in ("6e3fe1eea1881ebbd96dad4d704e0ab426b938cb", "2fd6b624a52fd5261c838f85a9c637f29b94e098297beb37f58c76b71705ca4a",
                   "82 lines, CRLF", "80 data rows", "ERTDITM,10,16,865.4084,6,4.493726", "VEKELGK,278,284,802.4669,6,2.585667",
                   "76 match UniProt P00519 exactly", "4,488,501", "1,468,412", "Is it sufficient for an endpoint? No"):
        assert needle in DOC, needle
    assert REC["csv"]["inspection"]["file_shape"]["last_data_row"] == "VEKELGK,278,284,802.4669,6,2.585667"
    low = DOC.lower()
    assert not [w for w in ("authentic", "pristine", "revolutionary", "obliterate", "definitive", "zero hyperparameters") if w in low]
