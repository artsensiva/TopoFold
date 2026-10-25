"""Oncology search (2026-10-25): helper functions on SYNTHETIC input; the decision record against the generated source record."""

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "benchmarks/audit/oncology_search_2026-10-25"
sys.path.insert(0, str(D))
import collect_sources as S  # noqa: E402

REC = json.loads((D / "results/sources_verified.json").read_text())
DOC = (ROOT / "docs/ONCOLOGY_EVALUATION_SEARCH_2026-10-25.md").read_text()


def test_plain_text_strips_tags_and_collapses_space():
    assert S.plain_text("<p>Data  <b>availability</b>\n x &amp; y</p>") == "Data availability x & y"


def test_excerpt_starts_before_the_match_and_is_capped():
    text = "a" * 300 + " TARGET " + "b" * 2000
    e = S.excerpt(text, r"TARGET", 50, 1500)
    assert "TARGET" in e and len(e) <= S.EXCERPT_MAX and S.excerpt("none here", r"TARGET") is None


def test_hdx_file_summary_counts_raw_files_per_protein_and_time():
    names = [("1-x_ABL_WT_10minD2O.raw.zip", 10), ("1-y_ABL_WT_10minD2Orep.raw.zip", 20), ("2-x_ABL_M472I_4hD2O.raw.zip", 5),
             ("3-x_ABL_M472I_ctrl2.raw.zip", 7), ("table.csv", 3)]
    s = S.hdx_file_summary(names)
    assert s["raw_files_per_protein"] == {"M472I": 2, "WT": 2} and s["total_bytes"] == 45
    assert s["labels_per_protein"]["M472I"] == ["4h", "ctrl2"]


def test_record_has_every_source_and_no_submitter_email():
    assert set(REC["crossref"]) == {"ayaz_2023", "pantsar_2018", "hunter_2015", "rennella_2025", "parker_2022"}
    assert all(v is not None for v in REC["statements"].values())
    assert "@" not in json.dumps(REC["pride"])
    assert REC["hunter_2015_europepmc"]["pmcid"] is None and REC["hunter_2015_europepmc"]["isOpenAccess"] == "N"
    assert REC["crossref"]["hunter_2015"]["licence_urls"] == []
    assert any("by/4.0" in u for u in REC["crossref"]["ayaz_2023"]["licence_urls"])
    assert any("by-nc-nd" in u for u in REC["crossref"]["parker_2022"]["licence_urls"])


def test_decision_record_numbers_match_the_generated_record():
    p = REC["pride"]["PXD034008"]
    assert p["files_by_category"] == {"RAW": 43, "SEARCH": 1}
    assert p["raw_files_per_protein"] == {"A344V": 13, "M472I": 14, "WT": 12}
    assert f"{p['total_bytes'] / 1e9:.2f} GB" in DOC and "43 raw files" in DOC and "WT 12, A344V 13, M472I 14" in DOC
    assert dict(p["non_raw_files"])["REI_ABL_KB_dynamX_export_18May22.csv"] == 3386 and "3,386 B" in DOC
    assert "non-commercial use" in REC["statements"]["ayaz_data_availability"]
    assert "PXD034008" in REC["statements"]["ayaz_data_availability"]
    assert "336,034,779 B" in DOC


def test_decision_record_states_the_required_conclusions():
    must = ["No candidate passes a defensible feasibility gate today", "Status correction", "not \"two negative method results\"",
            "Structural agreement is not functional or clinical evidence", "Retained options", "Mpro M1", "IDP", "RNA",
            "Rascore", "replicate-aware structural comparison", "Why now rather than Mpro M1"]
    assert not [m for m in must if m not in DOC]
    low = DOC.lower()
    assert not [w for w in ("authentic", "pristine", "revolutionary", "obliterate", "definitive", "zero hyperparameters") if re.search(w, low)]


def test_logs_and_addendum_carry_the_correction():
    log = (ROOT / "docs/DECISION_LOG.md").read_text()
    assert "| D71 |" in log and "| D72 |" in log
    add = (ROOT / "benchmarks/audit/abl_myr_1b/ADDENDUM_2026-10-25_closeout_lineage.md").read_text()
    assert "Status wording" in add and "not \"two negative method results\"" in add
    road = (ROOT / "docs/ENGINE_MISSION_AND_ROADMAP.md").read_text()
    assert "Oncology evaluation search (D71, D72)" in road and "no candidate passes a feasibility gate today" in road
