"""Tests for the ABL preregistration input schemas and validator. SYNTHETIC documents only (no ABL data scored)."""

import copy
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
pytest.importorskip("jsonschema")
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_prereg"))
import validate_inputs as V  # noqa: E402

MAP = json.loads((ROOT / "benchmarks/audit/abl_prereg/hdx_peptide_mapping.json").read_text())


def _perm(status="granted"):
    return {"status": status, "evidence": "SYNTHETIC test value"}


def synthetic_inventory():
    systems = [{"system_id": f"{v}_{lig or 'apo'}", "variant": v, "ligands": [lig] if lig else [], "starting_structure_pdb": "1OPJ",
                "construct_first_cAbl": 229, "construct_last_cAbl": 500, "numbering": "cAbl_1a", "force_field": "a99SB-disp",
                "water_model": "a99SB-disp water", "temperature_K": 310, "ensemble": "NVT", "ionic_strength_mM": 150,
                "protonation_notes": "SYNTHETIC"} for v in ("WT", "M472I") for lig in ("imatinib", None)]
    runs = [{"run_id": f"{s['system_id']}_r{i}", "system_id": s["system_id"], "independent_start": "same structure, new velocities",
             "length_us": 10.0, "frame_interval_ns": 1.0, "n_frames": 10000, "atoms_included": "protein_heavy",
             "files": [{"path": f"SYNTHETIC_{s['system_id']}_{i}.dcd", "format": "dcd", "bytes": 1, "sha256": "0" * 64}]}
            for s in systems for i in range(3)]
    return {"schema_version": "abl-traj-inventory/1",
            "source": {"provider": "SYNTHETIC", "received_utc": "2026-10-24T00:00:00Z", "reference": "SYNTHETIC"},
            "permissions": {k: _perm() for k in ("internal_evaluation", "publication_of_derived_results",
                                                 "redistribution_of_derived_files", "later_commercial_use")},
            "systems": systems, "runs": runs}


def test_inventory_valid_and_pairs():
    r = V.validate_inventory(synthetic_inventory())
    assert r["valid"], r["errors"]
    assert {tuple(p["ligands"]) for p in r["comparable_WT_M472I_pairs"]} == {("imatinib",), ("none (apo)",)}
    assert r["apo_pair_available"] and all(p["WT_runs"] == 3 and p["M472I_runs"] == 3 for p in r["comparable_WT_M472I_pairs"])


def test_inventory_without_apo_mutant_is_flagged():
    doc = synthetic_inventory()
    doc["systems"] = [s for s in doc["systems"] if s["system_id"] != "M472I_apo"]
    doc["runs"] = [r for r in doc["runs"] if r["system_id"] != "M472I_apo"]
    r = V.validate_inventory(doc)
    assert r["valid"] and not r["apo_pair_available"]


@pytest.mark.parametrize("break_it, message", [
    (lambda d: d["runs"].append(copy.deepcopy(d["runs"][0])), "duplicate run_id"),
    (lambda d: d["runs"][0].update(n_frames=500), "inconsistent with length/interval"),
    (lambda d: d["runs"][0].update(system_id="nope"), "unknown system_id"),
    (lambda d: d["systems"][0].update(numbering="cAbl_1b"), "no offset"),
    (lambda d: d["permissions"].pop("later_commercial_use"), "later_commercial_use"),
])
def test_inventory_refusals(break_it, message):
    doc = synthetic_inventory()
    break_it(doc)
    r = V.validate_inventory(doc)
    assert not r["valid"] and any(message in e for e in r["errors"]), r["errors"]


def test_inventory_unknown_permission_blocks():
    doc = synthetic_inventory()
    doc["permissions"]["publication_of_derived_results"] = _perm("unknown")
    assert V.validate_inventory(doc)["permissions_blocking_evaluation"] == ["publication_of_derived_results"]


def synthetic_hdx():
    """SYNTHETIC uptake values on real peptide sequences (WT peptides and the M472I-specific ones)."""
    wt = [p for p in MAP["peptides"] if p["exact_match_to_P00519"]][:3]
    mut = {p["sequence"]: p for p in MAP["peptides"] if any(d["cAbl_residue"] == 472 for d in p["differences"])}
    rows = []
    for state, peps in (("WT", wt + [{"sequence": "LMRACWQWNPSDRPSF", "construct_start": 243, "construct_end": 258}]),
                        ("M472I", wt + [mut["LIRACWQWNPSDRPSF"]])):
        for p in peps:
            for t in (10.0, 60.0):
                for rep in (1, 2):
                    rows.append({"peptide_sequence": p["sequence"], "construct_start": str(p["construct_start"]),
                                 "construct_end": str(p["construct_end"]), "state": state, "ligands": "none",
                                 "labeling_time_s": str(t), "replicate": str(rep), "uptake_Da": "1.5",
                                 "charge_states_combined": "true"})
    return rows


def test_hdx_valid_including_mutant_peptide():
    r = V.validate_hdx(synthetic_hdx())
    assert r["valid"], r["errors"]
    assert r["states"] == ["M472I|none", "WT|none"] and r["labeling_times_s"] == [10.0, 60.0]


def test_hdx_refusals():
    rows = synthetic_hdx()
    one_rep = [r for r in rows if not (r["state"] == "WT" and r["replicate"] == "2" and r["labeling_time_s"] == "10.0")]
    assert any("fewer than 2 replicates" in e for e in V.validate_hdx(one_rep)["errors"])
    wrong = copy.deepcopy(rows)
    wrong[0]["peptide_sequence"] = wrong[0]["peptide_sequence"][::-1]
    assert any("does not match P00519" in e for e in V.validate_hdx(wrong)["errors"])
    gap = [r for r in rows if not (r["state"] == "M472I" and r["labeling_time_s"] == "60.0")]
    assert any("differ in peptide/time coverage" in e for e in V.validate_hdx(gap)["errors"])
    bad = copy.deepcopy(rows)
    bad[0]["state"] = "T315I"
    assert not V.validate_hdx(bad)["valid"]
