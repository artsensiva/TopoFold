"""Consistency of the committed MET export results (the git-ignored PDB files are not needed)."""
import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

D = Path(__file__).resolve().parent.parent / "benchmarks/audit/met_feature_export_2026-10-27"
R = D / "results"
spec = importlib.util.spec_from_file_location("met_export", D / "met_export.py")
me = importlib.util.module_from_spec(spec)
sys.modules["met_export"] = me
spec.loader.exec_module(me)
man = json.loads((R / "manifest.json").read_text())


def test_committed_tables_match_the_manifest_hashes():
    for n in ("features_structure.csv", "features_mutation.csv", "coverage.csv"):
        assert hashlib.sha256((R / n).read_bytes()).hexdigest() == man["tables"][n]["sha256"], n


def test_row_universe_regenerates_to_the_manifest_hash(tmp_path):
    c = me.load_contract()
    u = me.row_universe(me.read_fasta(D / c["reference"]["file"]), c)
    me.write_csv(tmp_path / "u.csv", u, ["row_id", "drug", "uniprot_position", "dms_position", "wt", "mut"])
    assert hashlib.sha256((tmp_path / "u.csv").read_bytes()).hexdigest() == man["tables"]["row_universe.csv"]["sha256"] and len(u) == man["row_universe"]["rows"] == 287 * 19 * 6


def test_structure_table_has_one_row_per_position_and_drug_and_no_outcome_columns():
    rows = list(csv.DictReader(open(R / "features_structure.csv")))
    assert len(rows) == 287 * 6 and len({(r["drug"], r["uniprot_position"]) for r in rows}) == len(rows)
    cols = set(rows[0])
    assert not [c for c in cols if any(t in c.lower() for t in ("effect", "score", "fitness", "gamma", "rosace", "dmso"))]
    assert len(me.KAPPA_COLS) == 7 and len(me.TAU_COLS) == 12 and len(me.DIST_COLS) == 36 and set(me.KAPPA_COLS + me.TAU_COLS + me.DIST_COLS) <= cols


def test_masked_rows_have_no_values_and_a_reason_and_the_two_angle_and_distance_blocks_share_support():
    for r in csv.DictReader(open(R / "features_structure.csv")):
        assert r["T_reason"] == r["C_reason"]
        t_vals = [r[c] for c in me.KAPPA_COLS + me.TAU_COLS + me.DIST_COLS]
        if r["T_reason"] == "ok":
            assert all(v != "" for v in t_vals)
        else:
            assert all(v == "" for v in t_vals) and r["T_reason"]
        if r["L_reason"] not in ("ok", "no_side_chain"):
            assert all(r[c] == "" for c in me.L_COLS)


def test_only_the_engineered_template_residues_are_masked_as_engineered():
    eng = {k: v["engineered_detail"] for k, v in man["entries"].items() if v["engineered_detail"]}
    assert eng == {"NVP": {"1194": "PHE", "1234": "PHE"}}
    assert man["wider_engineered_template_effects_not_removed_by_masking"] is True
