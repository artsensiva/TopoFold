"""SYNTHETIC test of the ABL-MYR-1 G1 member-audit diagnostic (benchmarks/audit/abl_myr_1/g1_member_audit.py)."""

import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_myr_1"))
import g1_member_audit as A  # noqa: E402
import run_abl_myr as R  # noqa: E402

CFG = {**R.CFG, "sets": {"Set1": [11], "Set2": [21]}, "conditions": {"Myr": "kinase_kinked_Myr"}}


def test_audit_reports_exactly_the_missing_member_from_listings_only(tmp_path):
    members = {"PCA.zip": ["PCA/KD_noIa_bb_ref.pdb", "PCA/average.pdb", "PCA/steps", "PCA/traj_kinked_Myr_11.xtc",
                           "PCA/traj_kinked_Myr_21.xtc"]}
    model = [f"kinase_kinked_Myr/{s}/{f}" for s in ("Set1", "Set2") for f in
             ("system_em.gro", "topol.top", "production_md.mdp", "nvt.mdp", "npt.mdp")]
    model += ["kinase_kinked_Myr/Set1/system.gro", "kinase_kinked_Myr/Set1/11/md.tpr", "kinase_kinked_Myr/Set2/21/md.tpr"]
    members["kinase_kinked_Myr.zip"] = model
    for name, items in members.items():
        with zipfile.ZipFile(tmp_path / name, "w") as zf:
            for n in items:
                zf.writestr(n, b"x")
    out = A.audit_archives(tmp_path, CFG)
    assert out["PCA.zip"]["missing"] == [] and out["PCA.zip"]["present"] == out["PCA.zip"]["required"] == 5
    m = out["kinase_kinked_Myr.zip"]
    assert [x["suffix"] for x in m["missing"]] == ["Set2/system.gro"] and m["present"] == m["required"] - 1


def test_readme_tree_parser_finds_set_level_gro_files():
    tree = "\n".join(["+---kinase_kinked_Myr", "|   \\---kinase_kinked_Myr", "|       +---Set1", "|       |   |   system.gro",
                      "|       |   |   system_em.gro", "|       |   \\---11", "|       |           md.tpr", "|       \\---Set2",
                      "|           |   system_em.gro", "|           |   topol.top", "+---kinase_kinked_noMyr",
                      "|   \\---kinase_kinked_noMyr", "|       \\---Set2", "|           |   system_em.gro"])
    got = A.readme_set_gro_files(tree)
    assert got["kinase_kinked_Myr"] == {"Set1": ["system.gro", "system_em.gro"], "Set2": ["system_em.gro"]}
    assert got["kinase_kinked_noMyr"] == {"Set2": ["system_em.gro"]}
