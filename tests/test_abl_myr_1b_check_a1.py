"""SYNTHETIC tests of amendment a1 of the ABL-MYR-1b check (repeated index group names). Needs MDAnalysis."""

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_myr_1b"))
if os.environ.get("TOPOFOLD_REQUIRE_MDANALYSIS_TESTS") == "1":
    import MDAnalysis  # noqa: F401
else:
    pytest.importorskip("MDAnalysis")
import run_comp_check as B  # noqa: E402
import run_comp_check_a1 as A1  # noqa: E402
import test_abl_myr_1b_check as T  # noqa: E402


def with_defaults(text):
    """Prepend make_ndx-style default groups that repeat residue-name groups (as in the published files)."""
    return "[ ATP ]\n9\n[ MG ]\n10\n[ ATP ]\n9\n[ MG ]\n10\n" + text


def test_repeated_identical_groups_are_accepted_and_conflicting_ones_marked():
    g = A1.parse_ndx(with_defaults("[ X ]\n1 2\n"))
    assert g["ATP"].tolist() == [8] and g["X"].tolist() == [0, 1]
    g = A1.parse_ndx("[ X ]\n1\n[ X ]\n2\n")
    assert g["X"] is A1.AMBIGUOUS
    with pytest.raises(ValueError):
        A1.parse_ndx("1\n[ X ]\n2\n")
    with pytest.raises(ValueError, match="duplicate"):                               # the frozen parser still refuses
        B.parse_ndx(with_defaults(""))


def test_ambiguous_tc_group_fails_k3_and_unique_groups_follow_the_frozen_rule():
    atoms, rec = T._sys(True)
    text = T.ndx(atoms, True)
    ok = A1.tc_group_check(A1.parse_ndx(with_defaults(text)), ["Protein_MYRP_ATP_MG", "Water_NA_CL"], rec["resnames"],
                           rec["resindex"], rec["protein_resnames"], "Myr")
    assert ok["ok"], ok["problems"]
    amb = A1.tc_group_check(A1.parse_ndx(text + "[ Water_NA_CL ]\n1\n"), ["Protein_MYRP_ATP_MG", "Water_NA_CL"],
                            rec["resnames"], rec["resindex"], rec["protein_resnames"], "Myr")
    assert not amb["ok"] and "ambiguous" in amb["problems"][0]


def test_patch_is_scoped_to_the_run():
    with A1.amended():
        assert B.parse_ndx is A1.parse_ndx
    assert B.parse_ndx is not A1.parse_ndx and B.tc_group_check is not A1.tc_group_check


def _run(tmp_path, **kw):
    orig, data = T.build(tmp_path, **kw)
    res, _ = A1.run_check(orig, data, tmp_path / "ext", T.BUDGET, loader=T.loader_factory(), log=lambda *_: None)
    return res


def test_end_to_end_with_repeated_default_groups(tmp_path):
    over = {}
    for s in B.SETS:
        for c in B.CONDITIONS:
            over[(s, c)] = with_defaults(T.ndx(T.system_atoms(c == "Myr"), c == "Myr"))
    r = _run(tmp_path, ndx_override=over)
    assert r["outcome"] == "comparable_for_1b", r.get("decision", r.get("stop"))
    orig, data = T.build(tmp_path / "frozen", ndx_override=over)
    frozen, _ = B.run_check(orig, data, tmp_path / "e2", T.BUDGET, loader=T.loader_factory())
    assert frozen["stop"]["code"] == "S3"                                               # the defect, reproduced


def test_end_to_end_ambiguous_tc_group_fails_k3(tmp_path):
    over = {("Set1", "Myr"): T.ndx(T.system_atoms(True), True) + "[ Water_NA_CL ]\n1\n"}
    r = _run(tmp_path, ndx_override=over)
    assert r["outcome"] == "not_comparable" and r["decision"]["failed"] == ["K3_tc_groups"]


def test_a1_freeze_preserves_attempt1_and_frozen_check():
    frz = ROOT / "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_a1_freeze.json"
    if not frz.exists():
        pytest.skip("a1 freeze not yet written")
    f = json.loads(frz.read_text())
    for p, h in f["attempt1_sha256"].items():
        assert B.R0.sha256_file(ROOT / p) == h, p
    assert B.R0.sha256_file(B.SELF) == f["frozen_check"]["code_sha256"]
