"""Tests of the post-run residue summary of ABL-MYR-1b (reads the results JSON only)."""

import copy
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_myr_1b"))
import summarize_results_1b as S  # noqa: E402


def test_segments():
    assert S.segments([5, 3, 4, 9, 11, 10, 3]) == [[3, 5], [9, 11]] and S.segments([7]) == [[7]] and S.segments([]) == []


@pytest.mark.skipif(not S.RES.exists(), reason="no ABL-MYR-1b results yet")
def test_summary_matches_the_results_json_and_is_reproducible():
    res = json.loads(S.RES.read_text())
    s = S.summarize(res)
    idx = {r["index"]: r["resid"] for r in res["residues"]}
    assert s["R_T"]["residues"] == [idx[i] for i in res["outcome"]["R_T"]] and s["R_T"]["count"] == len(res["outcome"]["R_T"])
    assert s["R_T_outside_N2_R_B"]["count"] == len(res["outcome"]["R_T_outside_N2_R_B"])
    assert s["outcome_label"] == res["outcome"]["label"] and "reproduced" not in json.dumps(s)
    assert json.loads(S.OUT.read_text()) == s if S.OUT.exists() else True
    bad = copy.deepcopy(res)
    bad["gates"]["G3"]["alignment"]["author_numbering_diagnostic"]["reference_index_minus_author_resid"] = {"19": 251}
    with pytest.raises(ValueError, match="not the reference numbering"):
        S.summarize(bad)
