"""Tests of the ABL-MYR-1b supplement generator: SYNTHETIC helper tests, and consistency of the committed record."""

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / "benchmarks/audit/abl_myr_1b"
sys.path.insert(0, str(STUDY))
import plot_residue_coverage_supplement as P  # noqa: E402


def test_segments_merge_adjacent_and_keep_gaps():
    assert P.segments([5, 3, 4, 9, 11, 10, 3]) == [(3, 5), (9, 11)]
    assert P.segments([]) == []


def test_neighbourhood_is_plus_minus_two():
    assert P.neighbourhood({10}) == {8, 9, 10, 11, 12}
    assert P.neighbourhood({10, 13}) == set(range(8, 16))
    assert P.neighbourhood({10}, k=0) == {10}


def _synthetic(r_t, r_b, span=(254, 258)):
    """SYNTHETIC results/summary pair on a 5-residue topology (indices 0-4 = resid 254-258)."""
    idx = lambda xs: [x - span[0] for x in xs]  # noqa: E731
    n2b, n2t = P.neighbourhood(r_b), P.neighbourhood(r_t)
    ot, ob = sorted(set(r_t) - n2b), sorted(set(r_b) - n2t)
    res = {"residues": [{"index": i, "resid": span[0] + i} for i in range(span[1] - span[0] + 1)],
           "outcome": {"R_T": idx(r_t), "R_B": idx(r_b), "R_phi_psi": [], "R_T_outside_N2_R_B": idx(ot),
                       "R_B_outside_N2_R_T": idx(ob), "label": "SYNTHETIC", "decision_rule_step": 4}}
    summ = {"R_T": {"residues": r_t}, "R_B": {"residues": r_b}, "R_phi_psi": {"residues": []},
            "R_T_outside_N2_R_B": {"residues": ot}, "R_B_outside_N2_R_T": {"residues": ob},
            "alpha_I": {"represented": False, "not_represented": list(range(518, 532))},
            "topology_residue_range": list(span)}
    return res, summ


def test_load_check_accepts_consistent_synthetic_record(monkeypatch):
    monkeypatch.setattr(P, "SPAN", (254, 258))
    res, summ = _synthetic([255, 256], [256, 258])
    d = P.load_and_check(res, summ)
    assert d["R_T"] == [255, 256] and d["T_outside"] == [] and d["B_outside"] == []


def test_load_check_rejects_a_summary_that_disagrees_with_results(monkeypatch):
    monkeypatch.setattr(P, "SPAN", (254, 258))
    res, summ = _synthetic([255, 256], [256, 258])
    summ["R_T"]["residues"] = [255, 257]
    with pytest.raises(ValueError, match="disagree"):
        P.load_and_check(res, summ)


def test_load_check_rejects_a_wrong_outside_set(monkeypatch):
    monkeypatch.setattr(P, "SPAN", (254, 258))
    res, summ = _synthetic([255, 256], [256, 258])
    summ["R_B_outside_N2_R_T"]["residues"] = [258]
    with pytest.raises(ValueError, match="outside sets"):
        P.load_and_check(res, summ)


def test_committed_record_matches_the_recorded_numbers_and_files():
    res, summ = json.loads(P.RES.read_text()), json.loads(P.SUMMARY.read_text())
    d = P.load_and_check(res, summ)
    assert (len(d["R_T"]), len(d["R_B"]), len(d["R_phi_psi"])) == (34, 40, 15)
    assert d["T_outside"] == [] and d["B_outside"] == [317]
    prov = json.loads((P.FIG_DIR / f"{P.STEM}_provenance.json").read_text())
    for rel, h in prov["inputs_sha256"].items():
        assert P.sha256_file(STUDY / rel) == h
    assert P.sha256_file(P.FIG_DIR / f"{P.STEM}.png") == prov["output"]["sha256"]
    assert P.sha256_file(Path(P.__file__).resolve()) == prov["generator"]["sha256"]
    assert prov["output"]["bytes"] < 5 * 1024 * 1024
