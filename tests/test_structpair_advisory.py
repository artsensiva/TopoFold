"""structpair advisory output: regression cases for the correction of 2026-10-25 (output semantics of schema /3).

The three cases of the correction request, written to run against the committed 0.2.0 code (where they fail) and the corrected
code (where they pass). They read the advisory columns by either name: `separated` (schema /2) or `advisory_flag` (schema /3).
"""

import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs import chain_variation as cv  # noqa: E402
from topofold_inputs import structpair as sp  # noqa: E402


def flag(row):
    """The advisory flag of a row under either schema: True / False / None (not assessed)."""
    v = row["advisory_flag"] if "advisory_flag" in row else row["separated"]
    return {"True": True, "False": False, "": None}.get(v, v) if isinstance(v, str) else v


def boundary(row):
    v = row["advisory_at_boundary"] if "advisory_at_boundary" in row else row["at_boundary"]
    return {"True": True, "False": False, "": None}.get(v, v) if isinstance(v, str) else v


def _rows(value_fn, params=None, ca=("A", "B"), cb=("A", "B"), n=3):
    keys = cv.pair_keys(list(ca), list(cb))
    pv = {k: {ch: [value_fn(k, ch, i) for i in range(n)] for ch in cv.CHANNELS} for g in keys.values() for k in g}
    return cv.variation_rows(pv, n, [(f"A:{i}", f"B:{i}") for i in range(n)], [("w3", "w4")] * n, list(ca), list(cb),
                             ("A", "A"), cv.parse_params(params))


# --- case 1: identical chains, every value zero -------------------------------------------------------------------------------

def test_identical_chains_with_all_values_zero_are_never_flagged_by_default():
    rows = _rows(lambda k, ch, i: 0.0)
    assert all(r["between_max"] == 0.0 and r["within_max"] == 0.0 for r in rows)
    assert [r for r in rows if flag(r) is True] == []                    # never "separated"
    assert all(flag(r) is None for r in rows)                             # default: not assessed, with a reason
    assert all(r.get("advisory_not_assessed_reason", r.get("not_assessed_reason")) for r in rows)


def test_identical_chains_with_a_configured_floor_are_flagged_false_not_true():
    rows = _rows(lambda k, ch, i: 0.0, {"advisory": {"min_difference": {"disp": 0.5}}})
    disp = [r for r in rows if r["channel"] == "disp"]
    assert disp and all(flag(r) is False for r in disp)


# --- case 2: differences indistinguishable from zero --------------------------------------------------------------------------

@pytest.mark.parametrize("tiny", [1e-12, 5e-10, 9.99e-10])
def test_differences_indistinguishable_from_zero_do_not_become_a_finding(tiny):
    cross = lambda k, ch, i: tiny if k[0] == "cross" else 0.0   # noqa: E731
    rows = _rows(cross)
    assert [r for r in rows if flag(r) is True] == []
    rows = _rows(cross, {"advisory": {"min_difference": {c: 1e-6 for c in cv.CHANNELS}}})
    assert [r for r in rows if flag(r) is True] == []


def test_default_reports_the_numbers_but_leaves_the_flag_unassessed():
    big = lambda k, ch, i: 5.0 if k[0] == "cross" else 1.0   # noqa: E731
    rows = _rows(big)
    r = [x for x in rows if x["channel"] == "disp"][0]
    assert (r["between_min"], r["between_max"], r["within_max"]) == (5.0, 5.0, 1.0)      # numbers are always reported
    assert flag(r) is None                                                                # no configured minimum difference: unassessed
    assert "minimum difference" in (r.get("advisory_not_assessed_reason") or r.get("not_assessed_reason") or "")
    on = _rows(big, {"advisory": {"min_difference": {"disp": 2.0}}})
    d = [x for x in on if x["channel"] == "disp"][0]
    assert flag(d) is True                                                                # configured: 5.0 >= max(2.0, 1.0 x 1.0)
    assert flag([x for x in on if x["channel"] == "rmsd3"][0]) is None                    # other channels stay unassessed


def test_equality_at_a_positive_boundary_counts_and_a_hairline_below_does_not():
    p = {"advisory": {"statistic": "min", "within_factor": 1.0, "min_difference": {"disp": 1.0}}}
    at = _rows(lambda k, ch, i: 1.0 if k[0] == "cross" else 0.5, p)
    r = [x for x in at if x["channel"] == "disp"][0]
    assert flag(r) is True and boundary(r) is True                                        # equality meets the rule, marked at the boundary
    below = _rows(lambda k, ch, i: (1.0 - 1e-12) if k[0] == "cross" else 0.5, p)
    r = [x for x in below if x["channel"] == "disp"][0]
    assert flag(r) is False and boundary(r) is True                                       # no tolerance widening; the boundary is marked
    above = _rows(lambda k, ch, i: 1.1 if k[0] == "cross" else 0.5, p)
    assert flag([x for x in above if x["channel"] == "disp"][0]) is True
    factor = _rows(lambda k, ch, i: 1.0 if k[0] == "cross" else 0.8, dict(p, advisory=dict(p["advisory"], within_factor=2.0)))
    r = [x for x in factor if x["channel"] == "disp"][0]
    assert flag(r) is False                                                               # threshold max(1.0, 2 x 0.8) = 1.6 > 1.0


# --- case 3: non-finite numbers are rejected at spec parsing -----------------------------------------------------------------

NAN, INF = float("nan"), float("inf")


@pytest.mark.parametrize("block", [
    {"displacement_trim_fraction": NAN}, {"displacement_trim_fraction": INF}, {"displacement_trim_fraction": -INF},
    {"advisory": {"within_factor": NAN}}, {"advisory": {"within_factor": INF}}, {"advisory": {"within_factor": -INF}},
    {"advisory": {"min_difference": {"disp": NAN}}}, {"advisory": {"min_difference": {"disp": INF}}},
    {"advisory": {"min_difference": {"d_kappa": -INF}}}])
def test_nan_and_infinity_are_rejected_at_parsing(block):
    with pytest.raises(cv.VariationError):
        cv.parse_params(block)


def test_non_finite_values_are_rejected_when_a_spec_file_is_loaded(tmp_path):
    for name in ("a.pdb", "b.pdb"):
        (tmp_path / name).write_text("END\n")
    for bad in ("NaN", "Infinity", "-Infinity"):
        text = ('{"structures": {"A": {"file": "a.pdb", "chain": "A"}, "B": {"file": "b.pdb", "chain": "A"}}, '
                '"correspondence": {"same_numbering": "1-6"}, "chain_variation": {"advisory": {"within_factor": %s}}}' % bad)
        (tmp_path / "spec.json").write_text(text)
        with pytest.raises(sp.SpecError, match="within_factor"):
            sp.load_spec(tmp_path / "spec.json")
        text = text.replace('"within_factor": %s' % bad, '"min_difference": {"disp": %s}' % bad)
        (tmp_path / "spec.json").write_text(text)
        with pytest.raises(sp.SpecError, match="min_difference"):
            sp.load_spec(tmp_path / "spec.json")


def test_a_minimum_difference_must_be_positive_and_above_the_numerical_tolerance():
    for bad in (0.0, 1e-9, 5e-10, -1.0):
        with pytest.raises(cv.VariationError, match="min_difference"):
            cv.parse_params({"advisory": {"min_difference": {"disp": bad}}})
    ok = cv.parse_params({"advisory": {"min_difference": {"disp": 1e-6}}})
    assert ok["advisory"]["min_difference"] == {"disp": 1e-6}


# --- end to end: identical chains through the command -------------------------------------------------------------------------

def _adapter_python():
    cand = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
    if not Path(cand).is_file():
        return None
    return cand if subprocess.run([cand, "-c", "import MDAnalysis"], capture_output=True).returncode == 0 else None


def _pdb_two_chains(path, coords):
    lines = []
    for ch in "AB":
        for i, (x, y, z) in enumerate(coords):
            lines.append(f"ATOM  {len(lines) + 1:5d}  CA  ALA {ch}{i + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C")
        lines.append("TER")
    path.write_text("\n".join(lines) + "\nEND\n")


def _helix(n):
    t = np.arange(n) * math.radians(100)
    return np.c_[2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)]


def _run(tmp_path, name, extra):
    spec = {"name": "SYNTHETIC identical chains", "correspondence": {"same_numbering": "1-12"},
            "structures": {"A": {"label": "SYNTHETIC A", "file": "SYNTHETIC_A.pdb", "chains": ["A", "B"], "primary_chain": "A"},
                           "B": {"label": "SYNTHETIC B", "file": "SYNTHETIC_B.pdb", "chains": ["A", "B"], "primary_chain": "A"}}}
    spec.update(extra)
    (tmp_path / f"{name}.json").write_text(json.dumps(spec))
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=_adapter_python())
    return subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(tmp_path / f"{name}.json"), "--out", str(tmp_path / name)],
                          env=env, capture_output=True, text=True)


def test_end_to_end_identical_chains_in_both_entries_are_never_described_as_separated(tmp_path):
    if _adapter_python() is None:
        pytest.skip("no interpreter with MDAnalysis available")
    X = _helix(12)
    _pdb_two_chains(tmp_path / "SYNTHETIC_A.pdb", X)
    _pdb_two_chains(tmp_path / "SYNTHETIC_B.pdb", X)                         # the same coordinates in all four chains
    r = _run(tmp_path, "default", {})
    assert r.returncode == 0, r.stderr
    rows = list(csv.DictReader(open(tmp_path / "default/chain_variation.csv")))
    assert max(float(x["between_max"]) for x in rows if x["between_max"]) < 1e-6        # all values are zero
    assert max(float(x["within_max"]) for x in rows if x["within_max"]) < 1e-6
    assert [x for x in rows if flag(x) is True] == []
    report = (tmp_path / "default/report.md").read_text().lower()
    assert "separated" not in report                                                     # the word is not used for a zero comparison
    r = _run(tmp_path, "floor", {"chain_variation": {"advisory": {"min_difference": {"disp": 0.5}}}})
    assert r.returncode == 0, r.stderr
    rows = [x for x in csv.DictReader(open(tmp_path / "floor/chain_variation.csv")) if x["channel"] == "disp"]
    assert rows and all(flag(x) is False for x in rows)
