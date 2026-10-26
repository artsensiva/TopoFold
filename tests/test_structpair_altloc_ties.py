"""structpair: ties in the OPTIONAL altloc-sensitivity policy (KNOWN_ISSUES m13). SYNTHETIC inputs only.

A tie met only by the optional `altloc_sensitivity: highest_occupancy` must leave the primary comparison intact and report the
sensitivity as unavailable (never zero); a tie under the primary policy and every unrelated error must still fail."""

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
from topofold_inputs import structpair as sp  # noqa: E402


def _adapter_python():
    cand = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
    if not Path(cand).is_file():
        return None
    return cand if subprocess.run([cand, "-c", "import MDAnalysis"], capture_output=True).returncode == 0 else None


PY = _adapter_python()
needs_adapter = pytest.mark.skipif(PY is None, reason="no interpreter with MDAnalysis available")
N = 14


def _helix(shift=0.0):
    t = np.arange(N) * math.radians(100)
    return np.c_[2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(N)] + shift


def _pdb(path, alts=(), shift=0.0, disp=0.4):
    """CA-only chain A of N residues. alts: {resid: [(altloc, occupancy), ...]}; each alternate is displaced by `disp` Å x index."""
    X = _helix(shift)
    alts = dict(alts)
    lines = []
    for i in range(N):
        resid = i + 1
        for k, (alt, occ) in enumerate(alts.get(resid, [(" ", 1.0)])):
            p = X[i] + disp * k
            lines.append(f"ATOM  {len(lines) + 1:5d}  CA {alt}ALA A{resid:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}{occ:6.2f} 20.00           C")
    Path(path).write_text("\n".join(lines) + "\nTER\nEND\n")


def _spec(tmp, sens="highest_occupancy", altloc_a="blank_then_A", altloc_b="blank_then_A", name="spec"):
    s = {"name": "SYNTHETIC altloc tie", "correspondence": {"same_numbering": "1-12"},
         "structures": {"A": {"label": "SYNTHETIC A", "file": "SYNTHETIC_A.pdb", "chain": "A", "altloc": altloc_a},
                        "B": {"label": "SYNTHETIC B", "file": "SYNTHETIC_B.pdb", "chain": "A", "altloc": altloc_b}}}
    if sens:
        s["altloc_sensitivity"] = sens
    (tmp / f"{name}.json").write_text(json.dumps(s))
    return tmp / f"{name}.json"


def _cli(spec, out):
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=PY)
    return subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(spec), "--out", str(out)], env=env,
                          capture_output=True, text=True)


def _summary(out):
    return json.loads((Path(out) / "summary.json").read_text())


def _same_primary(a, b):
    for f in ("windows.csv", "residues.csv", "exclusions.csv"):
        assert (Path(a) / f).read_bytes() == (Path(b) / f).read_bytes(), f
    sa, sb = _summary(a), _summary(b)
    for k in sa:
        if k != "altloc_sensitivity":
            assert sa[k] == sb[k], k


@needs_adapter
def test_no_alternate_locations_is_not_applicable(tmp_path):
    _pdb(tmp_path / "SYNTHETIC_A.pdb")
    _pdb(tmp_path / "SYNTHETIC_B.pdb", shift=0.2)
    r = _cli(_spec(tmp_path), tmp_path / "out")
    assert r.returncode == 0, r.stderr
    assert _summary(tmp_path / "out")["altloc_sensitivity"]["status"].startswith("not applicable")


@needs_adapter
def test_unique_highest_occupancy_alternate_is_computed(tmp_path):
    _pdb(tmp_path / "SYNTHETIC_A.pdb")
    _pdb(tmp_path / "SYNTHETIC_B.pdb", alts={5: [(" ", 0.4), ("B", 0.6)]})            # B (0.6) wins under highest_occupancy; blank under blank_then_A
    r = _cli(_spec(tmp_path), tmp_path / "out")
    assert r.returncode == 0, r.stderr
    s = _summary(tmp_path / "out")["altloc_sensitivity"]
    assert s["status"] == "computed" and max(v for v in s["max_abs_change"].values() if v is not None) > 0


@needs_adapter
def test_a_tie_in_the_optional_policy_keeps_the_primary_comparison_and_reports_unavailable(tmp_path):
    _pdb(tmp_path / "SYNTHETIC_A.pdb", alts={5: [(" ", 0.5), ("B", 0.5)]})
    _pdb(tmp_path / "SYNTHETIC_B.pdb", shift=0.1)
    r = _cli(_spec(tmp_path), tmp_path / "out")
    assert r.returncode == 0, r.stderr                                                # was: exit 2 and no outputs
    s = _summary(tmp_path / "out")["altloc_sensitivity"]
    assert s["status"].startswith("unavailable") and s["available"] is False and s["values_reported"] is False
    assert "max_abs_change" not in s and "global_rmsd_angstrom" not in s             # no number, in particular no zero
    assert s["ties"] == [{"structure": "A", "chain": "A", "resid": 5, "icode": "", "in_correspondence": True,
                          "altlocs": ["", "B"], "occupancies": [0.5, 0.5], "policy": "highest_occupancy"}]
    assert "partial" in s["why_not_partial"].lower() or "whole" in s["why_not_partial"].lower()
    r2 = _cli(_spec(tmp_path, sens=None, name="plain"), tmp_path / "out_plain")
    assert r2.returncode == 0, r2.stderr
    _same_primary(tmp_path / "out", tmp_path / "out_plain")                           # primary values, mapping and exclusions identical
    man = json.loads((tmp_path / "out/manifest.json").read_text())
    assert man["policies"]["altloc_sensitivity_outcome"].startswith("unavailable")
    rep = (tmp_path / "out/report.md").read_text()
    assert "UNAVAILABLE" in rep and "not zero" in rep and "A:5" in rep and "0.5" in rep


@needs_adapter
def test_ties_in_both_structures_and_several_residues_are_all_listed_including_unmapped_ones(tmp_path):
    _pdb(tmp_path / "SYNTHETIC_A.pdb", alts={3: [("A", 0.5), ("B", 0.5)], 8: [("A", 0.3), ("B", 0.3), ("C", 0.3)]})
    _pdb(tmp_path / "SYNTHETIC_B.pdb", alts={6: [("A", 0.5), ("B", 0.5)], 14: [("A", 0.5), ("B", 0.5)]})   # residue 14 is outside 1-12
    r = _cli(_spec(tmp_path, altloc_a="blank_then_A"), tmp_path / "out")
    assert r.returncode == 0, r.stderr
    ties = _summary(tmp_path / "out")["altloc_sensitivity"]["ties"]
    assert [(t["structure"], t["resid"], t["in_correspondence"]) for t in ties] == [("A", 3, True), ("A", 8, True), ("B", 6, True), ("B", 14, False)]
    assert ties[1]["altlocs"] == ["A", "B", "C"] and ties[1]["occupancies"] == [0.3, 0.3, 0.3]


@needs_adapter
def test_a_tie_in_the_primary_policy_still_fails_explicitly(tmp_path):
    _pdb(tmp_path / "SYNTHETIC_A.pdb", alts={5: [(" ", 0.5), ("B", 0.5)]})
    _pdb(tmp_path / "SYNTHETIC_B.pdb")
    r = _cli(_spec(tmp_path, sens=None, altloc_a="highest_occupancy"), tmp_path / "out")
    assert r.returncode == 2 and "tie on occupancy 0.5" in r.stderr and "('A', 5, '')" in r.stderr
    assert not (tmp_path / "out/summary.json").exists()
    r = _cli(_spec(tmp_path, sens="highest_occupancy", altloc_a="highest_occupancy", name="both"), tmp_path / "out2")
    assert r.returncode == 2 and "tie on occupancy" in r.stderr                       # the optional entry does not rescue a primary tie


def test_an_unrelated_failure_in_the_optional_calculation_is_not_swallowed(tmp_path, monkeypatch):
    if PY is None:
        pytest.skip("no interpreter with MDAnalysis available")
    _pdb(tmp_path / "SYNTHETIC_A.pdb", alts={5: [(" ", 0.4), ("B", 0.6)]})
    _pdb(tmp_path / "SYNTHETIC_B.pdb")
    orig = sp.extract

    def broken(spec, key, altloc, work, adapter_python, tag, chain=None):
        if tag == "altloc_sens":
            raise RuntimeError("extraction failed: unrelated numerical problem")
        return orig(spec, key, altloc, work, adapter_python, tag, chain)
    monkeypatch.setattr(sp, "extract", broken)
    with pytest.raises(RuntimeError, match="unrelated numerical problem"):
        sp.run(_spec(tmp_path), tmp_path / "out", PY)


@needs_adapter
def test_adapter_marks_only_tie_errors_in_the_error_file(tmp_path):
    from topofold_inputs import ca_adapter as ca
    _pdb(tmp_path / "t.pdb", alts={3: [("A", 0.5), ("B", 0.5)], 7: [("A", 0.5), ("B", 0.5)]})
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"))
    base = [PY, "-m", "topofold_inputs", str(tmp_path / "t.pdb"), "--run-id", "x", "--out", str(tmp_path), "--altloc", "highest_occupancy"]
    r = subprocess.run(base + ["--error-json", str(tmp_path / "err.json")], env=env, capture_output=True, text=True)
    assert r.returncode == 2 and r.stderr.startswith("rejected: residue ('A', 3, ''): altlocs ['A', 'B'] tie on occupancy 0.5")
    rec = json.loads((tmp_path / "err.json").read_text())
    assert rec["kind"] == "altloc_tie" and [t["resid"] for t in rec["ties"]] == [3, 7] and issubclass(ca.AltlocTieError, ca.AdapterError)
    # an unrelated refusal (policy 'reject' on the same file) and a missing file write no tie record
    (tmp_path / "err.json").unlink()
    r = subprocess.run(base[:-1] + ["reject", "--error-json", str(tmp_path / "err.json")], env=env, capture_output=True, text=True)
    assert r.returncode == 2 and not (tmp_path / "err.json").exists()
    r = subprocess.run([PY, "-m", "topofold_inputs", str(tmp_path / "missing.pdb"), "--run-id", "x", "--out", str(tmp_path), "--error-json", str(tmp_path / "err.json")],
                       env=env, capture_output=True, text=True)
    assert r.returncode != 0 and not (tmp_path / "err.json").exists()


@needs_adapter
def test_report_and_viewer_distinguish_unavailable_from_zero(tmp_path):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    res = tmp_path / "results"
    _pdb(tmp_path / "SYNTHETIC_A.pdb", alts={5: [(" ", 0.5), ("B", 0.5)]})
    _pdb(tmp_path / "SYNTHETIC_B.pdb", shift=0.1)
    assert _cli(_spec(tmp_path), res / "tie").returncode == 0
    # a computed sensitivity whose alternates coincide: a real zero change, worded differently
    z = tmp_path / "zero_in"
    z.mkdir()
    _pdb(z / "SYNTHETIC_A.pdb", alts={5: [(" ", 0.4), ("B", 0.6)]}, disp=0.0)
    _pdb(z / "SYNTHETIC_B.pdb", shift=0.1)
    assert _cli(_spec(z), res / "zero").returncode == 0
    zs = _summary(res / "zero")["altloc_sensitivity"]
    assert zs["status"] == "computed" and max(v for v in zs["max_abs_change"].values() if v is not None) < 1e-6
    assert "UNAVAILABLE" not in (res / "zero/report.md").read_text()
    assert "largest change per metric" in (res / "zero/report.md").read_text()
    at = AppTest.from_file(str(ROOT / "apps/structpair_viewer/app.py"), default_timeout=90)
    at.session_state["results_root"] = str(res)
    at.run()
    at.selectbox(key="result_choice").set_value("tie").run()
    assert not at.exception
    text = " ".join(w.value for w in at.warning)
    assert "UNAVAILABLE" in text and "not zero" in text and "A:5" in text
    at.selectbox(key="result_choice").set_value("zero").run()
    assert "UNAVAILABLE" not in " ".join(w.value for w in at.warning)
