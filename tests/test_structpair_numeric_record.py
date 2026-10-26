"""structpair numeric_record.json (KNOWN_ISSUES m13b): computed values before six-figure display rounding. SYNTHETIC inputs."""

import csv
import json
import shutil
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tests"))
import test_structpair_altloc_ties as T  # noqa: E402  (SYNTHETIC helpers: _pdb, _spec, _cli, _summary, PY)
from topofold_inputs import numeric_record as nr  # noqa: E402
from topofold_inputs import structpair as sp  # noqa: E402
from topofold_inputs import structpair_verify as V  # noqa: E402

pytestmark = T.needs_adapter


def _inputs(tmp, alts_a=None):
    T._pdb(tmp / "SYNTHETIC_A.pdb", alts=alts_a or {})
    T._pdb(tmp / "SYNTHETIC_B.pdb", shift=0.137)


def _spec_with(tmp, name, numeric, sens=None, pairs=None):
    p = T._spec(tmp, sens=sens, name=name)
    s = json.loads(p.read_text())
    if numeric is not None:
        s["numeric_record"] = numeric
    if pairs is not None:
        s["correspondence"] = {"pairs": pairs}
    p.write_text(json.dumps(s))
    return p


def _run(tmp, name, **kw):
    r = T._cli(_spec_with(tmp, name, **kw), tmp / name.replace("spec", "out"))
    assert r.returncode == 0, r.stderr
    return tmp / name.replace("spec", "out")


def test_enabling_the_export_leaves_every_existing_output_unchanged(tmp_path):
    _inputs(tmp_path)
    off = _run(tmp_path, "spec_off", numeric=None)
    on = _run(tmp_path, "spec_on", numeric=True)
    for f in ("windows.csv", "residues.csv", "exclusions.csv"):
        assert (off / f).read_bytes() == (on / f).read_bytes(), f
    assert json.loads((off / "summary.json").read_text()) == json.loads((on / "summary.json").read_text())
    mo, mn = (json.loads((d / "manifest.json").read_text()) for d in (off, on))
    assert set(mn["outputs_sha256"]) - set(mo["outputs_sha256"]) == {"numeric_record.json"}
    assert all(mn["outputs_sha256"][k] == v for k, v in mo["outputs_sha256"].items() if k != "report.md")
    assert not (off / "numeric_record.json").exists() and mn["policies"]["numeric_record"] is True and mo["policies"]["numeric_record"] is False
    assert V.verify(on)["all_ok"] and V.verify(off)["all_ok"]


def test_record_equals_the_preformat_computed_values_and_rounds_to_the_csv(tmp_path):
    _inputs(tmp_path)
    out = _run(tmp_path, "spec_on", numeric=True)
    rec = nr.read(out)
    spec = sp.load_spec(tmp_path / "spec_on.json")
    work = tmp_path / "w"
    work.mkdir()
    c = sp.compare(spec, {"A": "blank_then_A", "B": "blank_then_A"}, work, T.PY, "check")
    n_vals = 0
    for k, vals in c["metrics"].items():
        assert rec["metrics"][k]["values"] == [None if v is None else float(v) for v in vals], k      # exact float64 equality, not a tolerance
        n_vals += sum(v is not None for v in vals)
    assert n_vals > 50
    g = c["global"]
    assert rec["global"]["rmsd_qcp_angstrom"] == g["rmsd_qcp_angstrom"] and rec["global"]["n_ca"] == g["n_ca"]
    rows = list(csv.DictReader(open(out / "windows.csv")))
    for k, unit in sp.COLUMNS:
        for r, v in zip(rows, rec["metrics"][k]["values"]):
            assert (r[f"{k}_{unit}"] == "") == (v is None)
            if v is not None:
                assert r[f"{k}_{unit}"] == f"{v:.6g}"                       # the CSV is the six-figure rendering of the record
    assert rec["metrics"]["d_kappa"]["unit"] == "rad" and rec["metrics"]["d_i_i2"]["unit"] == "angstrom"
    assert "not higher experimental accuracy" in rec["statement"] and "float32" in rec["provenance"]["coordinate_precision"]


def test_missing_values_are_null_with_a_reason_never_zero(tmp_path):
    _inputs(tmp_path)
    pairs = [[i, i] for i in range(1, 13) if i != 6]                       # residue 6 absent from the correspondence: a numbering gap
    out = _run(tmp_path, "spec_gap", numeric=True, pairs=pairs)
    rec = nr.read(out)
    d = rec["metrics"]["d_kappa"]
    assert d["values"][0] is None and d["status"][0].startswith("not defined")            # chain end
    assert d["values"][-1] is None and d["status"][-1].startswith("not defined")
    gap = [i for i, s in enumerate(d["status"]) if s.startswith("excluded")]
    assert gap and all(d["values"][i] is None for i in gap) and "numbering_discontinuity" in d["status"][gap[0]]
    assert not any(v == 0 for m in rec["metrics"].values() for v, s in zip(m["values"], m["status"]) if s != "value")
    ex = list(csv.DictReader(open(out / "exclusions.csv")))
    assert {i for i, s in enumerate(d["status"]) if s.startswith("excluded")} == {int(r["position"]) for r in ex if r["window_size"] == "3"}


def test_identifiers_align_with_the_residue_table_and_a_misaligned_record_is_refused(tmp_path):
    _inputs(tmp_path)
    out = _run(tmp_path, "spec_on", numeric=True)
    rec = nr.read(out)
    rows = list(csv.DictReader(open(out / "residues.csv")))
    assert [(p["position"], p["A_residue"], p["B_residue"]) for p in rec["positions"]] == [(int(r["position"]), r["A_residue"], r["B_residue"]) for r in rows]
    # a record whose identifiers or lengths are wrong is refused even if its hash is made to match
    for mutate in (lambda r: r["positions"][3].__setitem__("A_residue", "A:99"), lambda r: r["metrics"]["d_kappa"]["values"].pop(),
                   lambda r: r["metrics"]["d_kappa"]["status"].__setitem__(5, "value") if r["metrics"]["d_kappa"]["values"][5] is None else r["metrics"]["d_kappa"]["values"].__setitem__(5, None)):
        bad = tmp_path / "bad"
        shutil.rmtree(bad, ignore_errors=True)
        shutil.copytree(out, bad)
        r = json.loads((bad / "numeric_record.json").read_text())
        mutate(r)
        (bad / "numeric_record.json").write_text(json.dumps(r))
        m = json.loads((bad / "manifest.json").read_text())
        m["outputs_sha256"]["numeric_record.json"] = sp._sha256(bad / "numeric_record.json")
        (bad / "manifest.json").write_text(json.dumps(m))
        with pytest.raises(nr.NumericRecordError):
            nr.read(bad)


def test_declared_but_missing_modified_or_unlisted_record_fails_verification(tmp_path):
    _inputs(tmp_path)
    out = _run(tmp_path, "spec_on", numeric=True)
    assert V.verify(out)["all_ok"] and nr.read(out)
    mod = tmp_path / "mod"
    shutil.copytree(out, mod)
    f = mod / "numeric_record.json"
    f.write_text(f.read_text().replace("1.", "2.", 1))
    r = V.verify(mod)
    assert not r["all_ok"] and [c["path"] for c in r["checks"] if c["status"] == "mismatch"] == ["numeric_record.json"]
    with pytest.raises(nr.NumericRecordError, match="differs from its manifest hash"):
        nr.read(mod)
    gone = tmp_path / "gone"
    shutil.copytree(out, gone)
    (gone / "numeric_record.json").unlink()
    r = V.verify(gone)
    assert not r["all_ok"] and r["n_missing"] == 1
    with pytest.raises(nr.NumericRecordError, match="missing"):
        nr.read(gone)
    unl = tmp_path / "unlisted"
    shutil.copytree(out, unl)
    m = json.loads((unl / "manifest.json").read_text())
    del m["outputs_sha256"]["numeric_record.json"]
    (unl / "manifest.json").write_text(json.dumps(m))
    r = V.verify(unl)
    assert not r["all_ok"] and any("declared by policies.numeric_record" in str(c["expected"]) for c in r["checks"])
    with pytest.raises(nr.NumericRecordError, match="lists no hash"):
        nr.read(unl)


def test_a_result_written_without_the_export_has_no_record_and_the_option_must_be_boolean(tmp_path):
    _inputs(tmp_path)
    out = _run(tmp_path, "spec_off", numeric=None)
    with pytest.raises(nr.NumericRecordError, match="without a numeric record"):
        nr.read(out)
    p = _spec_with(tmp_path, "spec_bad", numeric="yes")
    with pytest.raises(sp.SpecError, match="numeric_record must be true or false"):
        sp.load_spec(p)


def test_unavailable_sensitivity_stays_distinct_from_a_computed_zero(tmp_path):
    ta = tmp_path / "tie"
    ta.mkdir()
    T._pdb(ta / "SYNTHETIC_A.pdb", alts={5: [(" ", 0.5), ("B", 0.5)]})
    T._pdb(ta / "SYNTHETIC_B.pdb", shift=0.1)
    tie = _run(ta, "spec_on", numeric=True, sens="highest_occupancy")
    s = nr.read(tie)["altloc_sensitivity"]
    assert s["available"] is False and "max_abs_change" not in s and s["status"].startswith("unavailable")
    za = tmp_path / "zero"
    za.mkdir()
    T._pdb(za / "SYNTHETIC_A.pdb", alts={5: [(" ", 0.4), ("B", 0.6)]}, disp=0.0)       # alternates at the same position: a real zero
    T._pdb(za / "SYNTHETIC_B.pdb", shift=0.1)
    zero = _run(za, "spec_on", numeric=True, sens="highest_occupancy")
    z = nr.read(zero)["altloc_sensitivity"]
    assert z["status"] == "computed" and all(v == 0.0 or v is None for v in z["max_abs_change"].values()) and any(v == 0.0 for v in z["max_abs_change"].values())
