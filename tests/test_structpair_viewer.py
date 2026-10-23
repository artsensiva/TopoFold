"""Tests for the TopoFold Structural Comparison viewer (apps/structpair_viewer). Read-only behaviour, displayed values,
residue highlighting, missing residues, insertion codes and invalid inputs. SYNTHETIC fixture plus the committed KRAS
structpair demonstration outputs (read only)."""

import csv
import hashlib
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
APP_DIR = ROOT / "apps/structpair_viewer"
sys.path.insert(0, str(APP_DIR))
import results as R  # noqa: E402

DEMO = ROOT / "benchmarks/audit/structpair_demo/results"
STUDIES = [DEMO / "kras1", DEMO / "kras2"]
streamlit = pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest  # noqa: E402


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# --- SYNTHETIC fixture: 8 positions, insertion code 52A in A, residue missing in B at position 4 -------------------------

A_RES = ["50", "51", "52", "52A", "53", "54", "55", "56"]
B_RES = ["150", "151", "152", "153", "154", "155", "156", "157"]
B_MISSING = 4


def _pdb(path, chain, labels, skip=()):
    lines = []
    for i, lab in enumerate(labels):
        if i in skip:
            continue
        num, ic = (lab[:-1], lab[-1]) if lab[-1].isalpha() else (lab, " ")
        x, y, z = 3.8 * i, 0.6 * (-1) ** i, 0.3 * i
        lines.append(f"ATOM  {len(lines) + 1:5d}  CA  ALA {chain}{int(num):4d}{ic}   {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C")
    path.write_text("\n".join(lines) + "\nEND\n")


def make_synthetic(tmp: Path) -> Path:
    out = tmp / "SYNTHETIC_result"
    out.mkdir()
    _pdb(tmp / "SYNTHETIC_A.pdb", "A", A_RES)
    _pdb(tmp / "SYNTHETIC_B.pdb", "B", B_RES, skip=(B_MISSING,))
    n = len(A_RES)
    res_rows, win_rows, exc_rows = [], [], []
    for i in range(n):
        res_rows.append({"position": i, "A_residue": f"A:{A_RES[i]}", "A_resname": "ALA", "A_status": "present",
                         "A_altloc_used": "", "B_residue": f"B:{B_RES[i]}", "B_resname": "ALA" if i != B_MISSING else "",
                         "B_status": "residue_not_found" if i == B_MISSING else "present", "B_altloc_used": "",
                         "resname_mismatch": "False"})
    for i in range(n):
        w3 = list(range(i - 1, i + 2)) if 1 <= i <= n - 2 else []
        w4 = list(range(i - 1, i + 3)) if 1 <= i <= n - 3 else []
        ok3, ok4 = w3 and B_MISSING not in w3, w4 and B_MISSING not in w4
        row = {"position": i, "A_residue": f"A:{A_RES[i]}", "B_residue": f"B:{B_RES[i]}", "resname_A": "ALA",
               "resname_B": res_rows[i]["B_resname"],
               "window3_A": f"A:{A_RES[w3[0]]}..A:{A_RES[w3[-1]]}" if w3 else "",
               "window4_A": f"A:{A_RES[w4[0]]}..A:{A_RES[w4[-1]]}" if w4 else "",
               "altloc_in_window3": "False", "altloc_in_window4": "False",
               "mismatch_in_window3": "False", "mismatch_in_window4": "False"}
        for k, unit in (("kappa_A", "rad"), ("kappa_B", "rad"), ("d2_A", "angstrom"), ("d2_B", "angstrom"),
                        ("d_kappa", "rad"), ("d_i_i2", "angstrom"), ("rmsd3", "angstrom")):
            row[f"{k}_{unit}"] = f"{0.1 * i + 0.01:.6g}" if ok3 else ""
        for k, unit in (("tau_A", "rad"), ("tau_B", "rad"), ("d3_A", "angstrom"), ("d3_B", "angstrom"),
                        ("d_tau", "rad"), ("d_i_i3", "angstrom"), ("rmsd4", "angstrom")):
            row[f"{k}_{unit}"] = f"{0.2 * i + 0.02:.6g}" if ok4 else ""
        win_rows.append(row)
        for size, ok, w in ((3, ok3, w3), (4, ok4, w4)):
            if w and not ok:
                exc_rows.append({"window_size": size, "position": i, "structure": "B",
                                 "window_A": f"A:{A_RES[w[0]]}..A:{A_RES[w[-1]]}", "cause": "missing_coordinates",
                                 "location": f"B:{B_RES[B_MISSING]}", "source": "adapter residue status",
                                 "detail": "residue status residue_not_found"})
    for name, rows in (("residues.csv", res_rows), ("windows.csv", win_rows), ("exclusions.csv", exc_rows)):
        with open(out / name, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    (out / "summary.json").write_text(json.dumps({"global": {"rmsd_qcp_angstrom": 0.5, "n_ca": 7,
                                                             "rmsd_numpy_kabsch_angstrom": 0.5}}))
    manifest = {"schema": R.SCHEMA, "name": "SYNTHETIC viewer fixture", "statement": "SYNTHETIC test data.",
                "command": "test", "spec": {"path": "n/a", "sha256": "0" * 64}, "created_utc": "n/a",
                "inputs": {k: {"label": f"SYNTHETIC {k}", "path": str(tmp / f"SYNTHETIC_{k}.pdb"),
                               "sha256": _sha(tmp / f"SYNTHETIC_{k}.pdb"), "chain": k, "altloc_policy": "reject"}
                           for k in "AB"},
                "correspondence": {"form": "pairs", "n_positions": n}, "breaks": {"A": [], "B": []},
                "policies": {"on_missing": "exclude", "resname_mismatch": "reject", "altloc_sensitivity": None},
                "software": {"structpair_version": "test", "module_sha256": "0" * 64, "python": "", "numpy": "",
                             "qcp": {"function": "n/a", "MDAnalysis": "n/a"}}}
    manifest["outputs_sha256"] = {f: _sha(out / f) for f in ("windows.csv", "residues.csv", "exclusions.csv",
                                                              "summary.json")}
    (out / "manifest.json").write_text(json.dumps(manifest))
    return out


# --- loader ------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("study", STUDIES, ids=lambda p: p.name)
def test_values_match_source_csv_and_missing_stays_missing(study):
    res = R.load_result(study)
    raw = list(csv.DictReader(open(study / "windows.csv")))
    n_cmp = 0
    for mk, (col, unit, size, _) in R.METRICS.items():
        df = R.metric_frame(res, mk)
        assert df.attrs["unit"] == unit
        for row, (_, d) in zip(raw, df.iterrows()):
            if row[col] == "":
                assert math.isnan(d["value"]) and d["status"] != "value"     # never zero
            else:
                assert d["value"] == float(row[col]) and d["status"] == "value"
                n_cmp += 1
    assert n_cmp > 900
    assert all(ok for ok in res.integrity.values())


def test_kras_values_agree_with_frozen_study_tables():
    """The viewer shows the structpair reproductions; spot-check them against the frozen KRAS study tables."""
    pairs = [(DEMO / "kras1", ROOT / "benchmarks/audit/kras_struct_cmp/results/KRAS-STRUCT-1_table.csv"),
             (DEMO / "kras2", ROOT / "benchmarks/audit/kras_struct_g12d/results/KRAS-STRUCT-2_table.csv")]
    for study, table in pairs:
        res = R.load_result(study)
        frozen = {r["resid"]: r for r in csv.DictReader(open(table))}
        for mk, (col, *_ ) in R.METRICS.items():
            df = R.metric_frame(res, mk)
            for _, d in df.iterrows():
                f = frozen[d["A_residue"].split(":")[1]][mk]
                assert (f == "" and math.isnan(d["value"])) or abs(float(f) - d["value"]) <= 1e-5 * max(1, abs(float(f)))


def test_synthetic_insertion_codes_missing_residue_and_highlight(tmp_path):
    res = R.load_result(make_synthetic(tmp_path))
    assert list(res.residues["A_residue"]) == [f"A:{x}" for x in A_RES]            # 52A preserved
    assert list(res.residues["B_residue"]) == [f"B:{x}" for x in B_RES]            # chain IDs preserved
    df = R.metric_frame(res, "rmsd3")
    for p in (3, 4, 5):                                                             # windows touching position 4
        assert math.isnan(df.loc[p, "value"]) and df.loc[p, "status"] == "excluded"
    assert df.loc[0, "status"] == "outside mapped range" and df.loc[7, "status"] == "outside mapped range"
    assert df.loc[2, "value"] == pytest.approx(0.21)
    assert df.loc[2, "segment"] != df.loc[6, "segment"]                             # no line across the gap
    assert R.window_members(res, 3, 3) == [2, 3, 4]
    assert list(res.residues.iloc[R.window_members(res, 3, 4)]["A_residue"]) == ["A:52", "A:52A", "A:53", "A:54"]
    ex = R.exclusion_explanations(res, 4, 3)
    assert list(ex["cause"]) == ["missing_coordinates"] and "no Cα coordinates" in ex["explanation"].iloc[0]
    coords = {}
    for k in "AB":
        f, msg = R.locate_input(res, k, ROOT)
        assert f is not None, msg
        coords[k] = R.read_ca(f, k, "reject")
    assert "A:52A" in coords["A"] and "B:154" not in coords["B"]
    A, B = R.display_coordinates(res, coords)
    assert np.isnan(B[B_MISSING]).all() and np.isfinite(A[B_MISSING]).all()


def test_invalid_inputs_are_refused(tmp_path):
    good = make_synthetic(tmp_path)
    with pytest.raises(R.ResultError, match="not a directory"):
        R.load_result(tmp_path / "nope")
    d = tmp_path / "missing_file"
    shutil.copytree(good, d)
    (d / "exclusions.csv").unlink()
    with pytest.raises(R.ResultError, match="missing exclusions.csv"):
        R.load_result(d)
    d = tmp_path / "tampered"
    shutil.copytree(good, d)
    with open(d / "windows.csv", "a") as fh:
        fh.write("\n")
    with pytest.raises(R.ResultError, match="differ from the manifest"):
        R.load_result(d)
    d = tmp_path / "schema"
    shutil.copytree(good, d)
    m = json.loads((d / "manifest.json").read_text())
    m["schema"] = "something/else"
    (d / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(R.ResultError, match="schema"):
        R.load_result(d)
    # input structure with a different hash: molecular view disabled, results still readable
    d = tmp_path / "pdb_changed"
    shutil.copytree(good, d)
    res = R.load_result(d)
    (tmp_path / "SYNTHETIC_B.pdb").write_text((tmp_path / "SYNTHETIC_B.pdb").read_text() + "REMARK changed\n")
    f, msg = R.locate_input(res, "B", ROOT)
    assert f is None and "differs" in msg


# --- app -----------------------------------------------------------------------------------------------------------------

def _df(at, column):
    """The first dataframe on the page having this column (read-only dataframes carry no key in AppTest)."""
    for d in at.dataframe:
        if column in d.value.columns:
            return d.value
    raise AssertionError(f"no dataframe with column {column}")


def _app(root):
    at = AppTest.from_file(str(APP_DIR / "app.py"), default_timeout=60)
    at.session_state["results_root"] = str(root)
    return at.run()


def test_app_kras_display_focus_and_sources_unchanged():
    watched = STUDIES + [ROOT / "benchmarks/audit/kras_struct_cmp/data_pdb", ROOT / "benchmarks/audit/kras_struct_g12d/data_pdb"]
    before = R.snapshot([p for p in watched if p.exists()])
    at = _app(DEMO)
    assert not at.exception and at.title[0].value == "TopoFold Structural Comparison"
    assert at.selectbox(key="result_choice").value == "kras1"                    # the KRAS studies come first
    at.selectbox(key="focus").set_value(30).run()
    rows = _df(at, "A_status")
    assert list(rows["A_residue"]) == ["A:30", "A:31", "A:32"] and list(rows["B_residue"]) == ["A:30", "A:31", "A:32"]
    raw = {r["position"]: r for r in csv.DictReader(open(DEMO / "kras1/windows.csv"))}
    shown = " ".join(m.value for m in at.markdown)
    assert f"d_kappa = {float(raw['30']['d_kappa_rad']):.4f} rad" in shown
    at.segmented_control(key="hl_size").set_value(4).run()
    assert len(_df(at, "A_status")) == 4
    for f in ("with values", "without values", "altloc in window", "all"):
        at.segmented_control(key="row_filter").set_value(f).run()
        assert not at.exception
    at.selectbox(key="sort_by").set_value("rmsd4").run()
    table = at.dataframe(key="win_table").value
    assert table["rmsd4_angstrom"].iloc[0] == max(float(r["rmsd4_angstrom"]) for r in raw.values() if r["rmsd4_angstrom"])
    at.pills(key="metrics").set_value(["d_tau", "d_i_i3"]).run()
    at.selectbox(key="result_choice").set_value("kras2").run()
    assert not at.exception
    ex = _df(at, "explanation")
    assert len(ex) == 10 and set(ex["cause"]) == {"numbering_discontinuity"}
    assert R.snapshot([p for p in watched if p.exists()]) == before                # viewing changed nothing


def test_app_synthetic_missing_residue_and_invalid_folder(tmp_path):
    out = make_synthetic(tmp_path)
    at = _app(out)
    assert not at.exception
    at.selectbox(key="focus").set_value(4).run()
    warn = " ".join(w.value for w in at.warning)
    assert "missing_coordinates" in warn and "B:154" in warn
    assert "unavailable" not in warn
    at = _app(tmp_path / "empty_dir_that_does_not_exist")
    assert not at.exception and "No structpair results" in at.error[0].value
    bad = tmp_path / "bad"
    shutil.copytree(out, bad)
    (bad / "summary.json").write_text("{}")                                          # hash no longer matches
    at = _app(bad)
    assert not at.exception and "cannot be displayed" in at.error[0].value


def test_viewer_does_not_import_core_or_mdanalysis():
    import subprocess
    code = ("import sys; from streamlit.testing.v1 import AppTest; "
            f"at = AppTest.from_file({str(APP_DIR / 'app.py')!r}, default_timeout=60); "
            f"at.session_state['results_root'] = {str(DEMO / 'kras2')!r}; at.run(); "
            "assert not at.exception; "
            "print(sorted(m for m in sys.modules if m.split('.')[0] in ('topofold', 'topofold_inputs', 'MDAnalysis')))")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr[-1500:]
    assert r.stdout.strip().splitlines()[-1] == "[]"
