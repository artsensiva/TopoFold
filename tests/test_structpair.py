"""Tests for topofold_inputs.structpair (two-structure comparison). SYNTHETIC inputs only."""

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


def _spec(tmp_path, **kw):
    for name in ("a.pdb", "b.pdb"):
        (tmp_path / name).write_text("END\n")
    s = {"structures": {"A": {"file": "a.pdb", "chain": "A"}, "B": {"file": "b.pdb", "chain": "A"}},
         "correspondence": {"same_numbering": "1-6"}}
    s.update(kw)
    p = tmp_path / "spec.json"
    p.write_text(json.dumps(s))
    return p


def test_spec_ranges_pairs_and_defaults(tmp_path):
    s = sp.load_spec(_spec(tmp_path, correspondence={"same_numbering": "1-3, 7, 9-10"}))
    assert [a for a, _ in s["_pairs"]] == [(1, ""), (2, ""), (3, ""), (7, ""), (9, ""), (10, "")]
    assert s["on_missing"] == "reject" and s["resname_mismatch"] == "reject"
    assert s["structures"]["A"]["altloc"] == "reject"          # no silent altloc choice
    s = sp.load_spec(_spec(tmp_path, correspondence={"pairs": [[1, 11], [2, "12A"], [3, {"resid": 13}], ["4", 14]]}))
    assert s["_pairs"][1] == ((2, ""), (12, "A"))


@pytest.mark.parametrize("change, message", [
    ({"structures": {"A": {"file": "a.pdb"}, "B": {"file": "b.pdb", "chain": "A"}}}, "chain is required"),
    ({"correspondence": {"align": "auto"}}, "never inferred"),
    ({"correspondence": {"pairs": [[1, 1], [2, 2], [3, 3], [1, 4]]}}, "repeats a residue"),
    ({"correspondence": {"same_numbering": "1-3"}}, "at least 4"),
    ({"on_missing": "guess"}, "on_missing"),
    ({"resname_mismatch": "ignore"}, "resname_mismatch"),
    ({"structures": {"A": {"file": "missing.pdb", "chain": "A"}, "B": {"file": "b.pdb", "chain": "A"}}}, "not found"),
])
def test_spec_refusals(tmp_path, change, message):
    with pytest.raises(sp.SpecError, match=message):
        sp.load_spec(_spec(tmp_path, **change))


def test_kabsch_proper_rotation_reference():
    rng = np.random.default_rng(1)
    P = rng.normal(size=(20, 3))
    a = 1.1
    R = np.array([[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]])
    assert sp.kabsch_rmsd(P, P @ R.T + [3, -2, 1]) < 1e-12
    assert sp.kabsch_rmsd(P, P * [1, 1, -1]) > 0.1      # a mirror image is not superposable by a proper rotation


def _side(coords, kappa, tau):
    n = len(coords)
    return {"coords": np.asarray(coords, float), "kappa": np.asarray(kappa, float), "tau": np.asarray(tau, float),
            "kappa_valid": np.isfinite(kappa), "tau_valid": np.isfinite(tau)}


def test_window_metrics_equal_coverage_and_circular_tau():
    n = 6
    X = np.c_[3.8 * np.arange(n), 0.5 * (-1) ** np.arange(n), np.zeros(n)]
    k = np.r_[np.nan, [1.5] * (n - 2), np.nan]
    ta = np.r_[np.nan, [math.pi - 0.1] * (n - 3), np.nan]
    tb = np.r_[np.nan, [-math.pi + 0.1] * (n - 3), np.nan]
    kb = k.copy()
    kb[2] = np.nan                                       # invalid in B only: the whole three-residue window is excluded
    rm = {"rmsd3": [None] + [0.0] * (n - 2) + [None], "rmsd4": [None] + [0.0] * (n - 3) + [None, None]}
    m = sp.window_metrics(_side(X, k, ta), _side(X, kb, tb), rm)
    assert all(m[x][2] is None for x in ("d_kappa", "d_i_i2", "rmsd3"))
    assert m["d_kappa"][1] == 0.0 and m["d_i_i2"][1] == 0.0
    assert all(abs(m["d_tau"][i] - 0.2) < 1e-12 for i in range(1, n - 2))   # wrapped across ±π, not 2π − 0.2


def _pdb(path, coords, chain="A", resnames=None, skip=()):
    lines = []
    for i, (x, y, z) in enumerate(coords):
        r = i + 1
        if r in skip:
            continue
        name = (resnames or {}).get(r, "ALA")
        lines.append(f"ATOM  {len(lines) + 1:5d}  CA  {name} {chain}{r:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C")
    path.write_text("\n".join(lines) + "\nEND\n")


def _helix(n):
    t = np.arange(n) * math.radians(100)
    return np.c_[2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)]


def _adapter_python():
    cand = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
    if not Path(cand).is_file():
        return None
    return cand if subprocess.run([cand, "-c", "import MDAnalysis"], capture_output=True).returncode == 0 else None


def test_end_to_end_synthetic_pair(tmp_path):
    py = _adapter_python()
    if py is None:
        pytest.skip("no interpreter with MDAnalysis available")
    n = 14
    X = _helix(n)
    a = 0.7
    R = np.array([[math.cos(a), -math.sin(a), 0], [math.sin(a), math.cos(a), 0], [0, 0, 1]])
    Y = X @ R.T + [5.0, -3.0, 2.0]
    Y[9] += [0.0, 0.0, 0.6]                       # one local perturbation at residue 10
    _pdb(tmp_path / "SYNTHETIC_A.pdb", X)
    _pdb(tmp_path / "SYNTHETIC_B.pdb", Y, chain="B", resnames={3: "GLY"}, skip=(6,))
    spec = {"name": "SYNTHETIC pair", "correspondence": {"same_numbering": f"1-{n}"},
            "structures": {"A": {"file": "SYNTHETIC_A.pdb", "chain": "A"}, "B": {"file": "SYNTHETIC_B.pdb", "chain": "B"}},
            "on_missing": "exclude", "resname_mismatch": "flag"}
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=py)
    run = lambda d: subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(tmp_path / "spec.json"),  # noqa: E731
                                    "--out", str(tmp_path / d)], env=env, capture_output=True, text=True)
    r = run("o1")
    assert r.returncode == 0, r.stderr
    assert run("o1").returncode == 2                     # refuses a non-empty output directory
    assert run("o2").returncode == 0
    assert (tmp_path / "o1/windows.csv").read_bytes() == (tmp_path / "o2/windows.csv").read_bytes()
    win = {int(w["A_residue"].split(":")[1]): w for w in csv.DictReader(open(tmp_path / "o1/windows.csv"))}
    # rigid motion removed: windows away from residues 6 and 10 show no difference
    assert float(win[3]["rmsd3_angstrom"]) < 1e-3 and float(win[13]["rmsd3_angstrom"]) < 1e-3
    assert float(win[10]["rmsd3_angstrom"]) > 0.1 and float(win[10]["d_kappa_rad"]) > 0.01
    # residue 6 missing in B: its windows are excluded for every metric, with the cause named
    for c in (5, 6, 7):
        assert win[c]["d_kappa_rad"] == "" and win[c]["rmsd3_angstrom"] == "" and win[c]["d_i_i2_angstrom"] == ""
    ex = list(csv.DictReader(open(tmp_path / "o1/exclusions.csv")))
    assert {(e["structure"], e["cause"], e["location"]) for e in ex} >= {("B", "missing_coordinates", "B:6")}
    assert win[3]["mismatch_in_window3"] == "True" and win[3]["resname_B"] == "GLY"
    summ = json.loads((tmp_path / "o1/summary.json").read_text())
    assert abs(summ["global"]["rmsd_qcp_angstrom"] - summ["global"]["rmsd_numpy_kabsch_angstrom"]) < 1e-6
    man = json.loads((tmp_path / "o1/manifest.json").read_text())
    assert set(man["outputs_sha256"]) >= {"windows.csv", "report.md", "summary.json"}
    assert man["inputs"]["B"]["chain"] == "B" and len(man["inputs"]["A"]["sha256"]) == 64
