"""structpair 0.2.0: several explicitly selected chains per entry (within-entry crystal-copy variation).

Pure-function tests and SYNTHETIC end-to-end tests. They check implementation consistency against analytic cases and
against the frozen K1 implementation, not biological validity.
"""

import csv
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/kras_k1"))
from topofold_inputs import chain_variation as cv  # noqa: E402
from topofold_inputs import structpair as sp  # noqa: E402

FIX = ROOT / "tests/fixtures/structpair_v0_1_0"


def _adapter_python():
    cand = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
    if not Path(cand).is_file():
        return None
    return cand if subprocess.run([cand, "-c", "import MDAnalysis"], capture_output=True).returncode == 0 else None


def _rot(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    return q if np.linalg.det(q) > 0 else q * np.array([1, 1, -1])


# --- pure functions ------------------------------------------------------------------------------------------------------

def test_aggregate_statistics_median_of_even_count_and_missing_values():
    a = cv.aggregate([3.0, None, 1.0, float("nan"), 2.0, 10.0])
    assert (a["n"], a["min"], a["median"], a["max"]) == (4, 1.0, 2.5, 10.0)       # even count: mean of the middle two
    assert a["mean"] == pytest.approx(4.0)
    e = cv.aggregate([None, float("nan")])
    assert e == {"n": 0, "min": None, "median": None, "mean": None, "max": None}   # empty, never zero


def test_pair_displacement_equals_the_frozen_k1_implementation_and_recovers_a_local_shift():
    import k1_endpoint as K
    rng = np.random.default_rng(3)
    n = 50
    ref = rng.normal(size=(n, 3)) * 8
    mov = ref.copy()
    mov[10:15] += np.array([3.0, 0.0, 0.0])
    mov = mov @ _rot(rng).T + np.array([1.0, 2.0, 3.0])
    d = cv.pair_displacement(ref, mov, 0.2)
    assert d[10:15] == pytest.approx(3.0, abs=1e-9) and np.delete(d, range(10, 15)).max() < 1e-9
    k1 = K.pair_displacement(ref, mov, list(range(n)), 0.2)
    assert np.abs(d - k1).max() < 1e-9                                  # implementation consistency with the K1 record's code
    noisy = ref + rng.normal(size=ref.shape) * 0.4
    assert np.abs(cv.pair_displacement(ref, noisy, 0.2) - K.pair_displacement(ref, noisy, list(range(n)), 0.2)).max() < 1e-9


def test_pair_displacement_missing_positions_and_too_few_positions():
    rng = np.random.default_rng(4)
    a = rng.normal(size=(12, 3)) * 5
    b = a + 0.1
    b[3] = np.nan
    d = cv.pair_displacement(a, b, 0.2)
    assert np.isnan(d[3]) and np.isfinite(np.delete(d, 3)).all()
    few = a.copy()
    few[3:] = np.nan
    assert np.isnan(cv.pair_displacement(few, a, 0.2)).all()              # three positions in both: fewer than 4, undefined
    four = a.copy()
    four[4:] = np.nan
    assert np.isfinite(cv.pair_displacement(four, a, 0.2)[:4]).all()      # exactly 4 positions are enough


def test_trim_ties_drop_the_lower_position_first():
    assert cv.trim_indices(np.array([1.0, 1.0, 1.0, 0.5]), 2) == [0, 1]
    assert cv.trim_indices(np.array([2.0, 3.0, 1.0]), 1) == [1]


def test_structpair_global_check_uses_the_same_kabsch_and_keeps_its_old_arithmetic():
    rng = np.random.default_rng(1)
    P = rng.normal(size=(20, 3))
    R = np.array([[1, 0, 0], [0, math.cos(1.1), -math.sin(1.1)], [0, math.sin(1.1), math.cos(1.1)]])
    Q = P @ R.T + [3, -2, 1]
    assert sp.kabsch_rmsd(P, Q) < 1e-12
    r, cp, cq = cv.kabsch_frame(P, Q)
    assert np.linalg.det(r) == pytest.approx(1.0)
    assert sp.kabsch_rmsd(P, P * [1, 1, -1]) > 0.1                         # a mirror image is not superposable


def test_parse_params_defaults_and_refusals():
    p = cv.parse_params(None)
    assert p == {"displacement_trim_fraction": 0.2, "advisory": {"statistic": "min", "within_factor": 1.0, "min_difference": {}}}
    k1 = cv.parse_params({"advisory": {"statistic": "mean", "within_factor": 2, "min_difference": {"disp": 1.0}}})
    assert k1["advisory"] == {"statistic": "mean", "within_factor": 2.0, "min_difference": {"disp": 1.0}}
    for bad in ({"displacement_trim_fraction": 0.6}, {"advisory": {"statistic": "mode"}}, {"advisory": {"within_factor": -1}},
                {"advisory": {"min_difference": {"nope": 1}}}, {"surprise": 1}, {"advisory": {"extra": 1}}):
        with pytest.raises(cv.VariationError):
            cv.parse_params(bad)


def test_advisory_rule_is_inclusive_at_the_boundary_and_never_false_when_not_assessable():
    p = cv.parse_params(None)["advisory"]
    eq = cv.advisory(1.0, 1.0, 4, "both", "disp", p)
    assert eq["separated"] is True and eq["at_boundary"] is True            # tie: equal values count as separated
    assert cv.advisory(1.0 - 1e-6, 1.0, 4, "both", "disp", p)["separated"] is False
    assert cv.advisory(1.0 - 1e-12, 1.0, 4, "both", "disp", p)["separated"] is True       # within 1e-9
    k1 = cv.parse_params({"advisory": {"statistic": "mean", "within_factor": 2.0, "min_difference": {"disp": 1.0}}})["advisory"]
    assert cv.advisory(1.5, 0.2, 12, "both", "disp", k1)["separated"] is True            # floor 1.0 governs (2 x 0.2 < 1.0)
    assert cv.advisory(1.5, 0.9, 12, "both", "disp", k1)["threshold"] == pytest.approx(1.8)
    na = cv.advisory(5.0, None, 4, "none", "disp", p)
    assert na["separated"] is None and "within-entry" in na["reason"]
    assert cv.advisory(None, 1.0, 0, "both", "disp", p)["separated"] is None


def _pv(chains_a, chains_b, n, fn):
    keys = cv.pair_keys(chains_a, chains_b)
    return {k: {ch: [fn(k, ch, i) for i in range(n)] for ch in cv.CHANNELS} for g in keys.values() for k in g}


def test_variation_rows_counts_scope_primary_and_missing_pairs():
    ca, cb, n = ["A", "B"], ["A", "B", "C"], 3
    # position 1: cross pairs with B:C undefined; position 2: no within-entry pair of B defined; values depend on the pair
    def fn(k, ch, i):
        kind, x, y = k
        if i == 1 and "C" in (x, y) and kind != "withinA":
            return None
        v = {"cross": 5.0, "withinA": 1.0, "withinB": 2.0}[kind] + (0.5 if (x, y) == ("A", "A") else 0.0)
        return v
    rows = cv.variation_rows(_pv(ca, cb, n, fn), n, [(f"A:{i}", f"B:{i}") for i in range(n)], [("w3", "w4")] * n, ca, cb, ("A", "A"),
                             cv.parse_params(None))
    r0 = [r for r in rows if r["position"] == 0 and r["channel"] == "disp"][0]
    assert (r0["between_n"], r0["between_min"], r0["between_max"], r0["primary"]) == (6, 5.0, 5.5, 5.5)
    assert (r0["withinA_n"], r0["withinB_n"], r0["within_scope"], r0["within_max"]) == (1, 3, "both", 2.0)
    r1 = [r for r in rows if r["position"] == 1 and r["channel"] == "disp"][0]
    assert r1["between_n"] == 4 and r1["withinB_n"] == 1 and r1["withinB_max"] == 2.0     # only the A-B pair of entry B remains
    assert r0["separated"] is True and r0["at_boundary"] is False
    w = [r for r in rows if r["position"] == 0 and r["channel"] == "rmsd4"][0]
    assert w["window_A"] == "w4" and [r for r in rows if r["channel"] == "disp"][0]["window_A"] == "A:0"


def test_variation_rows_one_entry_without_within_pairs_is_scoped_not_false():
    ca, cb, n = ["A", "B"], ["A"], 2
    rows = cv.variation_rows(_pv(ca, cb, n, lambda k, ch, i: 3.0 if k[0] == "cross" else 1.0), n, [("a", "b")] * n, [("w", "w")] * n,
                             ca, cb, ("A", "A"), cv.parse_params(None))
    r = rows[0]
    assert r["within_scope"] == "A_only" and r["withinB_n"] == 0 and r["withinB_max"] is None and r["separated"] is True
    only_b = cv.variation_rows(_pv(["A"], ["A", "B"], n, lambda k, ch, i: 3.0 if k[0] == "cross" else 1.0), n, [("a", "b")] * n,
                               [("w", "w")] * n, ["A"], ["A", "B"], ("A", "A"), cv.parse_params(None))
    assert only_b[0]["within_scope"] == "B_only"
    none = cv.variation_rows(_pv(["A"], ["A"], n, lambda k, ch, i: 3.0), n, [("a", "b")] * n, [("w", "w")] * n, ["A"], ["A"], ("A", "A"),
                             cv.parse_params(None))
    assert none[0]["within_scope"] == "none" and none[0]["separated"] is None and none[0]["not_assessed_reason"]


def test_summarize_status_and_block_text():
    rows = cv.variation_rows(_pv(["A", "B"], ["A", "B"], 2, lambda k, ch, i: 1.0 if k[0] != "cross" else 2.0), 2, [("a", "b")] * 2,
                             [("w", "w")] * 2, ["A", "B"], ["A", "B"], ("A", "A"), cv.parse_params(None))
    b = cv.summarize(rows, ["A", "B"], ["A", "B"], ("A", "A"), cv.parse_params(None))
    assert b["status"] == "available" and b["n_pairs"] == {"cross": 4, "withinA": 1, "withinB": 1}
    assert b["channels"]["disp"]["separated"] == 2 and "not independent experimental replicates" in b["statement"]
    assert "not a calibrated noise floor" in b["statement"] and "HDX replicate" in b["statement"]
    assert cv.summarize(rows, ["A", "B"], ["A"], ("A", "A"), cv.parse_params(None))["status"] == "partial"
    assert cv.unavailable_block()["status"] == "unavailable" and cv.unavailable_block()["reason"] == "one chain per entry"
    text = "\n".join(cv.report_lines(b, rows, 3))
    assert text.startswith("## Within-entry crystal-copy variation") and "descriptive" in text
    assert "Unavailable: one chain per entry" in "\n".join(cv.report_lines(cv.unavailable_block(), [], 3))


def _spec_file(tmp_path, a=None, b=None, **kw):
    for name in ("a.pdb", "b.pdb"):
        (tmp_path / name).write_text("END\n")
    s = {"structures": {"A": {"file": "a.pdb", "chain": "A"}, "B": {"file": "b.pdb", "chain": "A"}},
         "correspondence": {"same_numbering": "1-6"}}
    if a is not None:
        s["structures"]["A"] = dict({"file": "a.pdb"}, **a)
    if b is not None:
        s["structures"]["B"] = dict({"file": "b.pdb"}, **b)
    s.update(kw)
    p = tmp_path / "spec.json"
    p.write_text(json.dumps(s))
    return p


def test_spec_chains_are_explicit_and_the_primary_is_never_chosen_for_the_user(tmp_path):
    s = sp.load_spec(_spec_file(tmp_path, a={"chains": ["A", "B"], "primary_chain": "B"}, b={"chains": ["C"]}))
    assert s["structures"]["A"]["chain"] == "B" and s["structures"]["A"]["chains"] == ["A", "B"]
    assert s["structures"]["B"]["chain"] == "C" and s["_variation"]["displacement_trim_fraction"] == 0.2
    legacy = sp.load_spec(_spec_file(tmp_path))
    assert legacy["structures"]["A"]["chains"] == ["A"]
    for a, msg in (({"chains": ["A", "B"]}, "primary_chain is required"), ({"chains": ["A", "B"], "primary_chain": "Z"}, "primary_chain"),
                   ({"chains": ["A", "A"], "primary_chain": "A"}, "repeats"), ({"chains": [], "primary_chain": "A"}, "non-empty"),
                   ({"chain": "A", "chains": ["A", "B"], "primary_chain": "A"}, "either chain or chains"),
                   ({"chain": "A", "primary_chain": "A"}, "needs"), ({"chains": "AB", "primary_chain": "A"}, "non-empty list")):
        with pytest.raises(sp.SpecError, match=msg):
            sp.load_spec(_spec_file(tmp_path, a=a))
    with pytest.raises(sp.SpecError, match="chain is required"):
        sp.load_spec(_spec_file(tmp_path, a={}))
    with pytest.raises(sp.SpecError, match="advisory"):
        sp.load_spec(_spec_file(tmp_path, chain_variation={"advisory": {"statistic": "mode"}}))


# --- SYNTHETIC end-to-end ---------------------------------------------------------------------------------------------------

def _helix(n):
    t = np.arange(n) * math.radians(100)
    return np.c_[2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)]


def _multi_pdb(path, chains: dict):
    """chains: {chain: (coords, skip residue numbers)}; one file with several chains, residue names ALA."""
    lines = []
    for ch, (xyz, skip) in chains.items():
        for i, (x, y, z) in enumerate(xyz):
            if i + 1 in skip:
                continue
            lines.append(f"ATOM  {len(lines) + 1:5d}  CA  ALA {ch}{i + 1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C")
        lines.append("TER")
    path.write_text("\n".join(lines) + "\nEND\n")


def _scenario(tmp_path):
    """SYNTHETIC: entry A has chains A, B; entry B has chains A, B, C. Residues 8-10 of entry B are displaced by 2.0 A; chain C of
    entry B lacks residue 3; copies differ by small coordinate noise and rigid motions."""
    rng = np.random.default_rng(11)
    n = 16
    X = _helix(n)
    Y = X.copy()
    Y[7:10] += np.array([0.0, 0.0, 2.0])
    mv = lambda P, s: (P + rng.normal(size=P.shape) * s) @ _rot(rng).T + rng.normal(size=3) * 5   # noqa: E731
    _multi_pdb(tmp_path / "SYNTHETIC_A.pdb", {"A": (mv(X, 0.0), ()), "B": (mv(X, 0.05), ())})
    _multi_pdb(tmp_path / "SYNTHETIC_B.pdb", {"A": (mv(Y, 0.0), ()), "B": (mv(Y, 0.05), ()), "C": (mv(Y, 0.05), (3,))})
    spec = {"name": "SYNTHETIC chains", "correspondence": {"same_numbering": f"1-{n}"}, "on_missing": "exclude",
            "structures": {"A": {"label": "SYNTHETIC A", "file": "SYNTHETIC_A.pdb", "chains": ["A", "B"], "primary_chain": "A"},
                           "B": {"label": "SYNTHETIC B", "file": "SYNTHETIC_B.pdb", "chains": ["A", "B", "C"], "primary_chain": "A"}}}
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    return spec


def _run_cli(spec_path, out, py, extra_env=None):
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=py, **(extra_env or {}))
    return subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(spec_path), "--out", str(out)], env=env,
                          capture_output=True, text=True)


def test_single_chain_outputs_are_byte_identical_to_version_0_1_0(tmp_path):
    py = _adapter_python()
    if py is None:
        pytest.skip("no interpreter with MDAnalysis available")
    for f in FIX.glob("SYNTHETIC_[AB].pdb"):
        shutil.copy(f, tmp_path / f.name)
    shutil.copy(FIX / "SYNTHETIC_spec.json", tmp_path / "spec.json")
    r = _run_cli(tmp_path / "spec.json", tmp_path / "out", py)
    assert r.returncode == 0, r.stderr
    for name in ("windows", "residues", "exclusions"):
        assert (tmp_path / f"out/{name}.csv").read_bytes() == (FIX / f"SYNTHETIC_expected_{name}.csv").read_bytes()
    new, old = json.loads((tmp_path / "out/summary.json").read_text()), json.loads((FIX / "SYNTHETIC_expected_summary_numeric.json").read_text())
    assert all(new[k] == old[k] for k in old)                                  # every legacy number unchanged
    assert new["within_entry_variation"] == cv.unavailable_block()
    man = json.loads((tmp_path / "out/manifest.json").read_text())
    assert man["schema"] == "topofold_inputs.structpair/2" and man["software"]["structpair_version"] == "0.2.0"
    assert man["within_entry_variation"]["status"] == "unavailable"
    assert not (tmp_path / "out/chain_variation.csv").exists() and not (tmp_path / "out/chains.csv").exists()
    assert "Unavailable: one chain per entry" in (tmp_path / "out/report.md").read_text()


def test_multi_chain_end_to_end_synthetic(tmp_path):
    py = _adapter_python()
    if py is None:
        pytest.skip("no interpreter with MDAnalysis available")
    _scenario(tmp_path)
    r = _run_cli(tmp_path / "spec.json", tmp_path / "out", py)
    assert r.returncode == 0, r.stderr
    assert _run_cli(tmp_path / "spec.json", tmp_path / "out", py).returncode == 2           # refuses a non-empty directory
    assert _run_cli(tmp_path / "spec.json", tmp_path / "out2", py).returncode == 0
    for f in ("chain_variation.csv", "chains.csv", "windows.csv"):                          # deterministic reruns
        assert (tmp_path / "out" / f).read_bytes() == (tmp_path / "out2" / f).read_bytes()
    rows = list(csv.DictReader(open(tmp_path / "out/chain_variation.csv")))
    disp = {int(x["A_residue"].split(":")[1]): x for x in rows if x["channel"] == "disp"}
    # pair counts: 6 cross pairings, 1 pair inside A, 3 inside B
    assert (disp[1]["between_n"], disp[1]["withinA_n"], disp[1]["withinB_n"]) == ("6", "1", "3")
    # the injected 2.0 A shift of residues 8-10 appears in every cross pairing; the rest agrees within noise
    for r_ in (8, 9, 10):
        assert float(disp[r_]["between_min"]) > 1.8 and float(disp[r_]["between_max"]) < 2.2 and disp[r_]["separated"] == "True"
    assert max(float(disp[r_]["between_max"]) for r_ in (1, 2, 4, 5, 13, 14, 15, 16)) < 0.4
    assert max(float(x["within_max"]) for x in disp.values()) < 0.4                           # copies differ by noise only
    # chain C lacks residue 3: its pairings are undefined there, the other pairs keep their values and the counts say so
    assert (disp[3]["between_n"], disp[3]["withinB_n"], disp[3]["withinA_n"]) == ("4", "1", "1")
    chains = {(c["entry"], c["chain"]): c for c in csv.DictReader(open(tmp_path / "out/chains.csv"))}
    assert chains[("B", "C")]["n_missing"] == "1" and chains[("B", "C")]["missing"] == "C:3" and chains[("A", "A")]["primary"] == "True"
    # the legacy tables are the primary pair (A:A vs B:A) and the primary column repeats them
    win = {int(w["A_residue"].split(":")[1]): w for w in csv.DictReader(open(tmp_path / "out/windows.csv"))}
    prim = {int(x["A_residue"].split(":")[1]): x for x in rows if x["channel"] == "rmsd3"}
    for i in (2, 9, 12):
        assert float(prim[i]["primary"]) == pytest.approx(float(win[i]["rmsd3_angstrom"]), abs=1e-6)
    summ = json.loads((tmp_path / "out/summary.json").read_text())["within_entry_variation"]
    assert summ["status"] == "available" and summ["n_pairs"] == {"cross": 6, "withinA": 1, "withinB": 3}
    assert summ["channels"]["disp"]["separated"] >= 3 and "not independent experimental replicates" in summ["statement"]
    rep = (tmp_path / "out/report.md").read_text()
    assert rep.index("## Within-entry crystal-copy variation") < rep.index("## Mapping") and "## Exclusions" in rep
    man = json.loads((tmp_path / "out/manifest.json").read_text())
    assert {"chain_variation.csv", "chains.csv", "windows.csv"} <= set(man["outputs_sha256"])
    assert man["inputs"]["B"]["chains"] == ["A", "B", "C"] and man["inputs"]["B"]["chain"] == "A"
    # the batch QCP gives the single-pair QCP value for a non-primary pair (the pair A, B inside entry B)
    work = tmp_path / "out/work"
    single = work / "single.json"
    subprocess.run([py, "-m", "topofold_inputs.qcp_rmsd", str(work / "primary_B.npz"), str(work / "chain_B_B.npz"), str(single)],
                   check=True, env=dict(os.environ, PYTHONPATH=str(ROOT / "python")))
    batch = json.loads((work / "chain_pairs_qcp.json").read_text())["results"]["withinB|A|B"]
    one = json.loads(single.read_text())
    assert one["global"] == batch["global"] and one["rmsd3"] == batch["rmsd3"] and one["rmsd4"] == batch["rmsd4"]


def test_multi_chain_refuses_residue_name_disagreement_between_chains_of_an_entry(tmp_path):
    py = _adapter_python()
    if py is None:
        pytest.skip("no interpreter with MDAnalysis available")
    spec = _scenario(tmp_path)
    txt = (tmp_path / "SYNTHETIC_B.pdb").read_text().replace(" ALA C   5", " GLY C   5")
    (tmp_path / "SYNTHETIC_B.pdb").write_text(txt)
    r = _run_cli(tmp_path / "spec.json", tmp_path / "out", py)
    assert r.returncode == 2 and "chains of one entry must agree" in r.stderr
    spec["resname_mismatch"] = "flag"
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    assert _run_cli(tmp_path / "spec.json", tmp_path / "out_flag", py).returncode == 0


def test_k1_style_advisory_setting_flags_only_the_displaced_residues(tmp_path):
    py = _adapter_python()
    if py is None:
        pytest.skip("no interpreter with MDAnalysis available")
    spec = _scenario(tmp_path)
    spec["chain_variation"] = {"displacement_trim_fraction": 0.2,
                               "advisory": {"statistic": "mean", "within_factor": 2.0, "min_difference": {"disp": 1.0}}}
    (tmp_path / "spec.json").write_text(json.dumps(spec))
    r = _run_cli(tmp_path / "spec.json", tmp_path / "out", py)
    assert r.returncode == 0, r.stderr
    rows = [x for x in csv.DictReader(open(tmp_path / "out/chain_variation.csv")) if x["channel"] == "disp"]
    flagged = sorted(int(x["A_residue"].split(":")[1]) for x in rows if x["separated"] == "True")
    assert flagged == [8, 9, 10]
