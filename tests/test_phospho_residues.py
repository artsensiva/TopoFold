"""Opt-in support for documented phosphorylated residues (SEP, TPO, PTR) in the PDB-to-structpair workflow (KNOWN_ISSUES m14). SYNTHETIC inputs.

Design: benchmarks/audit/structpair_phospho_support_2026-10-26/DESIGN.md (written before the implementation)."""

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
sys.path.insert(0, str(ROOT / "benchmarks/audit/kin_pair_1"))
import kinpair as kp  # noqa: E402
from topofold_inputs import numeric_record as nr  # noqa: E402
from topofold_inputs import structpair as sp  # noqa: E402
from topofold_inputs import structpair_verify as V  # noqa: E402


def _adapter_python():
    cand = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
    if not Path(cand).is_file():
        return None
    return cand if subprocess.run([cand, "-c", "import MDAnalysis"], capture_output=True).returncode == 0 else None


PY = _adapter_python()
pytestmark = pytest.mark.skipif(PY is None, reason="no interpreter with MDAnalysis available")
PARENT = {"SEP": "SER", "TPO": "THR", "PTR": "TYR"}
BACKBONE = ("N", "CA", "C", "O")


def _res(resseq, name="ALA", icode="", record="ATOM  ", atoms=BACKBONE, alt=None):
    return {"resseq": resseq, "icode": icode, "name": name, "record": record, "atoms": atoms, "alt": alt}


def _chain(n=14, mods=None, over=None):
    """n residues numbered 1..n; mods {resseq: (name, record)}; over {resseq: dict of overrides of the residue spec}."""
    out = []
    for i in range(n):
        s = i + 1
        name, rec = (mods or {}).get(s, ("ALA", "ATOM  "))
        out.append(_res(s, name, record=rec, **(over or {}).get(s, {})))
    return out


def _helix(n, jump_after=None):
    t = np.arange(n) * math.radians(100)
    X = np.c_[2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)]
    if jump_after is not None:
        X[jump_after:] += np.array([0.0, 0.0, 6.0])
    return X


def _write(path, residues, chain="A", modres=(), shift=0.0, jump_after=None, extra=(), tag="1ABC"):
    n = len(residues)
    X = _helix(n, jump_after) + shift
    a = 0.326                                          # C(i)-N(i+1) = 3.83 * (1 - 2a) = 1.33 Å for equal Cα spacing (exact in the helix)
    lines = [f"MODRES {tag} {m[0]:>3s} {m[1]} {m[2]:4d}{m[3]:1s} {m[4]:>3s}  {m[5]}" for m in modres]
    k = 0
    for i, r in enumerate(residues):
        fwd = X[min(i + 1, n - 1)] - X[i] if i < n - 1 else X[i] - X[i - 1]
        back = X[i] - X[i - 1] if i else X[1] - X[0]
        xyz = {"N": X[i] - a * back, "CA": X[i], "C": X[i] + a * fwd, "O": X[i] + a * fwd + np.array([0, 0, 1.0])}
        for an in r["atoms"]:
            alts = [(" ", 1.0, 0.0)] if not r["alt"] else r["alt"]
            for lab, occ, off in alts:
                p = xyz[an] + off
                k += 1
                lines.append(f"{r['record']:<6s}{k:5d} {an:^4s}{lab}{r['name']:>3s} {chain}{r['resseq']:4d}{r['icode']:1s}   {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}{occ:6.2f} 20.00           {an[0]}")
    lines += list(extra) + ["TER", "END"]
    Path(path).write_text("\n".join(lines) + "\n")


def _mod(resseq, name, chain="A", icode="", parent=None, comment="MODIFIED RESIDUE"):
    return (name, chain, resseq, icode, parent or PARENT[name], comment)


def _adapter(path, out, extra_args=(), correspondence=None, altloc="blank_then_A", policy="phospho_pdb"):
    args = [PY, "-m", "topofold_inputs", str(path), "--run-id", "x", "--out", str(out), "--altloc", altloc, *extra_args]
    if policy:
        args += ["--modified-residues", policy]
    if correspondence is not None:
        cf = Path(out).parent / "corr.json"
        cf.write_text(json.dumps(correspondence))
        args += ["--correspondence", str(cf)]
    Path(out).mkdir(parents=True, exist_ok=True)
    r = subprocess.run(args, capture_output=True, text=True, env=dict(os.environ, PYTHONPATH=str(ROOT / "python")))
    rep = json.loads((Path(out) / "x.json").read_text()) if r.returncode == 0 else None
    return r, rep


def _by(rep, resid):
    return [x for x in rep["residues"] if x["resid"] == resid][0]


# --- adapter ----------------------------------------------------------------------------------------------------------

def test_documented_sep_tpo_ptr_are_admitted_with_identity_parent_source_and_original_coordinates(tmp_path):
    mods = {3: ("SEP", "HETATM"), 6: ("TPO", "HETATM"), 9: ("PTR", "ATOM  ")}
    _write(tmp_path / "a.pdb", _chain(mods=mods), modres=[_mod(3, "SEP"), _mod(6, "TPO"), _mod(9, "PTR", comment="PHOSPHORYLATED")])
    r, rep = _adapter(tmp_path / "a.pdb", tmp_path / "o")
    assert r.returncode == 0, r.stderr
    for resid, name in ((3, "SEP"), (6, "TPO"), (9, "PTR")):
        x = _by(rep, resid)
        assert x["resname"] == name and x["status"] == "present" and x["chain"] == "A" and x["icode"] == ""     # deposited name kept, never the parent
        m = x["modification"]
        assert m["deposited_name"] == name and m["parent"] == PARENT[name] and m["record_type"] in ("HETATM", "ATOM")
        assert "MODRES" in m["parent_source"] and "policy table" in m["parent_source"] and m["annotation"] in ("MODIFIED RESIDUE", "PHOSPHORYLATED")
    assert [x["resname"] for x in rep["residues"]].count("ALA") == 11
    blk = rep["modified_residues"]
    assert blk["policy"] == "phospho_pdb" and len(blk["admitted"]) == 3 and blk["excluded"] == []
    xyz = np.load(tmp_path / "o/x.npz")["coordinates_angstrom"][0]
    want = _helix(14)
    assert np.abs(xyz - want).max() < 1e-3                                              # original coordinates, as written in the file (3 decimals)


def test_default_policy_keeps_historical_behaviour_for_the_same_file(tmp_path):
    _write(tmp_path / "a.pdb", _chain(mods={3: ("SEP", "HETATM")}), modres=[_mod(3, "SEP")])
    r, rep = _adapter(tmp_path / "a.pdb", tmp_path / "o", policy=None)
    assert r.returncode == 0 and 3 not in [x["resid"] for x in rep["residues"]] and "modified_residues" not in rep
    assert all("modification" not in x for x in rep["residues"])
    _write(tmp_path / "s.pdb", _chain())
    r1, a = _adapter(tmp_path / "s.pdb", tmp_path / "o1", policy=None)
    r2, b = _adapter(tmp_path / "s.pdb", tmp_path / "o2", policy="phospho_pdb")
    assert r1.returncode == r2.returncode == 0
    assert np.array_equal(np.load(tmp_path / "o1/x.npz")["coordinates_angstrom"], np.load(tmp_path / "o2/x.npz")["coordinates_angstrom"])
    assert a["residues"] == b["residues"] and a["breaks"] == b["breaks"]               # standard residues: nothing changes


def test_missing_ca_and_ca_only_inputs_are_not_admitted_and_say_why(tmp_path):
    res = _chain(mods={3: ("SEP", "HETATM"), 7: ("TPO", "HETATM")}, over={3: {"atoms": ("N", "C", "O")}, 7: {"atoms": ("CA",)}})
    _write(tmp_path / "a.pdb", res, modres=[_mod(3, "SEP"), _mod(7, "TPO")])
    r, rep = _adapter(tmp_path / "a.pdb", tmp_path / "o")
    assert r.returncode == 0, r.stderr
    assert 3 not in [x["resid"] for x in rep["residues"]] and 7 not in [x["resid"] for x in rep["residues"]]
    ex = {e["resid"]: e for e in rep["modified_residues"]["excluded"]}
    assert "Cα" in ex[3]["reason"] and "N, CA and C" in ex[3]["reason"]
    assert "backbone" in ex[7]["reason"] and "Cα-only" in ex[7]["reason"] and ex[7]["name"] == "TPO"
    assert rep["modified_residues"]["admitted"] == []


def test_a_real_gap_and_a_numbering_gap_keep_their_break_records(tmp_path):
    res = _chain(mods={5: ("SEP", "HETATM")})
    _write(tmp_path / "a.pdb", res, modres=[_mod(5, "SEP")], jump_after=7)                # a 6 Å jump after residue 7: geometric break
    r, rep = _adapter(tmp_path / "a.pdb", tmp_path / "o")
    assert r.returncode == 0, r.stderr
    idx = [x["resid"] for x in rep["residues"]]
    assert 5 in idx and any(b["reasons"] == ["ca_ca_distance"] and idx[b["after_index"]] == 7 for b in rep["breaks"])
    res2 = _chain(mods={5: ("SEP", "HETATM")})
    for r_ in res2[5:]:
        r_["resseq"] += 1                                                               # residue number 6 absent: numbering discontinuity after 5
    _write(tmp_path / "b.pdb", res2, modres=[_mod(5, "SEP")])
    r, rep = _adapter(tmp_path / "b.pdb", tmp_path / "o2")
    idx = [x["resid"] for x in rep["residues"]]
    assert any(b["reasons"] == ["numbering_discontinuity"] and idx[b["after_index"]] == 5 for b in rep["breaks"])


def test_missing_or_contradictory_annotation_is_an_explicit_refusal_and_other_chains_are_out_of_scope(tmp_path):
    res = _chain(mods={4: ("SEP", "HETATM")})
    _write(tmp_path / "none.pdb", res)                                                  # bonded SEP, no MODRES
    r, _ = _adapter(tmp_path / "none.pdb", tmp_path / "o1")
    assert r.returncode == 2 and "SEP" in r.stderr and "no MODRES record" in r.stderr and "A:4" in r.stderr
    _write(tmp_path / "parent.pdb", res, modres=[_mod(4, "SEP", parent="THR")])          # parent contradicts the policy table
    r, _ = _adapter(tmp_path / "parent.pdb", tmp_path / "o2")
    assert r.returncode == 2 and "contradict" in r.stderr and "THR" in r.stderr
    _write(tmp_path / "name.pdb", res, modres=[("TPO", "A", 4, "", "THR", "MODIFIED RESIDUE")])   # MODRES names another residue
    r, _ = _adapter(tmp_path / "name.pdb", tmp_path / "o3")
    assert r.returncode == 2 and "contradict" in r.stderr
    # an unannotated bonded TPO in chain B is outside the requested correspondence: no error, not admitted, not output
    both = tmp_path / "two.pdb"
    _write(both, _chain(), extra=[])
    txt = both.read_text().replace("TER\nEND\n", "")
    _write(tmp_path / "b.pdb", _chain(mods={4: ("TPO", "HETATM")}), chain="B")
    both.write_text(txt + (tmp_path / "b.pdb").read_text())
    corr = [{"chain": "A", "resid": i, "icode": "", "label": f"p{i}"} for i in range(1, 15)]
    r, rep = _adapter(both, tmp_path / "o4", correspondence=corr)
    assert r.returncode == 0, r.stderr
    assert [x["chain"] for x in rep["residues"]] == ["A"] * 14 and rep["modified_residues"]["admitted"] == []


def test_insertion_codes_and_alternate_locations(tmp_path):
    res = [_res(1), _res(2), _res(3), _res(3, "SEP", icode="A", record="HETATM", alt=[(" ", 1.0, 0.0), ("B", 0.4, 0.3)]), _res(4), _res(5), _res(6)]
    _write(tmp_path / "a.pdb", res, modres=[_mod(3, "SEP", icode="A")])
    r, rep = _adapter(tmp_path / "a.pdb", tmp_path / "o")
    assert r.returncode == 0, r.stderr
    x = [x for x in rep["residues"] if x["resname"] == "SEP"][0]
    assert (x["resid"], x["icode"]) == (3, "A") and x["modification"]["icode"] == "A" and x["altloc"] == ""      # blank_then_A takes the unlabelled atom
    assert not rep["breaks"]                                                            # 3 -> 3A -> 4 is consecutive numbering
    r, rep = _adapter(tmp_path / "a.pdb", tmp_path / "o2", altloc="highest_occupancy")
    assert r.returncode == 0 and [x for x in rep["residues"] if x["resname"] == "SEP"][0]["altloc"] == ""
    tie = [_res(1), _res(2), _res(3, "SEP", record="HETATM", alt=[("A", 0.5, 0.0), ("B", 0.5, 0.3)]), _res(4), _res(5), _res(6)]
    _write(tmp_path / "t.pdb", tie, modres=[_mod(3, "SEP")])
    r, _ = _adapter(tmp_path / "t.pdb", tmp_path / "o3", altloc="highest_occupancy")
    assert r.returncode == 2 and "tie on occupancy" in r.stderr                         # the existing explicit tie refusal applies to the modified residue too
    r, rep = _adapter(tmp_path / "t.pdb", tmp_path / "o4", altloc="blank_then_A")
    assert r.returncode == 0 and [x for x in rep["residues"] if x["resname"] == "SEP"][0]["altloc"] == "A"


def test_a_free_component_and_a_calcium_ion_are_not_admitted(tmp_path):
    res = _chain(mods={3: ("SEP", "HETATM")})
    far = [f"HETATM 9001  N   SEP A 201      40.000  40.000  40.000  1.00 20.00           N",
           f"HETATM 9002  CA  SEP A 201      41.400  40.000  40.000  1.00 20.00           C",
           f"HETATM 9003  C   SEP A 201      42.000  41.200  40.000  1.00 20.00           C",
           f"HETATM 9004 CA    CA A 301      60.000  60.000  60.000  1.00 20.00          CA"]
    _write(tmp_path / "a.pdb", res, modres=[_mod(3, "SEP")], extra=far)
    r, rep = _adapter(tmp_path / "a.pdb", tmp_path / "o")
    assert r.returncode == 0, r.stderr
    assert [x["resid"] for x in rep["residues"]] == list(range(1, 15))                  # neither the free SEP 201 nor calcium 301 is output
    ex = rep["modified_residues"]["excluded"]
    assert [(e["resid"], e["name"]) for e in ex] == [(201, "SEP")] and "not peptide-bonded" in ex[0]["reason"] and "component" in ex[0]["reason"]
    assert [a["resid"] for a in rep["modified_residues"]["admitted"]] == [3]


def test_the_policy_refuses_a_non_pdb_topology(tmp_path):
    _write(tmp_path / "a.txt", _chain(mods={3: ("SEP", "HETATM")}), modres=[_mod(3, "SEP")])
    code = ("import sys, MDAnalysis as mda; sys.path.insert(0, %r); from topofold_inputs.ca_adapter import extract_ca, AdapterError\n"
            "u = mda.Universe(%r, topology_format='PDB')\n"
            "try:\n    extract_ca(u, 'x', altloc='blank_then_A', modified_residues='phospho_pdb')\nexcept AdapterError as e:\n    print('REFUSED', e)\n") % (str(ROOT / "python"), str(tmp_path / "a.txt"))
    r = subprocess.run([PY, "-c", code], capture_output=True, text=True)
    assert "REFUSED" in r.stdout and "PDB-format" in r.stdout, r.stdout + r.stderr


# --- structpair -------------------------------------------------------------------------------------------------------

def _spec(tmp, name, flag="flag", pairs=None, policy=True, numeric=True, sens=None):
    s = {"name": "SYNTHETIC phospho", "correspondence": {"pairs": pairs or [[i, i] for i in range(1, 15)]}, "on_missing": "exclude", "resname_mismatch": flag,
         "structures": {"A": {"label": "SYNTHETIC A", "file": "SYNTHETIC_A.pdb", "chain": "A", "altloc": "blank_then_A"},
                        "B": {"label": "SYNTHETIC B", "file": "SYNTHETIC_B.pdb", "chain": "A", "altloc": "blank_then_A"}}}
    if policy:
        s["modified_residues"] = "phospho_pdb"
    if numeric:
        s["numeric_record"] = True
    if sens:
        s["altloc_sensitivity"] = sens
    (tmp / f"{name}.json").write_text(json.dumps(s))
    return tmp / f"{name}.json"


def _cli(spec, out):
    return subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(spec), "--out", str(out)], capture_output=True, text=True,
                          env=dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=PY))


def _pair(tmp, mods_a, mods_b=None, jump_a=None):
    ma = [_mod(s, n) for s, (n, _) in mods_a.items()]
    mb = [_mod(s, n) for s, (n, _) in (mods_b or {}).items()]
    _write(tmp / "SYNTHETIC_A.pdb", _chain(mods=mods_a), modres=ma, jump_after=jump_a)
    _write(tmp / "SYNTHETIC_B.pdb", _chain(mods=mods_b or {}), modres=mb, shift=0.0)
    # B differs slightly from A so that every window has a nonzero difference
    t = (tmp / "SYNTHETIC_B.pdb").read_text().splitlines()
    out = []
    for ln in t:
        if ln[:6] in ("ATOM  ", "HETATM"):
            i = int(ln[22:26])
            x, y, z = float(ln[30:38]) + 0.05 * i, float(ln[38:46]) - 0.03 * i * (i % 3), float(ln[46:54]) + 0.02 * (i % 4)
            ln = ln[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + ln[54:]
        out.append(ln)
    (tmp / "SYNTHETIC_B.pdb").write_text("\n".join(out) + "\n")


def test_sep_versus_ser_stays_visible_and_the_modification_is_carried_to_every_output(tmp_path):
    _pair(tmp_path, {7: ("SEP", "HETATM"), 10: ("TPO", "HETATM")})
    r = _cli(_spec(tmp_path, "s_reject", flag="reject"), tmp_path / "o_rej")
    assert r.returncode == 2 and "SEP" in r.stderr and "ALA" in r.stderr                    # refused under the default: the name difference is not hidden by the parent
    r = _cli(_spec(tmp_path, "s"), tmp_path / "o")
    assert r.returncode == 0, r.stderr
    o = tmp_path / "o"
    rows = list(csv.DictReader(open(o / "residues.csv")))
    row = rows[6]
    assert row["A_resname"] == "SEP" and row["B_resname"] == "ALA" and row["resname_mismatch"] == "True"
    assert "SEP" in row["A_modification"] and "parent SER" in row["A_modification"] and "MODRES" in row["A_modification"] and row["B_modification"] == ""
    assert rows[9]["A_modification"].startswith("TPO") and "parent THR" in rows[9]["A_modification"] and rows[0]["A_modification"] == ""
    w = list(csv.DictReader(open(o / "windows.csv")))
    assert w[7]["modified_in_window3"] == "True" and w[1]["modified_in_window3"] == "False" and w[7]["mismatch_in_window3"] == "True"
    s = json.loads((o / "summary.json").read_text())["modified_residues"]
    assert s["policy"] == "phospho_pdb" and [(x["structure"], x["residue"], x["deposited_name"], x["parent"]) for x in s["retained"]] == [("A", "A:7", "SEP", "SER"), ("A", "A:10", "TPO", "THR")]
    assert "not evidence of functional state" in s["statement"] and s["excluded"] == []
    man = json.loads((o / "manifest.json").read_text())
    assert man["policies"]["modified_residues"] == "phospho_pdb" and man["software"]["structpair_version"] == sp.STRUCTPAIR_VERSION
    rec = nr.read(o)
    assert rec["positions"][6]["A_modification"].startswith("SEP") and rec["positions"][6]["B_modification"] is None and rec["positions"][0]["A_modification"] is None
    assert "Modified residues" in (o / "report.md").read_text() and "SEP" in (o / "report.md").read_text()
    assert V.verify(o)["all_ok"]


def test_default_structpair_output_is_unchanged_when_the_policy_is_off_and_standard_payload_is_identical_when_on(tmp_path):
    _pair(tmp_path, {})
    off = tmp_path / "off"
    on = tmp_path / "on"
    assert _cli(_spec(tmp_path, "s_off", policy=False), off).returncode == 0 and _cli(_spec(tmp_path, "s_on"), on).returncode == 0
    assert "A_modification" not in open(off / "residues.csv").readline() and "modified_in_window3" not in open(off / "windows.csv").readline()
    assert "modified_residues" not in json.loads((off / "summary.json").read_text())
    ro, rn = nr.read(off), nr.read(on)
    assert ro["metrics"] == rn["metrics"] and ro["global"] == rn["global"]                       # numerical payload identical
    wo, wn = list(csv.DictReader(open(off / "windows.csv"))), list(csv.DictReader(open(on / "windows.csv")))
    assert all({k: v for k, v in b.items() if k in a} == a for a, b in zip(wo, wn))               # old columns identical
    assert json.loads((on / "summary.json").read_text())["modified_residues"]["retained"] == []


def test_new_geometry_matches_independent_references_and_rigid_motion(tmp_path):
    _pair(tmp_path, {5: ("SEP", "HETATM"), 9: ("PTR", "HETATM")})
    o = tmp_path / "o"
    assert _cli(_spec(tmp_path, "s"), o).returncode == 0
    rec = nr.read(o)
    P = {k: kp.parse_pdb((tmp_path / f"SYNTHETIC_{k}.pdb").read_text()) for k in "AB"}
    ra, rb = P["A"].residues["A"], P["B"].residues["A"]
    assert [x.resname for x in ra].count("SEP") == 1 and not [x for x in ra if x.resname == "SEP"][0].read      # the independent parser alone does not read it ...
    for r in ra:
        if r.resname in ("SEP", "PTR") and r.ca is not None:
            r.read = True                                                                   # ... so the reference applies the same admission rule explicitly
    w = kp.independent_windows(ra, rb)
    for k in ("d_kappa", "d_i_i2", "rmsd3", "d_tau", "d_i_i3", "rmsd4"):
        got = rec["metrics"][k]["values"]
        assert [v is None for v in got] == [v is None for v in w[k]], k
        assert max(abs(a - b) for a, b in zip(got, w[k]) if a is not None) < 1e-4, k
    assert abs(rec["global"]["rmsd_qcp_angstrom"] - w["global_rmsd"]) < 1e-4 and rec["global"]["n_ca"] == 14
    # rigid motion of both structures by exact coordinate permutations and translations
    for k, R, t in (("A", np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], float), np.array([12.5, -7.25, 3.0])),
                    ("B", np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], float), np.array([-20.0, 5.0, 40.0]))):
        lines = []
        for ln in (tmp_path / f"SYNTHETIC_{k}.pdb").read_text().splitlines():
            if ln[:6] in ("ATOM  ", "HETATM"):
                v = R @ np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]) + t
                ln = ln[:30] + "".join(f"{x:8.3f}" for x in v) + ln[54:]
            lines.append(ln)
        (tmp_path / f"MOVED_{k}.pdb").write_text("\n".join(lines) + "\n")
    s = json.loads(_spec(tmp_path, "s").read_text())
    s["structures"]["A"]["file"], s["structures"]["B"]["file"] = "MOVED_A.pdb", "MOVED_B.pdb"
    (tmp_path / "s_moved.json").write_text(json.dumps(s))
    assert _cli(tmp_path / "s_moved.json", tmp_path / "om").returncode == 0
    rm = nr.read(tmp_path / "om")
    dev = max(abs(a - b) for k in rec["metrics"] for a, b in zip(rec["metrics"][k]["values"], rm["metrics"][k]["values"]) if a is not None and b is not None)
    assert dev < 1e-4 and abs(rm["global"]["rmsd_qcp_angstrom"] - rec["global"]["rmsd_qcp_angstrom"]) < 1e-4


def test_excluded_cases_are_reported_in_the_summary_and_the_exclusion_windows_stay(tmp_path):
    mods = {4: ("SEP", "HETATM"), 8: ("TPO", "HETATM")}
    _write(tmp_path / "SYNTHETIC_A.pdb", _chain(mods=mods, over={4: {"atoms": ("N", "C", "O")}}), modres=[_mod(4, "SEP"), _mod(8, "TPO")])
    _write(tmp_path / "SYNTHETIC_B.pdb", _chain(), shift=0.05)
    o = tmp_path / "o"
    r = _cli(_spec(tmp_path, "s"), o)
    assert r.returncode == 0, r.stderr
    s = json.loads((o / "summary.json").read_text())["modified_residues"]
    assert [x["residue"] for x in s["retained"]] == ["A:8"] and [(x["structure"], x["residue"], x["name"]) for x in s["excluded"]] == [("A", "A:4", "SEP")]
    assert "N, CA and C" in s["excluded"][0]["reason"]
    res = list(csv.DictReader(open(o / "residues.csv")))
    assert res[3]["A_status"] == "residue_not_found" and res[3]["A_resname"] == ""            # unobserved for the workflow: the position is a missing row
    ex = list(csv.DictReader(open(o / "exclusions.csv")))
    # position 3 is missing: windows centred on 2, 3, 4 need it; the run before it has only 3 positions (< 4), so centre 1 has no value either (existing rules)
    assert {int(x["position"]) for x in ex if x["window_size"] == "3"} == {1, 2, 3, 4}


def test_viewer_shows_retained_modifications_and_excluded_cases(tmp_path):
    pytest.importorskip("streamlit")
    from streamlit.testing.v1 import AppTest
    mods = {4: ("SEP", "HETATM"), 8: ("TPO", "HETATM")}
    _write(tmp_path / "SYNTHETIC_A.pdb", _chain(mods=mods, over={4: {"atoms": ("N", "C", "O")}}), modres=[_mod(4, "SEP"), _mod(8, "TPO")])
    _write(tmp_path / "SYNTHETIC_B.pdb", _chain(), shift=0.05)
    res = tmp_path / "results"
    assert _cli(_spec(tmp_path, "s"), res / "phospho").returncode == 0
    at = AppTest.from_file(str(ROOT / "apps/structpair_viewer/app.py"), default_timeout=90)
    at.session_state["results_root"] = str(res)
    at.run()
    assert not at.exception
    text = " ".join(m.value for m in list(at.info) + list(at.warning) + list(at.markdown))
    assert "TPO" in text and "retained" in text.lower() and "A:4" in text and "excluded" in text.lower() and "not renamed" in text
