"""SYNTHETIC tests of amendment a2 of the ABL-MYR-1b check (independent GROMACS reading of md.tpr). Needs MDAnalysis; the
real-GROMACS tests need TOPOFOLD_GMX (they generate their own water box with gmx and read nothing from Dryad)."""

import json
import os
import stat
import subprocess
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
import run_comp_check_a2 as A2  # noqa: E402
import test_abl_myr_1b_check as T  # noqa: E402

GMX = os.environ.get("TOPOFOLD_GMX")
needs_gmx = pytest.mark.skipif(not GMX, reason="TOPOFOLD_GMX not set")
SPEC = json.loads((ROOT / "benchmarks/audit/abl_myr_1b/check_amendment_a2.json").read_text())

DUMP = """inputrec:
   integrator                     = md
   dt                             = 0.002
   nsteps                         = 250000000
   ld-seed                        = -1234
   tcoupl                         = V-rescale
   nsttcouple                     = 10
   pcoupl                         = Parrinello-Rahman
   continuation                   = true
   constraints                    = h-bonds
grpopts:
   nrdf:       13911        90000
   ref-t:         300          300
   tau-t:           1            1
header:
   bIr    = present
   bBox   = present
   bTop   = present
   bX     = present
   bV     = present
   bF     = not present
   natoms = 3
topology:
   #atoms                         = 3
   molblock (0):
      moltype              = 0 "Protein_chain_A"
      #molecules           = 1
   molblock (1):
      moltype              = 1 "SOL"
      #molecules           = 7
box (3x3):
   box[    0]={ 4.00000e+00,  0.00000e+00,  0.00000e+00}
   box[    1]={ 0.00000e+00,  4.00000e+00,  0.00000e+00}
   box[    2]={ 2.00000e+00,  2.00000e+00,  2.82843e+00}
box_rel (3x3):
   box_rel[    0]={ 0.00000e+00,  0.00000e+00,  0.00000e+00}
x (3x3):
   x[    0]={ 1.00000e+00,  2.00000e+00, -3.00000e-01}
   x[    1]={ 1.10000e+00,  2.00000e+00, -3.00000e-01}
   x[    2]={-1.23456e-02,  3.14159e+00,  4.00000e+00}
v (3x3):
   v[    0]={ 5.52213e-01, -1.97867e+00,  3.63006e-01}
   v[    1]={ 0.00000e+00,  0.00000e+00,  0.00000e+00}
   v[    2]={ 1.00000e+00,  2.00000e+00,  3.00000e+00}
"""


def test_parse_dump_reads_box_coordinates_velocities_and_settings():
    d = A2.parse_dump(DUMP.splitlines())
    assert d["natoms"] == 3 and d["flags"] == {"bBox": True, "bX": True, "bV": True, "bF": False}
    assert d["box"].tolist() == [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [2.0, 2.0, 2.82843]]
    assert d["x"][2].tolist() == [-0.0123456, 3.14159, 4.0] and d["v"].shape == (3, 3)
    assert d["ir"]["tcoupl"] == "V-rescale" and d["ir"]["ld-seed"] == "-1234" and d["ir"]["continuation"] == "true"
    assert (d["ref_t"], d["tau_t"], d["nrdf"]) == ("300          300", "1            1", "13911        90000")
    assert d["molblocks"] == [{"moltype": 0, "name": "Protein_chain_A", "molecules": 1},
                              {"moltype": 1, "name": "SOL", "molecules": 7}]


@pytest.mark.parametrize("mutate,msg", [
    (lambda t: t.replace("x[    1]=", "x[    2]="), "out of order"),
    (lambda t: t.replace("   x[    2]={-1.23456e-02,  3.14159e+00,  4.00000e+00}\n", ""), "truncated x"),
    (lambda t: t.replace("box (3x3):", "box (3x4):"), "not N x 3"),
    (lambda t: t.replace("v (3x3):\n", "").replace("   v[    0]={ 5.52213e-01, -1.97867e+00,  3.63006e-01}\n", "")
     .replace("   v[    1]={ 0.00000e+00,  0.00000e+00,  0.00000e+00}\n", "").replace("   v[    2]={ 1.00000e+00,  2.00000e+00,  3.00000e+00}\n", "")
     .replace("box (3x3):", "nothing"), "no box"),
])
def test_parse_dump_refuses_malformed_input(mutate, msg):
    with pytest.raises(ValueError, match=msg):
        A2.parse_dump(mutate(DUMP).splitlines())


def test_parse_dump_without_velocities():
    t = DUMP.split("v (3x3):")[0].replace("bV     = present", "bV     = not present")
    d = A2.parse_dump(t.splitlines())
    assert d["v"] is None and d["flags"]["bV"] is False


def test_box_validity_and_conversion_to_dimensions():
    ok, _ = A2.box_valid([[4.0, 0, 0], [0, 4.0, 0], [2.0, 2.0, 2.82843]])
    assert ok
    assert not A2.box_valid([[4.0, 0.1, 0], [0, 4.0, 0], [0, 0, 4.0]])[0]
    assert not A2.box_valid([[4.0, 0, 0], [0, 0.0, 0], [0, 0, 4.0]])[0]
    assert not A2.box_valid([[np.nan, 0, 0], [0, 4, 0], [0, 0, 4.0]])[0]
    dims = A2.box_dims_A([[4.65918, 0, 0], [0, 4.65918, 0], [2.32959, 2.32959, 3.29454]])        # rhombic dodecahedron
    assert dims[:3] == pytest.approx([46.5918] * 3, abs=1e-3) and dims[3:] == pytest.approx([60.0, 60.0, 90.0], abs=1e-3)


def rec_and_dump(n=5, seed=0):
    rng = np.random.default_rng(seed)
    x = rng.uniform(0, 4, size=(n, 3))
    v = rng.normal(size=(n, 3))
    rec = {"n_atoms": n, "positions": x.copy(), "velocities": v.copy()}                         # the earlier reader: nm, nm/ps
    g = {"natoms": n, "x": x.copy(), "v": v.copy(), "box": np.array([[4.0, 0, 0], [0, 4.0, 0], [2.0, 2.0, 2.8]]),
         "flags": {"bBox": True, "bX": True, "bV": True}}
    return rec, g


def test_agreement_tolerances_and_failures():
    rec, g = rec_and_dump()
    ok = A2.agreement(rec, g, g["box"])
    assert ok["ok"] and ok["max_abs_dx_nm"] == 0.0 and ok["box_vs_mdanalysis_extract_box_info_nm"] == 0.0
    g["x"] = g["x"] + 0.0009
    assert A2.agreement(rec, g, g["box"])["checks"]["positions_agree"]                         # inside 0.001 nm
    g["x"] = g["x"] + 0.0005
    bad = A2.agreement(rec, g, g["box"])
    assert not bad["ok"] and "positions_agree" in bad["problems"]
    rec, g = rec_and_dump()
    g["v"] = g["v"] + 0.002
    assert "velocities_agree" in A2.agreement(rec, g, g["box"])["problems"]
    rec, g = rec_and_dump()
    g["v"] = None
    assert "velocities_agree" in A2.agreement(rec, g, g["box"])["problems"]                    # present in one reader only
    rec["velocities"] = None
    assert A2.agreement(rec, g, g["box"])["checks"]["velocities_agree"]                        # absent in both: fine
    rec, g = rec_and_dump()
    assert A2.agreement(rec, g, None)["ok"] and A2.agreement(rec, g, None)["box_vs_mdanalysis_extract_box_info_nm"] is None
    assert "box_agrees_with_extract_box_info" in A2.agreement(rec, g, g["box"] + 0.001)["problems"]
    g["box"] = np.array([[4.0, 1.0, 0], [0, 4.0, 0], [0, 0, 4.0]])
    assert not A2.agreement(rec, g, None)["ok"]
    rec, g = rec_and_dump()
    g["natoms"] = 6
    r = A2.agreement(rec, g, None)
    assert not r["ok"] and "natoms_equal" in r["problems"] and "max_abs_dx_nm" not in r
    rec, g = rec_and_dump()
    g["flags"]["bBox"] = False
    assert "header_box_and_x_present" in A2.agreement(rec, g, None)["problems"]
    rec, g = rec_and_dump()
    rec["positions"] = rec["positions"] * 10.0                                                  # a wrong unit assumption fails
    assert not A2.agreement(rec, g, None)["checks"]["positions_agree"]


def test_unwrap_by_bonds_makes_a_wrapped_molecule_whole_and_picks_the_main_component():
    dims = np.array([20.0, 20.0, 20.0, 90.0, 90.0, 90.0])
    chain = np.array([[18.0 + 1.5 * i, 5.0, 5.0] for i in range(6)])
    pos = np.vstack([chain, [[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]]])                             # + a second, 2-atom molecule
    pos[:6, 0] %= 20.0
    bonds = [(i, i + 1) for i in range(5)] + [(6, 7)]
    xyz, comps = A2.unwrap_by_bonds(pos, dims, bonds, heavy_idx=[0, 1, 2, 3, 4, 5, 6, 7])
    assert comps == 2 and len(xyz) == 6
    assert np.allclose(np.linalg.norm(np.diff(xyz, axis=0), axis=1), 1.5, atol=1e-4)
    tri = np.array([50.0, 50.0, 50.0, 60.0, 60.0, 90.0])
    w, _ = A2.unwrap_by_bonds(np.array([[1.0, 1.0, 1.0], [49.0, 1.0, 1.0]]), tri, [(0, 1)], [0, 1])
    assert np.linalg.norm(w[1] - w[0]) == pytest.approx(2.0, abs=1e-4)
    none, comps = A2.unwrap_by_bonds(pos, dims, [], [0, 1])
    assert none.shape == (0, 3) and comps == 0


def fake_gmx(tmp_path, dump_text, rc=0, version_ok=True):
    prefix = tmp_path / "env"
    (prefix / "conda-meta").mkdir(parents=True)
    (prefix / "conda-meta/gromacs-2024.5-x.json").write_text(json.dumps(
        {"name": "gromacs", "version": "2024.5", "build": "nompi_x", "channel": "conda-forge", "sha256": "ab" * 32,
         "url": f"{prefix}/pkg.conda"}))
    exe = prefix / "bin.AVX2_256/gmx"
    exe.parent.mkdir()
    exe.write_bytes(b"binary")
    dump = tmp_path / "dump.txt"
    dump.write_text(dump_text)
    sh = prefix / "bin/gmx"
    sh.parent.mkdir()
    sh.write_text(f"""#!/bin/sh
if [ "$1" = "--version" ]; then
  {'echo "GROMACS version:     2024.5-conda_forge"; echo "Precision:           mixed"; echo "MPI library:         thread_mpi"; echo "Executable:   ' + str(exe) + '"; echo "Data prefix:  ' + str(prefix) + '"' if version_ok else 'exit 3'}
  exit 0
fi
if [ "$1" = "dump" ]; then
  if [ {rc} -ne 0 ]; then echo "boom" >&2; exit {rc}; fi
  cat {dump}
  exit 0
fi
exit 9
""")
    sh.chmod(sh.stat().st_mode | stat.S_IXUSR)
    return str(sh)


def test_gmx_executable_discovery(tmp_path, monkeypatch):
    assert A2.gmx_executable({"TOPOFOLD_GMX": str(tmp_path / "nothing")}) is None
    g = fake_gmx(tmp_path, DUMP)
    assert A2.gmx_executable({"TOPOFOLD_GMX": g}) == g
    monkeypatch.setenv("PATH", "")
    assert A2.gmx_executable({}) is None


def test_provenance_and_dump_reader_with_a_fake_gmx(tmp_path):
    g = fake_gmx(tmp_path, DUMP)
    p = A2.gmx_provenance(g)
    assert p["version"] == "2024.5-conda_forge" and p["precision"] == "mixed" and p["package"]["build"] == "nompi_x"
    assert str(tmp_path) not in json.dumps(p) and "<env>" in p["package"]["url"]
    assert len(p["executable_sha256"]) == 64 and p["launcher_sha256"] != p["executable_sha256"]
    d = A2.read_gmx(tmp_path / "x.tpr", g)
    assert d["natoms"] == 3 and d["returncode"] == 0
    with pytest.raises(ValueError, match="exit code 5"):
        A2.read_gmx(tmp_path / "x.tpr", fake_gmx(tmp_path / "b", DUMP, rc=5))
    with pytest.raises(A2.GmxUnavailable):
        A2.gmx_provenance(fake_gmx(tmp_path / "c", DUMP, version_ok=False))


def make_loaders(**kw):
    base = T.loader_factory()

    def gmx_reader(path):
        rec = base(path)
        n = rec["n_atoms"]
        pos_nm = np.asarray(rec["positions"]) / 10.0
        v = None if rec["velocities"] is None else np.asarray(rec["velocities"]) / 10.0
        g = {"natoms": n, "x": pos_nm + kw.get("dx", {}).get(Path(path).read_bytes().decode().split()[3], 0.0), "v": v,
             "box": np.array([[T.BOX_NM, 0, 0], [0, T.BOX_NM, 0], [0, 0, T.BOX_NM]]),
             "flags": {"bBox": True, "bX": True, "bV": v is not None}, "ir": {"tcoupl": "V-rescale"}, "ref_t": "300 300",
             "tau_t": "1 1", "nrdf": "1 1", "molblocks": [{"name": "SOL", "molecules": 3}]}
        return g

    def mda_loader(path):                                                                   # the earlier reader: nm, nm/ps
        rec = dict(base(path))
        rec["positions"] = np.asarray(rec["positions"]) / 10.0
        rec["velocities"] = None if rec["velocities"] is None else np.asarray(rec["velocities"]) / 10.0
        rec["dims"] = None
        return rec
    return mda_loader, gmx_reader


def test_loader_replaces_reader_values_with_gromacs_values_in_angstrom():
    A2.SIDE.clear()
    mda_loader, gmx_reader = make_loaders()
    p = "/x/kinase_kinked_Myr/kinase_kinked_Myr/Set1/11/md.tpr"
    class Fake:
        pass
    import tempfile
    d = Path(tempfile.mkdtemp())
    f = d / "kinase_kinked_Myr/kinase_kinked_Myr/Set1/11/md.tpr"
    f.parent.mkdir(parents=True)
    f.write_bytes(b"SYNTHETIC Myr Set1 11")
    out = A2.load_tpr_gmx(f, "gmx", mda_loader=mda_loader, gmx_reader=gmx_reader, box_reader=lambda p: gmx_reader(p)["box"],
                          heavy_fn=lambda rec, pos, dims, idx: pos[idx])
    assert out["dims"].tolist() == [T.BOX_NM * 10, T.BOX_NM * 10, T.BOX_NM * 10, 90.0, 90.0, 90.0]
    assert out["positions"].shape == out["velocities"].shape and len(out["ca_positions"]) == 126
    assert list(A2.SIDE) == ["Set1:Myr:11"] and A2.SIDE["Set1:Myr:11"]["agreement"]["ok"]
    assert out["heavy_whole"].shape[1] == 3 and p
    f.write_bytes(b"SYNTHETIC Myr Set1 11")
    mda_loader2, gmx_reader2 = make_loaders(dx={"11": 0.01})
    with pytest.raises(ValueError, match="disagree.*positions_agree"):
        A2.load_tpr_gmx(f, "gmx", mda_loader=mda_loader2, gmx_reader=gmx_reader2, box_reader=lambda p: None,
                        heavy_fn=lambda rec, pos, dims, idx: pos[idx])
    assert A2.SIDE["Set1:Myr:11"]["agreement"]["ok"] is False                               # recorded even when refused


def run_a2(tmp_path, dx=None):
    orig, data = T.build(tmp_path)
    mda_loader, gmx_reader = make_loaders(dx=dx or {})
    loader = lambda p: A2.load_tpr_gmx(p, "gmx", mda_loader=mda_loader, gmx_reader=gmx_reader,  # noqa: E731
                                       box_reader=lambda q: gmx_reader(q)["box"], heavy_fn=lambda rec, pos, dims, idx: pos[idx])
    return A2.run_check(orig, data, tmp_path / "ext", T.BUDGET, "gmx", loader=loader, log=lambda *_: None)[0]


def test_run_check_passes_and_reports_geometry_separately_by_condition(tmp_path):
    r = run_a2(tmp_path)
    assert r["outcome"] == "comparable_for_1b", r.get("decision")
    a = r["a2"]
    assert a["runs"] == 40 and a["agreement_ok"] == 40 and a["velocities_in_both"] == 40
    g = r["R2_production_start_geometry"]["Set1"]
    assert {"within_Myr", "within_noMyr", "between_Myr_noMyr"} <= set(g) and g["between_Myr_noMyr"]["ca_rmsd_A"]["n"] == 100
    assert set(a["run_input_settings"]) == {f"{s}:{c}" for s in B.SETS for c in B.CONDITIONS}
    assert a["myr_vs_nomyr_run_input_differences_descriptive"] == {"Set1": {}, "Set2": {}}
    assert "Amendment a2" in A2.render_report(r) and "independent reader" in A2.render_report(r)


def test_run_check_k1_fails_when_one_run_disagrees(tmp_path):
    r = run_a2(tmp_path, dx={"13": 0.05})
    assert r["outcome"] == "not_comparable" and r["decision"]["failed"] == ["K1_usable_coordinates"]
    bad = [k for k, v in r["R1_usable_coordinates"]["per_run"].items() if not v["ok"]]
    assert sorted(bad) == ["Set1:Myr:13", "Set1:noMyr:13"] and r["a2"]["agreement_ok"] == 38


def test_side_key_and_stop_s5_message(tmp_path):
    assert A2.side_key("/e/kinase_kinked_noMyr/kinase_kinked_noMyr/Set2/25/md.tpr") == "Set2:noMyr:25"


def test_a2_freeze_preserves_earlier_records():
    frz = ROOT / "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_a2_freeze.json"
    if not frz.exists():
        pytest.skip("a2 freeze not yet written")
    f = json.loads(frz.read_text())
    for p, h in f["a1_records_sha256"].items():
        assert B.R0.sha256_file(ROOT / p) == h, p
    assert B.R0.sha256_file(B.SELF) == f["frozen_check"]["code_sha256"]
    assert f["gromacs_at_freeze"]["version"] == "2024.5-conda_forge"


def test_spec_constants_match_the_code():
    t = SPEC["reader_agreement_added_to_K1"]["tolerances"]
    assert (A2.POS_TOL_NM, A2.VEL_TOL_NM_PS, A2.BOX_TOL_NM) == (t["positions_nm"], t["velocities_nm_per_ps"], t["box_nm"])
    assert SPEC["budget"]["cpu_seconds_max"] == 1800 and "S5" in SPEC["stop_conditions_added"]
    assert SPEC["single_attempt"].startswith("This is the only further attempt")


# --- real GROMACS on a SYNTHETIC water box generated by the test ----------------------------------------------------------

@pytest.fixture(scope="module")
def synth_tpr(tmp_path_factory):
    if not GMX:
        pytest.skip("TOPOFOLD_GMX not set")
    d = tmp_path_factory.mktemp("synth")
    top = Path(GMX).resolve().parents[1] / "share/gromacs/top"
    run = lambda *a: subprocess.run([GMX, *a], cwd=d, capture_output=True, text=True, check=True)  # noqa: E731
    run("editconf", "-f", str(top / "spc216.gro"), "-bt", "dodecahedron", "-d", "0.8", "-o", "box.gro")
    last = (d / "box.gro").read_text().splitlines()[-1]
    (d / "empty.gro").write_text(f"SYNTHETIC empty box\n0\n{last}\n")
    (d / "topol.top").write_text('#include "oplsaa.ff/forcefield.itp"\n#include "oplsaa.ff/spc.itp"\n[ system ]\nSYNTHETIC water\n[ molecules ]\n')
    run("solvate", "-cp", "empty.gro", "-cs", str(top / "spc216.gro"), "-o", "w.gro", "-p", "topol.top")
    (d / "md.mdp").write_text("integrator = md\nnsteps = 0\ngen_vel = yes\ngen_temp = 300\ngen_seed = 17\nconstraints = h-bonds\n")
    run("grompp", "-f", "md.mdp", "-c", "w.gro", "-p", "topol.top", "-o", "synth.tpr", "-maxwarn", "2")
    return d


@needs_gmx
def test_real_gmx_dump_matches_the_input_structure_and_box(synth_tpr):
    g = A2.read_gmx(synth_tpr / "synth.tpr", GMX)
    lines = (synth_tpr / "w.gro").read_text().splitlines()
    x = np.array([[float(ln[20:28]), float(ln[28:36]), float(ln[36:44])] for ln in lines[2:-1]])
    assert g["natoms"] == len(x) and np.abs(g["x"] - x).max() < 1e-9 and g["v"] is not None and g["flags"]["bV"]
    dims = A2.box_dims_A(g["box"])
    assert dims == pytest.approx(B.gro_box_dims(synth_tpr / "w.gro"), abs=1e-3)               # box from gro text vs from the tpr
    assert g["ir"]["integrator"] == "md" and g["molblocks"][0]["name"] == "SOL"
    prov = A2.gmx_provenance(GMX)
    assert prov["version"].startswith("2024.5") and prov["package"]["name"] == "gromacs"


@needs_gmx
def test_real_gmx_and_mdanalysis_agree_on_the_synthetic_run_input(synth_tpr):
    A2.SIDE.clear()
    out = A2.load_tpr_gmx(synth_tpr / "synth.tpr", GMX)
    ag = list(A2.SIDE.values())[0]["agreement"]
    assert ag["ok"] and ag["max_abs_dx_nm"] < 1e-6 and ag["max_abs_dv_nm_ps"] < 1e-4
    assert ag["box_vs_mdanalysis_extract_box_info_nm"] < 1e-6
    assert out["dims"][:3] == pytest.approx([46.592] * 3, abs=1e-3) and out["dims"][3:] == pytest.approx([60, 60, 90], abs=1e-3)
