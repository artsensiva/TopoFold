"""Tests for topofold_inputs.mdcath_bridge on a SYNTHETIC HDF5 fixture in the mdCATH layout (plus an optional
mapping check on the local development file). Needs h5py and MDAnalysis, e.g.:
    PYTHONPATH=python benchmarks/audit/wpd_local_retrieval/build/venv/bin/python -m pytest tests/test_mdcath_bridge.py
"""

import sys
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("MDAnalysis")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from topofold_inputs import AdapterError  # noqa: E402
from topofold_inputs.mdcath_bridge import extract_ca_mdcath, mdcath_universe  # noqa: E402

RES = ["GLY", "ALA", "SER", "LEU", "VAL"]


def pdb_text(n_res):
    lines, k = ["MODEL        1"], 1
    for r in range(n_res):
        for name in ("N", "CA", "C"):
            lines.append(f"ATOM  {k:5d}  {name:<3s} {RES[r % 5]:>3s} 0{r + 1:4d}    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}"
                         f"  0.00  0.00      P0   {name[0]}")
            k += 1
    return "\n".join(lines + ["ENDMDL", "END"]) + "\n"


def make_h5(path, n_res=6, n_frames=12, unit="Angstrom", resid_shift=0, extra_atom=False, box_nm=4.0):
    """SYNTHETIC mdCATH-like file: one domain, temperatures 320/348, replicas 0/1."""
    n_atoms = 3 * n_res
    rng = np.random.default_rng(0)
    with h5py.File(path, "w") as f:
        f.attrs["layout"] = "mdcath-only-protein-v1.0"
        g = f.create_group("SYN0A00")
        g.create_dataset("pdbProteinAtoms", data=pdb_text(n_res))
        resid = np.repeat(np.arange(1, n_res + 1), 3) + resid_shift
        g.create_dataset("resid", data=resid)
        g.create_dataset("resname", data=np.array([RES[r % 5] for r in range(n_res) for _ in range(3)], dtype="S3"))
        g.create_dataset("chain", data=np.array(["0"] * n_atoms, dtype="S1"))
        for T in ("320", "348"):
            for rep in ("0", "1"):
                grp = g.create_group(f"{T}/{rep}")
                X = np.zeros((n_frames, n_atoms + int(extra_atom), 3), np.float32)
                for fr in range(n_frames):
                    for r in range(n_res):
                        base = np.array([3.8 * r, 0.6 * (-1) ** r, 0.01 * fr + 0.1 * int(rep)])
                        X[fr, 3 * r:3 * r + 3] = base + rng.normal(scale=0.05, size=(3, 3))
                ds = grp.create_dataset("coords", data=X)
                ds.attrs["unit"] = unit
                grp.attrs["numFrames"] = n_frames
                grp.create_dataset("box", data=np.diag([box_nm] * 3).astype(np.float32))
    return path


def test_mapping_units_frames_and_provenance(tmp_path):
    p = make_h5(tmp_path / "SYNTHETIC.h5")
    res = extract_ca_mdcath(p, "SYN0A00", "348", "1", "SYNTHETIC_bridge", start=3, stop=11, step=2)
    with h5py.File(p) as f:
        X = f["SYN0A00/348/1/coords"][()]
    ca_cols = [3 * r + 1 for r in range(6)]                 # independent: CA is the 2nd atom of each residue
    np.testing.assert_array_equal(res.coordinates, X[[3, 5, 7, 9]][:, ca_cols].astype(np.float64))
    assert res.frame_index.tolist() == [3, 5, 7, 9]
    assert [r["resid"] for r in res.residues] == [1, 2, 3, 4, 5, 6]
    assert [r["resname"] for r in res.residues] == RES + ["GLY"]
    rep = res.report
    assert rep["frames_inferred_time"]["inferred_time_ns"] == [3.0, 5.0, 7.0, 9.0]
    assert "INFERRED" in rep["frames_inferred_time"]["basis"] and np.isnan(res.time_ps).all()
    assert rep["frame_identity"]["status"] == "not_verifiable"
    assert rep["provenance"]["bridge"]["units"]["box_angstrom"] == [40.0, 40.0, 40.0]
    assert rep["provenance"]["source_files"][0]["role"] == "hdf5"
    assert rep["provenance"]["bridge"]["dataset"]["coords"] == "/SYN0A00/348/1/coords"


def test_rejects_wrong_units_and_mapping(tmp_path):
    with pytest.raises(AdapterError, match="unit"):
        mdcath_universe(make_h5(tmp_path / "a.h5", unit="nm"), "SYN0A00", "320", "0")
    with pytest.raises(AdapterError, match="mapping differs"):
        mdcath_universe(make_h5(tmp_path / "b.h5", resid_shift=1), "SYN0A00", "320", "0")
    with pytest.raises(AdapterError, match="atoms but coords"):
        mdcath_universe(make_h5(tmp_path / "c.h5", extra_atom=True), "SYN0A00", "320", "0")
    with pytest.raises(AdapterError, match="no frames"):
        mdcath_universe(make_h5(tmp_path / "d.h5"), "SYN0A00", "320", "0", start=50)


def test_pbc_uses_converted_box(tmp_path):
    p = make_h5(tmp_path / "SYNTHETIC.h5", box_nm=1.0)       # 10 Å box: the 3.8 Å steps stay below half a box
    a = extract_ca_mdcath(p, "SYN0A00", "320", "0", "SYNTHETIC_pbc", pbc="make_whole_sequential")
    assert a.report["pbc"]["box_present_in_any_frame"] and a.report["pbc"]["policy"] == "make_whole_sequential"


def test_real_development_file_mapping():
    p = ROOT / "benchmarks/audit/wpd_local_retrieval/data/mdcath/mdcath_dataset_1ha8A00.h5"
    if not p.is_file():
        pytest.skip("local development file not present")
    res = extract_ca_mdcath(p, "1ha8A00", "320", "0", "map_check", start=43, stop=46)
    with h5py.File(p) as f:
        g = f["1ha8A00"]
        X = g["320/0/coords"][43:46]
        resid = g["resid"][()]
        text = g["pdbProteinAtoms"][()].decode()
    # independent Cα column lookup from the PDB text in the file (ATOM records in coords order)
    atoms = [ln for ln in text.splitlines() if ln.startswith("ATOM")]
    ca = [i for i, ln in enumerate(atoms) if ln[12:16].strip() == "CA"]
    assert len(ca) == 51 and res.coordinates.shape == (3, 51, 3)
    np.testing.assert_array_equal(res.coordinates, X[:, ca].astype(np.float64))
    assert [r["resid"] for r in res.residues] == sorted(set(resid.tolist()))
    assert res.breaks == []
