"""Engineering tests for topofold_inputs.ca_adapter. All structures here are SYNTHETIC fixtures.

Requires MDAnalysis (skipped otherwise; the core ``topofold`` package does not depend on it). Run with an
environment that has it, e.g.:
    PYTHONPATH=python benchmarks/audit/wpd_local_retrieval/build/venv/bin/python -m pytest tests/test_ca_adapter.py
"""

import hashlib
import json
import random
import sys
from pathlib import Path

import numpy as np
import pytest

mda = pytest.importorskip("MDAnalysis")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from topofold_inputs import AdapterError, ResidueSpec, extract_ca, load_universe  # noqa: E402


# --- SYNTHETIC fixture helpers ---------------------------------------------------------------------------------------

def pdb_line(serial, name, alt, resname, chain, resid, icode, xyz, occ=1.0):
    x, y, z = xyz
    return (f"ATOM  {serial:5d} {(' ' + name):<4s}{alt or ' ':1s}{resname:>3s} {chain:1s}{resid:4d}{icode or ' ':1s}"
            f"   {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{0.0:6.2f}          {name[0]:>2s}")


def write_pdb(path, models):
    """models: list of lists of (name, alt, resname, chain, resid, icode, xyz[, occ])."""
    lines = ["REMARK   1 SYNTHETIC test fixture"]
    for m, atoms in enumerate(models):
        lines.append(f"MODEL     {m + 1:4d}")
        lines += [pdb_line(i + 1, *a) for i, a in enumerate(atoms)]
        lines += ["TER", "ENDMDL"]
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n")
    return path


def zigzag(k, frame=0, dy=0.0):
    return (3.8 * k, 0.5 * (-1) ** k + dy + 0.1 * frame * ((k % 3) - 1), 0.05 * frame)


def chain_atoms(chain, resids, frame=0, dy=0.0, start_k=0):
    return [("CA", "", "ALA", chain, r, ic, zigzag(start_k + i, frame, dy)) for i, (r, ic) in enumerate(resids)]


def load(tmp_path, models, name="SYNTHETIC.pdb"):
    return load_universe(write_pdb(tmp_path / name, models))


# --- tests ---------------------------------------------------------------------------------------------------------

def test_duplicate_resids_across_chains_keep_distinct_identity(tmp_path):
    res = [(1, ""), (2, ""), (3, "")]
    u = load(tmp_path, [chain_atoms("A", res) + chain_atoms("B", res, dy=20.0)])
    out = extract_ca(u, "SYNTHETIC_dup")
    assert [(r["chain"], r["resid"]) for r in out.residues] == [("A", 1), ("A", 2), ("A", 3),
                                                               ("B", 1), ("B", 2), ("B", 3)]
    assert {"after_index": 2, "reasons": ["chain_boundary"]} in out.breaks
    assert out.segments == [[0, 3], [3, 6]]
    pick = extract_ca(u, "SYNTHETIC_dup", correspondence=[ResidueSpec("B", 2, label="b2"),
                                                          ResidueSpec("A", 2, label="a2")])
    np.testing.assert_array_equal(pick.coordinates[0, 0], np.array(zigzag(1, dy=20.0), dtype=np.float32))
    np.testing.assert_array_equal(pick.coordinates[0, 1], np.array(zigzag(1), dtype=np.float32))
    assert [r["label"] for r in pick.residues] == ["b2", "a2"]


def test_insertion_codes_are_distinct_and_contiguous(tmp_path):
    u = load(tmp_path, [chain_atoms("A", [(51, ""), (52, ""), (52, "A"), (52, "B"), (53, "")])])
    out = extract_ca(u, "SYNTHETIC_icode")
    assert [(r["resid"], r["icode"]) for r in out.residues] == [(51, ""), (52, ""), (52, "A"), (52, "B"), (53, "")]
    assert out.breaks == [] and out.segments == [[0, 5]]
    one = extract_ca(u, "SYNTHETIC_icode", correspondence=[{"chain": "A", "resid": 52, "icode": "A"}])
    np.testing.assert_array_equal(one.coordinates[0, 0], np.array(zigzag(2), dtype=np.float32))


def test_missing_residues_and_atoms_are_reported_never_joined(tmp_path):
    atoms = chain_atoms("A", [(1, ""), (2, ""), (3, "")])
    atoms += [("N", "", "ALA", "A", 4, "", zigzag(3))]                     # residue 4 has no Cα
    atoms += chain_atoms("A", [(5, ""), (9, ""), (10, "")], start_k=4)      # numbering gap 5 -> 9
    u = load(tmp_path, [atoms])
    out = extract_ca(u, "SYNTHETIC_gaps", max_ca_ca_angstrom=10.0)
    st = {r["resid"]: r["status"] for r in out.residues}
    assert st[4] == "no_CA_atom" and not out.present[3] and np.isnan(out.coordinates[0, 3]).all()
    reasons = {b["after_index"]: b["reasons"] for b in out.breaks}
    assert reasons[2] == ["missing_CA"] and reasons[3] == ["missing_CA"]
    assert reasons[4] == ["numbering_discontinuity"]
    assert out.segments == [[0, 3], [4, 5], [5, 7]]
    marked = extract_ca(u, "SYNTHETIC_gaps", correspondence=[{"chain": "A", "resid": 7}])
    assert marked.residues[0]["status"] == "residue_not_found" and not marked.present[0]
    with pytest.raises(AdapterError, match="not found"):
        extract_ca(u, "SYNTHETIC_gaps", correspondence=[{"chain": "A", "resid": 7}], on_missing="reject")


def test_geometric_break_without_numbering_gap(tmp_path):
    atoms = chain_atoms("A", [(1, ""), (2, "")]) + [("CA", "", "ALA", "A", 3, "", (30.0, 0.0, 0.0))]
    out = extract_ca(load(tmp_path, [atoms]), "SYNTHETIC_geo")
    assert out.breaks[0]["after_index"] == 1 and out.breaks[0]["reasons"] == ["ca_ca_distance"]
    assert out.segments == [[0, 2], [2, 3]]


def test_reordered_atoms_give_identical_output(tmp_path):
    resids = [(i, "") for i in range(1, 7)]
    ordered = []
    for i, (r, _) in enumerate(resids):
        x, y, z = zigzag(i)
        ordered += [("N", "", "ALA", "A", r, "", (x - 1.0, y, z)), ("CA", "", "ALA", "A", r, "", (x, y, z)),
                    ("C", "", "ALA", "A", r, "", (x + 1.0, y, z))]
    shuffled = []
    rng = random.Random(7)
    for i in range(len(resids)):
        block = ordered[3 * i: 3 * i + 3]
        rng.shuffle(block)
        shuffled += block
    a = extract_ca(load(tmp_path, [ordered], "SYNTHETIC_a.pdb"), "SYNTHETIC_order")
    b = extract_ca(load(tmp_path, [shuffled], "SYNTHETIC_b.pdb"), "SYNTHETIC_order")
    np.testing.assert_array_equal(a.coordinates, b.coordinates)
    assert [r["resid"] for r in a.residues] == [r["resid"] for r in b.residues]


def test_correspondence_preserved_across_frames_and_reordering(tmp_path):
    resids = [(i, "") for i in range(1, 6)]
    models = [chain_atoms("A", resids, frame=f) for f in range(4)]
    u = load(tmp_path, models)
    order = [3, 1, 5, 2, 4]
    out = extract_ca(u, "SYNTHETIC_frames", correspondence=[{"chain": "A", "resid": r, "label": f"L{r}"}
                                                              for r in order])
    assert out.coordinates.shape == (4, 5, 3) and out.frame_index.tolist() == [0, 1, 2, 3]
    for f in range(4):
        for j, r in enumerate(order):
            np.testing.assert_array_equal(out.coordinates[f, j], np.array(zigzag(r - 1, frame=f), dtype=np.float32))
    assert len({r["atom_ix"] for r in out.residues}) == 5
    sub = extract_ca(u, "SYNTHETIC_frames", start=1, step=2)
    assert sub.frame_index.tolist() == [1, 3]
    np.testing.assert_array_equal(sub.coordinates[1], extract_ca(u, "x").coordinates[3])


def test_altloc_policies(tmp_path):
    atoms = chain_atoms("A", [(1, "")]) + [
        ("CA", "A", "ALA", "A", 2, "", zigzag(1), 0.4), ("CA", "B", "ALA", "A", 2, "", (3.8, -0.5, 0.3), 0.6)]
    u = load(tmp_path, [atoms])
    with pytest.raises(AdapterError, match="altloc policy is 'reject'"):
        extract_ca(u, "SYNTHETIC_alt")
    a = extract_ca(u, "SYNTHETIC_alt", altloc="blank_then_A")
    assert a.residues[1]["altloc"] == "A"
    occ = extract_ca(u, "SYNTHETIC_alt", altloc="highest_occupancy")
    assert occ.residues[1]["altloc"] == "B" and occ.report["altloc_decisions"][0]["used"] == "B"
    tie = chain_atoms("A", [(1, "")]) + [
        ("CA", "A", "ALA", "A", 2, "", zigzag(1), 0.5), ("CA", "B", "ALA", "A", 2, "", (3.8, -0.5, 0.3), 0.5)]
    with pytest.raises(AdapterError, match="tie"):
        extract_ca(load(tmp_path, [tie], "SYNTHETIC_tie.pdb"), "SYNTHETIC_alt", altloc="highest_occupancy")


def test_gro_loses_chain_and_icode_and_is_rejected_not_guessed(tmp_path):
    u = load(tmp_path, [chain_atoms("A", [(1, ""), (2, ""), (2, "A")]) + chain_atoms("B", [(1, "")], dy=20.0)])
    gro = tmp_path / "SYNTHETIC.gro"
    u.atoms.write(str(gro))
    with pytest.raises(AdapterError, match="without distinguishing altloc|duplicate residue key"):
        extract_ca(load_universe(gro), "SYNTHETIC_gro")


def test_unit_conversion_nm_to_angstrom(tmp_path):
    gro = tmp_path / "SYNTHETIC.gro"
    gro.write_text("SYNTHETIC unit fixture\n    3\n"
                   "    1ALA     CA    1   0.100   0.200   0.300\n"
                   "    2ALA     CA    2   0.480   0.200   0.300\n"
                   "    3ALA     CA    3   0.860   0.250   0.300\n"
                   "   5.00000   5.00000   5.00000\n")
    out = extract_ca(load_universe(gro), "SYNTHETIC_units")
    np.testing.assert_allclose(out.coordinates[0], [[1.0, 2.0, 3.0], [4.8, 2.0, 3.0], [8.6, 2.5, 3.0]], atol=1e-5)
    assert out.report["units"]["source_length_unit"] == "nm" and out.report["units"]["coordinates"] == "angstrom"
    xtc = tmp_path / "SYNTHETIC.xtc"
    u = load_universe(gro)
    with mda.Writer(str(xtc), n_atoms=3) as w:
        for _ in range(2):
            w.write(u.atoms)
    t = extract_ca(load_universe(gro, xtc), "SYNTHETIC_units")
    np.testing.assert_allclose(t.coordinates[1], out.coordinates[0], atol=0.01)   # XTC is lossy (writer precision)
    assert np.isfinite(t.time_ps).all() and t.report["units"]["time"] == "ps"
    assert np.isnan(extract_ca(load(tmp_path, [chain_atoms("A", [(1, ""), (2, "")])]), "p").time_ps).all()


def _memory_universe(positions, dimensions=None):
    from MDAnalysis.coordinates.memory import MemoryReader
    n = positions.shape[1]
    u = mda.Universe.empty(n, n_residues=n, atom_resindex=np.arange(n), trajectory=True)
    u.add_TopologyAttr("names", ["CA"] * n)
    u.add_TopologyAttr("resnames", ["ALA"] * n)
    u.add_TopologyAttr("resids", np.arange(1, n + 1))
    u.add_TopologyAttr("segids", ["SYN"])
    u.load_new(positions.astype(np.float32), format=MemoryReader, dimensions=dimensions)
    return u


def test_non_finite_coordinates_rejected():
    pos = np.array([[zigzag(k) for k in range(4)]] * 2, dtype=np.float64)
    pos[1, 2, 0] = np.nan
    with pytest.raises(AdapterError, match="non-finite.*frame 1"):
        extract_ca(_memory_universe(pos), "SYNTHETIC_nan", assume_angstrom_if_undeclared=True)


def test_pbc_policies():
    L = 10.0
    pos = np.array([[zigzag(k) for k in range(5)]], dtype=np.float64)
    wrapped = pos.copy()
    wrapped[0, 3:, 0] -= L                     # chain split by an orthorhombic boundary
    u = _memory_universe(wrapped, dimensions=[L, L, L, 90, 90, 90])
    raw = extract_ca(u, "SYNTHETIC_pbc", assume_angstrom_if_undeclared=True)
    assert any("ca_ca_distance" in b["reasons"] for b in raw.breaks)
    whole = extract_ca(u, "SYNTHETIC_pbc", pbc="make_whole_sequential", assume_angstrom_if_undeclared=True)
    np.testing.assert_allclose(whole.coordinates[0], pos[0], atol=1e-5)
    assert whole.breaks == [] and whole.report["pbc"]["frames_changed"] == 1
    tri = _memory_universe(wrapped, dimensions=[L, L, L, 90, 90, 60])
    with pytest.raises(AdapterError, match="not a valid orthorhombic"):
        extract_ca(tri, "SYNTHETIC_pbc", pbc="make_whole_sequential", assume_angstrom_if_undeclared=True)
    with pytest.raises(AdapterError, match="no box"):
        extract_ca(_memory_universe(wrapped), "SYNTHETIC_pbc", pbc="make_whole_sequential",
                   assume_angstrom_if_undeclared=True)
    assert raw.report["provenance"]["source_files"] == [] or all(
        f["sha256"] is None for f in raw.report["provenance"]["source_files"])


def test_pbc_unwrapping_never_crosses_chain_boundary():
    from MDAnalysis.coordinates.memory import MemoryReader
    L = 10.0
    a = [zigzag(k) for k in range(3)]                       # chain A ends near x = 7.6
    b = [(x + 13.6, y, z) for x, y, z in a]                 # chain B starts 6 Å past A's end: round(0.6) = 1 box
    pos = np.array([a + b], dtype=np.float32)
    u = mda.Universe.empty(6, n_residues=6, n_segments=2, atom_resindex=np.arange(6),
                           residue_segindex=[0, 0, 0, 1, 1, 1], trajectory=True)
    u.add_TopologyAttr("names", ["CA"] * 6)
    u.add_TopologyAttr("resnames", ["ALA"] * 6)
    u.add_TopologyAttr("resids", [1, 2, 3, 1, 2, 3])
    u.add_TopologyAttr("segids", ["SYA", "SYB"])
    u.load_new(pos, format=MemoryReader, dimensions=[L, L, L, 90, 90, 90])
    out = extract_ca(u, "SYNTHETIC_pbc_chains", pbc="make_whole_sequential", assume_angstrom_if_undeclared=True)
    assert {"after_index": 2, "reasons": ["chain_boundary"]} in out.breaks
    np.testing.assert_array_equal(out.coordinates[0], pos[0].astype(np.float64))


def test_precision_recorded_and_float32_values_kept_exactly(tmp_path):
    u = load(tmp_path, [chain_atoms("A", [(1, ""), (2, "")])])
    out = extract_ca(u, "SYNTHETIC_prec")
    assert out.coordinates.dtype == np.float64
    assert out.report["precision"]["reader_position_dtype"] == "float32"
    np.testing.assert_array_equal(out.coordinates[0], u.atoms.positions.astype(np.float64))


def test_end_to_end_local_file_to_export(tmp_path):
    src = write_pdb(tmp_path / "SYNTHETIC_e2e.pdb",
                    [chain_atoms("A", [(1, ""), (2, ""), (3, "")], frame=f) for f in range(3)])
    corr = tmp_path / "SYNTHETIC_corr.json"
    corr.write_text(json.dumps([{"chain": "A", "resid": r, "label": f"ref{r}"} for r in (1, 2, 3)]))
    from topofold_inputs.__main__ import main as cli
    assert cli([str(src), "--run-id", "SYNTHETIC_e2e", "--out", str(tmp_path / "out"),
                "--correspondence", str(corr)]) == 0
    rep = json.loads((tmp_path / "out" / "SYNTHETIC_e2e.json").read_text())
    npz_path = tmp_path / "out" / "SYNTHETIC_e2e.npz"
    assert rep["export"]["npz_sha256"] == hashlib.sha256(npz_path.read_bytes()).hexdigest()
    assert rep["provenance"]["source_files"][0]["sha256"] == hashlib.sha256(src.read_bytes()).hexdigest()
    data = np.load(npz_path)                    # no pickled objects
    assert data["coordinates_angstrom"].shape == (3, 3, 3)
    assert data["label"].tolist() == ["ref1", "ref2", "ref3"] and data["resid"].tolist() == [1, 2, 3]
    for f in range(3):
        for k in range(3):
            np.testing.assert_array_equal(data["coordinates_angstrom"][f, k],
                                          np.array(zigzag(k, frame=f), dtype=np.float32))
    assert rep["residues"][1]["label"] == "ref2" and rep["units"]["coordinates"] == "angstrom"


def test_invalid_configuration_rejected(tmp_path):
    u = load(tmp_path, [chain_atoms("A", [(1, ""), (2, "")])])
    for kw in ({"altloc": "first"}, {"pbc": "wrap"}, {"on_missing": "skip"}, {"max_ca_ca_angstrom": -1.0}):
        with pytest.raises(AdapterError):
            extract_ca(u, "SYNTHETIC_cfg", **kw)
    with pytest.raises(AdapterError, match="more than once"):
        extract_ca(u, "SYNTHETIC_cfg", correspondence=[{"chain": "A", "resid": 1}, {"chain": "A", "resid": 1}])
    with pytest.raises(FileNotFoundError):
        load_universe(tmp_path / "absent.pdb")


# --- regression tests: frame atom identity, insertion-code continuity, unit enforcement (2026-10-23) ----------------

def test_multimodel_pdb_with_reordered_atoms_in_later_model_is_rejected(tmp_path):
    m1 = chain_atoms("A", [(1, ""), (2, ""), (3, ""), (4, "")])
    m2 = chain_atoms("A", [(1, ""), (2, ""), (3, ""), (4, "")], frame=1)
    m2[1], m2[2] = m2[2], m2[1]                 # model 2 lists residue 3 before residue 2
    u = load(tmp_path, [m1, m2])
    with pytest.raises(AdapterError, match="atom identity differs"):
        extract_ca(u, "SYNTHETIC_reorder")


def test_frame_identity_check_is_reported(tmp_path):
    u = load(tmp_path, [chain_atoms("A", [(1, ""), (2, "")], frame=f) for f in range(2)])
    rep = extract_ca(u, "SYNTHETIC_ident").report["frame_identity"]
    assert rep["status"] == "verified" and rep["models_checked"] == 2
    gro = tmp_path / "SYNTHETIC.gro"
    u.atoms.write(str(gro))
    xtc = tmp_path / "SYNTHETIC.xtc"
    with mda.Writer(str(xtc), n_atoms=u.atoms.n_atoms) as w:
        for _ in u.trajectory:
            w.write(u.atoms)
    rep = extract_ca(load_universe(gro, xtc), "SYNTHETIC_ident", chain_attribute="segid").report["frame_identity"]
    assert rep["status"] == "not_verifiable"


@pytest.mark.parametrize("seq,expect_break", [
    ([(52, ""), (52, "A"), (52, "B"), (53, "")], False),   # ordinary insertion run
    ([(52, ""), (52, "B")], True),                          # skipped insertion code
    ([(52, "B"), (52, "A")], True),                         # reversed insertion codes
    ([(52, "A"), (52, "")], True),                          # insertion before its base residue
    ([(52, ""), (53, "A")], True),                          # next number starting with an insertion code
    ([(52, "A"), (53, "")], False),                         # insertion run then next number
])
def test_insertion_code_continuity_is_conservative(tmp_path, seq, expect_break):
    out = extract_ca(load(tmp_path, [chain_atoms("A", seq)]), "SYNTHETIC_icode_cont")
    has = any("numbering_discontinuity" in b["reasons"] for b in out.breaks)
    assert has == expect_break


def test_distance_never_establishes_continuity(tmp_path):
    # Numbering jump with physically adjacent Cα atoms: still a break.
    out = extract_ca(load(tmp_path, [chain_atoms("A", [(1, ""), (2, ""), (7, ""), (8, "")])]), "SYNTHETIC_dist")
    assert {"after_index": 1, "reasons": ["numbering_discontinuity"]} in out.breaks


def test_unconverted_units_are_rejected(tmp_path):
    gro = tmp_path / "SYNTHETIC.gro"
    gro.write_text("SYNTHETIC unit fixture\n    2\n"
                   "    1ALA     CA    1   0.100   0.200   0.300\n"
                   "    2ALA     CA    2   0.480   0.200   0.300\n"
                   "   5.00000   5.00000   5.00000\n")
    u = mda.Universe(str(gro), convert_units=False)
    with pytest.raises(AdapterError, match="unit"):
        extract_ca(u, "SYNTHETIC_units_off")


def test_undeclared_units_require_explicit_opt_in():
    pos = np.array([[zigzag(k) for k in range(3)]], dtype=np.float64)
    with pytest.raises(AdapterError, match="declares no length unit"):
        extract_ca(_memory_universe(pos), "SYNTHETIC_mem")
    out = extract_ca(_memory_universe(pos), "SYNTHETIC_mem", assume_angstrom_if_undeclared=True)
    assert "assumed by caller" in out.report["units"]["basis"]
