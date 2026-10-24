"""Tests for topofold_inputs.atom_adapter and kinase_geometry (capability "B").

SYNTHETIC fixtures for policies; public PDB structures (KIN-DESC-1 and KRAS inputs, git-ignored, skipped if absent) for
reference checks against independently computed recorded values and against the Cα adapter.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs.kinase_geometry import DESCRIPTOR_ATOMS  # noqa: E402

# --- extraction (MDAnalysis needed) ----------------------------------------------------------------------------------

mda = pytest.importorskip("MDAnalysis")
from topofold_inputs.atom_adapter import KINASE_DESCRIPTOR_ROLES, extract_atoms, role_request  # noqa: E402
from topofold_inputs.ca_adapter import AdapterError, ResidueSpec, extract_ca, load_universe  # noqa: E402
from topofold_inputs.kinase_geometry import kinase_descriptors  # noqa: E402


def _atom_line(serial, name, resname, chain, resid, xyz, icode=" ", altloc=" ", occ=1.0):
    nm = f" {name:<3}" if len(name) < 4 else name
    return (f"ATOM  {serial:5d} {nm}{altloc}{resname:>3} {chain}{resid:4d}{icode}   "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{occ:6.2f} 20.00           {name[0]}")


def synthetic_pdb(path, *, drop=(), altloc_occ=(0.4, 0.6)):
    """SYNTHETIC chain A: GLY 10, PHE 11 (CB/CG/CZ in altlocs A/B), GLY 11A (insertion code), LYS 12."""
    lines, n = [], 0
    rng = np.random.default_rng(0)
    for resid, icode, resname, names in ((10, " ", "GLY", ("N", "CA", "C", "O")),
                                         (11, " ", "PHE", ("N", "CA", "C", "O", "CB", "CG", "CZ")),
                                         (11, "A", "GLY", ("N", "CA", "C", "O")),
                                         (12, " ", "LYS", ("N", "CA", "C", "O", "CB", "NZ"))):
        base = np.array([3.8 * len(lines) / 6.0, 0.0, 0.0])
        for name in names:
            if (resid, icode.strip(), name) in drop:
                continue
            xyz = base + rng.normal(scale=1.0, size=3)
            if resname == "PHE" and name in ("CB", "CG", "CZ"):
                for alt, occ, shift in (("A", altloc_occ[0], 0.0), ("B", altloc_occ[1], 1.5)):
                    n += 1
                    lines.append(_atom_line(n, name, resname, "A", resid, xyz + shift, icode, alt, occ))
            else:
                n += 1
                lines.append(_atom_line(n, name, resname, "A", resid, xyz, icode))
    Path(path).write_text("\n".join(lines) + "\nEND\n")
    return path


def _req(*items):
    return [(ResidueSpec(chain="A", resid=r, icode=ic), names) for r, ic, names in items]


def test_missing_atom_mark_and_reject(tmp_path):
    u = load_universe(synthetic_pdb(tmp_path / "SYNTHETIC.pdb", drop={(12, "", "NZ")}))
    req = _req((12, "", ["CA", "NZ"]), (99, "", ["CA"]))
    t = extract_atoms(u, "SYNTHETIC", req)
    st = {(a["resid"], a["name"]): a["status"] for a in t.atoms}
    assert st == {(12, "CA"): "present", (12, "NZ"): "atom_not_found", (99, "CA"): "residue_not_found"}
    assert np.isnan(t.xyz("A", 12, "NZ")).all() and np.isfinite(t.xyz("A", 12, "CA")).all()
    assert t.report["summary"]["missing_by_status"] == {"residue_not_found": 1, "atom_not_found": 1}
    with pytest.raises(AdapterError, match="atom NZ not found"):
        extract_atoms(u, "SYNTHETIC", _req((12, "", ["CA", "NZ"])), on_missing="reject")
    with pytest.raises(AdapterError, match="not found in selection"):
        extract_atoms(u, "SYNTHETIC", _req((99, "", ["CA"])), on_missing="reject")


def test_altloc_policies_are_residue_consistent(tmp_path):
    u = load_universe(synthetic_pdb(tmp_path / "SYNTHETIC.pdb"))
    req = _req((11, "", ["CA", "CB", "CG", "CZ"]))
    with pytest.raises(AdapterError, match="altloc policy is 'reject'"):
        extract_atoms(u, "S", req, altloc="reject")
    a = extract_atoms(u, "S", req, altloc="blank_then_A")
    assert [x["altloc"] for x in a.atoms] == ["", "A", "A", "A"]
    b = extract_atoms(u, "S", req, altloc="highest_occupancy")
    assert [x["altloc"] for x in b.atoms] == ["", "B", "B", "B"]            # one conformer for the whole residue
    assert b.report["altloc_decisions"][0]["used"] == "B"
    assert np.allclose(b.xyz("A", 11, "CZ") - a.xyz("A", 11, "CZ"), 1.5, atol=1e-3)
    u2 = load_universe(synthetic_pdb(tmp_path / "SYNTHETIC_tie.pdb", altloc_occ=(0.5, 0.5)))
    with pytest.raises(AdapterError, match="tie on mean occupancy"):
        extract_atoms(u2, "S", req, altloc="highest_occupancy")


def test_insertion_code_and_request_refusals(tmp_path):
    u = load_universe(synthetic_pdb(tmp_path / "SYNTHETIC.pdb"))
    t = extract_atoms(u, "S", _req((11, "A", ["N", "CA"]), (11, "", ["CA"])), altloc="blank_then_A")
    assert [(a["resid"], a["icode"], a["resname"]) for a in t.atoms] == [(11, "A", "GLY"), (11, "A", "GLY"),
                                                                        (11, "", "PHE")]
    with pytest.raises(AdapterError, match="more than once"):
        extract_atoms(u, "S", _req((10, "", ["CA"]), (10, "", ["N"])))
    with pytest.raises(AdapterError, match="unique"):
        extract_atoms(u, "S", _req((10, "", ["CA", "CA"])))
    with pytest.raises(AdapterError, match="roles not mapped"):
        role_request({"beta3_Lys": ResidueSpec("A", 12)})


def test_descriptor_missing_atom_reports_reason_not_zero(tmp_path):
    u = load_universe(synthetic_pdb(tmp_path / "SYNTHETIC.pdb", drop={(11, "", "CZ")}))
    roles = {r: ResidueSpec("A", 11) for r in KINASE_DESCRIPTOR_ROLES}  # SYNTHETIC: every role on one residue
    roles.update(beta3_Lys=ResidueSpec("A", 12), alphaC_Glu=ResidueSpec("A", 12))
    t = extract_atoms(u, "S", role_request(roles), altloc="blank_then_A")
    out = kinase_descriptors(t, roles)
    assert np.isnan(out["values"]["D1"]).all() and np.isnan(out["values"]["D2"]).all()
    assert out["missing"]["D1"] == [("DFG_Phe", "CZ", "atom_not_found")]
    assert np.isfinite(out["values"]["alphaC_CbCb"]).all()


# --- public structures --------------------------------------------------------------------------------------------

KIN = ROOT / "benchmarks/audit/kinase_descriptors"
KIN_MAP = {"D1": "D1", "D2": "D2", "alphaC_CbCb": "alphaC_CbCb", "Aloop_NO": "Aloop_Hbond_NO",
           "chi1_Phe": "Phe382_chi1", "phi_XDFG": "phi380", "psi_XDFG": "psi380", "phi_Asp": "phi381",
           "psi_Asp": "psi381", "phi_Phe": "phi382", "psi_Phe": "psi382"}
ABL_ROLES_1A = {"XDFG_minus1": 379, "beta3_Lys": 271, "alphaC_Glu": 286, "alphaC_Glu_plus4": 290, "XDFG": 380,
                "DFG_Asp": 381, "DFG_Phe": 382, "DFG_Gly": 383, "DFGxxX": 386, "XHRD": 360}


@pytest.mark.skipif(not (KIN / "data_pdb/1OPJ.pdb").is_file(), reason="KIN-DESC-1 inputs not present (git-ignored)")
def test_abl_descriptors_match_independent_kin_desc_1_values():
    """MDAnalysis-based extraction vs KIN-DESC-1's frozen values from an independent fixed-column parser."""
    spec = json.loads((KIN / "run_spec.json").read_text())
    rec = json.loads((KIN / "results/KIN-DESC-1_results.json").read_text())["established_values"]
    n_cmp = 0
    for pid, s in spec["structures"].items():
        u = load_universe(ROOT / s["file"])
        off = -s["numbering_offset_to_cAbl_1a"]      # file numbering = 1a numbering + off
        for ch in s["chains"]:
            roles = {r: ResidueSpec(ch, n + off) for r, n in ABL_ROLES_1A.items()}
            t = extract_atoms(u, f"{pid}{ch}", role_request(roles), altloc="blank_then_A", on_missing="reject")
            out = kinase_descriptors(t, roles)
            assert not out["missing"]
            for ours, theirs in KIN_MAP.items():
                v, ref = float(out["values"][ours][0]), rec[f"{pid}:{ch}"][theirs]
                angular = len(DESCRIPTOR_ATOMS[ours]) == 4          # dihedrals: compare on the circle (rad)
                diff = abs(math.atan2(math.sin(v - ref), math.cos(v - ref))) if angular else abs(v - ref)
                assert diff < 2e-4, (pid, ch, ours, v, ref)              # Å or rad; PDB coordinates have 0.001 Å steps
                n_cmp += 1
    assert n_cmp == 4 * len(KIN_MAP)


KRAS = ROOT / "benchmarks/audit/kras_struct_cmp/data_pdb/9O0N.pdb"


@pytest.mark.skipif(not KRAS.is_file(), reason="KRAS inputs not present (git-ignored)")
@pytest.mark.parametrize("policy", ["blank_then_A", "highest_occupancy"])
def test_calpha_identical_to_ca_adapter_on_public_structure(policy):
    """Cα coordinates from the atom adapter equal the Cα adapter's, bit for bit, on 9O0N (which has altlocs)."""
    u = load_universe(KRAS)
    corr = [ResidueSpec("A", r) for r in range(1, 169)]
    ca = extract_ca(u, "ca", correspondence=corr, on_missing="reject", altloc=policy)
    at = extract_atoms(u, "at", [(s, ["CA"]) for s in corr], on_missing="reject", altloc=policy)
    assert np.array_equal(ca.coordinates, at.coordinates)
    assert [d["resid"] for d in ca.report["altloc_decisions"]] == [d["resid"] for d in at.report["altloc_decisions"]]
