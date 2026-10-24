"""Structure-side tests for topofold_inputs.kinase_roles (MDAnalysis needed).

SYNTHETIC PDB files written by the tests and the committed SYNTHETIC altloc fixture run everywhere MDAnalysis is
installed; the comparison with the independent parser on the KIN-ROB-1 panel PDB files (git-ignored) is skipped when
those files are absent.
"""

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

# With TOPOFOLD_REQUIRE_MDANALYSIS_TESTS=1 (the CI adapter job) a missing MDAnalysis fails instead of skipping.
if os.environ.get("TOPOFOLD_REQUIRE_MDANALYSIS_TESTS") == "1":
    import MDAnalysis as mda  # noqa: F401
else:
    mda = pytest.importorskip("MDAnalysis")
from topofold_inputs.atom_adapter import extract_atoms  # noqa: E402
from topofold_inputs.ca_adapter import load_universe  # noqa: E402
from topofold_inputs.kinase_geometry import kinase_descriptors  # noqa: E402
from topofold_inputs.kinase_roles import (STANDARD, Numbering, ObservedResidue, map_kinase_roles_in_universe,  # noqa: E402
                                          observed_chain)
from test_kinase_roles import POS, synthetic_reference  # noqa: E402

THREE = {v: k for k, v in STANDARD.items()}


def write_synthetic_pdb(path, ref, *, drop=(), icodes=None, altloc_resid=None):
    """SYNTHETIC straight-chain backbone: N, CA, C, O, CB per residue (no CB on Gly), C(i)–N(i+1) = 1.33 Å; dropped
    residues leave their space empty, so the flanking C–N distance exceeds 2 Å. Author resid = reference index + 500."""
    lines, n = [], 0
    for i, aa in enumerate(ref.sequence, start=1):
        if i in drop:
            continue
        resid, icode = (icodes or {}).get(i, (i + 500, " "))
        x0 = 3.82 * i
        atoms = [("N", (x0, 0, 0)), ("CA", (x0 + 1.46, 0, 0)), ("C", (x0 + 2.49, 0, 0)), ("O", (x0 + 2.49, 1.23, 0))]
        if aa != "G":
            atoms.append(("CB", (x0 + 1.46, -1.53, 0)))
        for name, xyz in atoms:
            alts = [("A", 0.6, 0.0), ("B", 0.4, 0.5)] if (i == altloc_resid and name == "CB") else [(" ", 1.0, 0.0)]
            for alt, occ, dz in alts:
                n += 1
                nm = f" {name:<3}"
                lines.append(f"ATOM  {n:5d} {nm}{alt}{THREE[aa]:>3} A{resid:4d}{icode}   "
                             f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2] + dz:8.3f}{occ:6.2f} 20.00           {name[0]}")
    Path(path).write_text("REMARK   1 SYNTHETIC kinase-like chain for kinase_roles tests\n" + "\n".join(lines) + "\nEND\n")
    return path


def test_observed_chain_breaks_icodes_and_altlocs_on_synthetic_pdb(tmp_path):
    ref = synthetic_reference()
    pdb = write_synthetic_pdb(tmp_path / "SYNTHETIC_kin.pdb", ref, drop=range(102, 108),
                              icodes={100: (599, "A"), 101: (599, "B")}, altloc_resid=60)
    obs = observed_chain(load_universe(pdb), "A")
    assert len(obs) == len(ref.sequence) - 6
    assert [o for o in obs if o.bonded_to_previous is False] == [ObservedResidue(608, "", THREE[ref.sequence[107]],
                                                                                 False)]
    assert sum(o.bonded_to_previous is None for o in obs) == 1   # the first residue only
    assert (599, "A") in [(o.resid, o.icode) for o in obs] and [o.resid for o in obs].count(560) == 1


def test_mapping_feeds_extraction_and_descriptors(tmp_path):
    ref = synthetic_reference()
    pdb = write_synthetic_pdb(tmp_path / "SYNTHETIC_kin.pdb", ref, drop=range(102, 108),
                              icodes={100: (599, "A"), 101: (599, "B")})
    u = load_universe(pdb)
    m = map_kinase_roles_in_universe(u, "A", ref, Numbering.reference_index(), POS)
    assert m.roles["DFG_Gly"] is None and m.roles["DFGxxX"] is None
    assert (m.roles["DFG_Asp"].resid, m.roles["DFG_Asp"].icode) == (599, "A")
    traj = extract_atoms(u, "SYNTHETIC", m.atom_request(), on_missing="mark", altloc="blank_then_A")
    d = kinase_descriptors(traj, m.residue_specs())
    assert d["missing"]["psi_Phe"] == [("DFG_Gly", "N", "role_not_mapped")]
    assert np.isnan(d["values"]["Aloop_NO"][0])
    # alphaC Cβ–Cβ against the written coordinates: Cβ of residue i is at x = 3.82 i + 1.46
    assert d["values"]["alphaC_CbCb"][0] == pytest.approx(3.82 * (POS["alphaC_Glu"] - POS["beta3_Lys"]), abs=2e-3)


def test_observed_chain_on_committed_altloc_fixture():
    u = load_universe(ROOT / "tests/fixtures/SYNTHETIC_kinase_like_altloc.pdb")
    obs = observed_chain(u, "A")
    keys = [(o.resid, o.icode) for o in obs]
    assert len(keys) == len(set(keys)) and (6, "A") in keys and keys.count((4, "")) == 1


PANEL = ROOT / "benchmarks/audit/kin_rob_1/data_pdb"
FX = ROOT / "tests/fixtures/kin_rob_1_panel_observed_chains.json"


@pytest.mark.skipif(not (PANEL / "1ATP.pdb").is_file(), reason="KIN-ROB-1 inputs not present (git-ignored)")
def test_observed_chain_matches_independent_parser_on_kin_rob_1_panel():
    for c in json.loads(FX.read_text())["chains"]:
        obs = observed_chain(load_universe(PANEL / f"{c['pdb']}.pdb"), c["chain"])
        assert [[o.resid, o.icode, o.resname, o.bonded_to_previous] for o in obs] == c["residues"], c["pdb"]
