"""Modified-residue inventory helpers (benchmarks/audit/modified_residue_inventory_2026-10-26): SYNTHETIC coordinates with known truth."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
INV = ROOT / "benchmarks/audit/modified_residue_inventory_2026-10-26"
sys.path.insert(0, str(INV))
import inventory as inv  # noqa: E402


def _chain(chain="A", n=8, mod=None, free=None):
    """CA/N/C/O helix-like chain. mod: {resseq: (resname, record)} replaces the residue name; free: resseq of a residue moved 30 Å away (not bonded)."""
    lines = []
    t = np.arange(n) * np.radians(100)
    X = np.c_[2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)]
    k = 0
    for i in range(n):
        resseq = i + 1
        name, rec = (mod or {}).get(resseq, ("ALA", "ATOM  "))
        a = 0.326                                            # C(i)-N(i+1) = 3.83 * (1 - 2a) = 1.33 Å for equal Cα spacing
        fwd = X[min(i + 1, n - 1)] - X[i] if i < n - 1 else X[i] - X[i - 1]
        back = X[i] - X[i - 1] if i else X[1] - X[0]
        u = fwd
        shift = np.array([30.0, 0, 0]) if free == resseq else 0.0
        for an, p in (("N", X[i] - a * back), ("CA", X[i]), ("C", X[i] + a * fwd), ("O", X[i] + a * fwd + np.array([0, 0, 1.0]))):
            p = p + shift
            k += 1
            lines.append(f"{rec:<6s}{k:5d} {an:^4s} {name:>3s} {chain}{resseq:4d}    {p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00 20.00           {an[0]}")
    return lines


def test_scan_finds_modified_residues_from_coordinates_even_without_modres():
    text = "\n".join(_chain(mod={3: ("TPO", "HETATM"), 6: ("XYZ", "HETATM")}) + ["END"]) + "\n"
    chains, n_models = inv.scan(text)
    assert n_models == 0
    res = chains["A"]
    nonstd = [(r["resseq"], r["resname"], r["record"]) for r in res if r["resname"] not in inv.STD20]
    assert nonstd == [(3, "TPO", "HETATM"), (6, "XYZ", "HETATM")]                  # found without any MODRES record
    assert all(inv.pick(res[2]["atoms"].get(a)) is not None for a in ("N", "CA", "C", "O"))
    assert inv.bonded(res[1], res[2]) < inv.PEPTIDE_BOND_MAX and inv.bonded(res[2], res[3]) < inv.PEPTIDE_BOND_MAX


def test_a_backbone_residue_that_is_not_bonded_is_not_called_a_polymer_member():
    text = "\n".join(_chain(mod={4: ("ABC", "HETATM")}, free=4) + ["END"]) + "\n"
    res = inv.scan(text)[0]["A"]
    r = res[3]
    assert r["resname"] == "ABC" and inv.bonded(res[2], r) > inv.PEPTIDE_BOND_MAX and inv.bonded(r, res[4]) > inv.PEPTIDE_BOND_MAX


def test_modres_is_read_as_annotation_including_a_residue_without_coordinates():
    head = "MODRES 1ABC SEP A    9  SER  MODIFIED RESIDUE\nMODRES 1ABC TPO A    3  THR  MODIFIED RESIDUE\nREMARK 465   M RES C SSSEQI\nREMARK 465     GLY A    12\n"
    modres, r465 = inv.parse_header(head + "\n".join(_chain(mod={3: ("TPO", "HETATM")})) + "\n")
    assert modres[("A", 9, "")] == {"resname": "SEP", "parent": "SER", "comment": "MODIFIED RESIDUE"}
    assert ("A", 12, "") in r465
    chains, _ = inv.scan(head + "\n".join(_chain(mod={3: ("TPO", "HETATM")})) + "\n")
    assert ("A", 9, "") not in {(r["chain"], r["resseq"], r["icode"]) for r in chains["A"]}      # annotated, but no coordinate record (category 1)


def test_only_the_first_model_is_read():
    one = "\n".join(_chain())
    text = "MODEL        1\n" + one + "\nENDMDL\nMODEL        2\n" + one.replace(" ALA ", " TPO ") + "\nENDMDL\nEND\n"
    chains, n_models = inv.scan(text)
    assert n_models == 2 and all(r["resname"] == "ALA" for r in chains["A"])


def _adapter_python():
    cand = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
    if not Path(cand).is_file():
        return None
    return cand if subprocess.run([cand, "-c", "import MDAnalysis"], capture_output=True).returncode == 0 else None


def test_the_real_adapter_drops_phosphorylated_names_but_keeps_standard_and_selenomethionine(tmp_path, monkeypatch):
    py = _adapter_python()
    if py is None:
        pytest.skip("no interpreter with MDAnalysis available")
    monkeypatch.setattr(inv, "ADAPTER", py)
    monkeypatch.setattr(inv, "WORK", tmp_path / "work")
    f = tmp_path / "SYNTHETIC.pdb"
    f.write_text("\n".join(_chain(mod={2: ("SEP", "HETATM"), 3: ("TPO", "HETATM"), 4: ("PTR", "HETATM"), 5: ("MSE", "HETATM")})) + "\nEND\n")
    got, _, _ = inv.adapter_residues(f, "diag")
    present = {k[1]: v["resname"] for k, v in got.items()}
    assert 2 not in present and 3 not in present and 4 not in present                  # SEP, TPO, PTR are absent from the adapter output
    assert present[5] == "MSE" and present[1] == "ALA" and present[6] == "ALA"          # selenomethionine and standard residues are read
    assert set(inv.STD20) >= {"ALA", "SER"} and "SEP" not in inv.STD20 and "SEP" not in kp_prot_res()


def kp_prot_res():
    import kinpair
    return kinpair.PROT_RES
