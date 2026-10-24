"""Build tests/fixtures/kin_rob_1_panel_observed_chains.json from the KIN-ROB-1 inputs (public data, not SYNTHETIC).

    python3 benchmarks/audit/kin_map_1/build_panel_fixture.py

Reads the git-ignored KIN-ROB-1 downloads (benchmarks/audit/kin_rob_1/data_pdb, data_uniprot), refuses to run if any
SHA-256 differs from KIN-ROB-1_results.json, and writes, per panel chain, the observed residues in file order with
peptide-bond flags. Parsing is an independent fixed-column reader (no MDAnalysis, no TopoFold adapter): first model,
ATOM records plus HETATM records whose residue name is a standard or listed variant amino acid; a residue is bonded to
the previous one when the shortest C(i-1)-N(i) distance over alternate locations is <= 2.0 Å. The UniProt sequences
(CC BY 4.0) are included with their headers and hashes so the committed tests run without the downloads.
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ROB = ROOT / "benchmarks/audit/kin_rob_1"
OUT = ROOT / "tests/fixtures/kin_rob_1_panel_observed_chains.json"
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs.kinase_roles import ONE_LETTER, PEPTIDE_BOND_MAX  # noqa: E402  (residue-name table only)


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def observed(path, chain):
    res, order = {}, []
    for line in Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")) or line[21] != chain:
            continue
        name3 = line[17:20].strip()
        if line.startswith("HETATM") and name3 not in ONE_LETTER:
            continue
        key = (int(line[22:26]), line[26].strip())
        if key not in res:
            res[key] = {"resname": name3, "C": [], "N": []}
            order.append(key)
        atom = line[12:16].strip()
        if atom in ("C", "N"):
            res[key][atom].append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    out = []
    for k, key in enumerate(order):
        bonded = None
        if k:
            c, n = np.array(res[order[k - 1]]["C"]), np.array(res[key]["N"])
            if len(c) and len(n):
                bonded = bool(np.linalg.norm(c[:, None] - n[None], axis=-1).min() <= PEPTIDE_BOND_MAX)
        out.append([key[0], key[1], res[key]["resname"], bonded])
    return out


def main():
    spec = json.loads((ROB / "run_spec.json").read_text())
    recorded = json.loads((ROB / "results/KIN-ROB-1_results.json").read_text())["inputs"]
    bad = [k for k, v in recorded.items() if sha256(ROOT / v["path"]) != v["sha256"]]
    if bad:
        raise SystemExit(f"STOP: inputs differ from KIN-ROB-1 records: {bad}")
    fx = {"description": "Observed residues of the KIN-ROB-1 panel chains and their UniProt sequences. Public data "
                         "(PDB archive CC0; UniProt CC BY 4.0), not SYNTHETIC. Built by "
                         "benchmarks/audit/kin_map_1/build_panel_fixture.py; residue rows are "
                         "[resid, icode, resname, bonded_to_previous].",
          "peptide_bond_max_angstrom": PEPTIDE_BOND_MAX, "sources": {}, "uniprot": {}, "chains": []}
    for entry in spec["panel"]:
        pdb, acc = entry["pdb"], entry["uniprot"]
        for key in (pdb, acc):
            fx["sources"][key] = {"url": recorded[key]["url"], "sha256": recorded[key]["sha256"]}
        text = (ROB / f"data_uniprot/{acc}.fasta").read_text()
        fx["uniprot"][acc] = {"header": text.splitlines()[0],
                              "sequence": "".join(ln.strip() for ln in text.splitlines()[1:])}
        for chain in entry["chains"]:
            fx["chains"].append({"pdb": pdb, "chain": chain, "uniprot": acc,
                                 "residues": observed(ROB / f"data_pdb/{pdb}.pdb", chain)})
    OUT.write_text(json.dumps(fx, separators=(",", ":")) + "\n")
    print(OUT.relative_to(ROOT), OUT.stat().st_size, "bytes", len(fx["chains"]), "chains")


if __name__ == "__main__":
    main()
