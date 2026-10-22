"""Engineering example for the Cα input adapter on a SYNTHETIC structure (not biological validation).

Writes a SYNTHETIC three-model PDB and a residue-correspondence JSON into --out, then runs the adapter's command-line
path on them (local file -> MDAnalysis -> adapter -> exported .npz + .json). The structure contains, on purpose:
duplicate residue numbers across chains A and B, an insertion code (12A), alternate Cα locations (A13), a numbering
gap (13 -> 16), a residue without a Cα atom (A18), and a requested residue absent from the file (A14).

Run (from the repository root, in an environment that has MDAnalysis):
    PYTHONPATH=python benchmarks/audit/wpd_local_retrieval/build/venv/bin/python \
        python/examples/ca_adapter_example.py --out <output directory>
"""

import argparse
import json
from pathlib import Path

from topofold_inputs.__main__ import main as adapter_cli

# (chain, resid, icode, resname, atoms) with atoms as (name, altloc, occupancy, offset index)
RESIDUES = [
    ("A", 10, "", "ALA", [("CA", "", 1.00, 0)]),
    ("A", 11, "", "GLY", [("N", "", 1.00, 1), ("CA", "", 1.00, 1)]),
    ("A", 12, "", "SER", [("CA", "", 1.00, 2)]),
    ("A", 12, "A", "SER", [("CA", "", 1.00, 3)]),
    ("A", 13, "", "LEU", [("CA", "A", 0.60, 4), ("CA", "B", 0.40, 4)]),
    ("A", 16, "", "VAL", [("CA", "", 1.00, 7)]),
    ("A", 17, "", "ILE", [("CA", "", 1.00, 8)]),
    ("A", 18, "", "LYS", [("N", "", 1.00, 9)]),
    ("B", 10, "", "ALA", [("CA", "", 1.00, 0)]),
    ("B", 11, "", "GLY", [("CA", "", 1.00, 1)]),
    ("B", 12, "", "SER", [("CA", "", 1.00, 2)]),
]


def _xyz(chain: str, k: int, alt: str, model: int):
    """Zig-zag trace (consecutive Cα–Cα 3.93 Å); chain B offset by 20 Å; small per-model displacement."""
    x, y, z = 3.8 * k, 0.5 * (-1) ** k, 0.0
    if chain == "B":
        y += 20.0
    if alt == "B":
        z += 0.3
    return x, y + 0.1 * model * ((k % 3) - 1), z + 0.05 * model


def write_synthetic_pdb(path: Path, n_models: int = 3) -> None:
    lines = ["REMARK   1 SYNTHETIC engineering fixture for topofold_inputs; not experimental data"]
    for m in range(n_models):
        lines.append(f"MODEL     {m + 1:4d}")
        serial = 1
        for chain, resid, icode, resname, atoms in RESIDUES:
            for name, alt, occ, k in atoms:
                x, y, z = _xyz(chain, k, alt, m)
                lines.append(f"ATOM  {serial:5d} {(' ' + name):<4s}{alt or ' ':1s}{resname:>3s} {chain:1s}"
                             f"{resid:4d}{icode or ' ':1s}   {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{0.0:6.2f}"
                             f"          {name[0]:>2s}")
                serial += 1
        lines += ["TER", "ENDMDL"]
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="output directory (created if needed)")
    out = Path(ap.parse_args().out)
    out.mkdir(parents=True, exist_ok=True)
    pdb = out / "SYNTHETIC_ca_adapter_fixture.pdb"
    write_synthetic_pdb(pdb)
    corr = [{"chain": "A", "resid": r, "icode": ic, "label": f"SYNTH_A:{r}{ic}"}
            for r, ic in [(10, ""), (11, ""), (12, ""), (12, "A"), (13, ""), (14, ""), (16, ""), (17, ""), (18, "")]]
    corr += [{"chain": "B", "resid": r, "label": f"SYNTH_B:{r}"} for r in (10, 11, 12)]
    corr_path = out / "SYNTHETIC_correspondence.json"
    corr_path.write_text(json.dumps(corr, indent=2) + "\n")
    return adapter_cli([str(pdb), "--run-id", "SYNTHETIC_example_run", "--out", str(out),
                        "--correspondence", str(corr_path), "--altloc", "highest_occupancy"])


if __name__ == "__main__":
    raise SystemExit(main())
