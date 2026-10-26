"""KIN-GEN-2 gates G3 and G8 against MDAnalysis (run with the adapter interpreter).

    python gate_mda.py JOB.json OUT.json
JOB: pairs [[windowA, windowB], ...] (9x3 Cα arrays) and d1d2 [{file, chain, anchors, D1, D2}]. OUT: MDAnalysis QCP RMSD per pair and the largest D1/D2 deviation.
"""

import json
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.rms import rmsd

job = json.loads(Path(sys.argv[1]).read_text())
out = {"rmsd": [float(rmsd(np.array(a, float), np.array(b, float), center=True, superposition=True)) for a, b in job["pairs"]]}
worst = 0.0
for e in job["d1d2"]:
    u = mda.Universe(e["file"])
    def atom(role, name):
        ch, num = e["anchors"][role].split(":")
        resid = int("".join(c for c in num if c.isdigit() or c == "-"))
        icode = "".join(c for c in num if c.isalpha())
        sel = u.select_atoms(f"segid {ch} or chainID {ch}").select_atoms(f"resid {resid} and name {name}")
        sel = [a for a in sel if str(getattr(a, "altLoc", "") or "") in ("", "A")]
        return sel[0].position
    d1 = float(np.linalg.norm(atom("alphaC_Glu_plus4", "CA") - atom("DFG_Phe", "CZ")))
    d2 = float(np.linalg.norm(atom("beta3_Lys", "CA") - atom("DFG_Phe", "CZ")))
    worst = max(worst, abs(d1 - e["D1"]), abs(d2 - e["D2"]))
out["max_d1d2_deviation"] = worst
Path(sys.argv[2]).write_text(json.dumps(out))
