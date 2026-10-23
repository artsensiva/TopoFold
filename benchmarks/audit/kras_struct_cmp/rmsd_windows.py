"""Per-window proper-rotation Calpha RMSD with MDAnalysis QCP (run in an interpreter with MDAnalysis).

    python rmsd_windows.py A.npz B.npz OUT.json
Reads adapter exports (coordinates_angstrom, one frame each); windows of 3 and 4 consecutive output residues, plus a
global value over all residues present in both.
"""
import json
import sys

import numpy as np
from MDAnalysis.analysis.rms import rmsd

a = np.load(sys.argv[1])["coordinates_angstrom"][0]
b = np.load(sys.argv[2])["coordinates_angstrom"][0]
n = len(a)
ok = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
out = {"rmsd3": [None] * n, "rmsd4": [None] * n,
       "global": float(rmsd(a[ok], b[ok], center=True, superposition=True)), "n_global": int(ok.sum())}
for i in range(n):
    if 1 <= i <= n - 2 and ok[i - 1: i + 2].all():
        out["rmsd3"][i] = float(rmsd(a[i - 1: i + 2], b[i - 1: i + 2], center=True, superposition=True))
    if 1 <= i <= n - 3 and ok[i - 1: i + 3].all():
        out["rmsd4"][i] = float(rmsd(a[i - 1: i + 3], b[i - 1: i + 3], center=True, superposition=True))
json.dump(out, open(sys.argv[3], "w"))
