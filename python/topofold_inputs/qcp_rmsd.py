"""Per-window and global proper-rotation Cα RMSD with MDAnalysis QCP (run in an interpreter with MDAnalysis).

    python -m topofold_inputs.qcp_rmsd A.npz B.npz OUT.json

Reads two single-frame adapter exports with the same residue order. Windows of 3 and 4 consecutive output rows, indexed
by the κ / τ-hinge convention (window i = rows i−1 .. i+1 or i−1 .. i+2); a window is computed only where all its rows
have finite coordinates in both. The global value uses every row finite in both. Same calculation as
``benchmarks/audit/kras_struct_cmp/rmsd_windows.py`` (MDAnalysis.analysis.rms.rmsd, center and superposition).
"""

import json
import sys

import numpy as np


def main(argv=None) -> int:
    import MDAnalysis
    from MDAnalysis.analysis.rms import rmsd

    a_path, b_path, out = (argv or sys.argv[1:])
    a = np.load(a_path)["coordinates_angstrom"]
    b = np.load(b_path)["coordinates_angstrom"]
    if a.shape[0] != 1 or b.shape[0] != 1 or a.shape != b.shape:
        raise SystemExit(f"expected two single-frame exports of equal shape, got {a.shape} and {b.shape}")
    a, b = a[0].astype(np.float64), b[0].astype(np.float64)
    n = len(a)
    ok = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
    res = {"rmsd3": [None] * n, "rmsd4": [None] * n, "n_global": int(ok.sum()),
           "global": float(rmsd(a[ok], b[ok], center=True, superposition=True)) if ok.sum() >= 3 else None,
           "software": {"MDAnalysis": MDAnalysis.__version__, "numpy": np.__version__,
                        "function": "MDAnalysis.analysis.rms.rmsd(center=True, superposition=True)"}}
    for i in range(n):
        if 1 <= i <= n - 2 and ok[i - 1: i + 2].all():
            res["rmsd3"][i] = float(rmsd(a[i - 1: i + 2], b[i - 1: i + 2], center=True, superposition=True))
        if 1 <= i <= n - 3 and ok[i - 1: i + 3].all():
            res["rmsd4"][i] = float(rmsd(a[i - 1: i + 3], b[i - 1: i + 3], center=True, superposition=True))
    with open(out, "w") as fh:
        json.dump(res, fh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
