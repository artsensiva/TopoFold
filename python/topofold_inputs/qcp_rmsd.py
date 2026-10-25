"""Per-window and global proper-rotation Cα RMSD with MDAnalysis QCP (run in an interpreter with MDAnalysis).

    python -m topofold_inputs.qcp_rmsd A.npz B.npz OUT.json
    python -m topofold_inputs.qcp_rmsd --batch PAIRS.json OUT.json     # PAIRS.json: {"pairs": [{"id", "a", "b"}, ...]}

Reads two single-frame adapter exports with the same residue order. Windows of 3 and 4 consecutive output rows, indexed
by the κ / τ-hinge convention (window i = rows i−1 .. i+1 or i−1 .. i+2); a window is computed only where all its rows
have finite coordinates in both. The global value uses every row finite in both. Same calculation as
``benchmarks/audit/kras_struct_cmp/rmsd_windows.py`` (MDAnalysis.analysis.rms.rmsd, center and superposition).
"""

import json
import sys

import numpy as np


def _software():
    import MDAnalysis
    return {"MDAnalysis": MDAnalysis.__version__, "numpy": np.__version__,
            "function": "MDAnalysis.analysis.rms.rmsd(center=True, superposition=True)"}


def _pair_result(a, b, software) -> dict:
    """Window and global values of one pair of (n, 3) coordinate arrays (the single calculation of this module)."""
    from MDAnalysis.analysis.rms import rmsd

    a, b = a.astype(np.float64), b.astype(np.float64)
    n = len(a)
    ok = np.isfinite(a).all(axis=1) & np.isfinite(b).all(axis=1)
    res = {"rmsd3": [None] * n, "rmsd4": [None] * n, "n_global": int(ok.sum()),
           "global": float(rmsd(a[ok], b[ok], center=True, superposition=True)) if ok.sum() >= 3 else None,
           "software": software}
    for i in range(n):
        if 1 <= i <= n - 2 and ok[i - 1: i + 2].all():
            res["rmsd3"][i] = float(rmsd(a[i - 1: i + 2], b[i - 1: i + 2], center=True, superposition=True))
        if 1 <= i <= n - 3 and ok[i - 1: i + 3].all():
            res["rmsd4"][i] = float(rmsd(a[i - 1: i + 3], b[i - 1: i + 3], center=True, superposition=True))
    return res


def main(argv=None) -> int:
    args = argv or sys.argv[1:]
    if args and args[0] == "--batch":
        spec = json.load(open(args[1]))
        software = _software()
        results = {}
        for p in spec["pairs"]:
            a, b = np.load(p["a"])["coordinates_angstrom"], np.load(p["b"])["coordinates_angstrom"]
            if a.shape[0] != 1 or b.shape[0] != 1 or a.shape != b.shape:
                raise SystemExit(f"pair {p['id']}: expected two single-frame exports of equal shape, got {a.shape} and {b.shape}")
            results[p["id"]] = _pair_result(a[0], b[0], software)
        with open(args[2], "w") as fh:
            json.dump({"results": results, "software": software}, fh)
        return 0
    a_path, b_path, out = args
    a = np.load(a_path)["coordinates_angstrom"]
    b = np.load(b_path)["coordinates_angstrom"]
    if a.shape[0] != 1 or b.shape[0] != 1 or a.shape != b.shape:
        raise SystemExit(f"expected two single-frame exports of equal shape, got {a.shape} and {b.shape}")
    res = _pair_result(a[0], b[0], _software())
    with open(out, "w") as fh:
        json.dump(res, fh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
