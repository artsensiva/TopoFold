"""Independent re-computation of the KRAS K1 crystal endpoint, compared with the recorded endpoint (reference check).

    .venv/bin/python benchmarks/audit/kras_k1/verify_k1_endpoint.py

Uses a different superposition algorithm (Horn's quaternion method, 4x4 eigenproblem) and its own minimal PDB reader instead of
the frozen implementation's SVD-based Kabsch, and compares per-residue between-state, within-state and threshold values and
the endpoint set with results/KRAS-K1_endpoint.json. Reads only the recorded record and the two local crystal files; it does not
modify the record and opens no MD file. Writes results/KRAS-K1_endpoint_verification.json.
"""

from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RECORD = HERE / "results/KRAS-K1_endpoint.json"
PARAMS = HERE / "protocol_params.json"
OUT = HERE / "results/KRAS-K1_endpoint_verification.json"


def read_ca(path: Path) -> dict:
    out = {}
    for ln in path.read_text().splitlines():
        if ln[:6] == "ATOM  " and ln[12:16] == " CA ":
            out.setdefault(ln[21], {})[int(ln[22:26])] = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return out


def horn_rotation(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Rotation R minimising sum |R p_i - q_i|^2 for centred p, q (Horn 1987, quaternion form)."""
    s = p.T @ q
    sxx, sxy, sxz, syx, syy, syz, szx, szy, szz = s.ravel()
    n = np.array([[sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
                  [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
                  [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
                  [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz]])
    vals, vecs = np.linalg.eigh(n)
    w, x, y, z = vecs[:, np.argmax(vals)]
    return np.array([[w*w + x*x - y*y - z*z, 2*(x*y - w*z), 2*(x*z + w*y)],
                     [2*(x*y + w*z), w*w - x*x + y*y - z*z, 2*(y*z - w*x)],
                     [2*(x*z - w*y), 2*(y*z + w*x), w*w - x*x - y*y + z*z]])


def pair_disp(ref: np.ndarray, mov: np.ndarray, resids: np.ndarray, trim_fraction: float) -> np.ndarray:
    def fit(rows):
        pc, qc = mov[rows].mean(0), ref[rows].mean(0)
        r = horn_rotation(mov[rows] - pc, ref[rows] - qc)
        return lambda x: (x - pc) @ r.T + qc
    allrows = np.arange(len(resids))
    d1 = np.linalg.norm(ref - fit(allrows)(mov), axis=1)
    m = int(math.floor(trim_fraction * len(resids) + 1e-9))
    order = np.lexsort((resids, -d1))                 # most displaced first; ties: lower residue number first
    keep = np.setdiff1d(allrows, order[:m])
    return np.linalg.norm(ref - fit(keep)(mov), axis=1)


def main() -> int:
    rec, p = json.loads(RECORD.read_text()), json.loads(PARAMS.read_text())["endpoint"]
    scored = np.array(rec["scored_residues"])
    a = read_ca(ROOT / rec["inputs"]["5US4"]["path"])
    b = read_ca(ROOT / rec["inputs"]["5USJ"]["path"])
    gdp = {c: np.stack([a[c][r] for r in scored]) for c in p["structures"]["GDP"]["chains"]}
    gnp = {c: np.stack([b[c][r] for r in scored]) for c in p["structures"]["GNP"]["chains"]}
    tf = p["pair_fit"]["trim_fraction"]
    between = np.stack([pair_disp(gdp[x], gnp[y], scored, tf) for x in sorted(gdp) for y in sorted(gnp)])
    within = np.stack([pair_disp(g[x], g[y], scored, tf) for g in (gdp, gnp) for x, y in itertools.combinations(sorted(g), 2)])
    bm, wm = between.mean(0), within.max(0)
    thr = np.maximum(p["min_displacement_A"], p["within_state_factor"] * wm)
    pos = scored[bm >= thr - 1e-9].tolist()
    pr = rec["per_residue"]
    diffs = {"b_mean": float(np.abs(bm - np.array(pr["b_mean"])).max()), "w_max": float(np.abs(wm - np.array(pr["w_max"])).max()),
             "threshold": float(np.abs(thr - np.array(pr["threshold"])).max())}
    out = {"id": "KRAS-K1-S1-endpoint-verification", "method": "Horn quaternion superposition and an independent PDB reader",
           "n_pairs": [between.shape[0], within.shape[0]], "max_abs_difference_A": diffs, "E_independent": pos, "E_recorded": rec["E"],
           "E_equal": pos == rec["E"], "tolerance_A": 1e-6,
           "agrees": pos == rec["E"] and max(diffs.values()) < 1e-6, "md_opened": False,
           "numpy": np.__version__, "python": sys.version.split()[0]}
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps({k: out[k] for k in ("n_pairs", "max_abs_difference_A", "E_equal", "agrees")}))
    return 0 if out["agrees"] else 1


if __name__ == "__main__":
    sys.exit(main())
