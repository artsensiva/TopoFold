"""Kinase-state descriptors from named atoms (engineering capability; definitions only, no classification).

Definitions follow Modi & Dunbrack, *PNAS* 116:6818–6827 (2019), doi:10.1073/pnas.1814279116:

* ``D1`` = dist(Cα of αC-Glu+4, Cζ of DFG-Phe); ``D2`` = dist(Cα of β3-Lys, Cζ of DFG-Phe)   (Å)
* ``alphaC_CbCb`` = dist(Cβ of αC-Glu, Cβ of β3-Lys); the paper calls ≤ 10 Å "C-helix-in"      (Å)
* ``Aloop_NO`` = dist(N of DFGxxX, O of X-HRD), the extended-activation-loop backbone H-bond      (Å)
* ``chi1_Phe`` = N–CA–CB–CG of DFG-Phe; ``phi``/``psi`` of X-DFG, DFG-Asp, DFG-Phe               (rad, IUPAC sign)

Physics (AGENTS rule 7): D1/D2 and χ1 locate the Phe ring and need side-chain atoms; Cα-only geometry cannot see a
χ1 rotation. φ/ψ need backbone N/C. The αC Cβ–Cβ distance reports whether the Glu–Lys salt bridge can form, not that it
does. No threshold or state label is applied here; the module returns values only.

Every function returns an array over frames with NaN where a needed atom is missing (never zero), and
:func:`kinase_descriptors` reports which atom caused it.
"""

from __future__ import annotations

from typing import Any, Mapping

import numpy as np

from .atom_adapter import AtomTrajectory
from .ca_adapter import ResidueSpec

DESCRIPTOR_UNITS = {"D1": "angstrom", "D2": "angstrom", "alphaC_CbCb": "angstrom", "Aloop_NO": "angstrom",
                    "chi1_Phe": "rad", "phi_XDFG": "rad", "psi_XDFG": "rad", "phi_Asp": "rad", "psi_Asp": "rad",
                    "phi_Phe": "rad", "psi_Phe": "rad"}

# descriptor -> list of (role, atom) in the order the formula uses them
DESCRIPTOR_ATOMS: dict[str, list[tuple[str, str]]] = {
    "D1": [("alphaC_Glu_plus4", "CA"), ("DFG_Phe", "CZ")],
    "D2": [("beta3_Lys", "CA"), ("DFG_Phe", "CZ")],
    "alphaC_CbCb": [("alphaC_Glu", "CB"), ("beta3_Lys", "CB")],
    "Aloop_NO": [("DFGxxX", "N"), ("XHRD", "O")],
    "chi1_Phe": [("DFG_Phe", "N"), ("DFG_Phe", "CA"), ("DFG_Phe", "CB"), ("DFG_Phe", "CG")],
    "phi_XDFG": [("XDFG_minus1", "C"), ("XDFG", "N"), ("XDFG", "CA"), ("XDFG", "C")],
    "psi_XDFG": [("XDFG", "N"), ("XDFG", "CA"), ("XDFG", "C"), ("DFG_Asp", "N")],
    "phi_Asp": [("XDFG", "C"), ("DFG_Asp", "N"), ("DFG_Asp", "CA"), ("DFG_Asp", "C")],
    "psi_Asp": [("DFG_Asp", "N"), ("DFG_Asp", "CA"), ("DFG_Asp", "C"), ("DFG_Phe", "N")],
    "phi_Phe": [("DFG_Asp", "C"), ("DFG_Phe", "N"), ("DFG_Phe", "CA"), ("DFG_Phe", "C")],
    "psi_Phe": [("DFG_Phe", "N"), ("DFG_Phe", "CA"), ("DFG_Phe", "C"), ("DFG_Gly", "N")],
}


def distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Euclidean distance per frame; NaN propagates."""
    return np.linalg.norm(np.asarray(b, float) - np.asarray(a, float), axis=-1)


def dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    """IUPAC torsion angle p0-p1-p2-p3 in (-π, π] per frame (positive = clockwise looking from p1 to p2).
    Degenerate (collinear) geometry or missing atoms give NaN."""
    p0, p1, p2, p3 = (np.atleast_2d(np.asarray(x, float)) for x in (p0, p1, p2, p3))
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    n1 = np.linalg.norm(b1, axis=-1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        b1n = b1 / n1
        v = b0 - np.sum(b0 * b1n, axis=-1, keepdims=True) * b1n
        w = b2 - np.sum(b2 * b1n, axis=-1, keepdims=True) * b1n
        x = np.sum(v * w, axis=-1)
        y = np.sum(np.cross(b1n, v) * w, axis=-1)
        ang = np.arctan2(y, x)
    degenerate = (np.linalg.norm(v, axis=-1) < 1e-9) | (np.linalg.norm(w, axis=-1) < 1e-9) | (n1[..., 0] < 1e-9)
    ang = np.where(degenerate, np.nan, ang)
    return np.where(ang <= -np.pi, np.pi, ang)


def kinase_descriptors(traj: AtomTrajectory, roles: Mapping[str, ResidueSpec | Mapping[str, Any]]) -> dict[str, Any]:
    """All descriptors per frame for one role -> residue mapping. Returns
    {"values": {name: (F,) array}, "units": ..., "missing": {name: [(role, atom, status)]}}."""
    spec = {r: (s if isinstance(s, ResidueSpec) else ResidueSpec.from_mapping(s)) for r, s in roles.items()}
    status = {(a["chain"], a["resid"], a["icode"], a["name"]): a["status"] for a in traj.atoms}

    def get(role: str, atom: str):
        s = spec[role]
        return traj.xyz(s.chain, s.resid, atom, s.icode), status.get((s.chain, s.resid, s.icode, atom), "not_requested")

    values, missing = {}, {}
    for name, needs in DESCRIPTOR_ATOMS.items():
        if any(r not in spec for r, _ in needs):
            values[name] = np.full(len(traj.frame_index), np.nan)
            missing[name] = [(r, a, "role_not_mapped") for r, a in needs if r not in spec]
            continue
        pts, bad = [], []
        for r, a in needs:
            try:
                xyz, st = get(r, a)
            except KeyError:
                xyz, st = np.full((len(traj.frame_index), 3), np.nan), "not_requested"
            pts.append(xyz)
            if st != "present":
                bad.append((r, a, st))
        values[name] = distance(*pts) if len(pts) == 2 else dihedral(*pts)
        if bad:
            values[name] = np.full(len(traj.frame_index), np.nan)
            missing[name] = bad
    return {"values": values, "units": DESCRIPTOR_UNITS, "missing": missing,
            "definitions": "Modi & Dunbrack, PNAS 116:6818 (2019); see module docstring"}
