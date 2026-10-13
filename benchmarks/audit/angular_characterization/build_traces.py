"""Harness trace construction (not production code).

Natural-extension placement of a Cα polyline with constant bond length, prescribed turning
angles κ (angle between consecutive bond vectors) and torsions τ. The torsion convention is
the IUPAC dihedral of four consecutive points; agreement with production `compute_torsions`
is established by the recovery check, not assumed (protocol §3).

SYNTHETIC: these traces are geometric test objects, not proteins.
"""
import numpy as np

from common import BOND


def build_traces(kappa, tau, bond=BOND):
    """kappa: (F, W-2), tau: (F, W-3) radians -> coordinates (F, W, 3) float64."""
    kappa = np.asarray(kappa, dtype=np.float64)
    tau = np.asarray(tau, dtype=np.float64)
    F, nk = kappa.shape
    W = nk + 2
    assert tau.shape == (F, W - 3)
    X = np.zeros((F, W, 3))
    X[:, 1, 0] = bond
    k0 = kappa[:, 0]
    X[:, 2, 0] = bond + bond * np.cos(k0)
    X[:, 2, 1] = bond * np.sin(k0)
    for j in range(3, W):
        a, b, c = X[:, j - 3], X[:, j - 2], X[:, j - 1]
        bc = c - b
        bc /= np.linalg.norm(bc, axis=1, keepdims=True)
        n = np.cross(b - a, bc)
        n /= np.linalg.norm(n, axis=1, keepdims=True)
        m = np.cross(n, bc)
        k = kappa[:, j - 2]          # turning angle at vertex j-1
        t = tau[:, j - 3]            # dihedral (j-3, j-2, j-1, j)
        # new bond direction: turn by k away from bc, rotated by t about bc
        d = (np.cos(k)[:, None] * bc
             + np.sin(k)[:, None] * (np.cos(t)[:, None] * m + np.sin(t)[:, None] * n))
        X[:, j] = c + bond * d
    return X


def dihedral(p0, p1, p2, p3):
    """IUPAC dihedral, independent NumPy reference (radians, (−π, π])."""
    b0 = p1 - p0
    b1 = p2 - p1
    b2 = p3 - p2
    n1 = np.cross(b0, b1)
    n2 = np.cross(b1, b2)
    b1_hat = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    x = np.sum(n1 * n2, axis=-1)
    y = np.sum(np.cross(n1, n2) * b1_hat, axis=-1)
    return np.arctan2(y, x)


def turning_angle(p0, p1, p2):
    t0 = p1 - p0
    t1 = p2 - p1
    c = np.sum(t0 * t1, axis=-1) / (np.linalg.norm(t0, axis=-1) * np.linalg.norm(t1, axis=-1))
    return np.arccos(np.clip(c, -1.0, 1.0))


def recover(X):
    """Independent NumPy recovery of κ (F, W-2) and τ (F, W-3) from coordinates (F, W, 3)."""
    W = X.shape[1]
    kap = np.stack([turning_angle(X[:, j], X[:, j + 1], X[:, j + 2]) for j in range(W - 2)], 1)
    tau = np.stack([dihedral(X[:, j], X[:, j + 1], X[:, j + 2], X[:, j + 3])
                    for j in range(W - 3)], 1)
    return kap, tau
