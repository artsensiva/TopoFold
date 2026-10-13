"""Research implementations for the legacy κ–τ engine evaluation (docs/LEGACY_ENGINE_EVALUATION_PROTOCOL.md).

Harness code, not production. Radians and Å. Distances take query features (Q, …) and database features
(D, …) and return (Q, D). The legacy Fréchet mirrors topofold-core `discrete_frechet_invariants`,
including its NaN semantics (Rust f64::max/min ignore a NaN operand; numpy fmax/fmin do the same) and its
pairing (sequence length min(#κ, #τ)); it is checked against the Rust path, not assumed equal to it.
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WPD = HERE.parent / 'wpd_local_retrieval'
sys.path.insert(0, str(WPD))
from wpd_lib import (ca_pair_vectors, d_chord, d_circ, d_rms, horn_rmsd, kabsch_rmsd,  # noqa: E402,F401
                     query_metrics, wrap)


def d_kappa(kq, kd):
    """K: RMS |Δκ| over the W − 2 turning angles."""
    return d_rms(kq, kd)


def d_aligned(kq, tq, kd, td):
    """A: √(½·mean Δκ² + ½·mean wrap(Δτ)²) over all W − 2 κ and W − 3 τ (channel-count normalized)."""
    dk = kq[:, None, :] - kd[None, :, :]
    dt = wrap(tq[:, None, :] - td[None, :, :])
    return np.sqrt(0.5 * np.mean(dk * dk, axis=-1) + 0.5 * np.mean(dt * dt, axis=-1))


def legacy_pairs(k, t):
    """Legacy element sequence: (κ[j], τ[j]) for j < min(#κ, #τ)."""
    m = min(k.shape[-1], t.shape[-1])
    return k[..., :m], t[..., :m]


def pointwise(kq, tq, kd, td):
    """P[q, d, i, j] = √(Δκ² + wrap(Δτ)²) between query element i and database element j."""
    dk = np.abs(kq[:, None, :, None] - kd[None, :, None, :])
    dt = np.abs(wrap(tq[:, None, :, None] - td[None, :, None, :]))
    return np.sqrt(dk * dk + dt * dt)


def frechet_legacy(kq, tq, kd, td):
    """Discrete Fréchet over legacy pairs, same recurrence and NaN semantics as the Rust code.
    Empty sequences return 0.0 (as the Rust code does)."""
    kq, tq = legacy_pairs(kq, tq)
    kd, td = legacy_pairs(kd, td)
    Q, D = kq.shape[0], kd.shape[0]
    la, lb = kq.shape[-1], kd.shape[-1]
    if la == 0 or lb == 0:
        return np.zeros((Q, D))
    P = pointwise(kq, tq, kd, td)
    ca = np.empty((Q, D, la, lb))
    ca[..., 0, 0] = P[..., 0, 0]
    for j in range(1, lb):
        ca[..., 0, j] = np.fmax(ca[..., 0, j - 1], P[..., 0, j])
    for i in range(1, la):
        ca[..., i, 0] = np.fmax(ca[..., i - 1, 0], P[..., i, 0])
        for j in range(1, lb):
            pred = np.fmin(np.fmin(ca[..., i - 1, j], ca[..., i, j - 1]), ca[..., i - 1, j - 1])
            ca[..., i, j] = np.fmax(pred, P[..., i, j])
    return ca[..., la - 1, lb - 1]


def frechet_lockstep(kq, tq, kd, td):
    """F_lock: max over the legacy pairs of the pointwise distance on the diagonal (no warping)."""
    kq, tq = legacy_pairs(kq, tq)
    kd, td = legacy_pairs(kd, td)
    dk = kq[:, None, :] - kd[None, :, :]
    dt = wrap(tq[:, None, :] - td[None, :, :])
    return np.fmax.reduce(np.sqrt(dk * dk + dt * dt), axis=-1)


def peptide_plane_atoms(BB, O, s, W):
    """T_pp atoms for window start s (0-based): C_r, O_r, N_{r+1} for r = s … s+W−2 → (F, 3(W−1), 3)."""
    r = np.arange(s, s + W - 1)
    C = BB[:, r, 2]
    Ox = O[:, r]
    Nn = BB[:, r + 1, 0]
    return np.stack([C, Ox, Nn], axis=2).reshape(BB.shape[0], 3 * (W - 1), 3)
