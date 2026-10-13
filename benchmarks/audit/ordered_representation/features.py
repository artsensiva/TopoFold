"""WP-C candidate: residue-ordered (cos τ, sin τ) features and same-sequence distances.

Harness code (SYNTHETIC evaluation), not production. Radians throughout; Cα distances in Å.
Protocol: docs/ORDERED_ANGULAR_REPRESENTATION_PROTOCOL.md §2–§3, §5a.
"""
import math

import numpy as np


def window_torsions(chain_torsions, i, W):
    """Torsions τ_i … τ_{i+W-4} of window i (W residues); τ_j uses Cα j…j+3."""
    return np.asarray(chain_torsions)[..., i:i + W - 3]


def bond_labels(i, W):
    """Virtual bond (j+1, j+2) to which each window torsion τ_j is attributed."""
    return [(j + 1, j + 2) for j in range(i, i + W - 3)]


def ordered_features(t):
    """(…, m) torsions -> (…, 2m) features (cos τ_1, sin τ_1, …, cos τ_m, sin τ_m)."""
    t = np.asarray(t, dtype=np.float64)
    out = np.empty(t.shape[:-1] + (2 * t.shape[-1],))
    out[..., 0::2] = np.cos(t)
    out[..., 1::2] = np.sin(t)
    return out


def wrap(a):
    w = np.mod(np.asarray(a, dtype=np.float64) + np.pi, 2 * np.pi) - np.pi
    return np.where(w <= -np.pi, w + 2 * np.pi, w)


def _pairwise_rms(A, B, diff, chunk=256):
    """RMS over the last axis of diff(A[i], B[j]) for all i, j, computed directly (no Gram trick)."""
    out = np.empty((A.shape[0], B.shape[0]))
    for s in range(0, A.shape[0], chunk):
        d = diff(A[s:s + chunk, None, :], B[None, :, :])
        out[s:s + chunk] = np.sqrt(np.mean(d * d, axis=-1))
    return out


def d_chord(ta, tb):
    """D_ch = ||x_a - x_b|| / sqrt(m) for torsion arrays (Fa, m), (Fb, m)."""
    m = ta.shape[-1]
    xa, xb = ordered_features(ta), ordered_features(tb)
    return _pairwise_rms(xa, xb, lambda a, b: a - b) * math.sqrt(2.0)   # mean over 2m -> /m


def d_circ(ta, tb):
    """D_circ = sqrt((1/m) Σ wrap(τa - τb)^2)."""
    return _pairwise_rms(np.asarray(ta, float), np.asarray(tb, float), lambda a, b: wrap(a - b))


def d_scalar(ma, mb):
    """Historical comparator: |wrap(m_a - m_b)| of per-frame circular means."""
    return np.abs(wrap(np.asarray(ma)[:, None] - np.asarray(mb)[None, :]))


def ca_distance_vectors(X):
    """(F, W, 3) -> (F, P) Cα–Cα distances for all pairs q - p >= 2 within the window (Å)."""
    W = X.shape[1]
    pairs = [(p, q) for p in range(W) for q in range(p + 2, W)]
    return np.stack([np.linalg.norm(X[:, q] - X[:, p], axis=-1) for p, q in pairs], axis=1), pairs


def d_euclid_rms(A, B):
    return _pairwise_rms(np.asarray(A, float), np.asarray(B, float), lambda a, b: a - b)


# --- feature-contract checks (protocol §5a) -------------------------------------------------------

FC1_TAU = (0.1, 0.7, -2.0, 3.0, -0.4)


def feature_contracts(feature_fn=ordered_features, window_fn=window_torsions, bond_fn=bond_labels):
    """Direct convention checks. Functions are parameters so corruptions can be substituted."""
    checks = []
    x = np.asarray(feature_fn(np.array([FC1_TAU])))[0]
    expect = []
    for t in FC1_TAU:
        expect += [math.cos(t), math.sin(t)]
    checks.append(('FC1 values (asymmetric known torsions, order, cos/sin entries)',
                   x.shape == (10,) and bool(np.all(np.abs(x - np.array(expect)) <= 1e-15)),
                   {'got': [float(v) for v in x], 'expected': expect}))
    chain = np.array([0.1 * (j + 1) for j in range(9)])          # 12-residue chain: 9 torsions
    xw = np.asarray(feature_fn(window_fn(chain, 3, 6)[None, :]))[0]
    exp2 = [math.cos(0.4), math.sin(0.4), math.cos(0.5), math.sin(0.5), math.cos(0.6), math.sin(0.6)]
    checks.append(('FC2 window offset (i=3, W=6 uses τ_3, τ_4, τ_5)',
                   xw.shape == (6,) and bool(np.all(np.abs(xw - np.array(exp2)) <= 1e-15)),
                   {'got': [float(v) for v in xw], 'expected': exp2}))
    bl = bond_fn(3, 6)
    checks.append(('FC3 residue mapping (bonds (4,5), (5,6), (6,7))', bl == [(4, 5), (5, 6), (6, 7)],
                   {'got': bl}))
    xs = np.asarray(feature_fn(np.array([[-0.7]])))[0]
    checks.append(('FC4 chirality (sine entry of τ = -0.7 is -sin 0.7)',
                   abs(xs[1] + math.sin(0.7)) <= 1e-15 and abs(xs[0] - math.cos(0.7)) <= 1e-15,
                   {'got': [float(v) for v in xs]}))
    return checks
