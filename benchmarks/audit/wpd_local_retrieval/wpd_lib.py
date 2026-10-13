"""WP-D harness library (docs/WPD_LOCAL_RETRIEVAL_PROTOCOL.md §3–§6, §10). Radians, Å.

Harness code, not production. Distances are computed in bounded-memory query batches by callers.
"""
import numpy as np

TIE_REL = 1e-12            # A3
REGRET_ROUND = 1e-10       # A3: [-1e-10, 0) -> 0; below -> numerical failure
NORM_DENOM_MIN = 1e-12     # A3: normalized regret undefined below this
SIN_DEGENERATE = 1e-3      # A1
R_ILL = 0.1                # A1: C1 mean ill-conditioned below this resultant length


def wrap(a):
    w = np.mod(np.asarray(a, float) + np.pi, 2 * np.pi) - np.pi
    return np.where(w <= -np.pi, w + 2 * np.pi, w)


def dihedral(p0, p1, p2, p3):
    """IUPAC dihedral (same formula as WP-A build_traces.dihedral), radians in (−π, π]."""
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
    b1h = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    return np.arctan2(np.sum(np.cross(n1, n2) * b1h, axis=-1), np.sum(n1 * n2, axis=-1))


def sin_angle(p0, p1, p2):
    """sin of the angle at p1 between bonds p1→p0 and p1→p2 (0 for collinear triples)."""
    a, b = p0 - p1, p2 - p1
    return (np.linalg.norm(np.cross(a, b), axis=-1)
            / (np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1)))


def turning_angle(p0, p1, p2):
    t0, t1 = p1 - p0, p2 - p1
    c = np.sum(t0 * t1, axis=-1) / (np.linalg.norm(t0, axis=-1) * np.linalg.norm(t1, axis=-1))
    return np.arccos(np.clip(c, -1.0, 1.0))


# --- tie-aware membership and per-query metrics (§6, A3) -----------------------------------

def membership(X, k):
    """Row-wise expected k-NN membership with uniform selection from the k-th boundary tie set.
    Non-finite entries are never members unless needed; rows with fewer than k finite values are
    flagged invalid. Returns (P, valid_rows, tie_rows)."""
    X = np.asarray(X, float)
    finite = np.isfinite(X)
    valid = finite.sum(1) >= k
    Xf = np.where(finite, X, np.inf)
    dk = np.partition(Xf, k - 1, axis=1)[:, k - 1]
    dk_f = np.where(np.isfinite(dk), dk, 0.0)
    tol = TIE_REL * (1.0 + np.abs(dk_f))
    less = finite & (Xf < (dk_f - tol)[:, None])
    eq = finite & (np.abs(Xf - dk_f[:, None]) <= tol[:, None])
    n_less, n_eq = less.sum(1), eq.sum(1)
    need = k - n_less
    ok = valid & (n_eq >= need) & (need >= 1)
    P = less.astype(float) + eq * np.where(n_eq > 0, need / np.maximum(n_eq, 1), 0.0)[:, None]
    P[~ok] = 0.0
    return P, ok, ok & (n_eq > need)


def query_metrics(R, T, k):
    """Per-query overlap@k, raw regret@k, normalized regret@k (§6) with A3 handling."""
    PR, okR, tieR = membership(R, k)
    PT, okT, _ = membership(T, k)
    ok = okR & okT
    overlap = (PR * PT).sum(1) / k
    retrieved = (PR * T).sum(1) / k
    oracle = (PT * T).sum(1) / k
    random_mean = T.mean(1)
    regret = retrieved - oracle
    rounding = (regret < 0) & (regret >= -REGRET_ROUND)
    failure = ok & (regret < -REGRET_ROUND)
    regret = np.where(rounding, 0.0, regret)
    denom = random_mean - oracle
    norm_def = ok & (denom > NORM_DENOM_MIN)
    norm = np.where(norm_def, regret / np.where(norm_def, denom, 1.0), np.nan)
    return {'valid': ok, 'overlap': overlap, 'regret': regret, 'norm_regret': norm,
            'norm_defined': norm_def, 'random_regret': denom, 'ties': tieR & ok,
            'n_regret_rounded': int(np.sum(rounding & ok)), 'n_regret_failures': int(np.sum(failure)),
            'random_overlap': k / R.shape[1]}


# --- distances (§4) ------------------------------------------------------------------------------

def d_chord(tq, td):
    """M1: sqrt((1/m) Σ ||u(a)-u(b)||²) = sqrt((1/m) Σ 4 sin²(Δ/2)) for (Q,m), (D,m)."""
    d = wrap(tq[:, None, :] - td[None, :, :])
    return np.sqrt(np.mean(4.0 * np.sin(d / 2.0) ** 2, axis=-1))


def d_circ(tq, td):
    """M2 (and the target with φ/ψ stacked): circular RMS."""
    d = wrap(tq[:, None, :] - td[None, :, :])
    return np.sqrt(np.mean(d * d, axis=-1))


def d_rms(vq, vd):
    """B1, B3: RMS difference of feature vectors."""
    d = vq[:, None, :] - vd[None, :, :]
    return np.sqrt(np.mean(d * d, axis=-1))


def d_scalar(mq, md):
    """C1 (A1): |wrap(m_a - m_b)|, +inf for database frames whose mean is ill-conditioned (NaN)."""
    d = np.abs(wrap(mq[:, None] - md[None, :]))
    return np.where(np.isnan(md)[None, :] | np.isnan(mq)[:, None], np.inf, d)


def circular_mean_R(t):
    z = np.mean(np.exp(1j * t), axis=-1)
    return np.angle(z), np.abs(z)


def horn_rmsd(A, B):
    """B2: minimal RMSD over proper rotations (Horn quaternion), A (Q,W,3), B (D,W,3) -> (Q,D)."""
    A = A - A.mean(1, keepdims=True)
    B = B - B.mean(1, keepdims=True)
    W = A.shape[1]
    ga = np.sum(A * A, axis=(1, 2))
    gb = np.sum(B * B, axis=(1, 2))
    M = np.einsum('qwi,dwj->qdij', A, B)
    Sxx, Sxy, Sxz = M[..., 0, 0], M[..., 0, 1], M[..., 0, 2]
    Syx, Syy, Syz = M[..., 1, 0], M[..., 1, 1], M[..., 1, 2]
    Szx, Szy, Szz = M[..., 2, 0], M[..., 2, 1], M[..., 2, 2]
    K = np.stack([
        np.stack([Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx], -1),
        np.stack([Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz], -1),
        np.stack([Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy], -1),
        np.stack([Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz], -1)], -2)
    lam = np.linalg.eigvalsh(K)[..., -1]
    msd = (ga[:, None] + gb[None, :] - 2.0 * lam) / W
    return np.sqrt(np.maximum(msd, 0.0))


def kabsch_rmsd(a, b):
    """Reference (SVD) for one pair, proper rotation."""
    a, b = a - a.mean(0), b - b.mean(0)
    U, _, Vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt(np.mean(np.sum((a @ R.T - b) ** 2, axis=1))))


def ca_pair_vectors(X):
    """B1 features: (F, W, 3) -> (F, P) distances for pairs q - p >= 2."""
    W = X.shape[1]
    pairs = [(p, q) for p in range(W) for q in range(p + 2, W)]
    return np.stack([np.linalg.norm(X[:, q] - X[:, p], axis=-1) for p, q in pairs], axis=1)


# --- rank correlation ---------------------------------------------------------------------------

def rankdata_avg(x):
    x = np.asarray(x, float)
    order = np.argsort(x, kind='mergesort')
    xs = x[order]
    ranks = np.empty(len(x))
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and xs[j + 1] == xs[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks


def spearman(a, b):
    ra, rb = rankdata_avg(a), rankdata_avg(b)
    ra -= ra.mean()
    rb -= rb.mean()
    den = np.sqrt(np.sum(ra * ra) * np.sum(rb * rb))
    return float(np.sum(ra * rb) / den) if den > 0 else None
