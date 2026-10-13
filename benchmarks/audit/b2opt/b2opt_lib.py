"""BQ-1 harness library (docs/B2_SINGLE_QUERY_OPTIMIZATION_PROTOCOL.md). Harness code, not production.

Arms (same proper-rotation Cα RMSD objective, exhaustive search, float64, one thread):
- V0: EB-1 harness kernel (eb1_lib.window_dist → wpd_lib.horn_rmsd: centring per call, 4×4 key matrix, eigvalsh).
- V1: database windows prepared once (centred coordinates and G_B); per query: cross-covariance on the prepared
  array, key matrix, eigvalsh.
- V2: as V1, largest key-matrix eigenvalue by Newton–Raphson on the characteristic polynomial (QCP; Theobald 2005,
  Liu et al. 2010; start E0 = (G_A + G_B)/2), eigvalsh fallback for rows that fail the frozen conditions.
Shared: frozen top-k/tie policy (`topk_lex`), ranking-equivalence test (`rank_check`), independent SVD reference
(`ref_svd`).
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'eb1'))
import eb1_lib as E  # noqa: E402  (harness kernel, frozen)

K = 10
W = 8
TOL_D = 1e-6                    # Å, distance agreement
TAU_RANK = 2e-6                 # Å, ranking tolerance
NEWTON_TOL = 1e-11              # relative, on max(λ, 1 Å²)
NEWTON_MAX = 50
GAP_V2B = 1e-4                  # Amendment 1: V2b near-degeneracy fallback, |P'(λ)| ≤ GAP·max(λ, 1)³


# --- top-k and ranking ----------------------------------------------------------------------------------------------

def topk_lex(d, k=K):
    """First k positions in ascending (distance, position) order; exact ties broken by the lower position."""
    n = d.shape[0]
    if n <= k:
        return np.lexsort((np.arange(n), d))
    v = d[np.argpartition(d, k - 1)[:k]].max()
    cand = np.flatnonzero(d <= v)
    return cand[np.lexsort((cand, d[cand]))][:k]


def rank_check(idx, d_ref, k=K, tau=TAU_RANK):
    """(equivalent, strict, near_tie) of an arm's top-k positions against reference distances (protocol §1.2)."""
    kk = min(k + 1, d_ref.shape[0])
    s = np.sort(np.partition(d_ref, kk - 1)[:kk]) if kk < d_ref.shape[0] else np.sort(d_ref)
    kq = min(k, d_ref.shape[0])
    equiv = len(set(idx.tolist())) == kq and bool(np.all(np.abs(d_ref[idx] - s[:kq]) <= tau))
    strict = bool(np.array_equal(idx, topk_lex(d_ref, k)))
    near = bool(np.any(np.diff(s) <= tau))
    return equiv, strict, near


# --- key matrix and characteristic polynomial ------------------------------------------------------------------------

def key_matrix(M):
    """Horn key matrix (n, 4, 4) from cross-covariance M (n, 3, 3), M[i, j] = Σ_w a_wi b_wj (as wpd_lib.horn_rmsd)."""
    Sxx, Sxy, Sxz = M[:, 0, 0], M[:, 0, 1], M[:, 0, 2]
    Syx, Syy, Syz = M[:, 1, 0], M[:, 1, 1], M[:, 1, 2]
    Szx, Szy, Szz = M[:, 2, 0], M[:, 2, 1], M[:, 2, 2]
    return np.stack([
        np.stack([Sxx + Syy + Szz, Syz - Szy, Szx - Sxz, Sxy - Syx], -1),
        np.stack([Syz - Szy, Sxx - Syy - Szz, Sxy + Syx, Szx + Sxz], -1),
        np.stack([Szx - Sxz, Sxy + Syx, -Sxx + Syy - Szz, Syz + Szy], -1),
        np.stack([Sxy - Syx, Szx + Sxz, Syz + Szy, -Sxx - Syy + Szz], -1)], -2)


def det3(M):
    return (M[:, 0, 0] * (M[:, 1, 1] * M[:, 2, 2] - M[:, 1, 2] * M[:, 2, 1])
            - M[:, 0, 1] * (M[:, 1, 0] * M[:, 2, 2] - M[:, 1, 2] * M[:, 2, 0])
            + M[:, 0, 2] * (M[:, 1, 0] * M[:, 2, 1] - M[:, 1, 1] * M[:, 2, 0]))


def det4(A):
    """Determinant of (n, 4, 4) by the 2×2-minor (Laplace) expansion on rows 0–1 and 2–3."""
    def m(r0, r1, c0, c1):
        return A[:, r0, c0] * A[:, r1, c1] - A[:, r0, c1] * A[:, r1, c0]
    s0, s1, s2 = m(0, 1, 0, 1), m(0, 1, 0, 2), m(0, 1, 0, 3)
    s3, s4, s5 = m(0, 1, 1, 2), m(0, 1, 1, 3), m(0, 1, 2, 3)
    c5, c4, c3 = m(2, 3, 2, 3), m(2, 3, 1, 3), m(2, 3, 1, 2)
    c2, c1, c0 = m(2, 3, 0, 3), m(2, 3, 0, 2), m(2, 3, 0, 1)
    return s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0


def charpoly(M, Kmat=None):
    """Coefficients (c2, c1, c0) of P(λ) = λ⁴ + c2 λ² + c1 λ + c0 for the traceless key matrix of M.
    Derived from Newton's identities: c2 = −½ tr K² = −2 Σ M_ij², c1 = −⅓ tr K³ = −8 det M, c0 = det K."""
    Kmat = key_matrix(M) if Kmat is None else Kmat
    return -2.0 * np.einsum('nij,nij->n', M, M), -8.0 * det3(M), det4(Kmat)


def lam_max_qcp(M, E0, tol=NEWTON_TOL, max_iter=NEWTON_MAX, gap=None):
    """Largest key-matrix eigenvalue by Newton–Raphson from E0; returns (λ, fallback mask) with fallback rows
    already replaced by eigvalsh. gap (V2b only): also fall back where |P'(λ)| ≤ gap·max(λ, 1)³, i.e. where λ_max is
    (nearly) a multiple root and the polynomial root has a precision floor near √ε."""
    Kmat = key_matrix(M)
    c2, c1, c0 = charpoly(M, Kmat)
    lam = E0.copy()
    act = np.ones(lam.shape, bool)
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        for _ in range(max_iter):
            l2 = lam * lam
            p = (l2 + c2) * l2 + c1 * lam + c0
            dp = 4.0 * l2 * lam + 2.0 * c2 * lam + c1
            step = np.where(act, p / dp, 0.0)
            lam = lam - step
            act &= ~(np.abs(step) <= tol * np.maximum(np.abs(lam), 1.0))
            if not act.any():
                break
    bad = act | ~np.isfinite(lam) | (lam < 0.0) | (lam > E0 * (1.0 + 1e-12))
    if gap is not None:
        with np.errstate(invalid='ignore', over='ignore'):
            dp = 4.0 * lam ** 3 + 2.0 * c2 * lam + c1
            bad |= ~(np.abs(dp) > gap * np.maximum(np.abs(lam), 1.0) ** 3)
    if bad.any():
        lam[bad] = np.linalg.eigvalsh(Kmat[bad])[:, -1]
    return lam, bad


# --- arms ------------------------------------------------------------------------------------------------------------

def v0_window(X, q, di, s):
    """EB-1 harness kernel for one query position q, window s (database rows di)."""
    return E.window_dist('B2', W, s, X, None, None, np.array([q]), di)[0]


def prepare(X, di, windows):
    """Centred database windows (n_win, n_D, W, 3) and G_B (n_win, n_D), float64, owned and C-contiguous."""
    P = np.empty((len(windows), len(di), W, 3))
    G = np.empty((len(windows), len(di)))
    for i, s in enumerate(windows):
        B = X[di, s:s + W]
        B -= B.mean(1, keepdims=True)
        P[i] = B
        G[i] = np.sum(B * B, axis=(1, 2))
    return P, G


def prepared_window(arm, xq, P_i, G_i):
    """Distances of one query window xq (W, 3) to the prepared database window; returns (d, fallback rows)."""
    a = xq - xq.mean(0)
    ga = float(np.sum(a * a))
    M = np.einsum('wi,dwj->dij', a, P_i)
    if arm == 'V1':
        lam = np.linalg.eigvalsh(key_matrix(M))[:, -1]
        nfb = 0
    else:
        lam, bad = lam_max_qcp(M, 0.5 * (ga + G_i), gap=GAP_V2B if arm == 'V2b' else None)
        nfb = int(bad.sum())
    msd = (ga + G_i - 2.0 * lam) / W
    return np.sqrt(np.maximum(msd, 0.0)), nfb


def search(arm, X, q, di, windows, prep=None, keep_dist=False):
    """Single-query exhaustive top-k over windows. Returns (top-k (n_win, k) positions in di order, distances or None,
    fallback row count)."""
    out = np.empty((len(windows), min(K, len(di))), np.int64)
    D = np.empty((len(windows), len(di))) if keep_dist else None
    nfb = 0
    for i, s in enumerate(windows):
        if arm == 'V0':
            d = v0_window(X, q, di, s)
        else:
            d, f = prepared_window(arm, X[q, s:s + W], prep[0][i], prep[1][i])
            nfb += f
        out[i] = topk_lex(d)
        if keep_dist:
            D[i] = d
    return out, D, nfb


# --- independent reference -------------------------------------------------------------------------------------------

def ref_svd(xq, B):
    """Kabsch via SVD of the cross-covariance with the proper-rotation sign correction. xq (W, 3), B (n, W, 3)."""
    a = xq - xq.mean(0)
    b = B - B.mean(1, keepdims=True)
    H = np.einsum('wi,dwj->dij', a, b)
    U, S, Vt = np.linalg.svd(H)
    sgn = np.sign(np.linalg.det(U) * np.linalg.det(Vt))
    msd = (np.sum(a * a) + np.sum(b * b, axis=(1, 2)) - 2.0 * (S[:, 0] + S[:, 1] + sgn * S[:, 2])) / a.shape[0]
    return np.sqrt(np.maximum(msd, 0.0))


def ref_kabsch_pair(a, b):
    """wpd_lib.kabsch_rmsd: explicit rotation and residuals for one pair."""
    return E.kabsch_rmsd(a, b)


def ref_qcprot_pair(a, b):
    """MDAnalysis lib.qcprot (inputs centred, as that routine requires)."""
    from MDAnalysis.lib import qcprot
    a = np.ascontiguousarray(a - a.mean(0), np.float64)
    b = np.ascontiguousarray(b - b.mean(0), np.float64)
    return float(qcprot.CalcRMSDRotationalMatrix(a, b, a.shape[0], None, None))


# --- gate helpers (SYNTHETIC adversarial suite and transient-memory calibration) ------------------------------------

def _chains(rng, n, L=20, step=3.8):
    v = rng.normal(size=(n, L, 3))
    v *= step / np.linalg.norm(v, axis=-1, keepdims=True)
    return np.cumsum(v, axis=1)


def _rot(rng):
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return Q * np.sign(np.linalg.det(Q))


def _dist(arm, X, q, di, s):
    if arm == 'V0':
        return v0_window(X, q, di, s), 0
    P, G = prepare(X, di, [s])
    return prepared_window(arm, X[q, s:s + W], P[0], G[0])


def synthetic_suite(arm, seed=20261015):
    """Protocol §5.3 on SYNTHETIC chains. Returns {check: {'pass', 'max_err', ...}}."""
    rng = np.random.default_rng(seed)
    out = {}
    base = _chains(rng, 1)[0]
    t = np.linspace(0, 20, 20)[:, None]
    line = t * np.array([[1.0, 2.0, 2.0]]) / 3.0
    plane = np.c_[rng.normal(size=(20, 2)) * 5, np.zeros(20)]
    mirror = base * np.array([1.0, 1.0, -1.0])
    X = np.stack([base, line, plane, mirror, base @ _rot(rng).T + 7.0, base.copy(), base.copy()])
    di = np.arange(len(X))
    err, nfb = 0.0, 0
    for q in range(4):
        for s in range(13):
            d, f = _dist(arm, X, q, di, s)
            err = max(err, float(np.max(np.abs(d - ref_svd(X[q, s:s + W], X[di, s:s + W])))))
            nfb += f
    out['degenerate_collinear_planar_mirror'] = {'pass': err <= TOL_D, 'max_err': err, 'fallback_rows': nfb}
    d, _ = _dist(arm, X, 0, di, 0)
    out['identical_and_exact_duplicates'] = {
        'pass': bool(d[0] <= TOL_D and d[5] == d[6] == d[0]
                     and np.array_equal(topk_lex(d, 4), np.lexsort((np.arange(len(d)), d))[:4])
                     and topk_lex(d, 1)[0] == int(np.flatnonzero(d == d.min()).min())),
        'd_identical': float(d[0])}
    Y = _chains(rng, 300)
    dY = np.arange(1, 300)
    ok_perm, e_perm, e_rig, ok_near = True, 0.0, 0.0, True
    for s in (0, 6, 12):
        d, _ = _dist(arm, Y, 0, dY, s)
        perm = rng.permutation(len(dY))
        dp, _ = _dist(arm, Y, 0, dY[perm], s)
        e_perm = max(e_perm, float(np.max(np.abs(dp - d[perm]))))
        ok_perm &= rank_check(perm[topk_lex(dp)], d)[0]
        Yr = Y.copy()
        Yr[0] = Y[0] @ _rot(rng).T + np.array([3.0, -40.0, 11.0])
        e_rig = max(e_rig, float(np.max(np.abs(_dist(arm, Yr, 0, dY, s)[0] - d))))
        Yn = np.concatenate([Y, Y[topk_lex(d) + 1] + 1e-9 * rng.normal(size=(K, Y.shape[1], 3))])
        dn, _ = _dist(arm, Yn, 0, np.arange(1, len(Yn)), s)
        ok_near &= rank_check(topk_lex(dn), ref_svd(Yn[0, s:s + W], Yn[1:, s:s + W]))[0]
    out['permutation'] = {'pass': bool(ok_perm and e_perm <= 1e-9), 'max_err': e_perm}
    out['rigid_transform'] = {'pass': e_rig <= 1e-8, 'max_err': e_rig}
    out['near_ties'] = {'pass': bool(ok_near)}
    out['pass'] = all(v['pass'] for v in out.values())
    return out


def calibrate_alpha(arm, n_d=(1000, 2000), seed=20261016):
    """Peak traced transient bytes per database row of one query window (tracemalloc), SYNTHETIC."""
    import gc
    import tracemalloc
    rng = np.random.default_rng(seed)
    per = []
    for n in n_d:
        X = _chains(rng, n + 1, L=12)
        di = np.arange(1, n + 1)
        prep = None if arm == 'V0' else prepare(X, di, [0])
        search(arm, X, 0, di, [0], prep)
        gc.collect()
        tracemalloc.start()
        o, _, _ = search(arm, X, 0, di, [0], prep)
        per.append((tracemalloc.get_traced_memory()[1] - o.nbytes) / n)
        tracemalloc.stop()
    return {'alpha_bytes_per_row': float(max(per)), 'per_n_d': dict(zip(map(str, n_d), per))}
