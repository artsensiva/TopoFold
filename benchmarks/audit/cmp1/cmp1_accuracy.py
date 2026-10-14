"""CMP-1 accuracy tracks (protocol §4). Pure functions plus the adapter-level degeneracy suite.

Reference: b2opt_lib.ref_svd on the canonical float64 Å input (independent SVD Kabsch, proper rotation).

Track F64 (strict): |d − d_ref| ≤ 1e-6 Å for every pair; ranking equivalent (BQ-1 §1.2) with τ = 2e-6 Å.
Track F32 (common relaxed): |d − d_ref| ≤ 1e-3 Å for every pair; ranking equivalent with τ = 2e-3 Å. 1e-3 Å is one
tenth of the coordinate quantization of the inputs (AdK XTC precision 1000 = 0.001 nm; mdCATH float32 Å values on a
0.01 Å grid), so it is set by input resolution, not by any implementation's error.
Identical requirements apply to every implementation within a track. Native-precision characterization (descriptive
only): max |d − d_ref| and the msd ratio |d² − d_ref²| / (u·(G_q + G_j)/W), u = 2⁻²⁴.
"""
import numpy as np

import cmp1_adapters as A

B = A.B
W = A.W
K = A.K
U32 = 2.0 ** -24
TOL = {'F64': 1e-6, 'F32': 1e-3}
TAU = {'F64': 2e-6, 'F32': 2e-3}


def window_G(X, rows, s):
    b = X[np.asarray(rows), s:s + W]
    b = b - b.mean(1, keepdims=True)
    return np.sum(b * b, axis=(1, 2))


def ref_distances(X, di, wins, q):
    return np.stack([B.ref_svd(X[q, s:s + W], X[np.asarray(di), s:s + W]) for s in wins])


def native_ratio(d, d_ref, Gq, Gdb):
    """Descriptive: max |d² − d_ref²| in units of u·(G_q + G_j)/W."""
    return float(np.max(np.abs(d ** 2 - d_ref ** 2) / (U32 * (Gq + Gdb) / W)))


def check_window(track, d, d_ref, idx, Gq=None, Gdb=None):
    """Returns (distance_ok, max_abs_err_A, rank_equivalent, strict_identical, near_tie)."""
    finite = bool(np.isfinite(d).all())
    err = float(np.max(np.abs(d - d_ref))) if finite else float('inf')
    eq, st, nr = B.rank_check(idx, d_ref, tau=TAU[track])
    return finite and err <= TOL[track], err, eq, st, nr


def evaluate(track, D, D_ref, TI, Gq_list, Gdb_list):
    """D, D_ref: (n_win, n_D); TI: (n_win, k). Aggregates check_window over windows."""
    ok, mx, n_bad, n_st, n_nr, nat = True, 0.0, 0, 0, 0, 0.0
    for i in range(D.shape[0]):
        dok, err, eq, st, nr = check_window(track, D[i], D_ref[i], TI[i], Gq_list[i], Gdb_list[i])
        if np.isfinite(D[i]).all():
            nat = max(nat, native_ratio(D[i], D_ref[i], Gq_list[i], Gdb_list[i]))
        ok &= dok
        mx = max(mx, err)
        n_bad += not eq
        n_st += st
        n_nr += nr
    return {'distance_ok': bool(ok), 'max_abs_err_A': mx, 'native_msd_ratio_descriptive': nat,
            'rank_nonequivalent': n_bad, 'strict_identical': n_st, 'near_tie': n_nr, 'windows': D.shape[0],
            'pass': bool(ok and n_bad == 0)}


def _run(impl, X, di, wins, q):
    h = impl.prepare(X, di, wins)
    D = impl.distances(h, X, q)
    TI = np.stack([B.topk_lex(D[i]) for i in range(D.shape[0])])
    return D, TI


def _eval_case(impl, X, di, wins, q):
    D, TI = _run(impl, X, di, wins, q)
    Dr = ref_distances(X, di, wins, q)
    Gq = [window_G(X, [q], s)[0] for s in wins]
    Gd = [window_G(X, di, s) for s in wins]
    return D, TI, {t: evaluate(t, D, Dr, TI, Gq, Gd) for t in ('F64', 'F32')}


def degeneracy_suite(impl, seed):
    """SYNTHETIC: collinear, planar, mirror, identical, exact duplicates, permutation, rigid transform, near ties."""
    rng = np.random.default_rng(seed)

    def chains(n, L=20):
        v = rng.normal(size=(n, L, 3))
        v *= 3.8 / np.linalg.norm(v, axis=-1, keepdims=True)
        return np.cumsum(v, axis=1)

    def rot():
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        return Q * np.sign(np.linalg.det(Q))

    def f32(Z):                                                   # canonical: float32-representable float64
        Z = np.ascontiguousarray(np.asarray(Z, np.float32).astype(np.float64))
        Z.setflags(write=False)
        return Z
    out = {t: {} for t in ('F64', 'F32')}
    base = chains(1)[0]
    line = np.linspace(0, 20, 20)[:, None] * np.array([[1.0, 2.0, 2.0]]) / 3.0
    plane = np.c_[rng.normal(size=(20, 2)) * 5, np.zeros(20)]
    X = f32(np.stack([base, line, plane, base * np.array([1.0, 1.0, -1.0]), base @ rot().T + 7.0, base, base]))
    di, wins = np.arange(7), list(range(13))
    names = {0: 'query_identical_to_database_rows', 1: 'collinear_query', 2: 'planar_query', 3: 'mirror_image_query'}
    for q, nm in names.items():
        _, _, ev = _eval_case(impl, X, di, wins, q)
        for t in out:
            out[t][nm] = bool(ev[t]['pass'])
            out[t][nm + '_max_abs_err_A'] = ev[t]['max_abs_err_A']
    D, TI, ev = _eval_case(impl, X, di, wins, 0)
    for t in out:
        out[t]['exact_duplicates_bitwise_equal_and_position_order'] = bool(
            np.array_equal(D[:, 5], D[:, 6]) and all(TI[i].tolist() == np.lexsort((np.arange(7), D[i]))[:7].tolist()
                                                    for i in range(len(wins))))
    Y = f32(chains(301))
    dY = np.arange(1, 301)
    perm = rng.permutation(300)
    D1, _, ev1 = _eval_case(impl, Y, dY, [0, 6, 12], 0)
    D2, TI2, ev2 = _eval_case(impl, Y, dY[perm], [0, 6, 12], 0)
    Yr = np.array(Y)
    Yr[0] = Y[0] @ rot().T + np.array([3.0, -40.0, 11.0])
    Yr = f32(Yr)
    D3, _, ev3 = _eval_case(impl, Yr, dY, [0, 6, 12], 0)
    tw = np.array(Y)
    tw = f32(np.concatenate([tw, tw[1:11] + 1e-9 * rng.normal(size=(10, 20, 3))]))
    _, _, ev4 = _eval_case(impl, tw, np.arange(1, len(tw)), [0, 6, 12], 0)
    for t in out:
        out[t]['permutation'] = bool(ev1[t]['pass'] and ev2[t]['pass'])
        out[t]['rigid_transform'] = bool(ev3[t]['pass'])
        out[t]['near_ties'] = bool(ev4[t]['pass'])
        out[t]['pass'] = all(v for k, v in out[t].items() if k != 'pass' and isinstance(v, bool))
    out['descriptive'] = {'max_abs_diff_permuted_vs_original_A': float(np.max(np.abs(D2 - D1[:, perm]))),
                          'max_abs_diff_rigid_vs_original_A': float(np.max(np.abs(D3 - D1)))}
    return out
