"""SYNTHETIC verification of the legacy κ–τ engine and the research methods (protocol §6).

Fixtures are constructed geometric objects, not proteins; they verify mechanisms, not biological usefulness.
Tests marked `rust` call the actual Rust path through the `topofold` binding (ConformationalIndex.query_subcurve,
compute_invariants).
Run: ../wpd_local_retrieval/build/venv/bin/python -B -m pytest -p no:cacheprovider test_legacy.py
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / 'angular_characterization'))
from build_traces import build_traces, recover  # noqa: E402
from legacy_lib import (ca_pair_vectors, d_aligned, d_chord, d_circ, d_kappa, d_rms,  # noqa: E402
                        frechet_legacy, frechet_lockstep, horn_rmsd, kabsch_rmsd)

tf = pytest.importorskip('topofold')
BOND = 3.8


def rot(seed):
    q = np.random.default_rng(seed).normal(size=4)
    q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])


def chain(kappa, tau):
    return build_traces(np.atleast_2d(kappa), np.atleast_2d(tau), bond=BOND)


def all_distances(Xq, Xd):
    """Every Cα method on windows Xq (Q, W, 3), Xd (D, W, 3); features from the NumPy reference."""
    kq, tq = recover(Xq)
    kd, td = recover(Xd)
    return {'M1': d_chord(tq, td), 'M2': d_circ(tq, td), 'K': d_kappa(kq, kd), 'A': d_aligned(kq, tq, kd, td),
            'F': frechet_legacy(kq, tq, kd, td), 'F_lock': frechet_lockstep(kq, tq, kd, td),
            'B1': d_rms(ca_pair_vectors(Xq), ca_pair_vectors(Xd)), 'B2': horn_rmsd(Xq, Xd)}


def rust_frechet(db_frames, start, end, query):
    """Actual legacy path: exhaustive query_subcurve over all indexed frames; returns distances by frame id."""
    idx = tf.ConformationalIndex(window_size=4, coarse_radius=0.5)
    idx.fit(np.ascontiguousarray(db_frames, dtype=np.float32))
    hits = idx.query_subcurve(start, end, np.ascontiguousarray(query, dtype=np.float32), k=len(db_frames))
    out = np.empty(len(db_frames))
    for fid, d in hits:
        out[fid] = d
    return out


R = np.random.default_rng(11)
K8 = R.uniform(1.0, 1.9, 6)            # W = 8: 6 turning angles
T8 = R.uniform(-np.pi, np.pi, 5)       # 5 torsions


# --- rigid motion, independent per structure ---------------------------------------------------

def test_independent_rigid_motions_leave_every_method_unchanged():
    Xa = chain(K8, T8)[0]
    Xb = chain(K8 + 0.05, T8 + 0.3)[0]
    ref = all_distances(Xa[None], Xb[None])
    Ya = Xa @ rot(1).T + np.array([10.0, -3.0, 7.0])
    Yb = Xb @ rot(2).T + np.array([-40.0, 2.0, 0.5])
    moved = all_distances(Ya[None], Yb[None])
    for m in ref:
        assert abs(ref[m][0, 0] - moved[m][0, 0]) < 1e-9, m


# --- torsion wrap --------------------------------------------------------------------------------

def test_torsion_wrap_boundary():
    ta, tb = T8.copy(), T8.copy()
    ta[2], tb[2] = math.pi - 0.01, -math.pi + 0.01
    d = all_distances(chain(K8, ta), chain(K8, tb))
    assert d['M2'][0, 0] == pytest.approx(0.02 / math.sqrt(5), abs=1e-9)
    assert d['F'][0, 0] == pytest.approx(0.02, abs=1e-9)          # pair 2 is a legacy pair
    assert d['A'][0, 0] == pytest.approx(math.sqrt(0.5 * 0.02 ** 2 / 5), abs=1e-9)


# --- channel-specific perturbations; boundary κ omission ---------------------------------------

def test_kappa_only_change_and_boundary_omission():
    k_first, k_last = K8.copy(), K8.copy()
    k_first[0] += 0.2
    k_last[-1] += 0.2                      # vertex W−2: omitted by the legacy pairing
    base = chain(K8, T8)
    d_first = all_distances(base, chain(k_first, T8))
    d_last = all_distances(base, chain(k_last, T8))
    for d in (d_first, d_last):
        assert d['M1'][0, 0] < 1e-9 and d['M2'][0, 0] < 1e-9
        assert d['K'][0, 0] == pytest.approx(0.2 / math.sqrt(6), abs=1e-9)
        assert d['A'][0, 0] == pytest.approx(math.sqrt(0.5 * 0.04 / 6), abs=1e-9)
    assert d_first['F'][0, 0] == pytest.approx(0.2, abs=1e-9)
    assert d_last['F'][0, 0] < 1e-12 and d_last['F_lock'][0, 0] < 1e-12


def test_rust_path_ignores_last_window_kappa():
    """Real engine: a change confined to the last interior turning angle of the window is invisible to F."""
    N, W, s = 14, 8, 3
    k = R.uniform(1.0, 1.9, N - 2)
    t = R.uniform(-np.pi, np.pi, N - 3)
    k_last, k_first = k.copy(), k.copy()
    k_last[s + W - 3] += 0.3          # production κ index of vertex s+W−2
    k_first[s] += 0.3                 # vertex s+1
    db = np.concatenate([chain(k, t), chain(k_last, t), chain(k_first, t)])
    f = rust_frechet(db, s, s + W - 1, db[0])
    assert f[0] < 1e-5 and f[1] < 1e-5 and f[2] == pytest.approx(0.3, abs=1e-4)


def test_tau_only_change():
    t2 = T8.copy()
    t2[1] += 0.4
    d = all_distances(chain(K8, T8), chain(K8, t2))
    assert d['K'][0, 0] < 1e-9 and d['M2'][0, 0] == pytest.approx(0.4 / math.sqrt(5), abs=1e-9)
    assert d['F'][0, 0] == pytest.approx(0.4, abs=1e-9)


# --- scalar-mean collision -------------------------------------------------------------------------

def test_scalar_mean_collision_resolved_by_ordered_methods():
    ta = np.array([0.5, -1.0, 2.0, 0.1, -0.7])
    tb = ta[::-1].copy()                  # same multiset → identical circular mean
    za, zb = np.mean(np.exp(1j * ta)), np.mean(np.exp(1j * tb))
    assert abs(np.angle(za) - np.angle(zb)) < 1e-12
    d = all_distances(chain(K8, ta), chain(K8, tb))
    assert d['M1'][0, 0] > 0.5 and d['A'][0, 0] > 0.3 and d['F'][0, 0] > 0.3


# --- reflection ----------------------------------------------------------------------------------

def test_reflection_behaviour():
    Xa = chain(K8, T8)
    Xm = Xa * np.array([1.0, 1.0, -1.0])
    ka, ta = recover(Xa)
    km, tm = recover(Xm)
    assert np.allclose(km, ka, atol=1e-12) and np.allclose(tm, -ta, atol=1e-12)
    d = all_distances(Xa, Xm)
    assert d['B1'][0, 0] < 1e-9                       # distances cannot see chirality
    assert d['K'][0, 0] < 1e-9                        # κ cannot see chirality
    assert d['B2'][0, 0] > 0.5 and d['M1'][0, 0] > 0.5 and d['F'][0, 0] > 0.5


# --- unequal lengths ---------------------------------------------------------------------------------

def test_unequal_lengths_follow_min_pairing():
    kq, tq = K8[None, :], T8[None, :3]                 # 6 κ, 3 τ → 3 legacy pairs
    kd = K8[None, :] + np.array([0, 0, 0, 1.0, 1.0, 1.0])
    assert frechet_legacy(kq, tq, kd, T8[None, :3])[0, 0] == 0.0


# --- repeated, collinear, near-collinear -------------------------------------------------------------

def test_degenerate_inputs_in_binding_and_fit():
    X = chain(K8, T8)[0].astype(np.float32)
    rep = X.copy()
    rep[4] = rep[3]
    with pytest.raises(Exception):
        tf.compute_invariants(rep)
    col = X.copy()
    col[3] = np.round(col[3] * 64) / 64               # exactly representable, exactly collinear 3,4,5
    col[4] = col[3] + np.array([2.0, 0.0, 0.0], np.float32)
    col[5] = col[3] + np.array([4.0, 0.0, 0.0], np.float32)
    with pytest.raises(Exception):
        tf.compute_invariants(col)
    idx = tf.ConformationalIndex(window_size=4, coarse_radius=0.5)
    with pytest.raises(Exception):                    # one collinear frame aborts the whole fit (m1)
        idx.fit(np.stack([X, col]))


def test_float32_rounded_collinear_passes_binding_check_silently():
    """Documented defect (LE-D6): GEOMETRY_EPSILON = 1e-9 on |T_i × T_{i+1}| is below float32 coordinate
    rounding, so a collinear triple built in float32 arithmetic is accepted and yields arbitrary torsions."""
    X = chain(K8, T8)[0].astype(np.float32)
    col = X.copy()
    col[5] = col[4] + (col[4] - col[3])
    k, t, _ = tf.compute_invariants(col)              # no error
    assert np.asarray(k)[3] < 1e-6                    # vertex 4 turning angle ~ rounding noise
    assert np.all(np.isfinite(np.asarray(t)))         # finite, but τ about bonds 3–4 and 4–5 is noise


def test_near_collinear_torsion_is_unstable():
    k = K8.copy()
    k[2] = 1e-7                                       # turning angle ~1e-7 rad at vertex 3
    X = chain(k, T8)[0]
    t0 = recover(X[None])[1][0]
    Xp = X + np.random.default_rng(3).normal(scale=1e-6, size=X.shape)
    t1 = recover(Xp[None])[1][0]
    jump = np.abs(np.angle(np.exp(1j * (t1 - t0))))
    assert jump[1] > 0.1 or jump[2] > 0.1             # τ about the neighbouring bonds is ill-defined
    assert jump[4] < 1e-3                             # far torsion unaffected


# --- float32 vs float64 ----------------------------------------------------------------------------

def test_float32_rounding_effect_on_features_is_small_for_regular_geometry():
    X = chain(K8, T8)[0] + 20.0
    k64, t64 = recover(X[None])
    k32, t32 = recover(X.astype(np.float32).astype(np.float64)[None])
    dt = np.abs(np.angle(np.exp(1j * (t32 - t64)))).max()
    dk = np.abs(k32 - k64).max()
    assert dt < 1e-5 and dk < 1e-5
    kr, tr, _ = tf.compute_invariants(X.astype(np.float32))
    assert np.abs(np.asarray(kr) - k32[0]).max() < 1e-9
    assert np.abs(np.angle(np.exp(1j * (np.asarray(tr) - t32[0])))).max() < 1e-9


# --- Fréchet misregistration -------------------------------------------------------------------------

def test_frechet_misregistration_numpy_and_rust():
    """Query τ pattern [x, y, y, y, y] vs candidate [x, x, y, y, y]: different conformations at residue
    position 1, but the minimax coupling maps x→(x, x) and returns 0."""
    k = np.full(6, 1.5)
    tq = np.array([0.3, 2.5, 2.5, 2.5, 2.5])
    tc = np.array([0.3, 0.3, 2.5, 2.5, 2.5])
    d = all_distances(chain(k, tq), chain(k, tc))
    assert d['F'][0, 0] < 1e-9
    assert d['F_lock'][0, 0] > 2.0 and d['M2'][0, 0] > 0.9 and d['A'][0, 0] > 0.6
    db = np.concatenate([chain(k, tq), chain(k, tc)])
    f = rust_frechet(db, 0, 7, db[0])
    assert f[1] < 1e-5                                # the real engine treats them as identical


# --- empty and invalid inputs ------------------------------------------------------------------------

def test_empty_and_nan_semantics_mirror_rust():
    z = np.zeros((1, 0))
    assert frechet_legacy(z, z, K8[None], T8[None])[0, 0] == 0.0          # empty matches everything
    tn = T8.copy()
    tn[1] = np.nan
    assert frechet_legacy(K8[None], T8[None], K8[None], tn[None])[0, 0] == 0.0  # NaN pair ignored
    assert np.isnan(frechet_legacy(K8[None], T8[None], K8[None], np.full((1, 5), np.nan))[0, 0])
    with pytest.raises(Exception):
        tf.ConformationalIndex(window_size=4).fit(np.zeros((1, 3, 3), np.float32))


# --- research Fréchet vs Rust path on random windows -----------------------------------------------

def test_numpy_frechet_matches_rust_path_on_random_chains():
    F, N, W, s = 40, 16, 8, 4
    r = np.random.default_rng(21)
    db = np.concatenate([chain(r.uniform(0.6, 2.2, N - 2), r.uniform(-np.pi, np.pi, N - 3)) for _ in range(F)])
    db32 = db.astype(np.float32)
    f_rust = rust_frechet(db32, s, s + W - 1, db32[7])
    kk, tt = [], []
    for f in range(F):
        kf, tf_, _ = tf.compute_invariants(db32[f])
        kk.append(np.asarray(kf)[s:s + W - 2])
        tt.append(np.asarray(tf_)[s:s + W - 3])
    kk, tt = np.array(kk), np.array(tt)
    f_np = frechet_legacy(kk[7:8], tt[7:8], kk, tt)[0]
    assert np.max(np.abs(f_np - f_rust)) < 1e-9


def test_horn_kabsch_and_frechet_symmetry():
    r = np.random.default_rng(5)
    A, B = r.normal(size=(4, 8, 3)) * 3, r.normal(size=(6, 8, 3)) * 3
    H = horn_rmsd(A, B)
    assert max(abs(H[i, j] - kabsch_rmsd(A[i], B[j])) for i in range(4) for j in range(6)) < 1e-6
    ka, ta = r.uniform(0.5, 2, (5, 6)), r.uniform(-3, 3, (5, 5))
    assert np.allclose(frechet_legacy(ka, ta, ka, ta), frechet_legacy(ka, ta, ka, ta).T, atol=1e-15)
