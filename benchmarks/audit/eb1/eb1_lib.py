"""EB-1 stage implementations (docs/EB1_ENGINEERING_BENCHMARK_PROTOCOL.md). Harness code, not production.

Only existing implementations are used: MDAnalysis XTC reader, frozen mdcath_reader_r1 (h5py), WP-B make_whole,
the Rust binding compute_invariants, legacy_lib.d_aligned (A), wpd_lib.horn_rmsd (B2), wpd_lib.membership.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
AUDIT = HERE.parent
REPO = AUDIT.parents[1]
WPD = AUDIT / 'wpd_local_retrieval'
sys.path.insert(0, str(WPD))
sys.path.insert(0, str(AUDIT / 'legacy_engine'))
from legacy_lib import d_aligned  # noqa: E402
from wpd_lib import dihedral, horn_rmsd, kabsch_rmsd, membership, turning_angle, wrap  # noqa: E402

_spec = importlib.util.spec_from_file_location('wpb_adk_lib', AUDIT / 'adk_feasibility' / 'adk_lib.py')
ADK = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ADK)

RAW = AUDIT / 'adk_feasibility' / 'data' / 'raw'
D1_GRO, D1_XTC = RAW / 'ORADD_open_state.gro', RAW / 'ORADD_open_state_traj.xtc'
D1_DERIVED = AUDIT / 'adk_feasibility' / 'data' / 'derived' / 'open_start_ca_whole.npy'
D2_H5 = WPD / 'data' / 'mdcath' / 'mdcath_dataset_1ha8A00.h5'
D2_TRAJ = (('320', '0'), ('320', '1'), ('348', '0'), ('348', '1'))
D1_DB = (1238, 2253)
D1_Q = (0, 1013)
K = 10
BATCH = 128
TIE_REL = 1e-9


# --- S1 ingestion -------------------------------------------------------------------------------------------------

def ingest_d1(n_frames):
    """XTC → stored float32 Cα (n, 214, 3) in ascending resid order, plus per-frame box (n, 3) Å."""
    import MDAnalysis as mda
    u = mda.Universe(str(D1_GRO), str(D1_XTC))
    ca = u.select_atoms('protein and name CA')
    ca = ca[np.argsort(ca.resids, kind='stable')]
    out = np.empty((n_frames, ca.n_atoms, 3), np.float32)
    box = np.empty((n_frames, 3))
    for k, ts in enumerate(u.trajectory[:n_frames]):
        out[k] = ca.positions
        box[k] = ts.dimensions[:3]
    return out, box, {'atoms_decoded_per_frame': int(u.trajectory.n_atoms), 'ca_indices': ca.indices.tolist()}


def ingest_d2(n_traj):
    """HDF5 → stored float32 Cα per trajectory (frozen reader; reads N, CA, C, O), plus box (3,) Å."""
    from mdcath_reader_r1 import MdcathDomainR1
    d = MdcathDomainR1(D2_H5)
    out = []
    for t, r in D2_TRAJ[:n_traj]:
        bb, _O, _op = d.backbone_r1(t, r)
        out.append((bb[:, :, 1].astype(np.float32), d.box_A(t, r)))
    d.close()
    return out


# --- S2 reconstruction --------------------------------------------------------------------------------------------

def reconstruct(ca32, box, burn_fraction=0.0):
    """Per-frame sequential minimum image along the Cα chain (WP-B make_whole), burn-in slice, float32 rounding.
    Returns (whole float64 before rounding, rounded float64)."""
    F = ca32.shape[0]
    whole = np.empty(ca32.shape, np.float64)
    for f in range(F):
        b = box[f] if np.ndim(box) == 2 else box
        whole[f] = ADK.make_whole(ca32[f].astype(np.float64), b)
    b0 = int(np.floor(burn_fraction * F))
    whole = whole[b0:]
    return whole, whole.astype(np.float32).astype(np.float64)


# --- S3 features -------------------------------------------------------------------------------------------------

def features(ca_rounded):
    """Binding κ (F, L−2) and τ (F, L−3) on float32 Cα."""
    import topofold as tf
    F, L = ca_rounded.shape[:2]
    kap, tau = np.empty((F, L - 2)), np.empty((F, L - 3))
    x32 = ca_rounded.astype(np.float32)
    for f in range(F):
        k_, t_, _ = tf.compute_invariants(x32[f])
        kap[f], tau[f] = k_, t_
    return kap, tau


def features_numpy(ca_rounded):
    X = ca_rounded
    return (turning_angle(X[:, :-2], X[:, 1:-1], X[:, 2:]),
            dihedral(X[:, :-3], X[:, 1:-2], X[:, 2:-1], X[:, 3:]))


# --- S4 search ---------------------------------------------------------------------------------------------------

def window_dist(method, W, s, X, kap, tau, qi, di):
    if method == 'A':
        return d_aligned(kap[qi, s:s + W - 2], tau[qi, s:s + W - 3], kap[di, s:s + W - 2], tau[di, s:s + W - 3])
    return horn_rmsd(X[qi, s:s + W], X[di, s:s + W])


def topk(Dm, k=K):
    idx = np.argpartition(Dm, k - 1, axis=1)[:, :k]
    order = np.take_along_axis(Dm, idx, 1).argsort(1, kind='stable')
    return np.take_along_axis(idx, order, 1)


def search(method, W, windows, X, kap, tau, qi, di, batch=BATCH):
    """Exhaustive top-k for every query in qi, every window start in `windows`. Returns {s: (n_q, k) indices}."""
    out = {}
    for s in windows:
        res = []
        for b in range(0, len(qi), batch):
            res.append(topk(window_dist(method, W, s, X, kap, tau, qi[b:b + batch], di)))
        out[s] = np.vstack(res)
    return out


# --- correctness helpers -----------------------------------------------------------------------------------------

def reference_distances(method, W, s, X, kap, tau, q, di):
    """Per-pair recomputation: A via d_aligned on single pairs; B2 via Kabsch SVD per pair."""
    if method == 'A':
        return np.array([d_aligned(kap[[q], s:s + W - 2], tau[[q], s:s + W - 3], kap[[d], s:s + W - 2],
                                   tau[[d], s:s + W - 3])[0, 0] for d in di])
    return np.array([kabsch_rmsd(X[q, s:s + W], X[d, s:s + W]) for d in di])


def compare_topk(op_dist, ref_dist, k=K):
    """Returns (ok, near_tie, max_abs_diff). Near tie: d_(k+1) − d_(k) ≤ 1e-9 (1 + d_(k)) in the reference."""
    srt = np.sort(ref_dist)
    near_tie = bool(srt[k] - srt[k - 1] <= TIE_REL * (1 + srt[k - 1]))
    diff = float(np.max(np.abs(op_dist - ref_dist)))
    if near_tie:
        Po, _, _ = membership(op_dist[None, :], k)
        Pr, _, _ = membership(ref_dist[None, :], k)
        ok = bool(np.max(np.abs(Po - Pr)) <= 1e-12)
    else:
        ok = set(np.argpartition(op_dist, k - 1)[:k].tolist()) == set(np.argpartition(ref_dist, k - 1)[:k].tolist())
    return ok, near_tie, diff
