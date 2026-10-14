"""CMP-1 adapters (docs/CMP1_COMPETITIVE_EVALUATION_PROTOCOL.md §3–§5). Harness code, not production.

Canonical physical input: a read-only float64 array X (rows, L, 3) in Å (S2 rounded values or SYNTHETIC tiled rows
recast once to float32-representable values). Each adapter converts explicitly, never mutates X, and exposes:

    prepare(X, di, windows) -> handle          database-side preparation (timed as preparation)
    query(handle, X, q) -> (topk_idx (n_win, k) int64, topk_d (n_win, k) float64 Å, kernel_s, nfb)
    distances(handle, X, q) -> (n_win, n_D) float64 Å     (verification only; not timed)

Every adapter passes its per-window distance vector (cast to float64, in Å) through the same `topk_lex`.
"""
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'b2opt'))
import b2opt_lib as B  # noqa: E402  (V2b, V0 harness, topk_lex, ref_svd; unchanged)

W = 8
K = 10
NM_PER_A = 0.1


def _finish(d64):
    idx = B.topk_lex(d64)
    return idx, d64[idx]


# --- T-V2b: unchanged BQ harness ---------------------------------------------------------------------------------------

class V2b:
    name, precision = 'T-V2b', 'float64'

    def prepare(self, X, di, wins):
        return {'prep': B.prepare(X, di, wins), 'wins': wins}

    def _dist(self, h, X, q, i):
        return B.prepared_window('V2b', X[q, h['wins'][i]:h['wins'][i] + W], h['prep'][0][i], h['prep'][1][i])

    def query(self, h, X, q):
        n = len(h['wins'])
        ti, td = np.empty((n, K), np.int64), np.empty((n, K))
        ks, nfb = 0.0, 0
        for i in range(n):
            t0 = time.perf_counter()
            d, f = self._dist(h, X, q, i)
            ks += time.perf_counter() - t0
            nfb += f
            ti[i], td[i] = _finish(d)
        return ti, td, ks, nfb

    def distances(self, h, X, q):
        return np.stack([self._dist(h, X, q, i)[0] for i in range(len(h['wins']))])

    def prepared_bytes(self, h):
        return int(h['prep'][0].nbytes + h['prep'][1].nbytes)


# --- T-V0: EB-1 harness kernel (no preparation) ------------------------------------------------------------------------

class V0:
    name, precision = 'T-V0', 'float64'

    def prepare(self, X, di, wins):
        return {'di': np.asarray(di), 'wins': wins}

    def query(self, h, X, q):
        n = len(h['wins'])
        ti, td = np.empty((n, K), np.int64), np.empty((n, K))
        ks = 0.0
        for i, s in enumerate(h['wins']):
            t0 = time.perf_counter()
            d = B.v0_window(X, q, h['di'], s)
            ks += time.perf_counter() - t0
            ti[i], td[i] = _finish(d)
        return ti, td, ks, 0

    def distances(self, h, X, q):
        return np.stack([B.v0_window(X, q, h['di'], s) for s in h['wins']])

    def prepared_bytes(self, h):
        return 0


# --- E-MDT: MDTraj 1.11.2, pre-sliced windows with cached centering -------------------------------------------------

class MDT:
    """Database: one 8-atom Trajectory per window (atom_slice of a full nm float32 trajectory), centred once with
    center_coordinates(); md.rmsd(..., precentered=True, parallel=False) is valid only with atom_indices=None
    (verified in the pinned source), which pre-slicing guarantees. Query: per window, an 8-atom 1-frame Trajectory on
    the prepared window topology, centred, then md.rmsd. Distances are float32 nm → float64 Å (×10)."""
    name, precision = 'E-MDT', 'float32'

    def __init__(self):
        import mdtraj as md
        self.md = md

    def to_nm32(self, A):
        return (np.asarray(A, np.float64) * NM_PER_A).astype(np.float32)

    def prepare(self, X, di, wins):
        md = self.md
        L = X.shape[1]
        top = md.Topology()
        ch = top.add_chain()
        for r in range(L):
            res = top.add_residue('ALA', ch)
            top.add_atom('CA', md.element.carbon, res)
        full = md.Trajectory(self.to_nm32(X[np.asarray(di)]), top)
        slices = []
        for s in wins:
            t = full.atom_slice(np.arange(s, s + W))
            t.center_coordinates()
            slices.append(t)
        del full
        return {'slices': slices, 'wins': wins}

    def _qtraj(self, h, X, q, i):
        s = h['wins'][i]
        t = self.md.Trajectory(self.to_nm32(X[q, s:s + W])[None], h['slices'][i].topology)
        t.center_coordinates()
        return t

    def _dist(self, h, X, q, i):
        qt = self._qtraj(h, X, q, i)
        t0 = time.perf_counter()
        d = self.md.rmsd(h['slices'][i], qt, 0, precentered=True, parallel=False)
        return d.astype(np.float64) * 10.0, time.perf_counter() - t0

    def query(self, h, X, q):
        n = len(h['wins'])
        ti, td = np.empty((n, K), np.int64), np.empty((n, K))
        ks = 0.0
        for i in range(n):
            d, k = self._dist(h, X, q, i)
            ks += k
            ti[i], td[i] = _finish(d)
        return ti, td, ks, 0

    def distances(self, h, X, q):
        return np.stack([self._dist(h, X, q, i)[0] for i in range(len(h['wins']))])

    def prepared_bytes(self, h):
        return int(sum(t.xyz.nbytes + (t._rmsd_traces.nbytes if t._rmsd_traces is not None else 0)
                       for t in h['slices']))


# --- E-MDA: MDAnalysis 2.10.0 qcprot, independently fitted windows, prepared float64 path --------------------------

class MDA:
    """Database: per window, centred float64 coordinates (n_D, 8, 3) in Å (qcprot requires centred inputs, verified in
    the pinned source). Query: per window, centre the query window, then CalcRMSDRotationalMatrix(query, db_j, 8,
    None, None) for every database row (independent proper-rotation fit per pair; QCP internally in double)."""
    name, precision = 'E-MDA', 'float64'

    def __init__(self):
        from MDAnalysis.lib import qcprot
        self.qcp = qcprot.CalcRMSDRotationalMatrix

    def prepare(self, X, di, wins):
        P = np.empty((len(wins), len(di), W, 3))
        for i, s in enumerate(wins):
            b = X[np.asarray(di), s:s + W].astype(np.float64)
            P[i] = b - b.mean(1, keepdims=True)
        return {'P': P, 'wins': wins}

    def _dist(self, h, X, q, i):
        s = h['wins'][i]
        a = np.array(X[q, s:s + W], np.float64)
        a -= a.mean(0)
        Pi = h['P'][i]
        f = self.qcp
        t0 = time.perf_counter()
        d = np.fromiter((f(a, Pi[j], W, None, None) for j in range(Pi.shape[0])), np.float64, Pi.shape[0])
        return d, time.perf_counter() - t0

    def query(self, h, X, q):
        n = len(h['wins'])
        ti, td = np.empty((n, K), np.int64), np.empty((n, K))
        ks = 0.0
        for i in range(n):
            d, k = self._dist(h, X, q, i)
            ks += k
            ti[i], td[i] = _finish(d)
        return ti, td, ks, 0

    def distances(self, h, X, q):
        return np.stack([self._dist(h, X, q, i)[0] for i in range(len(h['wins']))])

    def prepared_bytes(self, h):
        return int(h['P'].nbytes)


IMPLS = {'T-V2b': V2b, 'T-V0': V0, 'E-MDT': MDT, 'E-MDA': MDA}
