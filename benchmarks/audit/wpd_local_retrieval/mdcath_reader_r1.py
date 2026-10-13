"""mdCATH reader, revision 1 (Amendment 3 §2, §5). Supersedes mdcath_reader.py for the held-out runner; the original
file is preserved unmodified (hash in results/wpd_heldout_amendment2.json).

Changes:
- N, CA, C are mandatory; O is optional. A missing or duplicated O makes only that residue's O unavailable (NaN,
  `o_present[r] = False`); it never affects N/CA/C data or the primary φ/ψ target.
- Terminal caps are CHARMM patch atoms inside the terminal residues (verified on the schema sample).

Memory scaling (accurate statement): `backbone()` returns float64 arrays of size F × L × 3 × 3 (N, CA, C) plus
F × L × 3 (O), i.e. linear in frames × residues of ONE trajectory; trajectories are read one at a time. Reads go in
frame chunks over the selected atoms only, so the transient buffer is chunk × (4L) × 3 float32 in addition to the
output. The whole-file size (all atoms incl. hydrogens, forces, 25 trajectories) never enters memory.
"""
import numpy as np

from mdcath_reader import NM_TO_A, MdcathDomain as _Base, ReaderError, _s  # noqa: F401

MANDATORY = ('N', 'CA', 'C')


class MdcathDomainR1(_Base):
    def mandatory_ok(self):
        m = self.mapping
        return not any(m['missing'][a] or m['duplicated'][a] for a in MANDATORY)

    def backbone_r1(self, temp, rep, chunk=64):
        """Returns (bb, O, o_present): bb (F, L, 3, 3) N/CA/C float64 from stored float32; O (F, L, 3) with NaN rows
        for residues without a unique O; o_present (L,) bool. Raises if any residue lacks a unique N, CA or C."""
        if not self.mandatory_ok():
            raise ReaderError('missing or duplicated mandatory backbone atom (N, CA or C)')
        ds = self.g[temp][rep]['coords']
        if _s(ds.attrs.get('unit', '')) != 'Angstrom':
            raise ReaderError(f"coords unit {ds.attrs.get('unit')!r}")
        nF = int(self.g[temp][rep].attrs['numFrames'])
        if ds.shape[0] != nF:
            raise ReaderError('numFrames differs from coords length')
        idx = self.mapping['index']
        o_present = idx['O'] >= 0
        cols = np.concatenate([np.stack([idx[a] for a in MANDATORY], 1).ravel(), idx['O'][o_present]])
        order = np.argsort(cols)
        inv = np.empty_like(order)
        inv[order] = np.arange(len(order))
        srt = cols[order]
        if np.any(np.diff(srt) == 0):
            raise ReaderError('atom index used twice')
        L = len(o_present)
        bb = np.empty((nF, L, 3, 3))
        O = np.full((nF, L, 3), np.nan)
        nb = 3 * L
        for f0 in range(0, nF, chunk):
            blk = ds[f0:f0 + chunk, srt, :][:, inv, :].astype(np.float64)
            n = blk.shape[0]
            bb[f0:f0 + n] = blk[:, :nb].reshape(n, L, 3, 3)
            O[f0:f0 + n, o_present] = blk[:, nb:]
        return bb, O, o_present


def make_whole_r1(bb, O, box):
    """Sequential minimum image along N→CA→C→N(+1) (as WP-B make_whole); O follows its own C; NaN O stays NaN."""
    F, L = bb.shape[:2]
    chain = bb.reshape(F, 3 * L, 3).copy()
    for i in range(1, 3 * L):
        chain[:, i] -= box * np.round((chain[:, i] - chain[:, i - 1]) / box)
    whole = chain.reshape(F, L, 3, 3)
    C = whole[:, :, 2]
    Ow = O - box * np.round((O - C) / box)
    return whole, Ow
