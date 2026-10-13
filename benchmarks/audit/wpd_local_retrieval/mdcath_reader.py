"""Minimal mdCATH reader for the WP-D held-out tier (feasibility; Amendment 2 §5). Harness code, not production.

Conventions (verified on the schema sample, results/wpd_mdcath_schema_check.json):
- atom names come from the `pdbProteinAtoms` PDB text (columns 13–16); `resid`/`resname` datasets must agree;
- backbone residues = residues with exactly one N, CA, C (and O for T_pp); terminal caps are CHARMM patch atoms
  inside the terminal residues (CAY/CY/OY, NT/CAT), not separate residues;
- coordinates in Å (`coords.attrs['unit']`), float32 as stored; `box` is 3×3, read as diag × 10 (nm → Å); a
  non-diagonal box is an error;
- frame order is storage order; there are no timestamps;
- coordinates are read only for the selected backbone atoms, in bounded frame chunks.
"""
from pathlib import Path

import h5py
import numpy as np

BACKBONE = ('N', 'CA', 'C', 'O')
NM_TO_A = 10.0


class ReaderError(RuntimeError):
    pass


def _s(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def parse_pdb_atoms(text):
    rows = [ln for ln in text.split('\n') if ln.startswith(('ATOM', 'HETATM'))]
    return ([ln[12:16].strip() for ln in rows], [ln[17:21].strip() for ln in rows], [int(ln[22:26]) for ln in rows])


class MdcathDomain:
    def __init__(self, path, domain=None):
        self.path = Path(path)
        self.f = h5py.File(self.path, 'r')
        if _s(self.f.attrs.get('layout', '')) != 'mdcath-only-protein-v1.0':
            raise ReaderError(f"unexpected layout {self.f.attrs.get('layout')!r}")
        keys = list(self.f.keys())
        if domain is None:
            if len(keys) != 1:
                raise ReaderError(f'expected one domain group, found {keys}')
            domain = keys[0]
        self.domain = domain
        self.g = self.f[domain]
        names, rn, ri = parse_pdb_atoms(_s(self.g['pdbProteinAtoms'][()]))
        R = self.g['resid'][:]
        RN = [_s(x) for x in self.g['resname'][:]]
        if len(names) != len(R) or not np.array_equal(np.asarray(ri), R) or rn != RN:
            raise ReaderError('pdbProteinAtoms does not match resid/resname datasets')
        if len(names) != int(self.g.attrs['numProteinAtoms']):
            raise ReaderError('atom count differs from numProteinAtoms')
        self.atom_names, self.resid, self.resname = names, R, RN
        self.mapping = self._map()

    def _map(self):
        resids = sorted(set(self.resid.tolist()))
        idx, missing, dup = {a: [] for a in BACKBONE}, {a: [] for a in BACKBONE}, {a: [] for a in BACKBONE}
        for r in resids:
            sel = np.where(self.resid == r)[0]
            for a in BACKBONE:
                hit = [k for k in sel if self.atom_names[k] == a]
                if len(hit) == 0:
                    missing[a].append(r)
                    idx[a].append(-1)
                elif len(hit) > 1:
                    dup[a].append(r)
                    idx[a].append(-1)
                else:
                    idx[a].append(hit[0])
        bb = [r for k, r in enumerate(resids) if all(idx[a][k] >= 0 for a in ('N', 'CA', 'C'))]
        contiguous = bool(bb) and bb == list(range(bb[0], bb[0] + len(bb)))
        return {'resids': resids, 'index': {a: np.asarray(v) for a, v in idx.items()}, 'missing': missing,
                'duplicated': dup, 'backbone_resids': bb, 'backbone_contiguous': contiguous}

    def replicas(self):
        out = {}
        for t in self.g.keys():
            if t.isdigit():
                out[t] = {r: int(self.g[t][r].attrs['numFrames']) for r in self.g[t].keys()}
        return out

    def box_A(self, temp, rep):
        b = np.asarray(self.g[temp][rep]['box'][:], float)
        if np.any(np.abs(b - np.diag(np.diag(b))) > 0):
            raise ReaderError('non-diagonal box')
        return np.diag(b) * NM_TO_A

    def backbone(self, temp, rep, atoms=BACKBONE, chunk=64):
        """(F, L, len(atoms), 3) float64 holding the stored float32 values, read in frame chunks of `chunk`
        and only for the requested atoms. Residues lacking an atom raise unless the atom is 'O'."""
        ds = self.g[temp][rep]['coords']
        if _s(ds.attrs.get('unit', '')) != 'Angstrom':
            raise ReaderError(f"coords unit {ds.attrs.get('unit')!r}")
        nF = int(self.g[temp][rep].attrs['numFrames'])
        if ds.shape[0] != nF:
            raise ReaderError('numFrames differs from coords length')
        cols = np.stack([self.mapping['index'][a] for a in atoms], axis=1)       # (L, A)
        if np.any(cols < 0):
            raise ReaderError('missing or duplicated backbone atom')
        flat = cols.ravel()
        order = np.argsort(flat)
        inv = np.empty_like(order)
        inv[order] = np.arange(len(order))
        sorted_idx = flat[order]
        if np.any(np.diff(sorted_idx) == 0):
            raise ReaderError('atom index used twice')
        out = np.empty((nF, cols.shape[0], cols.shape[1], 3))
        for f0 in range(0, nF, chunk):
            blk = ds[f0:f0 + chunk, sorted_idx, :]            # h5py needs increasing indices
            out[f0:f0 + blk.shape[0]] = blk[:, inv, :].astype(np.float64).reshape(blk.shape[0], *cols.shape, 3)
        return out

    def close(self):
        self.f.close()


def make_whole_backbone(bb, box):
    """Sequential minimum image along N→CA→C→N(+1); O follows its own C. bb: (F, L, 4, 3) N, CA, C, O."""
    F, L = bb.shape[:2]
    chain = bb[:, :, :3].reshape(F, 3 * L, 3).copy()
    for i in range(1, 3 * L):
        chain[:, i] -= box * np.round((chain[:, i] - chain[:, i - 1]) / box)
    out = bb.copy()
    out[:, :, :3] = chain.reshape(F, L, 3, 3)
    C = out[:, :, 2]
    out[:, :, 3] = bb[:, :, 3] - box * np.round((bb[:, :, 3] - C) / box)
    return out
