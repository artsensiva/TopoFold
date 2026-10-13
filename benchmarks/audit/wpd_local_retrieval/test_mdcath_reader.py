"""SYNTHETIC tests for the minimal mdCATH reader (constructed HDF5 files; no mdCATH data).

Run: build/venv/bin/python -B -m pytest -p no:cacheprovider test_mdcath_reader.py
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mdcath_reader import MdcathDomain, ReaderError, make_whole_backbone  # noqa: E402


def write_fixture(path, drop=None, offdiag=False, unit='Angstrom'):
    """SYNTHETIC 3-residue 'domain' with a cap atom on residue 1, shuffled atom order within residues."""
    atoms = [(1, 'GLY', 'CAY'), (1, 'GLY', 'N'), (1, 'GLY', 'CA'), (1, 'GLY', 'O'), (1, 'GLY', 'C'),
             (2, 'ALA', 'N'), (2, 'ALA', 'CB'), (2, 'ALA', 'CA'), (2, 'ALA', 'C'), (2, 'ALA', 'O'),
             (3, 'CYS', 'C'), (3, 'CYS', 'N'), (3, 'CYS', 'CA'), (3, 'CYS', 'O'), (3, 'CYS', 'NT')]
    if drop:
        atoms = [a for a in atoms if a[:1] + a[2:] != drop]
    lines = ['MODEL        1'] + [f"ATOM  {k + 1:5d} {n:<4s} {rn:<4s}0{r:4d}    {0:8.3f}{0:8.3f}{0:8.3f}"
                                   for k, (r, rn, n) in enumerate(atoms)] + ['ENDMDL', 'END']
    nA, F = len(atoms), 5
    with h5py.File(path, 'w') as f:
        f.attrs['layout'] = 'mdcath-only-protein-v1.0'
        g = f.create_group('x00')
        g.attrs['numProteinAtoms'] = nA
        g['pdbProteinAtoms'] = '\n'.join(lines).encode()
        g['resid'] = np.array([a[0] for a in atoms])
        g['resname'] = np.array([a[1].encode() for a in atoms], dtype=object)
        rep = g.create_group('320').create_group('0')
        rep.attrs['numFrames'] = F
        X = np.arange(F * nA * 3, dtype=np.float32).reshape(F, nA, 3)
        d = rep.create_dataset('coords', data=X)
        d.attrs['unit'] = unit
        b = np.diag([5.0, 6.0, 7.0]).astype(np.float32)
        if offdiag:
            b[0, 1] = 0.1
        rep['box'] = b
    return atoms, X


def test_mapping_and_selected_reads(tmp_path):
    p = tmp_path / 'a.h5'
    atoms, X = write_fixture(p)
    d = MdcathDomain(p)
    m = d.mapping
    assert m['backbone_resids'] == [1, 2, 3] and m['backbone_contiguous']
    bb = d.backbone('320', '0')
    names = [(a[0], a[2]) for a in atoms]
    for r in (1, 2, 3):
        for j, a in enumerate(('N', 'CA', 'C', 'O')):
            k = names.index((r, a))
            assert np.array_equal(bb[:, r - 1, j], X[:, k].astype(np.float64))
    assert np.allclose(d.box_A('320', '0'), [50.0, 60.0, 70.0])
    assert bb.shape == (5, 3, 4, 3)
    bb2 = d.backbone('320', '0', chunk=2)
    assert np.array_equal(bb, bb2)


def test_missing_atom_and_invalid_metadata(tmp_path):
    p = tmp_path / 'b.h5'
    write_fixture(p, drop=(2, 'O'))
    d = MdcathDomain(p)
    assert d.mapping['missing']['O'] == [2]
    with pytest.raises(ReaderError):
        d.backbone('320', '0')
    assert d.backbone('320', '0', atoms=('N', 'CA', 'C')).shape == (5, 3, 3, 3)
    p2 = tmp_path / 'c.h5'
    write_fixture(p2, offdiag=True)
    with pytest.raises(ReaderError):
        MdcathDomain(p2).box_A('320', '0')
    p3 = tmp_path / 'd.h5'
    write_fixture(p3, unit='nm')
    with pytest.raises(ReaderError):
        MdcathDomain(p3).backbone('320', '0')


def test_make_whole_restores_split_chain_and_oxygen():
    box = np.array([30.0, 30.0, 30.0])
    L = 6
    t = np.arange(3 * L) * 1.4
    chain = np.stack([t, np.zeros_like(t), np.zeros_like(t)], 1).reshape(L, 3, 3)
    O = chain[:, 2] + np.array([0.0, 1.23, 0.0])
    bb = np.concatenate([chain, O[:, None]], 1)[None]
    wrapped = bb.copy()
    wrapped[0, 4:, :, 0] -= 30.0                 # residues 4–5 shifted by one box length
    wrapped[0, 2, 3, 1] += 30.0                  # one O shifted
    whole = make_whole_backbone(wrapped, box)
    assert np.allclose(whole, bb)
