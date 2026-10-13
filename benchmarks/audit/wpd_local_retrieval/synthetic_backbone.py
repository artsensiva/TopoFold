"""SYNTHETIC backbone generator and mdCATH-layout HDF5 fixture writer for runner tests (not proteins, not MD).

Backbone N, CA, C (and O) placed by natural extension from fixed bond lengths/angles and prescribed φ, ψ (ω = π).
Fixtures may be wrapped into a periodic box so that molecules cross the boundary.
"""
import h5py
import numpy as np

B_NCA, B_CAC, B_CN, B_CO = 1.458, 1.525, 1.329, 1.231
A_NCAC, A_CACN, A_CNCA, A_CACO = np.radians(111.2), np.radians(116.2), np.radians(121.7), np.radians(120.5)


def place(a, b, c, bond, angle, tors):
    """Place D so that |CD| = bond, angle(B, C, D) = angle, dihedral(A, B, C, D) = tors (IUPAC)."""
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d = np.array([-bond * np.cos(angle), bond * np.sin(angle) * np.cos(tors), bond * np.sin(angle) * np.sin(tors)])
    return c + d[0] * bc + d[1] * m + d[2] * n


def backbone(phi, psi):
    """phi, psi: (L,) radians → (L, 4, 3) N, CA, C, O. φ_0 and ψ_{L−1} are unused."""
    L = len(phi)
    X = np.zeros((L, 4, 3))
    X[0, 0] = [0.0, 0.0, 0.0]
    X[0, 1] = [B_NCA, 0.0, 0.0]
    X[0, 2] = X[0, 1] + B_CAC * np.array([-np.cos(A_NCAC), np.sin(A_NCAC), 0.0])
    for i in range(1, L):
        N = place(X[i - 1, 0], X[i - 1, 1], X[i - 1, 2], B_CN, A_CACN, psi[i - 1])
        CA = place(X[i - 1, 1], X[i - 1, 2], N, B_NCA, A_CNCA, np.pi)
        C = place(X[i - 1, 2], N, CA, B_CAC, A_NCAC, phi[i])
        X[i, :3] = N, CA, C
    for i in range(L):
        X[i, 3] = place(X[i, 0], X[i, 1], X[i, 2], B_CO, A_CACO, psi[i] + np.pi)
    return X


def trajectory(L, F, seed, base=None, noise=0.25, shift=(0.0, 0.0, 0.0)):
    r = np.random.default_rng(seed)
    if base is None:
        base = (np.full(L, -1.1) + r.normal(0, 0.3, L), np.full(L, -0.7) + r.normal(0, 0.3, L))
    out = np.empty((F, L, 4, 3))
    for f in range(F):
        out[f] = backbone(base[0] + r.normal(0, noise, L), base[1] + r.normal(0, noise, L)) + np.asarray(shift)
    return out, base


def write_domain(path, trajs, box_nm, missing_o=(), nonfinite=None, name='syn00A00', caps=True):
    """trajs: {(temp, rep): (F, L, 4, 3) Å}. Writes the mdCATH layout with CHARMM-style cap patch atoms."""
    L = next(iter(trajs.values())).shape[1]
    atoms = []
    for r in range(1, L + 1):
        if caps and r == 1:
            atoms += [(r, 'CAY'), (r, 'CY'), (r, 'OY')]
        atoms += [(r, 'N'), (r, 'CA'), (r, 'C')]
        if r not in missing_o:
            atoms.append((r, 'O'))
        atoms.append((r, 'CB'))
        if caps and r == L:
            atoms += [(r, 'NT'), (r, 'CAT')]
    lines = ['MODEL        1'] + [f"ATOM  {k + 1:5d} {n:<4s} ALA 0{r:4d}    {0:8.3f}{0:8.3f}{0:8.3f}"
                                   for k, (r, n) in enumerate(atoms)] + ['ENDMDL', 'END']
    with h5py.File(path, 'w') as f:
        f.attrs['layout'] = 'mdcath-only-protein-v1.0'
        g = f.create_group(name)
        g.attrs['numProteinAtoms'] = len(atoms)
        g.attrs['numResidues'] = L
        g.attrs['numChains'] = 1
        g['pdbProteinAtoms'] = '\n'.join(lines).encode()
        g['resid'] = np.array([a[0] for a in atoms])
        g['resname'] = np.array([b'ALA'] * len(atoms), dtype=object)
        g['chain'] = np.array([b'0'] * len(atoms), dtype=object)
        g['z'] = np.array([7 if n.startswith('N') else 8 if n.startswith('O') else 6 for _, n in atoms])
        for (t, rep), X in trajs.items():
            F = X.shape[0]
            C = np.zeros((F, len(atoms), 3), np.float32)
            for k, (r, n) in enumerate(atoms):
                j = {'N': 0, 'CA': 1, 'C': 2, 'O': 3}.get(n)
                if j is not None:
                    C[:, k] = X[:, r - 1, j]
                else:
                    C[:, k] = X[:, r - 1, 1] + np.array([0.5, 0.5, 0.5])
            box = np.asarray(box_nm, float) * 10.0
            C = (C - box * np.floor(C / box)).astype(np.float32)          # wrap every atom into the box
            if nonfinite is not None and (t, rep) == nonfinite[0]:
                C[nonfinite[1], 0, 0] = np.nan
            grp = g.require_group(t).create_group(rep)
            grp.attrs['numFrames'] = F
            d = grp.create_dataset('coords', data=C)
            d.attrs['unit'] = 'Angstrom'
            grp['box'] = np.diag(np.asarray(box_nm, np.float32))
    return atoms
