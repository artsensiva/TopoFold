#!/usr/bin/env python3
"""Verify the minimal reader on the mdCATH schema-inspection sample (development data; original protocol §3.1).

Checks mapping, units, frame ordering, periodic-boundary conventions, missing data, bounded memory and binding
compatibility. No retrieval, target or score is computed. Writes results/wpd_mdcath_schema_check.json.
Third-party MD derivative (mdCATH, CC BY 4.0). Run: build/venv/bin/python -B mdcath_schema_check.py
"""
import datetime
import hashlib
import json
import resource
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
SEL = HERE / 'results' / 'wpd_mdcath_selection.json'
OUT = HERE / 'results' / 'wpd_mdcath_schema_check.json'
WHEEL_SHA = 'e4360f15875a532bc7b98196c7592ed4fc92672a57c0a621355961cafb17a6dd'
TRAJ = [('320', '0'), ('320', '1'), ('348', '0'), ('348', '1')]
BOND_MAX = 2.0


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 22), b''):
            h.update(c)
    return h.hexdigest()


def memory_probe(path, mode):
    """Run in a subprocess: 'baseline' imports only; 'backbone' reads the four trajectories' backbone atoms in
    chunks; 'full' reads one full coords dataset. Prints peak RSS (kB)."""
    import tracemalloc
    import h5py  # noqa: F401
    from mdcath_reader import MdcathDomain
    d = MdcathDomain(path)
    n = 0
    tracemalloc.start()
    if mode == 'backbone':
        for t, r in TRAJ:
            n += d.backbone(t, r).nbytes                  # each result released before the next read
    elif mode == 'full':
        for t, r in TRAJ:
            n += d.g[t][r]['coords'][:].nbytes
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    print(json.dumps({'mode': mode, 'bytes_returned_total': n, 'tracemalloc_peak_bytes': peak,
                      'maxrss_kB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}))


def kabsch_rmsd(a, b):
    a, b = a - a.mean(0), b - b.mean(0)
    U, _, Vt = np.linalg.svd(a.T @ b)
    s = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, s]) @ U.T
    return float(np.sqrt(np.mean(np.sum((a @ R.T - b) ** 2, axis=1))))


def main():
    if len(sys.argv) == 4 and sys.argv[1] == 'probe':
        return memory_probe(sys.argv[2], sys.argv[3])
    import h5py
    import topofold as tf
    import topofold.topofold as ext
    from mdcath_reader import MdcathDomain, make_whole_backbone
    from wpd_lib import dihedral, sin_angle, turning_angle, wrap
    sel = json.loads(SEL.read_text())
    entry = sel['files'][0]
    assert entry['role'] == 'schema_sample_development'
    path = HERE / 'data' / 'mdcath' / f"mdcath_dataset_{entry['domain']}.h5"
    out = {'_meta': {'producer': 'mdcath_schema_check.py', 'label': 'THIRD-PARTY MD DERIVATIVE (mdCATH); '
                     'development data (schema sample)', 'checked_utc': datetime.datetime.now(datetime.UTC).isoformat()},
           'source': {'domain': entry['domain'], 'url': entry['url'], 'revision': sel['_meta']['revision'],
                      'expected_size': entry['size_bytes'], 'size': path.stat().st_size,
                      'expected_lfs_sha256': entry['lfs_sha256'], 'sha256': sha256(path)},
           'dependencies': {'h5py': h5py.__version__, 'hdf5': h5py.version.hdf5_version, 'numpy': np.__version__,
                            'python': sys.version.split()[0], 'h5py_wheel': 'h5py-3.16.0-cp314-cp314-manylinux_2_28_x86_64.whl',
                            'h5py_wheel_sha256': WHEEL_SHA, 'topofold_extension_sha256': sha256(ext.__file__)}}
    src = out['source']
    out['source']['verified'] = src['size'] == src['expected_size'] and src['sha256'] == src['expected_lfs_sha256']
    d = MdcathDomain(path)
    g = d.g
    m = d.mapping
    L = len(m['backbone_resids'])
    first, last = m['resids'][0], m['resids'][-1]
    out['mapping'] = {
        'layout': 'mdcath-only-protein-v1.0', 'domain_attrs': {k: int(v) for k, v in g.attrs.items()},
        'n_atoms': len(d.atom_names), 'resid_range': [int(first), int(last)], 'n_residues': len(m['resids']),
        'backbone_residues': L, 'backbone_contiguous': m['backbone_contiguous'],
        'missing': m['missing'], 'duplicated': m['duplicated'],
        'first_residue_atoms': [d.atom_names[k] for k in np.where(d.resid == first)[0]],
        'last_residue_atoms': [d.atom_names[k] for k in np.where(d.resid == last)[0]],
        'separate_cap_residues': sorted({rn for rn in d.resname if rn in ('ACE', 'NME', 'NMA', 'CT3')}),
        'chains': sorted({x.decode() if isinstance(x, bytes) else str(x) for x in g['chain'][:]}),
        'hydrogens_present': bool(np.any(g['z'][:] == 1)),
        'cath_length_from_selection': entry['cath_length']}
    reps = d.replicas()
    out['replicas'] = {'numFrames': reps, 'all_temps_reps_present': all(len(v) == 5 for v in reps.values())
                       and len(reps) == 5}
    traj = {}
    for t, r in TRAJ:
        rep = g[t][r]
        ds = rep['coords']
        bb = d.backbone(t, r)
        nF = bb.shape[0]
        box = d.box_A(t, r)
        raw_box = rep['box'][:]
        N, CA, C, O = (bb[:, :, k] for k in range(4))
        def bonds(N, CA, C, O):
            return {'N-CA': np.linalg.norm(CA - N, axis=-1), 'CA-C': np.linalg.norm(C - CA, axis=-1),
                    'C-N+1': np.linalg.norm(N[:, 1:] - C[:, :-1], axis=-1), 'C-O': np.linalg.norm(O - C, axis=-1),
                    'CA-CA+1': np.linalg.norm(CA[:, 1:] - CA[:, :-1], axis=-1)}
        braw = bonds(N, CA, C, O)
        w = make_whole_backbone(bb, box)
        Nw, CAw, Cw, Ow = (w[:, :, k] for k in range(4))
        bw = bonds(Nw, CAw, Cw, Ow)
        long_raw = int(np.sum(np.any(np.concatenate([braw[k].reshape(nF, -1) for k in ('N-CA', 'CA-C', 'C-N+1',
                                                                                         'C-O')], 1) > BOND_MAX, 1)))
        long_whole = int(np.sum(np.any(np.concatenate([bw[k].reshape(nF, -1) for k in ('N-CA', 'CA-C', 'C-N+1',
                                                                                        'C-O')], 1) > BOND_MAX, 1)))
        extent = float(np.max(np.ptp(CAw, axis=1)))
        # frame-order plausibility: consecutive-frame Cα RMSD vs random pairs (no timestamps exist)
        rng = np.random.default_rng(20261014)
        cons = [kabsch_rmsd(CAw[i], CAw[i + 1]) for i in range(nF - 1)]
        ii, jj = rng.integers(0, nF, 400), rng.integers(0, nF, 400)
        rand = [kabsch_rmsd(CAw[i], CAw[j]) for i, j in zip(ii, jj) if abs(i - j) > 1]
        # binding compatibility and degeneracy on Cα
        CA32 = CAw.astype(np.float32)
        rt = rk = 0.0
        binfail = 0
        tref = dihedral(CAw[:, :-3], CAw[:, 1:-2], CAw[:, 2:-1], CAw[:, 3:])
        kref = turning_angle(CAw[:, :-2], CAw[:, 1:-1], CAw[:, 2:])
        for fidx in range(nF):
            try:
                k_, t_, _ = tf.compute_invariants(CA32[fidx])
            except Exception:  # noqa: BLE001
                binfail += 1
                continue
            rt = max(rt, float(np.abs(wrap(np.asarray(t_) - tref[fidx])).max()))
            rk = max(rk, float(np.abs(np.asarray(k_) - kref[fidx]).max()))
        sCA = sin_angle(CAw[:, :-2], CAw[:, 1:-1], CAw[:, 2:])
        traj[f'{t}/{r}'] = {
            'numFrames': int(rep.attrs['numFrames']), 'coords_shape': list(ds.shape), 'coords_dtype': str(ds.dtype),
            'coords_unit_attr': ds.attrs.get('unit'), 'coords_chunks': ds.chunks, 'coords_compression': ds.compression,
            'box_raw_diag': [float(x) for x in np.diag(raw_box)], 'box_offdiag_max': float(np.max(np.abs(
                raw_box - np.diag(np.diag(raw_box))))), 'box_A': [float(x) for x in box],
            'ca_extent_max_A': extent,
            'nonfinite_backbone_values': int(np.sum(~np.isfinite(bb))),
            'frames_with_backbone_bond_over_2A_raw': long_raw, 'frames_with_backbone_bond_over_2A_after_whole': long_whole,
            'frames_changed_by_whole': int(np.sum(np.any(np.abs(w - bb) > 0, axis=(1, 2, 3)))),
            'bond_lengths_after_whole_A': {k: [float(v.min()), float(np.median(v)), float(v.max())] for k, v in bw.items()},
            'consecutive_frame_ca_rmsd_median_A': float(np.median(cons)),
            'random_pair_ca_rmsd_median_A': float(np.median(rand)),
            'binding_failures': binfail, 'binding_vs_numpy_max_tau_diff': rt, 'binding_vs_numpy_max_kappa_diff': rk,
            'min_sin_ca_angle': float(sCA.min()), 'post_burnin_frames': nF - int(np.floor(0.10 * nF)),
            'meets_numFrames_ge_200': nF >= 200,
            'cis_omega_frame_bonds': int(np.sum(np.abs(dihedral(CAw[:, :-1], Cw[:, :-1], Nw[:, 1:], CAw[:, 1:]))
                                                < np.pi / 2))}
    out['trajectories_used_by_protocol'] = traj
    out['windows'] = {f'W{W}': L - W + 1 for W in (5, 8)}
    # bounded memory: subprocess measurements
    mem = {}
    for mode in ('baseline', 'backbone', 'full'):
        p = subprocess.run([sys.executable, '-B', __file__, 'probe', str(path), mode], capture_output=True, text=True)
        mem[mode] = json.loads(p.stdout.strip().splitlines()[-1])
    ds = g['320']['0']['coords']
    out['memory'] = {**mem, 'note': 'tracemalloc peak covers Python/NumPy allocations during the reads (NumPy '
                     'buffers are traced); maxrss is dominated by imports for a file this small',
                     'per_trajectory_full_coords_bytes': int(np.prod(ds.shape) * ds.dtype.itemsize),
                     'per_trajectory_backbone_float64_bytes': int(ds.shape[0] * L * 4 * 3 * 8),
                     'largest_primary_domain_backbone_float64_bytes_estimate': int(
                         450 * max(f['cath_length'] for f in sel['files'] if f['role'] == 'primary') * 4 * 3 * 8)}
    out['absent_metadata'] = ['per-frame timestamps (none stored; 1 ns spacing is from the paper, not the file)',
                              'temperature or replica attributes beyond group names',
                              'explicit box unit attribute (nm inferred from the authors notebook ×10 conversion; '
                              'checked below against protein extent)']
    out['observed_anomalies'] = []
    if _s(g['320']['0']['dssp'].attrs.get('unit')) == 'nm':
        out['observed_anomalies'].append("dssp dataset carries unit attribute 'nm' (DSSP codes have no unit); "
                                         'file-metadata error, not used')
    out['convention_differences_vs_protocol'] = []
    if not out['mapping']['separate_cap_residues']:
        out['convention_differences_vs_protocol'].append(
            'Caps are CHARMM patch atoms inside the terminal residues (first residue: CAY/CY/OY; last: NT/CAT), not '
            'separate residues. Original protocol §4.3 says "caps lack CA and are therefore not window residues"; '
            'the rule "exactly one N, CA, C" still selects every residue, so the window set is unaffected, but the '
            'wording is inaccurate.')
    out['verdict'] = {
        'source_verified': out['source']['verified'],
        'mapping_ok': (m['backbone_contiguous'] and L == out['mapping']['n_residues']
                       and not any(m['missing'][a] for a in ('N', 'CA', 'C', 'O'))),
        'units_ok': all(v['coords_unit_attr'] == 'Angstrom' and 1.2 < v['bond_lengths_after_whole_A']['N-CA'][1] < 1.6
                        and all(b > v['ca_extent_max_A'] for b in v['box_A']) for v in traj.values()),
        'pbc_ok_after_whole': all(v['frames_with_backbone_bond_over_2A_after_whole'] == 0 for v in traj.values()),
        'no_missing_values': all(v['nonfinite_backbone_values'] == 0 for v in traj.values()),
        'binding_reference_ok': all(v['binding_failures'] == 0 and v['binding_vs_numpy_max_tau_diff'] <= 1e-9
                                    and v['binding_vs_numpy_max_kappa_diff'] <= 1e-9 for v in traj.values()),
        'frame_order_plausible': all(v['consecutive_frame_ca_rmsd_median_A'] < v['random_pair_ca_rmsd_median_A']
                                     for v in traj.values())}
    OUT.write_text(json.dumps(out, indent=1, default=str) + '\n')
    print(json.dumps(out['verdict'], indent=1))
    return 0


def _s(x):
    return x.decode() if isinstance(x, bytes) else (None if x is None else str(x))


if __name__ == '__main__':
    sys.exit(main())
