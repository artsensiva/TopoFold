#!/usr/bin/env python3
"""WP-B G0: data integrity, topology and mapping (SPEC.md + AMENDMENTS.md, Amendments 1–2).

Third-party MD data (Orädd et al., Zenodo 10.5281/zenodo.5583119, CC BY 4.0): not SYNTHETIC.
No structural observable is computed. Exits nonzero on any stop condition.
Run: build/venv/bin/python -B g0.py
"""
import datetime
import hashlib
import json
import struct
import sys
import urllib.request
from pathlib import Path

import numpy as np

from adk_lib import THREE, map_reference, parse_pdb_ca

HERE = Path(__file__).resolve().parent
RAW = HERE / 'data' / 'raw'
RES = HERE / 'results'
SPEC = json.loads((HERE / 'wpb_spec.json').read_text())
AMD = json.loads((HERE / 'wpb_spec_amendment2.json').read_text())
SPEC_FILES = ['SPEC.md', 'wpb_spec.json', 'AMENDMENTS.md', 'wpb_spec_amendment1.json',
              'wpb_spec_amendment2.json']
OUT = RES / 'adk_g0_provenance.json'
RESIDS = list(range(AMD['resid_range'][0], AMD['resid_range'][1] + 1))


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    return h.hexdigest()


def first_frame_header(path):
    with open(path, 'rb') as fh:
        b = fh.read(200)
    magic, natoms, step = struct.unpack('>iii', b[:12])
    return {'magic': magic, 'natoms': natoms, 'step': step,
            'time_ps': struct.unpack('>f', b[12:16])[0],
            'precision_per_nm': struct.unpack('>f', b[56:60])[0]}


def main():
    import MDAnalysis as mda
    stops, notes = [], []
    prov = {'_meta': {'producer': 'g0.py', 'data': 'third-party MD (not SYNTHETIC)',
                      'source_doi': SPEC['source_doi'], 'licence': 'CC BY 4.0 (Zenodo record)',
                      'governing_spec': 'SPEC.md + AMENDMENTS.md (Amendments 1–2)',
                      'started_utc': datetime.datetime.now(datetime.UTC).isoformat()},
            'spec_sha256': {f: sha256(HERE / f) for f in SPEC_FILES},
            'environment': {'python': sys.version.split()[0], 'numpy': np.__version__,
                            'MDAnalysis': mda.__version__,
                            'pip_freeze': (HERE / 'build' / 'pip_freeze.txt').read_text().split()}}
    freeze = json.loads((RES / 'adk_spec_amendment2_freeze.json').read_text())
    prov['amendment2_freeze_utc'] = freeze['frozen_utc']
    if prov['spec_sha256'] != freeze['files']:
        stops.append('specification files differ from the Amendment-2 freeze record')

    try:
        with urllib.request.urlopen('https://zenodo.org/api/records/5583119', timeout=60) as r:
            api = {f['key']: f for f in json.load(r)['files']}
        prov['zenodo_api_at_g0'] = {k: {'size': v['size'], 'checksum': v['checksum']} for k, v in api.items()}
    except Exception as e:
        prov['zenodo_api_at_g0'] = {'error': repr(e)}

    files = {}
    for name, exp in SPEC['files'].items():
        p = RAW / name
        if not p.exists():
            stops.append(f'missing file {name}')
            continue
        m, s, size = md5(p), sha256(p), p.stat().st_size
        files[name] = {'size': size, 'md5': m, 'sha256': s,
                       'md5_matches_source_record': m == exp['md5'], 'size_matches_source_record': size == exp['size']}
        if m != exp['md5'] or size != exp['size']:
            stops.append(f'checksum/size mismatch {name}')
    prov['files'] = files

    refs = {}
    for rid in ('4AKE', '1AKE'):
        p = RAW / f'{rid}.pdb'
        if not p.exists():
            stops.append(f'missing reference {rid}')
            continue
        chains, het, revdat = parse_pdb_ca(p)
        entry = {'sha256': sha256(p), 'size': p.stat().st_size, 'revdat_dates': revdat,
                 'hetero_residues': het, 'chains': {}}
        for ch, role in ((AMD['primary_chains'][rid], 'primary'), (AMD['sensitivity_chains'][rid], 'sensitivity')):
            _, rep = map_reference(chains.get(ch, {}), AMD['declared_sequence'], RESIDS)
            entry['chains'][ch] = dict(rep, role=role)
            if not rep['valid']:
                (stops if role == 'primary' else notes).append(
                    f'{rid}:{ch} ({role}) mapping invalid: missing {rep["missing"]}, '
                    f'mismatched {rep["mismatched"]}')
        refs[rid] = entry
    prov['references'] = refs
    if stops:
        return finish(prov, stops, notes)

    trajs = {}
    for label, (gro, xtc) in SPEC['pairs'].items():
        u = mda.Universe(str(RAW / gro), str(RAW / xtc))
        ca = u.select_atoms(SPEC['selection'])
        seq = ''.join(THREE.get(r, 'X') for r in ca.resnames)
        times, dims = [], []
        for ts in u.trajectory:                       # timestamps and boxes only
            times.append(ts.time)
            dims.append(ts.dimensions.copy())
        times, dims = np.array(times, float), np.array(dims, float)
        dt = np.diff(times)
        vals, counts = np.unique(np.round(dt, 3), return_counts=True)
        nominal = float(vals[np.argmax(counts)])
        irregular = int(np.sum(np.abs(dt - nominal) > AMD['cadence_tolerance_ps']))
        hdr = first_frame_header(RAW / xtc)
        info = {'gro_n_atoms': int(mda.Universe(str(RAW / gro)).atoms.n_atoms),
                'xtc_n_atoms': int(u.trajectory.n_atoms), 'n_ca': int(ca.n_atoms),
                'ca_resids_equal_declared': bool(np.array_equal(ca.resids, np.array(RESIDS))),
                'n_segments': int(len(ca.segments)),
                'sequence_equals_declared': seq == AMD['declared_sequence'],
                'n_frames': int(len(times)), 'time_first_ps': float(times[0]),
                'time_last_ps': float(times[-1]),
                'dt_values_ps_with_counts': {str(v): int(c) for v, c in zip(vals, counts)},
                'nominal_dt_ps': nominal, 'cadence_irregularities': irregular,
                'times_strictly_increasing': bool(np.all(dt > 0)),
                'box_angles_max_dev_deg': float(np.abs(dims[:, 3:] - 90.0).max()),
                'box_lengths_A_min_max': [float(dims[:, :3].min()), float(dims[:, :3].max())],
                'raw_first_frame_header': hdr,
                'coordinate_quantization_A': 10.0 / hdr['precision_per_nm'],
                'units': 'MDAnalysis: Å and ps (xtc stores nm; converted on read)'}
        stated = SPEC['source_statements'][f'{label.split("_")[0]}_start_ns']
        info['source_statement_ns'] = stated
        info['file_duration_ns'] = (times[-1] - times[0]) / 1000.0
        info['frames_expected_from_statement'] = int(round(stated * 1000 / nominal)) + 1
        trajs[label] = info
        if info['gro_n_atoms'] != info['xtc_n_atoms']:
            stops.append(f'{label}: atom count gro != xtc')
        if info['n_ca'] != len(RESIDS) or not info['ca_resids_equal_declared'] or info['n_segments'] != 1:
            stops.append(f'{label}: Cα mapping failure')
        if not info['sequence_equals_declared']:
            stops.append(f'{label}: sequence differs from declared')
        if not info['times_strictly_increasing']:
            stops.append(f'{label}: non-monotonic timestamps')
        if info['box_angles_max_dev_deg'] > SPEC['pbc']['orthorhombic_tol_deg']:
            stops.append(f'{label}: non-orthorhombic box')
        if abs(info['n_frames'] - info['frames_expected_from_statement']) > 1:
            notes.append(f'{label}: {info["n_frames"]} frames vs {info["frames_expected_from_statement"]} '
                         f'expected from the stated {stated} ns (recorded discrepancy, not a stop)')
        if irregular:
            notes.append(f'{label}: {irregular} cadence irregularities (runs will be split there)')
    prov['trajectories'] = trajs
    return finish(prov, stops, notes)


def finish(prov, stops, notes):
    prov['stops'], prov['discrepancies_and_notes'] = stops, notes
    prov['g0_passed'] = not stops
    prov['_meta']['finished_utc'] = datetime.datetime.now(datetime.UTC).isoformat()
    OUT.write_text(json.dumps(prov, indent=1) + '\n')
    print('G0', 'PASS' if not stops else 'STOP', stops, notes)
    return 0 if not stops else 1


if __name__ == '__main__':
    sys.exit(main())
