#!/usr/bin/env python3
"""WP-D: AdK backbone eligibility for a full-backbone φ/ψ reference target (read-only).

A separate artifact, not a rewrite of WP-B G0. It reads the WP-B raw files after verifying their
SHA-256 against the WP-B G0 record, and checks:
- atom mapping: N, CA, C present exactly once per residue 1–214;
- the neighbour atoms that φ (C_{r−1}, N, CA, C) and ψ (N, CA, C, N_{r+1}) require;
- finite coordinates in every frame;
- backbone covalent continuity (N–CA, CA–C, C–N_{+1}) in raw coordinates, to detect periodic
  splitting or chain breaks.
No dihedral, target or representation is computed.
Third-party MD (Zenodo 10.5281/zenodo.5583119, CC BY 4.0): not SYNTHETIC.
Run: ../adk_feasibility/build/venv/bin/python -B adk_backbone_eligibility.py
"""
import datetime
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WPB = HERE.parent / 'adk_feasibility'
RAW = WPB / 'data' / 'raw'
OUT = HERE / 'results' / 'wpd_adk_backbone_eligibility.json'
PAIRS = {'open_start': ('ORADD_open_state.gro', 'ORADD_open_state_traj.xtc'),
         'closed_start': ('ORADD_closed_state.gro', 'ORADD_closed_state_traj.xtc')}
RESIDS = list(range(1, 215))
BOND_MAX_A = 2.0      # covalent backbone bonds are ~1.3–1.5 Å; > 2.0 Å means a split or break


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 22), b''):
            h.update(c)
    return h.hexdigest()


def main():
    import MDAnalysis as mda
    g0 = json.loads((WPB / 'results' / 'adk_g0_provenance.json').read_text())
    out = {'_meta': {'producer': 'adk_backbone_eligibility.py', 'data': 'third-party MD (not SYNTHETIC)',
                     'source_doi': '10.5281/zenodo.5583119', 'mdanalysis': mda.__version__,
                     'started_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                     'relation_to_wpb': 'separate WP-D artifact; WP-B G0/G1 files are read, not modified',
                     'bond_max_A': BOND_MAX_A},
           'wpb_g0_sha256': sha256(WPB / 'results' / 'adk_g0_provenance.json')}
    bind = {n: sha256(RAW / n) == g0['files'][n]['sha256'] for p in PAIRS.values() for n in p}
    out['raw_files_match_wpb_g0'] = bind
    if not (g0.get('g0_passed') and all(bind.values())):
        out['eligible'] = False
        out['stop'] = 'raw files do not match the WP-B G0 record'
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(out, indent=1) + '\n')
        print('STOP', out['stop'])
        return 1

    out['trajectories'] = {}
    all_ok = True
    for lab, (gro, xtc) in PAIRS.items():
        u = mda.Universe(str(RAW / gro), str(RAW / xtc))
        prot = u.select_atoms('protein')
        mapping, idx = {}, {}
        for name in ('N', 'CA', 'C'):
            sel = prot.select_atoms(f'name {name}')
            counts = {r: int(np.sum(sel.resids == r)) for r in RESIDS}
            mapping[name] = {'missing': [r for r, c in counts.items() if c == 0],
                             'duplicated': [r for r, c in counts.items() if c > 1],
                             'outside_1_214': sorted(set(sel.resids.tolist()) - set(RESIDS))}
            if not mapping[name]['missing'] and not mapping[name]['duplicated']:
                order = np.argsort(sel.resids, kind='stable')
                idx[name] = sel[order]
        mapping_ok = all(not m['missing'] and not m['duplicated'] and not m['outside_1_214']
                         for m in mapping.values())
        # φ needs C_{r-1}, N_r, CA_r, C_r (r = 2..214); ψ needs N_r, CA_r, C_r, N_{r+1} (r = 1..213)
        phi_ok = [r for r in RESIDS[1:]] if mapping_ok else []
        psi_ok = [r for r in RESIDS[:-1]] if mapping_ok else []
        info = {'n_protein_atoms': int(prot.n_atoms), 'mapping': mapping, 'mapping_ok': mapping_ok,
                'residues_with_defined_phi': [min(phi_ok), max(phi_ok)] if phi_ok else None,
                'residues_with_defined_psi': [min(psi_ok), max(psi_ok)] if psi_ok else None,
                'terminal_atom_names': {
                    'first_residue': sorted(set(prot.select_atoms('resid 1').names.tolist())),
                    'last_residue': sorted(set(prot.select_atoms('resid 214').names.tolist()))}}
        if mapping_ok:
            N, CA, C = idx['N'], idx['CA'], idx['C']
            nonfinite = 0
            max_b = {'N-CA': 0.0, 'CA-C': 0.0, 'C-N+1': 0.0}
            min_b = {'N-CA': 9e9, 'CA-C': 9e9, 'C-N+1': 9e9}
            frames_with_long_bond, long_bond_examples = 0, []
            for k, ts in enumerate(u.trajectory):
                n, ca, c = N.positions.astype(float), CA.positions.astype(float), C.positions.astype(float)
                if not (np.isfinite(n).all() and np.isfinite(ca).all() and np.isfinite(c).all()):
                    nonfinite += 1
                    continue
                b = {'N-CA': np.linalg.norm(ca - n, axis=1), 'CA-C': np.linalg.norm(c - ca, axis=1),
                     'C-N+1': np.linalg.norm(n[1:] - c[:-1], axis=1)}
                long_here = False
                for key, v in b.items():
                    max_b[key] = max(max_b[key], float(v.max()))
                    min_b[key] = min(min_b[key], float(v.min()))
                    if v.max() > BOND_MAX_A:
                        long_here = True
                        if len(long_bond_examples) < 10:
                            long_bond_examples.append({'frame': k, 'bond': key,
                                                       'residue_index': int(np.argmax(v)) + 1,
                                                       'length_A': float(v.max())})
                frames_with_long_bond += long_here
            info.update({'n_frames': len(u.trajectory), 'frames_with_nonfinite_backbone': nonfinite,
                         'backbone_bond_max_A': max_b, 'backbone_bond_min_A': min_b,
                         'frames_with_bond_over_max': frames_with_long_bond,
                         'long_bond_examples': long_bond_examples})
            ok = nonfinite == 0 and frames_with_long_bond == 0
        else:
            ok = False
        info['eligible'] = bool(mapping_ok and ok)
        all_ok &= info['eligible']
        out['trajectories'][lab] = info
        print(lab, 'eligible' if info['eligible'] else 'NOT eligible')
    out['eligible'] = bool(all_ok)
    out['_meta']['finished_utc'] = datetime.datetime.now(datetime.UTC).isoformat()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    return 0 if all_ok else 1


if __name__ == '__main__':
    sys.exit(main())
