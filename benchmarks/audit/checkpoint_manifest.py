#!/usr/bin/env python3
"""Checkpoint provenance manifest (2026-10-13). Not a benchmark: no scientific computation is rerun.

It reads files and **writes** one file, `checkpoint_2026-10-13_manifest.json`. For every input or derived artifact
that is excluded from Git, it records the size and SHA-256 and compares them with the hash already recorded in the
tracked results files (where one exists). It also checks that the compressed held-out results file decompresses to
the protected SHA-256, and records tool versions.
Run from benchmarks/audit/: ../../benchmarks/audit/wpd_local_retrieval/build/venv/bin/python -B checkpoint_manifest.py
"""
import datetime
import hashlib
import json
import lzma
import platform
import subprocess
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
OUT = HERE / 'checkpoint_2026-10-13_manifest.json'
WPD = HERE / 'wpd_local_retrieval'
PROTECTED_HELDOUT = '8dbdcf8c813c4a49c84fa6458a50c3de21dc54acb08194d7c66f8e1e55d5f5ab'


def sha256(path, algo='sha256'):
    h = hashlib.new(algo)
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def rel(p):
    return str(Path(p).resolve().relative_to(REPO))


def entry(p, recorded=None, algo='sha256'):
    p = Path(p)
    if not p.exists():
        return {'path': rel(p) if p.is_absolute() else str(p), 'present': False, 'recorded': recorded, 'match': None}
    h = sha256(p, algo)
    return {'path': rel(p), 'present': True, 'size_bytes': p.stat().st_size, algo: h, 'recorded': recorded,
            'match': None if recorded is None else h == recorded}


def array_hash(p):
    return hashlib.sha256(np.load(p).tobytes()).hexdigest()


def tool(cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=60).stdout.strip().splitlines()[0]
    except Exception as e:                                        # recorded, not hidden
        return f'unavailable: {e!r}'


def main():
    out = {'recorded_utc': datetime.datetime.now(datetime.UTC).isoformat(), 'producer': 'checkpoint_manifest.py',
           'git_head': tool(['git', '-C', str(REPO), 'rev-parse', 'HEAD'])}
    # 1. compressed protected artifact
    raw = WPD / 'results' / 'wpd_heldout_THIRD_PARTY_MD_DERIVATIVE_results.json'
    xz = raw.with_suffix('.json.xz')
    h = hashlib.sha256()
    with lzma.open(xz) as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    out['compressed_in_git'] = [{
        'original': entry(raw, PROTECTED_HELDOUT), 'compressed': entry(xz),
        'decompressed_sha256': h.hexdigest(), 'decompressed_matches_protected': h.hexdigest() == PROTECTED_HELDOUT,
        'reason': 'original 131 MB exceeds the GitHub 100 MB per-file limit; xz -9 copy is lossless',
        'use': 'xz -dk wpd_heldout_THIRD_PARTY_MD_DERIVATIVE_results.json.xz before running make_report_heldout.py or '
               'make_closeout_addendum.py'}]
    # 2. third-party inputs (not redistributed in Git; CC BY 4.0 sources with pinned retrieval)
    spec_md5 = {'ORADD_open_state.gro': 'a73df3d6a60e69515bb7a7758af59f17',
                'ORADD_closed_state.gro': '541635b532a85a5e98ff14a7228c6a2b',
                'ORADD_open_state_traj.xtc': '05de6fa777dd5c86ba3e4c32b2c9080f',
                'ORADD_closed_state_traj.xtc': '014063067a9617dd89a39e35f7350d25'}   # from adk_feasibility/SPEC.md
    adk = []
    for name, md5 in spec_md5.items():
        p = HERE / 'adk_feasibility' / 'data' / 'raw' / name
        e = entry(p, md5, 'md5')
        if p.exists():
            e['sha256'] = sha256(p)
        adk.append(e)
    co = json.loads((WPD / 'results' / 'wpd_mdcath_cohort.json').read_text())
    mdc = [entry(WPD / 'data' / 'mdcath' / f"mdcath_dataset_{d['domain']}.h5", d['sha256']) for d in co['domains']]
    cohort_names = {f"mdcath_dataset_{d['domain']}.h5" for d in co['domains']}
    extra = sorted(p.name for p in (WPD / 'data' / 'mdcath').glob('*.h5') if p.name not in cohort_names)
    out['third_party_inputs_not_in_git'] = {
        'adk_orädd_2021': {'source': 'Zenodo 10.5281/zenodo.5583119, CC BY 4.0', 'retrieval': 'adk_feasibility/README.md',
                           'files': adk},
        'mdcath_heldout_cohort': {'source': 'mdCATH (doi:10.1038/s41597-024-04140-z), CC BY 4.0; pinned Hugging Face '
                                            'revision URLs in wpd_local_retrieval/results/wpd_mdcath_*ledger.json',
                                  'files': mdc, 'local_h5_not_in_cohort': extra}}
    # 3. derived caches and references (regenerable; hashes recorded by their producers)
    eb1m = json.loads((HERE / 'eb1' / 'results' / 'eb1_measurements.json').read_text())
    a1m = json.loads((HERE / 'eb1' / 'results' / 'eb1a1_measurements.json').read_text())
    der = [entry(HERE / 'eb1' / 'build' / 'cache' / n, h) for n, h in eb1m['prep_cache_sha256'].items()]
    der += [entry(HERE / 'eb1' / 'build' / 'cache_a1' / n, h) for n, h in a1m['cache_a1_sha256'].items()]
    refs = []
    b1 = json.loads((HERE / 'b2opt' / 'results' / 'b2opt_measurements.json').read_text())
    for wid, v in b1['refcheck']['workloads'].items():
        p = HERE / 'b2opt' / 'build' / 'ref' / f'{wid}.npy'
        refs.append({'path': rel(p), 'present': p.exists(), 'recorded_array_sha256': v['sha256'],
                     'match': (array_hash(p) == v['sha256']) if p.exists() else None})
    b2 = json.loads((HERE / 'bq2' / 'results' / 'bq2_measurements.json').read_text())
    for wid, v in b2['ref'].items():
        p = HERE / 'bq2' / 'build' / 'ref' / f'{wid}.npy'
        if v and v.get('status') == 'ok':
            refs.append({'path': rel(p), 'present': p.exists(), 'recorded_array_sha256': v['sha256'],
                         'match': (array_hash(p) == v['sha256']) if p.exists() else None})
    ck = sorted((WPD / 'build' / 'heldout_checkpoints').glob('*.json'))
    out['derived_not_in_git'] = {
        'eb1_stage_caches': {'regenerate': 'eb1_run.py / eb1a1_run.py prep steps from the AdK and mdCATH inputs',
                             'files': der},
        'reference_distance_arrays': {'regenerate': 'b2opt_run.py refcheck / bq2_run.py ref steps', 'files': refs},
        'heldout_checkpoints': {'note': 'per-domain checkpoints of the single authorized held-out run; they are not '
                                        'regenerable without a new authorized run',
                                'files': [entry(p) for p in ck]}}
    # 4. environment
    freeze = WPD / 'build' / 'pip_freeze_after_h5py.txt'
    out['environment'] = {'python': platform.python_version(), 'platform': platform.platform(),
                          'rustc': tool(['rustc', '--version']), 'cargo': tool(['cargo', '--version']),
                          'xz': tool(['xz', '--version']),
                          'pip_freeze_copy': 'benchmarks/audit/ENVIRONMENT_audit_venv_pip_freeze.txt',
                          'pip_freeze_source': entry(freeze)}
    checks = [e['match'] for e in adk + mdc + der] + [r['match'] for r in refs]
    out['summary'] = {'all_recorded_hashes_match': all(c is True for c in checks if c is not None),
                      'missing_local_files': sum(c is None for c in checks),
                      'compressed_heldout_ok': out['compressed_in_git'][0]['decompressed_matches_protected']}
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + '\n')
    print(json.dumps(out['summary']))


if __name__ == '__main__':
    main()
