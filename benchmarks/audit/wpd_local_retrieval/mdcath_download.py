#!/usr/bin/env python3
"""Download the 30 PRIMARY held-out mdCATH files (authorized 2026-10-05) from the pinned revision and verify each
file's recorded size and SHA-256 (selection: results/wpd_mdcath_selection.json). Reserves are never downloaded
here. Resumes partial files (curl -C -). Writes results/wpd_mdcath_download_manifest.json after every file.

Exit status: 0 if all 30 verified; 2 if any file failed verification or download (mandatory check).
Run: build/venv/bin/python -B mdcath_download.py
"""
import datetime
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEL = HERE / 'results' / 'wpd_mdcath_selection.json'
DATA = HERE / 'data' / 'mdcath'
OUT = HERE / 'results' / 'wpd_mdcath_download_manifest.json'
SELECTION_SHA256 = 'de0736ad3d832e4342f7b4522723ea2915d6ae39c43a9155bc22bc364477d16c'


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 22), b''):
            h.update(c)
    return h.hexdigest()


def now():
    return datetime.datetime.now(datetime.UTC).isoformat()


def main():
    if sha256(SEL) != SELECTION_SHA256:
        print('STOP: selection file differs from the frozen selection')
        return 2
    sel = json.loads(SEL.read_text())
    primary = [f for f in sel['files'] if f['role'] == 'primary']
    assert len(primary) == 30
    DATA.mkdir(parents=True, exist_ok=True)
    man = {'_meta': {'producer': 'mdcath_download.py', 'authorization': '2026-10-05 user authorization: 30 primary '
                     'domains from the pinned revision; no reserves', 'revision': sel['_meta']['revision'],
                     'selection_sha256': SELECTION_SHA256, 'started_utc': now()}, 'files': []}
    ok_all = True
    for f in primary:
        path = DATA / f"mdcath_dataset_{f['domain']}.h5"
        rec = {'order': f['order'], 'domain': f['domain'], 'url': f['url'], 'expected_size': f['size_bytes'],
               'expected_sha256': f['lfs_sha256'], 'path': str(path.relative_to(HERE)), 'started_utc': now()}
        attempts = 0
        while True:
            attempts += 1
            if path.exists() and path.stat().st_size == f['size_bytes']:
                break
            p = subprocess.run(['curl', '-sS', '-L', '--fail', '-C', '-', '--retry', '5', '--retry-delay', '10',
                                '-o', str(path), f['url']], capture_output=True, text=True)
            rec.setdefault('curl_exit_codes', []).append(p.returncode)
            if p.returncode not in (0, 33) and attempts >= 3:
                break
            if attempts >= 6:
                break
        rec['attempts'] = attempts
        rec['size'] = path.stat().st_size if path.exists() else None
        rec['sha256'] = sha256(path) if path.exists() else None
        rec['verified'] = rec['size'] == f['size_bytes'] and rec['sha256'] == f['lfs_sha256']
        rec['finished_utc'] = now()
        ok_all &= rec['verified']
        man['files'].append(rec)
        man['_meta']['n_verified'] = sum(r['verified'] for r in man['files'])
        OUT.write_text(json.dumps(man, indent=1) + '\n')
        print(f['order'], f['domain'], 'verified' if rec['verified'] else 'FAILED', flush=True)
    man['_meta']['finished_utc'] = now()
    man['_meta']['all_verified'] = bool(ok_all)
    man['_meta']['total_bytes'] = sum(r['size'] or 0 for r in man['files'])
    OUT.write_text(json.dumps(man, indent=1) + '\n')
    return 0 if ok_all else 2


if __name__ == '__main__':
    sys.exit(main())
