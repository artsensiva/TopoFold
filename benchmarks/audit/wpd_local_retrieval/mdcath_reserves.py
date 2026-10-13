#!/usr/bin/env python3
"""Complete the held-out cohort from the frozen reserve queue (authorized 2026-10-06: at most the 10 frozen reserves,
in order). For each reserve, in order: download from the pinned revision, verify recorded size and SHA-256
(mandatory), apply the unchanged structural eligibility criteria (heldout_runner.structural_eligibility). Stop as
soon as 30 eligible domains exist (24 eligible primaries + eligible reserves). Failing primaries are replaced by
eligible reserves in reserve order (original §4.3). Writes results/wpd_mdcath_reserve_ledger.json and
results/wpd_mdcath_cohort.json. Exit 2 on integrity failure; exit 4 if the queue is exhausted (insufficient evidence).
Run: build/venv/bin/python -B mdcath_reserves.py
"""
import datetime
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from heldout_runner import MandatoryCheckError, structural_eligibility  # noqa: E402

RES = HERE / 'results'
SEL = RES / 'wpd_mdcath_selection.json'
LEDGER = RES / 'wpd_mdcath_eligibility_ledger.json'
OUT_LEDGER = RES / 'wpd_mdcath_reserve_ledger.json'
COHORT = RES / 'wpd_mdcath_cohort.json'
DATA = HERE / 'data' / 'mdcath'
SELECTION_SHA256 = 'de0736ad3d832e4342f7b4522723ea2915d6ae39c43a9155bc22bc364477d16c'
TARGET_N = 30
MAX_RESERVES = 10


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
        print('STOP: selection changed')
        return 2
    sel = json.loads(SEL.read_text())
    led = json.loads(LEDGER.read_text())
    prim = {d['domain']: d for d in led['domains']}
    primaries = [f for f in sel['files'] if f['role'] == 'primary']
    failing = [f for f in primaries if prim[f['domain']].get('eligible') is False]
    eligible_primary = [f for f in primaries if prim[f['domain']].get('eligible') is True]
    if any(prim[f['domain']].get('eligible') is None for f in primaries):
        print('STOP: primary eligibility incomplete')
        return 2
    reserves = [f for f in sel['files'] if f['role'] == 'reserve'][:MAX_RESERVES]
    out = {'_meta': {'producer': 'mdcath_reserves.py', 'authorization': '2026-10-06 user authorization: frozen reserves '
                     'in order, at most 10, stop at 30 eligible', 'started_utc': now(),
                     'primary_ledger_sha256': sha256(LEDGER), 'selection_sha256': SELECTION_SHA256},
           'reserves': []}
    accepted = []
    for f in reserves:
        if len(eligible_primary) + len(accepted) >= TARGET_N:
            break
        path = DATA / f"mdcath_dataset_{f['domain']}.h5"
        rec = {'order': f['order'], 'domain': f['domain'], 'superfamily': f['superfamily'], 'url': f['url'],
               'expected_size': f['size_bytes'], 'expected_sha256': f['lfs_sha256'], 'download_started_utc': now()}
        for _ in range(6):
            if path.exists() and path.stat().st_size == f['size_bytes']:
                break
            p = subprocess.run(['curl', '-sS', '-L', '--fail', '-C', '-', '--retry', '5', '--retry-delay', '10',
                                '-o', str(path), f['url']], capture_output=True, text=True)
            rec.setdefault('curl_exit_codes', []).append(p.returncode)
        rec['download_finished_utc'] = now()
        try:
            e = structural_eligibility(path, f['size_bytes'], f['lfs_sha256'])
        except MandatoryCheckError as err:
            rec.update({'integrity_verified': False, 'error': str(err)})
            out['reserves'].append(rec)
            OUT_LEDGER.write_text(json.dumps(out, indent=1) + '\n')
            print('MANDATORY CHECK FAILED', f['domain'], err, flush=True)
            return 2
        rec.update({'integrity_verified': True, 'size': e['size'], 'sha256': e['sha256'], 'eligible': e['eligible'],
                    'reasons': e['reasons'], 'n_backbone_residues': e.get('n_backbone_residues'),
                    'numFrames': e.get('numFrames'), 'descriptive': e.get('descriptive')})
        if e['eligible']:
            rec['replaces_primary'] = failing[len(accepted)]['domain']
            accepted.append(f)
        out['reserves'].append(rec)
        OUT_LEDGER.write_text(json.dumps(out, indent=1) + '\n')
        print(f['order'], f['domain'], 'eligible' if e['eligible'] else f'INELIGIBLE {e["reasons"]}', flush=True)
    complete = len(eligible_primary) + len(accepted) >= TARGET_N
    out['_meta'].update({'finished_utc': now(), 'cohort_complete': complete,
                         'reserves_downloaded': len(out['reserves']), 'reserves_accepted': [f['domain'] for f in accepted]})
    OUT_LEDGER.write_text(json.dumps(out, indent=1) + '\n')
    rmap = {r['domain']: r for r in out['reserves']}
    cohort = []
    for f in eligible_primary + accepted:
        r = rmap.get(f['domain'])
        cohort.append({'order': f['order'], 'domain': f['domain'], 'role': f['role'], 'superfamily': f['superfamily'],
                       'size_bytes': f['size_bytes'], 'sha256': f['lfs_sha256'],
                       'replaces_primary': r.get('replaces_primary') if r else None})
    COHORT.write_text(json.dumps({'_meta': {'producer': 'mdcath_reserves.py', 'written_utc': now(),
                                            'complete': complete, 'n': len(cohort),
                                            'excluded_primaries': {f['domain']: prim[f['domain']]['reasons'] for f in failing},
                                            'primary_ledger_sha256': sha256(LEDGER), 'reserve_ledger_sha256': sha256(OUT_LEDGER)},
                                  'domains': cohort}, indent=1) + '\n')
    print('cohort', len(cohort), 'complete' if complete else 'INCOMPLETE (insufficient evidence rule)')
    return 0 if complete else 4


if __name__ == '__main__':
    sys.exit(main())
