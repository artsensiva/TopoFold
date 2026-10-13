#!/usr/bin/env python3
"""WP-D held-out tier: metadata-only eligibility, schema sample and seeded domain draw (protocol
docs/WPD_HELDOUT_MDCATH_PROTOCOL.md §3–§4, §8). Reads no trajectory and no dynamics summary.

Refuses to run unless the held-out freeze record matches the protocol and metadata hashes.
Outputs: results/wpd_mdcath_selection.json, HELDOUT_DOWNLOAD_REQUEST.md.
Run: build/venv/bin/python -B mdcath_select.py
"""
import datetime
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
META = HERE / 'build' / 'metadata'
RES = HERE / 'results'
PROTOCOL = REPO / 'docs' / 'WPD_HELDOUT_MDCATH_PROTOCOL.md'
FREEZE = RES / 'wpd_heldout_protocol_freeze.json'
OUT = RES / 'wpd_mdcath_selection.json'
REQ = HERE / 'HELDOUT_DOWNLOAD_REQUEST.md'
ADK = RES / 'wpd_adk_THIRD_PARTY_MD_DERIVATIVE_retrieval.json'
ADK_PROV = RES / 'wpd_adk_THIRD_PARTY_MD_DERIVATIVE_provenance.json'

REVISION = '5e3ed8aec62b689e01751db16275fdcdbc39e47f'
URL = 'https://huggingface.co/datasets/compsciencelab/mdCATH/resolve/{rev}/data/mdcath_dataset_{dom}.h5'
SEED = 20261007
N_PRIMARY, N_RESERVE = 30, 10
LEN_MIN, LEN_MAX = 50, 300
FIXED_EXCLUDED_SF = {'3.40.50.300': 'AdK (E. coli adenylate kinase, 4ake/1ake) development protein',
                     '3.40.190.100': 'superfamily of 1r9lA02, example domain in the mdCATH authors notebook'}
METADATA_FILES = ('cath-domain-list-v4_2_0.txt', 'cath-names-v4_2_0.txt', 'README-cath-list-file-format.txt',
                  'mdCATH_domains.txt', 'mdcath_tree_data.json', 'mdcath_tree_root.json', 'mdcath_api.json')
NOMINAL_FRAMES = 450
BURN = 0.10


def sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 22), b''):
            h.update(c)
    return h.hexdigest()


def read_cath(path):
    out = {}
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        p = line.split()
        out[p[0]] = {'superfamily': '.'.join(p[1:5]), 'length': int(p[10])}
    return out


def read_names(path):
    out = {}
    for line in open(path):
        if line.startswith('#') or not line.strip():
            continue
        node, _, rest = line.rstrip('\n').partition(' ')
        out[node] = rest.split(':', 1)[-1].strip() if ':' in rest else rest.strip()
    return out


def sf_key(sf):
    return tuple(int(x) for x in sf.split('.'))


def metadata_pool(domains, cath, tree):
    """§4.1 criteria 1–3. Returns (eligible dict, counts of failures)."""
    fail = {'not_in_tree': 0, 'not_in_cath': 0, 'length_out_of_range': 0}
    pool = {}
    for d in domains:
        if d not in tree or 'lfs' not in tree[d]:
            fail['not_in_tree'] += 1
            continue
        if d not in cath:
            fail['not_in_cath'] += 1
            continue
        if not LEN_MIN <= cath[d]['length'] <= LEN_MAX:
            fail['length_out_of_range'] += 1
            continue
        pool[d] = {'superfamily': cath[d]['superfamily'], 'length': cath[d]['length'],
                   'size_bytes': int(tree[d]['size']), 'lfs_sha256': tree[d]['lfs']['oid']}
    return pool, fail


def schema_sample(pool):
    """§3.1: smallest file in bytes, ties by name."""
    return min(pool, key=lambda d: (pool[d]['size_bytes'], d))


def draw(pool, excluded_sf, seed=SEED, n=N_PRIMARY + N_RESERVE):
    """§4.2 two-stage draw: superfamilies uniformly, then one domain uniformly within."""
    by_sf = {}
    for d in sorted(pool):
        if pool[d]['superfamily'] not in excluded_sf:
            by_sf.setdefault(pool[d]['superfamily'], []).append(d)
    sfs = sorted(by_sf, key=sf_key)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(sfs))
    chosen = []
    for j in perm:
        members = by_sf[sfs[j]]
        chosen.append(members[int(rng.integers(0, len(members)))])
        if len(chosen) == n:
            break
    return chosen, len(sfs), sum(len(v) for v in by_sf.values())


def adk_throughput():
    r = json.loads(ADK.read_text())
    pairs = sum(w['n_query_valid'] * w['n_db_valid'] for w in r['windows'] if 'methods' in w)
    return pairs, r['_meta']['wall_seconds'], json.loads(ADK_PROV.read_text())['versions']['workers']


def main():
    if not FREEZE.exists():
        print('STOP: no held-out freeze record')
        return 1
    fr = json.loads(FREEZE.read_text())
    if sha256(PROTOCOL) != fr['protocol_sha256']:
        print('STOP: protocol changed after freeze')
        return 1
    if sha256(Path(__file__)) != fr['selection_script_sha256']:
        print('STOP: selection script changed after freeze')
        return 1
    for n in METADATA_FILES:
        if sha256(META / n) != fr['metadata_sha256'][n]:
            print(f'STOP: metadata {n} changed after freeze')
            return 1
    cath = read_cath(META / 'cath-domain-list-v4_2_0.txt')
    names = read_names(META / 'cath-names-v4_2_0.txt')
    domains = [l.strip() for l in open(META / 'mdCATH_domains.txt') if l.strip()]
    tree = {x['path'].split('mdcath_dataset_')[1][:-3]: x for x in json.loads((META / 'mdcath_tree_data.json').read_text())
            if x['path'].startswith('data/mdcath_dataset_')}
    pool, fail = metadata_pool(domains, cath, tree)
    schema = schema_sample(pool)
    excluded = dict(FIXED_EXCLUDED_SF)
    excluded[pool[schema]['superfamily']] = f'superfamily of the schema-inspection sample {schema} (development data)'
    chosen, n_sf, n_dom = draw(pool, set(excluded))
    if len(chosen) < N_PRIMARY + N_RESERVE:
        print('STOP: too few superfamilies')
        return 1
    assert len({pool[d]['superfamily'] for d in chosen}) == len(chosen)
    assert not any(pool[d]['superfamily'] in excluded for d in chosen)

    def entry(d, role, order):
        sf = pool[d]['superfamily']
        return {'order': order, 'role': role, 'domain': d, 'superfamily': sf, 'superfamily_name': names.get(sf),
                'cath_length': pool[d]['length'], 'size_bytes': pool[d]['size_bytes'],
                'lfs_sha256': pool[d]['lfs_sha256'], 'url': URL.format(rev=REVISION, dom=d)}
    files = [entry(schema, 'schema_sample_development', 0)]
    files += [entry(d, 'primary', i + 1) for i, d in enumerate(chosen[:N_PRIMARY])]
    files += [entry(d, 'reserve', N_PRIMARY + i + 1) for i, d in enumerate(chosen[N_PRIMARY:])]
    gb = 1e9
    size = {role: sum(f['size_bytes'] for f in files if f['role'] == role)
            for role in ('schema_sample_development', 'primary', 'reserve')}
    pairs, wall, workers = adk_throughput()
    rate = pairs / wall
    est_pairs = 0
    for f in files:
        if f['role'] != 'primary':
            continue
        nf = NOMINAL_FRAMES - int(np.floor(BURN * NOMINAL_FRAMES))
        for W in (5, 8):
            est_pairs += 2 * (f['cath_length'] - W + 1) * nf * nf    # two temperatures
    out = {'_meta': {'producer': 'mdcath_select.py', 'generated_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                     'protocol': str(PROTOCOL.relative_to(REPO)), 'protocol_sha256': fr['protocol_sha256'],
                     'freeze_sha256': sha256(FREEZE), 'revision': REVISION, 'seed': SEED,
                     'reads': 'metadata only (CATH 4.2.0 lists, mdCATH domain list and file tree); no trajectory'},
           'pool': {'mdcath_domains': len(domains), 'metadata_failures': fail, 'metadata_eligible': len(pool),
                    'excluded_superfamilies': excluded,
                    'domains_in_excluded_superfamilies': sum(pool[d]['superfamily'] in excluded for d in pool),
                    'drawable_superfamilies': n_sf, 'drawable_domains': n_dom},
           'files': files,
           'budget': {'schema_sample_bytes': size['schema_sample_development'], 'primary_bytes': size['primary'],
                      'reserve_bytes_if_all_used': size['reserve'],
                      'requested_total_bytes_without_reserves': size['schema_sample_development'] + size['primary'],
                      'max_total_bytes_with_all_reserves': sum(size.values()),
                      'free_disk_bytes_at_selection': None,
                      'runtime_estimate': {'adk_pair_window_evaluations': pairs, 'adk_wall_seconds': wall,
                                           'adk_workers': workers, 'adk_rate_per_second': rate,
                                           'nominal_frames_per_replica': NOMINAL_FRAMES, 'burn_fraction': BURN,
                                           'heldout_pair_window_evaluations_primary': est_pairs,
                                           'heldout_seconds_estimate': est_pairs / rate,
                                           'note': 'estimate from AdK throughput; not a benchmark'}}}
    import shutil
    out['budget']['free_disk_bytes_at_selection'] = shutil.disk_usage(HERE).free
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    b = out['budget']
    L = ['# WP-D held-out tier: requested mdCATH download scope', '',
         f"Generated by `mdcath_select.py` from metadata only; protocol `{out['_meta']['protocol']}` "
         f"(SHA-256 `{fr['protocol_sha256'][:16]}…`), freeze `{FREEZE.name}`. Pinned revision `{REVISION}`; draw seed "
         f'{SEED}. **Nothing below has been downloaded.**', '',
         f"Pool: {len(domains)} mdCATH domains; {len(pool)} metadata-eligible (CATH length {LEN_MIN}–{LEN_MAX}); "
         f"{out['pool']['drawable_domains']} drawable in {n_sf} superfamilies after excluding "
         f"{len(excluded)} superfamilies ({out['pool']['domains_in_excluded_superfamilies']} domains).", '',
         '| Request | Files | Size (GB) |', '|---|---|---|',
         f"| 1. Schema sample (development data, download first) | 1 | {b['schema_sample_bytes'] / gb:.3f} |",
         f"| 2. Primary held-out sample | {N_PRIMARY} | {b['primary_bytes'] / gb:.3f} |",
         f"| Requested now (1 + 2) | {N_PRIMARY + 1} | {b['requested_total_bytes_without_reserves'] / gb:.3f} |",
         f"| 3. Reserves, only on a logged replacement (maximum) | {N_RESERVE} | {b['reserve_bytes_if_all_used'] / gb:.3f} |",
         f"| Upper bound (1 + 2 + 3) | {N_PRIMARY + N_RESERVE + 1} | {b['max_total_bytes_with_all_reserves'] / gb:.3f} |",
         '', f"Free disk at selection: {b['free_disk_bytes_at_selection'] / gb:.0f} GB. Storage: "
         '`benchmarks/audit/wpd_local_retrieval/data/mdcath/` (git-ignored). Whole files, verified against the LFS '
         'SHA-256 below.', '',
         f"Runtime estimate (primary sample, both temperatures, W = 5 and 8, {NOMINAL_FRAMES} nominal frames, "
         f"{int(BURN * 100)} % burn-in): {b['runtime_estimate']['heldout_pair_window_evaluations_primary']:.3g} "
         f"pair-window evaluations ≈ {b['runtime_estimate']['heldout_seconds_estimate'] / 60:.0f} min at the AdK rate "
         f"({workers} workers). An estimate, not a benchmark.", '',
         '| # | Role | Domain | CATH superfamily | Name | Length | Size (GB) | LFS SHA-256 |',
         '|---|---|---|---|---|---|---|---|']
    for f in files:
        L.append(f"| {f['order']} | {f['role']} | {f['domain']} | {f['superfamily']} | {f['superfamily_name'] or '—'} | "
                 f"{f['cath_length']} | {f['size_bytes'] / gb:.3f} | `{f['lfs_sha256']}` |")
    L += ['', 'URL pattern: `' + URL.format(rev=REVISION, dom='<domain>') + '`', '']
    REQ.write_text('\n'.join(L) + '\n')
    print('wrote', OUT.name, REQ.name)
    return 0


if __name__ == '__main__':
    sys.exit(main())
