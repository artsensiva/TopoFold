#!/usr/bin/env python3
"""C1 feasibility counts from the pinned Rascore label table (metadata only; no coordinates, no method run).
Reads sources/rascore_<commit>/src/rascore/util/data/entry.tsv; writes results/c1_counts.json.
Run: python3 c1_counts.py
"""
import collections
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMMIT = '8f77dc201cf9b1eb838d8e0d4d15bc0397dc3769'
ENTRY = HERE / 'sources' / f'rascore_{COMMIT}' / 'src' / 'rascore' / 'util' / 'data' / 'entry.tsv'
UNUSABLE = ('Outlier', 'Disordered')
CUTOFFS = ('2019-01-01', '2020-01-01', '2021-01-01', '2022-01-01')


def main():
    rows = list(csv.DictReader(open(ENTRY), delimiter='\t'))
    out = {'source': {'path': str(ENTRY.relative_to(HERE)), 'sha256': hashlib.sha256(ENTRY.read_bytes()).hexdigest(),
                      'commit': COMMIT, 'n_rows': len(rows)},
           'n_entries': len({r['pdb_code'] for r in rows}),
           'by_isoform': dict(collections.Counter(r['gene_class'] for r in rows)),
           'methods': dict(collections.Counter(r['method'] for r in rows)),
           'deposit_range': [min(r['deposit_data'] for r in rows), max(r['deposit_data'] for r in rows)]}
    for loop in ('SW1', 'SW2'):
        lab = collections.defaultdict(lambda: collections.Counter())
        groups = collections.defaultdict(lambda: {'entries': set(), 'crystal_forms': set(), 'ligands': set()})
        temporal = {c: collections.defaultdict(lambda: {'dev_entries': set(), 'test_entries': set()}) for c in CUTOFFS}
        for r in rows:
            c = r[loop]
            lab[c][r['gene_class']] += 1
            lab[c]['all_chains'] += 1
            g = groups[c]
            g['entries'].add(r['pdb_code'])
            g['crystal_forms'].add(r['crystal_form'])
            g['ligands'].add(r['pharmacological_ligand'])
            for cut in CUTOFFS:
                temporal[cut][c]['test_entries' if r['deposit_data'] >= cut else 'dev_entries'].add(r['pdb_code'])
        out[loop] = {c: {'chains': dict(lab[c]), 'distinct_entries': len(groups[c]['entries']),
                         'distinct_crystal_forms': len(groups[c]['crystal_forms']),
                         'distinct_pharmacological_ligands': len(groups[c]['ligands'] - {'None'}),
                         'usable_label': c not in UNUSABLE}
                     for c in sorted(lab)}
        out[loop + '_usable_chains'] = sum(v['chains']['all_chains'] for c, v in out[loop].items() if c not in UNUSABLE)
        out[loop + '_temporal_entries_by_cutoff'] = {
            cut: {c: {'dev': len(v['dev_entries']), 'test': len(v['test_entries'])} for c, v in sorted(temporal[cut].items())
                  if c not in UNUSABLE} for cut in CUTOFFS}
        iso = collections.defaultdict(lambda: collections.defaultdict(set))
        for r in rows:
            if r[loop] not in UNUSABLE:
                iso[r[loop]]['KRAS' if r['gene_class'] == 'KRAS' else 'HRAS+NRAS'].add(r['pdb_code'])
        out[loop + '_isoform_entries'] = {c: {k: len(v) for k, v in d.items()} for c, d in sorted(iso.items())}
    both = [r for r in rows if r['SW1'] not in UNUSABLE and r['SW2'] not in UNUSABLE]
    out['both_loops_usable_chains'] = len(both)
    out['both_loops_usable_entries'] = len({r['pdb_code'] for r in both})
    per_entry = collections.Counter(r['pdb_code'] for r in rows)
    out['chains_per_entry'] = {'max': max(per_entry.values()), 'entries_with_gt1_chain': sum(v > 1 for v in per_entry.values())}
    mixed = collections.defaultdict(set)
    for r in rows:
        if r['SW1'] not in UNUSABLE:
            mixed[r['pdb_code']].add(r['SW1'])
    out['entries_with_multiple_SW1_labels'] = sum(len(v) > 1 for v in mixed.values())
    (HERE / 'results' / 'c1_counts.json').write_text(json.dumps(out, indent=1) + '\n')
    for k in ('n_entries', 'by_isoform', 'methods', 'deposit_range', 'SW1_usable_chains', 'SW2_usable_chains',
              'both_loops_usable_chains', 'both_loops_usable_entries', 'chains_per_entry', 'entries_with_multiple_SW1_labels'):
        print(k, out[k])
    for loop in ('SW1', 'SW2'):
        print(loop)
        for c, v in out[loop].items():
            print(f"  {c:20s} chains {v['chains']}  entries {v['distinct_entries']}  CF {v['distinct_crystal_forms']}")
        print('  isoform entries', out[loop + '_isoform_entries'])
        print('  temporal 2021', out[loop + '_temporal_entries_by_cutoff']['2021-01-01'])
        print('  temporal 2020', out[loop + '_temporal_entries_by_cutoff']['2020-01-01'])


if __name__ == '__main__':
    main()
