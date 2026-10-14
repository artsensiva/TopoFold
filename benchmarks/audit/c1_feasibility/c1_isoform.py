#!/usr/bin/env python3
"""C1 within-isoform counts (existing local metadata only; no download, no method run).
SW1 ON (Y32in.3P-ON) vs OFF (Y32out.2P-OFF), X-ray chains, grouped by Rascore crystal form (CF).
Writes results/c1_isoform.json. Run: python3 c1_isoform.py
"""
import collections
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENTRY = next((HERE / 'sources').glob('rascore_*/src/rascore/util/data/entry.tsv'))
ON, OFF = 'Y32in.3P-ON', 'Y32out.2P-OFF'


def main():
    rows = [r for r in csv.DictReader(open(ENTRY), delimiter='\t') if r['SW1'] in (ON, OFF) and r['method'] == 'XRAY']
    out = {'source_entry_tsv': str(ENTRY.relative_to(HERE)), 'by_isoform': {}}
    for iso in sorted({r['gene_class'] for r in rows}):
        rs = [r for r in rows if r['gene_class'] == iso]
        cf = collections.defaultdict(set)
        for r in rs:
            cf[r['crystal_form']].add(r['SW1'])
        d = {}
        for lab in (ON, OFF):
            sub = [r for r in rs if r['SW1'] == lab]
            cfs = {r['crystal_form'] for r in sub}
            d[lab] = {'chains': len(sub), 'entries': len({r['pdb_code'] for r in sub}), 'crystal_forms': len(cfs),
                      'crystal_forms_single_class': sum(1 for c in cfs if cf[c] == {lab}),
                      'largest_crystal_form_chains': max((sum(1 for r in sub if r['crystal_form'] == c) for c in cfs),
                                                         default=0),
                      'mutation_status': dict(collections.Counter(r['mutation_status'] for r in sub)),
                      'bound_ligand_classes': dict(collections.Counter(r['pharmacological_class'] for r in sub))}
        d['crystal_forms_with_both_classes'] = sum(1 for v in cf.values() if len(v) > 1)
        out['by_isoform'][iso] = d
    shared = collections.defaultdict(set)
    for r in rows:
        shared[r['crystal_form']].add(r['gene_class'])
    out['crystal_forms_spanning_isoforms'] = sum(1 for v in shared.values() if len(v) > 1)
    (HERE / 'results' / 'c1_isoform.json').write_text(json.dumps(out, indent=1) + '\n')
    for iso, d in out['by_isoform'].items():
        print(iso, 'mixed CF', d['crystal_forms_with_both_classes'])
        for lab in (ON, OFF):
            x = d[lab]
            print(f"  {lab:14s} chains {x['chains']:3d} entries {x['entries']:3d} CF {x['crystal_forms']:2d} "
                  f"(single-class {x['crystal_forms_single_class']}, largest {x['largest_crystal_form_chains']}) "
                  f"mut {x['mutation_status']} pharm {x['bound_ligand_classes']}")
    print('CF spanning isoforms', out['crystal_forms_spanning_isoforms'])


if __name__ == '__main__':
    main()
