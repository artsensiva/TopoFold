#!/usr/bin/env python3
"""C1 grouping and precision inputs for the proposed SW1 ON-vs-OFF task (metadata only; no method run).
Unit of independence: Rascore crystal form (CF) for X-ray chains; entries never split. Writes results/c1_groups.json.
Run: python3 c1_groups.py
"""
import collections
import csv
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ENTRY = next((HERE / 'sources').glob('rascore_*/src/rascore/util/data/entry.tsv'))
CLASSES = ('Y32in.3P-ON', 'Y32out.2P-OFF')
DISCORDANCE = (0.10, 0.20, 0.30)          # assumed paired disagreement rates between two methods (stated, not measured)


def main():
    rows = [r for r in csv.DictReader(open(ENTRY), delimiter='\t') if r['SW1'] in CLASSES]
    xr = [r for r in rows if r['method'] == 'XRAY']
    cf = collections.defaultdict(lambda: collections.Counter())
    for r in xr:
        cf[r['crystal_form']][r['SW1']] += 1
    mixed = [k for k, v in cf.items() if len(v) > 1]
    pure = {c: [k for k, v in cf.items() if len(v) == 1 and c in v] for c in CLASSES}
    out = {'task': 'SW1 binary: Y32in.3P-ON vs Y32out.2P-OFF; X-ray chains; group = crystal form',
           'chains_all_methods': dict(collections.Counter(r['SW1'] for r in rows)),
           'chains_xray': dict(collections.Counter(r['SW1'] for r in xr)),
           'entries_xray': {c: len({r['pdb_code'] for r in xr if r['SW1'] == c}) for c in CLASSES},
           'crystal_forms_xray_total': len(cf), 'crystal_forms_mixed_classes': len(mixed),
           'crystal_forms_single_class': {c: len(v) for c, v in pure.items()},
           'largest_crystal_form_chains': max(sum(v.values()) for v in cf.values()),
           'isoforms_xray': {c: dict(collections.Counter(r['gene_class'] for r in xr if r['SW1'] == c)) for c in CLASSES}}
    n = len(cf)                                # leave-one-crystal-form-out: every CF is a held-out unit once
    out['precision_95ci_halfwidth_for_paired_difference'] = {
        'unit': 'crystal form (conservative; chains within a CF are not independent)', 'n_units': n,
        'formula': '1.96 * sqrt(p_d / n) (paired difference near 0, McNemar-type, units = CFs)',
        'by_assumed_discordance': {str(p): round(1.96 * math.sqrt(p / n), 3) for p in DISCORDANCE}}
    (HERE / 'results' / 'c1_groups.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
