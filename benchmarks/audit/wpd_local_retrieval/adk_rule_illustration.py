#!/usr/bin/env python3
"""Apply the proposed held-out criterion (docs/WPD_HELDOUT_MDCATH_PROTOCOL.md §7) to the AdK development
results as a point-estimate illustration. Reads the frozen AdK results; writes a separate file.

AdK has no domain sample, so windows stand in for domains and no interval is computed: this is not a
decision and not evidence about other proteins. rho was fixed before this script was run; the AdK medians
were already known when it was chosen (disclosed in the protocol).
Third-party MD derivative (Zenodo 10.5281/zenodo.5583119).
Run: build/venv/bin/python -B adk_rule_illustration.py
"""
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
IN = HERE / 'results' / 'wpd_adk_THIRD_PARTY_MD_DERIVATIVE_retrieval.json'
OUT = HERE / 'results' / 'wpd_adk_THIRD_PARTY_MD_DERIVATIVE_rule_illustration.json'
METHODS = ('M1', 'M2', 'B1', 'B2', 'B3', 'C1')
NON_C1 = ('M1', 'M2', 'B1', 'B2', 'B3')
RHO = 0.10


def main():
    r = json.loads(IN.read_text())
    out = {'_meta': {'producer': 'adk_rule_illustration.py', 'label': 'THIRD-PARTY MD DERIVATIVE; illustration only',
                     'input_sha256': hashlib.sha256(IN.read_bytes()).hexdigest(), 'rho': RHO,
                     'statistics': 'G = 1 - median_w NR@10; L = median_w(NR_m - NR_b) / G_b, b = argmin median NR '
                                   'among M1, M2, B1, B2, B3 (C1 excluded from best selection: conditional coverage)',
                     'caveat': 'windows are not independent units; no interval; C1 values are on its own coverage'},
           'cells': {}}
    for traj in ('open_start', 'closed_start'):
        for W in (5, 8):
            ws = [w for w in r['windows'] if w['trajectory'] == traj and w['W'] == W and 'methods' in w]
            NR = {m: np.array([w['methods'][m]['k10']['mean_norm_regret'] for w in ws], float) for m in METHODS}
            med = {m: float(np.nanmedian(NR[m])) for m in METHODS}
            b = min(NON_C1, key=lambda m: med[m])
            Gb = 1.0 - med[b]
            cell = {'best': b, 'G': {m: 1.0 - med[m] for m in METHODS}, 'L': {}, 'L_le_rho': {}}
            for m in METHODS:
                d = NR[m] - NR[b]
                L = float(np.nanmedian(d) / Gb)
                cell['L'][m] = L
                cell['L_le_rho'][m] = bool(L <= RHO)
            cell['G_gt_rho'] = {m: bool(cell['G'][m] > RHO) for m in METHODS}
            out['cells'][f'{traj}/W{W}'] = cell
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    for k, c in out['cells'].items():
        print(k, 'best', c['best'], {m: round(c['L'][m], 3) for m in METHODS}, {m: round(c['G'][m], 3) for m in METHODS})


if __name__ == '__main__':
    main()
