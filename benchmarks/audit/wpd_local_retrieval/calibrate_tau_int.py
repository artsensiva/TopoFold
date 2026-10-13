#!/usr/bin/env python3
"""SYNTHETIC calibration of the Sokal τ_int estimator (WP-B adk_lib.tau_int) on AR(1) series.

Purpose: find for which ratios n / τ̂ the estimate can be trusted, since a short series of a slowly decorrelating
process yields a small τ̂ (systematic underestimation). Writes results/SYNTHETIC_tau_int_calibration.json.
Run: build/venv/bin/python -B calibrate_tau_int.py
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from heldout_dependence import TAU_INT  # noqa: E402

OUT = HERE / 'results' / 'SYNTHETIC_tau_int_calibration.json'
PHIS = (0.0, 0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.999)
NS = (400, 750, 1000, 5000)
REPS = 40
SEED = 20261013


def ar1(phi, n, rng):
    x = np.empty(n)
    x[0] = rng.normal() / np.sqrt(max(1 - phi * phi, 1e-12))     # stationary start
    e = rng.normal(size=n)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    return x


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    for n in NS:
        for phi in PHIS:
            true = (1 + phi) / (1 - phi)
            est = np.array([TAU_INT(ar1(phi, n, rng)) for _ in range(REPS)], float)
            rows.append({'n': n, 'phi': phi, 'tau_true': true, 'est_median': float(np.median(est)),
                         'est_p90': float(np.quantile(est, 0.9)), 'ratio_median_to_true': float(np.median(est) / true),
                         'n_over_est_median': float(n / np.median(est))})
    # Largest τ̂/n among cases with severe underestimation (median ratio < 0.8): any τ̂ at or above this fraction
    # could come from such a case.
    severe = [r for r in rows if r['ratio_median_to_true'] < 0.8]
    frac_min_severe = min(r['est_median'] / r['n'] for r in severe) if severe else None
    out = {'_meta': {'label': 'SYNTHETIC AR(1) calibration', 'reps': REPS, 'seed': SEED,
                     'estimator': 'WP-B adk_lib.tau_int (Sokal window c = 5)'},
           'rows': rows,
           'min_est_fraction_among_severely_biased_cases': frac_min_severe}
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    for r in rows:
        print(r['n'], r['phi'], round(r['tau_true'], 1), round(r['est_median'], 1), round(r['ratio_median_to_true'], 2))
    print('min τ̂/n among severely biased (ratio < 0.8):', frac_min_severe)


if __name__ == '__main__':
    main()
