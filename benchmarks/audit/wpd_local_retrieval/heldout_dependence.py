"""Within-domain temporal-dependence diagnostics for the held-out tier (Amendment 2 §4).

Each paired contrast series (per domain, temperature, window size and target: the window-averaged per-query
difference NR_X − NR_Y in query-frame order) is classified as exactly one of:
- 'structural_identity'  — the two methods are declared identical for this configuration (e.g. a method compared
                           with itself); the series must be exactly zero, otherwise an error is raised;
- 'constant_observed'    — not a declared identity, but the observed series is constant within tolerance;
- 'estimation_failure'   — too short, no finite values, estimator returned nothing, or τ_int too large relative
                           to the series length for the Sokal window estimate to be trusted;
- 'estimated'            — τ_int estimated (frames).
Diagnostics are reported; they do not change the domain-level bootstrap, whose unit (the domain) is unchanged.
Harness code, not production. SYNTHETIC tests in test_heldout_dependence.py.
"""
import importlib.util
from pathlib import Path

import numpy as np

CONST_TOL = 1e-12
MIN_N = 50
MAX_TAU_FRACTION = 0.015      # calibrated: results/SYNTHETIC_tau_int_calibration.json (severe underestimation from 0.017·n)


def _tau_int():
    spec = importlib.util.spec_from_file_location(
        'wpb_adk_lib', Path(__file__).resolve().parent.parent / 'adk_feasibility' / 'adk_lib.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.tau_int


TAU_INT = _tau_int()


def is_structural_identity(x, y, W):
    """Declared identities only. In the held-out method set the only one is a method compared with itself.
    (Legacy F vs F_lock at W = 5 is an identity too, but F is not in the held-out set.)"""
    return x == y or ({x, y} == {'F', 'F_lock'} and W == 5)


def classify_series(series, identity=False):
    s = np.asarray(series, float)
    fin = s[np.isfinite(s)]
    out = {'n': int(s.size), 'n_finite': int(fin.size)}
    if identity:
        if fin.size and np.max(np.abs(fin)) > CONST_TOL:
            raise ValueError('declared structural identity violated: series is not zero')
        out['class'] = 'structural_identity'
        return out
    if fin.size < MIN_N:
        out.update({'class': 'estimation_failure', 'reason': 'too_short_or_nonfinite'})
        return out
    scale = 1.0 + abs(float(np.mean(fin)))
    if np.ptp(fin) <= CONST_TOL * scale:
        out.update({'class': 'constant_observed', 'value': float(np.mean(fin))})
        return out
    tau = TAU_INT(fin)
    if tau is None or not np.isfinite(tau):
        out.update({'class': 'estimation_failure', 'reason': 'estimator_returned_none'})
        return out
    out['tau_int_frames'] = float(tau)
    if tau > MAX_TAU_FRACTION * fin.size:
        out.update({'class': 'estimation_failure', 'reason': 'tau_not_resolvable_from_series_length'})
        return out
    out['class'] = 'estimated'
    return out


def comparator_for_dependence(x, best, runner_up):
    """Amendment 2 §4: when X is itself the best method, the dependence diagnostic for X uses the runner-up,
    never X against itself."""
    return runner_up if x == best else best
