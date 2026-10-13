"""Dependence diagnostics, revision 1 (Amendment 3 §4). Supersedes heldout_dependence.py (preserved unmodified).

Corrections:
- The τ̂ ≤ 0.015·n cutoff is a LIMITED HEURISTIC derived from an AR(1) simulation of this estimator
  (results/SYNTHETIC_tau_int_calibration.json). It is not a reliability guarantee: processes with drift or several
  timescales can give a small τ̂ while true dependence is long, and AR(1) underestimation already starts below the
  cutoff for short series. 'estimated' therefore means "τ̂ is below the heuristic cutoff", nothing more.
- Non-finite observations are never dropped by compressing time. The estimator runs only on the longest contiguous
  finite segment, and gaps are reported; with too few contiguous values the class is 'insufficient_observations'.
- A declared identity with unavailable (non-finite) observations is distinguished from an identity violation
  (a finite non-zero value).
Classes: structural_identity | identity_violation (error) | constant_observed | insufficient_observations |
         estimation_failure | above_heuristic_cutoff | estimated_below_heuristic_cutoff.
"""
import numpy as np

from heldout_dependence import TAU_INT, comparator_for_dependence, is_structural_identity  # noqa: F401

CONST_TOL = 1e-12
MIN_N = 50
HEURISTIC_TAU_FRACTION = 0.015


class IdentityViolation(ValueError):
    pass


def longest_finite_run(s):
    ok = np.isfinite(s)
    best, start, cur = (0, 0), None, 0
    for i, v in enumerate(ok):
        if v:
            if start is None:
                start = i
            cur = i - start + 1
            if cur > best[1] - best[0]:
                best = (start, i + 1)
        else:
            start = None
    return best


def classify_series(series, identity=False):
    s = np.asarray(series, float)
    fin = np.isfinite(s)
    out = {'n': int(s.size), 'n_finite': int(fin.sum()), 'n_unavailable': int((~fin).sum())}
    if identity:
        if np.any(np.abs(s[fin]) > CONST_TOL):
            raise IdentityViolation('declared structural identity violated by a finite non-zero value')
        out['class'] = 'structural_identity'
        return out
    a, b = longest_finite_run(s)
    seg = s[a:b]
    out['segment'] = [int(a), int(b)]
    out['n_gaps'] = int(np.sum(np.diff(fin.astype(int)) == -1)) if s.size else 0
    if seg.size < MIN_N:
        out['class'] = 'insufficient_observations'
        return out
    if np.ptp(seg) <= CONST_TOL * (1.0 + abs(float(np.mean(seg)))):
        out.update({'class': 'constant_observed', 'value': float(np.mean(seg))})
        return out
    tau = TAU_INT(seg)
    if tau is None or not np.isfinite(tau):
        out['class'] = 'estimation_failure'
        return out
    out['tau_int_frames'] = float(tau)
    out['heuristic_cutoff_frames'] = HEURISTIC_TAU_FRACTION * seg.size
    out['class'] = ('above_heuristic_cutoff' if tau > out['heuristic_cutoff_frames']
                    else 'estimated_below_heuristic_cutoff')
    return out
