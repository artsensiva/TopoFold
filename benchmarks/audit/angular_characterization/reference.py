"""REF: independent NumPy reference for the scan's documented formulas (protocol §3).

Written from the documented definitions, not by calling or transcribing the Rust code paths:
circular mean as arg(mean exp(iθ)), wrap to (−π, π], two-pass central moments, Sarle's BC
with the documented guards (n < 4, variance ≤ 1e-14, denominator ≤ 1e-12 → 0; clamp [0, 1]).
Agreement with Rust shows the code computes these formulas; it does not validate them.
"""
import numpy as np

from common import GUARD_VAR, wrap


def mean_and_resultant(x, axis=-1):
    z = np.mean(np.exp(1j * np.asarray(x, dtype=np.float64)), axis=axis)
    return np.angle(z), np.abs(z)


def moments(d):
    d = np.asarray(d, dtype=np.float64)
    n = d.size
    c = d - d.mean()
    m2 = np.mean(c ** 2)
    m3 = np.mean(c ** 3)
    m4 = np.mean(c ** 4)
    return n, m2, m3, m4


def bc_from_values(d):
    n, m2, m3, m4 = moments(d)
    out = {'n': n, 'variance': m2, 'skewness': 0.0, 'excess_kurtosis': 0.0, 'bc': 0.0,
           'guard': None}
    if n < 4:
        out['guard'] = 'n<4'
        return out
    if m2 <= GUARD_VAR:
        out['guard'] = 'variance<=1e-14'
        return out
    g = m3 / m2 ** 1.5
    k = m4 / m2 ** 2 - 3.0
    c = 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    out['skewness'], out['excess_kurtosis'] = g, k
    if k + c <= 1e-12:
        out['guard'] = 'denominator<=1e-12'
        return out
    out['bc'] = float(np.clip((g * g + 1.0) / (k + c), 0.0, 1.0))
    return out


def angular_pipeline(torsions):
    """torsions (F, n) -> per-frame means/R, across-frame mean/R, deviations, BC dict."""
    means, R = mean_and_resultant(torsions, axis=1)
    c, Ra = mean_and_resultant(means)
    dev = wrap(means - c)
    return {'means': means, 'R': R, 'across_mean': float(c), 'across_R': float(Ra),
            'deviations': dev, 'bc': bc_from_values(dev)}
