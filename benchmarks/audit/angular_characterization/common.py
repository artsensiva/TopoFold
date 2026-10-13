"""Shared helpers for the SYNTHETIC angular-characterization harness (WP-A).

All angles are radians. Raw arrays are little-endian float64, C order.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
BUILD = HERE / 'build'                       # regenerable binaries (not for version control)
FIXTURE_BIN = BUILD / 'fixtures'
ARRAY_BIN = BUILD / 'arrays'
FIXTURES_JSON = HERE / 'fixtures' / 'SYNTHETIC_fixtures.json'
FIXTURES_TSV = FIXTURE_BIN / 'SYNTHETIC_fixtures.tsv'   # flat manifest read by the Rust example
RESULTS = HERE / 'results'
RUST_JSON = RESULTS / 'SYNTHETIC_angular_characterization_rust.json'
PY_JSON = RESULTS / 'SYNTHETIC_angular_characterization_python.json'
SUMMARY_JSON = RESULTS / 'SYNTHETIC_angular_characterization_summary.json'
PROVENANCE_JSON = RESULTS / 'SYNTHETIC_angular_characterization_provenance.json'

U = 2.0 ** -53                # unit roundoff, float64
BOND = 3.8                    # virtual Cα–Cα bond length, Å
KAPPA = 1.2                   # constant turning angle, rad
UNDEFINED_BOUND = 0.1         # rad; conditioning diagnostic (protocol §4), not a production value
REL_PERTURB_MAX = 1e-6        # protocol §4 BC conditioning criterion
GUARD_VAR = 1e-14             # production variance guard (bimodality.rs), recorded for comparison
R2_RECOVERY_TOL = 1e-9        # protocol §4
P1_RECOVERY_FLAG = 1e-5       # protocol §4 (flag only)


def wrap(a):
    """Wrap to (−π, π]."""
    a = np.asarray(a, dtype=np.float64)
    w = np.mod(a + np.pi, 2 * np.pi) - np.pi
    return np.where(w <= -np.pi, w + 2 * np.pi, w)


def circ_dist(a, b):
    return np.abs(wrap(np.asarray(a) - np.asarray(b)))


def direction_bound(n, R, delta_in=0.0):
    """Protocol §4 conditioning bound b = ((n+1)u + δ_in) / R (inf if R == 0)."""
    R = np.asarray(R, dtype=np.float64)
    num = (np.asarray(n, dtype=np.float64) + 1.0) * U + delta_in
    with np.errstate(divide='ignore'):
        return np.where(R > 0, num / np.where(R > 0, R, 1.0), np.inf)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def write_f64(path, arr):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.ascontiguousarray(arr, dtype='<f8').tofile(path)
    return sha256_file(path)


def read_f64(path, shape=None):
    a = np.fromfile(path, dtype='<f8')
    return a.reshape(shape) if shape is not None else a


def load_json(path):
    with open(path) as f:
        return json.load(f)


def dump_json(obj, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w') as f:
        json.dump(obj, f, indent=1, sort_keys=False, allow_nan=False, default=_json_default)
        f.write('\n')


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return finite_or_none(float(o))
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return [_json_default(x) if isinstance(x, (np.generic,)) else x for x in o.tolist()]
    raise TypeError(type(o))


def finite_or_none(x):
    """JSON has no NaN/inf; non-finite values are stored as null with an explicit flag."""
    return x if math.isfinite(x) else None
