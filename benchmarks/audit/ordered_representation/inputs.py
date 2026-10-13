"""Read-only access to the frozen WP-A fixtures, with hash verification (protocol §4).

Files are read from the WP-A build/ directory and checked against the frozen manifest. If a file
is absent, the case is regenerated IN MEMORY with the WP-A generator functions and must reproduce
the manifest hash. Nothing is written into the WP-A directory.
"""
import hashlib
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
WPA = REPO / 'benchmarks' / 'audit' / 'angular_characterization'
sys.path.insert(0, str(WPA))

import build_traces as wpa_build  # noqa: E402  (WP-A, read-only use)
import common as wpa_common  # noqa: E402
import generate_fixtures as wpa_gen  # noqa: E402

WPA_MANIFEST_PATH = WPA / 'fixtures' / 'SYNTHETIC_fixtures.json'
WPA_FIXTURE_DIR = WPA / 'build' / 'fixtures'


class InputHashError(RuntimeError):
    pass


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


def arr_sha(a):
    return sha_bytes(np.ascontiguousarray(a, dtype='<f8').tobytes())


@lru_cache(maxsize=1)
def manifest():
    return {c['id']: c for c in wpa_common.load_json(WPA_MANIFEST_PATH)['cases']}


@lru_cache(maxsize=1)
def _regenerated():
    m = wpa_gen.Matrix()
    b1 = wpa_gen.f1(m)
    wpa_gen.f2(m, b1)
    wpa_gen.f3(m)
    wpa_gen.f4(m)
    wpa_gen.f6(m)
    return {c['id']: c for c in m.cases}


def read_verified(path, expected_sha, shape):
    path = Path(path)
    data = path.read_bytes()
    if sha_bytes(data) != expected_sha:
        raise InputHashError(f'{path}: SHA-256 does not match frozen manifest')
    return np.frombuffer(data, dtype='<f8').reshape(shape).copy()


def load_case(cid, fixture_dir=WPA_FIXTURE_DIR, allow_regenerate=True):
    """Return (torsions (F, m), coords (F, W, 3) or None, source) for a WP-A case."""
    c = manifest()[cid]
    F, m, W = c['F'], c['ncols'], c['W']
    tpath = Path(fixture_dir) / c['torsion_file']
    if tpath.exists():
        tors = read_verified(tpath, c['torsion_sha256'], (F, m))
        coords = None
        if c['coords_file']:
            coords = read_verified(Path(fixture_dir) / c['coords_file'], c['coords_sha256'],
                                   (F, W, 3))
        return tors, coords, 'wpa_build_file'
    if not allow_regenerate:
        raise FileNotFoundError(tpath)
    tors = np.asarray(_regenerated()[cid]['tors'], dtype=np.float64)
    if arr_sha(tors) != c['torsion_sha256']:
        raise InputHashError(f'{cid}: in-memory regeneration does not reproduce manifest hash')
    coords = None
    if c['coords_file']:
        coords = wpa_build.build_traces(np.full((F, W - 2), wpa_common.KAPPA), tors)
        if arr_sha(coords) != c['coords_sha256']:
            raise InputHashError(f'{cid}: regenerated coordinates do not match manifest hash')
    return tors, coords, 'in_memory_regeneration'


def f1_labels(cid, fixture_dir=WPA_FIXTURE_DIR):
    c = manifest()[cid]
    path = Path(fixture_dir) / c['labels_file']
    if path.exists():
        return read_verified(path, c['labels_sha256'], (c['F'],)).astype(int)
    labels = np.arange(c['F']) % 2
    if arr_sha(labels.astype(np.float64)) != c['labels_sha256']:
        raise InputHashError(f'{cid}: regenerated labels do not match manifest hash')
    return labels


def f4_labels(cid, tors):
    """Reproduce the φ₂ assignment from the recorded seed and WP-A draw order; verified by
    reproducing the torsion hash (protocol §4)."""
    c = manifest()[cid]
    p = c['params']
    rng = np.random.default_rng(c['seed'])
    perm = rng.permutation(c['F'])
    state = np.zeros(c['F'], dtype=bool)
    state[perm[:p['n_phi2']]] = True
    nu = rng.normal(0.0, p['sigma_rad'], size=c['F']) if p['sigma_rad'] > 0 else np.zeros(c['F'])
    x = np.where(state, p['phi2_rad'], 0.0) + nu
    rebuilt = wpa_common.wrap(np.stack([x, x], axis=1))
    if c['alpha_deg'] != 0:
        rebuilt = wpa_common.wrap(rebuilt + np.deg2rad(c['alpha_deg']))
    if arr_sha(rebuilt) != c['torsion_sha256'] or not np.array_equal(rebuilt, tors):
        raise InputHashError(f'{cid}: label reproduction failed (torsion hash mismatch)')
    return state.astype(int)
