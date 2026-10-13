"""Tests for the metadata-only held-out selection on SYNTHETIC metadata (no mdCATH content).

Run: build/venv/bin/python -B -m pytest -p no:cacheprovider test_mdcath_select.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mdcath_select import draw, metadata_pool, schema_sample  # noqa: E402


def synthetic():
    cath, tree, doms = {}, {}, []
    for s in range(12):
        for j in range(s % 3 + 1):
            d = f'x{s:02d}{j}A00'
            doms.append(d)
            cath[d] = {'superfamily': f'1.{s}.1.{j % 2 + 1}' if s < 6 else f'2.{s}.5.1', 'length': 40 + 30 * s}
            tree[d] = {'size': 1000 + 7 * s + j, 'lfs': {'oid': f'h{d}'}}
    doms.append('missingA00')
    return doms, cath, tree


def test_metadata_pool_filters_length_and_presence():
    doms, cath, tree = synthetic()
    pool, fail = metadata_pool(doms, cath, tree)
    assert fail['not_in_tree'] == 1
    assert all(50 <= pool[d]['length'] <= 300 for d in pool)
    assert fail['length_out_of_range'] == sum(1 for d in doms if d in cath and not 50 <= cath[d]['length'] <= 300)


def test_schema_sample_is_smallest_then_name():
    doms, cath, tree = synthetic()
    pool, _ = metadata_pool(doms, cath, tree)
    s = schema_sample(pool)
    assert pool[s]['size_bytes'] == min(p['size_bytes'] for p in pool.values())


def test_draw_deterministic_one_per_superfamily_and_excludes():
    doms, cath, tree = synthetic()
    pool, _ = metadata_pool(doms, cath, tree)
    excl = {pool[schema_sample(pool)]['superfamily']}
    a, n_sf, _ = draw(pool, excl, seed=1, n=5)
    b, _, _ = draw(pool, excl, seed=1, n=5)
    c, _, _ = draw(pool, excl, seed=2, n=5)
    assert a == b and len(a) == 5
    sfs = [pool[d]['superfamily'] for d in a]
    assert len(set(sfs)) == 5 and not set(sfs) & excl
    assert a != c
    full, n_sf2, _ = draw(pool, excl, seed=1, n=10 ** 6)
    assert len(full) == n_sf2 == n_sf


def test_draw_order_independent_of_input_order():
    doms, cath, tree = synthetic()
    pool, _ = metadata_pool(doms, cath, tree)
    rev = dict(reversed(list(pool.items())))
    assert draw(pool, set(), seed=3, n=6)[0] == draw(rev, set(), seed=3, n=6)[0]
