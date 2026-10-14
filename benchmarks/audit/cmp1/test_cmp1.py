"""CMP-1 tests: decision-logic failure paths and accuracy-track functions (SYNTHETIC; no timing).
Run: build/venv/bin/python -B -m pytest -p no:cacheprovider test_cmp1.py
"""
import copy
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cmp1_accuracy as ACC  # noqa: E402
import cmp1_criteria as C  # noqa: E402

H = ['h'] * 3


def cell(q, incr=100_000, peak=2_000_000, tracks=('F64', 'F32'), status='ok', n=5, hashes=H):
    g = {'status': status, 'tracks': {t: {'pass': t in tracks} for t in ('F64', 'F32')}, 'topk_sha256': H}
    reps = [{'status': 'ok', 'round': r, 'query_wall_s': [q, q, q], 'peak_rss_kB': peak,
             'rss_after_setup_kB': peak - incr, 'topk_sha256': list(hashes)} for r in range(n)]
    return {'gate': g, 'reps': reps}


def test_advantage_requires_speed_paired_rounds_and_memory():
    o = C.outcome({'T-V2b': cell(1.0), 'E-MDT': cell(2.0), 'E-MDA': cell(3.0)}, 'F32')
    assert o['E_star'] == 'E-MDT' and o['outcome'] == 'V2b_advantage'
    o = C.outcome({'T-V2b': cell(1.0, incr=300_000), 'E-MDT': cell(2.0)}, 'F32')
    assert o['outcome'] == 'no_decisive_difference' and not o['memory_ok']
    o = C.outcome({'T-V2b': cell(1.0), 'E-MDT': cell(1.2)}, 'F32')
    assert o['outcome'] == 'no_decisive_difference'
    o = C.outcome({'T-V2b': cell(1.0), 'E-MDT': cell(0.8)}, 'F32')
    assert o['outcome'] == 'comparator_preferable'


def test_competitor_failure_alone_never_gives_advantage():
    o = C.outcome({'T-V2b': cell(1.0), 'E-MDT': cell(0.5, tracks=('F32',))}, 'F64')
    assert o['E_star'] is None and o['outcome'] == 'insufficient_evidence'


def test_estar_is_fastest_eligible_external_per_track():
    cells = {'T-V2b': cell(1.0), 'E-MDT': cell(0.5, tracks=('F32',)), 'E-MDA': cell(3.0)}
    assert C.outcome(cells, 'F32')['E_star'] == 'E-MDT'
    assert C.outcome(cells, 'F64')['E_star'] == 'E-MDA'


def test_missing_failed_or_mismatched_evidence():
    for bad in (cell(1.0, status='failed'), cell(1.0, n=4), cell(1.0, hashes=['x'] * 3), cell(1.0, tracks=())):
        assert C.outcome({'T-V2b': bad, 'E-MDT': cell(2.0)}, 'F32')['outcome'] == 'insufficient_evidence'
    c = cell(1.0)
    c['reps'][2]['query_wall_s'] = [1.0, 1.0]
    assert C.outcome({'T-V2b': c, 'E-MDT': cell(2.0)}, 'F32')['outcome'] == 'insufficient_evidence'


def test_conclusion_categories():
    adv = {'outcome': 'V2b_advantage'}
    pref = {'outcome': 'comparator_preferable'}
    assert C.conclusion({'a': adv, 'b': adv}, ['a', 'b'], True)[0] == 'INTEGRATION_SUPPORTED'
    assert C.conclusion({'a': adv, 'b': adv}, ['a', 'b'], False)[0] == 'NO_DECISIVE_ADVANTAGE'
    assert C.conclusion({'a': pref, 'b': pref}, ['a', 'b'], True)[0] == 'EXISTING_LIBRARY_PREFERABLE'
    assert C.conclusion({'a': adv}, ['a', 'b'], True)[0] == 'INSUFFICIENT_EVIDENCE'
    assert C.conclusion({'a': adv, 'b': pref}, ['a', 'b'], True)[0] == 'NO_DECISIVE_ADVANTAGE'


def test_track_tolerances_identical_for_all_and_rank_check():
    d_ref = np.arange(20.0)
    d = d_ref + 5e-4
    idx = np.argsort(d_ref)[:10]
    ok32 = ACC.check_window('F32', d, d_ref, idx)
    ok64 = ACC.check_window('F64', d, d_ref, idx)
    assert ok32[0] and not ok64[0] and ok32[2]
    bad = d_ref.copy()
    bad[3] = np.nan
    assert not ACC.check_window('F32', bad, d_ref, idx)[0]
    swapped = idx.copy()
    swapped[9] = 15
    assert not ACC.check_window('F32', d_ref, d_ref, swapped)[2]
