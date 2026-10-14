"""BQ-3 failure-path tests of the acceptance function (SYNTHETIC inputs; no data, no timing).
Run: ../wpd_local_retrieval/build/venv/bin/python -B -m pytest -p no:cacheprovider test_bq3.py
"""
import copy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bq3_lib as Q  # noqa: E402

H = ['a' * 64, 'b' * 64, 'c' * 64]
REF = {'status': 'ok', 'finite': True}
GATE = {'status': 'ok', 'candidate_finite': True, 'max_abs_dist_err_A': 1e-9, 'rank_nonequivalent': 0,
        'topk_sha256': H}


def rep(i, q=(5.0, 5.0, 5.0), s=9.0, rss=3_000_000, hashes=H):
    return {'rep': i, 'status': 'ok', 'query_wall_s': list(q), 'first_query_span_s': s, 'peak_rss_kB': rss,
            'topk_sha256': list(hashes)}


REPS = [rep(i) for i in range(5)]


def test_all_pass():
    d = Q.decide(REF, GATE, REPS, 5)
    assert d['prepared_database'] == 'qualified' and d['first_query'] == 'qualified'


def test_missing_or_failed_gate_and_reference_block_qualification():
    for ref, gate in ((None, GATE), (REF, None), ({'status': 'failed'}, GATE), (REF, {'status': 'failed'}),
                      ({'status': 'ok', 'finite': False}, GATE)):
        d = Q.decide(ref, gate, REPS, 5)
        assert d['prepared_database'] == d['first_query'] == 'not_qualified_verification'


def test_gate_content_failures_block_qualification():
    for k, v in (('candidate_finite', False), ('max_abs_dist_err_A', 2e-6), ('max_abs_dist_err_A', float('nan')),
                 ('max_abs_dist_err_A', None), ('rank_nonequivalent', 1), ('rank_nonequivalent', None),
                 ('topk_sha256', H[:2]), ('topk_sha256', None)):
        g = copy.deepcopy(GATE)
        g[k] = v
        assert Q.decide(REF, g, REPS, 5)['prepared_database'] == 'not_qualified_verification', (k, v)


def test_timed_output_mismatch_blocks_qualification():
    reps = copy.deepcopy(REPS)
    reps[3]['topk_sha256'][1] = 'd' * 64
    d = Q.decide(REF, GATE, reps, 5)
    assert d['prepared_database'] == 'not_qualified_verification' and 'repetition 3' in d['verification_reasons'][0]


def test_missing_or_failed_repetition_is_insufficient_evidence():
    assert Q.decide(REF, GATE, REPS[:4], 5)['prepared_database'] == 'not_qualified_insufficient_evidence'
    reps = copy.deepcopy(REPS)
    reps[2] = {'rep': 2, 'status': 'skipped_memory_rule'}
    assert Q.decide(REF, GATE, reps, 5)['first_query'] == 'not_qualified_insufficient_evidence'


def test_limits_evaluated_separately_and_per_observation():
    reps = copy.deepcopy(REPS)
    reps[4]['first_query_span_s'] = 10.5                     # one cold sequence over the limit
    d = Q.decide(REF, GATE, reps, 5)
    assert d['prepared_database'] == 'qualified' and d['first_query'] == 'not_qualified_measured_limit'
    reps = copy.deepcopy(REPS)
    reps[0]['query_wall_s'][2] = 10.01                       # one query over the limit
    d = Q.decide(REF, GATE, reps, 5)
    assert d['prepared_database'] == 'not_qualified_measured_limit' and d['first_query'] == 'qualified'
    reps = copy.deepcopy(REPS)
    reps[1]['peak_rss_kB'] = Q.R4_KB + 1
    d = Q.decide(REF, GATE, reps, 5)
    assert d['prepared_database'] == d['first_query'] == 'not_qualified_measured_limit'
    reps = copy.deepcopy(REPS)
    reps[1]['query_wall_s'] = [5.0, 5.0]                     # a missing query latency: incomplete evidence
    assert Q.decide(REF, GATE, reps, 5)['prepared_database'] == 'not_qualified_insufficient_evidence'
