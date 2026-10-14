"""BQ-3 acceptance logic (docs/B2_DIRECT_P2_QUALIFICATION_PROTOCOL.md §3). Pure functions; tested before execution.

decide(ref, gate, reps, n_reps) returns separate verdicts for prepared-database and first-query operation. Any
missing, failed or mismatched verification check, or an incomplete set of repetitions, prevents qualification.
"""
R2_S = 10.0
R4_KB = 8 * 1024 ** 2
TOL_D = 1e-6
N_QUERIES = 3


def verification(ref, gate, reps):
    """Returns (ok, reasons). ref/gate: child result dicts or None; reps: list of repetition dicts."""
    why = []
    if not ref or ref.get('status') != 'ok':
        why.append('reference missing or failed')
    elif not ref.get('finite'):
        why.append('reference not finite')
    if not gate or gate.get('status') != 'ok':
        why.append('gate missing or failed')
    else:
        if not gate.get('candidate_finite'):
            why.append('candidate distances not finite')
        err = gate.get('max_abs_dist_err_A')
        if err is None or not err <= TOL_D:
            why.append(f'distance error {err} exceeds {TOL_D} Å')
        if gate.get('rank_nonequivalent') != 0:
            why.append(f"rank non-equivalent query-windows: {gate.get('rank_nonequivalent')}")
        hashes = gate.get('topk_sha256') or []
        if len(hashes) != N_QUERIES:
            why.append('gate top-k hashes missing')
        for r in reps:
            if r.get('status') == 'ok' and r.get('topk_sha256') != hashes:
                why.append(f"repetition {r.get('rep')}: timed outputs differ from verified outputs")
    return not why, why


def decide(ref, gate, reps, n_reps):
    ok_v, why_v = verification(ref, gate, reps)
    done = [r for r in reps if r.get('status') == 'ok' and len(r.get('query_wall_s') or []) == N_QUERIES
            and r.get('first_query_span_s') is not None and r.get('peak_rss_kB') is not None]
    complete = len(done) == n_reps and len(reps) == n_reps
    why_c = [] if complete else [f'{len(done)} of {n_reps} repetitions completed']
    mem_ok = bool(done) and all(r['peak_rss_kB'] <= R4_KB for r in done)
    q_all = [q for r in done for q in r['query_wall_s']]
    s_all = [r['first_query_span_s'] for r in done]
    prepared_pass = bool(q_all) and all(q <= R2_S for q in q_all)
    first_pass = bool(s_all) and all(s <= R2_S for s in s_all)

    def verdict(passed):
        if not ok_v:
            return 'not_qualified_verification'
        if not complete:
            return 'not_qualified_insufficient_evidence'
        return 'qualified' if (passed and mem_ok) else 'not_qualified_measured_limit'
    return {'verification_ok': ok_v, 'verification_reasons': why_v, 'complete': complete,
            'completeness_reasons': why_c, 'memory_within_R4_all_reps': mem_ok,
            'prepared_database': verdict(prepared_pass), 'first_query': verdict(first_pass),
            'n_query_latencies': len(q_all), 'n_first_query_spans': len(s_all)}
