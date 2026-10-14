"""CMP-1 decision logic (protocol §7). Pure functions; failure paths tested in test_cmp1.py.

Per workload w and accuracy track T:
  eligible(i, w, T)  ⇔ gate on w exists and passes in T (all queries) ∧ 5 complete timing repetitions ∧ every timed
                       top-k (indices and distances) equals the gate's.
  E*(w, T)           = eligible EXTERNAL implementation with the lowest median repeated-query latency (Q2–Q3 pooled);
                       ties: lower median incremental peak RSS, then name. Undefined if none is eligible.
  outcome(w, T)      ∈ {V2b_advantage, comparator_preferable, no_decisive_difference, insufficient_evidence}.
A competitor's failure (accuracy track, degeneracy) is a separately listed numerical limitation; it never by itself
produces V2b_advantage or justifies integration.
"""
import numpy as np

ADV_RATIO = 0.67          # V2b median ≤ 0.67 × E* median (≥ 1.5× faster)
PAIRED_MIN = 1.25         # every per-round paired ratio E*/V2b ≥ 1.25
MEM_RATIO = 1.25          # V2b incremental peak RSS ≤ 1.25 × E*'s
R4_KB = 8 * 1024 ** 2
N_REPS = 5
EXTERNAL = ('E-MDT', 'E-MDA')
CANDIDATE = 'T-V2b'


def med(x):
    return float(np.median(x))


def eligible(cell):
    """cell: {'gate': {'status', 'tracks': {T: {'pass'}}, 'topk_sha256'}, 'reps': [...]} for one (impl, workload)."""
    g = cell.get('gate') or {}
    reps = [r for r in cell.get('reps', []) if r.get('status') == 'ok']
    out = {}
    for T in ('F64', 'F32'):
        why = []
        if g.get('status') != 'ok':
            why.append('gate missing or failed')
        elif not (g.get('tracks') or {}).get(T, {}).get('pass'):
            why.append(f'gate fails track {T}')
        if len(reps) != N_REPS or len(cell.get('reps', [])) != N_REPS:
            why.append(f'{len(reps)} of {N_REPS} timing repetitions complete')
        if g.get('status') == 'ok' and any(r.get('topk_sha256') != g.get('topk_sha256') for r in reps):
            why.append('timed outputs differ from gate outputs')
        if any(len(r.get('query_wall_s') or []) != 3 for r in reps):
            why.append('missing query latencies')
        out[T] = (not why, why)
    return out


def summary(cell):
    reps = [r for r in cell['reps'] if r.get('status') == 'ok']
    return {'rep_q': {r['round']: med(r['query_wall_s'][1:]) for r in reps},
            'median_repeated_query_s': med([q for r in reps for q in r['query_wall_s'][1:]]),
            'median_incr_peak_kB': med([r['peak_rss_kB'] - r['rss_after_setup_kB'] for r in reps]),
            'max_peak_kB': max(r['peak_rss_kB'] for r in reps)}


def outcome(cells, T):
    """cells: {impl: cell} for one workload. Returns dict with E*, outcome and reasons."""
    el = {i: eligible(c)[T] for i, c in cells.items()}
    res = {'eligibility': {i: {'eligible': e[0], 'reasons': e[1]} for i, e in el.items()}}
    ext = [i for i in EXTERNAL if i in cells and el[i][0]]
    if not ext:
        res.update({'E_star': None, 'outcome': 'insufficient_evidence', 'reason': 'no eligible external comparator'})
        return res
    sm = {i: summary(cells[i]) for i in ext}
    estar = sorted(ext, key=lambda i: (sm[i]['median_repeated_query_s'], sm[i]['median_incr_peak_kB'], i))[0]
    res['E_star'] = estar
    if CANDIDATE not in cells or not el[CANDIDATE][0]:
        res.update({'outcome': 'insufficient_evidence', 'reason': 'candidate not eligible'})
        return res
    v, e = summary(cells[CANDIDATE]), sm[estar]
    rounds = sorted(set(v['rep_q']) & set(e['rep_q']))
    paired = [e['rep_q'][r] / v['rep_q'][r] for r in rounds]
    ratio = v['median_repeated_query_s'] / e['median_repeated_query_s']
    mem_ok = v['max_peak_kB'] <= R4_KB and v['median_incr_peak_kB'] <= MEM_RATIO * e['median_incr_peak_kB']
    res.update({'median_ratio_V2b_over_Estar': ratio, 'paired_ratios_Estar_over_V2b': paired,
                'memory_ok': bool(mem_ok), 'V2b': v, 'Estar_summary': e})
    if len(rounds) != N_REPS:
        res.update({'outcome': 'insufficient_evidence', 'reason': 'unpaired rounds'})
    elif ratio <= ADV_RATIO and min(paired) >= PAIRED_MIN and mem_ok:
        res['outcome'] = 'V2b_advantage'
    elif e['median_repeated_query_s'] <= v['median_repeated_query_s']:
        res['outcome'] = 'comparator_preferable'
    else:
        res['outcome'] = 'no_decisive_difference'
    return res


def conclusion(per_workload, required, v2b_degeneracy_pass_F64):
    """per_workload: {wid: outcome dict for the PRIMARY track F32}."""
    outs = {w: per_workload.get(w, {}).get('outcome', 'insufficient_evidence') for w in required}
    if any(o == 'insufficient_evidence' for o in outs.values()):
        return 'INSUFFICIENT_EVIDENCE', outs
    if all(o == 'V2b_advantage' for o in outs.values()) and v2b_degeneracy_pass_F64:
        return 'INTEGRATION_SUPPORTED', outs
    if all(o == 'comparator_preferable' for o in outs.values()):
        return 'EXISTING_LIBRARY_PREFERABLE', outs
    return 'NO_DECISIVE_ADVANTAGE', outs
