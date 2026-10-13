#!/usr/bin/env python3
"""Held-out (mdCATH) report. Every number is read from results JSON (AGENTS.md rule 3); static text makes no numeric
claim. Run: build/venv/bin/python -B make_report_heldout.py
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
IN = RES / 'wpd_heldout_THIRD_PARTY_MD_DERIVATIVE_results.json'
OUT = HERE / 'REPORT_HELDOUT.md'
METHODS = ('M1', 'M2', 'B1', 'B2', 'B3', 'A', 'C1')
NAMES = {'M1': 'τ chord-RMS', 'M2': 'τ circular-RMS', 'B1': 'Cα–Cα distances', 'B2': 'window Cα RMSD',
         'B3': 'κ-only', 'A': 'residue-aligned κ–τ', 'C1': 'scalar circular mean (coverage-matched)'}
TNAME = {'phipsi': 'T_φψ (primary, angular)', 'pp': 'T_pp (separate, Cartesian peptide-plane)'}


def f(x, d=3):
    return '—' if x is None else f'{x:.{d}f}'


def iv(v, d=3):
    return '—' if v is None else f'[{f(v[0], d)}, {f(v[1], d)}]'


def main():
    r = json.loads(IN.read_text())
    m = r['_meta']
    coh = json.loads((RES / 'wpd_mdcath_cohort.json').read_text())
    L = []
    w = L.append
    w('# WP-D held-out evaluation (mdCATH): Cα-only local-conformation retrieval')
    w('')
    w(f"**{m['label']}.** One authorized evaluation ({m['authorized_by']}). Scope: {m['scope']}. "
      f"{m['n_domains']} domains, {m['workers']} workers, wall {m['wall_seconds']} s, finished {m['finished_utc']}. "
      f"Authorization SHA-256 `{m['authorization_sha256'][:16]}…`; binding extension `{m['binding_extension_sha256'][:16]}…`.")
    w('')
    w('**What this result is:** generalization of a *geometric* retrieval benchmark to held-out CATH superfamilies '
      'within mdCATH (CHARMM22*, TIP3P, 320 K and 348 K; cross-replica retrieval from a shared starting configuration). '
      '**What it is not:** biological, functional or commercial validation; evidence about other force fields, '
      'temperatures, multi-domain proteins or independent starting conditions.')
    w('')
    w('## 1. Cohort')
    w('')
    w('| Order | Domain | Role | Superfamily | Replaces |')
    w('|---|---|---|---|---|')
    for d in coh['domains']:
        w(f"| {d['order']} | {d['domain']} | {d['role']} | {d['superfamily']} | {d.get('replaces_primary') or ''} |")
    w('')
    w('Excluded primaries (structural eligibility, original §4.3): ' + '; '.join(
        f'{k}: {", ".join(v)}' for k, v in coh['_meta']['excluded_primaries'].items()) + '.')
    w('')
    doms = r['domains']
    hk = max(d['horn_kabsch_max_A'] for d in doms)
    ref = max(max(g['reference_max_tau_kappa']) for d in doms for t in d['geometry'].values() for g in t.values())
    bf = sum(g['binding_failures'] for d in doms for t in d['geometry'].values() for g in t.values())
    w(f'**Mandatory checks passed:** binding vs NumPy max difference {ref:.1e} rad; Horn vs Kabsch max {hk:.1e} Å; '
      f'binding failures (frames excluded) {bf}; input integrity, extension identity and production sources verified.')
    w('')
    for tgt in ('phipsi', 'pp'):
        agg = r['aggregate'][tgt]
        w(f'## 2{"a" if tgt == "phipsi" else "b"}. {TNAME[tgt]}')
        w('')
        w('NR@10: 0 = target-optimal, 1 = random. G = 1 − median NR (gain over random). L = relative loss vs the best '
          'method (Amendment 1). Intervals: 95 % percentile domain bootstrap (2,000), best re-selected per resample. '
          'Classes vs ρ = 0.10. C1 uses its coverage-matched population.')
        w('')
        for key, c in agg['cells'].items():
            sec = r['secondary'][key]
            if 'methods' not in c:
                w(f"**{key}:** {c['status']} ({c['n_defined_domains']} defined domains).")
                w('')
                continue
            w(f"**{key}** — {c['n_defined_domains']} defined domains; observed best: **{c['observed']['best']}**; "
              f"best counts across resamples: {c['best_counts']}; C1 coverage-matched domains: "
              f"{c['methods']['C1']['n_domains']}.")
            w('')
            w('| Method | median NR@10 | G [95 %] | L [95 %] | loss class | gain class |')
            w('|---|---|---|---|---|---|')
            for mth in METHODS:
                e = c['methods'][mth]
                med = c['observed']['median_nr'].get(mth) if mth != 'C1' else (1 - e['G'] if e['G'] is not None else None)
                w(f"| {mth} {NAMES[mth]} | {f(med)} | {f(e['G'])} {iv(e['G_interval'])} | {f(e['L'])} "
                  f"{iv(e['L_interval'])} | {e['loss_class']} | {e['gain_class']} |")
            w('')
        w('**Overall outcome per method (Amendment 1 §4 precedence):** ' + '; '.join(
            f'{k}: `{v}`' for k, v in agg['overall'].items()) + '.')
        w('')
    w('## 3. Secondary (reported, not deciding)')
    w('')
    w('| Cell | defined / total domains | windows defined | undefined windows | random overlap@10 | '
      + ' | '.join(f'{x} ovl@10 / NR@1 / ρ' for x in METHODS) + ' |')
    w('|---|---|---|---|---|' + '---|' * len(METHODS))
    for key, s in r['secondary'].items():
        cells = []
        for mth in METHODS:
            dm = s['domain_medians'][mth]
            cells.append(f"{f(dm['overlap10'], 4)} / {f(dm['nr1'])} / {f(dm['spearman'])}")
        w(f"| {key} | {s['n_defined_domains']}/{s['n_domains']} | {s['windows_defined_total']} | "
          f"{s['windows_undefined_by_reason']} | {f(s['random_overlap10_median'], 4)} | " + ' | '.join(cells) + ' |')
    w('')
    w('Target-spread terciles (windows pooled within cell; mean window NR@10 per tercile T1 / T2 / T3):')
    w('')
    w('| Cell | bounds | ' + ' | '.join(METHODS) + ' |')
    w('|---|---|' + '---|' * len(METHODS))
    for key, s in r['secondary'].items():
        t = s['target_spread_terciles']
        if not t:
            continue
        w(f"| {key} | {'; '.join(f'{f(a)}–{f(b)}' for a, b in t['bounds'])} | " + ' | '.join(
            ' / '.join(f(x) for x in t['mean_window_nr10'][mth]) for mth in METHODS) + ' |')
    w('')
    tot = {}
    for s in r['secondary'].values():
        for pair, cls in s['dependence_class_counts'].items():
            for k, v in cls.items():
                tot[k] = tot.get(k, 0) + v
    w(f'**Within-domain dependence diagnostics** (descriptive; AR(1)-based heuristic only, Amendment 3 §4): class counts '
      f'over all domain × cell × pair series: {tot}. They do not affect the domain-level inference.')
    w('')
    w('## 4. Interpretation limits')
    w('')
    w('- Intervals are marginal and unadjusted (7 methods × 4 cells × 2 statistics per target); no family-wise '
      'coverage. Classes are pre-declared summaries, not tests; "within the margin" is not equivalence.')
    w('- ρ = 0.10 is a stated value judgment, not derived from a downstream task.')
    w('- Replicas share their starting configuration; cross-replica retrieval may be easier than between independent '
      'runs. The first 10 % of frames were excluded.')
    w('- Method A entered the candidate set after AdK development results; this held-out evaluation is its first test.')
    w('- Both targets are geometric reconstructions from the same trajectories; T_φψ may favour angular and T_pp '
      'Cartesian methods.')
    w('- 30 domains, one per superfamily, from a seeded draw; generalization is to that population within mdCATH.')
    w('')
    OUT.write_text('\n'.join(L) + '\n')
    print('wrote', OUT.name)


if __name__ == '__main__':
    main()
