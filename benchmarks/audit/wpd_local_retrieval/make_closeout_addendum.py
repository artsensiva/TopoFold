#!/usr/bin/env python3
"""Held-out closeout addendum (2026-10-12): interpretation corrections only. Reads the frozen held-out results; writes
results/wpd_heldout_closeout_addendum.json and docs/WPD_HELDOUT_MDCATH_CLOSEOUT_ADDENDUM_2026-10-12.md.
No result, protocol, implementation or checkpoint is modified. Every number in the addendum comes from the results JSON.
Run: build/venv/bin/python -B make_closeout_addendum.py
"""
import datetime
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
RES = HERE / 'results'
IN = RES / 'wpd_heldout_THIRD_PARTY_MD_DERIVATIVE_results.json'
RUN = RES / 'wpd_heldout_run_record.json'
OUT_JSON = RES / 'wpd_heldout_closeout_addendum.json'
OUT_MD = REPO / 'docs' / 'WPD_HELDOUT_MDCATH_CLOSEOUT_ADDENDUM_2026-10-12.md'
RHO = 0.10


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def f(x, d=3):
    return '—' if x is None else f'{x:.{d}f}'


def iv(v):
    return f'[{f(v[0])}, {f(v[1])}]'


def main():
    run = json.loads(RUN.read_text())
    if sha256(IN) != run['results_sha256']:
        raise SystemExit('STOP: results file differs from the run record')
    r = json.loads(IN.read_text())
    ph = r['aggregate']['phipsi']['cells']
    pp = r['aggregate']['pp']['cells']
    rec = {'_meta': {'producer': 'make_closeout_addendum.py', 'written_utc': datetime.datetime.now(datetime.UTC).isoformat(),
                     'results_sha256': sha256(IN), 'run_record_sha256': sha256(RUN), 'rho': RHO,
                     'nature': 'interpretation corrections; no rerun, no change to results or declared outcomes'},
           'B2_phipsi_W8': {k: {'L': ph[k]['methods']['B2']['L'], 'L_interval': ph[k]['methods']['B2']['L_interval'],
                                'loss_class': ph[k]['methods']['B2']['loss_class'],
                                'interval_contains_rho': ph[k]['methods']['B2']['L_interval'][0] <= RHO
                                <= ph[k]['methods']['B2']['L_interval'][1]}
                            for k in ('320/W8/phipsi', '348/W8/phipsi')},
           'W5_phipsi_A_B2': {k: {m: {'L': ph[k]['methods'][m]['L'], 'L_interval': ph[k]['methods'][m]['L_interval'],
                                      'loss_class': ph[k]['methods'][m]['loss_class']} for m in ('A', 'B2')}
                              | {'observed_best': ph[k]['observed']['best'], 'best_counts': ph[k]['best_counts']}
                              for k in ('320/W5/phipsi', '348/W5/phipsi')},
           'C1_gain': {k: {'G': c['methods']['C1']['G'], 'G_interval': c['methods']['C1']['G_interval'],
                           'gain_class': c['methods']['C1']['gain_class'], 'loss_class': c['methods']['C1']['loss_class'],
                           'G_lower_bound_above_0': c['methods']['C1']['G_interval'][0] > 0}
                       for cells in (ph, pp) for k, c in cells.items()},
           'declared_overall': {t: r['aggregate'][t]['overall'] for t in ('phipsi', 'pp')}}
    OUT_JSON.write_text(json.dumps(rec, indent=1) + '\n')
    b2 = rec['B2_phipsi_W8']
    w5 = rec['W5_phipsi_A_B2']
    c1 = rec['C1_gain']
    L = []
    w = L.append
    w('# WP-D Held-out (mdCATH): Closeout Addendum, 2026-10-12')
    w('')
    w('**Nature.** These are interpretation corrections only.')
    w('- No rerun, no tuning.')
    w('- No change to protocols, the scoring implementation, raw results, checkpoints or declared outcomes.')
    w(f"- Values are read by `make_closeout_addendum.py` from the frozen results (SHA-256 "
      f"`{rec['_meta']['results_sha256'][:16]}…`, matching `results/wpd_heldout_run_record.json`).")
    w('- Machine-readable record: `benchmarks/audit/wpd_local_retrieval/results/wpd_heldout_closeout_addendum.json`.')
    w('')
    w('The corrected statements were made in the assistant\'s chat summary of the run (2026-10-06). The generated report '
      '`REPORT_HELDOUT.md` contains the classes and intervals and is unchanged.')
    w('')
    w('## C1. B2 at W = 8 on T_φψ: inconclusive, not "beyond the margin"')
    w('')
    w('The summary said that at W = 8 method A beats B2 "with relative-loss intervals that clear the margin". That is '
      'wrong for B2. Its relative-loss intervals cross ρ = 0.10, so its declared class there is **inconclusive**:')
    w('')
    w('| Cell | L(B2) | 95 % interval | Contains ρ | Declared class |')
    w('|---|---|---|---|---|')
    for k, e in b2.items():
        w(f"| {k} | {f(e['L'])} | {iv(e['L_interval'])} | {e['interval_contains_rho']} | {e['loss_class']} |")
    w('')
    w('What is supported is narrower:')
    w('- A was the observed best in these cells and is within the margin.')
    w('- Whether B2 loses more than the declared acceptable fraction of A\'s gain is **not resolved**.')
    w('')
    w('## C2. W = 5 on T_φψ does not establish equivalence of A and B2')
    w('')
    w('The summary called A and B2 at W = 5 "statistically indistinguishable". The classes say only that both relative-loss '
      'intervals lie below ρ, i.e. both are within the declared acceptable-loss margin (Amendment 1 §4).')
    w('- That is neither an equivalence nor a non-inferiority claim.')
    w('- The intervals are marginal and unadjusted.')
    w('- The best method switched between A and B2 across resamples:')
    w('')
    w('| Cell | Observed best | Best counts over resamples | L(A) [95 %] | L(B2) [95 %] |')
    w('|---|---|---|---|---|')
    for k, e in w5.items():
        w(f"| {k} | {e['observed_best']} | {e['best_counts']} | {f(e['A']['L'])} {iv(e['A']['L_interval'])} | "
          f"{f(e['B2']['L'])} {iv(e['B2']['L_interval'])} |")
    w('')
    w('## C3. C1 failing the competitive criterion is not absence of gain over random')
    w('')
    w('The summary said that "every method except C1 is clearly better than random". That is wrong.')
    w('- C1\'s declared outcome *fails the criterion in all cells* refers to its relative loss against the best method '
      '(loss class *beyond the margin* everywhere).')
    w('- Its gain over random retrieval, G = 1 − median NR (coverage-matched), has a 95 % interval **above 0 in every '
      'cell**.')
    w('- On T_φψ its gain class is inconclusive or uninformative relative to ρ, not "no gain":')
    w('')
    w('| Cell | G(C1) | 95 % interval | Lower bound > 0 | Gain class vs ρ | Loss class |')
    w('|---|---|---|---|---|---|')
    for k, e in c1.items():
        w(f"| {k} | {f(e['G'])} | {iv(e['G_interval'])} | {e['G_lower_bound_above_0']} | {e['gain_class']} | "
          f"{e['loss_class']} |")
    w('')
    w('"Uninformative" in Amendment 1 means the upper bound of G is at or below ρ. That is a gain smaller than the '
      'declared threshold, not a gain of zero.')
    w('')
    w('## C4. Standing interpretation')
    w('')
    w('**Candidates for their respective geometric tasks:**')
    w('- **A** (residue-aligned κ–τ): Cα-only retrieval judged by an angular φ/ψ target. Declared T_φψ outcome: `'
      + rec['declared_overall']['phipsi']['A'] + '`.')
    w('- **B2** (window Cα RMSD): Cα-only retrieval judged by a Cartesian peptide-plane target. Declared T_pp outcome: `'
      + rec['declared_overall']['pp']['B2'] + '`.')
    w('')
    w('**Scope.** These are within-dataset generalization results:')
    w('- mdCATH, CHARMM22*, 320/348 K;')
    w('- 30 superfamilies;')
    w('- cross-replica retrieval from a shared starting configuration.')
    w('')
    w('**Not evidence of:**')
    w('- the legacy κ–τ Fréchet index or its writhe cascade (neither was part of this evaluation);')
    w('- any detector (bimodality scan, pocket, allostery);')
    w('- biological function;')
    w('- clinical or cancer-treatment efficacy.')
    w('')
    OUT_MD.write_text('\n'.join(L) + '\n')
    print('wrote', OUT_JSON.name, OUT_MD.name)
    print(json.dumps({'B2_W8': {k: (round(v['L'], 3), [round(x, 3) for x in v['L_interval']], v['loss_class'])
                                for k, v in b2.items()},
                      'C1_G_lower>0_all': all(e['G_lower_bound_above_0'] for e in c1.values())}))


if __name__ == '__main__':
    main()
