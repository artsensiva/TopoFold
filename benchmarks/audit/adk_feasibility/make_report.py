#!/usr/bin/env python3
"""Render REPORT.md from the WP-B G0/G1 results (Amendments 1–2). Every number is read from JSON."""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'


def f(x, p=2):
    return 'n/a' if x is None else f'{x:.{p}f}'


def main():
    fr0 = json.loads((RES / 'adk_spec_freeze.json').read_text())
    fr1 = json.loads((RES / 'adk_spec_amendment1_freeze.json').read_text())
    fr2 = json.loads((RES / 'adk_spec_amendment2_freeze.json').read_text())
    g0 = json.loads((RES / 'adk_g0_provenance.json').read_text())
    g1p = RES / 'adk_g1_summary.json'
    g1 = json.loads(g1p.read_text()) if g1p.exists() else None
    L = []
    w = L.append
    w('# WP-B G0–G1 Feasibility Report: AdK (third-party MD, not SYNTHETIC)\n')
    w('Data: Orädd et al., *Sci. Adv.* 2021, Zenodo 10.5281/zenodo.5583119 (CC BY 4.0). References: '
      'RCSB 4AKE (open) and 1AKE (closed). The report covers data integrity (G0) and an independent, '
      'reference-based assessment (G1). **No TopoFold or candidate-feature scoring, retrieval, '
      'localization or threshold tuning was performed.** AdK is a development dataset, not held-out '
      'validation.\n')
    w('**Specification history.** Original freeze ' + fr0['frozen_utc'] + '; Amendment 1 ' +
      fr1['frozen_utc'] + '; Amendment 2 ' + fr2['frozen_utc'] + '. All three were recorded before '
      'G0 or G1 ran and before any trajectory-derived structural quantity was computed. Amendment 1 '
      'answered review findings; Amendment 2 replaced a core algorithm that failed on a synthetic '
      'test (see `AMENDMENTS.md`). The original freeze record is preserved unchanged.\n')

    w('## G0: integrity and mapping\n')
    w(f"Result: **{'PASS' if g0['g0_passed'] else 'STOP'}** (finished {g0['_meta']['finished_utc']}). "
      f"Stops: {g0['stops'] or 'none'}.\n")
    w('| File | Size (B) | md5 and size match the source record | SHA-256 |\n|---|---:|---|---|')
    for n, x in g0['files'].items():
        w(f"| {n} | {x['size']} | {x['md5_matches_source_record'] and x['size_matches_source_record']} | `{x['sha256']}` |")
    for rid, r in g0['references'].items():
        w(f"\n**{rid}** (SHA-256 `{r['sha256']}`; REVDAT dates {r['revdat_dates'][:2]}; hetero residues "
          f"{r['hetero_residues']}):")
        for c, ch in r['chains'].items():
            w(f"- chain {c} ({ch['role']}): mapping valid **{ch['valid']}**; missing {ch['missing'] or 'none'}; "
              f"mismatched {ch['mismatched'] or 'none'}; residues with altlocs {len(ch['altlocs'])}"
              + (f" (used: {sorted({v['used'] for v in ch['altlocs'].values()})})" if ch['altlocs'] else ''))
    w('\n| Trajectory | atoms gro / xtc | Cα | frames | nominal Δt (ps) | cadence irregularities | first–last (ps) | file duration (ns) | stated (ns) | box angle dev (°) | quantization (Å) |')
    w('|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|')
    for lab, t in g0.get('trajectories', {}).items():
        w(f"| {lab} | {t['gro_n_atoms']} / {t['xtc_n_atoms']} | {t['n_ca']} | {t['n_frames']} | "
          f"{t['nominal_dt_ps']} | {t['cadence_irregularities']} | {t['time_first_ps']:.0f}–{t['time_last_ps']:.0f} | "
          f"{t['file_duration_ns']:.1f} | {t['source_statement_ns']} | {t['box_angles_max_dev_deg']:.1e} | "
          f"{t['coordinate_quantization_A']:.3f} |")
    w('\nVerified facts are those measured from the files above. The stated durations come from the '
      'paper and are compared with the files, not assumed. Discrepancies and notes: '
      + ('; '.join(g0['discrepancies_and_notes']) or 'none') + '. Coordinates are stored at the '
      'quantization shown; float64 arithmetic does not recover precision lost in the file.\n')
    if g1 is None:
        w('G1 not run.')
        (HERE / 'REPORT.md').write_text('\n'.join(L) + '\n')
        return 0
    if 'stopped' in g1:
        w(f"## G1 stopped\n\n{g1['stopped']}\n")
        (HERE / 'REPORT.md').write_text('\n'.join(L) + '\n')
        return 0

    w('## G1: independent reference-based assessment\n')
    b = g1['g0_binding']
    w(f"**Binding to G0:** {'all checks passed' if all(v for k, v in b.items() if k != 'g0_provenance_sha256') else b} "
      f"(spec, raw-file and reference hashes match the G0 record; G0 ran after the Amendment-2 freeze).\n")
    w('**Alignment cores** (crystal Cα only; two-phase algorithm of Amendment 2):\n')
    w('| Core | size | converged | valid | reason if invalid |\n|---|---:|---|---|---|')
    for k, c in g1['cores'].items():
        w(f"| {k} | {c.get('size', 'n/a')} | {c.get('converged', 'n/a')} | {c['valid']} | {c.get('invalid_reason') or ''} |")
    w(f"\nVariants used: {g1['variants_used']}. Excluded (invalid cores): {g1['variants_excluded'] or 'none'}.")
    w(f"d_ref (Å): {', '.join(f'{k} {v:.2f}' for k, v in g1['d_ref_A'].items())}. Non-core segments "
      f"(≥ 5 residues): {g1['segments_resids']}.")
    xc = g1['kabsch_crosscheck']
    w(f"\n**Kabsch vs MDAnalysis (QCP):** {xc['n_comparisons']} comparisons; max |diff| "
      f"{xc['max_abs_diff']:.1e} Å; predeclared tolerance {xc['tolerance_A']:.0e} Å → "
      f"**{'PASS' if xc['passed'] else 'FAIL'}**.\n")
    for lab, t in g1['trajectories'].items():
        w(f"### {lab} ({t['n_frames']} frames)\n")
        w(f"Preprocessing: frames with the raw Cα chain split across the periodic box: "
          f"{t['frames_chain_split_in_raw']}. Sequential minimum-image reconstruction was applied to "
          f"every frame; max consecutive Cα distance afterwards {t['max_consecutive_ca_after_whole_A']:.2f} Å. "
          'Raw files are unchanged.\n')
        w('| Observable (Å, primary core) | min | q05 | median | q95 | max | first 10 % mean | last 10 % mean | τ_int (frames, est.) | ESS (est.) |')
        w('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
        for k, x in t['observables'].items():
            w(f"| {k} | {f(x['min'])} | {f(x['q05'])} | {f(x['median'])} | {f(x['q95'])} | {f(x['max'])} | "
              f"{f(x['first10pct_mean'])} | {f(x['last10pct_mean'])} | {f(x['tau_int_frames_est'], 1)} | "
              f"{f(x['ess_est'], 1)} |")
        w('')
    import csv
    w('### Non-core segment RMSDs (pre-registered observable; from the per-frame CSVs)\n')
    w('Mean segment Cα RMSD (Å) after primary-core superposition, first vs last 10 % of frames.\n')
    for lab in g1['trajectories']:
        rows = list(csv.DictReader(open(RES / f"adk_g1_frames_{lab.split('_')[0]}.csv")))
        m = max(1, len(rows) // 10)
        segcols = [c for c in rows[0] if c.startswith('seg') and c.endswith('_rmsd_open')]
        w(f'**{lab}**\n')
        w('| Segment | to open: first → last | to closed: first → last |\n|---|---|---|')
        for c in segcols:
            cc = c.replace('_rmsd_open', '_rmsd_closed')
            mean = lambda col, part: sum(float(r[col]) for r in part) / len(part)
            w(f"| {c[3:].replace('_rmsd_open', '').replace('_', '–')} | {mean(c, rows[:m]):.2f} → {mean(c, rows[-m:]):.2f} | "
              f"{mean(cc, rows[:m]):.2f} → {mean(cc, rows[-m:]):.2f} |")
        w('')
    w('### Descriptive proximity classes by variant\n')
    w('A sustained run has an elapsed span ≥ 5.0 ns. Each cell gives frames / sustained runs / brief '
      'runs / longest span (ns). The classes are provisional and descriptive, not validated state '
      'definitions.\n')
    w('| Variant (radius Å) | category | trajectory | starts near own ref | near open | near closed | ambiguous |')
    w('|---|---|---|---|---|---|---|')
    for v, x in g1['classification'].items():
        for lab, t in x['trajectories'].items():
            s = t['summary']
            cell = lambda c: f"{s[c]['frames']} / {s[c]['sustained_runs']} / {s[c]['brief_runs']} / {s[c]['longest_span_ns']:.1f}"
            w(f"| {v} ({x['radius_A']:.2f}) | {x['category']} | {lab} | {t['starts_near_own_reference']} | "
              f"{cell('open')} | {cell('closed')} | {cell('ambiguous')} |")
    w('\n**Observed cross visits vs claimed transition evidence:**')
    for lab, e in g1['cross_visits_and_transition_evidence'].items():
        n_obs = sum(e['observed_cross_visit_by_variant'].values())
        w(f"- {lab}: sustained visit near the other reference observed in {n_obs} of "
          f"{len(e['observed_cross_visit_by_variant'])} variants; transition evidence claimed "
          f"(all variants, starting near own reference): **{e['transition_evidence_claimed']}**.")
    o = g1['outcome']
    w('\n## Outcome\n')
    w(f"- Categories by variant: {o['categories_by_variant']}; variants agree: {o['variants_agree']}.")
    w(f"- **Outcome: {o['final']}.**\n")
    w('**Evidence that would resolve the outcome** (any revision needs a new freeze before re-analysis; '
      'the current criteria are not changed after inspection):')
    w('1. State criteria fixed independently of these trajectories, e.g. segment-level criteria derived '
      'from the crystal structures, or verified literature domain observables.')
    w('2. An independent scale for thermal fluctuation around each reference state, to justify a '
      'proximity radius; the current 0.25·d_ref is a pre-declared choice, not a calibrated one.')
    w('3. Replicate trajectories per starting state.')
    w('4. Clarification of the closed-start system: it was simulated without the Ap5A ligand '
      'present in 1AKE, so its stability near the closed reference is not expected a priori.\n')
    w('**Interpretation limits.**')
    w('- Proximity to a reference by RMSD is not proof of a conformational transition.')
    w('- With one trajectory per starting state and autocorrelated frames, no equilibrium '
      'populations, rates or free energies can be inferred. The τ_int and ESS values are exploratory '
      'diagnostics, not evidence of equilibrium sampling; the observables drift over the runs.')
    w('- Segment observations are descriptive changes in pre-registered observables, not claimed '
      'transitions or assigned conformational states.')
    w('- The synthetic core-algorithm test (Amendment 2) is evidence about that construction only, '
      'not proof about AdK. On AdK the primary core converged with the size shown above.')
    w('- Pooling the two trajectories mixes starting states by construction and cannot establish '
      'spontaneous kinetics.')
    w('- If no transition occurred, that is a limitation of this dataset for the proposed task, not '
      'a failure of any representation; none was evaluated here.')
    w('- G2 label definitions and the held-out selection remain prerequisites for any scoring '
      'protocol.')
    (HERE / 'REPORT.md').write_text('\n'.join(L) + '\n')
    print('REPORT.md written')
    return 0


if __name__ == '__main__':
    sys.exit(main())
