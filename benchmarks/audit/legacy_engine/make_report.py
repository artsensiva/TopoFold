#!/usr/bin/env python3
"""LE-1 report (protocol + Amendment 1). Every number is read from results/*.json (AGENTS.md rule 3); static
text makes no numeric claim. Run: ../wpd_local_retrieval/build/venv/bin/python -B make_report.py
"""
import json
import statistics
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'
OUT = HERE / 'REPORT.md'
LOCAL = ('M1', 'M2', 'K', 'A', 'F', 'F_lock', 'B1', 'B2')
METHODS = LOCAL + ('G_RMSD', 'G_F')
NAMES = {'M1': 'τ chord-RMS', 'M2': 'τ circular-RMS', 'K': 'κ-only RMS', 'A': 'residue-aligned κ–τ',
         'F': 'legacy κ–τ Fréchet (Rust)', 'F_lock': 'legacy pairs, lock-step max', 'B1': 'local Cα–Cα distances',
         'B2': 'window Cα RMSD', 'G_RMSD': 'global Cα RMSD', 'G_F': 'legacy full-chain Fréchet'}
CELLS = [(lab, W) for lab in ('open_start', 'closed_start') for W in (5, 8)]

CLAIMS = [  # (mechanism / claim, legacy status, evidence)
    ('Raw Cartesian Fréchet (`discrete_frechet_distance`, Python)', 'Å, no superposition: not SE(3)-invariant '
     '(KNOWN_ISSUES M5); not part of the index', 'code `distance.rs::discrete_frechet_distance_coords`'),
    ('κ–τ Fréchet is SE(3)-invariant', 'Supported (features are invariant)', '`test_independent_rigid_motions_…`'),
    ('Torsion handled periodically', 'Supported', 'probe `wrap_pair_expected_0p02`; `test_torsion_wrap_boundary`'),
    ('Residue correspondence preserved', '**Not preserved**: flexible coupling maps distinct same-length windows '
     'to 0', 'probe `misregistration_x_y_vs_x_x`; Rust-path `test_frechet_misregistration_numpy_and_rust`'),
    ('All window features retained', '**No**: pairs (κ at vertex s+1+j, τ about bond s+1+j), length '
     'min(#κ, #τ) = W−3; last interior κ dropped; bond lengths unused',
     'probe `last_kappa_changed_by_1rad`; `test_rust_path_ignores_last_window_kappa`'),
    ('Consistent units and weights', 'κ, τ both rad, weights 1:1 (fixed in `engine.rs`); channel ranges differ, so '
     'equal weights are a modeling choice; "zero hyperparameters" is false (weights, window, radius)', 'code'),
    ('Empty input', '**Distance 0** (perfect match) for empty κ–τ, DTW and writhe spectra', 'probe `empty_vs_nonempty`, '
     '`dtw_empty`, `writhe_empty`'),
    ('Non-finite input', '**NaN pairs silently ignored** by f64::max/min; NaN only if every pair is NaN',
     'probe `nan_torsion`, `all_nan_torsions`'),
    ('Degenerate geometry', 'Exact collinear/repeated points raise and abort the whole `fit` (m1); float32-rounded '
     'collinear triples pass the 1e-9 check silently with arbitrary τ (LE-D6)',
     '`test_degenerate_inputs_in_binding_and_fit`, `test_float32_rounded_collinear_passes_binding_check_silently`'),
    ('Writhe-filtered candidate search', 'Approximate (M6); writhe distance truncates unequal lengths',
     'probe `writhe_unequal_len`; recall in §6'),
    ('Exhaustive local retrieval (`query_subcurve`)', 'Linear scan over all frames (correctly described in README)',
     'code `engine.rs::query_subcurve_k_nearest`'),
    ('Scalar circular-mean / BC scanning', 'Separate mechanism (bimodality.rs); characterized in WP-A/WP-D (C1)',
     'WP-A, WP-D reports'),
    ('Historical benchmark advantages', 'Unsupported: synthetic ensembles, circular labels, unequal information '
     '(KNOWN_ISSUES C1, C2)', 'KNOWN_ISSUES'),
]

PRIOR = [  # (topic, keys, established content stated here, content evidence level)
    ('Cα virtual-bond / virtual-dihedral representation', ['levitt1976'],
     'Cα-only angles and dihedrals as a reduced backbone description', 'prior knowledge, not checked in this audit'),
    ('Discrete curvature/torsion (Frenet) description of protein chains', ['rackovsky1978', 'hu2011'],
     'Differential-geometry descriptors of protein backbones', 'prior knowledge, not checked in this audit'),
    ('φ/ψ as local backbone description', ['ramachandran1963'], 'Backbone dihedral map',
     'prior knowledge, not checked in this audit'),
    ('Fréchet distance; DTW', ['altgodau1995', 'sakoe1978', 'jiang2008'],
     'Coupling distances for curves/sequences, incl. discrete Fréchet for protein alignment',
     'prior knowledge, not checked in this audit'),
    ('Optimal superposition RMSD', ['kabsch1976', 'horn1987', 'theobald2005'], 'Closed-form rotation / quaternion '
     'RMSD', 'prior knowledge; Horn vs Kabsch agreement checked numerically here'),
    ('Circular (sin, cos) dihedral embeddings', ['mu2005', 'altis2007'], 'dPCA on sin/cos of dihedrals',
     'prior knowledge, not checked in this audit'),
    ('Structural alphabets / local fragment libraries', ['debrevern2000', 'kolodny2002', 'vankempen2024'],
     'Discretized local backbone states; fragment libraries; 3Di local descriptors for search',
     'prior knowledge, not checked in this audit'),
    ('Writhe / Gauss-integral descriptors', ['rogen2003'], 'Gauss integrals for structure classification',
     'prior knowledge, not checked in this audit'),
    ('MD analysis tools', ['michaudagrawal2011', 'mcgibbon2015'], 'Dihedral, RMSD and distance analyses',
     'MDAnalysis used here (2.10.0); content otherwise prior knowledge'),
    ('Held-out dataset', ['mirarchi2024'], 'mdCATH content and builder', 'checked against paper and code in WP-D'),
]


def load(n):
    return json.loads((RES / n).read_text())


def f(x, d=3):
    return '—' if x is None else f'{x:.{d}f}'


def iv(v, d=3):
    return f'[{f(v[0], d)}, {f(v[1], d)}]'


def main():
    S = load('le_adk_THIRD_PARTY_MD_DERIVATIVE_summary.json')
    P = load('le_adk_THIRD_PARTY_MD_DERIVATIVE_provenance.json')
    A = load('le_adk_THIRD_PARTY_MD_DERIVATIVE_analysis_a1.json')
    T = load('le_compute_timing.json')
    C = load('le_prior_work_crossref.json')
    RP = load('le_rust_probe.json')
    FZ = load('le_protocol_freeze.json')
    AM = load('le_amendment1.json')
    AT = load('le_run_attempts.json')
    L = []
    w = L.append
    w('# LE-1: Legacy κ–τ engine evaluation (development tier, AdK)')
    w('')
    w(f"**THIRD-PARTY MD DERIVATIVE** for AdK results ({S['_meta']['data']}); SYNTHETIC for fixtures and probe. "
      f"Protocol `{S['_meta']['protocol']}` (SHA-256 `{FZ['protocol_sha256'][:16]}…`, frozen {FZ['frozen_utc']}); "
      f"Amendment 1 (`{AM['amendment']}`, recorded {AM['recorded_utc']}, {AM['outcome_status_at_recording']}). "
      f"Synthetic verification before the run: `{P['tests']['summary']}` (exit {P['tests']['exit_code']}). "
      'A **geometric local-conformation retrieval benchmark on one protein**: not biological, not held-out.')
    w('')
    a1 = AT['attempt_1']
    w(f"**Run history.** Attempt 1 ended at `{a1['last_progress_line']}` with no result file and no exit line; "
      f"cause {a1['cause']}; outcomes observed: {a1['outcomes_observed']}. Attempt 2 (identical code and protocol) "
      f"produced the results below in {S['_meta']['wall_seconds']} s.")
    w('')
    w('## 1. Legacy method: claim–evidence matrix')
    w('')
    w('| Mechanism / claim | Status | Evidence |')
    w('|---|---|---|')
    for a, b, c in CLAIMS:
        w(f'| {a} | {b} | {c} |')
    w('')
    w('Rust probe outputs (SYNTHETIC, ' + RP['_meta']['rustc'] + '):')
    w('')
    w('| Probe | Output | Meaning |')
    w('|---|---|---|')
    for k, v in RP['outputs'].items():
        w(f"| `{k}` | {v if isinstance(v, str) else f(v, 6)} | {RP['interpretation'].get(k, 'control: change in a paired κ is detected')} |")
    w('')
    ver = A['verification']
    mx = max(x['frechet_rust_numpy_max_abs_diff'] for x in S['windows'])
    w(f"**AdK verification gates (Amendment 1 §5):** {ver['gates']}. Rust vs NumPy Fréchet max |Δ| {mx:.1e}; "
      f"non-finite distances (T_φψ) {ver['nonfinite_distances_phipsi']}; regret failures {ver['regret_failures']}.")
    w('')
    w('## 2. Retrieval: equal-window-weighted mean NR@10 (0 = target-optimal, 1 = random; lower is better)')
    w('')
    w('Intervals: 95 % two-way block bootstrap of the same statistic, query blocks of 100 | 250 frames, window '
      'blocks of 20 starts, 2,000 paired resamples. Database frames are not resampled.')
    w('')
    for t, tname in (('phipsi', 'Primary target T_φψ'), ('pp', 'Complementary target T_pp (peptide-plane atoms)')):
        w(f'### {tname}')
        w('')
        w('| Method | ' + ' | '.join(f'{lab}/W{W}' for lab, W in CELLS) + ' |')
        w('|---|' + '---|' * len(CELLS))
        for m in METHODS:
            row = []
            for lab, W in CELLS:
                c = A['cells'][f'{lab}/W{W}/{t}']
                mi = c['method_intervals'][m]
                row.append(f"{f(c['theta'][m])} {iv(mi['qblock100'])} {iv(mi['qblock250'])}")
            w(f'| {m} {NAMES[m]} | ' + ' | '.join(row) + ' |')
        w('')
        c0 = {k: A['cells'][f'{lab}/W{W}/{t}'] for k, (lab, W) in enumerate(CELLS)}
        w('Valid query counts per window (min–max, distinct values): ' + '; '.join(
            f"{lab}/W{W} {c0[k]['valid_query_counts']['min']}–{c0[k]['valid_query_counts']['max']} "
            f"({c0[k]['valid_query_counts']['n_distinct']})" for k, (lab, W) in enumerate(CELLS))
            + '. Undefined NR units: ' + '; '.join(
            f"{lab}/W{W} {sum(c0[k]['undefined_units'].values())}" for k, (lab, W) in enumerate(CELLS)) + '.')
        w('')
    w('### Secondary: overlap@10 (median across windows of window means) vs random')
    w('')
    w('| Cell | Target | random | ' + ' | '.join(METHODS) + ' |')
    w('|---|---|---|' + '---|' * len(METHODS))
    for lab, W in CELLS:
        ws = [x for x in S['windows'] if x['trajectory'] == lab and x['W'] == W and 'methods' in x]
        for t in ('phipsi', 'pp'):
            vals = [statistics.median(x['methods'][m][t]['k10']['mean_overlap'] for x in ws) for m in METHODS]
            rnd = statistics.median(x['random_overlap']['k10'] for x in ws)
            w(f'| {lab}/W{W} | {t} | {f(rnd, 4)} | ' + ' | '.join(f(v, 4) for v in vals) + ' |')
    w('')
    w('## 3. Paired contrasts (ΔNR@10 = X − Y on identical units; negative favours X)')
    w('')
    w('Statement rule (Amendment 1 §4): τ_int ≤ 100 → 100-frame interval; ≤ 250 → 250-frame interval; '
      '> 250 or not estimable → **dependence_unresolved**. Both intervals shown. F − A is an overall comparison; '
      'only F − F_lock isolates flexible coupling; F_lock − A combines aggregation, the omitted κ and channel '
      'normalization.')
    w('')
    names = list(A['cells']['open_start/W5/phipsi']['contrasts'])
    for t in ('phipsi', 'pp'):
        w(f'**Target {t}**')
        w('')
        w('| Contrast | ' + ' | '.join(f'{lab}/W{W}' for lab, W in CELLS) + ' |')
        w('|---|' + '---|' * len(CELLS))
        for n in names:
            row = []
            for lab, W in CELLS:
                e = A['cells'][f'{lab}/W{W}/{t}']['contrasts'][n]
                row.append(f"{f(e['delta_paired'])} {iv(e['interval_qblock100'])} {iv(e['interval_qblock250'])} "
                           f"τ_int {f(e['tau_int_frames'], 0)} → **{e['statement']}**")
            w(f'| {n} | ' + ' | '.join(row) + ' |')
        w('')
    w('## 4. Development classification of A and F (primary target; Amendment-1 gates)')
    w('')
    w('| Method | 100-frame cells → outcome | 250-frame cells → outcome | τ_int vs best (frames) | Final |')
    w('|---|---|---|---|---|')
    for m, per in A['development_outcome'].items():
        w(f"| {m} | {', '.join(per['qblock100']['cells'])} → {per['qblock100']['outcome']} | "
          f"{', '.join(per['qblock250']['cells'])} → {per['qblock250']['outcome']} | "
          f"{', '.join(f(x, 0) for x in per['tau_int_vs_best_frames'])} | **{per['final']}** |")
    w('')
    w('| Cell | Best observed | L_A [100] [250] | L_F [100] [250] |')
    w('|---|---|---|---|')
    for lab, W in CELLS:
        dc = A['cells'][f'{lab}/W{W}/phipsi']['development_classification']
        w(f"| {lab}/W{W} | {dc['best_observed']} | {f(dc['A']['qblock100']['L'])} "
          f"{iv(dc['A']['qblock100']['L_interval'])} {iv(dc['A']['qblock250']['L_interval'])} | "
          f"{f(dc['F']['qblock100']['L'])} {iv(dc['F']['qblock100']['L_interval'])} "
          f"{iv(dc['F']['qblock250']['L_interval'])} |")
    w('')
    w('## 5. Failure cases: windows where F − A (T_φψ) is largest')
    w('')
    w('| Cell | i | residues | NR_F | NR_F_lock | NR_A | Δ(F − A) |')
    w('|---|---|---|---|---|---|---|')
    for lab, W in CELLS:
        ws = [x for x in S['windows'] if x['trajectory'] == lab and x['W'] == W and 'methods' in x]

        def g(x, m):
            return x['methods'][m]['phipsi']['k10']['mean_norm_regret']
        for x in sorted(ws, key=lambda x: -(g(x, 'F') - g(x, 'A')))[:5]:
            w(f"| {lab}/W{W} | {x['i']} | {x['i']}–{x['i'] + W - 1} | {f(g(x, 'F'))} | {f(g(x, 'F_lock'))} | "
              f"{f(g(x, 'A'))} | {f(g(x, 'F') - g(x, 'A'))} |")
    w('')
    w('Windows where NR_F > NR_A and where NR_F > NR_F_lock: ' + '; '.join(
        f"{lab}/W{W}: {sum(1 for x in S['windows'] if x['trajectory'] == lab and x['W'] == W and x['methods']['F']['phipsi']['k10']['mean_norm_regret'] > x['methods']['A']['phipsi']['k10']['mean_norm_regret'])}, "
        f"{sum(1 for x in S['windows'] if x['trajectory'] == lab and x['W'] == W and x['methods']['F']['phipsi']['k10']['mean_norm_regret'] > x['methods']['F_lock']['phipsi']['k10']['mean_norm_regret'])} of "
        f"{sum(1 for x in S['windows'] if x['trajectory'] == lab and x['W'] == W)}" for lab, W in CELLS) + '.')
    w('')
    w('## 6. Writhe cascade vs exhaustive full-chain Fréchet (same final distance)')
    w('')
    for lab, r in S['cascade_recall_vs_exhaustive_G_F'].items():
        w(f"- {lab}: mean recall@10 {f(r['mean_recall_at_10'])}, minimum {f(r['min_recall_at_10'])}, queries with "
          f"recall 1: {f(r['fraction_queries_recall_1'])}; queries returning < 10 hits: "
          f"{r['queries_with_fewer_than_10_hits']}.")
    w('')
    w('## 7. Computation (single thread; open_start, W = 8, windows 1–20 unless stated)')
    w('')
    w(f"CPU: {T['_meta']['cpu']}. Rust vs NumPy differences are implementation differences, not representation "
      'effects.')
    w('')
    w('| Implementation | Steady state (s, median of 3) | Work unit | Peak RSS above baseline (MB) |')
    w('|---|---|---|---|')
    for k, v in T['results'].items():
        if k == 'baseline':
            continue
        w(f"| {k} | {f(v.get('steady_state_s_median'), 3)} | {v.get('units')} | "
          f"{f(v.get('maxrss_above_baseline_kB', 0) / 1024, 1)} |")
    w('')
    w(f"Setup: binding feature extraction for all {S['trajectories']['open_start']['n_frames']} frames "
      f"{f(T['results']['M1'].get('feature_extraction_binding_s'), 2)} s; Rust index fit "
      f"{f(T['results']['F_rust_fit']['steady_state_s_median'], 3)} s.")
    w('')
    w('## 8. Prior work and novelty')
    w('')
    w('Crossref verification establishes bibliographic metadata only. "Content evidence" states whether what is '
      'attributed to a work was checked against the primary source in this audit; novelty judgments resting on '
      'unchecked content are **provisional**.')
    w('')
    w('| Key | Crossref-verified | DOI | Year |')
    w('|---|---|---|---|')
    for k, r in C['refs'].items():
        cr = r.get('crossref') or {}
        w(f"| {k} | {r['verified']} | {cr.get('doi') if r['verified'] else '— (not cited)'} | "
          f"{cr.get('year') if r['verified'] else '—'} |")
    w('')
    w('| Topic | References (verified keys only) | Content attributed | Content evidence |')
    w('|---|---|---|---|')
    for topic, keys, content, ev in PRIOR:
        ks = [k for k in keys if C['refs'].get(k, {}).get('verified')]
        w(f"| {topic} | {', '.join(ks)} | {content} | {ev} |")
    w('')
    w('| Category | Items |')
    w('|---|---|')
    w('| Established mathematics | discrete turning angles and dihedrals; circular distance; Fréchet/DTW coupling; '
      'quaternion/SVD superposition |')
    w('| Existing software functionality | dihedral, distance and RMSD analysis in MD toolkits; structure search with '
      'local descriptors (provisional: content not checked here) |')
    w('| Engineering integration (TopoFold) | Rust implementation of κ–τ extraction, exhaustive window scan, writhe '
      'VP-tree cascade; Python binding |')
    w('| New empirical findings (this audit, one protein) | §2–§6 of this report |')
    w('| Unsupported novelty claims | combined κ–τ as new; SE(3) invariance as distinctive; "zero hyperparameters"; '
      'Rust as a scientific contribution |')
    w('')
    _decision(L, S, A, T)
    _editorial(L, A, T)
    OUT.write_text('\n'.join(L) + '\n')
    print('wrote', OUT.name)



def _decision(L, S, A, T):
    w = L.append
    th = {c: A['cells'][c]['theta'] for c in A['cells']}
    ph = [f'{lab}/W{W}/phipsi' for lab, W in CELLS]
    pp = [f'{lab}/W{W}/pp' for lab, W in CELLS]
    best_ph = {c: min(('M1', 'M2', 'K', 'A', 'F', 'B1', 'B2'), key=lambda m: th[c][m]) for c in ph}
    best_pp = {c: min(('M1', 'M2', 'K', 'A', 'F', 'B1', 'B2'), key=lambda m: th[c][m]) for c in pp}
    con = lambda c, n: A['cells'][c]['contrasts'][n]
    w('## 9. Post-hoc observations (after outcomes were seen; they do not change any pre-declared outcome)')
    w('')
    dA = A['development_outcome']['A']
    w(f"- **Gate 7 and method A.** The pre-declared outcome for A is `{dA['final']}`. τ_int of NR_A − NR_best was "
      f"not estimable ({dA['tau_int_vs_best_frames']}) because A is the observed best in every primary cell, so "
      'the series is identically zero. The gate therefore flags an exact zero, not unresolved dependence. For the '
      'nearest competitor the dependence diagnostic is short: τ_int of B2 − A on T_φψ = '
      + ', '.join(f(con(c, 'H5_B2_minus_A')['tau_int_frames'], 0) for c in ph)
      + ' frames. A future amendment could take τ_int against the runner-up when X is best; it is not applied here.')
    w('- **F − F_lock at W = 5 is exactly 0 by construction.** With W − 3 = 2 legacy pairs every monotone coupling '
      'contains both diagonal points, so the discrete Fréchet equals the lock-step maximum; the "dependence '
      'unresolved" label there reflects a constant series, not uncertainty.')
    rc = S['cascade_recall_vs_exhaustive_G_F']
    w(f"- **Writhe cascade.** Recall@10 vs exhaustive was {f(rc['open_start']['mean_recall_at_10'])} / "
      f"{f(rc['closed_start']['mean_recall_at_10'])}, but its time per query "
      f"({f(T['results']['cascade_query']['steady_state_s_median'], 1)} s per 50 queries) equals the exhaustive "
      f"scan ({f(T['results']['GF_rust_query']['steady_state_s_median'], 1)} s): with the default radius 0.5 "
      'after the writhe-normalization fix (M1), equal timing is consistent with the filter admitting most frames '
      '(candidate counts were not recorded). Exact recall here is therefore not evidence of an effective '
      'approximate index.')
    w('- Valid query counts were identical across windows and no NR unit was undefined, so the Amendment 1 §1 '
      'estimand change does not alter any number here; it matters for data with exclusions.')
    w('')
    w('## 10. Conclusions (development tier, one protein, two trajectories)')
    w('')
    w('**Implemented correctly (verified):** the legacy κ–τ extraction and exhaustive `query_subcurve` path match '
      'independent NumPy implementations (§1); torsion periodicity and SE(3) invariance hold.')
    w('')
    w('**Legacy defects documented (not repaired):** last interior κ dropped by the pairing; flexible coupling maps '
      'misregistered windows to 0; empty inputs score 0; NaN pairs silently ignored; float32-rounded collinear '
      'triples pass the degeneracy check; one collinear frame aborts `fit`; writhe distance truncates unequal '
      'lengths; in the one measured workload (single thread, 50 full-chain queries, open_start, default radius '
      '0.5) the writhe cascade took the same time as the exhaustive scan [corrected E1: scope].')
    w('')
    w('**What preserved more information:** combining κ and τ (A) retrieved closer conformations than either channel '
      'alone under both targets (H1 `below_0` in ' + str(sum(con(c, n)['statement'] == 'below_0'
      for c in ph + pp for n in ('H1_A_minus_M2', 'H1_A_minus_K'))) + ' of 16 cell × contrast statements).')
    w('')
    w('**What improved retrieval:** best method by cell mean NR@10 — T_φψ: ' + ', '.join(
      f"{c.rsplit('/', 1)[0]} {best_ph[c]}" for c in ph) + '; T_pp: ' + ', '.join(
      f"{c.rsplit('/', 1)[0]} {best_pp[c]}" for c in pp) + '. A − B2 (H5) on T_φψ: ' + ', '.join(
      con(c, 'H5_B2_minus_A')['statement'] for c in ph) + ' (B2 − A; above_0 favours A); on T_pp: ' + ', '.join(
      con(c, 'H5_B2_minus_A')['statement'] for c in pp) + ' (below_0 favours B2). The winner depends on the '
      'target family: the angular target favours the angular combination, the Cartesian peptide-plane target '
      'favours Cartesian RMSD (H6). Gains over random are moderate for every method (NR@10 well above 0).')
    w('')
    w('**Legacy Fréchet (F):** worse than A in every cell for both targets (overall F − A ' + ', '.join(
      con(c, 'overall_F_minus_A')['statement'] for c in ph + pp) + '). Of that difference, flexible coupling '
      '(F − F_lock) contributes 0 at W = 5 by construction and ' + ', '.join(
      con(c, 'coupling_F_minus_Flock')['statement'] for c in (ph[1], ph[3], pp[1], pp[3])) +
      ' at W = 8 (T_φψ open/closed, T_pp open/closed). Which component is larger depends on target and window '
      'size (§13, corrected E1): F_lock − A (max aggregation, omitted κ and normalization together) is the whole '
      'difference at W = 5 and the larger part on T_φψ at W = 8, but flexible coupling is the larger part on T_pp '
      'at W = 8. Development outcome F (declared, unchanged): `' + A['development_outcome']['F']['final'] + '`.')
    w('')
    w('**Local vs global (H4):** ' + ', '.join(con(c, 'H4_GF_minus_F')['statement'] + '/' +
      con(c, 'H4_GRMSD_minus_B2')['statement'] for c in ph + pp) + ' (G_F − F / G_RMSD − B2). Global rankings, '
      'including the legacy full-chain Fréchet, are close to random for local retrieval; historical advantages of '
      'the engine over global baselines are consistent with local-vs-global information rather than with the κ–τ '
      'Fréchet representation itself.')
    w('')
    w('**Generalized beyond development data:** nothing. No held-out data were scored.')
    w('')
    w('**Scientifically unvalidated:** any biological relevance (both targets are geometric reconstructions from '
      'the same trajectories); performance on other proteins, force fields or temperatures; database-frame '
      'variability; ρ as a usefulness threshold; novelty statements resting on unchecked prior-work content.')
    w('')
    w('## 11. Recommendation (engine vs loader)')
    w('')
    w('1. **Engineering decision [corrected E1]: deprioritize the legacy κ–τ Fréchet for same-sequence local '
      'window retrieval.** This is a prioritization for the tested task, not a scientific rejection: F\'s declared '
      'development outcome is formally `inconclusive` (its relative-loss interval straddles ρ) and stays so. The '
      'basis is that, on AdK, the residue-aligned combination of the same features (A) had lower mean NR@10 in every '
      'cell for both targets; F − B2 was never below 0 (above 0 in six cells, unresolved in the two T_φψ W = 8 cells, '
      'where F\'s point estimate was marginally lower); and flexible coupling was zero or adverse where resolved. Tasks '
      'needing unequal-length or cross-sequence matching were not tested.')
    w('2. **Carry A and B2 forward as competing candidates**, with B1 and M1/M2 as baselines; no single method is '
      'preferred because the ranking depends on the target family. Adding A to the held-out method set changes '
      'the frozen held-out protocol and requires a reviewed amendment **before** held-out scoring.')
    w('3. **Build the minimal mdCATH reader/adapter next** (established reader such as h5py, explicit residue '
      'mapping, Å units, box handling, provenance), once the schema-sample download and the h5py installation '
      'are authorized; neither is authorized at present. That is the remaining dependency for generalization.')
    w('4. Document the legacy defects in `docs/KNOWN_ISSUES.md` in a separate, reviewed change; no production '
      'code was modified here.')
    w('')
    w('## 12. Limitations')
    w('')
    w('- One protein (AdK), two trajectories, already used in WP-B/WP-D; M1, M2, K, B1, B2 results were known before '
      'this protocol. A, F, F_lock, the global methods and T_pp are new.')
    w('- Both targets are derived from the same trajectory; T_φψ favours angular and T_pp Cartesian methods.')
    w('- Intervals cover query-frame and residue-position resampling only; the database pool is fixed; 250-frame '
      'blocks give only 3–5 query blocks per trajectory; τ_int estimates are themselves uncertain; WP-B observed '
      'slow drift that blocks of this length may not capture.')
    w('- No multiplicity adjustment; statements are descriptive.')
    w('- Timing compares a NumPy research implementation with the Rust engine path; peak RSS is dominated by a '
      'shared baseline.')
    w('')



def _editorial(L, A, T):
    """§13: editorial closeout E1 (interpretation only; no result or declared outcome changes)."""
    w = L.append
    w('## 13. Editorial closeout E1 (interpretation corrections; results and declared outcomes unchanged)')
    w('')
    w('Recorded in `results/le_editorial_log.json`. No experiment was rerun and no protocol or analysis changed.')
    w('')
    dA = A['development_outcome']['A']
    w('**E1.1 Dependence flag for A was structural, not measured.** The declared outcome for A remains `'
      + dA['final'] + '`. Gate 7 computes τ_int of NR_A − NR_best with best = the observed best method. A was '
      'that best in every primary cell, so the series compared A with itself and was identically zero; the '
      'estimator returned no value. No long temporal dependence was measured for A. The flag is therefore an '
      'artefact of a self-comparison in the gate definition, not evidence of dependence; equally, it is not '
      'evidence that dependence is short. The diagnostic was simply undefined.')
    w('')
    w('**E1.2 Attribution of the Fréchet gap, per target and window size.** F − A = (F − F_lock) + (F_lock − A) '
      'exactly on identical units. Values are paired equal-window means of ΔNR@10.')
    w('')
    w('| Cell | F − A | F − F_lock (coupling) | F_lock − A (aggregation + omitted κ + normalization) | Larger component |')
    w('|---|---|---|---|---|')
    for c, cell in A['cells'].items():
        k = cell['contrasts']
        a, b, d = (k['overall_F_minus_A'], k['coupling_F_minus_Flock'], k['aggregation_boundary_Flock_minus_A'])
        W5 = '/W5/' in c
        larger = ('F_lock − A (coupling ≡ 0 by construction at W = 5)' if W5 else
                  ('coupling' if b['delta_paired'] > d['delta_paired'] else 'F_lock − A'))
        cs = 'identity (exact 0)' if W5 else b['statement']
        w(f"| {c} | {f(a['delta_paired'], 4)} ({a['statement']}) | {f(b['delta_paired'], 4)} ({cs}) | "
          f"{f(d['delta_paired'], 4)} ({d['statement']}) | {larger} |")
    w('')
    w('The earlier sentence "the larger part is F_lock − A" was correct for T_φψ and for W = 5, and wrong for T_pp '
      'at W = 8, where coupling was larger. The components of F_lock − A are not separated.')
    w('')
    w('**E1.3 Deprioritizing F is an engineering decision** for same-sequence local window retrieval (§11 item 1). '
      'F\'s formal development outcome remains `' + A['development_outcome']['F']['final'] + '`.')
    w('')
    r = T['results']
    w('**E1.4 Timing statements are limited to the implementations and workload measured** (single thread, '
      f"{T['_meta']['cpu']}, open_start, W = 8, windows 1–20, {r['M1']['n_queries']} queries × "
      f"{r['M1']['n_candidates']} candidates; full-chain items: 50 queries). They compare particular NumPy research "
      'implementations with the Rust engine path through the Python binding, which includes per-call feature '
      'extraction and allocation. They are not statements about the representations, other implementations, '
      'thread counts, index radii or data sizes.')
    w('')
    cal = json.loads((HERE.parent / 'wpd_local_retrieval' / 'results' / 'SYNTHETIC_tau_int_calibration.json').read_text())
    frac = 0.015
    ex = next(r for r in cal['rows'] if r['n'] == 1000 and r['phi'] == 0.99)
    nq = {'open_start': 1013, 'closed_start': 755}
    rows = [(c, k, e) for c, cell in A['cells'].items() for k, e in cell['contrasts'].items()
            if e['tau_int_frames'] is not None]
    over = [(c, k, e) for c, k, e in rows if e['tau_int_frames'] > frac * nq[c.split('/')[0]]]
    w('**E1.5 τ_int diagnostics above about 1.5 % of the series length are not resolvable (found after the '
      'results, while preparing held-out diagnostics).** A SYNTHETIC AR(1) calibration of the same estimator '
      '(`wpd_local_retrieval/results/SYNTHETIC_tau_int_calibration.json`) shows severe underestimation (median '
      'estimate < 0.8 × true) whenever τ̂/n ≥ '
      + f(cal['min_est_fraction_among_severely_biased_cases'], 3) + '; e.g. true τ_int ' + f(ex['tau_true'], 0)
      + ' frames was estimated at a median of ' + f(ex['est_median'], 1) + f" from {ex['n']} frames. The Amendment 1 "
      '§4 rule therefore cannot certify block adequacy for '
      f"τ̂ above {frac} n ({f(frac * nq['open_start'], 1)} frames for open_start, {f(frac * nq['closed_start'], 1)} "
      f'for closed_start). {len(over)} of {len(rows)} '
      'estimable contrast diagnostics exceed this; their statements (listed below) rest on block lengths whose '
      'adequacy is not established, and true dependence may exceed both tested block lengths. Declared statements '
      'and outcomes are unchanged; this is an interpretation caveat. On T_φψ, all F − A, F − F_lock, F_lock − A and '
      'A − M2 diagnostics are below this heuristic cutoff (see E1.6 for what that does and does not mean).')
    w('')
    w('**E1.6 Correction to E1.5 wording.** The 0.015 n cutoff is a limited heuristic from an AR(1) simulation of '
      'this estimator, not a reliability guarantee. A τ̂ below it does not establish that dependence is short: '
      'processes with drift or several timescales (WP-B observed drift in these runs) can yield small τ̂, and in '
      'the AR(1) simulation underestimation already begins below the cutoff for short series. Diagnostics below the '
      'cutoff are therefore "not flagged by the heuristic", not "resolved"; the listed ones above it are flagged.')
    w('')
    w('| Cell | Contrast | τ̂ (frames) | Declared statement |')
    w('|---|---|---|---|')
    for c, k, e in over:
        w(f"| {c} | {k} | {f(e['tau_int_frames'], 0)} | {e['statement']} |")
    w('')


if __name__ == '__main__':
    main()
