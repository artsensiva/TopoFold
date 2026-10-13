#!/usr/bin/env python3
"""Render REPORT.md from the WP-C results (SYNTHETIC). Every number is read from JSON."""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = json.loads((HERE / 'results' / 'SYNTHETIC_ordered_representation_results.json').read_text())
P = json.loads((HERE / 'results' / 'SYNTHETIC_ordered_representation_provenance.json').read_text())
REPS = ['D_ch', 'D_circ', 'scalar_comparator', 'ca_distance', 'D_kappa']


def f(x, p=3):
    return 'null' if x is None else f'{x:.{p}f}'


def main():
    L = []
    w = L.append
    c, v = R['checks'], R['verdicts']
    w('# WP-C Ordered Angular Representation Report (SYNTHETIC)\n')
    w('Evaluation of a candidate representation (residue-ordered (cos τ, sin τ) features with '
      'same-sequence distances) on SYNTHETIC fixtures. It covers representation fidelity and '
      'fixed-fixture retrieval only. **No mode or state detection is performed or claimed, no '
      'production behaviour is changed, and no real-MD data are used.** Protocol: '
      '`docs/ORDERED_ANGULAR_REPRESENTATION_PROTOCOL.md` (amendments §11, recorded before '
      'execution). Every number below is read from `results/SYNTHETIC_ordered_representation_*.json`.\n')

    w('## Provenance\n')
    b = R['inputs']['binding']
    w(f"- Git HEAD `{P['git']['head']}` (working tree with {len(P['git']['status_porcelain'])} "
      f"uncommitted entries at run time). Python {P['versions']['python']}, NumPy "
      f"{P['versions']['numpy']}, {P['versions']['rustc']}.")
    w(f"- Production binding (T1): `{b['extension_path']}`, SHA-256 `{b['extension_sha256']}`. "
      f"Identical to the WP-A freshly built extension: **{b['equals_wpa_freshly_built_extension']}**. "
      f"Production sources unchanged since that build: "
      f"**{b['production_sources_unchanged_since_wpa_build']}** (compared: "
      f"{', '.join('`' + s + '`' for s in b['compared_sources'])}).")
    w(f"- WP-A inputs read from: {R['inputs']['wpa_input_source']}. Each was verified against the "
      f"frozen WP-A manifest (SHA-256 `{R['inputs']['wpa_manifest_sha256']}`); nothing was written "
      'into the WP-A directory, and the WP-A pipeline was not rerun.')
    w(f"- Historical scalar comparator: the unchanged WP-A Rust example (production functions), "
      f"writing into the WP-C `build/`. Outputs: JSON `{R['inputs']['comparator_output_sha256']['json']}`, "
      f"arrays `{R['inputs']['comparator_output_sha256']['arrays']}`.")
    w('- Commands and exit codes:\n')
    w('| Exit | Command |\n|---:|---|')
    for x in P['commands']:
        w(f"| {x['exit']} | `{x['cmd'][:140]}` |")
    w('\nPaths per case:')
    w('- Every case: T0 (prescribed), T2 (NumPy recovery from f64 coordinates) and T1 (production '
      'binding, float32 input).')
    w('- Retrieval fixtures: additionally the scalar comparator and the Cα-distance baseline.')
    w('- Case lists and families are in `results` → `inputs.case_families`.\n')

    w('## 1. Harness integrity (V1)\n')
    w('| Component | Result |\n|---|---|')
    for k, ok in v['V1_harness']['components'].items():
        w(f"| {k} | {'PASS' if ok else 'FAIL'} |")
    w(f"\nRecovery: max T2 error {c['recovery']['max_T2_err']:.2e} rad; max T1 (float32) error "
      f"{c['recovery']['max_T1_err']:.2e} rad.\n")

    w('## 2. Fidelity and continuity (V2)\n')
    w('**Feature-contract checks** (§5a; distances alone cannot detect these conventions):\n')
    w('| Check | Result |\n|---|---|')
    for x in c['feature_contracts']:
        w(f"| {x['name']} | {'PASS' if x['passed'] else 'FAIL'} |")
    w('\n**Distance-geometry checks** (§5b):\n')
    w('| Check | Measured | Result |\n|---|---|---|')
    w(f"| AD1 chord identity | max error {c['AD1_chord_identity']['max_abs_err']:.1e} | "
      f"{'PASS' if c['AD1_chord_identity']['passed'] else 'FAIL'} |")
    w(f"| AD2 (2/π)·D_circ ≤ D_ch ≤ D_circ on every retrieval matrix | failures "
      f"{len(c['AD2_bounds_on_retrieval_matrices']['failures'])} | "
      f"{'PASS' if c['AD2_bounds_on_retrieval_matrices']['passed'] else 'FAIL'} |")
    for cid, x in c['AD3_F1_analytic']['cases'].items():
        w(f"| AD3 {cid} A vs B | D_ch {x['D_ch']:.15f}, D_circ {x['D_circ']:.15f}, scalar {x['scalar']:.3e} | "
          f"{'PASS' if x['passed'] else 'FAIL'} |")
    w(f"| AD4 paired-shift invariance (F2) | max diff {c['AD4_paired_shift']['max_abs_diff']:.1e} | "
      f"{'PASS' if c['AD4_paired_shift']['passed'] else 'FAIL'} |")
    w(f"| AD5 Lipschitz continuity (F3) | max ‖Δx‖/δ {c['AD5_continuity']['max_feature_change_per_unit_delta']:.9f} | "
      f"{'PASS' if c['AD5_continuity']['lipschitz_passed'] else 'FAIL'} |")
    w(f"| AD6 mirror (N2) | D_ch error {c['AD6_mirror']['max_abs_err_D_ch']:.1e}; Cα-distance difference "
      f"{c['AD6_mirror']['max_ca_distance_diff_A']:.1e} Å | {'PASS' if c['AD6_mirror']['passed'] else 'FAIL'} |")
    w(f"| AD7 branch | D_ch {c['AD7_branch']['D_ch']:.2e} | {'PASS' if c['AD7_branch']['passed'] else 'FAIL'} |")
    w(f"| AD8 F6b equal torsions | max error {c['AD8_F6b']['max_abs_err']:.1e} | "
      f"{'PASS' if c['AD8_F6b']['passed'] else 'FAIL'} |")
    w(f"| Recovery-aware distances (T2, T1) | failures {len(c['distance_recovery_aware']['failures'])} | "
      f"{'PASS' if c['distance_recovery_aware']['passed'] else 'FAIL'} |")
    w('\n**Continuity across ±ε** (F3 noiseless, frame 0). The candidate varies continuously '
      'where the historical scalar summary jumps:\n')
    w('| ε | D_ch(+ε, −ε) | expected √2·sin ε | scalar comparator (+ε vs −ε) |\n|---|---:|---:|---:|')
    for e, x in c['AD5_continuity']['plus_minus_eps'].items():
        w(f"| {e} | {x['D_ch_plus_vs_minus']:.6e} | {x['expected']:.6e} | {x['scalar_comparator_plus_vs_minus']:.6f} |")

    w('\n## 3. Fixed-fixture retrieval (V3, descriptive)\n')
    w('Expected neighbour purity at k = 10 on T0, with label-independent tie handling. "ties" is '
      'the fraction of queries whose k-th boundary was tied. These are scores on fixed fixtures, '
      'with no confidence intervals and no population inference.\n')
    w('| Fixture | chance | D_ch | D_circ | scalar | Cα-dist | D_κ | ties D_ch / scalar / Cα / κ | D_ch (T1) |')
    w('|---|---:|---:|---:|---:|---:|---:|---|---:|')

    def row(name, r):
        t = r['T0']
        w(f"| {name} | {f(t['D_ch']['k10']['chance'])} | "
          + ' | '.join(f(t[x]['k10']['purity']) for x in REPS) + ' | '
          + ' / '.join(f(t[x]['k10']['tie_frequency'], 2) for x in ['D_ch', 'scalar_comparator',
                                                                    'ca_distance', 'D_kappa'])
          + f" | {f(r['T1']['D_ch']['k10']['purity'])} |")
    for name, r in R['retrieval'].items():
        if not name.startswith('F4'):
            row(name, r)
    f4 = {k: r for k, r in R['retrieval'].items() if k.startswith('F4')}
    w(f"\nF4 (α = 0; {len(f4)} fixtures), range across fixtures:\n")
    w('| Representation | min purity (k=10) | max purity (k=10) | min / max tie frequency |\n|---|---:|---:|---|')
    for x in REPS:
        ps = [r['T0'][x]['k10']['purity'] for r in f4.values()]
        ts = [r['T0'][x]['k10']['tie_frequency'] for r in f4.values()]
        w(f"| {x} | {f(min(ps))} | {f(max(ps))} | {f(min(ts), 2)} / {f(max(ts), 2)} |")
    k1 = [r['T0']['D_ch']['k1']['purity'] for r in R['retrieval'].values()]
    w(f"\nk = 1 (D_ch, T0): min purity over all retrieval fixtures {f(min(k1))}.\n")

    v3 = v['V3_retrieval_k10_T0']
    w(f"**Differences from comparators** (margin {v3['margin']}, pre-declared):")
    w(f"- Material deficits vs the scalar comparator: {v3['deficits_vs_scalar'] or 'none'}.")
    w(f"- Material deficits vs the Cα-distance baseline: {v3['deficits_vs_ca_distance'] or 'none'}.")
    big = sorted((n, x['diff_vs_scalar']) for n, x in v3['fixtures'].items() if x['diff_vs_scalar'] > 0.05)
    big_ca = sorted((n, x['diff_vs_ca']) for n, x in v3['fixtures'].items() if x['diff_vs_ca'] > 0.05)
    w(f"- D_ch above the scalar comparator by more than the margin: "
      + ', '.join(f'{n} (+{d:.3f})' for n, d in big) + '. These are the designed coordinate-wise '
      'collisions (F1a/F1b) and the mirror set.')
    w(f"- D_ch above the Cα-distance baseline by more than the margin: "
      + (', '.join(f'{n} (+{d:.3f})' for n, d in big_ca) or 'none')
      + '. This is the pre-stated chirality difference: the baseline is reflection-invariant.')
    n2 = R['retrieval']['N2_pooled_F1a_s0_mirror']['T0']
    w(f"- N2 pooled set: {R['retrieval']['N2_pooled_F1a_s0_mirror']['n_labels']} construction "
      f"labels (state × chirality); overall chance purity {n2['ca_distance']['k10']['chance']:.3f}. "
      f"The Cα-distance and scalar scores of {n2['ca_distance']['k10']['purity']:.3f} and "
      f"{n2['scalar_comparator']['k10']['purity']:.3f} are **not** chance-level overall. They are "
      'consistent with resolving one label factor and leaving pairs of classes unresolved: the '
      'Cα-distance baseline cannot separate mirror images, and the scalar summary cannot separate '
      'the two states (the per-frame mean is 60° in both, with sign reversed by the mirror).')
    w('- On all other fixtures, D_ch, D_circ and the Cα-distance baseline reach the same purity. '
      '**The fixtures saturate**: they confirm the designed differences but cannot rank the '
      'candidate against the Cα-distance baseline.\n')

    w('**Historical scan on the paired ±ε control N1** (production BC; shared assignment and '
      'noise, so only ε differs):\n')
    w('| ε | scan BC (R2) |\n|---|---:|')
    for k, x in sorted(R['n1_scan_bc_tau'].items(), key=lambda kv: float(kv[0].split('_e')[1])):
        w(f"| {k.split('_e')[1]} | {x:.4f} |")
    w('\nWith draws paired, the scan output is the same for ε = ±1e-6 and for ±1e-1, but differs '
      'between ε = +1e-3 and −1e-3. This one paired draw did **not reproduce** the ±1e-6 '
      'difference seen in WP-A. Because the WP-A configurations were unpaired, that earlier '
      'comparison cannot isolate a sign effect from draw variation, and a single paired draw does '
      'not establish that the earlier difference was caused entirely by draw variation. The ±1e-3 '
      'difference is likewise observed for this one paired draw only; its mechanism is not '
      'established. The ordered representation is not involved in these values.\n')

    w('## 4. Limitations\n')
    w('- **Constant curvature (κ = 1.2 rad).** D_κ carries no information here (purity = chance, '
      'ties everywhere). The κ block is untested.')
    w('- **Fixed bond lengths (3.8 Å)** and single-window traces. Chain breaks and inter-window '
      'context are not represented.')
    w('- **Constructed labels and simple geometries.** Two well-separated states, mostly W = 5 '
      '(F1b: W = 8), so retrieval saturates for every representation except where an information '
      'difference was designed in.')
    w('- Descriptive scores on fixed fixtures only. No generalization, no superiority claim, no '
      'statement about real proteins.\n')

    w('## 5. Decision record\n')
    w(f"- V1 harness integrity: **{'PASS' if v['V1_harness']['passed'] else 'FAIL'}**")
    w(f"- V2 fidelity: **{'PASS' if v['V2_fidelity']['passed'] else 'FAIL'}**")
    w(f"- V3 deficits vs scalar comparator: {v3['deficits_vs_scalar'] or 'none'}; vs Cα-distance "
      f"baseline: {v3['deficits_vs_ca_distance'] or 'none'}")
    w(f"- **Decision (protocol §9 ordered rules): {v['decision']}.**\n")
    w('Pilot eligibility means only the following:')
    w('- The candidate removed the demonstrated collisions and branch sensitivity.')
    w('- It lost nothing relative to the component it would replace, on these fixtures.')
    w('\nIt is **not** evidence that the candidate outperforms local Cα-distance features. The '
      'fixtures could not discriminate the two outside the designed chirality case. Conditions '
      'for the real-MD pilot:')
    w('- the Cα-distance baseline and the historical scan are included as comparators;')
    w('- a task is used that can separate representations (e.g. retrieval against independent G1 '
      'states);')
    w('- WP-B gates G0–G1 pass independently of this result.\n')
    w('## 6. Test corrections during execution (recorded, contracts unchanged)\n')
    w('Two harness tests failed on first execution because of test-design errors, not harness '
      'defects. Both were corrected without weakening their contracts:')
    w('1. `test_distance_checks_cannot_see_consistent_permutation` demanded bitwise equality of '
      'distances computed in different summation orders. It now requires equality to 1e-15, which '
      'is still the invariance that motivates amendment A.')
    w('2. `test_fully_tied_set_gives_chance_not_index_order` used an example where index-based '
      'tie-breaking happened to land near chance. It now checks the contract directly: the '
      'candidate scores exactly chance under two frame orderings, while index-based tie-breaking '
      'changes with the ordering.')
    w('## 7. Post-run editorial changes\n')
    w('Wording corrections made after the run (no rerun, results unchanged) are logged in '
      '`results/SYNTHETIC_ordered_representation_editorial_log.json`.\n')
    (HERE / 'REPORT.md').write_text('\n'.join(L) + '\n')
    print('REPORT.md written')
    return 0


if __name__ == '__main__':
    sys.exit(main())
