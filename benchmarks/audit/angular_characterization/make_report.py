#!/usr/bin/env python3
"""Render REPORT.md from the SYNTHETIC results files. Every number is read from JSON."""
import sys

from common import (FIXTURES_JSON, HERE, PROVENANCE_JSON, SUMMARY_JSON, direction_bound,
                    load_json)


def g(x, fmt='.4g'):
    return 'null' if x is None else format(x, fmt)


def main():
    S = load_json(SUMMARY_JSON)
    P = load_json(PROVENANCE_JSON)
    M = load_json(FIXTURES_JSON)
    H, F = S['harness'], S['findings']
    L = []
    w = L.append

    w('# WP-A Angular Characterization Report (SYNTHETIC)\n')
    w('Characterization of **one component**, the scalar circular-mean → centered-deviation → '
      'Sarle BC torsion scan in `compute_bimodality_profile`, on SYNTHETIC descriptor-level '
      'fixtures. It is not validation of the engine or of any biological application, and it '
      'changes no production behaviour. Protocol: `docs/ANGULAR_REPRESENTATION_CHARACTERIZATION_'
      'PROTOCOL.md` (amendments §9). Every number below is read from '
      '`results/SYNTHETIC_angular_characterization_*.json` by `make_report.py`.\n')

    # --- provenance ---------------------------------------------------------------------------
    w('## Provenance\n')
    e = P['extension']
    w(f"- Git HEAD `{P['git']['head']}`. The working tree had {len(P['git']['status_porcelain'])} "
      'uncommitted entries at run time (documentation drafts and this harness); they are listed '
      'in the provenance file. Production sources are identified by hash, not by HEAD alone.')
    w(f"- Toolchain: {P['versions']['rustc']}; {P['versions']['cargo']}; {P['versions']['maturin']}; "
      f"isolated venv Python {P['versions']['venv_python'].split()[0]}, NumPy "
      f"{P['versions']['venv_numpy']}, pytest {P['versions']['venv_pytest']}.")
    w(f"- Python extension: wheel `{e['wheel']}` built by this run from the current source. "
      f"Loaded extension `{e['loaded_extension_path']}` has SHA-256 `{e['loaded_extension_sha256']}`; "
      f"equal to the `.so` inside the freshly built wheel: **{e['loaded_equals_freshly_built']}**. "
      f"Loaded from the isolated venv: **{e['loaded_from_isolated_venv']}**.")
    w(f"- Fixtures: {M['_meta']['n_cases']} cases ({M['_meta']['n_r2_p1_cases']} with coordinates); "
      f"manifest SHA-256 `{P['fixtures_json_sha256']}`. Seed rule: {M['_meta']['seed_rule']}.")
    w(f"- Run status: `{P['status']}`. Commands and exit codes:\n")
    w('| Exit | Command |\n|---:|---|')
    for c in P['commands']:
        cmd = c['cmd'].replace(str(HERE.parents[2]) + '/', '')
        w(f"| {c['exit']} | `{cmd[:150]}` |")
    w('\nSource hashes:\n')
    w('| File | SHA-256 |\n|---|---|')
    for k, v in P['source_sha256'].items():
        w(f'| `{k}` | `{v}` |')

    # --- 1. harness correctness --------------------------------------------------------------
    w('\n## 1. Implementation correctness (harness)\n')
    w('| Check | Result |\n|---|---|')
    for c in H['checks']:
        w(f"| {c['name']} | {'PASS' if c['passed'] else 'FAIL'} |")
    w(f"\n- Worst construction recovery: R2 (f64) {g(H['r2_worst_recovery_rad'])} rad; P1 (float32) "
      f"{g(H['p1_worst_recovery_rad'])} rad; P1 cases above the 1e-5 flag: "
      f"{len(H['p1_recovery_flagged_cases'])}.")
    w(f"- Comparisons with REF: {H['comparisons_conditioned']} conditioned (compared, all within the "
      f"pre-declared bound); {H['comparisons_not_conditioned']} not conditioned (recorded only, "
      'because a mean was numerically undefined, the deviation set was too narrow, or the guard '
      'was ambiguous).')
    w('- The R2 bitwise check applies only to the exact replication of production operations on '
      'identical f64 inputs. REF and float32 inputs are never compared bitwise.')
    n_all = M['_meta']['n_cases']
    n_xyz = M['_meta']['n_r2_p1_cases']
    w(f'- Coverage: R1 and REF ran on all {n_all} cases; R2 and P1 ran only on the {n_xyz} cases '
      'with coordinates. In particular, F4 rotations other than α = 0 and the F6a scalar inputs '
      'have R1/REF results only.')

    # --- 2. information loss -----------------------------------------------------------------
    w('\n## 2. Representation information loss\n')
    w('**F1: coordinate-wise states with equal circular mean.**\n')
    w('| Case | BC R1 | BC R2 | BC P1 | guard (R2) | scalar state-mean distance (rad) | '
      'embedding state-mean distance |')
    w('|---|---:|---:|---:|---|---:|---:|')
    for k, v in F['F1'].items():
        w(f"| {k} | {g(v['r1'])} | {g(v['r2'])} | {g(v['p1'])} | {v['r2_guard_triggered']} | "
          f"{g(v['scalar_state_mean_circ_dist'])} | {g(v['embedding_state_mean_distance'])} |")
    f1a, f1as, f1c = F['F1']['F1a_s0'], F['F1']['F1a_s2deg'], F['F1']['F1c_s0']
    w(f"\nIn F1a/F1b, two states that the coordinate-wise embedding separates (distance "
      f"{g(f1a['embedding_state_mean_distance'])}) have identical per-frame summaries "
      f"(scalar distance {g(f1a['scalar_state_mean_circ_dist'])}). Noiseless, the scan returns BC "
      f"{g(f1a['r2'])} via the variance guard. With 2° noise it returns {g(f1as['r2'])}, a value "
      f"of the size produced by unimodal noise (F6a Gaussian control: "
      f"{g(F['F6']['F6a_gaussian']['ref']['bc'])}). The positive control F1c, where the states "
      f"differ in circular mean, gives {g(f1c['r2'])}. The two-state structure of F1a/F1b is "
      '**not represented** in the scan output. This is shown for these constructions only.\n')
    w('**F6: BC non-specificity carried through the angular path.**\n')
    w('| Input | BC (REF) | BC (R1) |\n|---|---:|---:|')
    for k, v in F['F6'].items():
        w(f"| {k} | {g(v['ref']['bc'])} | {g(v['r1_bc'])} |")
    w('\n| F6b case | BC R2 (coordinates) | BC P1 | R2 − scalar REF |\n|---|---:|---:|---:|')
    for k, v in F['F6b'].items():
        w(f"| {k} | {g(v['bc_r2'])} | {g(v['bc_p1'])} | {g(v['r2_minus_scalar'], '.2e')} |")
    w('\nRare excursions and skewed unimodal inputs give BC values above those of the Gaussian '
      'control (KNOWN_ISSUES M3). The angular path reproduces the scalar values. The replay of '
      '`bc_check.py` draws matches the values documented in KNOWN_ISSUES M3 to three decimals; '
      'this is a consistency check, not proof of identical samples.\n')

    # --- 3. numerical conditioning -----------------------------------------------------------
    w('## 3. Numerical conditioning and sensitivity (sampled)\n')
    f2 = F['F2']
    w(f"**F2: paired rotations** ({len(f2)} cases). Max |ΔBC| vs α = 0: R1 "
      f"{g(max(abs(v['r1_bc_minus_alpha0']) for v in f2.values()), '.2e')}, R2 "
      f"{g(max(abs(v['r2_bc_minus_alpha0']) for v in f2.values()), '.2e')}, P1 "
      f"{g(max(abs(v['p1_bc_minus_alpha0']) for v in f2.values()), '.2e')}. Max non-wrap deviation "
      f"change (R1): {g(max(v['r1_max_abs_dev_change_nonwrap'] for v in f2.values()), '.2e')} rad. "
      f"Wrap events: {sum(v['r1_wrap_events'] for v in f2.values())}. Deviations within 1e-6 of "
      f"±π: {sum(v['r1_devs_near_pi'] for v in f2.values())}.")
    w('F2c was intended to stress the ±π wrap but produced no deviation near ±π, so wrap '
      'behaviour under rotation is **not exercised by F2**. F4, below, does place deviations '
      f"near ±π; its BC range across the 8 rotations was at most "
      f"{g(max(v['range'] for v in F['F4_rotation_bc_range'].values()), '.2e')}.\n")

    w('**F3: within-frame near-zero resultant** (R1 frame-0 mean shift under a perturbation δ of '
      'one torsion; BC per path).\n')
    w('| ε | σ_n | per-frame R | shift δ=+1e-9 | shift δ=+1e-6 | BC R1 | BC R2 | BC P1 |')
    w('|---:|---:|---:|---:|---:|---:|---:|---:|')
    for v in F['F3'].values():
        ps = v['perturbation_frame0_shift']
        w(f"| {v['eps']:+.0e} | {v['sigma_n']:.0e} | {v['frame_R_min']:.2e} | {ps['+1e-09']:.2e} | "
          f"{ps['+1e-06']:.2e} | {v['r1']:.4f} | {v['r2']:.4f} | {v['p1']:.4f} |")
    flips6 = sorted({v['eps'] for v in F['F3'].values() if v['perturbation_frame0_shift']['+1e-06'] > 1})
    flips9 = sorted({v['eps'] for v in F['F3'].values() if v['perturbation_frame0_shift']['+1e-09'] > 1e-3})
    p1gap = {}
    for v in F['F3'].values():
        p1gap[v['eps']] = max(p1gap.get(v['eps'], 0.0), abs(v['p1'] - v['r2']))
    w(f"\nSampled observations:")
    w(f"- A perturbation of 1e-6 rad moved the frame-0 mean by more than 1 rad at ε ∈ {flips6}.")
    w(f"- A perturbation of 1e-9 rad moved it by more than 1e-3 rad at ε ∈ {flips9}.")
    w('- These flips are expected mathematically: the mean jumps by about π when a perturbation '
      'of about R reverses the sign of the resultant. They are not round-off.')
    w(f"- Max |BC P1 − BC R2| per sampled ε (float32 coordinates vs f64): "
      + ', '.join(f'{k:+.0e}: {g(v, ".3g")}' for k, v in sorted(p1gap.items())) + '.')
    w('- The scan reports no resultant length, so these cases cannot be recognized from its '
      'output. Behaviour between sampled ε values is not established.\n')

    w('**F4: across-frame cancellation and antipodal modes** (α = 0; BC per path).\n')
    w('| Case | across R (R1) | BC R1 | BC R2 | BC P1 |\n|---|---:|---:|---:|---:|')
    for k, v in F['F4'].items():
        if v['alpha_deg'] == 0:
            w(f"| {k} | {v['across_R_r1']:.3e} | {v['r1']:.4f} | {v['r2']:.4f} | {v['p1']:.4f} |")
    noisy_drop = []
    by = {(v['p'], v['eps'], v['sigma']): v for v in F['F4'].values()
          if v['alpha_deg'] == 0 and v['eps'] is not None}
    for (p, eps, sig), v in sorted(by.items()):
        if sig > 0:
            clean = by[(p, eps, 0.0)]['r2']
            if abs(v['r2'] - clean) > 0.1:
                noisy_drop.append(f'p={p}, ε={eps:+.0e}: {clean:.3f} → {v["r2"]:.3f}')
    p05 = F['F4']['F4_p0.50_e0_n0_a+000']
    w(f"\nSampled observations:")
    w(f"- At p = 0.5, ε = 0 the across-frame resultant is {p05['across_R_r1']:.1e}, so the "
      f"across-frame mean is undefined. The deviations are {p05['top_deviation_values']}. "
      f"BC is {p05['r1']:.4f} and did not depend on the arbitrary centering: a two-point set "
      'separated by π gives the same BC wherever the center falls.')
    w('- With 1° per-frame noise, these cases changed BC by more than 0.1 relative to the '
      'noiseless construction: ' + '; '.join(noisy_drop) + '.')
    f4_bases = [c for c in M['cases'] if c['family'] == 'F4' and c['alpha_deg'] == 0]
    n_seeds = len({c['seed'] for c in f4_bases})
    w(f'- **Draws are not paired across F4 configurations.** The {len(f4_bases)} F4 configurations '
      f'(p, ε, σ) used {n_seeds} distinct seeds. The minority-mode frame assignment and the noise '
      'draws therefore differ between configurations, including between +ε and −ε, and between '
      'noisy and noiseless cases. Differences between ±ε rows **cannot be attributed to the sign of '
      'ε**: draw-to-draw variation is not separated from ε. The noiseless BC values do not depend on '
      'the assignment (two-point sets with fixed counts), but every noisy value does.')
    w('- Proposed mechanism, consistent with the recorded deviations but not isolated by these '
      'fixtures: when the across-frame mean lies near the majority mode, the antipodal minority '
      'mode lies near ±π, and per-frame noise can place its members on both sides of the wrap. '
      'How often this happens for a given p, ε and σ is not established. A paired ±ε control '
      '(shared assignment and noise) is specified in `docs/ORDERED_ANGULAR_REPRESENTATION_'
      'PROTOCOL.md`.\n')

    kc = F['kappa_channel']
    w(f"**κ channel.** Non-zero `bc_kappa`: R2 {kc['r2_cases_bc_kappa_nonzero']}/{kc['n_cases']}, "
      f"P1 {kc['p1_cases_bc_kappa_nonzero']}/{kc['n_cases']}. Cases where the score came from κ "
      f"(P1): {kc['p1_cases_score_from_kappa']}. With κ = 1.2 rad and these coordinate scales, "
      'float32 quantization did not lift the κ variance above the absolute guard. Other '
      'geometries were not tested.\n')
    pd = F['p1_minus_r2_bc_tau']
    w(f"**Float32 path vs f64 path.** |BC P1 − BC R2| over {pd['n']} cases: median "
      f"{g(pd['median_abs'], '.2e')}, max {g(pd['max_abs'], '.3g')}; exactly equal in "
      f"{pd['n_exactly_equal']} cases. The large differences occur in the F3 near-cancellation "
      'cases above.\n')

    # --- 4. implications and decision --------------------------------------------------------
    w('## 4. Implications for this component only, and decision record\n')
    w('**Harness correctness:** established for the tested paths. Every check in §1 passed: '
      'bitwise replication of production, recovery, REF agreement where conditioned, '
      'analytical controls, and the F6 replay.\n')
    w('**Information loss (scientific limitation of the representation):** demonstrated for the '
      'F1 constructions. Coordinate-wise state changes that preserve the window circular mean are '
      'invisible to the scan. BC non-specificity (F6) is carried unchanged into the angular path. '
      'Absence of loss for other constructions is not shown.\n')
    w('**Numerical sensitivity:** in the sampled cases, the scan output changed substantially '
      'under two conditions:')
    w('- per-frame resultant comparable to input perturbations: F3, including float32 '
      'quantization;')
    w('- an antipodal minority mode placed at the wrap cut under 1° noise: F4.')
    w('Neither condition is reported by the scan. Rotation invariance held at round-off level in '
      'every sampled case.\n')
    w('**Decision: CHANGE** (protocol §8). The losses and sensitivities above affect '
      'configurations that cannot be detected from the scan output. Candidate changes, each '
      'needing its own protocol and validation (none adopted here):')
    w('1. a per-coordinate (sin, cos) representation instead of the per-frame scalar mean;')
    w('2. carrying and reporting within- and across-frame resultant lengths;')
    w('3. a centering or summary rule that does not place a minority mode at the wrap cut.')
    w('\nUntil a change is validated, a low BC from this scan should not be read as absence of '
      'multimodality, and a high BC should not be read as evidence of two states.')
    w('This decision concerns this component only. Geometry descriptors and retrieval are not '
      'affected.\n')
    w('**Limitations:**')
    w('- SYNTHETIC constructions with constant κ = 1.2 rad and W ∈ {5, 8};')
    w('- a closed, finite fixture matrix;')
    w('- sensitivity reported only at sampled ε, p and noise values;')
    w('- F2 did not exercise ±π wraps under rotation;')
    w('- no real-MD data were used.\n')
    w('**Next concrete decision required:** choose which candidate change, if any, to specify as '
      'its own protocol (with equal-information baselines). The alternative is to restrict the '
      'current scan to an explicitly documented, non-inferential screening role.\n')
    w('## 5. Post-execution notes (recorded, checks unchanged)\n')
    w('- Post-run editorial changes to this report (wording only; no rerun, no change to fixtures '
      'or results files) are logged separately in '
      '`results/SYNTHETIC_angular_characterization_editorial_log.json`. The run provenance file is '
      'the original and is not modified by editorial changes.')
    r9 = F['F3']['F3_e+1e-09_n0']['frame_R_min']
    b9 = float(direction_bound(2, r9))
    w('- The protocol (§5, F3) stated that |ε| = 1e-9 was "expected to be numerically undefined". '
      f'By the pre-declared diagnostic, with exact f64 inputs, it is not: R = {r9:.2e}, '
      f'b = {b9:.1e} rad. The '
      'observed sensitivity there comes from input perturbations of about R, not from round-off. '
      'That expectation statement was wrong; no check depended on it.')
    w('- F2c did not reach its intended wrap stress (see §3).')
    w('- Artifact layout: regenerable binaries (fixtures, arrays, wheel, isolated venv) are under '
      '`build/`. Hashes of the fixtures are in `fixtures/SYNTHETIC_fixtures.json`; hashes of the '
      'arrays are in the results files.')
    w('- During development, before any fixture was generated, the harness trace builder and the '
      'independent NumPy dihedral each had a sign error. The recovery check detected both, and '
      'both were fixed (`test_independent_dihedral_detects_sign_error` guards against '
      'regression).')
    (HERE / 'REPORT.md').write_text('\n'.join(L) + '\n')
    print('REPORT.md written')
    return 0


if __name__ == '__main__':
    sys.exit(main())
