#!/usr/bin/env python3
"""
M3-A2F ANALYSIS FRAMEWORK
=========================

Corrected collapse analysis and comprehensive estimator comparison.

Key corrections from external review:
1. Analyze ALL K2 fits, not just K2-selected cases
2. Report continuous quantities (pi_min, kappa_max) not arbitrary thresholds
3. Track optimization convergence properly
4. Compare U vs C vs P estimators systematically
"""

import numpy as np
import csv
from product_vm_three_estimators import *

#=============================================================================
# CORRECTED COLLAPSE ANALYSIS
#=============================================================================

def analyze_K2_fit_degeneracy(pi, kappa, n_samples, q,
                              kappa_ceiling=None, kappa1_for_reference=None):
    """
    Comprehensive K2 fit analysis (not just "collapse" binary).

    Reports continuous diagnostic quantities:
    - Minimum weight
    - Maximum concentration
    - Effective counts per component
    - Boundary hit (if ceiling specified)
    - Concentration relative to K1 fit

    Args:
        pi: (2,) weights
        kappa: (2, q) concentrations
        n_samples: data sample size
        q: dimensionality
        kappa_ceiling: concentration bound (if applicable)
        kappa1_for_reference: K=1 fitted kappas for comparison

    Returns:
        dict of diagnostic quantities
    """
    K = 2

    # Component weights
    pi_min = np.min(pi)
    pi_max = np.max(pi)

    # Effective counts
    n_eff = pi * n_samples
    n_eff_min = np.min(n_eff)

    # Concentrations
    kappa_min = np.min(kappa)
    kappa_max = np.max(kappa)
    kappa_mean = np.mean(kappa)

    # Boundary hit (if ceiling exists)
    boundary_hit = False
    if kappa_ceiling is not None:
        boundary_hit = np.any(kappa >= kappa_ceiling * 0.99)

    # Concentration ratio relative to K1
    kappa_ratio_to_K1 = None
    if kappa1_for_reference is not None:
        kappa1_mean = np.mean(kappa1_for_reference)
        if kappa1_mean > 0:
            kappa_ratio_to_K1 = kappa_max / kappa1_mean

    # Diagnostic indicators (continuous, not binary)
    diagnostics = {
        'pi_min': pi_min,
        'pi_max': pi_max,
        'n_eff_min': n_eff_min,
        'kappa_min': kappa_min,
        'kappa_max': kappa_max,
        'kappa_mean': kappa_mean,
        'boundary_hit': boundary_hit,
        'kappa_ratio_to_K1': kappa_ratio_to_K1
    }

    # Traditional "collapse" indicators (for reference, not classification)
    diagnostics['tiny_component'] = (n_eff_min < 5)  # Effective count near few obs
    diagnostics['extreme_concentration'] = (kappa_max > 100)  # Strong departure
    diagnostics['traditional_collapse'] = (pi_min < 0.05 or kappa_max > 100)

    return diagnostics

def stratified_collapse_analysis(results_list):
    """
    Stratified collapse analysis over multiple K2 fits.

    Separates:
    - All K2 fits
    - K2-selected fits (where BIC chose K=2)
    - K1-selected fits (where BIC chose K=1 despite fitting K=2)

    Reports continuous distributions, not just rates above arbitrary thresholds.

    Args:
        results_list: list of dicts, each containing:
            - selected_K
            - pi2, kappa2
            - converged2
            - details2 (with diagnostics)

    Returns:
        dict with stratified statistics
    """
    all_fits = []
    k2_selected = []
    k1_selected = []

    for res in results_list:
        # Extract diagnostics
        diag = analyze_K2_fit_degeneracy(
            res['pi2'], res['kappa2'],
            n_samples=res.get('n_samples', None),
            q=res.get('q', None),
            kappa_ceiling=res['details2'].get('kappa_constraint', None),
            kappa1_for_reference=res.get('kappa1', None)
        )
        diag['converged'] = res['converged2']
        diag['selected_K'] = res['selected_K']

        all_fits.append(diag)

        if res['selected_K'] == 2:
            k2_selected.append(diag)
        else:
            k1_selected.append(diag)

    def summarize_diagnostics(diag_list, label):
        """Summarize diagnostic distributions"""
        if not diag_list:
            return {f'{label}_count': 0}

        pi_mins = [d['pi_min'] for d in diag_list]
        kappa_maxs = [d['kappa_max'] for d in diag_list]
        n_eff_mins = [d['n_eff_min'] for d in diag_list if d['n_eff_min'] is not None]

        convergence_rate = np.mean([d['converged'] for d in diag_list])
        boundary_rate = np.mean([d['boundary_hit'] for d in diag_list])
        traditional_collapse_rate = np.mean([d['traditional_collapse'] for d in diag_list])

        return {
            f'{label}_count': len(diag_list),
            f'{label}_pi_min_median': np.median(pi_mins),
            f'{label}_pi_min_q05': np.percentile(pi_mins, 5),
            f'{label}_pi_min_q95': np.percentile(pi_mins, 95),
            f'{label}_kappa_max_median': np.median(kappa_maxs),
            f'{label}_kappa_max_q95': np.percentile(kappa_maxs, 95),
            f'{label}_convergence_rate': convergence_rate,
            f'{label}_boundary_hit_rate': boundary_rate,
            f'{label}_traditional_collapse_rate': traditional_collapse_rate,
        }

    summary = {}
    summary.update(summarize_diagnostics(all_fits, 'all_K2_fits'))
    summary.update(summarize_diagnostics(k2_selected, 'K2_selected'))
    summary.update(summarize_diagnostics(k1_selected, 'K1_selected'))

    return summary

#=============================================================================
# M3-A2E RESULTS POST-PROCESSOR
#=============================================================================

def reanalyze_m3a2e_checkpoint6(csv_path='/tmp/m3_a2e_checkpoint6_exact_family.csv'):
    """
    Post-process M3-A2e Checkpoint 6 results with corrected analysis.

    Reads the existing CSV (which has aggregate FPR) and adds:
    - Detailed collapse analysis
    - Convergence statistics
    - Continuous diagnostic distributions
    """
    import pandas as pd

    print("M3-A2e Checkpoint 6 Post-Processing")
    print("=" * 60)
    print()

    # Load existing results
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Checkpoint 6 results not found at {csv_path}")
        print("Waiting for M3-A2e benchmark to complete...")
        return None

    print(f"Loaded {len(df)} configurations")
    print()

    # Overall statistics
    total_reps = df['n_reps'].sum()
    total_K2 = df['K2_count'].sum()
    overall_fpr = total_K2 / total_reps

    print(f"Overall FPR: {overall_fpr:.1%} ({total_K2}/{total_reps})")
    print()

    # Stratified analysis
    print("FPR by dimensionality:")
    for q in sorted(df['q'].unique()):
        subset = df[df['q'] == q]
        q_K2 = subset['K2_count'].sum()
        q_total = subset['n_reps'].sum()
        q_fpr = q_K2 / q_total
        print(f"  q={q}: {q_fpr:.1%} ({q_K2}/{q_total})")
    print()

    print("FPR by concentration:")
    for kappa in sorted(df['kappa'].unique()):
        subset = df[df['kappa'] == kappa]
        k_K2 = subset['K2_count'].sum()
        k_total = subset['n_reps'].sum()
        k_fpr = k_K2 / k_total
        print(f"  κ={kappa}: {k_fpr:.1%} ({k_K2}/{k_total})")
    print()

    print("FPR by sample size:")
    for F in sorted(df['F'].unique()):
        subset = df[df['F'] == F]
        f_K2 = subset['K2_count'].sum()
        f_total = subset['n_reps'].sum()
        f_fpr = f_K2 / f_total
        print(f"  F={F}: {f_fpr:.1%} ({f_K2}/{f_total})")
    print()

    # Collapse analysis (from existing data)
    if 'collapse_count' in df.columns:
        total_collapse = df['collapse_count'].sum()
        print(f"Traditional collapse rate (π<5% OR κ>100):")
        print(f"  Among all K2 selections: {total_collapse}/{total_K2} = {total_collapse/total_K2:.1%}")
    print()

    # High-FPR configurations (descriptive, not threshold-based)
    high_configs = df.nlargest(10, 'FPR')
    print("Top 10 configurations by FPR:")
    for _, row in high_configs.iterrows():
        print(f"  q={row['q']}, κ={row['kappa']}, F={row['F']}: "
              f"FPR={row['FPR']:.1%}")
    print()

    return df

def reanalyze_m3a2e_checkpoint9(csv_path='/tmp/m3_a2e_checkpoint9_k2_recovery.csv'):
    """
    Post-process M3-A2e Checkpoint 9 (K2 recovery power).

    Reports power stratified by separation, NOT averaged over arbitrary >=60° threshold.
    """
    import pandas as pd

    print("M3-A2e Checkpoint 9 Post-Processing")
    print("=" * 60)
    print()

    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        print(f"Checkpoint 9 results not found at {csv_path}")
        return None

    print(f"Loaded {len(df)} configurations")
    print()

    # Power by separation
    print("Power by separation angle:")
    for sep in sorted(df['separation_deg'].unique()):
        subset = df[df['separation_deg'] == sep]
        avg_power = subset['power'].mean()
        print(f"  {sep}°: {avg_power:.1%} (mean over occupancy/κ/F)")
    print()

    # Power by concentration
    print("Power by concentration (averaged over separation/occupancy/F):")
    for kappa in sorted(df['kappa'].unique()):
        subset = df[df['kappa'] == kappa]
        avg_power = subset['power'].mean()
        print(f"  κ={kappa}: {avg_power:.1%}")
    print()

    # Power by occupancy
    print("Power by occupancy balance:")
    for _, group in df.groupby(['pi1', 'pi2']):
        pi1 = group['pi1'].iloc[0]
        pi2 = group['pi2'].iloc[0]
        avg_power = group['power'].mean()
        print(f"  ({pi1:.1f}, {pi2:.1f}): {avg_power:.1%}")
    print()

    # Detailed stratification (no arbitrary thresholds)
    print("Power stratification (no arbitrary >=60° grouping):")
    print()
    for sep in sorted(df['separation_deg'].unique()):
        print(f"Separation {sep}°:")
        subset = df[df['separation_deg'] == sep]
        for kappa in sorted(subset['kappa'].unique()):
            k_subset = subset[subset['kappa'] == kappa]
            avg_power = k_subset['power'].mean()
            print(f"  κ={kappa}: {avg_power:.1%}")
        print()

    return df

#=============================================================================
# ESTIMATOR COMPARISON VALIDATION
#=============================================================================

def compare_estimators_exact_family(q_values=[1, 2, 5],
                                    kappa_values=[1, 5, 20],
                                    F_values=[100, 500, 1000],
                                    n_reps=30,
                                    kappa_max_C=100,
                                    prior_shape_P=2.0,
                                    prior_rate_P=0.1):
    """
    Compare U, C, P estimators on exact-family K=1 nulls.

    For each configuration, generate n_reps datasets and fit with all three estimators.

    Returns:
        list of result dicts
    """
    print("Comparing Estimators U vs C vs P on Exact-Family K=1 Nulls")
    print("=" * 60)
    print()

    results = []
    config_idx = 0
    total_configs = len(q_values) * len(kappa_values) * len(F_values)

    for q in q_values:
        for kappa_true in kappa_values:
            for F in F_values:
                config_idx += 1
                print(f"Config {config_idx}/{total_configs}: q={q}, κ={kappa_true}, F={F}")

                mu_true = np.zeros(q)
                kappa_vec = np.full(q, kappa_true, dtype=float)

                config_results = {
                    'q': q,
                    'kappa': kappa_true,
                    'F': F,
                    'n_reps': n_reps,
                    'U': {'K2_count': 0, 'converged_count': 0},
                    'C': {'K2_count': 0, 'converged_count': 0},
                    'P': {'K2_count': 0, 'converged_count': 0}
                }

                for rep in range(n_reps):
                    seed = 50000 + config_idx * 1000 + rep

                    # Generate data
                    samples = generate_product_von_mises_K1(mu_true, kappa_vec, F, seed=seed)

                    # Compare estimators
                    comparison = compare_estimators_on_sample(
                        samples,
                        kappa_max_C=kappa_max_C,
                        prior_shape_P=prior_shape_P,
                        prior_rate_P=prior_rate_P,
                        seed=seed+100000
                    )

                    # Record results
                    for est in ['U', 'C', 'P']:
                        if comparison[est]['selected_K'] == 2:
                            config_results[est]['K2_count'] += 1
                        if comparison[est]['converged2']:
                            config_results[est]['converged_count'] += 1

                # Compute rates
                for est in ['U', 'C', 'P']:
                    config_results[est]['FPR'] = config_results[est]['K2_count'] / n_reps
                    config_results[est]['convergence_rate'] = config_results[est]['converged_count'] / n_reps

                print(f"  U: FPR={config_results['U']['FPR']:.1%}, "
                      f"C: FPR={config_results['C']['FPR']:.1%}, "
                      f"P: FPR={config_results['P']['FPR']:.1%}")

                results.append(config_results)

    return results

def compare_estimators_K2_power(separations_deg=[10, 30, 60, 90],
                                occupancies=[(0.5, 0.5), (0.8, 0.2)],
                                kappa_vals=[5, 20],
                                F_vals=[250, 1000],
                                q_vals=[2],  # Can extend to [1, 2, 5]
                                n_reps=20,
                                kappa_max_C=100):
    """
    Compare U, C, P estimators on K=2 recovery power.
    """
    print("Comparing Estimators U vs C vs P on K=2 Recovery Power")
    print("=" * 60)
    print()

    results = []
    config_idx = 0
    total_configs = (len(separations_deg) * len(occupancies) *
                     len(kappa_vals) * len(F_vals) * len(q_vals))

    for q in q_vals:
        for sep_deg in separations_deg:
            sep_rad = sep_deg * np.pi / 180

            for occ in occupancies:
                pi1_true, pi2_true = occ

                for kappa_true in kappa_vals:
                    for F in F_vals:
                        config_idx += 1
                        print(f"Config {config_idx}/{total_configs}: "
                              f"q={q}, sep={sep_deg}°, occ={occ}, κ={kappa_true}, F={F}")

                        mu1_true = np.zeros(q)
                        mu2_true = np.full(q, sep_rad)
                        kappa1_vec = np.full(q, kappa_true, dtype=float)
                        kappa2_vec = np.full(q, kappa_true, dtype=float)

                        config_results = {
                            'q': q,
                            'separation_deg': sep_deg,
                            'pi1': pi1_true,
                            'pi2': pi2_true,
                            'kappa': kappa_true,
                            'F': F,
                            'n_reps': n_reps,
                            'U': {'K2_recovered': 0},
                            'C': {'K2_recovered': 0},
                            'P': {'K2_recovered': 0}
                        }

                        for rep in range(n_reps):
                            seed = 60000 + config_idx * 1000 + rep

                            # Generate K=2 mixture
                            n1 = int(F * pi1_true)
                            n2 = F - n1

                            samples = np.zeros((F, q))
                            samples[:n1] = generate_product_von_mises_K1(mu1_true, kappa1_vec, n1, seed=seed)
                            samples[n1:] = generate_product_von_mises_K1(mu2_true, kappa2_vec, n2, seed=seed+1)

                            np.random.seed(seed+2)
                            np.random.shuffle(samples)

                            # Compare estimators
                            comparison = compare_estimators_on_sample(
                                samples, kappa_max_C=kappa_max_C, seed=seed+100000
                            )

                            for est in ['U', 'C', 'P']:
                                if comparison[est]['selected_K'] == 2:
                                    config_results[est]['K2_recovered'] += 1

                        # Compute power
                        for est in ['U', 'C', 'P']:
                            config_results[est]['power'] = config_results[est]['K2_recovered'] / n_reps

                        print(f"  U: {config_results['U']['power']:.1%}, "
                              f"C: {config_results['C']['power']:.1%}, "
                              f"P: {config_results['P']['power']:.1%}")

                        results.append(config_results)

    return results

#=============================================================================
# CONSTRAINT SENSITIVITY ANALYSIS
#=============================================================================

def kappa_max_sensitivity_analysis(kappa_max_grid=[20, 50, 100, 200, 500],
                                   q=2, kappa_true=5, F=500, n_reps=30):
    """
    Test sensitivity to κ_max choice in constrained estimator.

    For each κ_max value, report:
    - Exact-family K1 selection rate
    - Boundary hit frequency
    - Fitted parameter distributions
    """
    print("Constrained Estimator Sensitivity Analysis")
    print("=" * 60)
    print()
    print(f"Configuration: q={q}, κ_true={kappa_true}, F={F}")
    print(f"Testing κ_max: {kappa_max_grid}")
    print()

    mu_true = np.zeros(q)
    kappa_vec = np.full(q, kappa_true, dtype=float)

    results = []

    for kappa_max in kappa_max_grid:
        print(f"κ_max = {kappa_max}")

        K2_count = 0
        boundary_hits = 0
        all_kappa_max_fitted = []

        for rep in range(n_reps):
            seed = 70000 + rep

            samples = generate_product_von_mises_K1(mu_true, kappa_vec, F, seed=seed)

            res_C = select_K_via_bic_C(samples, kappa_max=kappa_max, seed=seed+100000)

            if res_C['selected_K'] == 2:
                K2_count += 1

            if res_C['details2']['boundary_hit']:
                boundary_hits += 1

            all_kappa_max_fitted.append(res_C['details2']['kappa_max_fitted'])

        fpr = K2_count / n_reps
        boundary_rate = boundary_hits / n_reps

        results.append({
            'kappa_max': kappa_max,
            'FPR': fpr,
            'boundary_hit_rate': boundary_rate,
            'kappa_fitted_median': np.median(all_kappa_max_fitted),
            'kappa_fitted_q95': np.percentile(all_kappa_max_fitted, 95)
        })

        print(f"  FPR: {fpr:.1%}, Boundary hit: {boundary_rate:.1%}, "
              f"κ_fitted median: {results[-1]['kappa_fitted_median']:.1f}")

    print()
    return results

if __name__ == "__main__":
    print("M3-A2f Analysis Framework")
    print("=" * 60)
    print()
    print("Use reanalyze_m3a2e_checkpoint6() to post-process M3-A2e results")
    print("Use compare_estimators_*() to compare U vs C vs P")
