#!/usr/bin/env python3
"""
M3-A2g CONVERGENCE GRID: Systematic max_iter and tolerance testing

Based on trajectory findings, systematically test:
1. max_iter sensitivity (100, 250, 500, 1000)
2. Tolerance criteria (absolute vs relative)
3. Wrapped-normal stress regime (σ=0.6-1.2)
4. Initialization strategies
5. Exact-family and K2 convergence rates

Saves results to CSV for analysis.
"""

import numpy as np
import csv
import sys
from pathlib import Path

# Import from sibling script
sys.path.insert(0, str(Path(__file__).parent))
from m3_a2g_convergence_audit import (
    generate_wrapped_normal,
    generate_product_von_mises_K1,
    generate_product_vm_K2,
    fit_K2_with_full_diagnostics,
    initialize_responsibilities
)
from product_vm_fast_accurate import select_K_via_bic_corrected

#=============================================================================
# MAX_ITER SENSITIVITY
#=============================================================================

def test_max_iter_sensitivity():
    """
    Test convergence rate vs max_iter on representative cases.

    max_iter ∈ {100, 250, 500, 1000}
    Cases: exact-family easy, wrapped σ=0.9, wrapped σ=1.0, clear K2
    """
    print("=" * 80)
    print("MAX_ITER SENSITIVITY TEST")
    print("=" * 80)

    max_iter_values = [100, 250, 500, 1000]

    results = []

    # Generate test cases once
    # Case 1: Exact-family
    mu1 = np.array([0.0, np.pi/4])
    kappa1 = np.array([5.0, 5.0])
    samples1 = generate_product_von_mises_K1(mu1, kappa1, 1000, seed=100)

    # Case 2: Wrapped σ=0.9
    samples2 = generate_wrapped_normal(np.array([0.0, 0.0]), 0.9, 1000, q=2, seed=101)

    # Case 3: Wrapped σ=1.0
    samples3 = generate_wrapped_normal(np.array([0.0, 0.0]), 1.0, 1000, q=2, seed=102)

    # Case 4: Clear K2
    pi4 = np.array([0.5, 0.5])
    mu4 = np.array([[0.0, 0.0], [np.pi/2, np.pi/2]])
    kappa4 = np.array([[20.0, 20.0], [20.0, 20.0]])
    samples4 = generate_product_vm_K2(pi4, mu4, kappa4, 1000, seed=103)

    cases = [
        ('exact_easy', samples1),
        ('wrapped_0.9', samples2),
        ('wrapped_1.0', samples3),
        ('k2_clear', samples4)
    ]

    for case_name, samples in cases:
        print(f"\nCase: {case_name}")
        for max_iter in max_iter_values:
            result = fit_K2_with_full_diagnostics(
                samples, initialization='split', n_init=1,
                max_iter=max_iter, tol_abs=1e-6, seed=100
            )

            best = result['best']
            hist = best['history']

            # Compute improvement after iter 100
            if len(hist) > 100:
                loglik_100 = hist[100]['loglik']
                loglik_final = best['loglik']
                improvement_after_100 = loglik_final - loglik_100
            else:
                improvement_after_100 = np.nan

            # Final delta
            if len(hist) > 1:
                final_delta_abs = hist[-1]['delta_loglik_abs']
                final_delta_rel = hist[-1]['delta_loglik_rel']
            else:
                final_delta_abs = np.nan
                final_delta_rel = np.nan

            results.append({
                'case': case_name,
                'max_iter': max_iter,
                'converged': best['converged'],
                'n_iters': best['n_iters'],
                'final_loglik': best['loglik'],
                'improvement_after_100': improvement_after_100,
                'final_delta_abs': final_delta_abs,
                'final_delta_rel': final_delta_rel,
                'max_kappa': hist[-1]['max_kappa'] if hist else np.nan
            })

            print(f"  max_iter={max_iter:4d}: converged={best['converged']}, "
                  f"iters={best['n_iters']:4d}, final_Δ={final_delta_abs:.2e}")

    # Save to CSV
    output_dir = Path(__file__).parent.parent / "results"
    output_dir.mkdir(exist_ok=True)
    csv_file = output_dir / "m3_a2g_max_iter_sensitivity.csv"

    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved to: {csv_file}")
    return results

#=============================================================================
# TOLERANCE CRITERION COMPARISON
#=============================================================================

def test_tolerance_criteria():
    """
    Compare absolute vs relative tolerance on same cases.

    Absolute: |Δloglik| < tol
    Relative: |Δloglik| / max(1, |loglik|) < tol
    """
    print("\n" + "=" * 80)
    print("TOLERANCE CRITERION COMPARISON")
    print("=" * 80)

    # Generate cases
    mu1 = np.array([0.0, np.pi/4])
    kappa1 = np.array([5.0, 5.0])
    samples1 = generate_product_von_mises_K1(mu1, kappa1, 1000, seed=200)

    samples2 = generate_wrapped_normal(np.array([0.0, 0.0]), 0.9, 1000, q=2, seed=201)

    cases = [
        ('exact_easy', samples1),
        ('wrapped_0.9', samples2)
    ]

    # Test tolerances
    tol_abs_values = [1e-4, 1e-5, 1e-6, 1e-7]
    tol_rel_values = [1e-5, 1e-6, 1e-7, 1e-8]

    results = []

    for case_name, samples in cases:
        print(f"\nCase: {case_name}")

        # Absolute tolerance
        for tol_abs in tol_abs_values:
            result = fit_K2_with_full_diagnostics(
                samples, initialization='split', n_init=1,
                max_iter=1000, tol_abs=tol_abs, tol_rel=None, seed=200
            )
            best = result['best']

            results.append({
                'case': case_name,
                'criterion': 'absolute',
                'tolerance': tol_abs,
                'converged': best['converged'],
                'n_iters': best['n_iters'],
                'final_loglik': best['loglik']
            })

            print(f"  Absolute tol={tol_abs:.0e}: converged={best['converged']}, iters={best['n_iters']:4d}")

        # Relative tolerance
        for tol_rel in tol_rel_values:
            result = fit_K2_with_full_diagnostics(
                samples, initialization='split', n_init=1,
                max_iter=1000, tol_abs=1e-20, tol_rel=tol_rel, seed=200
            )
            best = result['best']

            results.append({
                'case': case_name,
                'criterion': 'relative',
                'tolerance': tol_rel,
                'converged': best['converged'],
                'n_iters': best['n_iters'],
                'final_loglik': best['loglik']
            })

            print(f"  Relative tol={tol_rel:.0e}: converged={best['converged']}, iters={best['n_iters']:4d}")

    # Save
    output_dir = Path(__file__).parent.parent / "results"
    csv_file = output_dir / "m3_a2g_tolerance_comparison.csv"

    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved to: {csv_file}")
    return results

#=============================================================================
# WRAPPED-NORMAL STRESS REGIME
#=============================================================================

def test_wrapped_normal_stress():
    """
    Critical wrapped-normal stress map: σ=0.6-1.2, F=250-2500

    Focus on the convergence-critical σ=0.8-1.1 region.
    """
    print("\n" + "=" * 80)
    print("WRAPPED-NORMAL STRESS REGIME")
    print("=" * 80)

    sigma_values = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
    F_values = [250, 500, 1000, 2500]
    n_reps = 10  # Per (sigma, F) configuration

    results = []

    for sigma in sigma_values:
        for F in F_values:
            print(f"\nσ={sigma:.1f}, F={F}")

            conv_count = 0
            k1_conv_count = 0
            k2_conv_count = 0
            non_conv_count = 0

            for rep in range(n_reps):
                seed = 300 + len(results)
                samples = generate_wrapped_normal(
                    np.array([0.0, 0.0]), sigma, F, q=2, seed=seed
                )

                # Fit with reasonable max_iter=500
                selected_K, bic1, bic2, delta_bic, mu1, kappa1, loglik1, \
                    pi2, mu2, kappa2, loglik2, converged = select_K_via_bic_corrected(
                        samples, n_init=3, max_iter=500, tol=1e-6, seed=seed
                    )

                if converged:
                    conv_count += 1
                    if selected_K == 1:
                        k1_conv_count += 1
                    else:
                        k2_conv_count += 1
                else:
                    non_conv_count += 1

                results.append({
                    'sigma': sigma,
                    'F': F,
                    'rep': rep,
                    'converged': converged,
                    'selected_K': selected_K,
                    'loglik': loglik2,
                    'delta_bic': delta_bic,
                    'max_kappa': np.max(kappa2)
                })

            print(f"  Converged: {conv_count}/{n_reps}, K1: {k1_conv_count}, "
                  f"K2: {k2_conv_count}, Non-conv: {non_conv_count}")

    # Save
    output_dir = Path(__file__).parent.parent / "results"
    csv_file = output_dir / "m3_a2g_stress_regime.csv"

    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved to: {csv_file}")
    return results

#=============================================================================
# INITIALIZATION STRATEGY COMPARISON
#=============================================================================

def test_initialization_strategies():
    """
    Compare 5 initialization strategies on representative cases.

    Strategies: split, random, farthest, embedding, perturb_k1
    """
    print("\n" + "=" * 80)
    print("INITIALIZATION STRATEGY COMPARISON")
    print("=" * 80)

    strategies = ['split', 'random', 'farthest', 'embedding', 'perturb_k1']

    # Generate cases
    mu1 = np.array([0.0, np.pi/4])
    kappa1 = np.array([5.0, 5.0])
    samples1 = generate_product_von_mises_K1(mu1, kappa1, 1000, seed=400)

    samples2 = generate_wrapped_normal(np.array([0.0, 0.0]), 1.0, 1000, q=2, seed=401)

    # Close K2
    sep_rad = 20 * np.pi / 180
    pi3 = np.array([0.8, 0.2])
    mu3 = np.array([[0.0, 0.0], [sep_rad, sep_rad]])
    kappa3 = np.array([[20.0, 20.0], [20.0, 20.0]])
    samples3 = generate_product_vm_K2(pi3, mu3, kappa3, 1000, seed=402)

    cases = [
        ('exact_easy', samples1),
        ('wrapped_1.0', samples2),
        ('k2_close', samples3)
    ]

    results = []

    for case_name, samples in cases:
        print(f"\nCase: {case_name}")

        for strategy in strategies:
            # Run with n_init=5 for robust comparison
            result = fit_K2_with_full_diagnostics(
                samples, initialization=strategy, n_init=5,
                max_iter=500, tol_abs=1e-6, seed=400
            )

            best = result['best']

            results.append({
                'case': case_name,
                'strategy': strategy,
                'best_loglik': best['loglik'],
                'loglik_spread': result['loglik_spread'],
                'n_converged': result['n_converged'],
                'best_n_iters': best['n_iters'],
                'best_converged': best['converged']
            })

            print(f"  {strategy:12s}: loglik={best['loglik']:8.2f}, "
                  f"spread={result['loglik_spread']:.2e}, "
                  f"converged={result['n_converged']}/5")

    # Save
    output_dir = Path(__file__).parent.parent / "results"
    csv_file = output_dir / "m3_a2g_initialization_comparison.csv"

    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved to: {csv_file}")
    return results

#=============================================================================
# MAIN
#=============================================================================

def main():
    """Run all convergence grid tests."""
    print("M3-A2g CONVERGENCE GRID: Systematic Testing")
    print("=" * 80)
    print()

    # Test 1: max_iter sensitivity
    maxiter_results = test_max_iter_sensitivity()

    # Test 2: Tolerance criteria
    tol_results = test_tolerance_criteria()

    # Test 3: Wrapped-normal stress
    stress_results = test_wrapped_normal_stress()

    # Test 4: Initialization strategies
    init_results = test_initialization_strategies()

    print("\n" + "=" * 80)
    print("M3-A2g CONVERGENCE GRID: ALL TESTS COMPLETE")
    print("=" * 80)
    print("\nResults saved to benchmarks/audit/m3_multimodality/results/")

if __name__ == "__main__":
    main()
