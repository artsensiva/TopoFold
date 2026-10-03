#!/usr/bin/env python3
"""
M3-A2E FULL BENCHMARK VALIDATION
=================================

Long-running validation checkpoints for Product von Mises fitter.

Estimated runtime: 8-16 hours total
- Checkpoint 6: Exact-family K=1 benchmark (3-6 hours)
- Checkpoint 7: Likelihood collapse analysis (embedded in 6)
- Checkpoint 8: Constraint sensitivity (1.5-3 hours)
- Checkpoint 9: Controlled K=2 recovery (2-4 hours)
- Checkpoint 10: Independent wrapped-normal robustness (30 min)
- Checkpoint 11: Dependent toroidal robustness (30 min)
"""

import sys
import time
import csv
import numpy as np
from product_vm_fast_accurate import *

print("=" * 80)
print("M3-A2E FULL BENCHMARK VALIDATION")
print("=" * 80)
print()
print("WARNING: This validation suite will run for approximately 8-16 hours")
print()

start_time = time.time()

#=============================================================================
# CHECKPOINT 6: EXACT-FAMILY K=1 BENCHMARK
#=============================================================================

print("CHECKPOINT 6: Exact-Family K=1 Benchmark")
print("-" * 80)
print("Testing Product von Mises FPR on exact-family K=1 nulls")
print()

# Configuration grid
q_values = [1, 2, 5]
kappa_values = [1, 5, 20]
F_values = [100, 250, 500, 1000, 2500]
n_reps_per_config = 50

results_ch6 = []
total_configs = len(q_values) * len(kappa_values) * len(F_values)
config_idx = 0

for q in q_values:
    for kappa_true in kappa_values:
        for F in F_values:
            config_idx += 1
            print(f"Config {config_idx}/{total_configs}: q={q}, κ={kappa_true}, F={F}")

            # Prepare true parameters
            mu_true = np.zeros(q)
            kappa_vec = np.full(q, kappa_true, dtype=float)

            # Run replicates
            K2_count = 0
            collapse_count = 0

            for rep in range(n_reps_per_config):
                seed = 10000 + config_idx * 1000 + rep

                # Generate exact K=1 Product von Mises data
                samples = generate_product_von_mises_K1(mu_true, kappa_vec, F, seed=seed)

                # Fit and select
                result = select_K_via_bic_corrected(samples, seed=seed+100000)
                selected_K, bic1, bic2, delta_bic, mu1, kappa1, loglik1, pi2, mu2, kappa2, loglik2, converged = result

                # Count K=2 selections (false positives)
                if selected_K == 2:
                    K2_count += 1

                    # Analyze for collapse (Checkpoint 7)
                    min_weight = np.min(pi2)
                    max_kappa = np.max(kappa2)

                    # Collapse criteria: minority weight < 5% OR max κ > 100
                    if min_weight < 0.05 or max_kappa > 100:
                        collapse_count += 1

            # Compute FPR for this configuration
            fpr = K2_count / n_reps_per_config
            collapse_rate = collapse_count / n_reps_per_config

            results_ch6.append({
                'q': q,
                'kappa': kappa_true,
                'F': F,
                'n_reps': n_reps_per_config,
                'K2_count': K2_count,
                'FPR': fpr,
                'collapse_count': collapse_count,
                'collapse_rate': collapse_rate
            })

            print(f"  FPR: {fpr:.1%} ({K2_count}/{n_reps_per_config}), Collapse: {collapse_rate:.1%}")

# Save Checkpoint 6 results
with open('/tmp/m3_a2e_checkpoint6_exact_family.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['q', 'kappa', 'F', 'n_reps', 'K2_count', 'FPR', 'collapse_count', 'collapse_rate'])
    writer.writeheader()
    writer.writerows(results_ch6)

print()
print(f"Checkpoint 6 complete. Results saved to: /tmp/m3_a2e_checkpoint6_exact_family.csv")
print()

# Summary statistics
overall_K2 = sum(r['K2_count'] for r in results_ch6)
overall_total = sum(r['n_reps'] for r in results_ch6)
overall_fpr = overall_K2 / overall_total

print(f"OVERALL FPR: {overall_fpr:.1%} ({overall_K2}/{overall_total})")
print()

#=============================================================================
# CHECKPOINT 7: LIKELIHOOD COLLAPSE ANALYSIS
#=============================================================================

print("CHECKPOINT 7: Likelihood Collapse Analysis")
print("-" * 80)

# Detailed analysis of false K=2 selections from Checkpoint 6
if overall_K2 > 0:
    print(f"Analyzing {overall_K2} false K=2 selections from Checkpoint 6")
    print(f"Collapse rate (π_min < 5% OR κ_max > 100): {sum(r['collapse_count'] for r in results_ch6) / overall_K2:.1%}")
    print()

    # Group by configuration
    high_fpr_configs = [r for r in results_ch6 if r['FPR'] > 0.15]
    if high_fpr_configs:
        print(f"High-FPR configurations (>15%): {len(high_fpr_configs)}")
        for r in high_fpr_configs[:10]:
            print(f"  q={r['q']}, κ={r['kappa']}, F={r['F']}: FPR={r['FPR']:.1%}, Collapse={r['collapse_rate']:.1%}")
    else:
        print("No configurations exceeded 15% FPR")
else:
    print("No false K=2 selections in Checkpoint 6 - collapse analysis not applicable")

print()

#=============================================================================
# CHECKPOINT 9: CONTROLLED K=2 RECOVERY
#=============================================================================

print("CHECKPOINT 9: Controlled K=2 Recovery")
print("-" * 80)
print("Testing Product von Mises power to detect true K=2 mixtures")
print()

# Configuration grid
separations_deg = [10, 20, 30, 60, 90, 180]
occupancies = [(0.5, 0.5), (0.8, 0.2)]
kappa_vals = [5, 20, 100]
F_vals = [100, 500, 1000]
n_reps_k2 = 30

results_ch9 = []
total_configs_k2 = len(separations_deg) * len(occupancies) * len(kappa_vals) * len(F_vals)
config_idx_k2 = 0

for sep_deg in separations_deg:
    sep_rad = sep_deg * np.pi / 180

    for occ in occupancies:
        pi1_true, pi2_true = occ

        for kappa_true in kappa_vals:
            for F in F_vals:
                config_idx_k2 += 1
                print(f"Config {config_idx_k2}/{total_configs_k2}: sep={sep_deg}°, occ={occ}, κ={kappa_true}, F={F}")

                # Fixed q=2 for K=2 recovery
                q = 2

                # True component parameters
                mu1_true = np.array([0.0, 0.0])
                mu2_true = np.array([sep_rad, sep_rad])
                kappa1_vec = np.full(q, kappa_true, dtype=float)
                kappa2_vec = np.full(q, kappa_true, dtype=float)

                # Run replicates
                K2_recovered = 0

                for rep in range(n_reps_k2):
                    seed = 20000 + config_idx_k2 * 1000 + rep

                    # Generate K=2 mixture data
                    n1 = int(F * pi1_true)
                    n2 = F - n1

                    samples = np.zeros((F, q))
                    samples[:n1] = generate_product_von_mises_K1(mu1_true, kappa1_vec, n1, seed=seed)
                    samples[n1:] = generate_product_von_mises_K1(mu2_true, kappa2_vec, n2, seed=seed+1)

                    # Shuffle
                    np.random.seed(seed+2)
                    np.random.shuffle(samples)

                    # Fit and select
                    result = select_K_via_bic_corrected(samples, seed=seed+100000)
                    selected_K = result[0]

                    if selected_K == 2:
                        K2_recovered += 1

                # Compute power
                power = K2_recovered / n_reps_k2

                results_ch9.append({
                    'separation_deg': sep_deg,
                    'pi1': pi1_true,
                    'pi2': pi2_true,
                    'kappa': kappa_true,
                    'F': F,
                    'n_reps': n_reps_k2,
                    'K2_recovered': K2_recovered,
                    'power': power
                })

                print(f"  Power: {power:.1%} ({K2_recovered}/{n_reps_k2})")

# Save Checkpoint 9 results
with open('/tmp/m3_a2e_checkpoint9_k2_recovery.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['separation_deg', 'pi1', 'pi2', 'kappa', 'F', 'n_reps', 'K2_recovered', 'power'])
    writer.writeheader()
    writer.writerows(results_ch9)

print()
print(f"Checkpoint 9 complete. Results saved to: /tmp/m3_a2e_checkpoint9_k2_recovery.csv")
print()

# Summary: Power at clear separations
clear_sep_results = [r for r in results_ch9 if r['separation_deg'] >= 60]
if clear_sep_results:
    avg_power_clear = np.mean([r['power'] for r in clear_sep_results])
    print(f"Average power at clear separations (≥60°): {avg_power_clear:.1%}")
else:
    print("No clear separation results")

print()

#=============================================================================
# CHECKPOINT 10: INDEPENDENT WRAPPED-NORMAL ROBUSTNESS
#=============================================================================

print("CHECKPOINT 10: Independent Wrapped-Normal Robustness")
print("-" * 80)
print("Testing on independent wrapped-normal data (shape misspecification)")
print()

# Configuration
q_vals_wn = [2, 5]
sigma_vals = [0.5, 1.0, 2.0]  # Wrapped normal std dev
F_vals_wn = [250, 1000]
n_reps_wn = 30

results_ch10 = []

for q in q_vals_wn:
    for sigma in sigma_vals:
        for F in F_vals_wn:
            print(f"Config: q={q}, σ={sigma}, F={F}")

            K2_count = 0

            for rep in range(n_reps_wn):
                seed = 30000 + rep

                # Generate independent wrapped normal
                np.random.seed(seed)
                Z = np.random.normal(0, sigma, (F, q))
                samples = (Z % (2*np.pi)) - np.pi  # Wrap to [-π, π]

                # Fit and select
                result = select_K_via_bic_corrected(samples, seed=seed+100000)
                selected_K = result[0]

                if selected_K == 2:
                    K2_count += 1

            fpr = K2_count / n_reps_wn

            results_ch10.append({
                'q': q,
                'sigma': sigma,
                'F': F,
                'n_reps': n_reps_wn,
                'K2_count': K2_count,
                'FPR': fpr
            })

            print(f"  FPR: {fpr:.1%} ({K2_count}/{n_reps_wn})")

# Save results
with open('/tmp/m3_a2e_checkpoint10_wrapped_normal.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['q', 'sigma', 'F', 'n_reps', 'K2_count', 'FPR'])
    writer.writeheader()
    writer.writerows(results_ch10)

print()
print(f"Checkpoint 10 complete. Results saved to: /tmp/m3_a2e_checkpoint10_wrapped_normal.csv")
print()

#=============================================================================
# FINAL SUMMARY
#=============================================================================

elapsed_time = time.time() - start_time
hours = int(elapsed_time // 3600)
minutes = int((elapsed_time % 3600) // 60)

print("=" * 80)
print("M3-A2E FULL BENCHMARK COMPLETE")
print("=" * 80)
print()
print(f"Total runtime: {hours}h {minutes}m")
print()
print("Results files:")
print("  - /tmp/m3_a2e_checkpoint6_exact_family.csv")
print("  - /tmp/m3_a2e_checkpoint9_k2_recovery.csv")
print("  - /tmp/m3_a2e_checkpoint10_wrapped_normal.csv")
print()
print("Summary:")
print(f"  Checkpoint 6 (Exact-family FPR): {overall_fpr:.1%}")
if clear_sep_results:
    print(f"  Checkpoint 9 (K=2 recovery power): {avg_power_clear:.1%} (clear separations)")
overall_fpr_wn = sum(r['K2_count'] for r in results_ch10) / sum(r['n_reps'] for r in results_ch10)
print(f"  Checkpoint 10 (Wrapped-normal FPR): {overall_fpr_wn:.1%}")
print()
print("Next step: Analyze results and determine final classification")
print("Classification options:")
print("  1. IMPLEMENTATION STILL INVALID")
print("  2. VALID FITTER / MODEL ROBUSTNESS UNRESOLVED")
print("  3. REJECTED AS PRIMARY DETECTOR")
print("  4. ACCEPTABLE SCREENING MODEL")
print("  5. PRIMARY CANDIDATE FOR NEXT CALIBRATION STAGE")
print()
