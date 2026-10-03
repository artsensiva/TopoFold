#!/usr/bin/env python3
"""
M3-A2g EXACT-FAMILY AND K2 CONVERGENCE AUDIT

Reduced-scale convergence audit across the parameter regimes from M3-A2e:
- Exact Product-vM K1 nulls: q={1,2,5}, κ={1,5,20}, F={100,500,1000,2500}
- True K2 alternatives: sep={10,30,60,90}°, occupancy={50/50, 80/20}, κ={5,20}, F={250,1000}

Focus: convergence behavior, not full statistical validation.
Report separately: converged K1, converged K2, non-converged
"""

import numpy as np
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from m3_a2g_convergence_audit import (
    generate_product_von_mises_K1,
    generate_product_vm_K2
)
from product_vm_fast_accurate import select_K_via_bic_corrected

#=============================================================================
# EXACT-FAMILY CONVERGENCE AUDIT
#=============================================================================

def audit_exact_family_convergence():
    """
    Reduced exact-family convergence audit.

    Full grid would be 3×3×4 = 36 configurations × 50 reps = 1800 fits.
    Reduced: 3×3×4 = 36 configs × 10 reps = 360 fits (manageable).
    """
    print("=" * 80)
    print("EXACT PRODUCT-VM K1 CONVERGENCE AUDIT")
    print("=" * 80)

    q_values = [1, 2, 5]
    kappa_values = [1, 5, 20]
    F_values = [100, 500, 1000, 2500]
    n_reps = 10  # Reduced from 50 for this convergence audit

    results = []

    total_configs = len(q_values) * len(kappa_values) * len(F_values)
    config_idx = 0

    for q in q_values:
        for kappa_val in kappa_values:
            for F in F_values:
                config_idx += 1
                print(f"\n[{config_idx}/{total_configs}] q={q}, κ={kappa_val}, F={F}")

                conv_k1 = 0
                conv_k2 = 0
                non_conv = 0
                total_iters = []

                for rep in range(n_reps):
                    seed = 500 + len(results)

                    # Generate K=1 Product-vM
                    mu = np.zeros(q)
                    kappa = np.full(q, kappa_val, dtype=float)
                    samples = generate_product_von_mises_K1(mu, kappa, F, seed=seed)

                    # Fit with improved max_iter=500
                    selected_K, bic1, bic2, delta_bic, mu1, kappa1, loglik1, \
                        pi2, mu2, kappa2, loglik2, converged = select_K_via_bic_corrected(
                            samples, n_init=3, max_iter=500, tol=1e-6, seed=seed
                        )

                    # Classify outcome
                    if converged:
                        if selected_K == 1:
                            conv_k1 += 1
                        else:
                            conv_k2 += 1
                    else:
                        non_conv += 1

                    results.append({
                        'q': q,
                        'kappa': kappa_val,
                        'F': F,
                        'rep': rep,
                        'converged': converged,
                        'selected_K': selected_K,
                        'delta_bic': delta_bic,
                        'loglik': loglik2
                    })

                conv_rate = (conv_k1 + conv_k2) / n_reps * 100
                print(f"  Converged: {conv_k1+conv_k2}/{n_reps} ({conv_rate:.0f}%), "
                      f"K1: {conv_k1}, K2: {conv_k2}, Non-conv: {non_conv}")

    # Save
    output_dir = Path(__file__).parent.parent / "results"
    output_dir.mkdir(exist_ok=True)
    csv_file = output_dir / "m3_a2g_exact_family_convergence.csv"

    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved to: {csv_file}")
    return results

#=============================================================================
# TRUE K2 CONVERGENCE AUDIT
#=============================================================================

def audit_k2_convergence():
    """
    True K=2 convergence audit.

    Grid: sep={10,30,60,90}°, occ={50/50, 80/20}, κ={5,20}, F={250,1000}
    Total: 4 × 2 × 2 × 2 = 32 configs × 10 reps = 320 fits
    """
    print("\n" + "=" * 80)
    print("TRUE K=2 CONVERGENCE AUDIT")
    print("=" * 80)

    separation_deg = [10, 30, 60, 90]
    occupancies = [(0.5, 0.5), (0.8, 0.2)]
    kappa_values = [5, 20]
    F_values = [250, 1000]
    n_reps = 10

    results = []

    total_configs = len(separation_deg) * len(occupancies) * len(kappa_values) * len(F_values)
    config_idx = 0

    for sep_deg in separation_deg:
        for occ in occupancies:
            for kappa_val in kappa_values:
                for F in F_values:
                    config_idx += 1
                    sep_rad = sep_deg * np.pi / 180

                    print(f"\n[{config_idx}/{total_configs}] sep={sep_deg}°, "
                          f"occ={occ[0]:.1f}/{occ[1]:.1f}, κ={kappa_val}, F={F}")

                    conv_k1 = 0
                    conv_k2 = 0
                    non_conv = 0

                    for rep in range(n_reps):
                        seed = 600 + len(results)

                        # Generate bimodal Product-vM
                        pi = np.array([occ[0], occ[1]])
                        mu = np.array([[0.0, 0.0], [sep_rad, sep_rad]])
                        kappa = np.array([[kappa_val, kappa_val], [kappa_val, kappa_val]])
                        samples = generate_product_vm_K2(pi, mu, kappa, F, seed=seed)

                        # Fit
                        selected_K, bic1, bic2, delta_bic, mu1, kappa1, loglik1, \
                            pi2, mu2, kappa2, loglik2, converged = select_K_via_bic_corrected(
                                samples, n_init=3, max_iter=500, tol=1e-6, seed=seed
                            )

                        # Classify
                        if converged:
                            if selected_K == 1:
                                conv_k1 += 1
                            else:
                                conv_k2 += 1
                        else:
                            non_conv += 1

                        results.append({
                            'separation_deg': sep_deg,
                            'occupancy_1': occ[0],
                            'occupancy_2': occ[1],
                            'kappa': kappa_val,
                            'F': F,
                            'rep': rep,
                            'converged': converged,
                            'selected_K': selected_K,
                            'delta_bic': delta_bic,
                            'loglik': loglik2
                        })

                    conv_rate = (conv_k1 + conv_k2) / n_reps * 100
                    print(f"  Converged: {conv_k1+conv_k2}/{n_reps} ({conv_rate:.0f}%), "
                          f"K1: {conv_k1}, K2: {conv_k2}, Non-conv: {non_conv}")

    # Save
    output_dir = Path(__file__).parent.parent / "results"
    csv_file = output_dir / "m3_a2g_k2_convergence.csv"

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
    """Run exact-family and K2 convergence audits."""
    print("M3-A2g EXACT-FAMILY AND K2 CONVERGENCE AUDIT")
    print("=" * 80)
    print()

    # Audit 1: Exact-family
    exact_results = audit_exact_family_convergence()

    # Audit 2: True K2
    k2_results = audit_k2_convergence()

    # Summary
    print("\n" + "=" * 80)
    print("CONVERGENCE AUDIT COMPLETE")
    print("=" * 80)

    # Exact-family summary
    exact_conv = sum(1 for r in exact_results if r['converged'])
    exact_total = len(exact_results)
    print(f"\nExact-family: {exact_conv}/{exact_total} converged ({exact_conv/exact_total*100:.1f}%)")

    # K2 summary
    k2_conv = sum(1 for r in k2_results if r['converged'])
    k2_total = len(k2_results)
    print(f"True K2: {k2_conv}/{k2_total} converged ({k2_conv/k2_total*100:.1f}%)")

    print("\nResults saved to benchmarks/audit/m3_multimodality/results/")

if __name__ == "__main__":
    main()
