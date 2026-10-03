#!/usr/bin/env python3
"""
M3-A2g NEAR-DUPLICATE COMPONENT INVESTIGATION

Investigate whether slow convergence corresponds to two components
becoming nearly identical (weakly identifiable / duplicate-component solution).

For K1 and wrapped-null fits, quantify:
- Circular distance between K2 component centers
- Difference in kappa
- Component weights
- Change over late EM iterations
"""

import numpy as np
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from m3_a2g_convergence_audit import (
    generate_product_von_mises_K1,
    generate_wrapped_normal,
    fit_K2_with_full_diagnostics
)

def circular_distance(theta1, theta2):
    """Circular distance between two angles."""
    diff = np.abs(theta1 - theta2)
    return min(diff, 2*np.pi - diff)

def component_separation(mu1, mu2, q):
    """
    Toroidal distance between two component centers.

    Returns average circular distance across dimensions.
    """
    dists = []
    for dim in range(q):
        dists.append(circular_distance(mu1[dim], mu2[dim]))
    return np.mean(dists)

def investigate_duplicates():
    """Investigate near-duplicate components in null cases."""
    print("=" * 80)
    print("NEAR-DUPLICATE COMPONENT INVESTIGATION")
    print("=" * 80)

    cases = []

    # K1 cases
    for kappa_val in [1, 5, 20]:
        mu = np.array([0.0, np.pi/4])
        kappa = np.array([kappa_val, kappa_val])
        samples = generate_product_von_mises_K1(mu, kappa, 1000, seed=1100+kappa_val)
        cases.append((f'k1_kappa{kappa_val}', samples, 2))

    # Wrapped-normal cases
    for sigma in [0.6, 0.8, 0.9, 1.0, 1.2]:
        samples = generate_wrapped_normal(np.array([0.0, 0.0]), sigma, 1000, q=2, seed=1110+int(sigma*10))
        cases.append((f'wrapped_{sigma:.1f}', samples, 2))

    results = []

    for case_name, samples, q in cases:
        print(f"\nCase: {case_name}")

        # Fit K2 with tracking
        result = fit_K2_with_full_diagnostics(
            samples, initialization='split', n_init=1,
            max_iter=1000, tol_abs=1e-6, seed=1100
        )

        best = result['best']
        hist = best['history']

        # Extract final state
        mu_final = best['mu']
        kappa_final = best['kappa']
        pi_final = best['pi']

        # Component separation (average across dimensions)
        sep_final = component_separation(mu_final[0], mu_final[1], q)

        # Kappa difference
        kappa_diff = np.abs(kappa_final[0] - kappa_final[1]).mean()

        # Weight ratio
        weight_ratio = min(pi_final) / max(pi_final)

        # Trajectory of separation over iterations
        if len(hist) > 100:
            # Early (first 10 after iter 10)
            hist_early = hist[10:20] if len(hist) > 20 else hist[:10]
            # Late (last 10)
            hist_late = hist[-10:]

            # This would require tracking mu in history, which we don't have
            # For now, just report final state
            sep_trajectory = "not_tracked"
        else:
            sep_trajectory = "too_few_iters"

        print(f"  Final separation: {sep_final:.4f} rad ({sep_final*180/np.pi:.1f}°)")
        print(f"  Final Δκ: {kappa_diff:.2f}")
        print(f"  Final weight ratio: {weight_ratio:.4f}")
        print(f"  Converged: {best['converged']}, iters: {best['n_iters']}")

        results.append({
            'case': case_name,
            'converged': best['converged'],
            'n_iters': best['n_iters'],
            'final_loglik': best['loglik'],
            'final_separation_rad': sep_final,
            'final_separation_deg': sep_final * 180 / np.pi,
            'final_kappa_diff': kappa_diff,
            'final_weight_ratio': weight_ratio,
            'final_min_weight': min(pi_final),
            'final_max_kappa': np.max(kappa_final),
            'effective_n_components': 1.0 / np.sum(pi_final**2)
        })

    # Save
    output_dir = Path(__file__).parent.parent / "results"
    csv_file = output_dir / "m3_a2g_duplicate_components.csv"

    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved to: {csv_file}")

    # Summary
    print("\n" + "=" * 80)
    print("DUPLICATE COMPONENT SUMMARY")
    print("=" * 80)

    # Categorize by separation
    very_close = [r for r in results if r['final_separation_deg'] < 10]
    close = [r for r in results if 10 <= r['final_separation_deg'] < 30]
    moderate = [r for r in results if 30 <= r['final_separation_deg'] < 90]
    separated = [r for r in results if r['final_separation_deg'] >= 90]

    print(f"\nVery close (<10°): {len(very_close)}/{len(results)}")
    if very_close:
        for r in very_close:
            print(f"  {r['case']:20s}: sep={r['final_separation_deg']:5.1f}°, "
                  f"weight_ratio={r['final_weight_ratio']:.3f}, "
                  f"iters={r['n_iters']:4d}, conv={r['converged']}")

    print(f"\nClose (10-30°): {len(close)}/{len(results)}")
    if close:
        for r in close:
            print(f"  {r['case']:20s}: sep={r['final_separation_deg']:5.1f}°, "
                  f"weight_ratio={r['final_weight_ratio']:.3f}, "
                  f"iters={r['n_iters']:4d}, conv={r['converged']}")

    print(f"\nModerate (30-90°): {len(moderate)}/{len(results)}")
    print(f"Separated (≥90°): {len(separated)}/{len(results)}")

    return results

if __name__ == "__main__":
    investigate_duplicates()
