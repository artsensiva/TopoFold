#!/usr/bin/env python3
"""
M3-A2g EXTENDED max_iter TEST

For difficult stress fixtures, test max_iter ∈ {500, 1000, 2000}
to determine behavior and practical limits.

Focus on:
- Wrapped σ=1.0, F ∈ {1000, 2500}
- Wrapped σ=0.6-0.7, F ∈ {1000, 2500}
"""

import numpy as np
import csv
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from m3_a2g_convergence_audit import generate_wrapped_normal
from product_vm_fast_accurate import select_K_via_bic_corrected

def test_extended_maxiter():
    """Test extended max_iter on stress cases."""
    print("=" * 80)
    print("EXTENDED MAX_ITER TEST (STRESS CASES)")
    print("=" * 80)

    max_iter_values = [500, 1000, 2000]

    # Stress configurations
    stress_configs = [
        (0.6, 1000),
        (0.6, 2500),
        (0.7, 1000),
        (0.7, 2500),
        (1.0, 1000),
        (1.0, 2500),
    ]

    n_reps = 10
    results = []

    for sigma, F in stress_configs:
        print(f"\nσ={sigma:.1f}, F={F}")

        for max_iter in max_iter_values:
            conv_count = 0
            k1_conv = 0
            k2_conv = 0
            non_conv = 0
            total_iters = []
            runtimes = []

            for rep in range(n_reps):
                seed = 1000 + len(results)
                samples = generate_wrapped_normal(
                    np.array([0.0, 0.0]), sigma, F, q=2, seed=seed
                )

                start_time = time.time()

                selected_K, bic1, bic2, delta_bic, mu1, kappa1, loglik1, \
                    pi2, mu2, kappa2, loglik2, converged = select_K_via_bic_corrected(
                        samples, n_init=3, max_iter=max_iter, tol=1e-6, seed=seed
                    )

                runtime = time.time() - start_time
                runtimes.append(runtime)

                if converged:
                    conv_count += 1
                    if selected_K == 1:
                        k1_conv += 1
                    else:
                        k2_conv += 1
                else:
                    non_conv += 1

                # Record (simplified - would need to track iterations from fit)
                results.append({
                    'sigma': sigma,
                    'F': F,
                    'max_iter': max_iter,
                    'rep': rep,
                    'converged': converged,
                    'selected_K': selected_K,
                    'loglik': loglik2,
                    'delta_bic': delta_bic,
                    'max_kappa': np.max(kappa2),
                    'runtime_sec': runtime
                })

            conv_rate = conv_count / n_reps * 100
            median_runtime = np.median(runtimes)

            print(f"  max_iter={max_iter:4d}: conv={conv_count:2d}/{n_reps} ({conv_rate:4.0f}%), "
                  f"K1={k1_conv:2d}, K2={k2_conv:2d}, non-conv={non_conv:2d}, "
                  f"time={median_runtime:.1f}s")

    # Save
    output_dir = Path(__file__).parent.parent / "results"
    csv_file = output_dir / "m3_a2g_extended_maxiter.csv"

    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved to: {csv_file}")
    return results

if __name__ == "__main__":
    test_extended_maxiter()
