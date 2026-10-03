#!/usr/bin/env python3
"""
M3-A2g n_init SENSITIVITY TEST

Test whether split + 2 random (n_init=3) is sufficient.

Compare n_init ∈ {1, 3, 5, 10, 20} on representative cases:
- Exact K1
- Wrapped σ=1.0
- Close K2
- Clear K2
"""

import numpy as np
import csv
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from m3_a2g_convergence_audit import (
    generate_product_von_mises_K1,
    generate_wrapped_normal,
    generate_product_vm_K2,
    fit_K2_with_full_diagnostics
)

def test_ninit_sensitivity():
    """Test n_init sensitivity on 4 representative cases."""
    print("=" * 80)
    print("n_init SENSITIVITY TEST")
    print("=" * 80)

    n_init_values = [1, 3, 5, 10, 20]

    # Generate cases
    cases = []

    # Case 1: Exact K1
    mu1 = np.array([0.0, np.pi/4])
    kappa1 = np.array([5.0, 5.0])
    samples1 = generate_product_von_mises_K1(mu1, kappa1, 1000, seed=800)
    cases.append(('exact_k1', samples1))

    # Case 2: Wrapped σ=1.0 (stress case)
    samples2 = generate_wrapped_normal(np.array([0.0, 0.0]), 1.0, 1000, q=2, seed=801)
    cases.append(('wrapped_1.0', samples2))

    # Case 3: Close K2
    sep_rad = 20 * np.pi / 180
    pi3 = np.array([0.8, 0.2])
    mu3 = np.array([[0.0, 0.0], [sep_rad, sep_rad]])
    kappa3 = np.array([[20.0, 20.0], [20.0, 20.0]])
    samples3 = generate_product_vm_K2(pi3, mu3, kappa3, 1000, seed=802)
    cases.append(('k2_close', samples3))

    # Case 4: Clear K2
    pi4 = np.array([0.5, 0.5])
    mu4 = np.array([[0.0, 0.0], [np.pi/2, np.pi/2]])
    kappa4 = np.array([[20.0, 20.0], [20.0, 20.0]])
    samples4 = generate_product_vm_K2(pi4, mu4, kappa4, 1000, seed=803)
    cases.append(('k2_clear', samples4))

    results = []

    for case_name, samples in cases:
        print(f"\nCase: {case_name}")

        for n_init in n_init_values:
            start_time = time.time()

            # Use max_iter=1000 for adequate budget
            result = fit_K2_with_full_diagnostics(
                samples, initialization='split' if n_init == 1 else 'random',
                n_init=n_init, max_iter=1000, tol_abs=1e-6, seed=800
            )

            runtime = time.time() - start_time

            # Extract metrics
            best = result['best']
            all_logliks = [r['loglik'] for r in result['all_inits']]

            min_loglik = np.min(all_logliks)
            max_loglik = np.max(all_logliks)
            spread_range = max_loglik - min_loglik
            spread_std = np.std(all_logliks) if len(all_logliks) > 1 else 0.0

            # Check convergence
            n_converged = result['n_converged']
            conv_rate = n_converged / n_init if n_init > 0 else 0.0

            results.append({
                'case': case_name,
                'n_init': n_init,
                'best_loglik': best['loglik'],
                'min_loglik': min_loglik,
                'max_loglik': max_loglik,
                'spread_range': spread_range,
                'spread_std': spread_std,
                'n_converged': n_converged,
                'convergence_rate': conv_rate,
                'best_n_iters': best['n_iters'],
                'best_converged': best['converged'],
                'runtime_sec': runtime
            })

            print(f"  n_init={n_init:2d}: best_loglik={best['loglik']:8.2f}, "
                  f"spread_range={spread_range:.2e}, conv={n_converged}/{n_init}, "
                  f"time={runtime:.1f}s")

    # Save
    output_dir = Path(__file__).parent.parent / "results"
    csv_file = output_dir / "m3_a2g_ninit_sensitivity.csv"

    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved to: {csv_file}")
    return results

if __name__ == "__main__":
    test_ninit_sensitivity()
