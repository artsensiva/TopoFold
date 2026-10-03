#!/usr/bin/env python3
"""
M3-A2g EM MONOTONICITY VERIFICATION

Quick verification that EM maintains monotonicity across diverse cases.
Any systematic likelihood decrease indicates implementation bug.
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from m3_a2g_convergence_audit import (
    generate_wrapped_normal,
    generate_product_von_mises_K1,
    generate_product_vm_K2,
    fit_K2_with_full_diagnostics
)

def check_monotonicity():
    """Run diverse cases and check for monotonicity violations."""
    print("=" * 80)
    print("EM MONOTONICITY VERIFICATION")
    print("=" * 80)

    cases = []

    # Case 1: Exact-family
    mu1 = np.array([0.0, np.pi/4])
    kappa1 = np.array([5.0, 5.0])
    samples1 = generate_product_von_mises_K1(mu1, kappa1, 1000, seed=700)
    cases.append(('exact_easy', samples1))

    # Case 2: Wrapped-normal stress
    samples2 = generate_wrapped_normal(np.array([0.0, 0.0]), 1.0, 1000, q=2, seed=701)
    cases.append(('wrapped_1.0', samples2))

    # Case 3: Clear K2
    pi3 = np.array([0.5, 0.5])
    mu3 = np.array([[0.0, 0.0], [np.pi/2, np.pi/2]])
    kappa3 = np.array([[20.0, 20.0], [20.0, 20.0]])
    samples3 = generate_product_vm_K2(pi3, mu3, kappa3, 1000, seed=702)
    cases.append(('k2_clear', samples3))

    # Case 4: Close K2
    sep_rad = 20 * np.pi / 180
    pi4 = np.array([0.8, 0.2])
    mu4 = np.array([[0.0, 0.0], [sep_rad, sep_rad]])
    kappa4 = np.array([[20.0, 20.0], [20.0, 20.0]])
    samples4 = generate_product_vm_K2(pi4, mu4, kappa4, 1000, seed=703)
    cases.append(('k2_close', samples4))

    total_violations = 0
    total_iterations = 0

    for case_name, samples in cases:
        print(f"\nCase: {case_name}")

        result = fit_K2_with_full_diagnostics(
            samples, initialization='split', n_init=3,
            max_iter=500, tol_abs=1e-6, seed=700
        )

        # Check all initializations
        for init_result in result['all_inits']:
            violations = init_result['monotonicity_violations']
            n_iters = init_result['n_iters']

            total_iterations += n_iters
            total_violations += len(violations)

            if violations:
                print(f"  Init {init_result['init_idx']}: {len(violations)} violations in {n_iters} iters")
                for v in violations[:3]:
                    print(f"    Iter {v['iter']}: Δloglik = {v['decrease']:.2e}")

        if not any(r['monotonicity_violations'] for r in result['all_inits']):
            print(f"  ✓ No monotonicity violations ({sum(r['n_iters'] for r in result['all_inits'])} total iters)")

    # Summary
    print("\n" + "=" * 80)
    print("MONOTONICITY VERIFICATION SUMMARY")
    print("=" * 80)
    print(f"Total iterations checked: {total_iterations}")
    print(f"Total violations: {total_violations}")

    if total_violations == 0:
        print("✓ EM MONOTONICITY VERIFIED: No systematic likelihood decrease detected")
    else:
        violation_rate = total_violations / total_iterations * 100
        print(f"⚠ {violation_rate:.2f}% violation rate")
        print("  Note: Tiny violations (~1e-10) are numerical noise, not implementation bugs")

if __name__ == "__main__":
    check_monotonicity()
