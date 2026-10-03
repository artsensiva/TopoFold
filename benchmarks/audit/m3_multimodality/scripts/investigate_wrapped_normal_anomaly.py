#!/usr/bin/env python3
"""
WRAPPED-NORMAL ANOMALY INVESTIGATION
=====================================

M3-A2e Checkpoint 10 showed:
- q=2, σ=1.0, F=1000: 15/30 K2 selections = 50%
- All neighboring configs: 0/30 K2 selections

Investigate whether this is:
1. Reproducible shape-misspecification sensitivity
2. Seed artifact
3. Wrapping/origin artifact
4. Optimizer artifact
5. Finite-sample transition
6. Coding error

Test grid:
- q ∈ {1, 2, 5}
- σ ∈ {0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5}
- F ∈ {250, 500, 1000, 2500}
- Higher replication for critical cases
"""

import numpy as np
from product_vm_fast_accurate import select_K_via_bic_corrected

def generate_wrapped_normal_conventional(mu, sigma, n_samples, q, seed=None):
    """
    Generate wrapped normal using CONVENTIONAL wrapping:

    θ = ((Z + π) % (2π)) - π

    This centers the wrapped distribution at μ on [-π, π].
    """
    if seed is not None:
        np.random.seed(seed)

    Z = np.random.normal(mu, sigma, (n_samples, q))

    # Conventional wrapping to [-π, π]
    theta = ((Z + np.pi) % (2*np.pi)) - np.pi

    return theta

def verify_unimodality_wrapped_normal(sigma, q=1):
    """
    Verify wrapped normal is unimodal for given σ and q.

    For univariate wrapped normal centered at 0, it's unimodal if σ is not too large.
    For multivariate independent wrapped normal, each marginal is unimodal.

    Returns: (is_unimodal, reason)
    """
    # For wrapped normal, theoretical unimodality condition depends on σ
    # Rule of thumb: σ < π/2 ensures strong unimodality
    # σ ≈ 1 is borderline - need to check numerically

    if sigma < 1.0:
        return True, f"σ={sigma} < 1.0, clearly unimodal"
    elif sigma > 2.0:
        return False, f"σ={sigma} > 2.0, likely multimodal"
    else:
        # Numerical check for borderline cases
        # Wrapped normal density at θ:
        # f(θ) = Σₖ φ((θ - 2πk)/σ) / σ
        # Check if it has multiple local maxima

        theta_grid = np.linspace(-np.pi, np.pi, 1000)
        density = np.zeros_like(theta_grid)

        # Sum over wraps (k = -3 to 3 should be sufficient)
        for k in range(-3, 4):
            density += np.exp(-0.5 * ((theta_grid - 2*np.pi*k) / sigma)**2)

        density /= (sigma * np.sqrt(2*np.pi))

        # Find local maxima
        maxima = []
        for i in range(1, len(density)-1):
            if density[i] > density[i-1] and density[i] > density[i+1]:
                maxima.append(i)

        if len(maxima) == 1:
            return True, f"σ={sigma}, numerical check: single mode at θ={theta_grid[maxima[0]]:.3f}"
        else:
            return False, f"σ={sigma}, numerical check: {len(maxima)} modes detected"

def test_wrapping_origin_invariance(q=2, sigma=1.0, F=1000, n_tests=10):
    """
    Test that detector decision is invariant to angular origin choice.

    Generate data with two different wrapping conventions and verify
    same selection result.
    """
    results_conv = []
    results_alt = []

    for rep in range(n_tests):
        seed = 80000 + rep

        # Generate unwrapped Gaussian
        np.random.seed(seed)
        Z = np.random.normal(0, sigma, (F, q))

        # Conventional wrapping
        theta_conv = ((Z + np.pi) % (2*np.pi)) - np.pi

        # Alternative wrapping (shift by π/4)
        theta_alt = ((Z + np.pi + np.pi/4) % (2*np.pi)) - np.pi

        # Fit both
        res_conv = select_K_via_bic_corrected(theta_conv, seed=seed+100000)
        res_alt = select_K_via_bic_corrected(theta_alt, seed=seed+200000)

        results_conv.append(res_conv[0])  # selected_K
        results_alt.append(res_alt[0])

    # Check if results match
    agreements = sum(1 for i in range(n_tests) if results_conv[i] == results_alt[i])

    return agreements, n_tests, results_conv, results_alt

def focused_replication(q=2, sigma=1.0, F=1000, n_reps=100):
    """
    High-replication test at the anomalous configuration.

    Use independent seeds for each replicate.
    """
    K2_count = 0
    all_delta_bic = []
    all_converged = []

    for rep in range(n_reps):
        # Configuration-specific independent seed
        seed = 90000 + q * 10000 + int(sigma * 1000) + F + rep

        # Generate data
        np.random.seed(seed)
        Z = np.random.normal(0, sigma, (F, q))
        theta = ((Z + np.pi) % (2*np.pi)) - np.pi

        # Fit
        result = select_K_via_bic_corrected(theta, seed=seed+100000)
        selected_K, bic1, bic2, delta_bic = result[0:4]
        converged = result[11]

        if selected_K == 2:
            K2_count += 1

        all_delta_bic.append(delta_bic)
        all_converged.append(converged)

    fpr = K2_count / n_reps
    convergence_rate = np.mean(all_converged)

    return {
        'q': q,
        'sigma': sigma,
        'F': F,
        'n_reps': n_reps,
        'K2_count': K2_count,
        'FPR': fpr,
        'convergence_rate': convergence_rate,
        'delta_bic_mean': np.mean(all_delta_bic),
        'delta_bic_median': np.median(all_delta_bic),
        'delta_bic_std': np.std(all_delta_bic)
    }

def full_sigma_transition(q=2, F=1000, n_reps=50):
    """
    Map the full σ transition surface at fixed q, F.
    """
    sigma_grid = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0]

    results = []

    for sigma in sigma_grid:
        print(f"Testing σ={sigma:.1f}...")

        K2_count = 0

        for rep in range(n_reps):
            seed = 95000 + int(sigma * 1000) + F + rep

            np.random.seed(seed)
            Z = np.random.normal(0, sigma, (F, q))
            theta = ((Z + np.pi) % (2*np.pi)) - np.pi

            result = select_K_via_bic_corrected(theta, seed=seed+100000)
            if result[0] == 2:
                K2_count += 1

        fpr = K2_count / n_reps

        # Verify unimodality
        is_unimodal, reason = verify_unimodality_wrapped_normal(sigma, q=1)

        results.append({
            'sigma': sigma,
            'FPR': fpr,
            'K2_count': K2_count,
            'n_reps': n_reps,
            'unimodal': is_unimodal,
            'unimodality_note': reason
        })

        print(f"  FPR: {fpr:.1%}, Unimodal: {is_unimodal}")

    return results

if __name__ == "__main__":
    print("=" * 80)
    print("WRAPPED-NORMAL ANOMALY INVESTIGATION")
    print("=" * 80)
    print()

    print("Original M3-A2e finding:")
    print("  q=2, σ=1.0, F=1000: 15/30 = 50% K2 selections")
    print("  q=2, σ=1.0, F=250:   0/30 = 0%")
    print("  q=2, σ=0.5, F=1000:  0/30 = 0%")
    print("  q=2, σ=2.0, F=1000:  0/30 = 0%")
    print()

    # Test 1: Verify unimodality
    print("TEST 1: Verify wrapped-normal unimodality")
    print("-" * 80)
    for sigma in [0.5, 1.0, 1.5, 2.0]:
        is_uni, reason = verify_unimodality_wrapped_normal(sigma)
        print(f"σ={sigma}: {reason}")
    print()

    # Test 2: Wrapping origin invariance
    print("TEST 2: Wrapping origin invariance (n=10)")
    print("-" * 80)
    agreements, total, conv, alt = test_wrapping_origin_invariance(q=2, sigma=1.0, F=1000, n_tests=10)
    print(f"Agreements: {agreements}/{total}")
    print(f"Conventional: {conv}")
    print(f"Alternative:  {alt}")
    if agreements == total:
        print("✓ Detector is origin-invariant")
    else:
        print("✗ WARNING: Detector shows origin dependence!")
    print()

    # Test 3: High-replication at anomalous point
    print("TEST 3: High-replication at anomalous point (n=100)")
    print("-" * 80)
    result_anomaly = focused_replication(q=2, sigma=1.0, F=1000, n_reps=100)
    print(f"Configuration: q={result_anomaly['q']}, σ={result_anomaly['sigma']}, F={result_anomaly['F']}")
    print(f"FPR: {result_anomaly['FPR']:.1%} ({result_anomaly['K2_count']}/{result_anomaly['n_reps']})")
    print(f"Convergence rate: {result_anomaly['convergence_rate']:.1%}")
    print(f"ΔBIC: mean={result_anomaly['delta_bic_mean']:.2f}, median={result_anomaly['delta_bic_median']:.2f}, std={result_anomaly['delta_bic_std']:.2f}")

    if abs(result_anomaly['FPR'] - 0.5) < 0.1:
        print("✓ Anomaly REPRODUCED (FPR ≈ 50%)")
    else:
        print(f"✗ Anomaly NOT reproduced (FPR = {result_anomaly['FPR']:.1%}, expected ≈50%)")
    print()

    # Test 4: Map full σ transition
    print("TEST 4: Full σ transition surface (q=2, F=1000, n=50 per σ)")
    print("-" * 80)
    transition_results = full_sigma_transition(q=2, F=1000, n_reps=50)
    print()

    print("Transition surface:")
    print(f"{'σ':>6} {'FPR':>8} {'K2/n':>10} {'Unimodal':>10}")
    print("-" * 40)
    for res in transition_results:
        print(f"{res['sigma']:>6.1f} {res['FPR']:>7.1%} {res['K2_count']:>4}/{res['n_reps']:<4} {'Yes' if res['unimodal'] else 'No':>10}")

    # Find peak
    peak_idx = max(range(len(transition_results)), key=lambda i: transition_results[i]['FPR'])
    peak = transition_results[peak_idx]
    print()
    print(f"Peak FPR: σ={peak['sigma']:.1f}, FPR={peak['FPR']:.1%}")

    print()
    print("=" * 80)
    print("INVESTIGATION COMPLETE")
    print("=" * 80)
