#!/usr/bin/env python3
"""
M3-A2g OPTIMIZER CONVERGENCE & INITIALIZATION AUDIT

Systematic audit of Product von Mises EM optimizer to diagnose
the 97% non-convergence rate observed in M3-A2f.

Investigates:
1. Root causes of non-convergence
2. Convergence criterion adequacy (absolute vs relative tolerance)
3. max_iter sensitivity
4. Initialization strategy robustness
5. n_init sensitivity
6. Wrapped-normal stress regime behavior
7. Exact-family and true-K2 convergence rates
8. Degeneracy vs convergence correlation
9. EM monotonicity verification

Scientific autonomy applies: report findings objectively.
"""

import numpy as np
import sys
import os
from pathlib import Path

# Add scripts directory to path to import fitter
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir))

from product_vm_fast_accurate import (
    fit_product_von_mises_K1,
    fit_product_von_mises_K2_em_vectorized,
    generate_product_von_mises_K1,
    select_K_via_bic_corrected,
    von_mises_logpdf_stable,
    kappa_mle_fast
)

#=============================================================================
# DATA GENERATION
#=============================================================================

def generate_wrapped_normal(mu, sigma, n_samples, q=2, seed=None):
    """
    Generate from wrapped normal distribution.

    Wraps N(mu, sigma^2) onto [-π, π].
    Product independence assumed across dimensions.
    """
    if seed is not None:
        np.random.seed(seed)

    samples = np.zeros((n_samples, q))
    for dim in range(q):
        # Generate from normal and wrap
        x = np.random.normal(mu[dim], sigma, n_samples)
        samples[:, dim] = np.arctan2(np.sin(x), np.cos(x))

    return samples

def generate_product_vm_K2(pi, mu, kappa, n_samples, seed=None):
    """
    Generate from K=2 Product von Mises mixture.

    Args:
        pi: (2,) mixture weights
        mu: (2, q) component means
        kappa: (2, q) component concentrations
        n_samples: number of samples
        seed: random seed

    Returns:
        samples: (n_samples, q) array
    """
    if seed is not None:
        np.random.seed(seed)

    q = mu.shape[1]
    samples = np.zeros((n_samples, q))

    # Assign samples to components
    component_labels = np.random.choice(2, size=n_samples, p=pi)

    for i in range(n_samples):
        k = component_labels[i]
        for dim in range(q):
            # von Mises sampling via rejection
            angles = []
            while len(angles) < 1:
                theta_prop = np.random.uniform(-np.pi, np.pi)
                u = np.random.uniform()
                density = np.exp(kappa[k, dim] * np.cos(theta_prop - mu[k, dim]))
                max_density = np.exp(kappa[k, dim])
                if u < density / max_density:
                    angles.append(theta_prop)
            samples[i, dim] = angles[0]

    return samples

#=============================================================================
# ADVANCED EM FITTER WITH DIAGNOSTICS
#=============================================================================

def fit_K2_with_full_diagnostics(samples, initialization='split', n_init=1,
                                  max_iter=100, tol_abs=1e-6, tol_rel=None,
                                  seed=None, kappa_max=1e6):
    """
    Fit K=2 Product-vM with comprehensive diagnostic tracking.

    Args:
        samples: (n, q) data
        initialization: 'split', 'random', 'farthest', 'embedding', 'perturb_k1'
        n_init: number of independent initializations
        max_iter: maximum EM iterations
        tol_abs: absolute convergence tolerance
        tol_rel: relative convergence tolerance (if not None, overrides absolute)
        seed: random seed
        kappa_max: maximum concentration

    Returns:
        results: dict with fitted parameters and diagnostics
    """
    if seed is not None:
        np.random.seed(seed)

    n, q = samples.shape

    best_loglik = -np.inf
    best_result = None
    all_init_results = []

    for init_idx in range(n_init):
        # Initialize
        resp = initialize_responsibilities(samples, initialization, seed=seed)

        # Track iteration history
        history = []
        converged = False
        monotonicity_violations = []

        # Initialize parameters
        pi = np.mean(resp, axis=0)
        mu = np.zeros((2, q))
        kappa = np.zeros((2, q))

        prev_loglik = -np.inf

        for em_iter in range(max_iter):
            # M-step
            pi = np.mean(resp, axis=0)

            for k in range(2):
                for dim in range(q):
                    theta = samples[:, dim]
                    w = resp[:, k]
                    w_sum = w.sum()

                    if w_sum < 1e-10:
                        mu[k, dim] = 0.0
                        kappa[k, dim] = 0.0
                        continue

                    # Weighted circular mean
                    C = np.sum(w * np.cos(theta)) / w_sum
                    S = np.sum(w * np.sin(theta)) / w_sum
                    R_bar = np.sqrt(C**2 + S**2)
                    mu[k, dim] = np.arctan2(S, C)

                    kappa_val = kappa_mle_fast(R_bar)
                    kappa_val = min(kappa_val, kappa_max)
                    kappa[k, dim] = kappa_val

            # Compute log-likelihood
            log_lik_matrix = np.zeros((n, 2))
            for k in range(2):
                for dim in range(q):
                    log_lik_matrix[:, k] += von_mises_logpdf_stable(
                        samples[:, dim], mu[k, dim], kappa[k, dim]
                    )

            log_pi = np.log(pi)
            log_weighted = log_lik_matrix + log_pi[None, :]
            loglik = np.sum(np.logaddexp.reduce(log_weighted, axis=1))

            # Check monotonicity
            if em_iter > 0 and loglik < prev_loglik - 1e-10:
                monotonicity_violations.append({
                    'iter': em_iter,
                    'loglik': loglik,
                    'prev_loglik': prev_loglik,
                    'decrease': prev_loglik - loglik
                })

            # Compute diagnostics
            max_kappa = np.max(kappa)
            min_weight = np.min(pi)

            # Compute effective component count
            eff_n_components = 1.0 / np.sum(pi**2)  # Inverse Simpson index

            # Compute nearest observation distances
            min_dist_to_obs = []
            for k in range(2):
                dists = []
                for i in range(n):
                    # Toroidal distance
                    dist = 0.0
                    for dim in range(q):
                        angle_diff = np.abs(samples[i, dim] - mu[k, dim])
                        angle_diff = min(angle_diff, 2*np.pi - angle_diff)
                        dist += angle_diff**2
                    dists.append(np.sqrt(dist))
                min_dist_to_obs.append(np.min(dists))

            # Record history
            history.append({
                'iter': em_iter,
                'loglik': loglik,
                'delta_loglik_abs': loglik - prev_loglik if em_iter > 0 else np.nan,
                'delta_loglik_rel': (loglik - prev_loglik) / max(1.0, abs(prev_loglik)) if em_iter > 0 else np.nan,
                'pi': pi.copy(),
                'max_kappa': max_kappa,
                'min_weight': min_weight,
                'eff_n_components': eff_n_components,
                'min_dist_to_obs': min(min_dist_to_obs)
            })

            # Check convergence
            if em_iter > 0:
                delta_abs = abs(loglik - prev_loglik)

                # Absolute criterion
                converged_abs = delta_abs < tol_abs

                # Relative criterion (if specified)
                if tol_rel is not None:
                    delta_rel = delta_abs / max(1.0, abs(loglik))
                    converged_rel = delta_rel < tol_rel
                    converged = converged_rel
                else:
                    converged = converged_abs

                if converged:
                    break

            prev_loglik = loglik

            # E-step
            log_resp = log_weighted
            log_sum = np.logaddexp.reduce(log_resp, axis=1, keepdims=True)
            log_resp_normalized = log_resp - log_sum
            resp = np.exp(log_resp_normalized)

        # Record result for this initialization
        init_result = {
            'init_idx': init_idx,
            'pi': pi,
            'mu': mu,
            'kappa': kappa,
            'loglik': loglik,
            'converged': converged,
            'n_iters': len(history),
            'history': history,
            'monotonicity_violations': monotonicity_violations
        }

        all_init_results.append(init_result)

        if loglik > best_loglik:
            best_loglik = loglik
            best_result = init_result

    # Analyze spread across initializations
    logliks = [r['loglik'] for r in all_init_results]
    loglik_spread = np.max(logliks) - np.min(logliks) if len(logliks) > 1 else 0.0

    return {
        'best': best_result,
        'all_inits': all_init_results,
        'loglik_spread': loglik_spread,
        'n_converged': sum(r['converged'] for r in all_init_results),
        'n_init': n_init
    }

def initialize_responsibilities(samples, strategy='split', seed=None):
    """
    Initialize EM responsibilities using specified strategy.

    Strategies:
    - 'split': deterministic split at midpoint
    - 'random': random Dirichlet
    - 'farthest': farthest-point initialization on torus
    - 'embedding': sin/cos embedding + k-means
    - 'perturb_k1': perturb fitted K=1 solution
    """
    if seed is not None:
        np.random.seed(seed)

    n, q = samples.shape

    if strategy == 'split':
        # Split at midpoint (deterministic)
        resp = np.zeros((n, 2))
        resp[:n//2, 0] = 1.0
        resp[n//2:, 1] = 1.0

    elif strategy == 'random':
        # Random soft responsibilities
        resp = np.random.dirichlet([1.0, 1.0], size=n)

    elif strategy == 'farthest':
        # Farthest-point on torus
        # Pick first center randomly
        center1_idx = np.random.randint(n)

        # Find farthest point from center1
        dists = np.zeros(n)
        for i in range(n):
            dist = 0.0
            for dim in range(q):
                angle_diff = np.abs(samples[i, dim] - samples[center1_idx, dim])
                angle_diff = min(angle_diff, 2*np.pi - angle_diff)
                dist += angle_diff**2
            dists[i] = np.sqrt(dist)

        center2_idx = np.argmax(dists)

        # Assign to nearest center
        resp = np.zeros((n, 2))
        for i in range(n):
            dist1 = 0.0
            dist2 = 0.0
            for dim in range(q):
                angle_diff1 = np.abs(samples[i, dim] - samples[center1_idx, dim])
                angle_diff1 = min(angle_diff1, 2*np.pi - angle_diff1)
                dist1 += angle_diff1**2

                angle_diff2 = np.abs(samples[i, dim] - samples[center2_idx, dim])
                angle_diff2 = min(angle_diff2, 2*np.pi - angle_diff2)
                dist2 += angle_diff2**2

            if dist1 < dist2:
                resp[i, 0] = 1.0
            else:
                resp[i, 1] = 1.0

    elif strategy == 'embedding':
        # Transform to sin/cos embedding and use simple k-means
        embedding = np.zeros((n, 2*q))
        for dim in range(q):
            embedding[:, 2*dim] = np.cos(samples[:, dim])
            embedding[:, 2*dim + 1] = np.sin(samples[:, dim])

        # Simple k-means (2 iterations sufficient for initialization)
        # Random initialization of centers
        center_idxs = np.random.choice(n, size=2, replace=False)
        centers = embedding[center_idxs]

        for _ in range(2):
            # Assign to nearest
            labels = np.zeros(n, dtype=int)
            for i in range(n):
                dist1 = np.sum((embedding[i] - centers[0])**2)
                dist2 = np.sum((embedding[i] - centers[1])**2)
                labels[i] = 0 if dist1 < dist2 else 1

            # Update centers
            for k in range(2):
                if np.sum(labels == k) > 0:
                    centers[k] = np.mean(embedding[labels == k], axis=0)

        # Convert to responsibilities
        resp = np.zeros((n, 2))
        for i in range(n):
            resp[i, labels[i]] = 1.0

    elif strategy == 'perturb_k1':
        # Fit K=1, then perturb in opposite directions
        mu1, kappa1, _ = fit_product_von_mises_K1(samples)

        # Create two centers by perturbing π/3 in opposite directions
        perturbation = np.pi / 3

        mu_center1 = (mu1 + perturbation) % (2*np.pi) - np.pi
        mu_center2 = (mu1 - perturbation) % (2*np.pi) - np.pi

        # Assign to nearest (using first dimension as proxy)
        resp = np.zeros((n, 2))
        for i in range(n):
            dist1 = abs(samples[i, 0] - mu_center1[0])
            dist1 = min(dist1, 2*np.pi - dist1)
            dist2 = abs(samples[i, 0] - mu_center2[0])
            dist2 = min(dist2, 2*np.pi - dist2)

            if dist1 < dist2:
                resp[i, 0] = 1.0
            else:
                resp[i, 1] = 1.0

    else:
        raise ValueError(f"Unknown initialization strategy: {strategy}")

    return resp

#=============================================================================
# TRAJECTORY ANALYSIS
#=============================================================================

def analyze_trajectories():
    """
    Analyze iteration trajectories for representative cases.

    Cases:
    1. Exact-family easy null (Product-vM K1, q=2, κ=5, F=1000)
    2. Wrapped-normal stress null (q=2, σ=0.9, F=1000)
    3. Wrapped-normal extreme stress (q=2, σ=1.0, F=1000)
    4. Clear true-K2 (q=2, sep=90°, κ=20, 50/50, F=1000)
    5. Close true-K2 (q=2, sep=20°, κ=20, 80/20, F=1000)
    """
    print("=" * 80)
    print("TRAJECTORY ANALYSIS")
    print("=" * 80)

    results = []

    # Case 1: Exact-family easy null
    print("\n1. Exact-family easy null (Product-vM K1, q=2, κ=5, F=1000)")
    mu1 = np.array([0.0, np.pi/4])
    kappa1 = np.array([5.0, 5.0])
    samples1 = generate_product_von_mises_K1(mu1, kappa1, 1000, seed=42)

    result1 = fit_K2_with_full_diagnostics(samples1, initialization='split',
                                            n_init=1, max_iter=500, seed=42)
    results.append(('exact_easy', result1))
    print_trajectory_summary(result1['best'])

    # Case 2: Wrapped-normal stress σ=0.9
    print("\n2. Wrapped-normal stress (q=2, σ=0.9, F=1000)")
    mu2 = np.array([0.0, 0.0])
    samples2 = generate_wrapped_normal(mu2, 0.9, 1000, q=2, seed=43)

    result2 = fit_K2_with_full_diagnostics(samples2, initialization='split',
                                            n_init=1, max_iter=500, seed=43)
    results.append(('wrapped_0.9', result2))
    print_trajectory_summary(result2['best'])

    # Case 3: Wrapped-normal extreme stress σ=1.0
    print("\n3. Wrapped-normal extreme stress (q=2, σ=1.0, F=1000)")
    mu3 = np.array([0.0, 0.0])
    samples3 = generate_wrapped_normal(mu3, 1.0, 1000, q=2, seed=44)

    result3 = fit_K2_with_full_diagnostics(samples3, initialization='split',
                                            n_init=1, max_iter=500, seed=44)
    results.append(('wrapped_1.0', result3))
    print_trajectory_summary(result3['best'])

    # Case 4: Clear true-K2
    print("\n4. Clear true-K2 (sep=90°, κ=20, 50/50, F=1000)")
    pi4 = np.array([0.5, 0.5])
    mu4 = np.array([[0.0, 0.0], [np.pi/2, np.pi/2]])
    kappa4 = np.array([[20.0, 20.0], [20.0, 20.0]])
    samples4 = generate_product_vm_K2(pi4, mu4, kappa4, 1000, seed=45)

    result4 = fit_K2_with_full_diagnostics(samples4, initialization='split',
                                            n_init=1, max_iter=500, seed=45)
    results.append(('k2_clear', result4))
    print_trajectory_summary(result4['best'])

    # Case 5: Close true-K2
    print("\n5. Close true-K2 (sep=20°, κ=20, 80/20, F=1000)")
    sep_rad = 20 * np.pi / 180
    pi5 = np.array([0.8, 0.2])
    mu5 = np.array([[0.0, 0.0], [sep_rad, sep_rad]])
    kappa5 = np.array([[20.0, 20.0], [20.0, 20.0]])
    samples5 = generate_product_vm_K2(pi5, mu5, kappa5, 1000, seed=46)

    result5 = fit_K2_with_full_diagnostics(samples5, initialization='split',
                                            n_init=1, max_iter=500, seed=46)
    results.append(('k2_close', result5))
    print_trajectory_summary(result5['best'])

    return results

def print_trajectory_summary(result):
    """Print summary of iteration trajectory."""
    hist = result['history']

    print(f"  Converged: {result['converged']}")
    print(f"  Iterations: {result['n_iters']}")
    print(f"  Final loglik: {result['loglik']:.2f}")

    if len(hist) > 0:
        print(f"  First 10 Δloglik: {[h['delta_loglik_abs'] for h in hist[:10]]}")
        if result['n_iters'] >= 10:
            print(f"  Last 10 Δloglik:  {[h['delta_loglik_abs'] for h in hist[-10:]]}")

        print(f"  Final max κ: {hist[-1]['max_kappa']:.2f}")
        print(f"  Final min weight: {hist[-1]['min_weight']:.4f}")
        print(f"  Final eff. components: {hist[-1]['eff_n_components']:.2f}")

    if result['monotonicity_violations']:
        print(f"  ⚠ Monotonicity violations: {len(result['monotonicity_violations'])}")
        for v in result['monotonicity_violations'][:3]:
            print(f"    Iter {v['iter']}: decrease = {v['decrease']:.2e}")

#=============================================================================
# MAIN
#=============================================================================

def main():
    """Run M3-A2g convergence audit."""
    print("M3-A2g OPTIMIZER CONVERGENCE & INITIALIZATION AUDIT")
    print("=" * 80)
    print()

    # System information
    print("Environment:")
    print(f"  Python: {sys.version.split()[0]}")
    print(f"  NumPy: {np.__version__}")
    print(f"  Working directory: {os.getcwd()}")
    print()

    # Run trajectory analysis
    trajectory_results = analyze_trajectories()

    print("\n" + "=" * 80)
    print("TRAJECTORY ANALYSIS COMPLETE")
    print("=" * 80)

    # Save trajectory details to file
    output_dir = Path(__file__).parent.parent / "results"
    output_dir.mkdir(exist_ok=True)

    traj_file = output_dir / "m3_a2g_trajectories.txt"
    with open(traj_file, 'w') as f:
        f.write("M3-A2g ITERATION TRAJECTORIES\n")
        f.write("=" * 80 + "\n\n")

        for case_name, result in trajectory_results:
            f.write(f"\nCase: {case_name}\n")
            f.write("-" * 80 + "\n")

            best = result['best']
            f.write(f"Converged: {best['converged']}\n")
            f.write(f"Iterations: {best['n_iters']}\n")
            f.write(f"Final loglik: {best['loglik']:.6f}\n\n")

            f.write("iter,loglik,delta_abs,delta_rel,max_kappa,min_weight,eff_n_comp\n")
            for h in best['history']:
                f.write(f"{h['iter']},{h['loglik']:.6f},{h['delta_loglik_abs']:.2e},"
                       f"{h['delta_loglik_rel']:.2e},{h['max_kappa']:.2f},"
                       f"{h['min_weight']:.6f},{h['eff_n_components']:.4f}\n")

            f.write("\n")

    print(f"\nTrajectory details saved to: {traj_file}")

    print("\n" + "=" * 80)
    print("M3-A2g CHECKPOINT: Phase 1 (Trajectory Analysis) COMPLETE")
    print("=" * 80)

if __name__ == "__main__":
    main()
