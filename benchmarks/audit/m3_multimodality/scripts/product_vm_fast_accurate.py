#!/usr/bin/env python3
"""
FAST, ACCURATE PRODUCT VON MISES FITTER

Key improvements:
1. Canonical ΔBIC convention: delta_bic = BIC_K2 - BIC_K1
2. Fast Newton-refined κ MLE (not raw approximation)
3. Vectorized EM implementation
4. Independent scipy validation
5. EM monotonicity tracking
"""

import numpy as np
from scipy.special import i0, i1, ive
from scipy.optimize import brentq
from scipy.stats import vonmises as scipy_vm

#=============================================================================
# CANONICAL ΔBIC CONVENTION
#=============================================================================

DELTA_BIC_CONVENTION = """
CANONICAL ΔBIC CONVENTION (Project-Wide)
========================================

Definition:
    delta_bic = BIC_K2 - BIC_K1

Interpretation:
    delta_bic < 0  →  K=2 preferred (BIC_K2 < BIC_K1)
    delta_bic > 0  →  K=1 preferred (BIC_K1 < BIC_K2)

Selection:
    selects_K2 = (delta_bic < 0)
"""

def assert_delta_bic_convention():
    """
    Regression test: duplicate K=2 must select K=1.

    This validates the canonical ΔBIC convention is correctly implemented.
    """
    # Simple test data
    np.random.seed(42)
    samples = np.random.uniform(-np.pi, np.pi, (100, 2))

    # Fit K=1
    mu1, kappa1, loglik1 = fit_product_von_mises_K1(samples)
    bic1 = compute_bic(loglik1, 2*2, 100)

    # Construct duplicate K=2
    pi_dup = np.array([0.5, 0.5])
    mu_dup = np.stack([mu1, mu1])
    kappa_dup = np.stack([kappa1, kappa1])

    # Compute K=2 likelihood manually
    loglik2_dup = 0.0
    for i in range(100):
        log_mix = -np.inf
        for k in range(2):
            log_comp = np.sum([von_mises_logpdf_stable(samples[i, d], mu_dup[k, d], kappa_dup[k, d])
                              for d in range(2)])
            log_mix = np.logaddexp(log_mix, np.log(pi_dup[k]) + log_comp)
        loglik2_dup += log_mix

    bic2_dup = compute_bic(loglik2_dup, 1+4*2, 100)

    # CANONICAL CONVENTION
    delta_bic = bic2_dup - bic1
    selects_K2 = (delta_bic < 0)

    assert not selects_K2, f"ΔBIC convention test FAILED: duplicate K=2 selected (ΔBIC={delta_bic:.2f})"
    assert delta_bic > 0, f"ΔBIC convention test FAILED: ΔBIC should be positive, got {delta_bic:.2f}"

    return True

#=============================================================================
# STABLE VON MISES
#=============================================================================

def log_i0_stable(kappa):
    """Numerically stable log(I0(κ)) = κ + log(i0e(κ))"""
    kappa = np.asarray(kappa)
    result = np.zeros_like(kappa, dtype=float)

    small_mask = np.abs(kappa) < 700
    large_mask = ~small_mask

    if np.any(small_mask):
        result[small_mask] = np.log(i0(kappa[small_mask]))

    if np.any(large_mask):
        result[large_mask] = np.abs(kappa[large_mask]) + np.log(ive(0, kappa[large_mask]))

    return result

def von_mises_logpdf_stable(theta, mu, kappa):
    """Stable log p(θ|μ,κ) = κ*cos(θ-μ) - log(2π) - log(I0(κ))"""
    cos_term = kappa * np.cos(theta - mu)
    norm_term = np.log(2 * np.pi) + log_i0_stable(kappa)
    return cos_term - norm_term

#=============================================================================
# FAST ACCURATE KAPPA MLE (NEWTON-REFINED)
#=============================================================================

def A_function(kappa):
    """A(κ) = I1(κ)/I0(κ), using scaled Bessel functions"""
    if kappa < 1e-10:
        return 0.0
    i0e_val = ive(0, kappa)
    i1e_val = ive(1, kappa)
    if i0e_val < 1e-300:
        return 0.0
    return i1e_val / i0e_val

def A_derivative(kappa):
    """
    A'(κ) = 1 - A(κ)² - A(κ)/κ

    Used for Newton refinement.
    """
    if kappa < 1e-10:
        return 0.5  # Limit as κ→0

    A_val = A_function(kappa)
    return 1.0 - A_val**2 - A_val / kappa

def kappa_mle_approx_initial(R_bar):
    """
    Standard piecewise approximation as INITIAL GUESS ONLY.

    NOT used as final estimate.
    """
    if R_bar < 1e-10:
        return 0.0
    if R_bar >= 1.0 - 1e-10:
        return 1e6

    if R_bar < 0.53:
        return 2*R_bar + R_bar**3 + 5*R_bar**5/6
    elif R_bar < 0.85:
        return -0.4 + 1.39*R_bar + 0.43/(1-R_bar)
    else:
        return 1/(R_bar**3 - 4*R_bar**2 + 3*R_bar)

def kappa_mle_fast(R_bar, n_newton=2, tol=1e-10):
    """
    Fast high-accuracy κ MLE via approximation + Newton refinement.

    Args:
        R_bar: mean resultant length
        n_newton: number of Newton steps (default 2)
        tol: convergence tolerance

    Returns:
        κ estimate (high accuracy)
    """
    if R_bar < 1e-10:
        return 0.0
    if R_bar >= 1.0 - 1e-10:
        return 1e6  # Near-perfect concentration

    # Initial guess from approximation
    kappa = kappa_mle_approx_initial(R_bar)

    # Newton refinement: solve A(κ) = R_bar
    for _ in range(n_newton):
        A_val = A_function(kappa)
        A_prime = A_derivative(kappa)

        if abs(A_prime) < 1e-15:
            break

        residual = A_val - R_bar

        if abs(residual) < tol:
            break

        kappa_new = kappa - residual / A_prime

        # Ensure kappa stays positive
        if kappa_new <= 0:
            kappa_new = kappa / 2.0

        kappa = kappa_new

    return kappa

def kappa_mle_root_solve_reference(R_bar, tol=1e-10):
    """
    Reference root-solved κ MLE for validation.

    Slower but exact to numerical tolerance.
    """
    if R_bar < 1e-10:
        return 0.0
    if R_bar >= 1.0 - 1e-10:
        return 1e6

    # Find upper bracket
    kappa_upper = 1.0
    while A_function(kappa_upper) < R_bar and kappa_upper < 1e6:
        kappa_upper *= 2.0

    if kappa_upper >= 1e6:
        return 1e6

    try:
        kappa_hat = brentq(lambda k: A_function(k) - R_bar, 0, kappa_upper, xtol=tol)
        return kappa_hat
    except ValueError:
        return kappa_mle_approx_initial(R_bar)

#=============================================================================
# K=1 FITTER
#=============================================================================

def fit_product_von_mises_K1(samples, use_fast=True):
    """
    Fit K=1 Product von Mises.

    Args:
        samples: (n, q) array
        use_fast: if True, use fast Newton-refined MLE; if False, use root solve

    Returns: (mu, kappa, loglik)
    """
    n, q = samples.shape
    mu = np.zeros(q)
    kappa = np.zeros(q)

    for dim in range(q):
        theta = samples[:, dim]
        C = np.mean(np.cos(theta))
        S = np.mean(np.sin(theta))
        R_bar = np.sqrt(C**2 + S**2)
        mu[dim] = np.arctan2(S, C)

        if use_fast:
            kappa[dim] = kappa_mle_fast(R_bar)
        else:
            kappa[dim] = kappa_mle_root_solve_reference(R_bar)

    # Log-likelihood (vectorized)
    loglik = np.sum([von_mises_logpdf_stable(samples[:, d], mu[d], kappa[d])
                     for d in range(q)])

    return mu, kappa, loglik

#=============================================================================
# VECTORIZED K=2 EM FITTER
#=============================================================================

def fit_product_von_mises_K2_em_vectorized(samples, n_init=3, max_iter=100, tol=1e-6,
                                           seed=None, kappa_max=None, pi_min=None,
                                           use_fast=True, track_monotonicity=False):
    """
    Vectorized K=2 Product von Mises mixture EM.

    Args:
        samples: (n, q) array
        n_init: number of random initializations
        max_iter: maximum EM iterations
        tol: convergence tolerance
        seed: random seed
        kappa_max: maximum concentration (None = unconstrained, but numerically capped at 1e6)
        pi_min: minimum component weight (None = unconstrained)
        use_fast: if True, use fast Newton-refined κ MLE
        track_monotonicity: if True, return iteration history

    Returns:
        (pi, mu, kappa, loglik, converged, em_history)

        pi: (2,) weights
        mu: (2, q) means
        kappa: (2, q) concentrations
        loglik: final log-likelihood
        converged: bool
        em_history: list of (iter, loglik, pi, max_kappa) if track_monotonicity
    """
    if seed is not None:
        np.random.seed(seed)

    n, q = samples.shape
    K = 2

    # Apply numerical kappa bound if not specified
    if kappa_max is None:
        kappa_max_actual = 1e6
    else:
        kappa_max_actual = kappa_max

    best_loglik = -np.inf
    best_result = None

    for init in range(n_init):
        # Initialize responsibilities
        if init == 0:
            # Split at midpoint
            resp = np.zeros((n, K))
            resp[:n//2, 0] = 1.0
            resp[n//2:, 1] = 1.0
        else:
            # Random responsibilities
            resp = np.random.dirichlet([1.0] * K, size=n)

        em_history = [] if track_monotonicity else None

        for em_iter in range(max_iter):
            # M-step: update parameters
            pi = np.mean(resp, axis=0)

            # Apply weight constraint
            if pi_min is not None:
                pi = np.maximum(pi, pi_min)
                pi /= pi.sum()

            mu = np.zeros((K, q))
            kappa = np.zeros((K, q))

            for k in range(K):
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

                    # Estimate kappa
                    if use_fast:
                        kappa_val = kappa_mle_fast(R_bar)
                    else:
                        kappa_val = kappa_mle_root_solve_reference(R_bar)

                    # Apply concentration constraint
                    kappa_val = min(kappa_val, kappa_max_actual)
                    kappa[k, dim] = kappa_val

            # Compute log-likelihood (vectorized)
            log_lik_matrix = np.zeros((n, K))  # (n, K)

            for k in range(K):
                for dim in range(q):
                    log_lik_matrix[:, k] += von_mises_logpdf_stable(
                        samples[:, dim], mu[k, dim], kappa[k, dim]
                    )

            # Mixture log-likelihood
            log_pi = np.log(pi)  # (K,)
            log_weighted = log_lik_matrix + log_pi[None, :]  # (n, K)

            # Log-sum-exp over components
            loglik = np.sum(np.logaddexp.reduce(log_weighted, axis=1))

            # Track history
            if track_monotonicity:
                em_history.append({
                    'iter': em_iter,
                    'loglik': loglik,
                    'pi': pi.copy(),
                    'max_kappa': np.max(kappa),
                    'min_weight': np.min(pi)
                })

            # Check convergence
            if em_iter > 0:
                if track_monotonicity:
                    prev_loglik = em_history[em_iter - 1]['loglik']
                else:
                    prev_loglik = prev_loglik_stored

                if abs(loglik - prev_loglik) < tol:
                    converged = True
                    break
            else:
                converged = False

            if not track_monotonicity:
                prev_loglik_stored = loglik

            # E-step: update responsibilities (vectorized)
            # log_resp = log p(x_i | k) + log π_k
            log_resp = log_weighted  # (n, K)

            # Normalize: resp[i,k] = exp(log_resp[i,k]) / sum_k' exp(log_resp[i,k'])
            log_sum = np.logaddexp.reduce(log_resp, axis=1, keepdims=True)  # (n, 1)
            log_resp_normalized = log_resp - log_sum
            resp = np.exp(log_resp_normalized)

        else:
            converged = False

        # Keep best
        if loglik > best_loglik:
            best_loglik = loglik
            best_result = (pi, mu, kappa, loglik, converged, em_history)

    return best_result

#=============================================================================
# BIC AND SELECTION (CANONICAL CONVENTION)
#=============================================================================

def compute_bic(loglik, n_params, n_samples):
    """BIC = -2*logL + p*log(n)"""
    return -2 * loglik + n_params * np.log(n_samples)

def select_K_via_bic_corrected(samples, **kwargs):
    """
    Select K=1 or K=2 using BIC with CANONICAL CONVENTION.

    ΔBIC = BIC_K2 - BIC_K1
    selects_K2 = (ΔBIC < 0)

    Returns:
        selected_K, bic1, bic2, delta_bic,
        mu1, kappa1, loglik1,
        pi2, mu2, kappa2, loglik2, converged
    """
    n, q = samples.shape

    # Fit K=1
    mu1, kappa1, loglik1 = fit_product_von_mises_K1(samples)
    p1 = 2 * q
    bic1 = compute_bic(loglik1, p1, n)

    # Fit K=2
    pi2, mu2, kappa2, loglik2, converged, em_hist = fit_product_von_mises_K2_em_vectorized(
        samples, **kwargs
    )
    p2 = 1 + 4 * q
    bic2 = compute_bic(loglik2, p2, n)

    # CANONICAL CONVENTION
    delta_bic = bic2 - bic1  # BIC_K2 - BIC_K1
    selects_K2 = (delta_bic < 0)

    selected_K = 2 if selects_K2 else 1

    return (selected_K, bic1, bic2, delta_bic,
            mu1, kappa1, loglik1,
            pi2, mu2, kappa2, loglik2, converged)

#=============================================================================
# GENERATION
#=============================================================================

def generate_product_von_mises_K1(mu, kappa, n_samples, seed=None):
    """Generate from K=1 Product von Mises"""
    if seed is not None:
        np.random.seed(seed)

    q = len(mu)
    samples = np.zeros((n_samples, q))

    for dim in range(q):
        # Rejection sampling
        angles = []
        while len(angles) < n_samples:
            theta_prop = np.random.uniform(-np.pi, np.pi)
            u = np.random.uniform()
            density = np.exp(kappa[dim] * np.cos(theta_prop - mu[dim]))
            max_density = np.exp(kappa[dim])
            if u < density / max_density:
                angles.append(theta_prop)
        samples[:, dim] = angles[:n_samples]

    return samples

#=============================================================================
# SCIPY VALIDATION
#=============================================================================

def validate_against_scipy(samples, mu_fitted, kappa_fitted):
    """
    Independent validation using scipy.stats.vonmises.

    Returns: (max_logpdf_diff, total_loglik_diff)
    """
    n, q = samples.shape

    max_diff = 0.0
    total_ours = 0.0
    total_scipy = 0.0

    for dim in range(q):
        theta = samples[:, dim]
        mu = mu_fitted[dim]
        kappa = kappa_fitted[dim]

        # Our implementation
        logpdf_ours = von_mises_logpdf_stable(theta, mu, kappa)

        # Scipy implementation
        logpdf_scipy = scipy_vm.logpdf(theta, kappa, loc=mu)

        # Compare
        diff = np.abs(logpdf_ours - logpdf_scipy)
        max_diff = max(max_diff, np.max(diff))

        total_ours += np.sum(logpdf_ours)
        total_scipy += np.sum(logpdf_scipy)

    total_diff = abs(total_ours - total_scipy)

    return max_diff, total_diff
