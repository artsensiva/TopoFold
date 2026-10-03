#!/usr/bin/env python3
"""
THREE PRODUCT VON MISES MIXTURE ESTIMATORS
===========================================

Compares regularization strategies for von Mises mixture estimation:

U: Numerical-ceiling EM heuristic (κ ≤ 10⁶)
C: Constrained likelihood (finite justified κ_max)
P: Penalized likelihood (Gamma prior on κ)

All use the corrected Product von Mises infrastructure from M3-A2e.
"""

import numpy as np
from scipy.special import i0, i1, ive
from scipy.optimize import brentq
from scipy.stats import vonmises as scipy_vm

#=============================================================================
# SHARED INFRASTRUCTURE (from product_vm_fast_accurate.py)
#=============================================================================

def log_i0_stable(kappa):
    """Numerically stable log(I0(κ))"""
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
    """Stable von Mises log-density"""
    theta = np.asarray(theta)
    diff = theta - mu
    log_norm = np.log(2*np.pi) + log_i0_stable(kappa)
    return kappa * np.cos(diff) - log_norm

def A_function(kappa):
    """A(κ) = I₁(κ)/I₀(κ)"""
    if np.abs(kappa) < 700:
        return i1(kappa) / i0(kappa)
    else:
        # For large κ: A(κ) ≈ 1 - 1/(2κ)
        return 1 - 1/(2*kappa)

def A_derivative(kappa):
    """A'(κ) = 1 - A² - A/κ"""
    A_val = A_function(kappa)
    return 1 - A_val**2 - A_val/kappa

def kappa_mle_fast(R_bar, n_newton=2, tol=1e-10):
    """Fast Newton-refined κ MLE"""
    if R_bar < 1e-10:
        return 0.0

    # Initial approximation
    if R_bar < 0.53:
        kappa = 2 * R_bar + R_bar**3 + 5 * R_bar**5 / 6
    elif R_bar < 0.85:
        kappa = -0.4 + 1.39 * R_bar + 0.43 / (1 - R_bar)
    else:
        kappa = 1 / (2 * (1 - R_bar))

    # Newton refinement
    for _ in range(n_newton):
        A_val = A_function(kappa)
        A_prime = A_derivative(kappa)

        residual = A_val - R_bar
        if abs(residual) < tol:
            break

        kappa_new = kappa - residual / A_prime
        if kappa_new <= 0:
            kappa_new = kappa / 2.0

        kappa = kappa_new

    return kappa

def fit_product_von_mises_K1(samples):
    """K=1 Product von Mises MLE"""
    n, q = samples.shape

    mu = np.zeros(q)
    kappa = np.zeros(q)

    for dim in range(q):
        theta = samples[:, dim]

        # Circular mean
        S = np.mean(np.sin(theta))
        C = np.mean(np.cos(theta))
        R_bar = np.sqrt(S**2 + C**2)
        mu[dim] = np.arctan2(S, C)

        # MLE concentration
        kappa[dim] = kappa_mle_fast(R_bar)

    # Log-likelihood
    loglik = np.sum([von_mises_logpdf_stable(samples[:, dim], mu[dim], kappa[dim])
                     for dim in range(q)])

    return mu, kappa, loglik

def compute_bic(loglik, n_params, n_samples):
    """BIC = -2*logL + p*log(n)"""
    return -2 * loglik + n_params * np.log(n_samples)

def generate_product_von_mises_K1(mu, kappa, n_samples, seed=None):
    """Generate from K=1 Product von Mises"""
    if seed is not None:
        np.random.seed(seed)

    q = len(mu)
    samples = np.zeros((n_samples, q))

    for dim in range(q):
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
# ESTIMATOR U: NUMERICAL-CEILING EM HEURISTIC
#=============================================================================

def fit_U_estimator(samples, n_init=3, max_iter=100, tol=1e-6, seed=None,
                    track_detailed=False):
    """
    U: Numerical-ceiling EM heuristic

    Uses κ ≤ 10⁶ implicit numerical bound.
    3 random initializations, best local solution by likelihood.

    Returns:
        pi, mu, kappa, loglik, converged, details

        details: dict with optimization diagnostics if track_detailed=True
    """
    kappa_max = 1e6  # Explicit numerical ceiling

    if seed is not None:
        np.random.seed(seed)

    n, q = samples.shape
    K = 2

    best_loglik = -np.inf
    best_result = None
    all_init_logliks = []

    for init_idx in range(n_init):
        # Initialize responsibilities
        if init_idx == 0:
            resp = np.zeros((n, K))
            resp[:n//2, 0] = 1.0
            resp[n//2:, 1] = 1.0
        else:
            resp = np.random.dirichlet([1.0] * K, size=n)

        # EM iterations
        converged = False
        prev_loglik = -np.inf
        iter_history = []

        for em_iter in range(max_iter):
            # M-step
            pi = np.mean(resp, axis=0)
            mu = np.zeros((K, q))
            kappa = np.zeros((K, q))

            for k in range(K):
                for dim in range(q):
                    theta = samples[:, dim]
                    weights = resp[:, k]

                    S = np.sum(weights * np.sin(theta)) / np.sum(weights)
                    C = np.sum(weights * np.cos(theta)) / np.sum(weights)
                    R_bar = np.sqrt(S**2 + C**2)
                    mu[k, dim] = np.arctan2(S, C)

                    kappa_val = kappa_mle_fast(R_bar)
                    kappa[k, dim] = min(kappa_val, kappa_max)  # CEILING

            # E-step
            log_lik_matrix = np.zeros((n, K))
            for k in range(K):
                for dim in range(q):
                    log_lik_matrix[:, k] += von_mises_logpdf_stable(
                        samples[:, dim], mu[k, dim], kappa[k, dim]
                    )

            log_weighted = log_lik_matrix + np.log(pi)[None, :]
            log_sum = np.logaddexp.reduce(log_weighted, axis=1, keepdims=True)
            log_resp = log_weighted - log_sum
            resp = np.exp(log_resp)

            # Likelihood
            loglik = np.sum(log_sum.ravel())

            iter_history.append({
                'iter': em_iter,
                'loglik': loglik,
                'pi': pi.copy(),
                'max_kappa': np.max(kappa)
            })

            # Convergence check
            if abs(loglik - prev_loglik) < tol:
                converged = True
                break

            prev_loglik = loglik

        all_init_logliks.append(loglik)

        # Track best initialization
        if loglik > best_loglik:
            best_loglik = loglik
            best_result = {
                'pi': pi.copy(),
                'mu': mu.copy(),
                'kappa': kappa.copy(),
                'loglik': loglik,
                'converged': converged,
                'n_iter': em_iter + 1,
                'init_idx': init_idx,
                'iter_history': iter_history if track_detailed else None
            }

    details = {
        'estimator': 'U',
        'kappa_constraint': f'κ ≤ {kappa_max:.0e}',
        'n_init': n_init,
        'all_init_logliks': all_init_logliks,
        'loglik_spread': max(all_init_logliks) - min(all_init_logliks),
        'best_init': best_result['init_idx'],
        'boundary_hit': np.any(best_result['kappa'] >= kappa_max * 0.99),
        'pi_min': np.min(best_result['pi']),
        'kappa_max_fitted': np.max(best_result['kappa'])
    }

    if track_detailed:
        details['iter_history'] = best_result['iter_history']

    return (best_result['pi'], best_result['mu'], best_result['kappa'],
            best_result['loglik'], best_result['converged'], details)

#=============================================================================
# ESTIMATOR C: CONSTRAINED LIKELIHOOD
#=============================================================================

def fit_C_estimator(samples, kappa_max=100, n_init=3, max_iter=100, tol=1e-6,
                    seed=None, track_detailed=False):
    """
    C: Constrained likelihood estimator

    Explicit finite concentration domain: 0 ≤ κ ≤ κ_max.

    Args:
        kappa_max: finite concentration bound (scientifically justified)

    Returns:
        pi, mu, kappa, loglik, converged, details
    """
    if seed is not None:
        np.random.seed(seed)

    n, q = samples.shape
    K = 2

    best_loglik = -np.inf
    best_result = None
    all_init_logliks = []

    for init_idx in range(n_init):
        # Initialize responsibilities
        if init_idx == 0:
            resp = np.zeros((n, K))
            resp[:n//2, 0] = 1.0
            resp[n//2:, 1] = 1.0
        else:
            resp = np.random.dirichlet([1.0] * K, size=n)

        # EM iterations
        converged = False
        prev_loglik = -np.inf
        iter_history = []

        for em_iter in range(max_iter):
            # M-step
            pi = np.mean(resp, axis=0)
            mu = np.zeros((K, q))
            kappa = np.zeros((K, q))

            for k in range(K):
                for dim in range(q):
                    theta = samples[:, dim]
                    weights = resp[:, k]

                    S = np.sum(weights * np.sin(theta)) / np.sum(weights)
                    C = np.sum(weights * np.cos(theta)) / np.sum(weights)
                    R_bar = np.sqrt(S**2 + C**2)
                    mu[k, dim] = np.arctan2(S, C)

                    kappa_val = kappa_mle_fast(R_bar)
                    kappa[k, dim] = min(kappa_val, kappa_max)  # EXPLICIT CONSTRAINT

            # E-step
            log_lik_matrix = np.zeros((n, K))
            for k in range(K):
                for dim in range(q):
                    log_lik_matrix[:, k] += von_mises_logpdf_stable(
                        samples[:, dim], mu[k, dim], kappa[k, dim]
                    )

            log_weighted = log_lik_matrix + np.log(pi)[None, :]
            log_sum = np.logaddexp.reduce(log_weighted, axis=1, keepdims=True)
            log_resp = log_weighted - log_sum
            resp = np.exp(log_resp)

            # Likelihood
            loglik = np.sum(log_sum.ravel())

            iter_history.append({
                'iter': em_iter,
                'loglik': loglik,
                'pi': pi.copy(),
                'max_kappa': np.max(kappa)
            })

            # Convergence check
            if abs(loglik - prev_loglik) < tol:
                converged = True
                break

            prev_loglik = loglik

        all_init_logliks.append(loglik)

        # Track best initialization
        if loglik > best_loglik:
            best_loglik = loglik
            best_result = {
                'pi': pi.copy(),
                'mu': mu.copy(),
                'kappa': kappa.copy(),
                'loglik': loglik,
                'converged': converged,
                'n_iter': em_iter + 1,
                'init_idx': init_idx,
                'iter_history': iter_history if track_detailed else None
            }

    details = {
        'estimator': 'C',
        'kappa_constraint': f'κ ≤ {kappa_max}',
        'n_init': n_init,
        'all_init_logliks': all_init_logliks,
        'loglik_spread': max(all_init_logliks) - min(all_init_logliks),
        'best_init': best_result['init_idx'],
        'boundary_hit': np.any(best_result['kappa'] >= kappa_max * 0.99),
        'pi_min': np.min(best_result['pi']),
        'kappa_max_fitted': np.max(best_result['kappa'])
    }

    if track_detailed:
        details['iter_history'] = best_result['iter_history']

    return (best_result['pi'], best_result['mu'], best_result['kappa'],
            best_result['loglik'], best_result['converged'], details)

#=============================================================================
# ESTIMATOR P: PENALIZED LIKELIHOOD
#=============================================================================

def gamma_penalty_logpdf(kappa, shape=2.0, rate=0.1):
    """
    Gamma(α, β) prior on κ in log space.

    Standard literature choice (Hornik & Grün 2014):
    - shape α = 2 (weak prior favoring moderate concentration)
    - rate β = 0.1 (allows wide range but penalizes extreme values)

    log p(κ) = (α-1) log κ - β κ - log Γ(α) + α log β
    """
    from scipy.special import gammaln

    if kappa <= 0:
        return -np.inf

    return ((shape - 1) * np.log(kappa) - rate * kappa
            - gammaln(shape) + shape * np.log(rate))

def fit_P_estimator(samples, prior_shape=2.0, prior_rate=0.1,
                    n_init=3, max_iter=100, tol=1e-6, seed=None,
                    track_detailed=False):
    """
    P: Penalized likelihood estimator (Gamma prior on κ)

    Uses MAP estimation with Gamma(α, β) prior on each κ parameter.

    Penalized objective:
        Q(θ) = log L(θ) + Σ log p(κ)

    Args:
        prior_shape: Gamma shape parameter α (default 2.0)
        prior_rate: Gamma rate parameter β (default 0.1)

    Returns:
        pi, mu, kappa, loglik_raw, converged, details

        Note: loglik_raw is the DATA log-likelihood (unpenalized)
              details['penalized_objective'] contains Q(θ)
    """
    if seed is not None:
        np.random.seed(seed)

    n, q = samples.shape
    K = 2

    best_penalized_obj = -np.inf
    best_result = None
    all_init_objs = []

    for init_idx in range(n_init):
        # Initialize responsibilities
        if init_idx == 0:
            resp = np.zeros((n, K))
            resp[:n//2, 0] = 1.0
            resp[n//2:, 1] = 1.0
        else:
            resp = np.random.dirichlet([1.0] * K, size=n)

        # EM iterations (penalized M-step)
        converged = False
        prev_obj = -np.inf
        iter_history = []

        for em_iter in range(max_iter):
            # M-step (penalized)
            pi = np.mean(resp, axis=0)
            mu = np.zeros((K, q))
            kappa = np.zeros((K, q))

            for k in range(K):
                for dim in range(q):
                    theta = samples[:, dim]
                    weights = resp[:, k]

                    S = np.sum(weights * np.sin(theta)) / np.sum(weights)
                    C = np.sum(weights * np.cos(theta)) / np.sum(weights)
                    R_bar = np.sqrt(S**2 + C**2)
                    mu[k, dim] = np.arctan2(S, C)

                    # MAP estimate (penalized MLE)
                    # For Gamma(α, β) prior, MAP modifies the sufficient statistic:
                    # κ_MAP solves: n_eff * A(κ) = n_eff * R_bar + (α-1) - β·n_eff/κ
                    # Approximation: use MLE then adjust
                    kappa_mle = kappa_mle_fast(R_bar)

                    # Newton refinement for MAP
                    n_eff = np.sum(weights)
                    kappa_map = kappa_mle

                    for _ in range(5):  # MAP iterations
                        A_val = A_function(kappa_map)
                        A_prime = A_derivative(kappa_map)

                        # Penalized score
                        score = n_eff * (R_bar - A_val) + (prior_shape - 1)/kappa_map - prior_rate

                        # Penalized information
                        info = n_eff * A_prime - (prior_shape - 1)/(kappa_map**2)

                        if abs(score) < tol:
                            break

                        kappa_new = kappa_map + score / info
                        if kappa_new <= 0:
                            kappa_new = kappa_map / 2

                        kappa_map = kappa_new

                    kappa[k, dim] = max(kappa_map, 1e-6)  # Numerical floor only

            # E-step (unchanged)
            log_lik_matrix = np.zeros((n, K))
            for k in range(K):
                for dim in range(q):
                    log_lik_matrix[:, k] += von_mises_logpdf_stable(
                        samples[:, dim], mu[k, dim], kappa[k, dim]
                    )

            log_weighted = log_lik_matrix + np.log(pi)[None, :]
            log_sum = np.logaddexp.reduce(log_weighted, axis=1, keepdims=True)
            log_resp = log_weighted - log_sum
            resp = np.exp(log_resp)

            # Raw data log-likelihood
            loglik_raw = np.sum(log_sum.ravel())

            # Penalty term
            penalty = sum(gamma_penalty_logpdf(kappa[k, dim], prior_shape, prior_rate)
                         for k in range(K) for dim in range(q))

            # Penalized objective
            penalized_obj = loglik_raw + penalty

            iter_history.append({
                'iter': em_iter,
                'loglik_raw': loglik_raw,
                'penalty': penalty,
                'penalized_obj': penalized_obj,
                'pi': pi.copy(),
                'max_kappa': np.max(kappa)
            })

            # Convergence check on penalized objective
            if abs(penalized_obj - prev_obj) < tol:
                converged = True
                break

            prev_obj = penalized_obj

        all_init_objs.append(penalized_obj)

        # Track best initialization (by penalized objective)
        if penalized_obj > best_penalized_obj:
            best_penalized_obj = penalized_obj
            best_result = {
                'pi': pi.copy(),
                'mu': mu.copy(),
                'kappa': kappa.copy(),
                'loglik_raw': loglik_raw,
                'penalty': penalty,
                'penalized_obj': penalized_obj,
                'converged': converged,
                'n_iter': em_iter + 1,
                'init_idx': init_idx,
                'iter_history': iter_history if track_detailed else None
            }

    details = {
        'estimator': 'P',
        'penalty': f'Gamma({prior_shape}, {prior_rate}) prior on κ',
        'n_init': n_init,
        'all_init_objs': all_init_objs,
        'obj_spread': max(all_init_objs) - min(all_init_objs),
        'best_init': best_result['init_idx'],
        'penalty_value': best_result['penalty'],
        'penalized_objective': best_result['penalized_obj'],
        'pi_min': np.min(best_result['pi']),
        'kappa_max_fitted': np.max(best_result['kappa'])
    }

    if track_detailed:
        details['iter_history'] = best_result['iter_history']

    # Return RAW log-likelihood (not penalized)
    return (best_result['pi'], best_result['mu'], best_result['kappa'],
            best_result['loglik_raw'], best_result['converged'], details)

#=============================================================================
# MODEL SELECTION WRAPPERS
#=============================================================================

def select_K_via_bic_U(samples, **kwargs):
    """BIC model selection using estimator U"""
    n, q = samples.shape

    # Fit K=1
    mu1, kappa1, loglik1 = fit_product_von_mises_K1(samples)
    p1 = 2 * q
    bic1 = compute_bic(loglik1, p1, n)

    # Fit K=2
    pi2, mu2, kappa2, loglik2, conv2, details2 = fit_U_estimator(samples, **kwargs)
    p2 = 1 + 4 * q
    bic2 = compute_bic(loglik2, p2, n)

    delta_bic = bic2 - bic1
    selected_K = 2 if delta_bic < 0 else 1

    return {
        'selected_K': selected_K,
        'bic1': bic1,
        'bic2': bic2,
        'delta_bic': delta_bic,
        'loglik1': loglik1,
        'loglik2': loglik2,
        'mu1': mu1,
        'kappa1': kappa1,
        'pi2': pi2,
        'mu2': mu2,
        'kappa2': kappa2,
        'converged2': conv2,
        'details2': details2
    }

def select_K_via_bic_C(samples, kappa_max=100, **kwargs):
    """BIC model selection using estimator C"""
    n, q = samples.shape

    # Fit K=1 (constrained)
    mu1, kappa1_unconstrained, loglik1_unconstrained = fit_product_von_mises_K1(samples)
    kappa1 = np.minimum(kappa1_unconstrained, kappa_max)

    # Recompute likelihood with constrained kappa
    loglik1 = np.sum([von_mises_logpdf_stable(samples[:, dim], mu1[dim], kappa1[dim])
                      for dim in range(q)])

    p1 = 2 * q
    bic1 = compute_bic(loglik1, p1, n)

    # Fit K=2
    pi2, mu2, kappa2, loglik2, conv2, details2 = fit_C_estimator(
        samples, kappa_max=kappa_max, **kwargs
    )
    p2 = 1 + 4 * q
    bic2 = compute_bic(loglik2, p2, n)

    delta_bic = bic2 - bic1
    selected_K = 2 if delta_bic < 0 else 1

    return {
        'selected_K': selected_K,
        'bic1': bic1,
        'bic2': bic2,
        'delta_bic': delta_bic,
        'loglik1': loglik1,
        'loglik2': loglik2,
        'mu1': mu1,
        'kappa1': kappa1,
        'pi2': pi2,
        'mu2': mu2,
        'kappa2': kappa2,
        'converged2': conv2,
        'details2': details2
    }

def select_K_via_bic_P(samples, prior_shape=2.0, prior_rate=0.1, **kwargs):
    """
    Model selection using estimator P

    NOTE: Uses BIC on RAW log-likelihoods, NOT penalized objectives.
    This is heuristic - proper Bayesian model selection would use Bayes factors.
    """
    n, q = samples.shape

    # Fit K=1 (penalized)
    # For K=1, penalty is just on the q kappa parameters
    mu1, kappa1_mle, loglik1_raw_mle = fit_product_von_mises_K1(samples)

    # MAP refinement for K=1
    kappa1 = np.zeros(q)
    for dim in range(q):
        kappa_mle = kappa1_mle[dim]
        theta = samples[:, dim]
        R_bar = np.sqrt(np.mean(np.sin(theta))**2 + np.mean(np.cos(theta))**2)

        kappa_map = kappa_mle
        for _ in range(5):
            A_val = A_function(kappa_map)
            A_prime = A_derivative(kappa_map)
            score = n * (R_bar - A_val) + (prior_shape - 1)/kappa_map - prior_rate
            info = n * A_prime - (prior_shape - 1)/(kappa_map**2)

            if abs(score) < 1e-6:
                break

            kappa_new = kappa_map + score / info
            if kappa_new <= 0:
                kappa_new = kappa_map / 2
            kappa_map = kappa_new

        kappa1[dim] = max(kappa_map, 1e-6)

    # Raw likelihood with MAP kappas
    loglik1 = np.sum([von_mises_logpdf_stable(samples[:, dim], mu1[dim], kappa1[dim])
                      for dim in range(q)])

    p1 = 2 * q
    bic1 = compute_bic(loglik1, p1, n)

    # Fit K=2
    pi2, mu2, kappa2, loglik2, conv2, details2 = fit_P_estimator(
        samples, prior_shape=prior_shape, prior_rate=prior_rate, **kwargs
    )
    p2 = 1 + 4 * q
    bic2 = compute_bic(loglik2, p2, n)

    delta_bic = bic2 - bic1
    selected_K = 2 if delta_bic < 0 else 1

    return {
        'selected_K': selected_K,
        'bic1': bic1,
        'bic2': bic2,
        'delta_bic': delta_bic,
        'loglik1': loglik1,
        'loglik2': loglik2,
        'mu1': mu1,
        'kappa1': kappa1,
        'pi2': pi2,
        'mu2': mu2,
        'kappa2': kappa2,
        'converged2': conv2,
        'details2': details2
    }

#=============================================================================
# ESTIMATOR COMPARISON UTILITIES
#=============================================================================

def compare_estimators_on_sample(samples, kappa_max_C=100,
                                 prior_shape_P=2.0, prior_rate_P=0.1,
                                 seed=None):
    """
    Compare all three estimators on a single dataset.

    Returns:
        dict with keys: 'U', 'C', 'P', each containing selection results
    """
    results = {}

    # Estimator U
    results['U'] = select_K_via_bic_U(samples, seed=seed, track_detailed=False)

    # Estimator C
    results['C'] = select_K_via_bic_C(samples, kappa_max=kappa_max_C,
                                      seed=seed, track_detailed=False)

    # Estimator P
    results['P'] = select_K_via_bic_P(samples, prior_shape=prior_shape_P,
                                      prior_rate=prior_rate_P,
                                      seed=seed, track_detailed=False)

    return results

if __name__ == "__main__":
    print("Three Product von Mises Mixture Estimators")
    print("=" * 60)
    print()
    print("Estimators:")
    print("  U: Numerical-ceiling EM (κ ≤ 10⁶)")
    print("  C: Constrained likelihood (finite κ_max)")
    print("  P: Penalized likelihood (Gamma prior)")
    print()
    print("Use compare_estimators_on_sample() to test on data")
