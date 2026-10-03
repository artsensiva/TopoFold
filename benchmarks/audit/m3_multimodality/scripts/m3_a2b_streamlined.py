#!/usr/bin/env python3
"""M3-A2b STREAMLINED VALIDATION
Focus on essential corrected tests without the stuck issues.
"""

import sys
sys.path.insert(0, '/tmp')

import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal, norm
from scipy.optimize import minimize

from m3_a_diagnostic import (
    wrap, fit_product_von_mises_mixture_em, fit_product_von_mises_single,
    product_von_mises_logpdf
)

print("M3-A2b STREAMLINED VALIDATION", flush=True)
print("=" * 80, flush=True)
print()

# Wilson score CI
def wilson_ci(k, n, alpha=0.05):
    if n == 0:
        return (0.0, 1.0)
    z = norm.ppf(1 - alpha/2)
    p_hat = k / n
    denom = 1 + z**2 / n
    center = (p_hat + z**2 / (2*n)) / denom
    margin = z * np.sqrt((p_hat * (1 - p_hat) / n + z**2 / (4*n**2))) / denom
    return (max(0, center - margin), min(1, center + margin))

# Wrapped normal generator
def generate_wrapped_normal(mu, Sigma, n_samples, seed=None):
    if seed is not None:
        np.random.seed(seed)
    q = len(mu)
    Z = multivariate_normal.rvs(mean=mu, cov=Sigma, size=n_samples)
    if Z.ndim == 1:
        Z = Z.reshape(1, -1)
    return wrap(Z)

# Bivariate sine von Mises
def bivariate_sine_vm_log_density(theta1, theta2, mu1, mu2, kappa1, kappa2, lambda_param):
    return (kappa1 * np.cos(theta1 - mu1) +
            kappa2 * np.cos(theta2 - mu2) +
            lambda_param * np.sin(theta1 - mu1) * np.sin(theta2 - mu2))

def generate_bivariate_sine_vm(mu1, mu2, kappa1, kappa2, lambda_param, n_samples, seed=None):
    """Generate via rejection sampling"""
    if seed is not None:
        np.random.seed(seed)

    # Max density for rejection sampling
    grid = np.linspace(-np.pi, np.pi, 100)
    max_log_dens = -np.inf
    for t1 in grid:
        for t2 in grid:
            ld = bivariate_sine_vm_log_density(t1, t2, mu1, mu2, kappa1, kappa2, lambda_param)
            if ld > max_log_dens:
                max_log_dens = ld

    samples = []
    while len(samples) < n_samples:
        theta1 = np.random.uniform(-np.pi, np.pi)
        theta2 = np.random.uniform(-np.pi, np.pi)
        log_d = bivariate_sine_vm_log_density(theta1, theta2, mu1, mu2, kappa1, kappa2, lambda_param)
        log_u = np.log(np.random.uniform())
        if log_u < (log_d - max_log_dens):
            samples.append([theta1, theta2])

    return np.array(samples)

print("EXPERIMENT 1: WRAPPED NORMAL (SIMPLIFIED)", flush=True)
print("-" * 80, flush=True)

q = 5
mu = np.zeros(q)
Sigma = np.eye(q)  # Independent for simplicity
n_reps = 30  # Reduced for speed
F = 500

results_wn = []
n_false_pos = 0

print(f"Testing {n_reps} independent wrapped normal replicates (q={q}, F={F})...", flush=True)

for rep in range(n_reps):
    samples = generate_wrapped_normal(mu, Sigma, F, seed=42 + rep)

    # Fit K=1 and K=2
    mu1_fit, kappa1_fit, loglik1 = fit_product_von_mises_single(samples)
    pi_em, mu_em, kappa_em, loglik2, converged = fit_product_von_mises_mixture_em(
        samples, K=2, n_init=1, max_iter=20, seed=42+rep
    )

    n_params_1 = 2 * q
    n_params_2 = 2 * 2 * q + 1
    bic1 = float(-2 * loglik1 + n_params_1 * np.log(F))
    bic2 = float(-2 * loglik2 + n_params_2 * np.log(F))
    delta_bic = float(bic1 - bic2)

    selects_K2 = bool(delta_bic < 0)
    if selects_K2:
        n_false_pos += 1

    results_wn.append({
        'rep': rep,
        'delta_bic': delta_bic,
        'selects_K2': selects_K2
    })

    if (rep + 1) % 10 == 0:
        print(f"  Progress: {rep+1}/{n_reps}", flush=True)

fpr_wn = n_false_pos / n_reps
ci_wn = wilson_ci(n_false_pos, n_reps)

print(f"\nRESULTS:", flush=True)
print(f"  FPR: {fpr_wn:.1%} [{ci_wn[0]:.1%}, {ci_wn[1]:.1%}] (Wilson 95% CI)", flush=True)
print(f"  False positives: {n_false_pos}/{n_reps}", flush=True)
print()

print("EXPERIMENT 2: BIVARIATE SINE VON MISES (CORRECTED)", flush=True)
print("-" * 80, flush=True)

mu1, mu2 = 0.0, 0.0
kappa1, kappa2 = 5.0, 5.0
sqrt_kappa_product = np.sqrt(kappa1 * kappa2)  # 5.0

# CORRECTED eta parameterization
eta_values = [0.0, 0.5, 0.9]  # Reduced set for speed
n_reps_bsvm = 20  # Reduced for speed
F_bsvm = 500

results_bsvm = []

print(f"Unimodality boundary: |λ| < sqrt(κ₁·κ₂) = {sqrt_kappa_product:.2f}", flush=True)
print()

for eta in eta_values:
    lambda_param = eta * sqrt_kappa_product

    # CORRECTED unimodality check
    is_unimodal = (lambda_param**2 < kappa1 * kappa2)

    print(f"η={eta:.2f}, λ={lambda_param:.2f}, λ²={lambda_param**2:.1f} < κ·κ={kappa1*kappa2:.1f}? {is_unimodal}", flush=True)

    if not is_unimodal:
        print("  BIMODAL - SKIPPING", flush=True)
        continue

    n_false_pos_bsvm = 0

    for rep in range(n_reps_bsvm):
        # Generate bivariate samples
        samples_2d = generate_bivariate_sine_vm(mu1, mu2, kappa1, kappa2, lambda_param, F_bsvm, seed=100+rep)

        # Fit K=1 and K=2
        mu1_fit, kappa1_fit, loglik1 = fit_product_von_mises_single(samples_2d)
        pi_em, mu_em, kappa_em, loglik2, converged = fit_product_von_mises_mixture_em(
            samples_2d, K=2, n_init=1, max_iter=20, seed=100+rep
        )

        n_params_1 = 2 * 2
        n_params_2 = 2 * 2 * 2 + 1
        bic1 = float(-2 * loglik1 + n_params_1 * np.log(F_bsvm))
        bic2 = float(-2 * loglik2 + n_params_2 * np.log(F_bsvm))
        delta_bic = float(bic1 - bic2)

        selects_K2 = bool(delta_bic < 0)
        if selects_K2:
            n_false_pos_bsvm += 1

    fpr_bsvm = n_false_pos_bsvm / n_reps_bsvm
    ci_bsvm = wilson_ci(n_false_pos_bsvm, n_reps_bsvm)

    print(f"  FPR: {fpr_bsvm:.1%} [{ci_bsvm[0]:.1%}, {ci_bsvm[1]:.1%}]", flush=True)

    results_bsvm.append({
        'eta': eta,
        'lambda': lambda_param,
        'fpr': fpr_bsvm,
        'ci_lower': ci_bsvm[0],
        'ci_upper': ci_bsvm[1],
        'n_false_pos': n_false_pos_bsvm,
        'n_reps': n_reps_bsvm
    })

print()
print("=" * 80, flush=True)
print("M3-A2b STREAMLINED VALIDATION COMPLETE", flush=True)
print("=" * 80, flush=True)
print()

print("SUMMARY:", flush=True)
print(f"1. Wrapped Normal (independent): FPR = {fpr_wn:.1%} [{ci_wn[0]:.1%}, {ci_wn[1]:.1%}]", flush=True)
print(f"2. Bivariate Sine von Mises (CORRECTED λ² < κ₁·κ₂):", flush=True)
for r in results_bsvm:
    print(f"   η={r['eta']:.2f}: FPR = {r['fpr']:.1%} [{r['ci_lower']:.1%}, {r['ci_upper']:.1%}]", flush=True)
print()

# Save results
df_wn = pd.DataFrame(results_wn)
df_bsvm = pd.DataFrame(results_bsvm)
df_wn.to_csv('/tmp/m3_a2b_streamlined_wrapped_normal.csv', index=False)
df_bsvm.to_csv('/tmp/m3_a2b_streamlined_bivariate_sine.csv', index=False)

print("Files saved:")
print("  /tmp/m3_a2b_streamlined_wrapped_normal.csv")
print("  /tmp/m3_a2b_streamlined_bivariate_sine.csv")
