#!/usr/bin/env python3
"""
M3-A2c FITTING / LIKELIHOOD SANITY AUDIT

Critical issue: 100% FPR on eta=0 (lambda=0) bivariate sine von Mises,
which IS exactly a product von Mises distribution.

This cannot be explained by model misspecification - it indicates
implementation bugs or likelihood degeneracy.
"""

import sys
sys.path.insert(0, '/tmp')

import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.special import i0, i1

from m3_a_diagnostic import (
    wrap, fit_product_von_mises_mixture_em, fit_product_von_mises_single,
    product_von_mises_logpdf
)

print("M3-A2c FITTING / LIKELIHOOD SANITY AUDIT", flush=True)
print("=" * 80, flush=True)
print()

# Wilson CI
def wilson_ci(k, n, alpha=0.05):
    if n == 0:
        return (0.0, 1.0)
    z = norm.ppf(1 - alpha/2)
    p_hat = k / n
    denom = 1 + z**2 / n
    center = (p_hat + z**2 / (2*n)) / denom
    margin = z * np.sqrt((p_hat * (1 - p_hat) / n + z**2 / (4*n**2))) / denom
    return (max(0, center - margin), min(1, center + margin))

#=============================================================================
# 1. EXACT K=1 PRODUCT VON MISES NULLS
#=============================================================================

print("CHECKPOINT 1: EXACT K=1 PRODUCT VON MISES NULLS", flush=True)
print("-" * 80, flush=True)
print()
print("Testing on data generated from EXACTLY the model being fitted.", flush=True)
print("Expected: FPR ~5% (nominal alpha)", flush=True)
print()

def generate_product_von_mises_exact(mu, kappa, n_samples, seed=None):
    """Generate from product von Mises (exact null for the detector)"""
    if seed is not None:
        np.random.seed(seed)

    q = len(mu)
    samples = np.zeros((n_samples, q))

    for dim in range(q):
        # von Mises sampling via rejection
        angles = []
        while len(angles) < n_samples:
            theta_prop = np.random.uniform(-np.pi, np.pi)
            u = np.random.uniform()

            # von Mises density (unnormalized)
            density = np.exp(kappa[dim] * np.cos(theta_prop - mu[dim]))
            max_density = np.exp(kappa[dim])

            if u < density / max_density:
                angles.append(theta_prop)

        samples[:, dim] = np.array(angles[:n_samples])

    return samples

# Test configurations
configs = [
    {'q': 1, 'F': 500, 'kappa': 5.0, 'n_reps': 30},
    {'q': 2, 'F': 500, 'kappa': 5.0, 'n_reps': 30},
    {'q': 5, 'F': 500, 'kappa': 5.0, 'n_reps': 30},
]

results_exact = []

for config in configs:
    q = config['q']
    F = config['F']
    kappa_val = config['kappa']
    n_reps = config['n_reps']

    mu = np.zeros(q)
    kappa = np.ones(q) * kappa_val

    print(f"Testing q={q}, F={F}, kappa={kappa_val:.1f}...", flush=True)

    n_false_pos = 0

    for rep in range(n_reps):
        samples = generate_product_von_mises_exact(mu, kappa, F, seed=1000+rep)

        # Fit K=1
        mu1_fit, kappa1_fit, loglik1 = fit_product_von_mises_single(samples)

        # Fit K=2
        pi_em, mu_em, kappa_em, loglik2, converged = fit_product_von_mises_mixture_em(
            samples, K=2, n_init=1, max_iter=50, seed=1000+rep
        )

        # BIC
        n_params_1 = 2 * q
        n_params_2 = 4 * q + 1
        bic1 = -2 * loglik1 + n_params_1 * np.log(F)
        bic2 = -2 * loglik2 + n_params_2 * np.log(F)
        delta_bic = bic1 - bic2

        selects_K2 = (delta_bic < 0)
        if selects_K2:
            n_false_pos += 1

        results_exact.append({
            'q': q,
            'F': F,
            'kappa_true': kappa_val,
            'rep': rep,
            'loglik1': loglik1,
            'loglik2': loglik2,
            'bic1': bic1,
            'bic2': bic2,
            'delta_bic': delta_bic,
            'selects_K2': selects_K2,
            'converged': converged
        })

    fpr = n_false_pos / n_reps
    ci = wilson_ci(n_false_pos, n_reps)

    print(f"  FPR: {fpr:.1%} [{ci[0]:.1%}, {ci[1]:.1%}] ({n_false_pos}/{n_reps})", flush=True)
    print()

#=============================================================================
# 2. INDEPENDENT K=1 MLE VERIFICATION
#=============================================================================

print("CHECKPOINT 2: INDEPENDENT K=1 MLE VERIFICATION", flush=True)
print("-" * 80, flush=True)
print()

def compute_k1_mle_independent(samples):
    """Independent computation of K=1 MLE for verification"""
    n, q = samples.shape

    mu_indep = np.zeros(q)
    kappa_indep = np.zeros(q)

    for dim in range(q):
        theta = samples[:, dim]

        # Circular mean
        S = np.sum(np.sin(theta))
        C = np.sum(np.cos(theta))
        mu_indep[dim] = np.arctan2(S, C)

        # Mean resultant length
        R_bar = np.sqrt(S**2 + C**2) / n

        # Kappa via A^-1(R_bar) approximation
        if R_bar < 0.53:
            kappa_indep[dim] = 2 * R_bar + R_bar**3 + 5 * R_bar**5 / 6
        elif R_bar < 0.85:
            kappa_indep[dim] = -0.4 + 1.39 * R_bar + 0.43 / (1 - R_bar)
        else:
            kappa_indep[dim] = 1 / (R_bar**3 - 4 * R_bar**2 + 3 * R_bar)

    # Compute log-likelihood independently
    loglik_indep = 0.0
    for i in range(n):
        for dim in range(q):
            loglik_indep += (kappa_indep[dim] * np.cos(samples[i, dim] - mu_indep[dim])
                           - np.log(2 * np.pi) - np.log(i0(kappa_indep[dim])))

    return mu_indep, kappa_indep, loglik_indep

# Verify on one example
q = 5
F = 500
mu_true = np.zeros(q)
kappa_true = np.ones(q) * 5.0

samples_test = generate_product_von_mises_exact(mu_true, kappa_true, F, seed=42)

mu_code, kappa_code, loglik_code = fit_product_von_mises_single(samples_test)
mu_indep, kappa_indep, loglik_indep = compute_k1_mle_independent(samples_test)

print(f"Test dataset: q={q}, F={F}, kappa_true={kappa_true[0]:.1f}", flush=True)
print()
print("MLE Comparison (Code vs Independent):", flush=True)
print(f"  Max |mu diff|:    {np.max(np.abs(wrap(mu_code - mu_indep))):.6f}", flush=True)
print(f"  Max |kappa diff|: {np.max(np.abs(kappa_code - kappa_indep)):.6f}", flush=True)
print(f"  LogLik code:      {loglik_code:.4f}", flush=True)
print(f"  LogLik indep:     {loglik_indep:.4f}", flush=True)
print(f"  LogLik diff:      {abs(loglik_code - loglik_indep):.6f}", flush=True)
print()

#=============================================================================
# 3. BIC PARAMETER COUNTING VERIFICATION
#=============================================================================

print("CHECKPOINT 3: BIC PARAMETER COUNTING", flush=True)
print("-" * 80, flush=True)
print()

print("K=1 Product von Mises:", flush=True)
print("  Parameters: q means + q kappas", flush=True)
print("  p1 = 2*q", flush=True)
print(f"  Example (q=5): p1 = {2*5}", flush=True)
print()

print("K=2 Product von Mises Mixture:", flush=True)
print("  Parameters: 1 weight + 2*q means + 2*q kappas", flush=True)
print("  p2 = 1 + 4*q", flush=True)
print(f"  Example (q=5): p2 = {1 + 4*5}", flush=True)
print()

print("BIC Formula:", flush=True)
print("  BIC = -2*logL + p*log(n)", flush=True)
print("  where n = F (number of independent frames)", flush=True)
print("  NOT n = F*q (number of scalar coordinates)", flush=True)
print()

print("Code verification:")
q_test = 5
F_test = 500
p1_expected = 2 * q_test
p2_expected = 1 + 4 * q_test

print(f"  Expected p1 = {p1_expected}")
print(f"  Expected p2 = {p2_expected}")
print(f"  Sample size n = {F_test}")
print(f"  Penalty difference = (p2-p1)*log(n) = {(p2_expected - p1_expected) * np.log(F_test):.2f}")
print()

#=============================================================================
# 4. NESTED MODEL IDENTITY TEST
#=============================================================================

print("CHECKPOINT 4: NESTED MODEL IDENTITY TEST", flush=True)
print("-" * 80, flush=True)
print()

print("K=2 contains K=1 as special case: duplicate components with equal weights", flush=True)
print()

# Generate K=1 data
q = 5
F = 500
mu_true = np.zeros(q)
kappa_true = np.ones(q) * 5.0
samples_nested = generate_product_von_mises_exact(mu_true, kappa_true, F, seed=99)

# Fit K=1
mu1, kappa1, loglik1 = fit_product_von_mises_single(samples_nested)

# Construct K=2 with duplicate components (no optimization)
pi_dup = np.array([0.5, 0.5])
mu_dup = np.stack([mu1, mu1])
kappa_dup = np.stack([kappa1, kappa1])

# Compute K=2 log-likelihood manually
loglik2_dup = 0.0
for i in range(F):
    log_mix = -np.inf
    for k in range(2):
        log_comp = 0.0
        for dim in range(q):
            log_comp += (kappa_dup[k, dim] * np.cos(samples_nested[i, dim] - mu_dup[k, dim])
                        - np.log(2 * np.pi) - np.log(i0(kappa_dup[k, dim])))
        log_mix = np.logaddexp(log_mix, np.log(pi_dup[k]) + log_comp)
    loglik2_dup += log_mix

print(f"K=1 loglik:                    {loglik1:.6f}", flush=True)
print(f"K=2 loglik (duplicate):        {loglik2_dup:.6f}", flush=True)
print(f"Difference:                    {abs(loglik1 - loglik2_dup):.6f}", flush=True)
print()

p1 = 2 * q
p2 = 1 + 4 * q
bic1_nested = -2 * loglik1 + p1 * np.log(F)
bic2_dup = -2 * loglik2_dup + p2 * np.log(F)
delta_bic_dup = bic1_nested - bic2_dup

print(f"BIC(K=1):                      {bic1_nested:.2f}", flush=True)
print(f"BIC(K=2 duplicate):           {bic2_dup:.2f}", flush=True)
print(f"Δ BIC (K=1 - K=2):            {delta_bic_dup:.2f}", flush=True)
print(f"Expected: (p2-p1)*log(F) =    {(p2-p1)*np.log(F):.2f}", flush=True)
print()

if delta_bic_dup > 0:
    print("✓ PASS: BIC correctly prefers K=1 over duplicate K=2", flush=True)
else:
    print("✗ FAIL: BIC incorrectly prefers duplicate K=2", flush=True)
print()

#=============================================================================
# 5. FALSE-POSITIVE PARAMETER DIAGNOSTICS
#=============================================================================

print("CHECKPOINT 5: FALSE-POSITIVE PARAMETER DIAGNOSTICS", flush=True)
print("-" * 80, flush=True)
print()

print("Analyzing fitted parameters for exact K=1 false positives...", flush=True)
print()

q = 5
F = 500
mu_true = np.zeros(q)
kappa_true = np.ones(q) * 5.0

fp_diagnostics = []

for rep in range(10):
    samples_fp = generate_product_von_mises_exact(mu_true, kappa_true, F, seed=2000+rep)

    # Fit K=1 and K=2
    mu1, kappa1, loglik1 = fit_product_von_mises_single(samples_fp)
    pi_em, mu_em, kappa_em, loglik2, converged = fit_product_von_mises_mixture_em(
        samples_fp, K=2, n_init=3, max_iter=100, seed=2000+rep
    )

    # BIC
    p1 = 2 * q
    p2 = 1 + 4 * q
    bic1 = -2 * loglik1 + p1 * np.log(F)
    bic2 = -2 * loglik2 + p2 * np.log(F)
    delta_bic = bic1 - bic2

    selects_K2 = (delta_bic < 0)

    if selects_K2:
        # Diagnostic info
        min_weight = np.min(pi_em)
        max_kappa = np.max(kappa_em)

        fp_diagnostics.append({
            'rep': rep,
            'loglik1': loglik1,
            'loglik2': loglik2,
            'delta_bic': delta_bic,
            'weight1': pi_em[0],
            'weight2': pi_em[1],
            'min_weight': min_weight,
            'max_kappa': max_kappa,
            'kappa_mean_comp1': np.mean(kappa_em[0]),
            'kappa_mean_comp2': np.mean(kappa_em[1]),
            'converged': converged
        })

if fp_diagnostics:
    df_fp = pd.DataFrame(fp_diagnostics)
    print(f"Found {len(fp_diagnostics)} false positives out of 10 replicates", flush=True)
    print()
    print("False-Positive Component Diagnostics:", flush=True)
    print(df_fp[['rep', 'weight1', 'weight2', 'min_weight', 'max_kappa', 'delta_bic']].to_string(index=False))
    print()

    # Check for likelihood collapse pattern
    tiny_weight_count = np.sum(df_fp['min_weight'] < 0.05)
    large_kappa_count = np.sum(df_fp['max_kappa'] > 100)

    print(f"Likelihood collapse indicators:", flush=True)
    print(f"  Cases with minority weight < 5%:  {tiny_weight_count}/{len(fp_diagnostics)}", flush=True)
    print(f"  Cases with max kappa > 100:       {large_kappa_count}/{len(fp_diagnostics)}", flush=True)
    print()
else:
    print("No false positives found in 10 exact K=1 replicates", flush=True)
    print()

#=============================================================================
# 10. WILSON CI VALIDATION
#=============================================================================

print("CHECKPOINT 10: WILSON CI VALIDATION", flush=True)
print("-" * 80, flush=True)
print()

test_cases = [
    (0, 30),
    (30, 30),
    (0, 100),
    (100, 100),
]

print("Wilson score 95% confidence intervals:", flush=True)
print()
print("  x    n     lower    upper")
print("  " + "-" * 32)

for x, n in test_cases:
    ci = wilson_ci(x, n, alpha=0.05)
    print(f"  {x:3d}  {n:3d}   {ci[0]:6.1%}   {ci[1]:6.1%}")

print()

#=============================================================================
# SUMMARY
#=============================================================================

print("=" * 80, flush=True)
print("M3-A2c SANITY AUDIT COMPLETE", flush=True)
print("=" * 80, flush=True)
print()

# Save results
df_exact = pd.DataFrame(results_exact)
df_exact.to_csv('/tmp/m3_a2c_exact_null_results.csv', index=False)

print("Files saved:")
print("  /tmp/m3_a2c_exact_null_results.csv")
print()

print("CRITICAL FINDINGS:")
print()

# Check exact null FPR
for config in configs:
    q = config['q']
    subset = df_exact[df_exact['q'] == q]
    fpr = subset['selects_K2'].mean()
    n_fp = subset['selects_K2'].sum()
    n_total = len(subset)
    ci = wilson_ci(int(n_fp), n_total)

    print(f"  q={q}: FPR = {fpr:.1%} [{ci[0]:.1%}, {ci[1]:.1%}] ({int(n_fp)}/{n_total})")

print()
