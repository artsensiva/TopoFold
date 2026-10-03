# M3-A2 vs. M3-A2b: Methodological Corrections

This document explains the critical corrections made in M3-A2b to address flaws in the M3-A2 experimental design.

---

## Summary of Changes

| Aspect | M3-A2 (FLAWED) | M3-A2b (CORRECTED) |
|--------|----------------|-------------------|
| **Generators** | `θᵢ ~ vM(ρ·θᵢ₋₁, κ)` and `θᵢ ~ vM(α·ψ, κ)` | Wrapped normal, bivariate sine von Mises |
| **Validity** | Multiplication not branch-independent on S¹ | Mathematically rigorous, literature-cited |
| **Unimodality** | Assumed from code structure | Analytically verified or proven |
| **Correlation** | Reported generator parameter as "correlation" | Measured circular-circular correlation |
| **Confidence Intervals** | Normal approx: 0/n → [0%, 0%] | Wilson score: 0/n → [0%, 3.6%] |
| **Failure Modes** | Confounded | Oracle BIC separates power/optimization/misspecification |

---

## Flaw 1: Invalid Angular Multiplication

### M3-A2 Generator (FLAWED)

```python
# Chain correlation
theta_i | theta_{i-1} ~ vonMises(rho * theta_{i-1}, kappa)

# Common latent
psi ~ vonMises(0, kappa)
theta_i ~ vonMises(alpha * psi, kappa)
```

### Why This Is Invalid

For angular coordinates on S¹, multiplication by a scalar is **not branch-independent**:

**Example:** Consider `ψ = π/2` and `ψ' = π/2 + 2π`.

These represent the **same direction**:
```
ψ = π/2   (90°)
ψ' = 5π/2 (450° ≡ 90°)
```

But with `α = 0.7`:
```
α·ψ  = 0.7·(π/2)   = 0.35π ≈ 63°
α·ψ' = 0.7·(5π/2)  = 1.75π ≈ 315°
```

These are **different directions** (63° ≠ 315°).

**Consequence:** The operation `α·ψ` is not well-defined as a circular function unless α is an integer.

For non-integer `ρ` or `α`, the generator does not produce a valid distribution on T^q in the mathematical sense.

### M3-A2b Correction

Use **valid toroidal distributions** from the literature:

1. **Wrapped Normal:** Z ~ N(μ, Σ), θ = Z mod 2π (Mardia & Jupp 2000)
2. **Bivariate Sine von Mises:** Explicit joint density (Mardia & Sutton 1978)

Both have rigorously defined densities on T^q.

---

## Flaw 2: Unverified Unimodality

### M3-A2 Approach (FLAWED)

```python
# Generate with one code block
samples = generate_common_latent(n, q, mu, kappa, alpha, seed)

# Assume it's K=1 because we called it "one state"
# No verification
```

**Problem:** Some toroidal distributions can have:
- Unimodal marginals but multimodal joint density
- Parameter regimes where the distribution itself becomes multimodal

If the "null" is actually multimodal, detecting K=2 is CORRECT, not a false positive.

### M3-A2b Correction

**Wrapped Normal:**
- **Analytic guarantee:** Always unimodal for any positive-definite Σ
- **Proof:** Weighted sum of Gaussians centered at μ + 2πm has global mode at μ

**Bivariate Sine von Mises:**
- **Sufficient condition:** |λ| < κ₁·κ₂ ensures unimodality (Mardia & Sutton 1978)
- **Tested:** All λ ∈ {0, 5, 10, 15, 20} satisfy |λ| < 25
- **Enforcement:** Code raises error if condition violated

---

## Flaw 3: Generator Parameter ≠ Correlation

### M3-A2 Reporting (FLAWED)

```
Common latent, alpha=0.7 → 100% FPR
```

**Problem:** `α` is a **generator parameter**, not a circular correlation coefficient.

The relationship between `α` and actual circular dependence is:
1. Non-trivial
2. Depends on κ
3. Not validated

We cannot claim "biologically plausible correlation ≥0.5" based solely on `α≥0.5`.

### M3-A2b Correction

**Step 1:** Generate samples using valid distribution

**Step 2:** Measure actual circular correlation in generated samples:
```python
corr_matrix = measure_toroidal_dependence(samples)
# Uses circular-circular correlation coefficient
# (Jammalamadaka & SenGupta 2001)
```

**Step 3:** Report MEASURED correlation, not generator parameter:
```
Wrapped normal, rho_param=0.8 → measured_corr=0.63 → FPR=...
```

This allows comparison with empirically observed protein torsion correlation.

---

## Flaw 4: Improper Confidence Intervals

### M3-A2 Approach (FLAWED)

Used normal approximation:
```
FPR = k/n
SE = sqrt(FPR * (1-FPR) / n)
CI = [FPR - 1.96*SE, FPR + 1.96*SE]
```

**Failure at boundaries:**
```
k=0, n=30:
  FPR = 0/30 = 0.0
  SE = sqrt(0.0 * 1.0 / 30) = 0.0
  CI = [0.0, 0.0]  ← WRONG
```

This incorrectly implies **zero uncertainty** when observing 0/30.

### M3-A2b Correction

**Wilson Score Interval:**
```python
def wilson_score_interval(k, n, alpha=0.05):
    z = norm.ppf(1 - alpha/2)
    p_hat = k / n
    denominator = 1 + z**2 / n
    center = (p_hat + z**2 / (2*n)) / denominator
    margin = z * sqrt((p_hat * (1 - p_hat) / n + z**2 / (4*n**2))) / denominator
    return (max(0, center - margin), min(1, center + margin))
```

**Correct boundaries:**
```
k=0, n=100:   CI = [0.0%, 3.6%]    ← Retains uncertainty
k=100, n=100: CI = [96.4%, 100.0%] ← Retains uncertainty
```

**Reference:** Wilson, E.B. (1927), JASA, 22, 209-212

---

## Flaw 5: Confounded Failure Modes

### M3-A2 Conclusion (INCOMPLETE)

```
10° separation: ΔBIC > 0, gap_to_oracle = 350 log-units
→ Conclusion: "EM optimization failure"
```

**Problem:** Large oracle gap suggests EM issue, but we didn't verify whether **oracle itself** selects K=2.

If ΔBIC_Oracle > 0, then even perfect optimization would fail (insufficient power).

### M3-A2b Oracle BIC Analysis

For each case, compute:

1. BIC_K1_MLE: Best K=1 fit
2. BIC_K2_EM: Best K=2 EM fit
3. BIC_K2_Oracle: K=2 at TRUE parameters

**Decision tree:**
```
ΔBIC_Oracle = BIC_K2_Oracle - BIC_K1_MLE

if ΔBIC_Oracle < 0:
    # Oracle favors K=2
    if ΔBIC_EM > 0:
        → Optimizer failure (EM not finding oracle)
    else:
        → Success
else:
    → Insufficient power (even oracle favors K=1)
```

This **separates** power from optimization cleanly.

---

## Impact on Detector Classification

### If M3-A2b Shows Low FPR on Valid Nulls

**Conclusion:** M3-A2 rejection was based on invalid generators.

**Classification:** Product von Mises is **rehabilitated** as primary candidate (with EM improvements).

**M3-A3 Direction:** Focus on bootstrap calibration, multiple testing, temporal autocorrelation.

### If M3-A2b Confirms High FPR on Valid Nulls

**Conclusion:** Product von Mises genuinely misspecifies correlated states.

**Classification:** **Rejected** as primary detector.

**M3-A3 Direction:** Evaluate wrapped normal mixtures, multivariate von Mises, toroidal KDE.

---

## Why These Corrections Matter

### Scientific Rigor

M3-A2 would reject a detector based on potentially invalid nulls. This could:
1. Incorrectly eliminate a valid method
2. Waste effort developing unnecessary alternatives

M3-A2b ensures **any rejection is based on mathematically sound evidence**.

### Reproducibility

M3-A2 generators (`ρ * θ`) are not standard toroidal distributions. Another researcher could not:
1. Find them in the literature
2. Reproduce the exact FPR
3. Compare to their own methods

M3-A2b uses **literature-cited distributions** that can be independently verified.

### Biological Interpretation

M3-A2: "FPR=100% at α=0.7"
→ What does α=0.7 mean biologically?
→ How does it compare to real protein data?
→ Unknown.

M3-A2b: "FPR=X% at measured circular correlation=0.63"
→ Can be compared to empirically observed φ-ψ correlation
→ Biological plausibility can be assessed

---

## Summary

M3-A2b is not a minor refinement. It is a **fundamental correction** of the experimental design to ensure:

1. ✓ Mathematically valid null distributions
2. ✓ Verified unimodality (not assumed)
3. ✓ Measured circular correlation (not parameter-based)
4. ✓ Statistically proper confidence intervals
5. ✓ Separated failure modes (power/optimization/misspecification)

**Only after M3-A2b results can we confidently classify Product von Mises detector.**

---

**END OF DOCUMENT**
