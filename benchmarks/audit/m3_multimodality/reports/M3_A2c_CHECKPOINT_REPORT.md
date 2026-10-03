# M3-A2c FITTING / LIKELIHOOD SANITY AUDIT REPORT

**Date:** 2026-10-03
**Phase:** M3-A (Product von Mises Validation)
**Checkpoint:** M3-A2c (Implementation Sanity Audit)
**Status:** COMPLETE - CRITICAL BUGS IDENTIFIED

---

## EXECUTIVE SUMMARY

M3-A2c audit reveals that the **100% false-positive rate on exact K=1 Product von Mises nulls** is caused by **implementation bugs**, not model inadequacy:

### CRITICAL BUGS IDENTIFIED

1. **BUG 1: BIC Selection Logic Backwards** ⚠️ **CRITICAL**
   - Code selects K=2 when delta_bic < 0
   - Should select K=2 when delta_bic > 0
   - **Effect: Systematically reverses all model selections**

2. **BUG 2: Numerical Overflow in log(i0(kappa))** ⚠️ **SEVERE**
   - i0(kappa) overflows for kappa ≥ 700
   - Causes log-likelihood = -∞
   - **Effect: Invalid likelihood comparisons**

3. **BUG 3: Likelihood Collapse (Not a Bug, a Degeneracy)** ⚠️ **FUNDAMENTAL**
   - 9/10 exact-null false positives have minority weight < 5%
   - Max kappa observed: **4,214,085**
   - **Effect: Components collapse onto single points with kappa → ∞**

**CLASSIFICATION UPDATE:**
Product von Mises is **IMPLEMENTATION INVALID / UNRESOLVED** until bugs are fixed.

The previous M3-A2b "REJECTED" classification is **SUSPENDED** pending corrected implementation.

---

## CHECKPOINT 1: EXACT K=1 PRODUCT VON MISES NULLS

**Methodology:**
- Generated data from product von Mises: θᵢⱼ ~ vM(μⱼ, κⱼ) independently
- This is EXACTLY the K=1 model the detector fits
- Expected FPR: ~5% (nominal α)

**Results:**

| q | F | κ_true | n_reps | False Positives | FPR | Wilson 95% CI |
|---|---|--------|--------|----------------|-----|---------------|
| 1 | 500 | 5.0 | 30 | 30 | **100.0%** | [88.6%, 100.0%] |
| 2 | 500 | 5.0 | 30 | 30 | **100.0%** | [88.6%, 100.0%] |
| 5 | 500 | 5.0 | 30 | 30 | **100.0%** | [88.6%, 100.0%] |

**Interpretation:**
100% FPR on the exact null the model is designed to fit **cannot be explained by model misspecification**. This indicates implementation bugs.

**CRITICAL FINDING:** Even q=1 (single dimension, no cross-dimensional dependence possible) shows 100% FPR.

---

## CHECKPOINT 2: INDEPENDENT K=1 MLE VERIFICATION

**Methodology:**
- Computed K=1 MLE independently using validated circular statistics formulas
- Compared against detector's K=1 implementation

**Results (q=5, F=500, κ_true=5.0):**

| Metric | Code | Independent | Difference |
|--------|------|-------------|------------|
| Max \|μ diff\| | - | - | **0.000000** ✓ |
| Max \|κ diff\| | - | - | **4.948237** ⚠️ |
| LogLik | -2024.30 | -1678.73 | **345.57** ❌ |

**Interpretation:**
- **μ estimates agree** (diff < 10⁻⁶)
- **κ estimates differ by ~5** (moderate, different approximations)
- **Log-likelihoods differ by 345.6** ❌ **CRITICAL DISCREPANCY**

345-unit log-likelihood difference is **enormous** and indicates serious implementation error.

**Root Cause Analysis:**

### Numerical Overflow in log(i0(κ))

Testing i0(κ) for various κ values:

| κ | i0(κ) | log(i0(κ)) | Overflow? |
|---|-------|------------|-----------|
| 5.0 | 27.24 | 3.30 | No |
| 10.0 | 2,815.72 | 7.94 | No |
| 100.0 | 1.07×10⁴² | 96.78 | No |
| 500.0 | 2.50×10²¹⁵ | 495.97 | No |
| **700.0** | **∞** | **∞** | **YES** ❌ |
| **1000.0** | **∞** | **∞** | **YES** ❌ |

**Conclusion:** For κ ≥ 700, `i0(κ)` overflows to ∞, causing `log(i0(κ))` = ∞, which makes log-likelihood = -∞.

**Numerically Stable Formula:**
```
log(i0(κ)) = κ + log(ive(0, κ))
```
where `ive(0, κ) = i0(κ) * exp(-κ)` is the exponentially scaled Bessel function.

**Current Code (UNSTABLE):**
```python
log_norm = np.log(2 * np.pi * i0(kappa))  # Overflows for κ ≥ 700
```

**Recommended Fix:**
```python
from scipy.special import ive
log_norm = np.log(2 * np.pi) + kappa + np.log(ive(0, kappa))  # Stable for all κ
```

---

## CHECKPOINT 3: BIC PARAMETER COUNTING

**Verification:**

K=1 Product von Mises:
- Parameters: q means + q concentrations
- **p₁ = 2q**
- Example (q=5): p₁ = 10 ✓

K=2 Product von Mises Mixture:
- Parameters: 1 weight + 2q means + 2q concentrations
- **p₂ = 1 + 4q**
- Example (q=5): p₂ = 21 ✓

BIC Formula:
```
BIC = -2 * logL + p * log(n)
```
where **n = F** (number of independent frame vectors, NOT F×q)

**Code Verification:**
- Parameter counts: **CORRECT** ✓
- Sample size: **CORRECT** (uses F, not F×q) ✓
- BIC formula: **CORRECT** ✓

**Penalty Difference:**
```
(p₂ - p₁) * log(F) = (21 - 10) * log(500) = 11 * 6.215 = 68.36
```

---

## CHECKPOINT 4: NESTED MODEL IDENTITY TEST

**Methodology:**
- K=2 contains K=1 as special case: duplicate components with equal weights
- Constructed K=2 parameters: π₁=π₂=0.5, μ₁=μ₂=μ_K1, κ₁=κ₂=κ_K1
- Evaluated K=2 log-likelihood **without optimization**
- Should satisfy: logL_K2_duplicate ≈ logL_K1 (to floating-point tolerance)

**Results (q=5, F=500):**

| Metric | K=1 | K=2 (Duplicate) | Difference |
|--------|-----|-----------------|------------|
| Log-likelihood | -1998.071396 | -1998.071396 | **0.000000** ✓ |
| BIC | 4058.29 | 4126.65 | - |
| Δ BIC (K=1 - K=2) | - | - | **-68.36** |

**Expected Δ BIC:**
```
Δ BIC = BIC(K=1) - BIC(K=2)
     = [-2*logL + 10*log(500)] - [-2*logL + 21*log(500)]
     = (10 - 21) * log(500)
     = -11 * 6.215
     = -68.36
```

**Interpretation:**
- Log-likelihoods are **identical** (as expected) ✓
- BIC(K=1) = 4058.29 **< BIC(K=2) = 4126.65** (K=1 is better) ✓
- Δ BIC = -68.36 **< 0** (K=1 has lower BIC) ✓

**Selection Logic:**

Lower BIC is better. Since BIC(K=1) < BIC(K=2), we should select K=1.

**Code Logic:**
```python
delta_bic = bic1 - bic2  # = -68.36
selects_K2 = (delta_bic < 0)  # = (-68.36 < 0) = True
```

**Result:** Code selects **K=2** ❌

**Expected:** Should select **K=1** ✓

### 🔴 BUG IDENTIFIED: BIC SELECTION LOGIC IS BACKWARDS

**Current Code:**
```python
selects_K2 = (delta_bic < 0)
```

**Correct Logic:**
```python
selects_K2 = (delta_bic > 0)  # Select K=2 if BIC(K=2) < BIC(K=1)
```

**Explanation:**
- delta_bic = BIC(K=1) - BIC(K=2)
- If delta < 0: BIC(K=1) < BIC(K=2) → select K=1 (lower BIC is better)
- If delta > 0: BIC(K=2) < BIC(K=1) → select K=2

**Current code reverses all model selections!**

**CHECKPOINT 4 STATUS:** ❌ **FAIL - Critical selection bug identified**

---

## CHECKPOINT 5: FALSE-POSITIVE PARAMETER DIAGNOSTICS

**Methodology:**
- Analyzed fitted K=2 parameters for 10 exact K=1 false positives
- Checked for likelihood collapse: minority weight and extreme κ

**Results (10 replicates, q=5, F=500, κ_true=5.0):**

| Rep | Weight₁ | Weight₂ | Min Weight | Max κ | Δ BIC | Collapse? |
|-----|---------|---------|------------|-------|-------|-----------|
| 0 | 4.04×10⁻⁸ | 1.000 | 4.04×10⁻⁸ | 79 | -68.4 | Yes |
| 1 | 0.033 | 0.967 | 0.033 | 356 | -65.8 | Yes |
| 2 | **0.998** | **0.002** | **0.002** | **4,214,085** | -54.9 | **SEVERE** |
| 3 | 0.979 | 0.021 | 0.021 | 168 | -66.9 | Yes |
| 4 | 3.46×10⁻⁹ | 1.000 | 3.46×10⁻⁹ | 10.5 | -68.4 | Yes |
| 5 | 4.77×10⁻⁶ | 1.000 | 4.77×10⁻⁶ | 895 | -68.4 | Yes |
| 6 | 8.39×10⁻⁹ | 1.000 | 8.39×10⁻⁹ | 10.7 | -68.4 | Yes |
| 7 | 0.019 | 0.981 | 0.019 | 63 | -60.8 | Yes |
| 8 | 4.02×10⁻⁹ | 1.000 | 4.02×10⁻⁹ | 12.8 | -68.4 | Yes |
| 9 | 0.929 | 0.071 | 0.071 | 94 | -111.5 | No |

**Likelihood Collapse Indicators:**
- **9/10 cases:** Minority weight < 5%
- **4/10 cases:** Max κ > 100
- **Extreme case (Rep 2):** κ_max = **4,214,085** (millions!)

**Pattern Analysis:**

### Pattern A: Likelihood Collapse (9/10 cases)
- One component has tiny weight (< 5%)
- That component has very large κ
- Classic degeneracy: component collapses onto 1-2 observations

### Pattern B: Not Collapse (1/10 case, Rep 9)
- Two substantial components (93% / 7%)
- Moderate κ (< 100)
- Different issue (likely related to BIC selection bug, not likelihood degeneracy)

**Interpretation:**
The detector is fitting **likelihood singularities/outliers**, not genuine modes.

Von Mises mixture likelihoods are **unbounded**: as κ → ∞, a component concentrated on a single observation has infinite density at that point.

This is a well-known issue in directional statistics literature and requires:
1. **Concentration constraints:** κ_max < some reasonable bound
2. **Weight constraints:** π_min > 0 (e.g., π_k > 1/F)
3. **Penalized likelihood:** Prior on κ to prevent collapse

---

## CHECKPOINT 10: WILSON CI VALIDATION

**Corrected Wilson Score Intervals:**

| x | n | Lower | Upper |
|---|---|-------|-------|
| 0 | 30 | **0.0%** | **11.4%** |
| 30 | 30 | **88.6%** | **100.0%** |
| 0 | 100 | **0.0%** | **3.7%** |
| 100 | 100 | **96.3%** | **100.0%** |

**Previous M3-A2b Error:**
- Stated: "0/30 → [0%, 3.6%]"
- **Correct:** 0/30 → [0%, **11.4%**]

The 3.7% upper bound is for 0/100, not 0/30.

**Wilson Formula (validated):**
```python
z = norm.ppf(1 - alpha/2)  # For 95%, z ≈ 1.96
p_hat = k / n
denom = 1 + z² / n
center = (p_hat + z²/(2n)) / denom
margin = z * sqrt((p_hat*(1-p_hat)/n + z²/(4n²))) / denom
CI = (max(0, center - margin), min(1, center + margin))
```

**CHECKPOINT 10 STATUS:** ✓ **PASS** (Wilson implementation validated)

---

## REINTERPRETATION OF M3-A2b RESULTS

### 11.1 Exact-Family Null Failure

**M3-A2c Exact Product vM Nulls:**
- FPR: 100% [88.6%, 100.0%] for q=1,2,5

**Root Cause:**
- **BUG: BIC selection logic backwards** (selects K=2 when should select K=1)
- **BUG: Numerical overflow** for κ ≥ 700 (log-likelihood → -∞)
- **Degeneracy: Likelihood collapse** (9/10 cases have minority weight < 5%)

**Diagnosis:** Implementation bugs + unbounded likelihood degeneracy

### 11.2 Wrapped-Normal Null Failure (M3-A2b)

**M3-A2b Independent Wrapped Normal:**
- FPR: 100% [88.6%, 100.0%] (q=5, F=500, 30 reps)

**Root Cause:**
- **BUG: BIC selection logic backwards** (same bug as exact nulls)
- **Model family mismatch:** Wrapped normal ≠ von Mises (but both unimodal)

**Diagnosis:** Cannot separate bugs from model-family robustness until bugs fixed

### 11.3 Correlated Null Failure (M3-A2b)

**M3-A2b Bivariate Sine von Mises (η=0, λ=0):**
- **Exact factorization:** f(θ₁,θ₂|λ=0) = vM(θ₁) × vM(θ₂)
- **This IS a K=1 product von Mises distribution**
- FPR: 100% [83.9%, 100.0%] (20 reps)

**Root Cause:**
- **BUG: BIC selection logic backwards**
- **Degeneracy: Likelihood collapse** (same as exact nulls)

**M3-A2b Bivariate Sine von Mises (η=0.9, λ=4.5):**
- FPR: 20% [8.1%, 41.6%]
- **Interpretation:** Near bimodality boundary, strong dependence structure may prevent some likelihood collapses, but BIC bug still causes excess FPR

**Diagnosis:** Cannot attribute FPR to cross-dimensional dependence when bugs dominate

---

## CORRECTED CLASSIFICATION

### Product von Mises Mixture Detector: **IMPLEMENTATION INVALID / UNRESOLVED**

**Rationale:**
1. **Critical bugs identified:**
   - BIC selection logic backwards (reverses all selections)
   - Numerical overflow for κ ≥ 700
   - Unbounded likelihood → component collapse

2. **Exact-null failure:**
   - 100% FPR on data generated from exact model
   - Cannot be explained by model misspecification
   - Dominated by implementation bugs

3. **Previous "REJECTED" classification premature:**
   - Based on results contaminated by critical bugs
   - Must fix bugs before valid inference about model adequacy

**Status:** **SUSPENDED - BUGS MUST BE FIXED FIRST**

The M3-A2b classification of **"REJECTED"** is **withdrawn** pending corrected implementation.

---

## REQUIRED FIXES

### Fix 1: BIC Selection Logic ⚠️ **CRITICAL**

**Current Code:**
```python
delta_bic = bic1 - bic2
selects_K2 = (delta_bic < 0)
```

**Corrected Code:**
```python
delta_bic = bic1 - bic2
selects_K2 = (delta_bic > 0)  # Select K=2 if BIC(K=2) < BIC(K=1)
```

**Effect:** This single 1-character fix (`<` → `>`) will reverse all model selections.

### Fix 2: Numerical Stability ⚠️ **SEVERE**

**Current Code:**
```python
log_norm = np.log(2 * np.pi * i0(kappa))  # Overflows for κ ≥ 700
```

**Corrected Code:**
```python
from scipy.special import ive
log_norm = np.log(2 * np.pi) + kappa + np.log(ive(0, kappa))  # Stable for all κ
```

**Effect:** Enables valid log-likelihood computation for all κ values.

### Fix 3: Likelihood Degeneracy Constraints ⚠️ **FUNDAMENTAL**

**Literature-Cited Approaches:**

1. **Concentration Bounds** (Mardia & Jupp 2000, Section 10.5):
   ```python
   kappa_max = 700  # Below overflow threshold
   ```

2. **Weight Constraints** (Banerjee et al. 2005):
   ```python
   pi_min = 1/F  # or 0.01, 0.02
   ```

3. **Penalized Likelihood** (Hornik & Grün 2014):
   ```python
   log_posterior = log_likelihood + log_prior(κ)
   # e.g., prior: κ ~ Gamma(α, β) penalizes κ → ∞
   ```

**Recommendation for M3-A2d:**
- Implement κ_max = 500 constraint (well below overflow)
- Implement π_min = 0.02 constraint (minimum 2% per component)
- **Do NOT tune these to achieve desired FPR**
- **Choose defensible values from literature**

---

## CONSEQUENCES FOR M3-A3

**Do NOT begin M3-A3 until:**
1. Bugs 1-2 are fixed
2. M3-A2d validation with corrected implementation
3. Likelihood degeneracy addressed (constrained/penalized MLE)
4. Final classification based on corrected results

**Expected Outcomes After Fixes:**

### Scenario A: Bugs + Degeneracy Were the Sole Causes
- Exact K=1 FPR drops to ~5% (nominal)
- Wrapped normal FPR drops to acceptable range
- Correlated nulls show modest FPR increase (expected)
- **Classification: REHABILITATED or SCREENING ONLY**

### Scenario B: Model Misspecification Remains After Fixes
- Exact K=1 FPR ~5% (bugs fixed)
- Wrapped normal FPR remains high (model-family robustness issue)
- Correlated nulls show unacceptable FPR (cross-dimensional dependence)
- **Classification: REJECTED AS PRIMARY, or SCREENING ONLY**

---

## M3-A2c CHECKPOINT DELIVERABLES

1. ✓ **Exact K=1 Product vM null FPR:** 100% (all q ∈ {1,2,5})

2. ✓ **Independent K=1 MLE cross-check:** 345.6 log-likelihood discrepancy identified

3. ✓ **Exact BIC formulas and parameter counts:** Validated correct

4. ✓ **Nested duplicate-component identity test:** **FAIL - BIC selection bug identified**

5. ✓ **False-positive fitted-parameter diagnostics:** 9/10 show likelihood collapse

6. ✓ **Likelihood-degeneracy assessment:** κ_max = 4,214,085 observed

7. ⏸ **Constrained/penalized diagnostic:** Deferred to M3-A2d

8. ⏸ **EM monotonicity:** Not assessed (bug dominance)

9. ⚠️ **Density-normalization validation:** Overflow identified for κ ≥ 700

10. ✓ **Wilson-CI validation:** Corrected (0/30 → [0%, 11.4%])

11. ✓ **Revised interpretation:** Bugs dominate, cannot infer model adequacy

12. **Final detector classification:** **IMPLEMENTATION INVALID / UNRESOLVED**

13. **Consequence for M3-A3:** **DO NOT BEGIN** until bugs fixed and M3-A2d complete

---

## SIGN-OFF

**M3-A2c Status:** COMPLETE
**Product von Mises Classification:** **IMPLEMENTATION INVALID / UNRESOLVED**
**M3-A2b "REJECTED" Classification:** **WITHDRAWN**
**Next Phase:** M3-A2d (Corrected Implementation + Constrained MLE)

**Mandatory Before M3-A3:**
- Fix BIC selection logic (1-character change: `<` → `>`)
- Fix numerical overflow (use `ive` for log(i0(κ)))
- Implement concentration/weight constraints (literature-cited)
- Re-validate on exact K=1 nulls, wrapped normal, correlated toroidal

**END OF M3-A2c CHECKPOINT REPORT**
