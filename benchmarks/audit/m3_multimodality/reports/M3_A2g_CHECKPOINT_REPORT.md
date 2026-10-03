# M3-A2g CHECKPOINT REPORT
**Product von Mises EM Optimizer Convergence Audit**

**Date:** 2026-10-03
**Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c
**Supersedes:** M3-A2f (classification only)
**Status:** INCOMPLETE (initialization and stopping validation invalid)

---

## EXECUTIVE SUMMARY

### Investigation: Iteration Budget vs. Fundamental Failure

M3-A2f reported **97% EM convergence failure** in wrapped-normal stress regime (σ≈1.0). M3-A2g investigated whether this indicated fundamental optimizer failure (O1) or iteration budget exhaustion.

**Key findings:**
1. ⚠ Budget exhaustion confirmed (NOT fundamental failure) - non-converged fits still improving at cutoff
2. ⚠ **POOR convergence rates persist** even with extended budgets:
   - Exact-family (max_iter=500): 111/360 (**30.8%**), varies 8%-83% by F
   - Stress regime (max_iter=500): 102/280 (**36.4%**)
   - Extended stress (max_iter=2000, **unpaired** draws): 56/60 (93.3%)
3. ⚠ Limited monotonicity evidence (4 fixtures, 3384 states including triplicates, 1124 unique comparisons) - 0 violations observed, but **NOT comprehensive proof**
4. ❌ Initialization diagnostic **INVALID** (reseeding on inner starts)
5. ❌ Stopping diagnostic **INVALID** (relative tolerance unimplemented, other issues)

**Optimizer classification:** **O2 - FUNCTIONAL BUT FRAGILE** (provisional)

U — numerical-ceiling local-EM heuristic. O2 — FUNCTIONAL BUT FRAGILE, provisional. Optimization fragility documented; initialization and stopping validation incomplete; statistical robustness unresolved. Not production-ready; no validated production defaults.

**Product-vM detector classification:** Remains **VALID ALGORITHMIC HEURISTIC / STATISTICAL ROBUSTNESS UNRESOLVED**

This checkpoint documents optimization fragility but does NOT resolve optimization robustness. Initialization and stopping validation remain INCOMPLETE due to invalid diagnostics.

---

## 1. MOTIVATION

M3-A2f identified **97% EM convergence failure** in the wrapped-normal σ≈1.0 stress regime (max_iter=100, tol=1e-6).

**Critical question:** Is this fundamental optimizer failure (O1), or iteration budget exhaustion?

**M3-A2g scope:**
- Iteration trajectory analysis
- Convergence rates at extended budgets (max_iter=500, 1000, 2000)
- Initialization sensitivity (n_init=1,3,5,10,20)
- Stopping criterion validation (absolute vs. relative tolerance)
- EM monotonicity verification (limited fixtures)

**Out of scope:**
- Dependent-toroidal validation
- Regularization strategy (U/C/P)
- Bootstrap calibration
- Production deployment

---

## 2. CONVERGENCE AUDIT RESULTS

### 2.1 Exact-Family Convergence (max_iter=500)

**Data:** [results/m3_a2g_exact_family_convergence.csv](../results/m3_a2g_exact_family_convergence.csv)

**Grid:** q∈{1,2,5} × κ∈{1,5,20} × F∈{100,500,1000,2500} × 10 replicates = 360 fits

**Results:**
- **Overall:** 111/360 = **30.8%** converged
- **By F:**
  - F=100: 75/90 (83%)
  - F=500: 20/90 (22%)
  - F=1000: 9/90 (**10%**)
  - F=2500: 7/90 (**8%**)

**Interpretation:** Convergence DEGRADES severely with sample size. At F=2500 (realistic MD sample), only 8% of exact-family fits converge at max_iter=500.

### 2.2 True-K2 Convergence (max_iter=500)

**Data:** [results/m3_a2g_k2_convergence.csv](../results/m3_a2g_k2_convergence.csv)

**Grid:** sep∈{10°,30°,60°,90°} × κ∈{5,20} × occupancy∈{0.5,0.8} × F∈{250,1000} × 10 replicates = 320 fits

**Results:**
- **Overall:** 258/320 = 80.6% converged
- **By separation:**
  - sep=10°: 27/80 (**34%**)
  - sep=30°: 71/80 (89%)
  - sep=60°: 80/80 (100%)
  - sep=90°: 80/80 (100%)

**Interpretation:** Convergence degrades with closer separation (harder problem). At sep=10° (challenging but realistic), only 34% converge.

### 2.3 Stress Regime Convergence (max_iter=500)

**Data:** [results/m3_a2g_stress_regime.csv](../results/m3_a2g_stress_regime.csv)

**Grid:** 7 sigma values × 4 sample sizes × 10 replicates = 280 fits, q=2 fixed
- σ∈{0.6,0.7,0.8,0.9,1.0,1.1,1.2}
- F∈{250,500,1000,2500}

**Results:**
- **Overall:** 102/280 = **36.4%** converged
- **Key examples:**
  - σ=1.0, F=1000: 3/10 (30%)
  - σ=1.0, F=2500: 1/10 (**10%**)
  - σ=0.6, F=2500: 0/10 (**0%**)

**Interpretation:** Wrapped-normal stress regime has POOR convergence even at max_iter=500. This is the regime where Product-vM model is most misspecified.

### 2.4 Extended Budget Experiment (max_iter=2000, UNPAIRED)

**Data:** [results/m3_a2g_extended_maxiter.csv](../results/m3_a2g_extended_maxiter.csv)

**Grid:** σ∈{0.6,0.7,1.0} × F∈{1000,2500} × max_iter∈{500,1000,2000} × 10 replicates = 180 fits

**CRITICAL:** This experiment uses **UNPAIRED draws** - different random seeds per (sigma, F, max_iter, rep) configuration. Seed formula: `seed = 1000 + len(results)`, which increments across ALL prior results. Therefore, max_iter=500, 1000, 2000 for the same (sigma, F, rep) get DIFFERENT random samples.

**Results at max_iter=2000:**
- σ=0.6, F=1000: 10/10 (100%)
- σ=0.6, F=2500: 7/10 (70%)
- σ=0.7, F=1000: 9/10 (90%)
- σ=0.7, F=2500: 10/10 (100%)
- σ=1.0, F=1000: 10/10 (100%)
- σ=1.0, F=2500: 10/10 (100%)
- **Overall:** 56/60 = 93.3%

**Interpretation:** Extended budget helps, but comparison is CONFOUNDED by unpaired draws. The 56/60 (93.3%) at max_iter=2000 cannot be directly compared to the 102/280 (36.4%) at max_iter=500 from the main stress experiment, as they use different random samples.

**Example seed calculation for σ=1.0, F=1000, rep=0:**
- At max_iter=500: seed = 1120 (after 120 prior results in loop order)
- At max_iter=1000: seed = 1130 (after 130 prior results)
- At max_iter=2000: seed = 1140 (after 140 prior results)
- These are DIFFERENT random draws

### 2.5 EM Monotonicity Verification (Limited Scope)

**Script:** [scripts/m3_a2g_monotonicity_check.py](../scripts/m3_a2g_monotonicity_check.py)
**Note:** This is an in-memory verification, NOT derived from m3_a2g_initialization_comparison.csv

**Fixtures tested:**
1. exact_easy: Product-vM K=1, μ=[0, π/4], κ=[5,5], F=1000 (seed=700)
2. wrapped_1.0: Wrapped-normal, σ=1.0, F=1000 (seed=701)
3. k2_clear: Product-vM K=2, sep=90°, κ=20, balanced (seed=702)
4. k2_close: Product-vM K=2, sep=20°, κ=20, unbalanced (π=[0.8,0.2], seed=703)

**Configuration:**
- Uses `initialization='split', n_init=3`: split initialization repeated 3 times per fixture
- Split initialization is DETERMINISTIC, so all 3 runs produce identical trajectories (triplicates)
- 4 fixtures × 3 identical deterministic runs = 12 runs

**Results:**
- **Total iteration states:** 3,384 (including triplicates)
- **Total EM step comparisons:** 3,372 (including triplicates)
- **Unique comparisons:** 1,124 (after removing triplicates from identical trajectories)
- **Monotonicity violations:** 0

**Interpretation:** Zero monotonicity violations suggest EM implementation is not systematically broken. However, this is **LIMITED SCOPE** - only 4 deterministic fixtures, empirical verification only, NOT comprehensive mathematical proof of correctness.

---

## 3. INITIALIZATION AND STOPPING DIAGNOSTICS (BOTH INVALID)

### 3.1 Initialization Sensitivity (INVALID)

**Data:** [results/m3_a2g_ninit_sensitivity.csv](../results/m3_a2g_ninit_sensitivity.csv)
**Script:** [scripts/m3_a2g_ninit_sensitivity.py](../scripts/m3_a2g_ninit_sensitivity.py)

**Test:** Compare n_init∈{1,3,5,10,20} on 4 cases (exact_k1, wrapped_1.0, k2_close, k2_clear)

**CRITICAL BUG:** Script reseeds on every INNER start:
```python
# Line 74-75 in test_ninit_sensitivity():
result = fit_K2_with_full_diagnostics(
    samples, initialization='split' if n_init == 1 else 'random',
    n_init=n_init, max_iter=1000, tol_abs=1e-6, seed=800  # <-- FIXED SEED
)
```

The `seed=800` is passed to the fitter, which reseeds on every inner initialization. This means:
- n_init=1: uses split initialization (deterministic)
- n_init=3: uses 3 random initializations, ALL with same seed=800, producing IDENTICAL starts
- n_init=5: uses 5 random initializations, ALL with same seed=800, producing IDENTICAL starts
- etc.

**Observed result:** wrapped_1.0 shows loglik diff = 2.63 between n_init=1 (split) and n_init=3 (random). This indicates DIFFERENT SOLUTIONS reached, but the comparison is INVALID because of the seed bug.

**Required fix:** The defect is reseeding on every INNER start. Merely changing the outer seed to `800+n_init` or `800+rep` does NOT fix this defect - the inner reseeding must be removed or use independent per-start seeds.

**Status:** **INVALID FOR ROBUSTNESS CLAIMS**

### 3.2 Stopping Criterion Validation (INVALID)

**Data:** [results/m3_a2g_stopping_criterion_validation.csv](../results/m3_a2g_stopping_criterion_validation.csv)
**Script:** [scripts/m3_a2g_stopping_criterion_validation.py](../scripts/m3_a2g_stopping_criterion_validation.py)

**Test:** Compare absolute vs. relative tolerance stopping criteria using continuation test

**CRITICAL BUG:** Relative tolerance branch is unimplemented:
```python
# Line 69-73 in continuation_test():
if tol_rel is not None:
    # Re-run with monitoring to find when relative criterion would fire
    pass  # For simplicity, use tol_abs as proxy
```

When `tol_rel` is specified, the function does not implement relative tolerance logic - it falls through to using absolute criterion with tol=1e-20 (effectively no stopping).

**Configurations tested:**
- abs_1e-6: Uses absolute tolerance (IMPLEMENTED)
- rel_1e-6, rel_1e-7, rel_1e-8: All fall through to tol=1e-20 (UNIMPLEMENTED)

**Additional issues affecting ALL protocols:**
- Approximate stop iteration detection (finds nearest, not exact)
- Incorrect continuation endpoint (uses initial fit's final state, which may not match continuation start)
- Non-converged initial fits invalidate continuation test premise

**Status:** **INVALID FOR VALIDATION** - relative tolerance comparison unexecuted, other methodological issues present

---

## 4. ASSESSMENT

### 4.1 Budget Exhaustion vs. Fundamental Failure

**Finding:** Budget exhaustion confirmed (NOT fundamental EM failure).

**Evidence:**
- Non-converged fits at max_iter=100 still have material log-likelihood improvements in final iterations
- Extended budgets increase convergence rates (though remain poor in some regimes)
- Zero monotonicity violations in limited fixtures suggests EM implementation not systematically broken

**Classification implication:** Avoid O1 (systematic failure). Consider O2 (functional but fragile).

### 4.2 Convergence Fragility

**Finding:** Convergence rates remain POOR even with extended budgets.

**Evidence:**
- Exact-family at max_iter=500: 30.8% overall, only 8% at F=2500
- Stress regime at max_iter=500: 36.4% overall, only 10% at σ=1.0, F=2500
- True-K2 at max_iter=500: 80.6% overall, only 34% at sep=10°

**Interpretation:** Convergence is FRAGILE - heavily dependent on problem difficulty (sample size, separation, model misspecification). No single max_iter budget provides robust convergence across tested configurations.

### 4.3 Monotonicity Evidence (Limited)

**Finding:** Zero monotonicity violations in 1124 unique EM step comparisons.

**Evidence:** 4 deterministic fixtures, each run 3 times with identical split initialization (seed=700).

**Limitation:** LIMITED SCOPE - only 4 fixtures, empirical verification only, NOT comprehensive mathematical proof. Does not cover full range of configurations where convergence fails.

### 4.4 Diagnostics Status

**Initialization:** INVALID (reseeding on inner starts)
**Stopping:** INVALID (relative tolerance unimplemented, other issues)

**Implication:** Initialization and stopping validation INCOMPLETE. Cannot accept O3 or claim robust convergence until diagnostics are fixed and re-executed.

---

## 5. CLASSIFICATION

**Optimizer Classification:** **O2 - FUNCTIONAL BUT FRAGILE** (provisional)

**Rationale:**
- ✅ Budget exhaustion identified (NOT fundamental failure = avoid O1)
- ⚠ Poor convergence rates documented (30.8%, 36.4% at max_iter=500)
- ⚠ Limited monotonicity evidence (suggests not systematically broken, but NOT proof)
- ❌ Initialization validation INCOMPLETE (diagnostic invalid)
- ❌ Stopping validation INCOMPLETE (diagnostic invalid)
- ❌ No validated production defaults

**Product-vM Detector Classification:** Remains **VALID ALGORITHMIC HEURISTIC / STATISTICAL ROBUSTNESS UNRESOLVED**

U — numerical-ceiling local-EM heuristic. O2 — FUNCTIONAL BUT FRAGILE, provisional. Optimization fragility documented; initialization and stopping validation incomplete; statistical robustness unresolved. Not production-ready; no validated production defaults.

---

## 6. LIMITATIONS AND OPEN QUESTIONS

### 6.1 Incomplete Diagnostics

**Initialization:**
- Current diagnostic has reseeding bug - INVALID
- Required: Fix inner reseeding, re-execute with independent per-start seeds
- Required: Systematic comparison of advanced strategies (farthest-point, embedding-based, perturb-K1) on challenging cases

**Stopping:**
- Current diagnostic has unimplemented relative tolerance - INVALID
- Required: Implement actual relative tolerance criterion |Δloglik|/|loglik| < tol_rel
- Required: Address approximate stop iteration, continuation endpoint, non-converged fit issues

### 6.2 Unpaired Extended Budget Experiment

The extended max_iter experiment uses unpaired draws (different seeds per budget). This confounds interpretation:
- Cannot directly compare max_iter=500 convergence rate from main experiments to max_iter=2000 rate from extended experiment
- The 56/60 (93.3%) at max_iter=2000 may reflect easier random draws, not purely budget effect

**Required:** Paired draws experiment with fixed seeds across budgets for controlled comparison.

### 6.3 Limited Monotonicity Scope

Only 4 deterministic fixtures tested, all with identical split initialization. Does not cover:
- Full range of configurations where convergence fails
- Stochastic variation across random initializations
- Edge cases or pathological configurations

**Required:** Extended monotonicity verification on broader fixture set, including failing configurations.

### 6.4 No Production Defaults Established

M3-A2g documents fragility but provides NO validated production recommendations:
- No validated max_iter default (rates vary 8%-100% across configurations)
- No validated tolerance (absolute tol=1e-6 scales poorly, relative unvalidated)
- No validated n_init (comparison invalid due to seed bug)

**Required:** After fixing diagnostics, establish validated defaults with robustness guarantees.

### 6.5 Statistical Robustness Unresolved

M3-A2g addresses optimization fragility only. Statistical validation blockers remain:
- Dependent-toroidal misspecification (MANDATORY before production)
- Regularization strategy (U/C/P comparison incomplete)
- Extended dimensional validation (q=1, q=5)
- Temporal correlation sensitivity
- Calibrated significance testing

---

## 7. NEXT SCIENTIFIC BLOCKERS

**Priority 1: Complete Initialization and Stopping Validation**
- Fix reseeding bug in initialization diagnostic
- Implement relative tolerance in stopping diagnostic
- Re-execute corrected diagnostics
- Only then consider O3 upgrade if results support

**Priority 2: Dependent-Toroidal Validation (MANDATORY)**
- Bivariate sine von Mises with verified unimodality
- Correlated wrapped-normal with verified unimodality
- Measured circular association

**Priority 3: Regularization Strategy Decision**
- Systematic U/C/P comparison
- Literature-grounded penalty selection

**Priority 4: Production Default Validation**
- Establish validated max_iter, tolerance, n_init with robustness guarantees
- Only after fixing diagnostics and statistical validation

**Alternative:** Pivot to model-free multimodality detectors if fragility proves prohibitive.

---

## 8. CHECKPOINT CLOSURE

**M3-A2g Status:** INCOMPLETE

**Completed:**
- ✅ Budget exhaustion confirmed (avoid O1)
- ✅ Convergence fragility documented (30.8%, 36.4% at max_iter=500)
- ✅ Limited monotonicity evidence (0 violations, 1124 unique comparisons)
- ✅ Provisional O2 classification justified

**Incomplete:**
- ❌ Initialization validation (diagnostic invalid)
- ❌ Stopping validation (diagnostic invalid)
- ❌ Production defaults (not established)
- ❌ Statistical robustness (unresolved)

**Next Checkpoint Options:**
- **M3-A2h:** Fix initialization/stopping diagnostics, complete validation
- **M3-A3:** Dependent-toroidal validation (MANDATORY before production)
- **M3-B:** Regularization strategy systematic comparison

**External Review Required:** Yes - provisional O2 pending diagnostic completion and external validation.

---

**Repository Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c
**Report Date:** 2026-10-03
**Status:** Provisional O2, diagnostics incomplete, external review required
