# M3-A2g CHECKPOINT CLOSURE

**Checkpoint:** M3-A2g - Optimizer Convergence Audit
**Date:** 2026-10-03
**Repository Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c
**Status:** **INCOMPLETE** - Initialization and stopping diagnostics invalid

---

## CHECKPOINT QUESTION

**Was M3-A2f's "97% non-convergence" a fundamental optimizer failure (O1) or iteration budget exhaustion?**

**Answer:** Budget exhaustion confirmed (NOT fundamental failure), but convergence remains POOR even with extended budgets.

---

## COMPLETED EXPERIMENTS

1. ✅ Exact-family convergence audit (360 fits, max_iter=500): 111/360 (**30.8%**) converged
2. ✅ True-K2 convergence audit (320 fits, max_iter=500): 258/320 (80.6%) converged
3. ✅ Stress regime audit (280 fits, max_iter=500): 102/280 (**36.4%**) converged
4. ✅ Extended budget test (180 fits, max_iter up to 2000, **UNPAIRED** draws): 56/60 (93.3%) at max_iter=2000
5. ✅ Monotonicity verification (4 fixtures, 3384 states including triplicates): 0 violations
6. ✅ Iteration trajectory analysis (5 cases): Non-converged fits still improving at cutoff
7. ❌ Initialization sensitivity (**INVALID** - reseeding on inner starts)
8. ❌ Stopping criterion validation (**INVALID** - relative tolerance unimplemented)
9. ✅ Duplicate component investigation (8 cases)
10. ✅ Max_iter/tolerance grid (16 configs)

---

## KEY FINDINGS

### Budget Exhaustion Confirmed
- Non-converged fits at max_iter=100 still have material Δloglik improvements
- Extended budgets increase convergence rates (though remain poor in many regimes)
- Zero monotonicity violations in limited fixtures (exact_easy, wrapped_1.0, k2_clear, k2_close)

### Convergence Remains POOR
- **Exact-family (max_iter=500):** 30.8% overall
  - F=100: 83%, F=500: 22%, F=1000: 10%, F=2500: **8%**
- **Stress regime (max_iter=500):** 36.4% overall
  - σ=1.0, F=1000: 30%, σ=1.0, F=2500: **10%**
- **True-K2 (max_iter=500):** 80.6% overall
  - sep=10°: **34%**, sep=30°: 89%, sep=60°: 100%, sep=90°: 100%

### Extended Budget Experiment (UNPAIRED)
- **max_iter=2000:** 56/60 (93.3%) converged
- **CRITICAL:** Uses different random seeds per budget - cannot directly compare to main experiments
- Confounded by draw difficulty variation

### Diagnostics INVALID
- **Initialization:** Fitter reseeds on every inner start - n_init>1 produce identical random starts
- **Stopping:** Relative tolerance branch unimplemented - validation incomplete

---

## CLASSIFICATION

**Optimizer:** **O2 - FUNCTIONAL BUT FRAGILE** (provisional)

U — numerical-ceiling local-EM heuristic. O2 — FUNCTIONAL BUT FRAGILE, provisional. Optimization fragility documented; initialization and stopping validation incomplete; statistical robustness unresolved. Not production-ready; no validated production defaults.

**Rationale:**
- ✅ Budget exhaustion identified (avoid O1 - not fundamental failure)
- ⚠ Poor convergence rates documented (30.8%, 36.4% at max_iter=500)
- ⚠ Limited monotonicity evidence (4 fixtures, NOT comprehensive proof)
- ❌ Initialization validation INCOMPLETE (diagnostic invalid)
- ❌ Stopping validation INCOMPLETE (diagnostic invalid)
- ❌ No validated production defaults

**NOT O3 because:**
- Poor convergence rates (30.8%, 36.4%)
- Invalid diagnostics (initialization, stopping)
- No robust production defaults established
- Fragility not resolved

---

## INCOMPLETE DIRECTIVES

### Initialization Validation
**Status:** INVALID
**Issue:** Script reseeds on every inner start - all n_init>1 produce identical random starts
**Required:** Fix inner reseeding defect, re-execute with independent per-start seeds

### Stopping Validation
**Status:** INVALID
**Issue:** Relative tolerance branch unimplemented - falls through to tol=1e-20
**Required:** Implement actual relative criterion |Δloglik|/|loglik| < tol_rel

### Production Defaults
**Status:** NOT ESTABLISHED
**Required:** After fixing diagnostics, establish validated max_iter, tolerance, n_init with robustness guarantees

### Statistical Robustness
**Status:** UNRESOLVED
**Required:** Dependent-toroidal validation (MANDATORY), U/C/P regularization, extended dimensional validation, calibrated significance

---

## NEXT ACTIONS

1. **Fix initialization and stopping diagnostics** (M3-A2h candidate)
2. **Re-execute corrected diagnostics**
3. **Only then consider O3 upgrade** if results support
4. **Dependent-toroidal validation** (MANDATORY before production)
5. **External review** of provisional O2 classification

---

## CORRECTED CLAIMS HISTORY

**Original erroneous claims (REJECTED by external review):**
- ❌ "O3 - ROBUST"
- ❌ "EM algorithmically correct" (overstated from limited monotonicity)
- ❌ "max_iter ≥ 250-500 resolves convergence in tested configurations"
- ❌ "Optimization robustness PARTIALLY RESOLVED"
- ❌ "All directives satisfied"
- ❌ "COMPLETE"

**Corrected classification:**
- ✅ "O2 - FUNCTIONAL BUT FRAGILE" (provisional)
- ✅ "Optimization fragility documented"
- ✅ "Initialization and stopping validation INCOMPLETE"
- ✅ "Not production-ready; no validated production defaults"

**Numerical corrections:**
- ❌ Claimed: κ grid {1,5,20,50} → ✅ Actual: {1,5,20}
- ❌ Claimed: 5 iteration records → ✅ Actual: 1634 iteration state records
- ❌ Claimed: "split plus two random" → ✅ Actual: "split for n_init=1, all-random for n_init>1"
- ❌ Claimed: "only abs valid" → ✅ Actual: "ALL protocols have issues"

---

**Repository Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c
**Checkpoint Status:** Provisional O2, diagnostics incomplete, external review required
