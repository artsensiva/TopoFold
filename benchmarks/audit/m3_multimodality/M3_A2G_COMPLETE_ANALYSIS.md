# M3-A2g COMPLETE ANALYSIS
**Product von Mises EM Optimizer Convergence Audit - Full Data Analysis**

**Date:** 2026-10-03
**Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c
**Status:** INCOMPLETE (diagnostics invalid)

---

## SUMMARY

M3-A2g investigated whether M3-A2f's "97% non-convergence" indicated fundamental optimizer failure (O1) or iteration budget exhaustion.

**Finding:** Budget exhaustion confirmed (NOT fundamental failure), but convergence remains **POOR** even with extended budgets (exact: 30.8%, stress: 36.4% at max_iter=500). Limited monotonicity evidence (4 fixtures, 3384 states, 0 violations) suggests EM not systematically broken, but initialization and stopping diagnostics are **INVALID**.

**Classification:** **O2 - FUNCTIONAL BUT FRAGILE** (provisional)

U — numerical-ceiling local-EM heuristic. O2 — FUNCTIONAL BUT FRAGILE, provisional. Optimization fragility documented; initialization and stopping validation incomplete; statistical robustness unresolved. Not production-ready; no validated production defaults.

---

## CONVERGENCE AUDIT DATA (1235 CSV rows)

### Exact-Family (max_iter=500, 360 fits)
**Grid:** q∈{1,2,5} × κ∈{1,5,20} × F∈{100,500,1000,2500} × 10 reps

**Results:** 111/360 (30.8%) converged
- F=100: 75/90 (83%)
- F=500: 20/90 (22%)
- F=1000: 9/90 (10%)
- F=2500: 7/90 (**8%**)

**Interpretation:** POOR convergence, degrades severely with sample size.

### True-K2 (max_iter=500, 320 fits)
**Grid:** sep∈{10°,30°,60°,90°} × κ∈{5,20} × occupancy∈{0.5,0.8} × F∈{250,1000} × 10 reps

**Results:** 258/320 (80.6%) converged
- sep=10°: 27/80 (34%)
- sep=30°: 71/80 (89%)
- sep=60°: 80/80 (100%)
- sep=90°: 80/80 (100%)

**Interpretation:** Convergence degrades with close separation.

### Stress Regime (max_iter=500, 280 fits)
**Grid:** 7 sigma values × 4 sample sizes × 10 reps = 280 fits, q=2 fixed
- σ∈{0.6,0.7,0.8,0.9,1.0,1.1,1.2}, F∈{250,500,1000,2500}

**Results:** 102/280 (36.4%) converged
- σ=1.0, F=1000: 3/10 (30%)
- σ=1.0, F=2500: 1/10 (**10%**)

**Interpretation:** POOR convergence in misspecified regime.

### Extended Budget (max_iter=2000, 180 fits, UNPAIRED)
**Grid:** σ∈{0.6,0.7,1.0} × F∈{1000,2500} × max_iter∈{500,1000,2000} × 10 reps

**Results at max_iter=2000:** 56/60 (93.3%) converged

**CRITICAL:** Uses **UNPAIRED** random draws (different seeds per budget). Cannot directly compare to main experiments.

### Monotonicity (4 fixtures, 3384 states)
**Fixtures:** exact_easy, wrapped_1.0, k2_clear, k2_close
**Configuration:** Uses `initialization='split', n_init=3` - 3 repeated deterministic split starts per fixture (creates triplicates)

**Results:**
- 3384 iteration states (including triplicates)
- 1124 unique comparisons (after removing triplicates)
- 0 violations

**Interpretation:** LIMITED SCOPE evidence - suggests not systematically broken, but NOT comprehensive proof.

### Trajectory Analysis (1634 iteration records, 5 cases)
**Cases:** exact_easy, wrapped_0.9, wrapped_1.0, k2_clear, k2_close

**Finding:** Non-converged fits still improving at max_iter cutoff - confirms budget exhaustion, not fundamental failure.

---

## INVALID DIAGNOSTICS

### Initialization Sensitivity (INVALID)
**Data:** 20 fits (4 cases × 5 n_init values)

**Issue:** Script reseeds on every INNER start - all n_init>1 produce identical random starts.

**Observed:** wrapped_1.0 loglik diff = 2.63 between n_init=1 (split) and n_init=3 (random), but comparison INVALID.

**Status:** **INVALID FOR ROBUSTNESS CLAIMS**

### Stopping Criterion (INVALID)
**Data:** 20 fits (4 configs × 5 cases)

**Issue:** Relative tolerance branch unimplemented - falls through to tol=1e-20.

**Status:** **INVALID FOR VALIDATION**

---

## CLASSIFICATION RATIONALE

**O2 - FUNCTIONAL BUT FRAGILE** (provisional)

**Evidence AGAINST O1 (systematic failure):**
- Budget exhaustion confirmed (not algorithmic failure)
- Zero monotonicity violations in limited scope
- Extended budgets improve convergence (though remain poor)

**Evidence AGAINST O3 (robust):**
- Poor convergence rates (30.8%, 36.4% at max_iter=500)
- Invalid initialization diagnostic
- Invalid stopping diagnostic
- No validated production defaults

**Provisional O2 justified:**
- EM implementation not systematically broken
- Fragility documented
- Diagnostics incomplete
- Not production-ready

---

## REMAINING BLOCKERS

1. **Initialization validation:** Fix reseeding defect, re-execute
2. **Stopping validation:** Implement relative tolerance, re-execute
3. **Production defaults:** Establish validated max_iter, tolerance, n_init
4. **Statistical robustness:** Dependent-toroidal validation (MANDATORY), U/C/P, extended dimensional, calibrated significance

---

**Repository Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c
**Status:** Provisional O2, diagnostics incomplete, external review required
