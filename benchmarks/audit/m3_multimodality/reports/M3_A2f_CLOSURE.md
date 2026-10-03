# M3-A2f CHECKPOINT CLOSURE

**Date:** 2026-10-03
**Checkpoint:** M3-A2f (Regularized/Constrained Estimator Definition)
**Status:** CLOSED

---

## FINAL CLASSIFICATION

```
VALID ALGORITHMIC HEURISTIC / OPTIMIZATION ROBUSTNESS UNRESOLVED
```

---

## CHECKPOINT STATUS

### What Was Validated ✅

1. **Numerical implementation:** Correct
   - Likelihood agrees with scipy to machine precision
   - EM monotonicity confirmed
   - BIC convention correct

2. **Exact-family algorithmic behavior:** Promising
   - P(BIC selects K=2 | exact Product-vM K=1) = 0.13% (3/2250)
   - Stratified by q, κ, F
   - Rare K=2 selections at small F only

3. **K=2 power (q=2):** Excellent at clear separations
   - 100% power at separations ≥60°
   - Power stratified by separation/occupancy/κ/F
   - Close separations (10°) show expected κ/F dependence

4. **Mathematical foundation:** Corrected
   - Unbounded likelihood proven (valid diverging sequence)
   - 7 literature citations
   - Mixture model singularity documented

### What Remains Unresolved ⚠️

1. **Optimization convergence:** CRITICAL ISSUE
   - 97% non-convergence in wrapped-normal σ≈1.0 stress regime
   - Convergence rates across other conditions not audited
   - Most K=2 selections in stress regime are optimizer failures

2. **Shape-misspecification robustness:** Unclear
   - Raw 34% K=2 selection at σ≈1.0 (q=2, F=1000)
   - Cannot attribute to statistical sensitivity vs optimization failure
   - Stress regime performance unknown

3. **Initialization robustness:** Not tested
   - Current: split + 2 random initializations
   - Farthest-point, embedding k-means, perturbed K=1 not compared
   - Likelihood multimodality not explored

4. **Dimensional scaling:** Not validated
   - Only q=2 tested for K=2 power
   - q∈{1,5} required

5. **Dependent-toroidal robustness:** Not tested (MANDATORY for primary)
   - Bivariate sine von Mises not implemented
   - Correlated wrapped-normal not tested

6. **Regularization strategy:** Incomplete
   - Estimator U (numerical ceiling) tested
   - Estimator C (constrained) framework ready
   - Estimator P (penalized) literature source not implemented
   - Systematic U vs C vs P comparison not completed

### What Product-vM Is NOT ❌

- ❌ **NOT production-ready:** Convergence failures prevent deployment
- ❌ **NOT a validated screening detector:** Robustness unresolved
- ❌ **NOT a primary detector candidate:** Dependent-toroidal blocking

---

## SCIENTIFIC RECORD CORRECTIONS (Final)

### 1. Unbounded Likelihood Proof ✅ CORRECTED

**Original error:** π₁ → 0 and π₁√κ₁ → c > 0 (finite), then claimed divergence

**Corrected proof:**
- Keep π₁ = ε > 0 fixed
- Let κ₁ → ∞
- Then ε · √(κ₁/(2π)) → ∞

Alternative: If π₁ → 0, require π₁√κ₁ → ∞ (not finite constant)

### 2. Wrapped-Normal Unimodality ✅ CORRECTED

**Original error:** Claimed σ≈1 "theoretically in multimodality transition region" based on invalid discrete mode checker

**Corrected interpretation:**
- Independent wrapped-normal with σ<π is unimodal (analytical)
- σ≈0.7-1.2 tested range: all unimodal nulls
- Numerical mode checker had bug (equal neighboring values at symmetric peak)
- K=2 selections represent shape-misspecification sensitivity, BUT 97% non-convergence prevents clear attribution

### 3. Non-Converged Selections ✅ CORRECTED

**Original error:** Called raw K=2 selections "FPR"

**Corrected language:**
- "Raw K=2 selection frequency from local-EM heuristic"
- "Including non-converged fits"
- NOT "calibrated false-positive rate"
- Non-converged fits are optimizer failures, not valid detector decisions

### 4. Classification Language ✅ TIGHTENED

Explicitly states what Product-vM is NOT:
- NOT production-ready
- NOT validated screening detector
- NOT primary detector candidate

### 5. Phase Naming ✅ CORRECTED

Changed "M3-B (barrier analysis)" to "future M3 implementation/calibration phase"

Barrier/free-energy issues are separate audit blockers, not M3 scope.

---

## NEXT SCIENTIFIC BLOCKER

**Primary blocker:** Optimizer convergence and initialization robustness

**Before any deployment or calibration:**
1. Convergence audit across all validation conditions
2. Initialization robustness testing (multiple strategies)
3. Convergence diagnostics and failure modes
4. Stress-regime performance characterization

**Before primary candidate status:**
1. All above convergence items
2. Dependent-toroidal robustness test (MANDATORY)
3. Extended K=2 power (q>2)
4. U vs C vs P comparison with justified estimator selection

**Before production:**
- All above items
- Plus temporal correlation robustness
- Plus bootstrap calibration
- Plus independent external validation

---

## RECOMMENDATION FOR NEXT CHECKPOINT

**Two scientifically valid paths:**

### Path A: Complete M3-A2f Extended Validation

**Focus:** Resolve convergence and complete robustness testing

**Tasks:**
1. Convergence audit
2. Initialization robustness
3. Dependent toroidal (bivariate sine vM or correlated wrapped-normal)
4. Extended K=2 power (q∈{1,5})
5. U vs C vs P comparison
6. Literature-grounded penalized estimator

**Outcome:** Potential upgrade to "acceptable screening" or "primary candidate"

**Risk:** Convergence issue may not be easily resolved

### Path B: Proceed to M3-A3 (Alternative Detectors) ✅ PREFERRED

**Focus:** Comparative evaluation of alternative approaches

**Rationale:**
- Product-vM validated as functional heuristic
- Convergence issue documented
- Alternative detectors may address optimization failures
- M3-A3 not blocked by incomplete M3-A2f items
- Parallel evaluation is efficient

**Candidate alternatives:**
- Dirichlet Process mixture (nonparametric, avoids K selection)
- Hidden Markov Model (handles temporal correlation)
- Kernel density methods
- Other directional clustering approaches

**Outcome:** Comparative context for Product-vM; may inform refinement

---

## FILES DELIVERED

### Primary Deliverable
- **[M3_A2f_CHECKPOINT_REPORT.md](/tmp/M3_A2f_CHECKPOINT_REPORT.md)** - Complete scientific record

### Supporting Documentation
- [M3_A2f_LIKELIHOOD_DEGENERACY.md](/tmp/M3_A2f_LIKELIHOOD_DEGENERACY.md) - Mathematical foundation
- [M3_A2f_STATUS.md](/tmp/M3_A2f_STATUS.md) - Implementation framework
- [M3_A2f_SESSION_SUMMARY.md](/tmp/M3_A2f_SESSION_SUMMARY.md) - Session history
- [M3_A2f_CLOSURE.md](/tmp/M3_A2f_CLOSURE.md) - This closure document

### Implementation Artifacts
- [product_vm_three_estimators.py](/tmp/product_vm_three_estimators.py) - U, C, P framework
- [m3_a2f_analysis_framework.py](/tmp/m3_a2f_analysis_framework.py) - Analysis tools
- [investigate_wrapped_normal_anomaly.py](/tmp/investigate_wrapped_normal_anomaly.py) - Anomaly investigation
- [product_vm_fast_accurate.py](/tmp/product_vm_fast_accurate.py) - Estimator U (M3-A2e)

### Data Products
- [m3_a2e_checkpoint6_exact_family.csv](/tmp/m3_a2e_checkpoint6_exact_family.csv) - 3/2250 K=2
- [m3_a2e_checkpoint9_k2_recovery.csv](/tmp/m3_a2e_checkpoint9_k2_recovery.csv) - Power results
- [m3_a2e_checkpoint10_wrapped_normal.csv](/tmp/m3_a2e_checkpoint10_wrapped_normal.csv) - Shape sensitivity

---

## CLOSURE CHECKLIST

- ✅ Unbounded likelihood proof corrected
- ✅ Wrapped-normal unimodality interpretation corrected
- ✅ Non-converged selection language corrected
- ✅ Classification status tightened (explicit NOT statements)
- ✅ Phase naming corrected (not "barrier analysis")
- ✅ Final checkpoint report complete
- ✅ Scientific record accurate
- ✅ Convergence issue prominently documented
- ✅ Next blocker clearly stated
- ✅ Recommendation provided

---

## FINAL STATE

```
M3-A2f CHECKPOINT: CLOSED

Classification:
  VALID ALGORITHMIC HEURISTIC /
  OPTIMIZATION ROBUSTNESS UNRESOLVED

Next Scientific Blocker:
  Optimizer convergence and initialization robustness

Approved for Production:
  NO

Validated Screening Detector:
  NO

Primary Detector Candidate:
  NO

Recommended Next Checkpoint:
  M3-A3 (Alternative Detector Comparison)
```

---

**NO PRODUCTION CODE MODIFICATIONS**
**NO COMMITS**
**NO PUSH**
**NO M3-A3 STARTED**

**M3-A2f CHECKPOINT CLOSED**

---

**END OF M3-A2f**
