# M3-A2f CHECKPOINT REPORT: Product von Mises Regularization & Validation

**Date:** 2026-10-03
**Phase:** M3-A (Multimodality Detection Validation)
**Checkpoint:** M3-A2f (Regularized/Constrained Estimator Definition)
**Status:** COMPLETE

---

## EXECUTIVE SUMMARY

M3-A2f addresses the fundamental statistical degeneracy in von Mises mixture models and completes validation of the Product von Mises detector.

**Key Findings:**
1. ✅ Von Mises mixture likelihoods are mathematically unbounded (7 literature citations)
2. ✅ M3-A2e algorithmic benchmark shows 0.13% K=2 selection on exact Product-vM K=1 nulls
3. ⚠️ **Wrapped-normal shape-misspecification sensitivity:** Finite-sample transition at σ≈0.9-1.0 with 97% non-convergence
4. ✅ K=2 power: 100% at separations ≥60° (q=2 only, validated)
5. ⏸ Extended validation (q>2, dependent-toroidal, U vs C vs P) deferred due to convergence issue priority

**Classification:** VALID ALGORITHMIC HEURISTIC / OPTIMIZATION ROBUSTNESS UNRESOLVED

**Status:**
- ✅ Numerical implementation validated
- ✅ Exact-family algorithmic behavior promising (0.13% K=2 selection)
- ⚠️ Statistical robustness under model misspecification UNRESOLVED
- ⚠️ Optimization failures prevent reliable inference in stress regimes
- ❌ Product-vM is NOT production-ready
- ❌ Product-vM is NOT yet a validated screening detector
- ❌ Product-vM is NOT a primary detector candidate

**Recommendation:** Resolve optimizer convergence and initialization robustness before any deployment or further statistical calibration

---

## 1. CORRECTED MATHEMATICAL FOUNDATION

### 1.1 Likelihood Unboundedness

**Theorem:** The likelihood function of a finite von Mises mixture is unbounded when concentration parameters are unrestricted.

**Proof:** Consider n observations {θ₁, ..., θₙ} and K≥2 component mixture. Construct diverging sequence:

1. Fix component 1 center at observation θ₁: μ₁ = θ₁
2. Keep component 1 weight fixed: π₁ = ε > 0 (constant)
3. Let κ₁ → ∞
4. Keep remaining components with positive density on all observations

The likelihood contribution from θ₁:

```
L₁ = Σₖ πₖ f(θ₁ | μₖ, κₖ)
   ≥ π₁ f(θ₁ | μ₁, κ₁)
   ~ ε · √(κ₁/(2π))    [using I₀(κ) ~ exp(κ)/√(2πκ)]
   → ∞   as κ₁ → ∞
```

Therefore log L(θ₁) → ∞, proving unboundedness. ∎

### 1.2 Literature Support

**Primary sources:**
- Redner & Walker (1984) - Mixture degeneracy recognition
- Banerjee et al. (2005) - vMF mixtures require κ ≤ κ_max constraint
- Hornik & Grün (2014) - movMF package constrains κ to finite interval
- Mardia & Jupp (2000) - Directional Statistics (identifiability issues)

**Implication:** Unrestricted MLE does not exist. Any valid estimator must impose regularization.

### 1.3 Mixture Model Singularity

K=1 vs K=2 comparison is non-regular/singular:
- Component labels non-identifiable
- K=1 is on boundary of parameter space
- Standard BIC asymptotics do not apply

**Model-selection interpretation:** BIC is an empirical heuristic, NOT a calibrated significance test. Formal calibration requires bootstrap validation (future M3-B phase).

---

## 2. COMPLETED M3-A2E ALGORITHMIC BENCHMARK

M3-A2e tested estimator U (numerical-ceiling EM heuristic with κ ≤ 10⁶).

### 2.1 Exact Product-vM K=1 Nulls (Checkpoint 6)

**Configuration:** q ∈ {1,2,5}, κ ∈ {1,5,20}, F ∈ {100,250,500,1000,2500}, 50 reps/config

**Results:**
- **Total:** 3/2250 K=2 selections = 0.13%
- **Breakdown:**
  - q=1, κ=1, F=100: 1/50
  - q=1, κ=5, F=100: 1/50
  - q=1, κ=20, F=250: 1/50
  - All other 42 configs: 0/50

**Interpretation:** P(BIC selects K=2 | exact Product-vM K=1) = 0.13%

This is NOT "Type-I error at α=0.05" - it is empirical selection frequency under this algorithmic heuristic.

**Stratification:**
- By dimensionality: q=1 (3/750=0.4%), q=2 (0/750=0%), q=5 (0/750=0%)
- By concentration: κ=1 (1/750=0.13%), κ=5 (1/750=0.13%), κ=20 (1/750=0.13%)
- By sample size: F=100 (2/450=0.44%), F=250 (1/450=0.22%), F≥500 (0/1350=0%)

### 2.2 K=2 Recovery Power (Checkpoint 9, q=2 only)

**Configuration:** Separations {10°,20°,30°,60°,90°,180°}, occupancies {50/50, 80/20}, κ ∈ {5,20,100}, F ∈ {100,500,1000}, 30 reps

**Results by separation:**
- 10°: 0-70% (depends heavily on κ, F)
- 20°: 10-100%
- 30°: 60-100%
- 60°: 100% (all configs)
- 90°: 100% (all configs)
- 180°: 100% (all configs)

**Stratification critical:** "100% power at ≥60°" masks variation at close separations.

**Example detail (separation=10°):**
- κ=5, F=100: 0% power
- κ=20, F=1000: 13% power
- κ=100, F=500: 100% power

**Dimensionality:** ONLY q=2 tested. Extension to q∈{1,5} required before final classification.

### 2.3 Wrapped-Normal Shape Misspecification (Checkpoint 10)

**Original finding:** q=2, σ=1.0, F=1000: 15/30 = 50%

**Follow-up investigation (M3-A2f):**

**Reproducibility test (n=100 reps):**
- Raw K=2 selection frequency: 33% (including non-converged)
- **CRITICAL: Convergence rate: 3%** (97% non-convergence!)

**σ transition surface (q=2, F=1000, n=50):**

| σ | K=2 Selection Rate |
|---|-------------------|
| 0.4-0.6 | 0% |
| 0.7 | 6% |
| 0.8 | 12% |
| **0.9** | **36%** |
| **1.0** | **34%** |
| 1.1 | 18% |
| 1.2 | 6% |
| ≥1.3 | 0% |

Note: These are raw K=2 selection frequencies from the local-EM heuristic, including non-converged fits. With 97% non-convergence at peak, these do NOT represent calibrated false-positive rates.

**Scientific interpretation:**

Independent wrapped-normal with σ≈0.7-1.2 is unimodal (verified analytically for σ<π). The rise in K=2 selections represents shape-misspecification sensitivity, but:

1. **Shape misspecification:** σ≈0.9-1.0 wrapped-normal is diffuse enough that Product-vM attempts two-component approximation
2. **Sample-size dependent:** F=250→0%, F=1000→34% (larger sample enables finer approximation)
3. **Convergence failure dominates:** 97% non-convergence prevents attributing K=2 selections purely to statistical model misspecification
4. **Interpretation:** Most K=2 selections are optimizer failures, not valid detector decisions

**Wrapping origin invariance:** ✓ Verified (10/10 agreement across wrapping conventions)

**Recommendation:** Optimization convergence failures prevent reliable inference. Initialization robustness and convergence diagnostics required before any deployment.

---

## 3. ESTIMATOR DEFINITIONS

### 3.1 Estimator U: Numerical-Ceiling Algorithmic Heuristic

**Definition:**
```
Local EM optimization
Concentration bound: 0 ≤ κ ≤ 10⁶ (numerical ceiling)
3 random initializations, best local solution by likelihood
BIC model selection
```

**Status:** This is the M3-A2e benchmark estimator

**Label:** "3-start EM with numerical κ≤10⁶ ceiling, ordinary BIC"

**NOT claimed:** Unrestricted MLE (does not exist)

### 3.2 Estimator C: Constrained Likelihood Family

**Definition:**
```
Local EM optimization
Explicit finite concentration domain: 0 ≤ κ ≤ κ_max
κ_max subjected to sensitivity analysis
BIC model selection
```

**Status:** Framework implemented, sensitivity analysis incomplete due to convergence issue priority

**Sensitivity grid prepared:** κ_max ∈ {20, 50, 100, 200, 500, 1000}

**Note:** κ_max must be justified by interpretable angular-resolution scale, not chosen to produce desired classification

### 3.3 Estimator P: Literature-Grounded Penalized Likelihood

**Status:** Framework prepared but literature-supported penalty not yet implemented

**Original Gamma(2, 0.1) claim:** REMOVED - could not verify in Hornik & Grün (2014)

**Required:** Identify published penalty from:
- Chen, Li & Tan (2007)
- Ng (2023) - Penalized MLE for vMF mixtures
- Other primary sources

**Implementation requirements:**
- Exact penalty form from source
- Sample-size dependence
- Optimization objective
- Separate tracking: raw logL, penalty, penalized objective

### 3.4 Optional Estimator G: Exploratory Gamma-MAP

**If retained, requires mathematical correction:**

Penalized score:
```
score = n_eff * (R_bar - A(κ)) + (α-1)/κ - β
```

Information (CORRECTED sign):
```
information = n_eff * A'(κ) + (α-1)/κ²
```

NOT: `- (α-1)/κ²` (wrong sign in original)

**Status:** Deferred pending literature-grounded P implementation

---

## 4. CONVERGENCE ROBUSTNESS ISSUE

**Critical finding from wrapped-normal investigation:**

At shape-misspecification transition (q=2, σ=1.0, F=1000):
- Convergence rate: 3% (97% non-convergence!)
- FPR correlation with convergence unclear

**Implication:** The 0.13% exact-family K=2 rate may underestimate true sensitivity if:
1. Exact-family data has higher convergence rate
2. Non-convergence protects against false K=2 (selection defaults to K=1)

**Required before production deployment:**
1. Audit convergence rates across all validation conditions
2. Investigate initialization robustness
3. Compare multiple initialization strategies:
   - Random responsibilities
   - Toroidal farthest-point
   - Sin/cos embedding k-means (initialization only)
   - Split/perturbed K=1

**Status:** Convergence audit incomplete (deferred due to scope)

---

## 5. OUTSTANDING VALIDATION REQUIREMENTS

### 5.1 Extended K=2 Power (q>2)

**M3-A2e tested:** q=2 only

**Required:** q ∈ {1, 2, 5}

**Reduced grid:**
- Separations: {10°, 30°, 60°, 90°}
- κ: {5, 20}
- F: {250, 1000}
- Occupancy: {50/50, 80/20}

**Rationale:** Distinguish dimensional accumulation of evidence from genuine robustness

**Status:** Not completed (deferred due to convergence issue priority)

### 5.2 Dependent Toroidal Robustness (MANDATORY for "primary")

**Requirement:** At least one verified-unimodal dependent toroidal null

**Options:**
1. Corrected bivariate sine von Mises with λ² < κ₁·κ₂ (safe margin from boundary)
2. Verified-unimodal correlated wrapped-normal

**Report:**
- Exact generating model
- Theoretical/numerical unimodality verification
- Measured circular association
- K=2 selection behavior
- Convergence diagnostics

**Status:** Not completed (BLOCKS "primary candidate" classification)

### 5.3 Initialization Robustness

**M3-A2e used:** Split + 2 random initializations

**Required comparison:**
- Random responsibilities
- Toroidal farthest-point
- Sin/cos embedding clustering (initialization only)
- Split/perturbed K=1

**Report:** Likelihood multimodality, selected K stability, convergence

**Status:** Not completed

### 5.4 U vs C vs P Systematic Comparison

**Required:**
- Paired comparisons on same generated datasets
- Common random numbers
- All three valid estimators

**Grids:**
- Exact Product-vM K=1 (multiple q, κ, F)
- Product-vM K=2 (separations, occupancy, κ, F, q)
- Wrapped-normal (σ transition region)
- Dependent toroidal (when implemented)

**Report per estimator:**
- Selected K
- Raw likelihood
- Model criterion
- Convergence
- Parameters
- Boundary/regularization diagnostics
- Runtime

**Status:** Framework implemented, validation not completed

---

## 6. CLASSIFICATION

### 6.1 Classification Options

1. IMPLEMENTATION INVALID
2. VALID ALGORITHMIC HEURISTIC / STATISTICS UNRESOLVED
3. REJECTED AS PRIMARY DETECTOR
4. ACCEPTABLE SCREENING CANDIDATE
5. PRIMARY CANDIDATE FOR BOOTSTRAP CALIBRATION

### 6.2 Assessment

**Option 1: IMPLEMENTATION INVALID** ❌
- M3-A2e validated core numerical correctness
- Likelihood agrees with scipy
- EM monotonicity confirmed
- NOT applicable

**Option 2: VALID ALGORITHMIC HEURISTIC / STATISTICS UNRESOLVED** ✅ **SELECTED**

**Rationale:**
- ✅ Fitter is numerically correct
- ✅ Exact-family behavior excellent (0.13% K=2 rate)
- ✅ K=2 power excellent at clear separations (q=2)
- ⚠️ **Convergence robustness UNRESOLVED** (97% non-convergence at σ≈1.0 transition)
- ⏸ q>2 power not validated
- ⏸ Dependent-toroidal not tested
- ⏸ Initialization robustness not tested
- ⏸ U vs C vs P comparison incomplete

**Option 3: REJECTED AS PRIMARY** ❌
- No evidence of fundamental failure
- Convergence issue may be addressable
- NOT rejected, but not promoted yet

**Option 4: ACCEPTABLE SCREENING CANDIDATE** ⏸ **POSSIBLE after convergence resolution**
- If convergence failures are initialization-addressable
- If extended validation (q>2, dependent-toroidal) passes
- Requires completion of outstanding validation

**Option 5: PRIMARY CANDIDATE** ❌ **BLOCKED**
- Dependent-toroidal mandatory test not completed
- Convergence robustness not established
- Cannot proceed to M3-B without these

### 6.3 Final Classification

**Product von Mises Mixture Detector (Estimator U):**

**VALID ALGORITHMIC HEURISTIC / OPTIMIZATION ROBUSTNESS UNRESOLVED**

**Strengths:**
- Exact-family control: excellent (0.13%)
- K=2 power at clear separations: excellent (100% at ≥60°, q=2)
- Core implementation: validated

**Critical Unresolved Issues:**
1. **Convergence failures:** 97% non-convergence in shape-misspecification transition
2. **Initialization robustness:** Not tested
3. **Dimensional scaling:** q>2 not validated
4. **Dependent-toroidal robustness:** Not tested (MANDATORY)

**Not classified as:**
- "Production-ready" (convergence issue)
- "Primary candidate" (dependent-toroidal blocking)
- "Acceptable screening" (pending convergence resolution)

---

## 7. EXACT BLOCKERS

### 7.1 Blockers for "Acceptable Screening Candidate"

1. ⏸ **Convergence audit** across all validation conditions
2. ⏸ **Initialization robustness** demonstration
3. ⏸ **Extended K=2 power** (q>2) validation
4. ⏸ **Dependent-toroidal robustness** test

### 7.2 Blockers for "Primary Candidate for Bootstrap Calibration"

All of above, plus:

5. ⏸ **U vs C vs P comparison** showing U is competitive
6. ⏸ **Regularization justification** (if C or P selected)
7. ⏸ **Temporal correlation robustness** (future M3-C)

### 7.3 Blockers for M3-A3 (Alternative Detectors)

**M3-A3 is NOT blocked by M3-A2f incomplete items.**

Rationale:
- Product-vM algorithmic heuristic is validated as working detector
- Convergence issues do not prevent comparison with alternatives
- Alternative detectors may have better convergence
- Parallel evaluation is scientifically valid

**M3-A3 may proceed** while M3-A2f extended validation continues in parallel.

### 7.4 Blockers for Future M3 Implementation/Calibration Phase

**Future M3 phase requires PRIMARY detector candidate** - currently BLOCKED:

1. ⏸ Dependent-toroidal robustness (MANDATORY)
2. ⏸ Convergence robustness established
3. ⏸ Final estimator selection (U/C/P) with justification
4. ⏸ Classification upgrade to "primary candidate"

---

## 8. RECOMMENDATION FOR NEXT SCIENTIFIC CHECKPOINT

### 8.1 Option A: Complete M3-A2f Extended Validation

**Priority tasks:**
1. Convergence audit across all M3-A2e conditions
2. Initialization robustness testing
3. Dependent toroidal robustness (bivariate sine vM or correlated wrapped-normal)
4. Extended K=2 power (q∈{1,5})
5. U vs C vs P comparison

**Estimated time:** 4-8 hours computation + analysis

**Outcome:** Potential upgrade to "acceptable screening" or "primary candidate"

### 8.2 Option B: Proceed to M3-A3 (Alternative Detectors)

**Rationale:**
- Product-vM algorithmic heuristic validated as functional
- Alternative detectors may address convergence issues
- Parallel evaluation is efficient
- Does NOT require completing M3-A2f extended validation first

**Candidate alternatives:**
- Dirichlet Process mixture (nonparametric)
- Hidden Markov Model (temporal correlation)
- Kernel density ratio test
- Other directional clustering methods

**Estimated time:** Variable (depends on alternatives selected)

**Outcome:** Comparative evaluation may inform Product-vM refinement

### 8.3 Recommended Path

**Immediate:** Document convergence issue prominently

**Short-term:** Proceed to M3-A3 while leaving M3-A2f extended validation for parallel/future work

**Rationale:**
1. Convergence issue is documented
2. Product-vM is validated as functional heuristic
3. Alternative detector comparison is scientifically valuable
4. M3-A2f extended validation can continue in parallel
5. M3-B requires primary detector (blocked regardless)

**Do NOT:**
- Deploy Product-vM in production without convergence resolution
- Claim "primary candidate" status
- Proceed to M3-B barrier analysis yet

---

## 9. FILES CREATED (M3-A2f)

### Documentation
- [M3_A2f_LIKELIHOOD_DEGENERACY.md](/tmp/M3_A2f_LIKELIHOOD_DEGENERACY.md) - Corrected mathematics with 7 citations
- [M3_A2f_STATUS.md](/tmp/M3_A2f_STATUS.md) - Implementation framework status
- [M3_A2f_SESSION_SUMMARY.md](/tmp/M3_A2f_SESSION_SUMMARY.md) - Session accomplishments
- [M3_A2f_CHECKPOINT_REPORT.md](/tmp/M3_A2f_CHECKPOINT_REPORT.md) - This final report

### Implementation
- [product_vm_three_estimators.py](/tmp/product_vm_three_estimators.py) - U, C, P estimator framework (~800 lines)
- [m3_a2f_analysis_framework.py](/tmp/m3_a2f_analysis_framework.py) - Corrected analysis tools (~400 lines)
- [investigate_wrapped_normal_anomaly.py](/tmp/investigate_wrapped_normal_anomaly.py) - Anomaly investigation (~300 lines)

### From M3-A2e (Preserved)
- [product_vm_fast_accurate.py](/tmp/product_vm_fast_accurate.py) - Estimator U (516 lines)
- [m3_a2e_checkpoint6_exact_family.csv](/tmp/m3_a2e_checkpoint6_exact_family.csv) - Exact-family results
- [m3_a2e_checkpoint9_k2_recovery.csv](/tmp/m3_a2e_checkpoint9_k2_recovery.csv) - K=2 power results
- [m3_a2e_checkpoint10_wrapped_normal.csv](/tmp/m3_a2e_checkpoint10_wrapped_normal.csv) - Wrapped-normal results

---

## 10. CONCLUSION

M3-A2f establishes that:

1. ✅ Von Mises mixture degeneracy is mathematically proven and literature-supported
2. ✅ Product-vM estimator U (numerical-ceiling heuristic) is numerically correct
3. ✅ Exact-family control is excellent (0.13%)
4. ✅ K=2 power at clear separations is excellent (100% at ≥60°, q=2)
5. ⚠️ **Convergence robustness is UNRESOLVED** (97% failures in shape-misspecification transition)
6. ⏸ Extended validation incomplete (q>2, dependent-toroidal, initialization, U vs C vs P)

**Classification:** VALID ALGORITHMIC HEURISTIC / OPTIMIZATION ROBUSTNESS UNRESOLVED

**Next checkpoint options:**
- Complete M3-A2f extended validation (convergence, q>2, dependent-toroidal)
- OR proceed to M3-A3 alternative detector comparison

**Future M3 implementation/calibration:** BLOCKED (requires primary candidate with dependent-toroidal validation)

---

**NO PRODUCTION CODE MODIFICATIONS**
**NO COMMITS**
**NO PUSH**

**M3-A2f CHECKPOINT: COMPLETE**

---

**END OF M3-A2f CHECKPOINT REPORT**
