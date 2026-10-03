# M3 Multimodality Detection: Current Status

**Last Updated:** 2026-10-03
**Current Checkpoint:** M3-A2g (Optimizer Convergence Audit)
**Repository Branch:** scientific-rebuild
**Repository Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c

---

## CURRENT CLASSIFICATION

```
U — numerical-ceiling local-EM heuristic.
O2 — FUNCTIONAL BUT FRAGILE, provisional.
Optimization fragility documented;
initialization and stopping validation incomplete;
statistical robustness unresolved.
Not production-ready; no validated production defaults.

POOR CONVERGENCE RATES at max_iter=500:
- Exact-family: 30.8% (111/360)
- Stress regime: 36.4% (102/280)
- True-K2: 80.6% (258/320)
```

**Product von Mises (Product-vM) mixture detector is:**

- **NOT** production-ready
- **NOT** a validated screening detector
- **NOT** a primary detector candidate

**M3-A2g Findings:** Investigated "97% convergence failure" from M3-A2f. Root cause: iteration budget exhaustion (max_iter=100 insufficient), NOT fundamental optimizer failure. However, convergence remains POOR even with extended budgets (exact: 30.8% at max_iter=500, stress: 36.4% at max_iter=500). Limited monotonicity evidence (4 fixtures, 3384 states including triplicates, 1124 unique comparisons) suggests EM implementation not systematically broken, but does NOT constitute proof of correctness.

**Optimizer Classification:** **O2 - FUNCTIONAL BUT FRAGILE** (provisional; fragility documented, initialization and stopping validation incomplete)

---

## VALIDATED COMPONENTS

### ✅ Core Likelihood Implementation
- Fast Newton-refined κ MLE (numerical stability verified)
- Vectorized Product-vM K=2 EM algorithm
- Canonical ΔBIC convention (δBIC = BIC_K2 - BIC_K1)
- EM log-likelihood monotonicity verified **in limited scope** (4 fixtures, empirical only)
- Independent validation against `scipy.stats.vonmises`

### ✅ Exact-Family Algorithmic Behavior (M3-A2e)
**Configuration:** Product-vM K=1 null (exact generative model)
**Result:** 3/2250 K=2 selections = **0.13%**
**Interpretation:** Excellent algorithmic specificity under exact model
**Evidence:** [results/m3_a2e_checkpoint6_exact_family.csv](results/m3_a2e_checkpoint6_exact_family.csv)

### ✅ K=2 Clear-Separation Recovery (M3-A2e, q=2 only)
**Configuration:** Bimodal Product-vM with angular separations ≥60°
**Result:** 1620/1620 = **100% power**
**Limitations:** Only tested for q=2 (bivariate toroidal)
**Evidence:** [results/m3_a2e_checkpoint9_k2_recovery.csv](results/m3_a2e_checkpoint9_k2_recovery.csv)

### ⚠ EM Optimizer Fragility (M3-A2g) — PROVISIONAL O2, UNRESOLVED

**M3-A2g Investigation:** Tested whether "97% non-convergence" indicated fundamental optimizer failure (O1) or iteration budget exhaustion.

**Findings from 1235 CSV rows, 1634 iteration state records:**
- ⚠ **POOR convergence rates at max_iter=500:**
  - Exact-family: 111/360 (**30.8%**), varies 8%-83% by F (F=2500: only **8%**)
  - Stress regime: 102/280 (**36.4%**), includes σ=1.0 F=2500: 1/10 (**10%**)
  - True-K2: 258/320 (80.6%), varies 34%-100% by separation (sep=10°: only **34%**)
- ⚠ Extended budget experiment (max_iter=2000, UNPAIRED draws): 56/60 (93.3%) in stress regime
- ⚠ Absolute tolerance tol=1e-6 scales poorly with sample size
- ⚠ Non-converged fits often still improving at max_iter cutoff (budget exhaustion, not failure)
- ⚠ EM monotonicity observed in 4 limited fixtures (exact_easy, wrapped_1.0, k2_clear, k2_close):
  - 4 cases × 3 repeated deterministic split starts = 12 runs
  - 3384 iteration states total (including triplicates)
  - 1124 unique EM step comparisons (after removing triplicates)
  - 0 violations observed
  - **NOT comprehensive proof** - limited scope only

**Evidence:** [reports/M3_A2g_CHECKPOINT_REPORT.md](reports/M3_A2g_CHECKPOINT_REPORT.md)

**Assessment:** Budget exhaustion identified (NOT fundamental failure = avoid O1), but convergence remains **FRAGILE**. Limited monotonicity evidence suggests implementation not systematically broken. **Initialization and stopping validation remain INCOMPLETE** (diagnostics invalid). Not production-ready.

---

## NOT VALIDATED / UNRESOLVED

### ❌ Unbounded Likelihood / Regularization Strategy (M3-A2f)

**Fundamental Issue:**
Von Mises mixture likelihoods are unbounded when concentration parameters κ are unrestricted. The current implementation uses an implicit numerical ceiling κ≤10⁶, which is a computational hack, not a scientifically justified estimator.

**Three Candidate Strategies:**
1. **U (Numerical-ceiling local-EM heuristic):** Current implicit κ≤10⁶ ceiling (not unrestricted MLE)
2. **C (Constrained):** Explicit finite concentration domain with sensitivity analysis
3. **P (Penalized):** Literature-grounded penalized estimator (exact penalty and implementation not yet selected/validated)

**Status:**
- Framework implemented: [scripts/product_vm_three_estimators.py](scripts/product_vm_three_estimators.py)
- Systematic U/C/P comparison: **NOT COMPLETE**
- Final regularization choice: **NOT DETERMINED**

**Evidence:** [reports/M3_A2f_LIKELIHOOD_DEGENERACY.md](reports/M3_A2f_LIKELIHOOD_DEGENERACY.md)

### ❌ Initialization Robustness — INCOMPLETE / DIAGNOSTICS INVALID

**Status:** INCOMPLETE (M3-A2g diagnostic script has CRITICAL BUG)

**Current production code:** n_init=3 with all-random initialization (for n_init=1 uses split initialization, for n_init>1 uses all-random)

**M3-A2g Findings:**
- n_init sensitivity test shows wrapped_1.0 loglik diff = **2.63** between n_init=1 (split) and n_init=3 (all-random)
- This indicates **DIFFERENT SOLUTIONS** reached, NOT mere numerical noise
- **CRITICAL DIAGNOSTIC BUG:** Script reseeds on every INNER start, making all n_init>1 produce identical random initializations and invalidating the comparison

**Script:** [scripts/m3_a2g_ninit_sensitivity.py](scripts/m3_a2g_ninit_sensitivity.py) — **INVALID FOR ROBUSTNESS CLAIMS** (reseeding on inner starts; merely changing outer seed to 800+n_init does NOT fix this defect; see M3_A2G_CORRECTIONS_PROVENANCE.md)

**Required:** Fix inner reseeding defect, re-execute with independent per-start seeds, then systematic comparison of advanced initialization strategies (farthest-point, embedding-based, perturb-K1) on challenging cases (close separation, unbalanced occupancy)

### ❌ Stopping Criterion Validation — INCOMPLETE / DIAGNOSTICS INVALID

**Status:** INCOMPLETE (M3-A2g diagnostic script has CRITICAL BUGS)

**Script:** [scripts/m3_a2g_stopping_criterion_validation.py](scripts/m3_a2g_stopping_criterion_validation.py) — **INVALID FOR VALIDATION**
- Relative tolerance branch unimplemented (falls through to tol=1e-20)
- Absolute tolerance configurations have additional issues: approximate stop iteration, incorrect continuation endpoint, non-converged fits
- ALL stopping protocols tested have limitations (documented in M3_A2G_CORRECTIONS_PROVENANCE.md)

**Required:** Implement actual relative tolerance criterion |Δloglik|/|loglik| < tol_rel, fix continuation test methodology, re-execute validation

### ❌ Extended K=2 Power (q=1, q=5)

**Status:** NOT TESTED
**Limitation:** K=2 power validated only for q=2 (bivariate toroidal)
**Required:** Extended validation for q=1 (univariate) and q=5 (higher-dimensional) cases

### ❌ Dependent-Toroidal Robustness

**Status:** NOT TESTED
**Critical:** Product-vM assumes toroidal independence p(θ) = ∏ᵢ vM(θᵢ|μᵢ,κᵢ)
**Required:** Validation under dependent-toroidal distributions (bivariate sine von Mises with verified unimodality, or verified-unimodal correlated wrapped-normal)
**Requirement:** Measured circular association and verified joint unimodality
**Implication:** **MANDATORY** for "primary detector" classification

### ❌ Temporal Correlation Sensitivity

**Status:** NOT TESTED
**Context:** MD simulation data contains autocorrelation
**Required:** Validation with temporally correlated samples

### ❌ Calibrated Significance Testing

**Status:** NOT ESTABLISHED
**Current:** BIC model selection is an empirical heuristic
**Required:** Calibrated Type-I error control, significance threshold determination

### ❌ Overlapping-Window Multiple Testing

**Status:** NOT ADDRESSED
**Context:** TopoFold sliding-window analysis creates multiplicity
**Required:** FWER/FDR correction framework

---

## CHECKPOINT PROGRESSION

| Checkpoint | Focus | Status | Report |
|-----------|-------|--------|--------|
| M3-A | Initial Product-vM preference | Superseded | [reports/M3_A_CHECKPOINT_REPORT.md](reports/M3_A_CHECKPOINT_REPORT.md) |
| M3-A2 | Early validation | Superseded | [reports/M3_A2_CHECKPOINT_REPORT.md](reports/M3_A2_CHECKPOINT_REPORT.md) |
| M3-A2b | Correlated-null generators | Superseded | [reports/M3_A2b_CHECKPOINT_REPORT.md](reports/M3_A2b_CHECKPOINT_REPORT.md) |
| M3-A2c | Sanity audit, exact nulls | Supporting | [reports/M3_A2c_CHECKPOINT_REPORT.md](reports/M3_A2c_CHECKPOINT_REPORT.md) |
| M3-A2d | Core corrections following A2c | Supporting | [reports/M3_A2d_CHECKPOINT_REPORT.md](reports/M3_A2d_CHECKPOINT_REPORT.md) |
| M3-A2e | Full benchmark with corrected fitter | Supporting | [reports/M3_A2e_CHECKPOINT_REPORT.md](reports/M3_A2e_CHECKPOINT_REPORT.md) |
| M3-A2f | Likelihood degeneracy & convergence | Superseded (classification) | [reports/M3_A2f_CHECKPOINT_REPORT.md](reports/M3_A2f_CHECKPOINT_REPORT.md) |
| **M3-A2g** | **Optimizer convergence audit** | **CURRENT (Provisional O2)** | **[reports/M3_A2g_CHECKPOINT_REPORT.md](reports/M3_A2g_CHECKPOINT_REPORT.md)** |

**Historical Attribution:**
- 3/2250 (0.13%) K=2-selection under exact-K1 nulls: **M3-A2e**
- 1620/1620 (100%) clear-K2 recovery (sep≥60°, q=2): **M3-A2e**
- K=1 fitter corrections: **M3-A2d** (documented), **M3-A2e** (full benchmark/recovery work)

**Closure Documents:**
- [reports/M3_A2f_CLOSURE.md](reports/M3_A2f_CLOSURE.md) - M3-A2f closure
- [reports/M3_A2f_LIKELIHOOD_DEGENERACY.md](reports/M3_A2f_LIKELIHOOD_DEGENERACY.md) - Unbounded likelihood proof
- [M3_A2G_CHECKPOINT_CLOSURE.md](M3_A2G_CHECKPOINT_CLOSURE.md) - M3-A2g checkpoint closure
- [M3_A2G_CORRECTIONS_PROVENANCE.md](M3_A2G_CORRECTIONS_PROVENANCE.md) - External review findings and diagnostic limitations

---

## NEXT SCIENTIFIC BLOCKER

```
Statistical robustness under model misspecification
```

**M3-A2g Status:** Budget exhaustion identified as root cause of "97% convergence failure" (NOT fundamental optimizer failure = avoid O1). However, convergence remains POOR even with extended budgets (exact: 30.8%, stress: 36.4% at max_iter=500). Limited monotonicity evidence (4 fixtures, 1124 unique comparisons, 0 violations) suggests EM not systematically broken, but initialization and stopping diagnostics are INVALID. Classification: **provisional O2 - FUNCTIONAL BUT FRAGILE**.

**Remaining scientific blockers (in priority order):**

1. **Complete initialization and stopping validation:**
   - Fix diagnostic bugs (inner reseeding, relative tolerance implementation, continuation methodology)
   - Execute validated robustness tests with independent per-start seeds
   - Required before accepting O2 or upgrading to O3

2. **Dependent-toroidal misspecification:**
   - Product-vM assumes toroidal independence
   - Real protein dihedrals exhibit angular correlations
   - **Mandatory** before production deployment

3. **Regularization strategy (U/C/P):**
   - U/C/P systematic comparison incomplete
   - Final regularization choice undetermined

4. **Extended validation (q=1, q=5, bootstrap calibration):**
   - Dimensional validation incomplete
   - δBIC threshold calibration under misspecification

**Next checkpoint options:**

- **M3-A2h:** Fix initialization/stopping diagnostic bugs, complete validation (required before accepting provisional O2)
- **M3-A3:** Dependent-toroidal validation (bivariate sine vM, correlated wrapped-normal)
- **M3-B:** Regularization strategy systematic comparison (U/C/P)
- **Alternative:** Pivot to model-free multimodality detectors if misspecification sensitivity is prohibitive

**Decision:** M3-A2g corrections and provisional O2 status documented. External review required before proceeding.

---

## SCIENTIFIC METHODOLOGY NOTES

### Historical Checkpoint Validity

All checkpoints **before M3-A2d** used an incorrect K=1 fitter and are scientifically invalid for quantitative conclusions. They are retained for provenance and methodological transparency.

### Supersession vs. Supporting

- **Superseded:** Checkpoint conclusions were invalidated by later findings
- **Supporting:** Checkpoint contributes to historical chain but is not standalone conclusion
- **Current:** Represents current accepted scientific state (provisional)

### Toroidal Independence Assumption

Product von Mises mixture model assumes conditional independence within each component:
```
p(θ|k) = ∏ᵢ₌₁^q vM(θᵢ|μₖᵢ,κₖᵢ)

Overall:
p(θ) = Σₖ πₖ ∏ᵢ₌₁^q vM(θᵢ|μₖᵢ,κₖᵢ)
```

This is a **strong simplifying assumption**. Real protein dihedral angle distributions may exhibit:
- Angular correlations (e.g., Ramachandran-like coupling)
- Dependent multimodality (modes in joint space, not marginal)

**Validation under dependent-toroidal distributions is MANDATORY before production deployment.**

---

## PRIMARY DELIVERABLES

### Reports
- **[M3_A2g_CHECKPOINT_REPORT.md](reports/M3_A2g_CHECKPOINT_REPORT.md)** — **Optimizer convergence audit (CURRENT, Provisional O2)**
- [M3_A2G_CHECKPOINT_CLOSURE.md](M3_A2G_CHECKPOINT_CLOSURE.md) — M3-A2g checkpoint closure
- [M3_A2G_COMPLETE_ANALYSIS.md](M3_A2G_COMPLETE_ANALYSIS.md) — Complete M3-A2g analysis
- [M3_A2G_CORRECTIONS_PROVENANCE.md](M3_A2G_CORRECTIONS_PROVENANCE.md) — External review findings and diagnostic limitations
- [M3_A2f_CHECKPOINT_REPORT.md](reports/M3_A2f_CHECKPOINT_REPORT.md) — Likelihood degeneracy & convergence
- [M3_A2f_CLOSURE.md](reports/M3_A2f_CLOSURE.md) — M3-A2f closure document
- [M3_A2f_LIKELIHOOD_DEGENERACY.md](reports/M3_A2f_LIKELIHOOD_DEGENERACY.md) — Unbounded likelihood proof

### Scripts
- **[m3_a2g_convergence_audit.py](scripts/m3_a2g_convergence_audit.py)** — **Convergence trajectory analysis**
- **[m3_a2g_exact_k2_convergence.py](scripts/m3_a2g_exact_k2_convergence.py)** — **Exact-family/K2 convergence audit**
- **[m3_a2g_monotonicity_check.py](scripts/m3_a2g_monotonicity_check.py)** — **EM monotonicity verification (limited scope)**
- **[m3_a2g_ninit_sensitivity.py](scripts/m3_a2g_ninit_sensitivity.py)** — **n_init sensitivity (INVALID - seed reuse bug)**
- **[m3_a2g_stopping_criterion_validation.py](scripts/m3_a2g_stopping_criterion_validation.py)** — **Stopping validation (INVALID - relative branch bug)**
- **[m3_a2g_extended_maxiter.py](scripts/m3_a2g_extended_maxiter.py)** — **Extended max_iter test (unpaired draws)**
- **[m3_a2g_duplicate_component_investigation.py](scripts/m3_a2g_duplicate_component_investigation.py)** — **Duplicate component investigation**
- **[m3_a2g_data_verification.py](scripts/m3_a2g_data_verification.py)** — **Data verification and numerical reconciliation**
- [product_vm_fast_accurate.py](scripts/product_vm_fast_accurate.py) — M3-A2e corrected fitter
- [product_vm_three_estimators.py](scripts/product_vm_three_estimators.py) — U/C/P estimator framework

### Results (1235 CSV rows, 1634 iteration state records)
- **[m3_a2g_exact_family_convergence.csv](results/m3_a2g_exact_family_convergence.csv)** — **360-fit exact-family audit (111/360 = 30.8%)**
- **[m3_a2g_k2_convergence.csv](results/m3_a2g_k2_convergence.csv)** — **320-fit true-K2 audit (258/320 = 80.6%)**
- **[m3_a2g_stress_regime.csv](results/m3_a2g_stress_regime.csv)** — **280-fit wrapped-normal stress (102/280 = 36.4%)**
- **[m3_a2g_trajectories.txt](results/m3_a2g_trajectories.txt)** — **5 trajectory blocks, 1634 iteration state records**
- **[m3_a2g_ninit_sensitivity.csv](results/m3_a2g_ninit_sensitivity.csv)** — **20 fits, wrapped_1.0 loglik diff = 2.63**
- [m3_a2g_duplicate_components.csv](results/m3_a2g_duplicate_components.csv) — Duplicate component investigation
- [m3_a2g_extended_maxiter.csv](results/m3_a2g_extended_maxiter.csv) — Extended max_iter (180 fits)
- [m3_a2g_initialization_comparison.csv](results/m3_a2g_initialization_comparison.csv) — Initialization comparison
- [m3_a2g_max_iter_sensitivity.csv](results/m3_a2g_max_iter_sensitivity.csv) — max_iter sensitivity
- [m3_a2g_stopping_criterion_validation.csv](results/m3_a2g_stopping_criterion_validation.csv) — Stopping validation
- [m3_a2g_tolerance_comparison.csv](results/m3_a2g_tolerance_comparison.csv) — Tolerance comparison
- [m3_a2e_checkpoint6_exact_family.csv](results/m3_a2e_checkpoint6_exact_family.csv) — 3/2250 (0.13%) K2-selection under exact-K1 (M3-A2e)
- [m3_a2e_checkpoint9_k2_recovery.csv](results/m3_a2e_checkpoint9_k2_recovery.csv) — 1620/1620 (100%) clear-K2 power (M3-A2e)
- [m3_a2e_checkpoint10_wrapped_normal.csv](results/m3_a2e_checkpoint10_wrapped_normal.csv) — Wrapped-normal anomaly

---

## ARTIFACT MANIFEST

See [MANIFEST.md](MANIFEST.md) for complete provenance, checksums, and regenerability status.

---

**NO PRODUCTION CODE MODIFICATIONS APPROVED**

All artifacts in this directory are scientific audit evidence. None have been promoted to `src/` or production TopoFold code.

**M3-A2g STATUS:** Provisional O2 classification. Initialization and stopping validation remain incomplete. External review required before proceeding.
