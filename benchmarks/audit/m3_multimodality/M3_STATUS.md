# M3 Multimodality Detection: Current Status

**Last Updated:** 2026-10-03
**Current Checkpoint:** M3-A2f CLOSED
**Repository Branch:** scientific-rebuild
**Repository HEAD:** 95244fd486be1a96e29cc0cb37bb5cb244fa9407

---

## CURRENT CLASSIFICATION

```
VALID ALGORITHMIC HEURISTIC /
OPTIMIZATION ROBUSTNESS UNRESOLVED
```

**Product von Mises (Product-vM) mixture detector is:**

- **NOT** production-ready
- **NOT** a validated screening detector
- **NOT** a primary detector candidate

---

## VALIDATED COMPONENTS

### ✅ Core Likelihood Implementation
- Fast Newton-refined κ MLE (numerical stability verified)
- Vectorized Product-vM K=2 EM algorithm
- Canonical ΔBIC convention (δBIC = BIC_K2 - BIC_K1)
- EM log-likelihood monotonicity verified
- Independent validation against `scipy.stats.vonmises`

### ✅ Exact-Family Algorithmic Behavior
**Configuration:** Product-vM K=1 null (exact generative model)
**Result:** 3/2250 K=2 selections = **0.13%**
**Interpretation:** Excellent algorithmic specificity under exact model
**Evidence:** [results/m3_a2e_checkpoint6_exact_family.csv](results/m3_a2e_checkpoint6_exact_family.csv)

### ✅ K=2 Clear-Separation Recovery (q=2 only)
**Configuration:** Bimodal Product-vM with angular separations ≥60°
**Result:** 1620/1620 = **100% power**
**Limitations:** Only tested for q=2 (bivariate toroidal)
**Evidence:** [results/m3_a2e_checkpoint9_k2_recovery.csv](results/m3_a2e_checkpoint9_k2_recovery.csv)

---

## NOT VALIDATED / UNRESOLVED

### ⚠️ PRIMARY BLOCKER: Optimizer Convergence Robustness

**Critical Finding (M3-A2f):**
- Wrapped-normal σ≈1.0 transition: **97% EM convergence failure**
- Raw K=2 selection frequency: 33-36% (including non-converged fits)
- **Cannot distinguish:** statistical sensitivity vs. optimizer failure
- Most K=2 selections are non-converged local optima, not valid detector decisions

**Implication:**
Optimization robustness is the primary scientific blocker. Product-vM cannot be production-ready until convergence rates are characterized and initialization strategies are validated.

**Evidence:** [scripts/investigate_wrapped_normal_anomaly.py](scripts/investigate_wrapped_normal_anomaly.py)

### ❌ Unbounded Likelihood / Regularization Strategy

**Fundamental Issue (M3-A2f):**
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

### ❌ Initialization Robustness

**Status:** NOT TESTED
**Current:** Split initialization + 2 random initializations per EM fit
**Required:** Systematic evaluation of initialization strategies, convergence rate vs. n_init

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
| M3-A2d | Corrected fitter | Supporting | [reports/M3_A2d_CHECKPOINT_REPORT.md](reports/M3_A2d_CHECKPOINT_REPORT.md) |
| M3-A2e | Full benchmark with corrected fitter | Superseded | [reports/M3_A2e_CHECKPOINT_REPORT.md](reports/M3_A2e_CHECKPOINT_REPORT.md) |
| **M3-A2f** | **Likelihood degeneracy & convergence** | **CURRENT** | **[reports/M3_A2f_CHECKPOINT_REPORT.md](reports/M3_A2f_CHECKPOINT_REPORT.md)** |

**Closure Document:** [reports/M3_A2f_CLOSURE.md](reports/M3_A2f_CLOSURE.md)

---

## NEXT SCIENTIFIC BLOCKER

```
Optimizer convergence and initialization robustness
```

**Two possible next steps:**

1. **M3-A2g-style convergence work:**
   - Characterize EM convergence rates across parameter space
   - Evaluate initialization strategies (k-means, grid, mixture sampling)
   - Determine convergence diagnostics and retry policies
   - Complete U/C/P systematic comparison

2. **M3-A3 alternative detector comparison:**
   - Evaluate alternative multimodality detectors (Hartigan dip test, mixture model-free methods)
   - Compare Product-vM against alternatives under same validation framework
   - Determine if Product-vM convergence issues justify pivot to alternative approach

**Decision:** Requires external review of M3-A2f closure before proceeding.

---

## SCIENTIFIC METHODOLOGY NOTES

### Historical Checkpoint Validity

All checkpoints **before M3-A2d** used an incorrect K=1 fitter and are scientifically invalid for quantitative conclusions. They are retained for provenance and methodological transparency.

### Supersession vs. Supporting

- **Superseded:** Checkpoint conclusions were invalidated by later findings
- **Supporting:** Checkpoint contributes to historical chain but is not standalone conclusion
- **Current:** Represents current accepted scientific state

### Toroidal Independence Assumption

Product von Mises detector assumes:
```
p(θ₁, θ₂, ..., θ_q) = ∏ᵢ₌₁^q p(θᵢ)
```

This is a **strong simplifying assumption**. Real protein dihedral angle distributions may exhibit:
- Angular correlations (e.g., Ramachandran-like coupling)
- Dependent multimodality (modes in joint space, not marginal)

**Validation under dependent-toroidal distributions is MANDATORY before production deployment.**

---

## PRIMARY DELIVERABLES

### Reports
- [M3_A2f_CHECKPOINT_REPORT.md](reports/M3_A2f_CHECKPOINT_REPORT.md) — Final checkpoint with corrected degeneracy proof
- [M3_A2f_CLOSURE.md](reports/M3_A2f_CLOSURE.md) — Checkpoint closure document
- [M3_A2f_LIKELIHOOD_DEGENERACY.md](reports/M3_A2f_LIKELIHOOD_DEGENERACY.md) — Unbounded likelihood proof and literature

### Scripts
- [product_vm_fast_accurate.py](scripts/product_vm_fast_accurate.py) — M3-A2e corrected fitter
- [product_vm_three_estimators.py](scripts/product_vm_three_estimators.py) — U/C/P estimator framework
- [m3_a2f_analysis_framework.py](scripts/m3_a2f_analysis_framework.py) — Corrected collapse analysis tools
- [investigate_wrapped_normal_anomaly.py](scripts/investigate_wrapped_normal_anomaly.py) — Convergence failure investigation

### Results
- [m3_a2e_checkpoint6_exact_family.csv](results/m3_a2e_checkpoint6_exact_family.csv) — 0.13% exact-family K2 selections
- [m3_a2e_checkpoint9_k2_recovery.csv](results/m3_a2e_checkpoint9_k2_recovery.csv) — 100% K2 power (≥60°, q=2)
- [m3_a2e_checkpoint10_wrapped_normal.csv](results/m3_a2e_checkpoint10_wrapped_normal.csv) — Wrapped-normal anomaly

---

## ARTIFACT MANIFEST

See [MANIFEST.md](MANIFEST.md) for complete provenance, checksums, and regenerability status.

---

**NO PRODUCTION CODE MODIFICATIONS APPROVED**

All artifacts in this directory are scientific audit evidence. None have been promoted to `src/` or production TopoFold code.
