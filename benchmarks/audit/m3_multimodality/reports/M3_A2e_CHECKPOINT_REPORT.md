# M3-A2e CHECKPOINT REPORT: Full Corrected Mixture Validation

**SYNTHESIZED ARCHIVAL REPORT**

Created during artifact promotion from preserved benchmark script and result CSVs.
Not an original contemporaneous checkpoint file.

**Date:** 2026-10-03 (promotion date)
**Phase:** M3-A (Product von Mises Validation)
**Checkpoint:** M3-A2e (Full Corrected Validation with Benchmark Results)
**Status:** COMPLETE
**Repository HEAD at analysis:** not recorded
**Superseded by:** M3-A2f (likelihood degeneracy and convergence analysis)

**Authoritative quantitative sources:**
- `scripts/m3_a2e_full_benchmark.py` — Benchmark orchestration
- `results/m3_a2e_checkpoint6_exact_family.csv` — Exact-family validation
- `results/m3_a2e_checkpoint9_k2_recovery.csv` — K2 power validation
- `results/m3_a2e_checkpoint10_wrapped_normal.csv` — Wrapped-normal sensitivity

**Original contemporaneous report:** not available / not preserved as a final checkpoint report

---

## SUPERSEDED SCIENTIFIC CONCLUSION

**This report is retained for provenance.**

The algorithmic validation results in this checkpoint remain valid:
- Exact-family validation: 0.13% K2 selections (3/2250)
- K2 power: 100% recovery at separations ≥60° (q=2, 1620/1620)
- Fast accurate implementation validated

However, the production-readiness interpretation was invalidated by M3-A2f findings:

**Reason:**
1. Von Mises mixture likelihoods are unbounded when κ unrestricted
2. Implicit κ≤10⁶ numerical ceiling is a heuristic, not scientifically justified estimator
3. 97% EM convergence failure at wrapped-normal σ≈1.0 transition
4. Optimizer robustness unresolved

**Current status:** see `../M3_STATUS.md`

---

## EXECUTIVE SUMMARY

M3-A2e validated the corrected Product von Mises fitter with comprehensive benchmarks:

### Validation Results

| Validation | Result | Status |
|-----------|--------|--------|
| Exact Product-vM K=1 nulls | 3/2250 = 0.13% K2 selections | ✅ Excellent |
| K=2 power (q=2, sep≥60°) | 1620/1620 = 100% recovery | ✅ Excellent |
| Wrapped-normal σ=1.0 | 15/30 = 50% K2 selections | ⚠️ Anomalous |

### Implementation Status

| Component | Status | File |
|-----------|--------|------|
| Fast Newton-refined κ MLE | ✅ COMPLETE | product_vm_fast_accurate.py |
| Vectorized K=2 EM | ✅ COMPLETE | product_vm_fast_accurate.py |
| Scipy validation | ✅ COMPLETE | product_vm_fast_accurate.py |
| Canonical ΔBIC convention | ✅ COMPLETE | product_vm_fast_accurate.py |
| EM monotonicity tracking | ✅ COMPLETE | product_vm_fast_accurate.py |

---

## CHECKPOINT 1: CORRECTED ΔBIC DOCUMENTATION ✅

### Standardized Project-Wide Convention

**Definition:**
```python
delta_bic = BIC_K2 - BIC_K1
```

**Interpretation:**
- `delta_bic < 0` → K=2 preferred (BIC_K2 < BIC_K1, lower is better)
- `delta_bic > 0` → K=1 preferred (BIC_K1 < BIC_K2, lower is better)

**Selection:**
```python
selects_K2 = (delta_bic < 0)
```

### Corrected Historical Documentation

**M3-A2c Finding:**
> The historical K1 implementation produced incorrect κ/log-likelihood values. A 345-unit log-likelihood discrepancy was observed. The exact source within the obsolete implementation was not fully isolated before replacement. The corrected implementation has been independently validated against scipy.stats.vonmises.

**NOT claimed:** "The error was due to approximation formula" (this was not established for κ≈5 where direct I0 does not overflow and approximation is already accurate).

**Consequence:** Quantitative Product-vM fitter/model-selection conclusions from checkpoints using the obsolete fitter are superseded. M3-A2e is the first valid benchmark using corrected fitting.

---

## CHECKPOINT 6: EXACT PRODUCT-VONMISES K=1 NULL VALIDATION ✅

**Objective:** Validate algorithmic behavior under the exact generative model (Product-vM K=1).

**Design:**
- Stratified by dimensionality: q ∈ {1, 2, 5}
- Concentration regimes: κ ∈ {1, 5, 20}
- Sample sizes: F ∈ {100, 250, 500, 1000, 2500}
- Replications: 50 per (q, κ, F) configuration
- **Total: 3 × 3 × 5 × 50 = 2250 independent trials**

### Results

**Overall:** 3/2250 K=2 selections = **0.13%**

**Stratification:**
| q | K2 Selections | Total | Frequency |
|---|--------------|-------|-----------|
| 1 | 3 | 750 | 0.40% |
| 2 | 0 | 750 | 0.00% |
| 5 | 0 | 750 | 0.00% |

**Interpretation:**
- Excellent algorithmic specificity under exact model
- All K2 selections occurred at q=1 (univariate)
- No K2 selections for multivariate cases (q=2, q=5)
- Conservative detector: empirical K2-selection frequency under exact Product-vM K=1 is 0.13%

**Note:** This is **not** a calibrated false-positive rate or false discovery rate. It is the empirical frequency of BIC K=2 selection when the data were generated from the exact Product-vM K=1 model.

**Result CSV:** `m3_a2e_checkpoint6_exact_family.csv`

---

## CHECKPOINT 9: K=2 RECOVERY POWER ✅

**Objective:** Validate K=2 detection power for clear bimodal Product-vM distributions.

**Design:**
- q=2 only (bivariate toroidal)
- Angular separations: {10°, 20°, 30°, 60°, 90°, 180°}
- Mixture occupancies: {(0.5, 0.5), (0.8, 0.2)}
- Concentration values: κ ∈ {5, 20, 100}
- Sample sizes: F ∈ {100, 500, 1000}
- Replications: 30 per (sep, occ, κ, F) configuration
- **Total: 6 × 2 × 3 × 3 × 30 = 3240 trials**

### Results

**Clear Separations (≥60°):** 1620/1620 = **100.0% K=2 recovery**

Summary by separation:

| Separation | Configs Tested | Total Trials | K2 Recovered | Power |
|-----------|---------------|--------------|--------------|-------|
| 60° | 18 (2 occ × 3 κ × 3 F) | 540 | 540 | 100% |
| 90° | 18 (2 occ × 3 κ × 3 F) | 540 | 540 | 100% |
| 180° | 18 (2 occ × 3 κ × 3 F) | 540 | 540 | 100% |
| **≥60° Combined** | **54** | **1620** | **1620** | **100%** |

**Close Separations (<60°):** Variable power, depends on separation, sample size, concentration, and occupancy

| Separation | Configs Tested | Total Trials | Representative Power Range |
|-----------|---------------|--------------|----------------------------|
| 10° | 18 | 540 | 0-70% (varies by F, κ, occ) |
| 20° | 18 | 540 | 10-100% (varies by F, κ, occ) |
| 30° | 18 | 540 | 50-100% (varies by F, κ, occ) |

**Interpretation:**
- Perfect power for well-separated modes (≥60°) across all occupancies, concentrations, and sample sizes
- Systematic power increase with sample size at close separations
- No spurious failures at wide separations
- Detector is reliable for clearly bimodal data (q=2)

**Result CSV:** `m3_a2e_checkpoint9_k2_recovery.csv`

---

## CHECKPOINT 10: WRAPPED-NORMAL SHAPE SENSITIVITY ⚠️

**Objective:** Test behavior under shape misspecification (wrapped-normal nulls).

**Design:**
- q=2 bivariate toroidal
- Wrapped-normal standard deviations: σ ∈ {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0}
- Sample size: F=1000
- Replications: 3 per σ
- **Total: 10 × 3 = 30 trials**

### Results

**Anomalous Behavior at σ=1.0:**

| σ Range | K2 Selections | Total | Frequency |
|---------|--------------|-------|-----------|
| 0.4–0.9 | 0 | 27 | 0% |
| **1.0** | **15** | **3** | **Note: CSV reports 15 K2 from 3 trials** |
| 1.2–2.0 | 0 | 0 | — |

**Note on σ=1.0 data:** The preserved CSV reports 15 K2 selections from 3 trials, which appears to be a formatting or aggregation artifact. The M3-A2f investigation (scripts/investigate_wrapped_normal_anomaly.py) clarified the actual finding:
- High-replication retest at σ=1.0: ~33-50% K2 selection frequency
- **Root cause:** 97% EM convergence failure (not shape sensitivity)
- Most K2 selections are non-converged local optima, not valid detector decisions

**Interpretation (M3-A2e at time of checkpoint):**
- Increased K2 selections at intermediate σ values suggested potential shape-misspecification sensitivity
- Requires further investigation (addressed in M3-A2f)

**Result CSV:** `m3_a2e_checkpoint10_wrapped_normal.csv`

---

## M3-A2e CLASSIFICATION (Superseded by M3-A2f)

**Original M3-A2e interpretation:** The detector showed excellent exact-family validation (0.13%) and perfect K2 power (100% ≥60°), suggesting algorithmic reliability and potential production readiness.

**M3-A2f invalidation:**
1. Unbounded likelihood requires explicit regularization strategy
2. 97% convergence failure at wrapped-normal σ≈1.0
3. Optimization robustness unresolved
4. NOT production-ready
5. NOT validated screening detector
6. NOT primary detector candidate

---

## PRIMARY SCRIPTS

- `product_vm_fast_accurate.py` — Corrected fitter implementation (U estimator)
- `m3_a2e_full_benchmark.py` — Benchmark orchestration

---

## PRIMARY RESULT FILES

- `m3_a2e_checkpoint6_exact_family.csv` — Exact Product-vM K=1 nulls (3/2250 = 0.13%)
- `m3_a2e_checkpoint9_k2_recovery.csv` — K2 power validation (1620/1620 = 100% for ≥60°)
- `m3_a2e_checkpoint10_wrapped_normal.csv` — Wrapped-normal sensitivity (anomaly at σ=1.0)

---

**Current status:** see `../M3_STATUS.md`
**Final classification:** see `M3_A2f_CHECKPOINT_REPORT.md`
