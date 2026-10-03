# M3 Multimodality Detection Audit

**Scientific Audit:** Product von Mises (Product-vM) mixture multimodality detector
**Current Status:** **O2 - FUNCTIONAL BUT FRAGILE** (provisional; initialization and stopping validation incomplete)
**Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c (scientific-rebuild branch)
**Last Updated:** 2026-10-03

---

## QUICK START

**Primary Navigation:**
- **[M3_STATUS.md](M3_STATUS.md)** — Current scientific status, validated/unresolved components, next blockers
- **[MANIFEST.md](MANIFEST.md)** — Complete artifact inventory, checksums, regenerability status

**Current Checkpoint:**
- **[M3_A2G_CHECKPOINT_CLOSURE.md](M3_A2G_CHECKPOINT_CLOSURE.md)** — M3-A2g closure summary
- **[M3_A2G_COMPLETE_ANALYSIS.md](M3_A2G_COMPLETE_ANALYSIS.md)** — Complete M3-A2g analysis
- **[M3_A2G_CORRECTIONS_PROVENANCE.md](M3_A2G_CORRECTIONS_PROVENANCE.md)** — External review findings, diagnostic limitations
- **[reports/M3_A2g_CHECKPOINT_REPORT.md](reports/M3_A2g_CHECKPOINT_REPORT.md)** — Detailed M3-A2g checkpoint report

---

## CURRENT CLASSIFICATION

```
U — numerical-ceiling local-EM heuristic.
O2 — FUNCTIONAL BUT FRAGILE, provisional.
Optimization fragility documented;
initialization and stopping validation incomplete;
statistical robustness unresolved.
Not production-ready; no validated production defaults.
```

**NOT production-ready. NOT a validated screening detector. NOT a primary detector candidate.**

---

## KEY FINDINGS

### ✅ Validated (M3-A2e)
- 3/2250 (0.13%) K=2-selection under exact-K1 nulls — excellent algorithmic specificity
- 1620/1620 (100%) clear-K2 recovery (sep≥60°, q=2 only) — perfect power on well-separated modes

### ⚠ Provisional O2 (M3-A2g)
- Budget exhaustion identified as root cause of "97% non-convergence" (NOT fundamental optimizer failure)
- Limited monotonicity evidence (4 fixtures, 1124 unique comparisons, 0 violations) suggests EM not systematically broken
- Convergence remains POOR at max_iter=500 (exact: 30.8%, stress: 36.4%)
- Extended budget experiment (max_iter=2000, UNPAIRED draws): 56/60 (93.3%) in stress regime
- **Initialization and stopping validation INCOMPLETE** (diagnostic bugs documented)

### ❌ Unresolved / Critical Gaps
- Dependent-toroidal robustness: **NOT TESTED** (mandatory before production)
- Unbounded likelihood regularization: U/C/P comparison incomplete
- Initialization robustness: diagnostic has seed reuse bug (**INVALID**)
- Stopping criterion: diagnostic has unimplemented relative tolerance (**INVALID**)
- Extended dimensional validation (q=1, q=5): **NOT TESTED**
- Temporal correlation sensitivity: **NOT TESTED**
- Calibrated significance testing: **NOT ESTABLISHED**

---

## DIAGNOSTIC LIMITATIONS (CRITICAL)

**M3-A2g diagnostic scripts preserved with documented bugs:**

1. **[m3_a2g_ninit_sensitivity.py](scripts/m3_a2g_ninit_sensitivity.py)** — INVALID FOR ROBUSTNESS CLAIMS
   - Seed reuse bug: Reseeds on every INNER start, making all n_init>1 produce identical random initializations
   - wrapped_1.0 loglik diff = 2.63 between n_init=1 (split) and n_init=3 (all-random) indicates different solutions
   - Merely changing outer seed does NOT fix this defect

2. **[m3_a2g_stopping_criterion_validation.py](scripts/m3_a2g_stopping_criterion_validation.py)** — INVALID FOR VALIDATION
   - Relative tolerance branch unimplemented, falls through to tol=1e-20
   - Absolute tolerance configurations also have issues: approximate stop iteration, incorrect continuation endpoint, non-converged fits
   - ALL stopping protocols tested have limitations

3. **[m3_a2g_extended_maxiter.py](scripts/m3_a2g_extended_maxiter.py)** — VALID BUT UNPAIRED
   - Different seeds per max_iter value, not controlled comparison
   - Convergence differences partially reflect draw difficulty

4. **[m3_a2g_monotonicity_check.py](scripts/m3_a2g_monotonicity_check.py)** — VALID BUT LIMITED
   - Only 4 synthetic fixtures (exact_easy, wrapped_1.0, k2_clear, k2_close)
   - Each fixture × 3 identical inits (seed=700) = 12 runs
   - 3384 iteration states total (including triplicates)
   - 1124 unique comparisons (after removing triplicates)
   - 0 violations observed
   - Empirical verification, not comprehensive proof

**See [M3_A2G_CORRECTIONS_PROVENANCE.md](M3_A2G_CORRECTIONS_PROVENANCE.md) for detailed diagnostic analysis.**

---

## AUTHORITATIVE STATISTICS

**Source:** [scripts/m3_a2g_data_verification.py](scripts/m3_a2g_data_verification.py) (trajectory bug fixed)

- **Total CSV rows:** 1235
- **Trajectory records:** 1634 iteration state records
- **Exact-family convergence:** 111/360 (30.8%), varies 8%-83% by F
- **True-K2 convergence:** 258/320 (80.6%), varies 34%-100% by separation
- **Stress regime convergence:** 102/280 (36.4%)

**See [results/m3_a2g_data_verification.json](results/m3_a2g_data_verification.json) for complete statistics.**

---

## DIRECTORY STRUCTURE

```
benchmarks/audit/m3_multimodality/
├── M3_STATUS.md                           # Current status (PRIMARY)
├── README.md                              # This file
├── MANIFEST.md                            # Artifact inventory
├── M3_A2G_CHECKPOINT_CLOSURE.md           # M3-A2g closure
├── M3_A2G_COMPLETE_ANALYSIS.md            # M3-A2g complete analysis
├── M3_A2G_CORRECTIONS_PROVENANCE.md       # External review findings
├── reports/
│   ├── M3_A2g_CHECKPOINT_REPORT.md        # M3-A2g detailed report
│   ├── M3_A2f_CHECKPOINT_REPORT.md
│   ├── M3_A2f_CLOSURE.md
│   ├── M3_A2f_LIKELIHOOD_DEGENERACY.md
│   ├── M3_A2e_CHECKPOINT_REPORT.md
│   └── [earlier checkpoints...]
├── scripts/
│   ├── m3_a2g_convergence_audit.py        # Convergence trajectory analysis
│   ├── m3_a2g_exact_k2_convergence.py     # Exact/K2 convergence audit
│   ├── m3_a2g_monotonicity_check.py       # EM monotonicity (limited scope)
│   ├── m3_a2g_ninit_sensitivity.py        # n_init test (INVALID - seed bug)
│   ├── m3_a2g_stopping_criterion_validation.py  # Stopping test (INVALID - rel branch)
│   ├── m3_a2g_extended_maxiter.py         # Extended max_iter (unpaired)
│   ├── m3_a2g_duplicate_component_investigation.py
│   ├── m3_a2g_data_verification.py        # Data verification (trajectory bug FIXED)
│   ├── product_vm_fast_accurate.py        # M3-A2e corrected fitter
│   └── product_vm_three_estimators.py     # U/C/P framework
└── results/
    ├── m3_a2g_exact_family_convergence.csv      # 360 fits (111/360 = 30.8%)
    ├── m3_a2g_k2_convergence.csv                # 320 fits (258/320 = 80.6%)
    ├── m3_a2g_stress_regime.csv                 # 280 fits (102/280 = 36.4%)
    ├── m3_a2g_trajectories.txt                  # 5 blocks, 1634 iteration records
    ├── m3_a2g_ninit_sensitivity.csv             # 20 fits (loglik diff = 2.63)
    ├── [6 more M3-A2g CSVs...]
    ├── m3_a2e_checkpoint6_exact_family.csv      # 3/2250 (M3-A2e)
    ├── m3_a2e_checkpoint9_k2_recovery.csv       # 1620/1620 (M3-A2e)
    └── m3_a2g_data_verification.json            # Authoritative statistics
```

---

## HISTORICAL CONTEXT

**Valid Checkpoints:**
- **M3-A2d:** Core corrections following M3-A2c sanity audit
- **M3-A2e:** Full benchmark with corrected fitter (3/2250, 1620/1620 validated)
- **M3-A2f:** Likelihood degeneracy analysis, convergence investigation
- **M3-A2g:** Optimizer convergence audit (provisional O2, diagnostics incomplete)

**Pre-M3-A2d checkpoints used incorrect K=1 fitter** — scientifically invalid for quantitative conclusions, retained for provenance.

**See [M3_STATUS.md](M3_STATUS.md) checkpoint progression table for complete history.**

---

## NEXT STEPS

**Before Production Deployment:**

1. **Fix and re-execute initialization/stopping diagnostics** (M3-A2h candidate)
   - Independent seeds in n_init comparison
   - Implement relative tolerance stopping criterion
   - Only then consider O3 upgrade

2. **Dependent-toroidal validation** (M3-A3 candidate, **MANDATORY**)
   - Bivariate sine von Mises with verified unimodality
   - Correlated wrapped-normal with verified unimodality
   - Measured circular association

3. **Regularization strategy decision** (M3-B candidate)
   - Systematic U/C/P comparison
   - Literature-grounded penalty selection

4. **Extended validation**
   - q=1 (univariate), q=5 (higher-dimensional)
   - Temporal correlation sensitivity
   - Calibrated significance testing

**Alternative:** Pivot to model-free multimodality detectors if misspecification sensitivity is prohibitive.

---

## PRESERVATION POLICY

**Original M3-A2g experiment artifacts preserved unchanged** despite diagnostic limitations:
- Scripts with bugs preserved for provenance (bugs documented in M3_A2G_CORRECTIONS_PROVENANCE.md)
- CSV/TXT results unchanged (1235 rows, 1634 iteration records)
- Corrected verification utility: trajectory counting bug **FIXED** (generates new verification.json)

**No production code modifications.** All artifacts are scientific audit evidence only.

---

## EXTERNAL REVIEW STATUS

**Status:** Provisional O2 classification pending external review

**External Review Rejection (Initial O3 Claim):**
- Numerical reconciliation errors (corrected)
- Historical attribution errors (corrected)
- Diagnostic limitations not disclosed (now documented)
- Classification downgraded to provisional O2

**Required:** External review of corrected documentation package before proceeding.

---

**For questions or external review requests, see project maintainers.**

**Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c
**Last Updated:** 2026-10-03
