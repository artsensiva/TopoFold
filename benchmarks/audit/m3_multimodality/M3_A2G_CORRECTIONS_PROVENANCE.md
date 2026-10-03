# M3-A2g CORRECTIONS PROVENANCE

**Document Type:** External Review Findings and Diagnostic Limitations
**Date:** 2026-10-03
**Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c
**Scope:** Documentation corrections and verification tooling execution (corrected trajectory counting bug)

---

## EXTERNAL REVIEW REJECTION

**Initial Claim:** O3 - ROBUST (all directives satisfied)

**External Review Findings:**

1. **Numerical Reconciliation Errors:**
   - Claimed: "Total CSV rows: 1245"
   - Actual: **1235 rows** (authoritative from verification utility)

   - Claimed: "Exact-family F=2500 range: 0-10%"
   - Actual: **0-30%** (verified from CSV: q∈{1,2,5}, κ∈{1,5,20})

   - Claimed: "K2 sep=30° range: 70-100%"
   - Actual: **40-100%** (verified: varies by kappa and occupancy)

   - Claimed: "Stress σ=1.0, F=1000: 0/10"
   - Actual: **3/10** (verified from m3_a2g_stress_regime.csv)

   - Claimed: "Iteration records: 5"
   - Actual: **1634 iteration state records** (bug in verification utility, now fixed)

2. **Historical Attribution Errors:**
   - Incorrectly claimed M3-A2b "fixed K1 fitter" — requires verification
   - Incorrectly claimed M3-A2d "validated K1 with 1620 fits" — actually M3-A2e
   - Misattributed 3/2250 as "close-K2 recovery" — actually K2-selection under exact-K1 nulls (M3-A2e)

3. **Diagnostic Script Limitations (CRITICAL):**

**External Review Decision:** **O3 classification REJECTED as premature**

---

## DIAGNOSTIC SCRIPT LIMITATIONS

### 1. Initialization Sensitivity — INVALID FOR ROBUSTNESS CLAIMS

**Script:** [scripts/m3_a2g_ninit_sensitivity.py](scripts/m3_a2g_ninit_sensitivity.py)

**Bug:** Reseeding on every INNER start

**Code Analysis:**
```python
# Line 74-75 in test_ninit_sensitivity():
result = fit_K2_with_full_diagnostics(
    samples, initialization='split' if n_init == 1 else 'random',
    n_init=n_init, max_iter=1000, tol_abs=1e-6, seed=800  # <-- FIXED SEED
)
```

**Issue:** The fitter reseeds on every INNER start, making all n_init>1 produce **identical** random initializations. When n_init=3, all 3 random starts get the same seed, generating the same initial parameters. Only n_init=1 differs because it uses split initialization instead.

**Implication:** The observed wrapped_1.0 loglik diff = 2.63 between n_init=1 (split) and n_init=3 (all-random) indicates **DIFFERENT SOLUTIONS**, which could reflect either:
- True local optima sensitivity (concerning)
- Artifact of split vs random initialization comparison
- Combination of both

**Validity:** **INVALID FOR ROBUSTNESS CLAIMS** — preserved for provenance only

**Required Fix:** Fix inner reseeding defect in the fitter itself, then re-execute with independent per-start seeds. Merely changing the outer seed to `seed=800+n_init` does NOT fix this defect.

---

### 2. Stopping Criterion Validation — INVALID FOR VALIDATION

**Script:** [scripts/m3_a2g_stopping_criterion_validation.py](scripts/m3_a2g_stopping_criterion_validation.py)

**Bug:** Relative tolerance branch unimplemented; absolute protocols also have issues

**Code Analysis:**
```python
# Line 69-73 in continuation_test():
if tol_rel is not None:
    # Re-run with monitoring to find when relative criterion would fire
    pass  # For simplicity, use tol_abs as proxy
```

**Issue:** When `tol_rel` is specified, the function does not implement relative tolerance logic at all — it falls through to using the absolute criterion with tol=1e-20 (effectively no stopping).

**Configs Affected:**
- `rel_1e-6`, `rel_1e-7`, `rel_1e-8` all use tol=1e-20, NOT relative tolerance
- `abs_1e-6` and other absolute configs also have issues: approximate stop iteration, incorrect continuation endpoint, non-converged fits

**Implication:** ALL stopping protocols tested have limitations. The relative tolerance comparison is **unexecuted**, and the absolute tolerance methodology has multiple issues.

**Validity:** **INVALID FOR VALIDATION** — preserved for provenance only

**Required Fix:** Implement actual relative tolerance criterion `|Δloglik| / |loglik| < tol_rel`, fix continuation test methodology (correct stop iteration tracking, proper continuation endpoint, handle non-converged fits)

---

### 3. Monotonicity Verification — LIMITED SCOPE

**Script:** [scripts/m3_a2g_monotonicity_check.py](scripts/m3_a2g_monotonicity_check.py)

**Limitation:** Only 4 synthetic fixtures tested:
- exact_easy
- wrapped_1.0
- k2_clear
- k2_close

**Design:** Each fixture × 3 identical inits (seed=700) = 12 runs

**Scope:** Empirical verification on limited fixtures, NOT comprehensive mathematical proof

**Findings:**
- 3384 iteration states total (including triplicates)
- 3372 comparisons including repeats
- 1124 unique comparisons (after removing triplicates)
- 0 monotonicity violations

**Validity:** ✅ **VALID BUT LIMITED** — suggests EM not systematically broken but does not constitute proof of correctness

---

### 4. Extended max_iter Test — UNPAIRED DRAWS

**Script:** [scripts/m3_a2g_extended_maxiter.py](scripts/m3_a2g_extended_maxiter.py)

**Limitation:** Different seeds per (sigma, F, max_iter, rep) configuration

**Code Analysis:**
```python
# Line 56:
seed = 1000 + len(results)  # Increments across ALL prior results
```

**Issue:** Each max_iter value gets **different random draws**, so the comparison is not on identical data. For example:
- σ=1.0, F=1000, max_iter=500, rep=0: seed=1000
- σ=1.0, F=1000, max_iter=1000, rep=0: seed=1010 (different!)
- σ=1.0, F=1000, max_iter=2000, rep=0: seed=1020 (different!)

**Implication:** Convergence rate differences partially reflect **draw difficulty**, not purely max_iter sufficiency.

**Validity:** ⚠ **VALID BUT UNPAIRED** — shows max_iter trends but not controlled comparison on identical draws

**Better Design:** `seed = base_seed + stress_config_index * 100 + rep` to ensure paired draws

---

## CORRECTED AUTHORITATIVE STATISTICS

**Source:** [scripts/m3_a2g_data_verification.py](scripts/m3_a2g_data_verification.py) (trajectory counting bug fixed)

### CSV Data
- **Total CSV rows:** 1235 (not 1245)
- **Trajectory records:** 1634 iteration state records (not 5 headers)

### Exact-Family Convergence (360 fits)
- **Total:** 111/360 (30.8%)
- **By F:**
  - F=100: 75/90 (83%)
  - F=500: 20/90 (22%)
  - F=1000: 9/90 (10%)
  - F=2500: 7/90 (8%)
- **Configs <20%:** 20/36
- **Configs at 0%:** 9/36

**F=2500 range across q∈{1,2,5}, κ∈{1,5,20}:** 0-30% (not "0-10%")

### True-K2 Convergence (320 fits)
- **Total:** 258/320 (80.6%)
- **By separation:**
  - sep=10°: 27/80 (34%)
  - sep=30°: 71/80 (89%)
  - sep=60°: 80/80 (100%)
  - sep=90°: 80/80 (100%)

**sep=30° range across κ, occupancy:** 40-100% (not "70-100%")

### Stress Regime (280 fits)
- **Total:** 102/280 (36.4%)
- **Key examples:**
  - σ=1.0, F=1000: **3/10** (not 0/10)
  - σ=1.0, F=2500: 1/10
  - σ=0.6, F=2500: 0/10

### n_init Sensitivity
- **wrapped_1.0 loglik diff:** 2.6307315391 (between n_init=1 (split) and n_init=3 (all-random))
- **Interpretation:** Different solutions reached (NOT validation of robustness due to inner reseeding bug)

### Extended max_iter (180 fits, unpaired draws)
- σ=1.0, F=1000, max_iter=2000: 10/10 (100%)
- σ=1.0, F=2500, max_iter=2000: 10/10 (100%)

---

## CORRECTED HISTORICAL ATTRIBUTION

### M3-A2d (Core corrections following A2c)
**Documented:** K=1 fitter corrections identified in M3-A2c
**Deferred:** Full benchmark and recovery validation work

### M3-A2e (Full benchmark with corrected fitter)
**Validated:**
- **3/2250 (0.13%)** K=2-selection frequency under exact-K1 nulls (NOT "close-K2 recovery")
- **1620/1620 (100%)** clear-K2 power (sep≥60°, q=2 only)

**Source Files:**
- [results/m3_a2e_checkpoint6_exact_family.csv](results/m3_a2e_checkpoint6_exact_family.csv)
- [results/m3_a2e_checkpoint9_k2_recovery.csv](results/m3_a2e_checkpoint9_k2_recovery.csv)

---

## PROVISIONAL O2 CLASSIFICATION

**Corrected Classification:** **O2 - FUNCTIONAL BUT FRAGILE** (provisional working classification)

**Rationale:**
1. ✅ Budget exhaustion identified as root cause of "97% non-convergence" (not fundamental failure)
2. ✅ Limited monotonicity evidence (4 fixtures, 1124 unique comparisons, 0 violations) suggests EM not systematically broken
3. ⚠ Convergence remains POOR at max_iter=500 (exact: 30.8%, stress: 36.4%)
4. ⚠ Extended budget experiment (max_iter=2000, UNPAIRED draws): 56/60 (93.3%) in stress regime
5. ❌ Initialization robustness **INCOMPLETE** (inner reseeding bug invalidates diagnostic)
6. ❌ Stopping criterion **INCOMPLETE** (relative branch unimplemented, absolute protocols have issues)

**Not O3 because:**
- Initialization and stopping validation remain incomplete
- Diagnostic scripts have critical bugs that invalidate robustness claims
- Convergence rates remain poor at max_iter=500
- No validation under dependent-toroidal distributions

**Next Required Steps:**
1. Fix initialization diagnostic (fix inner reseeding defect, re-execute with independent per-start seeds)
2. Fix stopping validation (implement relative tolerance, fix continuation methodology)
3. Execute corrected diagnostics
4. Only then consider upgrading to O3 if results support

---

## PRESERVATION POLICY

**Original M3-A2g experiment artifacts PRESERVED UNCHANGED** despite diagnostic limitations:

**Preserved Scripts (with documented bugs):**
- [m3_a2g_ninit_sensitivity.py](scripts/m3_a2g_ninit_sensitivity.py) — inner reseeding bug
- [m3_a2g_stopping_criterion_validation.py](scripts/m3_a2g_stopping_criterion_validation.py) — relative branch unimplemented, absolute protocols have issues
- [m3_a2g_extended_maxiter.py](scripts/m3_a2g_extended_maxiter.py) — unpaired draws
- [m3_a2g_monotonicity_check.py](scripts/m3_a2g_monotonicity_check.py) — limited scope (4 fixtures)

**Preserved Results:**
- All 10 CSVs (1235 rows total) **UNCHANGED**
- [m3_a2g_trajectories.txt](results/m3_a2g_trajectories.txt) (1634 iteration records) **UNCHANGED**

**Corrected Tooling (Claude executed during correction pass):**
- [m3_a2g_data_verification.py](scripts/m3_a2g_data_verification.py) — trajectory counting bug **FIXED** (5→1634)
- [m3_a2g_data_verification.json](results/m3_a2g_data_verification.json) — **REGENERATED** derivative summary (not original experiment artifact)

**Rationale:** Raw experimental outputs are scientific provenance. Diagnostic bugs are documented in this provenance file and in corrected reports, but original artifacts remain as generated. The verification utility was corrected and re-executed to produce accurate summary statistics.

---

## POST-COMMIT CORRECTION: HISTORICAL DIAGNOSTIC MISINTERPRETATION

**Date:** 2026-10-03 (post-commit bf33346e52676f87925643681f3669879ef5fc18)
**Issue:** Erroneous interpretation of `/tmp/m3_a2_fast_validation.py` results during background process review

### WITHDRAWN CLAIMS

The following claims made during background process review are **WITHDRAWN**:

1. ❌ **"temporal correlation is fatal"** - WITHDRAWN
2. ❌ **"never a primary detector"** - WITHDRAWN
3. ❌ **"publication-quality evidence"** - WITHDRAWN
4. ❌ **"O2 final"** - WITHDRAWN

### CORRECT INTERPRETATION

**Historical M3-A2 evidence from `/tmp/` artifacts:**

The `/tmp/m3_a2_fast_validation.py` script and `/tmp/m3_a2_fast_correlated_nulls.csv` results are **historical M3-A2 evidence** from an old implementation with documented bugs, NOT new evidence about the corrected A2e/A2g implementation.

**Historical diagnostic bugs in `/tmp/m3_a_diagnostic.py`:**

1. **Likelihood calculation error (lines 528-545):**
   ```python
   # E-step: compute responsibilities
   log_resp = np.zeros((n, K))
   for k in range(K):
       for i in range(n):
           log_resp[i, k] = (np.log(weights[k] + 1e-10) +
                            product_von_mises_logpdf(samples[i], mus[k], kappas[k]))

   # Normalize responsibilities (numerical stability)
   log_resp_max = log_resp.max(axis=1, keepdims=True)
   resp = np.exp(log_resp - log_resp_max)
   resp /= resp.sum(axis=1, keepdims=True)  # ← Line 542: normalize to sum=1

   # Compute log-likelihood
   loglik = (log_resp_max.squeeze() + np.log(resp.sum(axis=1))).sum()  # ← Line 545: BUG
   ```

   **Bug:** After normalizing responsibilities to sum=1 (line 542), the log-likelihood calculation computes `log(resp.sum(axis=1))` which always gives log(1)=0, **losing the mixture normalization term** log(Σₖ πₖ p(x|θₖ)).

2. **Concentration approximation (lines 446-452):**
   ```python
   if R < 1e-6:
       kappa[dim] = 0.0
   elif R < 0.8:
       kappa[dim] = 2*R + R**3 + 5*R**5/6  # Taylor-like approximation
   else:
       # High R: use k ≈ 1/(1-R)
       kappa[dim] = 1.0 / (1.0 - R + 1e-10)
   ```

**CSV evidence from `/tmp/m3_a2_fast_correlated_nulls.csv`:**

Common_latent corr=0.7 results (lines 272-301):
- **30/30 K2 selections** (all delta_bic < 0, fp=True)
- **Only 5/30 converged runs** (reps 8, 13, 15, 23, 27)
- **25/30 non-converged** fits still produced K2 selections - invalid optimization

**Correct understanding:**
- The **common_latent parameter 0.7 is NOT a measured correlation coefficient**
- The generator **introduces dependence across coordinates within observations**, NOT temporal autocorrelation
- **Joint unimodality was NOT established** for these dependent-toroidal nulls
- The 30/30 K2 selections **cannot be interpreted as 100% FPR** when only 5/30 runs converged
- This historical evidence does **NOT** apply to the corrected M3-A2e/A2g implementation

**Preservation:** Historical `/tmp/` artifacts preserved unchanged for provenance. No modifications to original M3-A2g experiment artifacts.

---

## DOCUMENT STATUS

This corrections provenance document was created following external review rejection of the initial O3 classification.

**Corrections Applied:**
- ✅ Numerical reconciliation errors corrected
- ✅ Historical attribution corrected (M3-A2d/A2e)
- ✅ Diagnostic limitations documented
- ✅ Classification downgraded to provisional O2
- ✅ Verification utility trajectory bug fixed (5→1634)
- ✅ All canonical documents updated with corrected content
- ✅ **POST-COMMIT:** Historical diagnostic misinterpretation withdrawn (4 false claims)

**Scope:** Documentation corrections and verification tooling execution. Claude corrected the trajectory counting bug in m3_a2g_data_verification.py and regenerated m3_a2g_data_verification.json. Post-commit correction: Claude withdrew erroneous claims about historical `/tmp/` validation results and documented bugs in historical diagnostic implementation.

**Verification:** All statistics verified from preserved CSV/TXT artifacts using corrected verification utility.

**Next Action:** External review of corrected documentation package before proceeding to M3-A2h or alternative checkpoints.

---

**Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c
**Corrections Date:** 2026-10-03
**External Review Required:** Yes (provisional O2 pending validation completion)
