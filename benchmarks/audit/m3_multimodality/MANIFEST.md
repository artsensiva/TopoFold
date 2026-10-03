# M3 Multimodality Detection: Artifact Manifest

**Generated:** 2026-10-03 (Corrections Pass - FINAL)
**Repository Branch:** scientific-rebuild
**Repository Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c
**Current Status:** **O2 - FUNCTIONAL BUT FRAGILE** (provisional, diagnostics incomplete)

This manifest records provenance, checksums, and validity status of all M3-A2g scientific audit artifacts.

---

## STATUS AND NAVIGATION DOCUMENTS

| Path | Purpose | SHA-256 |
|------|---------|---------|
| `M3_STATUS.md` | Current M3 status (Provisional O2) | `a959947cc08a96ac71ebc29fff99a658e3b03e56984013f56928ed52e3fd713d` |
| `README.md` | Navigation and limitations | `c4eedfb63dad7faaef729b65238daf6f1e7cfc49c7b87c9b2efb7b4bc5f4b990` |
| `M3_A2G_CORRECTIONS_PROVENANCE.md` | External review findings | `eaed51a2ad72efb17b0014bd0a1b7ad33c5cd71d5e921951eb43173c6447e8fb` |
| `M3_A2G_CHECKPOINT_CLOSURE.md` | M3-A2g closure | `f40b7f561b05ac21a6985a840fb67aa80939d3ccffd9b65871c2e8c90c8a6d40` |
| `M3_A2G_COMPLETE_ANALYSIS.md` | M3-A2g analysis | `d51688da073539167254211c57cd3f98b1cb4f4c1daf3f4b9ff169656fa13c57` |
| `MANIFEST.md` | This provenance manifest | (self) |

---

## M3-A2g CHECKPOINT REPORTS (CURRENT)

| Path | Purpose | Status | SHA-256 |
|------|---------|--------|---------|
| `reports/M3_A2g_CHECKPOINT_REPORT.md` | M3-A2g detailed report | **CORRECTED (Provisional O2)** | `23837c646124bb2b065f24b74c83e1cc438d2c681085808babb3fefd4e7decc0` |
| `M3_A2G_CHECKPOINT_CLOSURE.md` | M3-A2g closure | **CORRECTED (Provisional O2)** | `f40b7f561b05ac21a6985a840fb67aa80939d3ccffd9b65871c2e8c90c8a6d40` |
| `M3_A2G_COMPLETE_ANALYSIS.md` | M3-A2g complete analysis | **CORRECTED (Provisional O2)** | `d51688da073539167254211c57cd3f98b1cb4f4c1daf3f4b9ff169656fa13c57` |
| `M3_A2G_CORRECTIONS_PROVENANCE.md` | External review findings | **CORRECTED** | `eaed51a2ad72efb17b0014bd0a1b7ad33c5cd71d5e921951eb43173c6447e8fb` |

---

## M3-A2g EXPERIMENT SCRIPTS (8 scripts, 1235 CSV rows)

| Script | Purpose | Validity | SHA-256 |
|--------|---------|----------|---------|
| `scripts/m3_a2g_convergence_audit.py` | Convergence trajectory analysis | ✅ VALID | `0e73e3c161f6ea5d9c7d0a78a384753651871c600a3e6be6f8ca47c25dcdbb46` |
| `scripts/m3_a2g_exact_k2_convergence.py` | Exact/K2 convergence (360+320 fits) | ✅ VALID | `bf81b22a4a2cf360b611133b38246ee797d66159daf0cce18d7372c07f310922` |
| `scripts/m3_a2g_monotonicity_check.py` | EM monotonicity (4 fixtures) | ⚠ VALID BUT LIMITED | `fe8745b607b5493060a6ee2a365a306e3d58eb60f8517632f068dd6addee04eb` |
| `scripts/m3_a2g_ninit_sensitivity.py` | n_init sensitivity (20 fits) | ❌ INVALID (inner reseeding) | `81be468a3a5b07817c438b67317f43d6afaac787474336dfc251165c425372a4` |
| `scripts/m3_a2g_stopping_criterion_validation.py` | Stopping validation | ❌ INVALID (ALL protocols) | `ab853368155185561cbfb5f1760c6791af851f5a1ff7f01a7b63651be8b5e9b7` |
| `scripts/m3_a2g_extended_maxiter.py` | Extended max_iter (180 fits) | ⚠ VALID BUT UNPAIRED | `3912af0a32251ec205d44aab59ce09f5e5e7145dd609dcbae51739e1f6175328` |
| `scripts/m3_a2g_duplicate_component_investigation.py` | Duplicate components | ✅ VALID | `131831e214e7d1ea3e240ca391482e5d0fbf1b4806628a570cc3a53147b7db60` |
| `scripts/m3_a2g_data_verification.py` | Data verification (trajectory bug FIXED) | ✅ VALID (corrected) | `a6fd0d49323cf8ccacbb577d2d7f9ca818afb3cefdb5615a0b5b5e417fa0168f` |

---

## M3-A2g RESULTS (10 CSVs, 1 TXT — PRESERVED UNCHANGED; 1 JSON — REGENERATED)

| File | Purpose | Rows/Records | SHA-256 |
|------|---------|--------------|---------|
| `results/m3_a2g_exact_family_convergence.csv` | Exact-family (111/360=30.8%) | 360 | `d72501064cc0622c7ae572d846b4bd77e112caf862f96086cfe0e5f879bd8430` |
| `results/m3_a2g_k2_convergence.csv` | True-K2 (258/320=80.6%) | 320 | `9897f9da806401e937aa5cd6bde8fcff3e01d03fe4bcaa64aa716c58583b23e4` |
| `results/m3_a2g_stress_regime.csv` | Wrapped stress (102/280=36.4%) | 280 | `efe8139cb048107e14f94534c86a54f47d560bd49e727073d1707d31381ba66e` |
| `results/m3_a2g_extended_maxiter.csv` | Extended max_iter (56/60 at 2000) | 180 | `e0fd53fcf5218f6124acf7137de6693a092ea251e0f547932473cecdf8d025ea` |
| `results/m3_a2g_ninit_sensitivity.csv` | n_init (loglik diff=2.63) | 20 | `8e8fffee27b874c2132aa4939f3ab34628675fc6e4fe79e2c20a18df626c0c0f` |
| `results/m3_a2g_stopping_criterion_validation.csv` | Stopping | 20 | `c6b11f820825bca1618290901ea61eda1ee48f4c48e00dd25f49a9f2799874d0` |
| `results/m3_a2g_tolerance_comparison.csv` | Tolerance | 16 | `e62e2280c52bc66f013ac4fb6c7a9915f308a4b0a77efbe9476c12e4a32a9367` |
| `results/m3_a2g_max_iter_sensitivity.csv` | Max_iter | 16 | `500fcd757f49ed0ceea802d6cbc51ccdddad5a97c7ac22fadd4e495f717f73fd` |
| `results/m3_a2g_initialization_comparison.csv` | Initialization | 15 | `9bdb14e031c44605e0b43ae9f9709bb067544df302fa7cfd2472d04b6f5e5882` |
| `results/m3_a2g_duplicate_components.csv` | Duplicates | 8 | `c3fbcc3b12c683619e70f69b3a38369b90daae82f90e5e865acd91c36dfcf8f8` |
| `results/m3_a2g_trajectories.txt` | 5 blocks, 1634 iterations | 1634 records | `fd298b3430602ff82b47b7e396021698ba6c70cd53923e2159f41527dd44e670` |
| `results/m3_a2g_data_verification.json` | **REGENERATED** derivative summary | (summary) | `af95fe40b277fb0f38d67bc50a3bc29e657933461ce2f418e2a428d8609bc5bb` |

**Total CSV rows:** 1235 (verified)

**Note:** All 10 CSVs and trajectories.txt are original experiment outputs preserved unchanged. The verification.json is a SYNTHETIC derivative regenerated by corrected verification utility (trajectory bug fixed: 5→1634).

---

## M3-A2e VALIDATION RESULTS (SUPPORTING)

| File | Purpose | Result | SHA-256 |
|------|---------|--------|---------|
| `results/m3_a2e_checkpoint6_exact_family.csv` | K2-selection under exact-K1 nulls | 3/2250 (0.13%) | `4b4e89e942c3ee58dfcf666d1f62978cfbd6146bc30f26fde3789bcaba20fa31` |
| `results/m3_a2e_checkpoint9_k2_recovery.csv` | Clear-K2 power (sep≥60°, q=2) | 1620/1620 (100%) | `28c73b63927707a8baf3d914409a81bb104ea23525e5ca35a4da6080b7dd1821` |
| `results/m3_a2e_checkpoint10_wrapped_normal.csv` | Wrapped-normal anomaly | (investigation) | `73eaabc4ef5dd830237902c22d20dfc315e1770dd1865ca7c28909d36767b267` |

---

## EARLIER RESULTS (HISTORICAL)

| File | Purpose | SHA-256 |
|------|---------|---------|
| `results/m3_a2c_exact_null_results.csv` | M3-A2c sanity audit | `4d28b21eb2e1d9c3e3e82831f22745718de084b4a170bf4b611d2e98dc4622c1` |

---

## CORE FITTERS AND FRAMEWORKS

| Script | Purpose | Checkpoint | SHA-256 |
|--------|---------|-----------|---------|
| `scripts/product_vm_fast_accurate.py` | M3-A2e corrected K=1/K=2 fitter | M3-A2e | `bb1a12f5b25f996d74e5aeaabb9f107cfe2b329c22391c8ef9e22fcac5260f47` |
| `scripts/product_vm_three_estimators.py` | U/C/P regularization framework | M3-A2f | `aee65e4f262851c72812bbdf811af5cb8be879c90187d0c3ab73bb461ddfac69` |
| `scripts/m3_a2g_data_verification.py` | Data verification (trajectory bug FIXED) | M3-A2g | `a6fd0d49323cf8ccacbb577d2d7f9ca818afb3cefdb5615a0b5b5e417fa0168f` |

---

## HISTORICAL CHECKPOINT REPORTS

| Path | Checkpoint | Status | SHA-256 |
|------|-----------|--------|---------|
| `reports/M3_A2f_CHECKPOINT_REPORT.md` | M3-A2f | Superseded (classification) | `7e2fb3bdb0e8d3fc3fb9f45ab188e879bfce904b81264824281df608d167bb07` |
| `reports/M3_A2f_CLOSURE.md` | M3-A2f | Supporting | `7122cb6783791b7a9cbffa9ca6e90ff5e80b2426094c5c9604397f5b35659ee9` |
| `reports/M3_A2f_LIKELIHOOD_DEGENERACY.md` | M3-A2f | Supporting | `dc249a3eba008bbd2a29131a414eb775663eaf32b29dcf2d1b75aa8326cf5dac` |
| `reports/M3_A2e_CHECKPOINT_REPORT.md` | M3-A2e | Supporting | `9b1d488d6cac511201a41c3fc942e57127cafd9bb8c3c9222002dfec409a5d21` |

**Pre-M3-A2d checkpoints:** Scientifically invalid (incorrect K=1 fitter), retained for provenance only.

---

## REGENERABILITY STATUS

### PRESERVED UNCHANGED (Authoritative)
**All 10 M3-A2g CSVs and trajectories.txt preserved exactly as generated** despite documented diagnostic bugs. Regenerability: **CANNOT REGENERATE** (experimental outputs are provenance artifacts).

### REGENERATED SYNTHETIC DERIVATIVES
**m3_a2g_data_verification.json:** REGENERATED by corrected verification utility (trajectory counting bug fixed: 5→1634). This is a SYNTHETIC derivative summary, not an original experiment artifact. Regenerability: **CAN REGENERATE** (run corrected verification utility on preserved CSVs/TXT).

### DIAGNOSTIC LIMITATIONS DOCUMENTED

1. **m3_a2g_ninit_sensitivity.py:** ❌ INVALID FOR ROBUSTNESS CLAIMS
   - Inner reseeding defect: reseeds on every INNER start, making all n_init>1 identical
   - Output CSV preserved for provenance only
   - Fix required: Fix inner reseeding defect in fitter, then re-execute with independent per-start seeds

2. **m3_a2g_stopping_criterion_validation.py:** ❌ INVALID FOR VALIDATION
   - Relative tolerance branch unimplemented (falls through to tol=1e-20)
   - Absolute tolerance protocols also have issues: approximate stop iteration, incorrect continuation endpoint, non-converged fits
   - ALL stopping protocols tested have limitations
   - Fix required: Implement actual relative tolerance, fix continuation methodology

3. **m3_a2g_extended_maxiter.py:** ⚠ VALID BUT UNPAIRED
   - Different seeds per max_iter value (seed = 1000 + len(results))
   - Not controlled comparison on identical draws
   - Convergence differences partially reflect draw difficulty
   - Better design: seed = base_seed + stress_config_index * 100 + rep

4. **m3_a2g_monotonicity_check.py:** ⚠ VALID BUT LIMITED
   - Only 4 synthetic fixtures tested (exact_easy, wrapped_1.0, k2_clear, k2_close)
   - Each fixture × 3 identical inits (seed=700) = 12 runs
   - 3384 iteration states total (including triplicates)
   - 1124 unique comparisons (after removing triplicates)
   - 0 violations observed
   - Empirical verification on limited fixtures, not comprehensive proof

**See [M3_A2G_CORRECTIONS_PROVENANCE.md](M3_A2G_CORRECTIONS_PROVENANCE.md) for detailed diagnostic analysis.**

---

## CORRECTED VERIFICATION UTILITY

**scripts/m3_a2g_data_verification.py** — Trajectory counting bug **FIXED** (now correctly counts 1634 iteration state records, not 5 block headers)

**Corrections applied by Claude during correction pass:**
- Path resolution from `__file__` (not current directory)
- CSV schema validation
- Grid coverage validation
- Input hash computation
- SYNTHETIC labeling with full provenance metadata

**Verification:** All statistics in corrected documents derived from this utility operating on preserved CSVs/TXT.

---

## PROVISIONAL O2 CLASSIFICATION

**Current:** O2 - FUNCTIONAL BUT FRAGILE (provisional)

**Wording:**
```
U — numerical-ceiling local-EM heuristic.
O2 — FUNCTIONAL BUT FRAGILE, provisional.
Optimization fragility documented;
initialization and stopping validation incomplete;
statistical robustness unresolved.
Not production-ready; no validated production defaults.
```

**Rationale:**
- ✅ Budget exhaustion identified as root cause (not fundamental O1 failure)
- ✅ Limited monotonicity evidence (4 fixtures, 1124 unique comparisons, 0 violations) suggests EM not systematically broken
- ⚠ Convergence remains POOR at max_iter=500 (exact: 30.8%, stress: 36.4%)
- ⚠ Extended budget experiment (max_iter=2000, UNPAIRED draws): 56/60 (93.3%) in stress regime
- ❌ Initialization validation INCOMPLETE (diagnostic invalid)
- ❌ Stopping validation INCOMPLETE (diagnostic invalid)

**NOT O3 because:**
- Initialization and stopping validation remain incomplete
- Diagnostic scripts have critical bugs that invalidate robustness claims
- Convergence rates remain poor at max_iter=500
- No validation under dependent-toroidal distributions

**Required before O3:**
1. Fix initialization diagnostic (fix inner reseeding defect in fitter, re-execute with independent per-start seeds)
2. Fix stopping validation (implement relative tolerance, fix continuation methodology)
3. Execute corrected diagnostics
4. Only then consider upgrading to O3 if results support

**Required before production:**
- Dependent-toroidal validation (MANDATORY)
- U/C/P regularization decision
- Extended dimensional validation (q=1, q=5)
- Temporal correlation sensitivity
- Calibrated significance testing

---

## DOCUMENT CORRECTION STATUS

**✅ Corrected (2026-10-03):**
- M3_STATUS.md (provisional O2, corrected numbers, budgets attached)
- README.md (diagnostic limitations documented)
- M3_A2G_CORRECTIONS_PROVENANCE.md (external review findings, diagnostic bugs)
- M3_A2G_CHECKPOINT_CLOSURE.md (corrected claims history, budget attachment)
- M3_A2G_COMPLETE_ANALYSIS.md (complete convergence data, UNPAIRED disclaimers)
- reports/M3_A2g_CHECKPOINT_REPORT.md (detailed report, corrected statistics, invalid diagnostic sections)
- scripts/m3_a2g_data_verification.py (trajectory bug fixed, schema/grid validation added)
- MANIFEST.md (this document - full SHA-256 hashes, regenerability status)

**All corrections completed.** All documents now consistent with provisional O2 classification, corrected statistics, and documented diagnostic limitations.

---

## PRESERVATION POLICY

**All M3-A2g experimental outputs (10 CSVs, 1 TXT) preserved unchanged** despite diagnostic script bugs. Bugs documented in M3_A2G_CORRECTIONS_PROVENANCE.md and validity status marked in this manifest.

**Regenerated derivative:** m3_a2g_data_verification.json is a SYNTHETIC summary regenerated by corrected verification utility (not original experiment artifact).

**No production code modifications.** All artifacts are scientific audit evidence only.

---

**Repository Baseline:** 37ce1350e934ea0e54d92c55ad83c9748e4ffc3c
**Corrections Date:** 2026-10-03
**External Review Status:** Provisional O2 pending review
