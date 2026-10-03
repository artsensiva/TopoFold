# M3-A2g POST-COMMIT CORRECTIONS SUMMARY

**Date:** 2026-10-04
**Session:** Continuation after bf33346e52676f87925643681f3669879ef5fc18
**Scope:** Bounded documentation/provenance correction pass (no production code changes)

---

## COMPLETED CORRECTIONS

### 1. Historical Diagnostic Misinterpretation Withdrawal ✅

**Location:** [M3_A2G_CORRECTIONS_PROVENANCE.md](M3_A2G_CORRECTIONS_PROVENANCE.md)

**Withdrawn claims:**
- ❌ "temporal correlation is fatal"
- ❌ "never a primary detector"
- ❌ "publication-quality evidence"
- ❌ "O2 final"

**Documented bugs in `/tmp/m3_a_diagnostic.py`:**
1. Likelihood calculation error (line 545): `log(resp.sum(axis=1))` computed AFTER normalization to sum=1, always gives log(1)=0, losing mixture normalization term
2. Concentration approximation (lines 446-452): Taylor-like approximation for R<0.8

**CSV evidence verified:**
- `/tmp/m3_a2_fast_correlated_nulls.csv` common_latent corr=0.7: 30/30 K2 selections, only 5/30 converged
- Cannot interpret as "100% FPR" when 25/30 runs failed to converge

### 2. Critical Regression Fixed ✅

**Issue:** Conflated two different fixture sets

**Corrected:**
- **Monotonicity** uses: `exact_easy`, `wrapped_1.0`, `k2_clear`, `k2_close`
  - Uses `initialization='split', n_init=3`: 3 repeated DETERMINISTIC split starts
  - NOT random initialization - deterministic algorithm creates triplicates
- **n_init sensitivity** uses: `exact_k1`, `wrapped_1.0`, `k2_close`, `k2_clear`

**Files corrected:**
- M3_STATUS.md
- README.md
- M3_A2G_COMPLETE_ANALYSIS.md
- reports/M3_A2g_CHECKPOINT_REPORT.md

### 3. Grid Specifications Corrected ✅

**True-K2 grid:**
- Was: F=1000 only
- Now: F∈{250,1000} ✓

**Stress regime grid:**
- Was: σ∈{0.6,0.7,0.8,0.9,1.0} × F∈{100,500,1000,2500} × "14 configs"
- Now: 7 sigma values × 4 sample sizes × 10 reps = 280 fits, q=2 fixed
  - σ∈{0.6,0.7,0.8,0.9,1.0,1.1,1.2}
  - F∈{250,500,1000,2500} ✓

**Files updated:**
- reports/M3_A2g_CHECKPOINT_REPORT.md
- M3_A2G_COMPLETE_ANALYSIS.md

### 4. Split Initialization Description Corrected ✅

**Was:** "3 identical random initializations" or "fitter seed=700 creates triplicates"
**Now:** "Uses `initialization='split', n_init=3` - 3 repeated deterministic split starts (creates triplicates)"

**Distinction preserved:** Deterministic repetition in monotonicity check vs. inner-reseeding defect in random-start n_init diagnostic

---

## REMAINING WORK (INCOMPLETE)

### 5. Additional Documentation Reconciliation ⚠️ PARTIAL

**Still needed in checkpoint report:**
- ✓ Provenance seeds: Example correctly shows 1120/1130/1140
- ❌ Mixture independence statements (need review for "conditional within components")
- ❌ Stopping diagnostic fabricated indices description
- ❌ Remove unsupported "strongest-misspecification" claims
- ❌ Clarify "partially resolved" review attribution
- ❌ Identify inner reseeding in audit wrapper (not production code)
- ❌ Enhance historical withdrawal explanation (combined implementation/null-validity/optimization limits)

### 6. Verification Utility Schema Fix ❌ NOT STARTED

**Current issue:**
- `results/m3_a2g_data_verification.json` shows `schema_validation.exact_family.valid: false`
- Expected columns: `final_loglik`, `n_iter`
- Actual columns: `loglik`, `selected_K`
- Schema mismatch for all 10 CSVs

**Required:**
- Fix expected schema to match actual producer CSVs
- Implement complete checks for all 10 CSVs + trajectory TXT
- Handle NaN values explicitly (7 improvement_after_100 entries)
- Return nonzero exit code on validation failures
- Regenerate JSON with SYNTHETIC labeling

### 7. Acceptance Checklist Replacement ❌ NOT STARTED

**Current issue:**
- Passes 25/25 checks despite schema valid=false
- Accepts "O2" or "provisional" anywhere in files (pattern matching)
- Counts hash strings without verifying targets

**Required:**
- Replace with substantive validation checking structured results
- Verify manifest coverage and actual digests
- Demonstrate invalid inputs fail checks
- No weakening requirements to obtain PASS

### 8. MANIFEST Completion ❌ NOT STARTED

**Current state:**
- 40 hash occurrences covering 36 unique files
- 13 tracked files omitted including:
  - Original convergence-grid script
  - Historical dependencies/reports
  - Correction utilities

**Required:**
- Build complete inventory from audit tree and Git
- Include all scientific artifacts
- Include retained correction utilities
- Purpose-specific limitations for multipurpose scripts
- Verify all SHA-256 digests match current files

### 9. Final Verification ❌ NOT STARTED

**Required checks:**
- Markdown link validity and anchor resolution
- Python syntax correctness
- Original artifact preservation (8 scripts, 10 CSVs, 1 TXT unchanged)
- Git state clean (no staged/uncommitted changes)
- File change inventory
- Limitation documentation

---

## FILES MODIFIED

**Documentation (6 files):**
1. M3_A2G_CORRECTIONS_PROVENANCE.md - Withdrawal section added
2. M3_STATUS.md - Fixture names, grid specs corrected
3. README.md - Monotonicity description, fixtures corrected
4. M3_A2G_COMPLETE_ANALYSIS.md - Grid specs, fixtures corrected
5. M3_A2G_CHECKPOINT_CLOSURE.md - Diagnostic descriptions clarified
6. reports/M3_A2g_CHECKPOINT_REPORT.md - Grid specs, fixtures, initialization corrected

**Unchanged (as required):**
- 8 experimental scripts (m3_a2g_*.py)
- 10 CSVs (m3_a2g_*.csv)
- 1 trajectory TXT (m3_a2g_trajectories.txt)
- 1 verification JSON (m3_a2g_data_verification.json) - requires regeneration
- scripts/m3_a2g_data_verification.py - requires fixes before regeneration
- scripts/acceptance_checklist.py - requires replacement
- MANIFEST.md - requires completion

---

## CONSTRAINTS PRESERVED

✅ Commit bf33346 and history preserved (no git operations)
✅ No production code modifications
✅ No changes to original experiment artifacts (scripts/CSVs/TXT)
✅ No new scientific experiments
✅ Provisional O2 classification maintained
✅ Canonical scientific wording retained

---

## NEXT STEPS FOR COMPLETION

1. **Verification utility** - Fix schema expectations, implement complete CSV checks, regenerate JSON
2. **Acceptance checklist** - Replace superficial checks with substantive validation
3. **MANIFEST** - Complete inventory with all 49 tracked files
4. **Final reconciliation** - Complete remaining checkpoint report corrections
5. **Validation** - Run complete verification suite, document all limitations

**Estimated remaining work:** 3-4 hours for completion of all tasks with thorough testing

---

## SESSION LIMITATIONS

**Token usage:** ~106K/200K at summary creation
**Time constraints:** Continuation session from context overflow
**Scope:** Bounded to documentation/provenance corrections only

This summary documents progress made in the post-commit correction pass. The session successfully corrected the critical historical diagnostic misinterpretation and fixture set regression, but verification utility and acceptance checklist fixes remain incomplete.

**Status:** PARTIAL COMPLETION - Core documentation corrections done, verification tooling requires additional session
