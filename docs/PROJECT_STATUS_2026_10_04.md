# TopoFold: Comprehensive Project Review and Restoration Plan

**Date:** 2026-10-04
**Author:** Claude Code (Sonnet 4.5)
**Purpose:** Full project review, context restoration, and next-steps roadmap

---

## EXECUTIVE SUMMARY

**TopoFold** is a research prototype computing SE(3)-invariant differential-geometry descriptors of biomolecular backbones (curvature κ, torsion τ, local writhe, Cβ orientation). The project is under **scientific re-validation** following discovery that v0.8.3 contained:

- Synthetic data presented as "real MD"
- Estimator artefacts reported as free energies
- Hard-coded statistics and physically impossible benchmarks
- A writhe implementation bug (−0.5× correct value)

**Current Status:**
- ✅ Core Rust tests: 18/18 passing
- ✅ Writhe bug: FIXED (M1, committed on scientific-rebuild)
- ✅ M3-A2g multimodality audit: COMPLETE with provisional O2 classification
- ⚠️ 10 critical/major issues remaining in KNOWN_ISSUES.md
- ⚠️ 11 uncommitted files in M3 audit directory

**This Session's Work:**
- Completed M3-A2g verification tooling (schema validation, acceptance checklist, MANIFEST)
- All 29 substantive acceptance criteria now PASS
- Documentation package ready for external review

---

## I. PROJECT STRUCTURE REVIEW

### A. Repository Organization

```
TopoFold/
├── crates/                          # Rust workspace
│   ├── topofold-core/              # Core geometry (κ, τ, writhe, bimodality)
│   ├── topofold-index/             # VP-tree conformational index
│   ├── topofold-io/                # PDB/DCD I/O
│   └── topofold-python/            # Python bindings (PyO3)
├── tests/                           # Python integration tests
├── benchmarks/
│   └── audit/                       # Scientific validation scripts
│       ├── writhe_check.py         # M1 reproduction (v0.8.3 bug)
│       ├── barrier_check.py        # C3 reproduction (PMF artefacts)
│       ├── bc_check.py             # M2 reproduction (circular stats)
│       ├── fair_baseline.py        # C2 reproduction (circular eval)
│       └── m3_multimodality/       # M3 multimodality detection audit
├── docs/
│   ├── KNOWN_ISSUES.md             # Development backlog (11 issues)
│   ├── MANUSCRIPT_DRAFT.md         # Research manuscript (needs update)
│   └── whitepaper.md               # Technical whitepaper (needs update)
├── AGENTS.md                        # Scientific integrity rules (hard)
├── CLAUDE.md                        # Claude-specific instructions
└── README.md                        # Corrected with withdrawal notice
```

### B. Scientific Integrity Framework

**Three-Layer Protection:**

1. **AGENTS.md** - 10 hard rules (data provenance, no silent fallbacks, physics first)
2. **Build/test protocol** - Rust tests → maturin rebuild → pytest → audit checks
3. **KNOWN_ISSUES.md** - Transparent backlog, issue closed only with regression test

**Key Principles:**
- Every synthetic dataset must say SYNTHETIC
- No hand-typed scientific results (code → JSON → figures/text)
- Equal information for baselines
- Independent ground truth labels
- Circular statistics for angles

### C. Current Branch Structure

```
main (5d609e5)
  └─ scientific-rebuild (bf33346) ← HEAD
       └─ docs/audit-corrections (merged, withdrawn branch)
```

**scientific-rebuild branch:**
- Committed: M1 writhe fix + M3-A2g corrections (bf33346)
- Uncommitted: M3 verification tooling completion (this session)

---

## II. DEVELOPMENT HISTORY ANALYSIS

### A. Git Timeline (Last 30 Commits)

**Phase 1: v0.8.3 Publication Push (May-Sept 2026)**
- `v0.8.3` - RNA ribonucleic geometry, Streamlit viewer
- `v0.8.2` - PROTAC benchmark (later withdrawn as invalid)
- `v0.8.1` - IDP spectral topological density
- `v0.8.0` - Metamorphic proteins, allosteric networks
- v0.7.x - Rotamer detection (later withdrawn), manuscript draft

**Phase 2: Audit Discovery (Sept-Oct 2026)**
- `bb87544` (v0.8.3-pre-audit tag) - Last commit before audit
- `5d609e5` - Withdrawal notice added to main
- Created `docs/audit-corrections` branch (withdrawn benchmarks)

**Phase 3: Scientific Rebuild (Oct 2026)**
- `418ca13` - CI checks added (Rust, Python, audit reproduction)
- `95244fd` - **M1 FIXED**: Writhe normalization/sign corrected
- `37ce135` - M3 multimodality artifacts preserved
- `bf33346` - **M3-A2g corrections** (provisional O2) ← HEAD

### B. Audit Milestones

**M3 Multimodality Detection Audit:**
- **M3-A2e** - Exact/wrapped-normal null validation (checkpoint 6-10)
- **M3-A2f** - Three estimators comparison (U/C/P regularization)
- **M3-A2g** - Stress regime, convergence, diagnostics (CURRENT)

**Classification Ladder:**
- O1 = Failure (fundamental error)
- O2 = Functional but fragile (limited validation) ← **M3-A2g STATUS**
- O3 = Robust (comprehensive validation, production-ready)

### C. v0.8.3 Problems Discovered

**Critical (5 issues):**
- C1: Synthetic data presented as real MD
- C2: Circular evaluation, unequal baselines
- C3: "Free energies" are estimator artefacts
- C4: θβ cannot detect rotamers (Cβ doesn't move with χ1)
- C5: PROTAC benchmark invalid (wrong PDB, hard-coded)

**Major (8 issues):**
- M1: ✅ **FIXED** - Writhe −0.5× bug
- M2: No circular statistics in bimodality scan
- M3: BC is not a test (no null model, no FDR)
- M4: RNA base vector cannot measure χ (glycosidic)
- M5: Fréchet not SE(3)-invariant, wrong units
- M6: Index cascade is approximate (not proven lower bound)
- M7: Mutual information inconsistent (3 descriptions)
- M8: Hard-coded results, unreproducible numbers

**Minor (11 documentation issues)**

---

## III. CONVERSATION HISTORY ANALYSIS

### Session Context

**Previous Session:** M3-A2g documentation corrections pass
- Started: User provided continuation context from previous session
- Completed: 28 files committed (bf33346), 25/25 acceptance criteria passed

**Current Session (This Conversation):**

**Phase 1: Completion of Deferred Work (Messages 1-11)**
- User: "Continue from where we left off"
- Task: Complete verification tooling deferred in deep review
- Accomplished:
  1. Fixed verification script schema validation
  2. Regenerated m3_a2g_data_verification.json (all schemas valid)
  3. Replaced acceptance_checklist.py with substantive checks
  4. Completed MANIFEST.md with updated hashes

**Phase 2: Final Validation (Messages 12-Present)**
- All 29 substantive checks PASS
- Created VERIFICATION_TOOLING_COMPLETION_REPORT.md
- User: "Comprehensive review and propose next steps"

### Key Insights from Our Work

1. **Verification tooling is critical** - Superficial pattern matching was insufficient
2. **Structured validation from data** - Checks must validate against actual CSVs/JSONs
3. **Complete provenance tracking** - MANIFEST with SHA-256 hashes for all artifacts
4. **Bounded correction passes** - Clear scope, constraints, completion criteria

---

## IV. CURRENT PROJECT STATE

### A. What Works ✅

**Core Rust Implementation:**
- 18/18 tests passing
- Writhe corrected with 8 regression tests
- SE(3) invariance validated
- No unsafe code

**M3-A2g Audit Package:**
- 8 experimental scripts (1235 CSV rows)
- 10 CSVs + 1 TXT trajectory (preserved unchanged)
- 6 corrected documentation files
- Verification tooling with 29 substantive checks
- Complete MANIFEST (49 files tracked)
- **Status:** Provisional O2, ready for external review

**Documentation:**
- Honest README with correction notice
- KNOWN_ISSUES.md backlog (transparent)
- Scientific integrity framework (AGENTS.md)

### B. What's Broken/Incomplete ⚠️

**Uncommitted Work:**
- 9 modified M3 files in working directory
- 5 new M3 files (reports, utilities)
- 1 patch file (topofold-v0.8.3-corrections.patch)

**KNOWN_ISSUES.md Backlog:**
- 4 critical issues (C1, C2, C3, C4, C5 withdrawn)
- 7 major issues (M2-M8, M1 resolved)
- 11 minor documentation issues

**No Real Validation:**
- All benchmarks are synthetic/interpolated
- No independently generated MD trajectories analyzed
- Baselines not fairly compared

### C. Git Working Directory Status

```
Modified (9 files):
  M M3_A2G_CHECKPOINT_CLOSURE.md
  M M3_A2G_COMPLETE_ANALYSIS.md
  M M3_A2G_CORRECTIONS_PROVENANCE.md
  M M3_STATUS.md
  M MANIFEST.md
  M README.md
  M reports/M3_A2g_CHECKPOINT_REPORT.md
  M results/m3_a2g_data_verification.json
  M scripts/m3_a2g_data_verification.py

Untracked (5 files):
  ?? DEEP_DOCUMENTATION_REVIEW_REPORT.md
  ?? POST_COMMIT_CORRECTIONS_SUMMARY.md
  ?? VERIFICATION_TOOLING_COMPLETION_REPORT.md
  ?? scripts/acceptance_checklist.py
  ?? scripts/compute_manifest_hashes.py
  ?? topofold-v0.8.3-corrections.patch (root dir)
```

**Total Changes:** 218 insertions, 45 deletions across 9 files

---

## V. CONTEXT RESTORATION PLAN

### A. Immediate Git Cleanup

**Step 1: Review and Commit Current M3 Work**

Rationale: The verification tooling completion is a discrete, bounded unit of work that passes all validation checks. Committing now establishes a clean checkpoint.

```bash
# 1. Review uncommitted changes
git diff benchmarks/audit/m3_multimodality/

# 2. Add verification tooling updates
git add benchmarks/audit/m3_multimodality/scripts/m3_a2g_data_verification.py
git add benchmarks/audit/m3_multimodality/scripts/acceptance_checklist.py
git add benchmarks/audit/m3_multimodality/scripts/compute_manifest_hashes.py
git add benchmarks/audit/m3_multimodality/results/m3_a2g_data_verification.json

# 3. Add documentation corrections
git add benchmarks/audit/m3_multimodality/M3_STATUS.md
git add benchmarks/audit/m3_multimodality/README.md
git add benchmarks/audit/m3_multimodality/M3_A2G_CORRECTIONS_PROVENANCE.md
git add benchmarks/audit/m3_multimodality/M3_A2G_CHECKPOINT_CLOSURE.md
git add benchmarks/audit/m3_multimodality/M3_A2G_COMPLETE_ANALYSIS.md
git add benchmarks/audit/m3_multimodality/MANIFEST.md
git add benchmarks/audit/m3_multimodality/reports/M3_A2g_CHECKPOINT_REPORT.md

# 4. Add session reports
git add benchmarks/audit/m3_multimodality/DEEP_DOCUMENTATION_REVIEW_REPORT.md
git add benchmarks/audit/m3_multimodality/POST_COMMIT_CORRECTIONS_SUMMARY.md
git add benchmarks/audit/m3_multimodality/VERIFICATION_TOOLING_COMPLETION_REPORT.md

# 5. Commit with descriptive message
git commit -m "feat(m3-audit): complete M3-A2g verification tooling

- Fix verification script schema validation (all 5 CSVs now valid)
- Replace superficial acceptance checks with 29 substantive criteria
- Complete MANIFEST with updated SHA-256 hashes (49 files tracked)
- Regenerate m3_a2g_data_verification.json (all schemas valid)
- Add acceptance_checklist.py and compute_manifest_hashes.py utilities

All substantive acceptance criteria PASS (29/29).
M3-A2g documentation package ready for external review.

Closes post-commit correction pass from bf33346."
```

**Step 2: Handle Remaining Files**

```bash
# Review patch file
cat topofold-v0.8.3-corrections.patch

# If it's historical/reference, move to benchmarks/audit/
git mv topofold-v0.8.3-corrections.patch benchmarks/audit/

# Or if superseded by commits, document and remove
git add -A  # Stage any moves
git commit -m "chore: archive v0.8.3 corrections patch"
# OR
rm topofold-v0.8.3-corrections.patch
```

### B. Project Documentation Updates

**Step 3: Add This Review Document to Project**

```bash
# Move to docs/ for permanent reference
mv COMPREHENSIVE_PROJECT_REVIEW_AND_RESTORATION_PLAN.md docs/PROJECT_STATUS_2026_10_04.md

git add docs/PROJECT_STATUS_2026_10_04.md
git commit -m "docs: add comprehensive project review and restoration plan (Oct 2026)"
```

### C. Session State Documentation

**Create Session Handoff Document:**

Create `docs/SESSION_HANDOFF.md` capturing:
- What was accomplished (verification tooling)
- Current state (provisional O2, ready for external review)
- Next steps (see Section VI)
- Known blockers (none for next phase)

This allows any future AI assistant or human to pick up where we left off.

---

## VI. PROPOSED NEXT WORK

### A. Priority Matrix

```
┌─────────────────────────────────────────────────────────────┐
│ IMMEDIATE (This Week)                                       │
├─────────────────────────────────────────────────────────────┤
│ 1. Commit current M3-A2g work (verification tooling)        │
│ 2. External scientific review of M3 package                 │
│ 3. Resolve git working directory state                      │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│ SHORT TERM (1-2 Weeks) - Continue Scientific Rebuild        │
├─────────────────────────────────────────────────────────────┤
│ 4. M2: Implement circular statistics in bimodality scan     │
│ 5. C1 Phase 1: Label all synthetic data as SYNTHETIC        │
│ 6. Add first real MD validation (BPTI or villin)           │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│ MEDIUM TERM (1 Month) - Core Fixes                          │
├─────────────────────────────────────────────────────────────┤
│ 7. M5: Fix Fréchet distance (SE(3) invariance, units)       │
│ 8. M4: Fix RNA χ descriptor (glycosidic angle)              │
│ 9. C4: Remove θβ rotamer claims (or add χ1/χ2)             │
│ 10. M8: Results pipeline (code → JSON → figures/text)       │
└─────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────┐
│ LONG TERM (2-3 Months) - Fair Validation                    │
├─────────────────────────────────────────────────────────────┤
│ 11. C2: Implement fair baseline comparisons                 │
│ 12. C1 Complete: Real MD benchmark suite                    │
│ 13. M3: Add null model + FDR to bimodality test            │
│ 14. Python test suite restoration (pytest currently fails)  │
└─────────────────────────────────────────────────────────────┘
```

### B. Detailed Next Steps

#### IMMEDIATE: Issue #1 - Commit Current Work

**Task:** Complete git cleanup from Section V.A
**Estimated Time:** 30 minutes
**Prerequisites:** Review uncommitted changes
**Deliverables:**
- Clean git working directory
- M3-A2g verification tooling committed
- Session reports preserved

**Commands:**
```bash
# Execute Section V.A steps 1-5
# Verify: git status should show clean working directory
```

---

#### IMMEDIATE: Issue #2 - External Review

**Task:** Send M3-A2g package for external scientific review
**Estimated Time:** Variable (external dependency)
**Prerequisites:** Issue #1 complete

**What to Send:**
- All files in `benchmarks/audit/m3_multimodality/`
- Specifically:
  - M3_STATUS.md (classification summary)
  - M3_A2G_COMPLETE_ANALYSIS.md (analysis)
  - M3_A2G_CORRECTIONS_PROVENANCE.md (all corrections)
  - reports/M3_A2g_CHECKPOINT_REPORT.md (detailed report)
  - All 10 CSVs + trajectories.txt (original data)
  - verification JSON + acceptance checklist

**Review Questions:**
1. Is provisional O2 classification justified?
2. Are diagnostic issues correctly characterized?
3. Are withdrawn claims comprehensive?
4. Any missing issues/concerns?

---

#### SHORT TERM: Issue #4 - M2 Circular Statistics

**KNOWN_ISSUES.md Reference:** M2
**Estimated Time:** 2-3 days
**Complexity:** Medium

**Problem:**
```rust
// Current (WRONG): bimodality.rs
let sum: f64 = window_values.iter().sum();
let mean = sum / window_values.len() as f64;
```

Angles near ±π produce false positives (BC 0.96 instead of 0.35).

**Solution Plan:**

1. **Add circular mean/variance utilities** (`crates/topofold-core/src/circular.rs`):
```rust
pub fn circular_mean(angles: &[f64]) -> f64 {
    let sum_sin: f64 = angles.iter().map(|&a| a.sin()).sum();
    let sum_cos: f64 = angles.iter().map(|&a| a.cos()).sum();
    sum_sin.atan2(sum_cos)
}

pub fn circular_variance(angles: &[f64]) -> f64 {
    // Mardia-Jupp circular variance
    let n = angles.len() as f64;
    let r = (/* resultant length */);
    1.0 - r
}
```

2. **Update bimodality computation** for τ, θβ (both circular)

3. **Regression test:**
```python
# benchmarks/audit/m2_circular_fix_validation.py
# Test: uniform around 180° → BC ≈ 0.33 (not 0.98)
# Test: β-strand near −170° → no false positive
```

4. **Update KNOWN_ISSUES.md:** M2 status → RESOLVED

---

#### SHORT TERM: Issue #5 - C1 Phase 1 Labeling

**KNOWN_ISSUES.md Reference:** C1
**Estimated Time:** 1 week
**Complexity:** High (many files)

**Problem:** Synthetic data presented as "real MD"

**Phase 1 Solution:** Add SYNTHETIC labels WITHOUT rewriting benchmarks

**Files to Update:**

1. **Benchmark scripts** - Add docstring/header:
```python
"""
BPTI Synthetic Test Ensemble

**SYNTHETIC BENCHMARK**
Data source: Interpolated between 5PTI (closed) and 3OTJ (open)
Generation: Linear interpolation + Gaussian noise (σ=0.5 Å)
States: Built-in by construction (t < 0.5 = closed, t ≥ 0.5 = open)
Purpose: Controlled validation of bimodality detection
NOT representative of: Real MD dynamics, free energies, kinetics
"""
```

2. **All figures** - Update titles:
```python
fig.suptitle("BPTI Synthetic Test: BC Profile (SYNTHETIC)")
plt.savefig("assets/SYNTHETIC_bpti_bc_profile.png")
```

3. **README/manuscript** - Already corrected with withdrawal notice

4. **Regression test:**
```bash
# Check all .py files in benchmarks/ for misleading claims
grep -r "authentic\|real MD\|experimental trajectory" benchmarks/
# Should return 0 matches (except in comments explaining what's NOT true)
```

---

#### SHORT TERM: Issue #6 - First Real MD Validation

**Prerequisites:** C1 Phase 1 complete
**Estimated Time:** 1-2 weeks
**Complexity:** High (new territory)

**Objective:** Analyze ONE independently generated MD trajectory

**Candidates:**
1. **BPTI** (Belyaeva et al., Zenodo 7347434) - Already cited
2. **Villin headpiece** (Shaw lab, D. E. Shaw Research)
3. **Trp-cage** (compact, well-studied)

**Requirements:**
- Independent source (not generated by our scripts)
- Provenance documented (who ran it, which force field, when)
- At least 100 ns equilibrated trajectory
- Known conformational transition (optional but helpful)

**Validation Plan:**

```python
# benchmarks/real_md/bpti_belyaeva_2023.py
"""
BPTI Trypsin-Inhibitor Complex Real MD Analysis

Data source: Belyaeva, Zlobin, Maslova, Golovin, PCCP 2023
DOI: 10.1039/d2cp05502c
Zenodo: 7347434
Force fields: AMBER14SB, CHARMM36m, OPLS-AA/L, GROMOS53a6, AMBER99SB-ILDN, AMBER03
System: Trypsin–BPTI complex, explicit solvent
Frames: [to be determined from actual trajectory]
Purpose: Independent validation of TopoFold bimodality detection
"""

# 1. Load trajectory (not generate)
# 2. Compute TopoFold descriptors
# 3. Compare with baseline (local RMSD, Cα-Cα distances)
# 4. Report results WITHOUT claiming superiority
# 5. Document limitations
```

**Success Criteria:**
- Script downloads real data (or errors if missing)
- No synthetic fallback
- Results written to JSON
- Limitations documented

---

#### MEDIUM TERM: Issue #7 - M5 Fréchet Fix

**KNOWN_ISSUES.md Reference:** M5
**Estimated Time:** 3-4 days
**Complexity:** Medium-High

**Problem:**
- `discrete_frechet_distance` on raw Cartesian coords (not SE(3)-invariant)
- Claimed as "exact SE(3)-invariant" in XCL1 benchmark
- Wrong units (claimed Å, actually in (κ,τ) radians)

**Solution Plan:**

1. **Separate APIs:**
```rust
// topofold-index/src/frechet.rs

/// Cartesian Fréchet distance (NOT SE(3)-invariant)
pub fn frechet_cartesian(
    traj1: &[Point3<f64>],
    traj2: &[Point3<f64>]
) -> f64 {
    // Current implementation
    // Returns distance in Å
}

/// Invariant Fréchet distance in (κ, τ) space
pub fn frechet_invariant(
    kappa_tau1: &[(f64, f64)],
    kappa_tau2: &[(f64, f64)]
) -> f64 {
    // Returns distance in radians
    // Already implemented, rename from discrete_frechet_distance
}

/// Per-residue RMSD (for same-sequence frames)
pub fn per_residue_distances(
    frame1: &[Point3<f64>],
    frame2: &[Point3<f64>]
) -> Vec<f64> {
    // Aligned per-residue distances in Å
}
```

2. **Fix XCL1 benchmark:**
```python
# Use per_residue_distances for same-sequence comparison
# OR align then use Cartesian Fréchet
# State units explicitly
```

3. **Regression test:**
```rust
#[test]
fn test_frechet_se3_invariance() {
    // Rotate/translate traj → invariant_frechet unchanged
    // Cartesian changes
}
```

---

### C. Long-Term Roadmap

#### **Q4 2026: Scientific Rebuild Complete**

**Goals:**
- ✅ All M-series issues resolved (M1-M8)
- ✅ C1 Phase 1 complete (SYNTHETIC labels)
- ✅ At least 1 real MD validation
- ✅ Python tests restored (pytest passing)

**Deliverables:**
- v0.9.0-rc1 (release candidate)
- Updated manuscript draft
- Fair baseline comparison on synthetic benchmark

#### **Q1 2027: Real Validation Suite**

**Goals:**
- C1 Complete (multiple real MD trajectories)
- C2 Complete (fair baselines: local RMSD, Cα-Cα, φ/ψ, TICA on same residues)
- M3 Complete (null model, FDR, significance testing)

**Deliverables:**
- v0.9.0 stable release
- Preprint submission with corrected claims
- Benchmark suite documented and reproducible

#### **Q2 2027: Production Readiness**

**Goals:**
- CI/CD pipeline (currently missing)
- Performance benchmarks (speed/memory)
- API stability guarantees
- Documentation complete

**Deliverables:**
- v1.0.0 (production-ready)
- Published paper
- Public benchmark repository

---

## VII. RISKS AND MITIGATION

### A. Technical Risks

**Risk 1: Circular Statistics Edge Cases**
- **Impact:** False positives in bimodality detection
- **Probability:** Medium
- **Mitigation:** Comprehensive test suite with angles near ±π, multiple wrappings

**Risk 2: Real MD Data Access**
- **Impact:** Cannot validate on independent trajectories
- **Probability:** Low (many public datasets)
- **Mitigation:** Use well-documented public archives (Zenodo, OSF, figshare)

**Risk 3: Fair Baselines Show TopoFold Equivalent**
- **Impact:** Method novelty questioned
- **Probability:** Medium-High (fair_baseline.py already showed local distances ≈ TopoFold)
- **Mitigation:**
  - Honest reporting (already committed to in AGENTS.md)
  - Focus on computational efficiency, alignment-free property
  - Identify specific use cases where TopoFold excels

### B. Scientific Risks

**Risk 4: External Review Finds More Issues**
- **Impact:** Additional corrections required, timeline delay
- **Probability:** Medium
- **Mitigation:**
  - Transparent KNOWN_ISSUES.md
  - Conservative claims
  - Quick incorporation of feedback

**Risk 5: Bimodality Without Statistical Testing**
- **Impact:** Cannot make definitive claims
- **Probability:** High (M3 already documented)
- **Mitigation:**
  - Call it "screening heuristic" (already done)
  - Add null model + FDR (M3 fix)
  - Compare with dip test, Gaussian mixture ΔBIC

### C. Project Management Risks

**Risk 6: Context Loss Between Sessions**
- **Impact:** Rework, inconsistency
- **Probability:** Medium
- **Mitigation:**
  - This restoration plan document
  - docs/SESSION_HANDOFF.md
  - Detailed git commit messages

**Risk 7: Scope Creep During Fixes**
- **Impact:** Never-ending corrections
- **Probability:** High
- **Mitigation:**
  - Work from KNOWN_ISSUES.md backlog only
  - One issue per branch
  - Mark resolved with regression test

---

## VIII. SUCCESS METRICS

### A. Technical Metrics

**Code Quality:**
- ✅ Rust tests: 18/18 passing → Maintain 100%
- ⚠️ Python tests: Currently broken → Restore to 100%
- ⚠️ Audit checks: 4/4 passing → Add more as issues fixed
- ⚠️ CI/CD: None → Implement (GitHub Actions)

**Documentation:**
- ✅ KNOWN_ISSUES.md: 11 issues documented → Resolve to 0
- ✅ README correction notice: Present → Maintain
- ⚠️ API docs: Incomplete → Coverage >80%

### B. Scientific Metrics

**Validation:**
- ✅ M3-A2g audit: Provisional O2 → Maintain or upgrade to O3
- ⚠️ Real MD validations: 0 → At least 3
- ⚠️ Fair baselines: 0 → At least 2 (local RMSD, Cα-Cα distances)

**Reproducibility:**
- ✅ Audit scripts reproduce bugs: 4/4 → Maintain
- ⚠️ Benchmark results in JSON: 0% → 100%
- ⚠️ Figures generated from JSON: 0% → 100%

### C. Project Health

**Velocity:**
- Current: ~1 issue/week (M1 took 1 week)
- Target: ~2 issues/week (as we build momentum)

**Quality Gates:**
- Each issue must have regression test
- KNOWN_ISSUES.md status updated
- Documentation updated if behavior changed

---

## IX. RECOMMENDATIONS

### For User (Artem)

1. **Immediate (Today):**
   - Review this restoration plan
   - Approve git commit strategy (Section V.A)
   - Decide on external review process for M3-A2g

2. **This Week:**
   - Execute git cleanup (or delegate to Claude)
   - Begin M2 circular statistics fix (next priority issue)
   - Consider pytest restoration (currently broken)

3. **This Month:**
   - Complete M2, C1 Phase 1, first real MD validation
   - Re-evaluate project timeline based on progress
   - Consider CI/CD setup (GitHub Actions)

### For Future AI Assistants

1. **Session Start Protocol:**
   - Read docs/KNOWN_ISSUES.md
   - Check docs/SESSION_HANDOFF.md (if exists)
   - Review this restoration plan (docs/PROJECT_STATUS_2026_10_04.md)
   - Ask user which issue to work on

2. **Work Protocol:**
   - One issue from KNOWN_ISSUES.md at a time
   - Create branch per issue (issue-M2-circular-stats)
   - Regression test mandatory
   - Update KNOWN_ISSUES.md status
   - Commit with descriptive message

3. **Scientific Integrity:**
   - AGENTS.md rules are NON-NEGOTIABLE
   - Never synthesize data when real files missing
   - Never hand-type scientific results
   - Always use circular statistics for angles
   - Document limitations honestly

---

## X. CONCLUSION

TopoFold is at a critical juncture:

**Where We Are:**
- Core geometry validated (writhe fixed, tests passing)
- M3-A2g audit complete (provisional O2)
- Scientific integrity framework established
- 11 known issues remaining

**Where We're Going:**
- Complete scientific rebuild (M2-M8)
- Real MD validation suite
- Fair baseline comparisons
- Production-ready v1.0.0

**How We Get There:**
- Work systematically from KNOWN_ISSUES.md backlog
- One issue → branch → test → commit → close
- Honest reporting (no overstated claims)
- Transparent progress tracking

**The project has already demonstrated the capacity for rigorous self-correction.** The hard work of discovering and documenting the issues is complete. Now we execute the fixes methodically.

**Next Immediate Action:** Execute git cleanup (Section V.A) to establish clean checkpoint, then begin M2 circular statistics fix.

---

**Prepared by:** Claude Code (Sonnet 4.5)
**Date:** 2026-10-04
**For:** TopoFold Scientific Rebuild Initiative
**Status:** Ready for review and execution
