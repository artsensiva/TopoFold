# M3-A2g Verification Correction: Implementation Report

**Date:** 2026-10-04
**Branch:** scientific-rebuild
**Starting HEAD:** 394ce0ccf8c46d3f37a2bf0e9a81e255886de5c0
**Status:** Second pass finished; NOT COMPLETE until independent review. Nothing staged or committed.

This report replaces the first-pass version, whose completion claims were wrong; they are listed in
[VERIFICATION_CORRECTION_RECORD.md](VERIFICATION_CORRECTION_RECORD.md) §1. The inputs are SYNTHETIC
experiment outputs. Passing these checks establishes format, bookkeeping and provenance integrity
only; the classification remains provisional O2.

## Working-tree inventory (relative to the starting HEAD)

Modified (tracked):
- `MANIFEST.md`
- `README.md`
- `results/m3_a2g_data_verification.json` (regenerated with `--write`)
- `scripts/acceptance_checklist.py`
- `scripts/compute_manifest_hashes.py`
- `scripts/m3_a2g_data_verification.py`

Untracked (new):
- `M3_VERIFICATION_IMPLEMENTATION_REPORT.md` (this file)
- `VERIFICATION_CORRECTION_RECORD.md`
- `scripts/test_m3_a2g_data_verification.py`

All paths are relative to `benchmarks/audit/m3_multimodality/`. No other file in the repository
was modified.

## Commands

Run from `benchmarks/audit/m3_multimodality/` with the repository virtual environment (numpy is
required; the system `python3` lacks it). `-B` and `-p no:cacheprovider` avoid writing caches.

```bash
python -B scripts/m3_a2g_data_verification.py
python -B scripts/acceptance_checklist.py
python -B scripts/compute_manifest_hashes.py
python -B -m pytest -p no:cacheprovider scripts/test_m3_a2g_data_verification.py
```

## Outputs (final run, after the last manifest update)

```text
$ python -B scripts/m3_a2g_data_verification.py
M3-A2g DATA VERIFICATION (SYNTHETIC inputs; read-only unless --write)
Datasets (header / fields / row checks / selection / coverage):
  OK   exact_family                   rows= 360 fields_err=0 row_err=0 sel_err=0 keys=360/360 missing=0 unexpected=0 dup=0
  OK   k2_convergence                 rows= 320 fields_err=0 row_err=0 sel_err=0 keys=320/320 missing=0 unexpected=0 dup=0
  OK   stress_regime                  rows= 280 fields_err=0 row_err=0 sel_err=0 keys=280/280 missing=0 unexpected=0 dup=0
  OK   extended_maxiter               rows= 180 fields_err=0 row_err=0 sel_err=0 keys=180/180 missing=0 unexpected=0 dup=0
  OK   ninit_sensitivity              rows=  20 fields_err=0 row_err=0 sel_err=0 keys=20/20 missing=0 unexpected=0 dup=0
  OK   stopping_criterion_validation  rows=  20 fields_err=0 row_err=0 sel_err=0 keys=20/20 missing=0 unexpected=0 dup=0
  OK   tolerance_comparison           rows=  16 fields_err=0 row_err=0 sel_err=0 keys=16/16 missing=0 unexpected=0 dup=0
  OK   max_iter_sensitivity           rows=  16 fields_err=0 row_err=0 sel_err=0 keys=16/16 missing=0 unexpected=0 dup=0
  OK   initialization_comparison      rows=  15 fields_err=0 row_err=0 sel_err=0 keys=15/15 missing=0 unexpected=0 dup=0
  OK   duplicate_components           rows=   8 fields_err=0 row_err=0 sel_err=0 keys=8/8 missing=0 unexpected=0 dup=0
  total CSV rows: 1235
Trajectories: OK blocks=5 records=1634 loglik decreases (tol 1e-06): 0
Stored derivative m3_a2g_data_verification.json: matches fresh result
PASS: integrity checks passed and stored derivative is current (format/bookkeeping only; not scientific validation)
[exit 0]

$ python -B scripts/acceptance_checklist.py
M3-A2g ACCEPTANCE CHECKLIST (fresh, non-writing; SYNTHETIC inputs)

A. Data integrity (fresh)
  PASS  exact_family: header
  PASS  exact_family: field values
  PASS  exact_family: row consistency
  PASS  exact_family: coverage 360/360 keys
  PASS  exact_family: selection bookkeeping selected_K
  PASS  k2_convergence: header
  PASS  k2_convergence: field values
  PASS  k2_convergence: row consistency
  PASS  k2_convergence: coverage 320/320 keys
  PASS  k2_convergence: selection bookkeeping selected_K
  PASS  stress_regime: header
  PASS  stress_regime: field values
  PASS  stress_regime: row consistency
  PASS  stress_regime: coverage 280/280 keys
  PASS  stress_regime: selection bookkeeping selected_K
  PASS  extended_maxiter: header
  PASS  extended_maxiter: field values
  PASS  extended_maxiter: row consistency
  PASS  extended_maxiter: coverage 180/180 keys
  PASS  extended_maxiter: selection bookkeeping selected_K
  PASS  ninit_sensitivity: header
  PASS  ninit_sensitivity: field values
  PASS  ninit_sensitivity: row consistency
  PASS  ninit_sensitivity: coverage 20/20 keys
  PASS  stopping_criterion_validation: header
  PASS  stopping_criterion_validation: field values
  PASS  stopping_criterion_validation: row consistency
  PASS  stopping_criterion_validation: coverage 20/20 keys
  PASS  stopping_criterion_validation: selection bookkeeping selected_K_conv+selected_K_cont
  PASS  tolerance_comparison: header
  PASS  tolerance_comparison: field values
  PASS  tolerance_comparison: row consistency
  PASS  tolerance_comparison: coverage 16/16 keys
  PASS  max_iter_sensitivity: header
  PASS  max_iter_sensitivity: field values
  PASS  max_iter_sensitivity: row consistency
  PASS  max_iter_sensitivity: coverage 16/16 keys
  PASS  initialization_comparison: header
  PASS  initialization_comparison: field values
  PASS  initialization_comparison: row consistency
  PASS  initialization_comparison: coverage 15/15 keys
  PASS  duplicate_components: header
  PASS  duplicate_components: field values
  PASS  duplicate_components: row consistency
  PASS  duplicate_components: coverage 8/8 keys
  PASS  trajectories: declarations, counts, sequences, records

B. Documented figures vs fresh data
  PASS  total CSV rows = 1235  [fresh 1235]
  PASS  trajectory blocks/records = 5/1634  [fresh 5/1634]
  PASS  exact_family = 111/360  [fresh 111/360]
  PASS  k2_convergence = 258/320  [fresh 258/320]
  PASS  stress_regime = 102/280  [fresh 102/280]
  PASS  stress sigma=1.0,F=1000 = 3/10  [fresh 3/10]
  PASS  extended_maxiter max_iter=2000 = 56/60  [fresh 56/60]
  PASS  ninit wrapped_1.0 loglik diff = 2.63  [fresh 2.6307315391318298]

C. Derivative freshness (read-only comparison)
  PASS  stored m3_a2g_data_verification.json exists
  PASS  stored input hashes match current input bytes
  PASS  stored derivative equals fresh result

D. Manifest and preservation
  PASS  manifest: no duplicate entries
  PASS  manifest: no omitted from manifest
  PASS  manifest: no listed but missing
  PASS  manifest: no checksum mismatches
  PASS  manifest: no unclassified files
  PASS  manifest: no classified but missing
  PASS  all tracked M3 files present
  PASS  protected artifacts identical to 394ce0ccf8c4 (50/50)

Scientific status (information only; not changed by this checklist)
  Provisional O2 - functional but fragile. n_init and stopping diagnostics INVALID;
  extended max_iter UNPAIRED; monotonicity evidence limited; dependent-toroidal
  validation not performed. Not production-ready. KNOWN_ISSUES M3 remains open.

Summary: 65/65 checks passed
RESULT: PASS (format, bookkeeping and provenance integrity only)
[exit 0]

$ python -B scripts/compute_manifest_hashes.py
Manifest: 57 entries, 57 files on disk (excluding MANIFEST.md and caches)
  duplicate_entries: 0
  omitted_from_manifest: 0
  listed_but_missing: 0
  checksum_mismatches: 0
  unclassified_files: 0
  classified_but_missing: 0
  tracked files missing on disk: 0
Protected vs 394ce0ccf8c4: 50/50 identical; changed=[] missing=[]
PASS
[exit 0]

$ python -B -m pytest -p no:cacheprovider -q scripts/test_m3_a2g_data_verification.py
37 passed
[exit 0]
```

The pytest timing is removed from the last line because it varies between runs.

## Not covered / unresolved

See [VERIFICATION_CORRECTION_RECORD.md](VERIFICATION_CORRECTION_RECORD.md) §4. In short: seeds are
not verified; no trajectory-to-CSV correspondence is asserted; the pre-M3-A2d scripts depend on
`/tmp/m3_a_diagnostic.py` outside the repository; historical correction reports with superseded
claims were not edited; Markdown anchors were not checked. Scientific limitations are unchanged:
the n_init and stopping diagnostics are invalid, extended max_iter is unpaired, monotonicity
evidence is limited, and dependent-toroidal validation has not been done.
