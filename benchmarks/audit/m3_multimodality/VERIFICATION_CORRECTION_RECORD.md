# M3-A2g Verification Correction Record

**Date:** 2026-10-04
**Branch:** scientific-rebuild
**Starting HEAD:** 394ce0ccf8c46d3f37a2bf0e9a81e255886de5c0
**Status:** NOT COMPLETE until independent review. Nothing staged or committed.

Scope: the M3-A2g verification tooling, its regression tests, the derivative JSON, `MANIFEST.md`,
`README.md` and these correction reports. No production code, no experiment scripts or outputs,
no new experiments. The inputs being verified are SYNTHETIC experiment outputs.

---

## 1. Superseded claims (first pass, 2026-10-04)

The first pass of this correction produced `VERIFICATION_CORRECTION_RECORD.md` and
`M3_VERIFICATION_IMPLEMENTATION_REPORT.md` with the claims below. An independent review found them
wrong. They are kept here as correction history and are withdrawn.

| Claim made | What was actually true |
|------------|------------------------|
| "Implementation COMPLETE", "All 10 requirements met", "29/29 acceptance" | Required checks were missing. Acceptance passed 29/29 on inputs the reviewer corrupted in memory: exact-family `rep=99` replacing an expected replicate, `selected_K_cont` contradicting `delta_bic_cont`, and trajectory iteration `99999` replacing iteration 0 |
| "MANIFEST.md updated with complete inventory" | `MANIFEST.md` was not modified: 36 unique files covered, 21 omitted (including the three new files), 4 duplicate entries, stale hashes for the verifier script and the regenerated JSON |
| "Working tree clean" | Four modified and three untracked files |
| "Fixed trajectory counting bug (5→1634)" as work of this pass | HEAD already reported 1634 records; the counting was corrected in commit bf33346 |
| "Input hashes ... all match Git checkpoint" | No comparison against the checkpoint, the manifest or the stored JSON was implemented; acceptance only counted that 11 hashes were computed, after the verifier had overwritten the JSON |
| "Markdown links verified" | No link check was run |
| Regression tests cover missing/duplicate replicates, trajectory errors and stale checksums | Those tests were written, failed against the implementation, and were removed rather than fixed. Several remaining tests re-implemented validation logic instead of calling the validator |
| `acceptance_checklist.py` "Status: NEW" | It was a modified existing file |
| Regenerated JSON reported n_init convergence 0/20 | The CSV records `best_converged`, not `converged`; actual bookkeeping count is 15/20 |

Earlier, preserved and unedited: `VERIFICATION_TOOLING_COMPLETION_REPORT.md` (bf33346 era) states
"Status: COMPLETE" for verification tooling that covered 5 of 10 CSVs. That claim is likewise
superseded; the file is classified as a historical correction report in `MANIFEST.md`.

---

## 2. Changes in this pass

### `scripts/m3_a2g_data_verification.py` (modified)
- **Read-only by default.** Recomputes everything from the input bytes, compares the stored JSON
  (input hashes and every section) with the fresh result, and exits nonzero on an integrity
  failure or a stale JSON. `--write` regenerates the JSON only when integrity passes, after
  reporting the comparison.
- **Generator-derived designs** for all 10 CSVs: expected key sets (configurations × replicates,
  or case × setting) from the loops in the generator scripts; missing, duplicate and unexpected
  keys are reported.
- **Headers:** missing, unexpected and duplicate columns (raw header row, so duplicates are not
  hidden by `DictReader`). Column order is not enforced.
- **Fields:** every field parsed as its type (int, float, bool `True`/`False`, str); empty required
  values, unparseable values and wrong row lengths are errors; non-finite values are errors except
  the documented `nan` sentinels in `max_iter_sensitivity`, which must appear exactly when the
  generator writes them (`n_iters <= 100`; `n_iters <= 1`).
- **Domains:** e.g. `selected_K` ∈ {1, 2}; κ-type columns ≥ 0 (κ = 0 permitted); counts within
  budgets; weights and effective component counts within their ranges.
- **Row consistency** derived from the generators: n_init rate and spread arithmetic; stopping
  diagnostic `K_changed`, `loglik_gain`, `bic_change` and its `convergence_iter` approximation;
  degree/radian agreement.
- **Selection bookkeeping:** `selected_K` vs the implemented rule `2 if delta_bic < 0 else 1`
  (equality selects K = 1), for every pair, including both `_conv` and `_cont` pairs of the
  stopping diagnostic. Kept separate from convergence; non-converged selections are counted
  separately and are not treated as validated.
- **Trajectories:** strict parser. Declarations (`Converged`, `Iterations`, `Final loglik`), record
  header, 7 fields per record, iteration sequence 0..N−1, declared count = record count, final
  loglik = last record, nan deltas exactly at iteration 0, expected case order. Malformed lines are
  errors, never skipped. Likelihood decreases are reported (tolerance 1e-6, the 6-decimal
  resolution of the file; scope: 5 single-fit trajectories) as observations, not integrity failures.
- **Summaries restored** from the HEAD version (exact-family by F, low/zero-convergence
  configurations, true-K2 by separation, extended max_iter by configuration, n_init by case and the
  wrapped_1.0 diff) and corrected (n_init uses `best_converged`).

### `scripts/acceptance_checklist.py` (modified)
Fresh and non-writing. Section A (data integrity), B (figures cited in M3 documents vs fresh
counts), C (stored JSON vs fresh result), D (manifest coverage, duplicates, checksums, tracked files
present, protected artifacts byte-identical to the starting HEAD). Nonzero exit on any failure.
Scientific status is printed as information only.

### `scripts/compute_manifest_hashes.py` (modified)
Classification of every file in the directory; read-only `--check` (default) for inventory coverage,
duplicate entries, stale checksums and byte comparison of protected artifacts against the starting
HEAD via `git show`; `--update` regenerates the inventory section of `MANIFEST.md`.

### `scripts/test_m3_a2g_data_verification.py` (created)
Tests call the actual validators and commands on copies of the inputs in pytest temporary
directories. They include the three corruptions the reviewer demonstrated, missing and duplicate
replicates, duplicate headers, malformed required fields, nan-sentinel misuse, trajectory count,
record and declaration errors, stale checksums, non-writing default behaviour, refusal to write on
failure, and manifest duplicate, omitted, stale and unclassified detection.

### `results/m3_a2g_data_verification.json` (regenerated with `--write`)
SYNTHETIC derivative summary. Provenance states that the 1634-record count dates from bf33346.

### `MANIFEST.md`, `README.md` (modified); this record and `M3_VERIFICATION_IMPLEMENTATION_REPORT.md` (rewritten)

---

## 3. Results

Actual command outputs are in
[M3_VERIFICATION_IMPLEMENTATION_REPORT.md](M3_VERIFICATION_IMPLEMENTATION_REPORT.md). Counts quoted
there are produced by the commands and stored in `results/m3_a2g_data_verification.json`.

---

## 4. Remaining limitations and unresolved items

Verification scope:
- Random seeds are not recorded in the CSVs and are not verified.
- No correspondence between trajectory blocks and CSV rows is asserted (none is recorded).
- `scripts/m3_a2b_streamlined.py` and `scripts/m3_a2c_sanity_audit.py` depend on
  `/tmp/m3_a_diagnostic.py`, which is outside the repository. Its hash is recorded in
  `MANIFEST.md`; it is not preserved by the repository.
- Historical correction reports (`VERIFICATION_TOOLING_COMPLETION_REPORT.md`,
  `POST_COMMIT_CORRECTIONS_SUMMARY.md`, `DEEP_DOCUMENTATION_REVIEW_REPORT.md`) were not edited
  (outside the authorised file set); their superseded claims are recorded here and in the manifest.
- Manifest and preservation checks need `git` and the baseline revision in the local repository.
- `scripts/__pycache__/` contains bytecode written during the first pass (git-ignored). This pass
  ran Python with `-B` and pytest with `-p no:cacheprovider`.

Scientific (unchanged by this pass):
- n_init diagnostic INVALID (inner reseeding defect); stopping diagnostic INVALID.
- Extended max_iter experiment UNPAIRED; monotonicity evidence limited to 4 fixtures.
- No dependent-toroidal validation (mandatory before production).
- Classification remains **provisional O2**. KNOWN_ISSUES M3 remains open.
