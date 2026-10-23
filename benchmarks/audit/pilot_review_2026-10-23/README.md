# Separate-Agent Usability Review of the Descriptive Tool (2026-10-23) — archived evidence

**What this is:** a usability review of `python -m topofold_inputs.ensemble` (v0.2.0, repository commit `389f2d7`),
carried out by a separate agent in `/tmp/topofold-reviewer-pilot-7oe2feoy/`.

**What it is not:**
- not independent human validation;
- not a fresh-install test (it used the existing interpreters and local data);
- not biological validation.

The reviewer had prior familiarity with the project from code reviews.

## Recorded outcome (from `REVIEW.md` / `REVIEW.json`)

- **Runs:** three CLI runs completed without implementation changes: the documented example (`1ha8A00`, 3 vs 3), a
  reviewer-chosen domain (`2yztA00`, 5 vs 5, source frames 13–129), and a gap check (`1s14B00`, 2 vs 2).
- **Criteria met:** preparing input, completing the run, and explaining the limitations.
- **Criterion partly met:** identifying the **exact exclusion causes from `tables/exclusions.csv`**. The 10 excluded
  channel–windows were identifiable, but the exact cause (a numbering discontinuity 2093 → 2117) required the
  per-member JSON. Tool v0.3.0 addresses this (`../dev_ensemble_report/results_v3/`).
- **Repository preservation:** verified by the reviewer. 570 tracked files were checked with no changes; `git status`
  and HEAD were unchanged.

## Contents and integrity

- **Archived files:**
  - `REVIEW.md`, `REVIEW.json`, `REVIEW_ARTIFACT_HASHES.json`;
  - `baseline.json` (the repository-preservation baseline);
  - `independent_mapping_check.json`, `reviewer_design.json`;
  - the two reviewer specs;
  - `*_run.json` (exact commands, working directory, environment overrides), with stdout/stderr;
  - per run: `report.md`, `manifest.json`, `tables/exclusions.csv`, `tables/between.csv`.
- **Integrity:** all 28 archived files match the SHA-256 values in `REVIEW_ARTIFACT_HASHES.json` (checked
  2026-10-23).
- **Left outside Git** (112 files, recorded by hash in `REVIEW_ARTIFACT_HASHES.json`): `distributions.npz`, the other
  tables and the per-member `members/` intermediates.
- **Regeneration:** run the command in each `*_run.json` from the repository root at commit `389f2d7`, with the same
  local mdCATH files (pinned SHA-256 in the specs' sources and in `manifest.json`).
