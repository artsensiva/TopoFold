# WP-D: Cα-only local retrieval vs a full-backbone φ/ψ target

Audit harness, not production code. Protocol: `docs/WPD_LOCAL_RETRIEVAL_PROTOCOL.md` (§10 amendments
A1–A5), frozen in `results/wpd_protocol_freeze.json` before any feature or score was computed.

## Data labelling

- AdK outputs (`results/wpd_adk_THIRD_PARTY_MD_DERIVATIVE_*.json`, `REPORT_ADK.md`) are **derived from
  third-party MD** (Orädd et al. 2021, Zenodo 10.5281/zenodo.5583119, CC BY 4.0). They are not SYNTHETIC
  and not TopoFold-generated trajectories.
- Test fixtures in `test_wpd.py` are **SYNTHETIC** (random matrices, constructed helices); no AdK data.
- The AdK tier is **descriptive only**; it takes no retain/reject decision and is not held-out validation.
  The mdCATH arm (§7–§8) is provisional and not executed.

## Files

| File | Role |
|------|------|
| `adk_backbone_eligibility.py` | Atom availability, finite coordinates, bond continuity (read-only; earlier task) |
| `wpd_lib.py` | Dihedrals, degeneracy sines, distances (M1, M2, B1, B2 Horn, B3, C1), tie-aware membership, overlap/regret, Spearman |
| `run_adk.py` | Input binding (protocol/eligibility/WP-B G0/binding identity), geometry, per-window retrieval in query batches, aggregation |
| `test_wpd.py` | SYNTHETIC correctness tests (A3 list plus degeneracy and binding indexing) |
| `make_report.py` | Writes `REPORT_ADK.md` from the results JSON |
| `adk_rule_illustration.py` | Applies the held-out criterion to AdK windows as a point-estimate illustration (not a decision) |
| `mdcath_select.py`, `test_mdcath_select.py` | Held-out tier: metadata-only eligibility, schema sample and seeded draw; writes `HELDOUT_DOWNLOAD_REQUEST.md` |

## Held-out tier (mdCATH)

Protocol `docs/WPD_HELDOUT_MDCATH_PROTOCOL.md`, frozen in `results/wpd_heldout_protocol_freeze.json` before the
draw. `results/wpd_mdcath_selection.json` lists the schema sample (development data), 30 primary and 10 reserve
domains with sizes and LFS SHA-256. *(Status 2026-10-12: 30 primary + 6 reserve files downloaded and scored; see below.)* Metadata files
(CATH 4.2.0 lists, mdCATH tree listing) are in `build/metadata/` with hashes in the freeze record.

Amendment 1 (`docs/WPD_HELDOUT_MDCATH_AMENDMENT_1.md`, record `results/wpd_heldout_amendment1.json`) replaces the
loss estimand, defines invalid-denominator and C1 coverage-matched handling, relabels classes and corrects the
frame-count arithmetic. Implemented in `heldout_decision.py`, tested in `test_heldout_decision.py` (SYNTHETIC).
The original protocol, freeze and selection are unchanged.

Amendment 2 (approved 2026-10-06, `docs/WPD_HELDOUT_MDCATH_AMENDMENT_2.md`, record `results/wpd_heldout_amendment2.json`) adds
method A and the separate peptide-plane target, and defines dependence diagnostics (`heldout_dependence.py`,
calibration `calibrate_tau_int.py` → `results/SYNTHETIC_tau_int_calibration.json`, tests
`test_heldout_dependence.py`; candidate set `heldout_decision_a2.py`). Minimal reader `mdcath_reader.py`
(tests `test_mdcath_reader.py`) verified on the schema sample only by `mdcath_schema_check.py` →
`results/wpd_mdcath_schema_check.json`. The schema sample (development data) is in `data/mdcath/` (git-ignored);
h5py 3.16.0 was added to `build/venv` (`build/pip_freeze_after_h5py.txt`).

Amendment 3 (approved 2026-10-06, `docs/WPD_HELDOUT_MDCATH_AMENDMENT_3.md`, record `results/wpd_heldout_amendment3.json`): mandatory checks exit 2, optional O, rounding after reconstruction, dependence r1 (AR(1)-based heuristic only), cap correction, memory scaling, runner `heldout_runner.py` (tests `test_heldout_runner.py`; dev validation `results/wpd_heldout_runner_dev_validation.json`). The 30 primary files were downloaded and verified (`mdcath_download.py` → `results/wpd_mdcath_download_manifest.json`); structural eligibility in `results/wpd_mdcath_eligibility_ledger.json` and `results/wpd_mdcath_eligibility_notes.json`. **Held-out evaluation completed 2026-10-06** under `docs/WPD_HELDOUT_MDCATH_EXECUTION_CORRECTIONS.md` (record `results/wpd_heldout_execution_corrections.json`; authorization `results/wpd_heldout_scoring_authorization.json`): cohort `results/wpd_mdcath_cohort.json` (24 primaries + 6 reserves; ledger `results/wpd_mdcath_reserve_ledger.json`), results `results/wpd_heldout_THIRD_PARTY_MD_DERIVATIVE_results.json`, report `REPORT_HELDOUT.md` (`make_report_heldout.py`), run record `results/wpd_heldout_run_record.json`. Execution tests: `test_heldout_exec.py`; dev command-path verification `results/wpd_heldout_exec_verification.json`.

## Run

```
build/venv/bin/python -B -m pytest -p no:cacheprovider test_wpd.py
build/venv/bin/python -B run_adk.py --workers 4      # runs the tests first; stops on any failed check
build/venv/bin/python -B make_report.py
```

`build/venv` is an isolated environment with the WP-A wheel (`topofold-0.8.3`, extension SHA-256 checked at
run time against the WP-A record) plus NumPy, MDAnalysis, SciPy and pytest (`build/pip_freeze.txt`). It is
git-ignored. Raw AdK files are read from `../adk_feasibility/data/raw/` and bound to the WP-B G0 hashes.

## Closeout (2026-10-12)

Interpretation corrections (no rerun, results unchanged): `docs/WPD_HELDOUT_MDCATH_CLOSEOUT_ADDENDUM_2026-10-12.md`,
generated by `make_closeout_addendum.py` → `results/wpd_heldout_closeout_addendum.json`. A and B2 are retained as
candidates for their respective geometric tasks only. Next-step proposals: `docs/ENGINEERING_BENCHMARK_EB1_PROPOSAL.md`,
`docs/ONCOLOGY_VALIDATION_REQUIREMENTS.md`.
