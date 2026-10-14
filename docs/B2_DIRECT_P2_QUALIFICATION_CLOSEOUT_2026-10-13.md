# BQ-3 Closeout Addendum (2026-10-13)

**Closes:** BQ-3, the direct P2 qualification of unchanged V2b (`docs/B2_DIRECT_P2_QUALIFICATION_PROTOCOL.md`;
report `benchmarks/audit/bq3/REPORT_SYNTHETIC_P2.md`; starting checkpoint `f55def1`). Nothing was rerun.

**Unchanged:** the protocol, acceptance code, measurements, analysis, report and the original run record. BQ-2's
historical **conclusion C** (insufficient evidence) also stands. BQ-3 is a separate direct measurement and does not
convert it.

## 1. Accepted outcome (scoped)

- **Prepared-database operation qualified** on the recorded SYNTHETIC_tiled, performance-only P2 workload
  (n_D = 10,000, L = 1,000, W = 8, 993 windows, 3 queries, 5 fresh processes). All 15 query latencies were within
  10 s, peak RSS was within 8 GiB, and verification passed.
- **Preparation plus the first query did not meet 10 s** in any repetition.
- **What the first-query span covers.** The span (`first_query_span_s`) runs from the start of preparation to the end
  of the first query. It **excludes**:
  - process start-up and imports;
  - the SYNTHETIC warm-up;
  - input loading and construction of the 10,000 + 3 tiled rows, recorded separately as `build_rows_wall_s` in
    `results/bq3_SYNTHETIC_measurements.json`.

  With real data, loading would add stage S1/S2 costs, which EB-1 evaluates separately under R3.
- **Scope:** the rows repeat real frames, so the ranking check has near ties everywhere (rankings were strictly
  identical to the reference). The results characterize one machine and session. They establish no biological
  capability, no novelty over existing RMSD software and no production readiness.

## 2. Editorial correction of the run record

- **Defect.** In the original `results/bq3_SYNTHETIC_run_record.json`:
  - `freeze_hashes_unchanged` is `true` although `make_report.py` differs from its pre-measurement hash;
  - `changed_after_pre_measurement` uses `true` to mean "hash matches", the opposite of its name.
- **What was preserved.** The original file is kept byte-for-byte. Its SHA-256 is recorded in the correction file.
- **Correction file:** `results/bq3_SYNTHETIC_run_record_correction_2026-10-13.json`, with explicitly named fields:
  - `all_pre_measurement_hashes_unchanged` (false);
  - `files_matching_pre_measurement_hash`;
  - `files_differing_from_pre_measurement_hash` (`make_report.py` only);
  - `decision_and_scientific_files_match_pre_measurement_hash` (true);
  - the before and after hashes of `make_report.py`;
  - output-file hashes compared with the original record.
- **History disclosure.** The run record was written twice before review:
  - The first version stated `freeze_hashes_unchanged = false` with `make_report.py` listed as changed. It was
    overwritten and is not preserved.
  - `make_report.py` was edited twice after measurement (report text only), and the intermediate hash was not
    recorded.

## 3. Next steps

- Neither integration option in `docs/B2_V2B_INTEGRATION_PROPOSAL.md` is implemented.
- The next step is review of the competitive-evaluation protocol draft,
  `docs/B2_COMPETITIVE_EVALUATION_PROTOCOL_DRAFT.md`.
