# BQ-3: direct P2 qualification of unchanged V2b

**Preflight (2026-10-13, first attempt).** The memory preflight stopped the attempt; no measurement was made.
- `bq3_memory_preflight.py` writes `results/bq3_memory_preflight.json`. That record is kept unchanged; running the
  script again would overwrite it.
- Details: [`docs/B2_DIRECT_P2_QUALIFICATION_PREFLIGHT.md`](../../../docs/B2_DIRECT_P2_QUALIFICATION_PREFLIGHT.md).

**Direct qualification (2026-10-13, from checkpoint `f55def1`).** Run after the memory recheck. It uses
SYNTHETIC_tiled rows and is performance-only.
- Protocol: [`docs/B2_DIRECT_P2_QUALIFICATION_PROTOCOL.md`](../../../docs/B2_DIRECT_P2_QUALIFICATION_PROTOCOL.md).
- Report: [`REPORT_SYNTHETIC_P2.md`](REPORT_SYNTHETIC_P2.md).
- Closeout and run-record correction: [`docs/B2_DIRECT_P2_QUALIFICATION_CLOSEOUT_2026-10-13.md`](../../../docs/B2_DIRECT_P2_QUALIFICATION_CLOSEOUT_2026-10-13.md); explicit-field correction record `results/bq3_SYNTHETIC_run_record_correction_2026-10-13.json` (the original run record is preserved unchanged).
- Integration proposal (review only): [`docs/B2_V2B_INTEGRATION_PROPOSAL.md`](../../../docs/B2_V2B_INTEGRATION_PROPOSAL.md).

| File | Role |
|---|---|
| `bq3_lib.py` | Frozen acceptance rule (verification, completeness, per-observation limits) |
| `test_bq3.py` | Failure-path tests of the acceptance rule |
| `bq3_run.py` | Tests → reference → gate → 5 timing repetitions (fresh processes, memory start rule) |
| `bq3_analyse.py`, `make_report.py` | Verdicts and report from the results files |

**Records**
- Freeze and pre-measurement: `results/bq3_protocol_freeze.json`, `results/bq3_pre_measurement_record.json`.
- Outputs: `results/bq3_SYNTHETIC_measurements.json`, `results/bq3_SYNTHETIC_analysis.json`,
  `results/bq3_SYNTHETIC_run.log`.
- Reference distances: `build/` (git-ignored).
