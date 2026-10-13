# EB-1: bounded engineering benchmark (development data only)

Frozen protocol: `docs/EB1_ENGINEERING_BENCHMARK_PROTOCOL.md` (`results/eb1_protocol_freeze.json`). Harness, analysis
and tests hashed before measurement (`results/eb1_pre_measurement_record.json`). Report: `REPORT.md`.

| File | Role |
|------|------|
| `eb1_lib.py` | Stage implementations (existing readers/functions only) and correctness helpers |
| `eb1_run.py` | Runner: freeze check → cache prep → mandatory correctness → measurement matrix (fresh subprocess per cell) |
| `eb1_analyse.py` | Fits, support rules, P2 projections, requirements and the frozen decision rule |
| `test_eb1.py` | SYNTHETIC tests of the analysis rules and top-k/near-tie comparison |
| `make_report.py` | `REPORT.md` from results (includes a labelled post-hoc section) |

Results: `results/eb1_measurements.json`, `results/eb1_analysis.json`, `results/eb1_run_record.json`, log
`results/eb1_run.log`. Stage-input caches are in `build/` (git-ignored). Run order:
`eb1_run.py` → `eb1_analyse.py` → `make_report.py` (venv `../wpd_local_retrieval/build/venv`).

## Amendment 1 (2026-10-12)

`docs/EB1_AMENDMENT_1.md` (freeze `results/eb1a1_freeze.json`); code `eb1a1_lib.py`, `eb1a1_run.py`,
`eb1a1_analyse.py`, tests `test_eb1a1.py`, SYNTHETIC calibration `results/eb1a1_temp_calibration.json`; outputs
`results/eb1a1_measurements.json`, `results/eb1a1_analysis.json`, `results/eb1a1_run_record.json`, `REPORT_A1.md`
(`make_report_a1.py`). Original EB-1 files are unchanged.

## Amendment 2 (2026-10-12)

`docs/EB1_AMENDMENT_2.md` (freeze `results/eb1a2_freeze.json`); read-only preflight `eb1a2_preflight.py` →
`results/eb1a2_preflight.json`; code `eb1a2_lib.py` (direct-read loaders, released repetitions), `eb1a2_run.py`,
`eb1a2_analyse.py`, tests `test_eb1a2.py`; outputs `results/eb1a2_measurements.json`, `results/eb1a2_analysis.json`,
`results/eb1a2_run_record.json`, log `results/eb1a2_run.log`, `REPORT_A2.md` (`make_report_a2.py`). EB-1 and
Amendment 1 files are unchanged. Run order: `eb1a2_preflight.py` → `eb1a2_run.py` → `eb1a2_analyse.py` →
`make_report_a2.py`.
