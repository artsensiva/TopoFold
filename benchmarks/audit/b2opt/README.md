# BQ-1: B2 single-query optimization experiment (harness only)

This experiment uses development data and SYNTHETIC tests. It changes no production code. EB-1 versions and outcomes
are unchanged.

**Protocol and freezes**
- `docs/B2_SINGLE_QUERY_OPTIMIZATION_PROTOCOL.md`; frozen in `results/b2opt_protocol_freeze.json`, together with the
  workloads.
- `docs/B2_SINGLE_QUERY_OPTIMIZATION_AMENDMENT_1.md` (adds arm V2b).
- `results/b2opt_pre_measurement_record.json`: code hashes recorded before any real-data computation.

| File | Role |
|------|------|
| `b2opt_workloads.py` | Reused grid vs fresh confirmation queries (read-only) → `results/b2opt_workloads.json` |
| `b2opt_lib.py` | Arms V0/V1/V2/V2b, frozen top-k tie policy, ranking check, SVD reference, SYNTHETIC gate suite, transient calibration |
| `b2opt_run.py` | Tests → reference check → gates → interleaved timing rounds (fresh subprocess, one thread) |
| `b2opt_analyse.py` | Direct measurements, speed-ups, models, confirmation, P2 projections, acceptance |
| `test_b2opt.py` | SYNTHETIC regression tests: numerical correctness, ties, degenerate geometry, allocation lifetimes |
| `make_report.py` | `REPORT.md` from the results files |

**Run order** (venv `../wpd_local_retrieval/build/venv`): `b2opt_run.py` → `b2opt_analyse.py` → `make_report.py`.

**Outputs**
- `results/b2opt_measurements.json`, `results/b2opt_analysis.json`, `results/b2opt_run.log`
- Reference distances in `build/ref/` (git-ignored by the local `.gitignore`)
