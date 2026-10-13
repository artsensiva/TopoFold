# BQ-2: scaling confirmation of unchanged V2b (harness only)

This study uses development data plus SYNTHETIC_tiled rows, which are performance-only. It changes no production code
and leaves BQ-1 (`../b2opt/`) unchanged.

**Protocol and freezes**
- `docs/B2_SCALING_CONFIRMATION_PROTOCOL.md`; frozen in `results/bq2_protocol_freeze.json`, together with the
  workloads.
- `results/bq2_pre_measurement_record.json`: code hashes recorded before measurement.

| File | Role |
|------|------|
| `bq2_workloads.py` | Fit grid, fresh real confirmations, SYNTHETIC_tiled stress cells, direct P2 cell → `results/bq2_workloads.json` |
| `bq2_lib.py` | Tiled stress rows, preparation repetitions with release, memory start rule |
| `bq2_run.py` | Tests → references → gates → interleaved timing rounds → P2D if memory permits |
| `bq2_analyse.py` | Direct measurements, models, confirmation, P2 projections, conclusion A/B/C |
| `test_bq2.py` | SYNTHETIC regression tests |
| `make_report.py` | `REPORT.md` from the results files |

**Run order** (venv `../wpd_local_retrieval/build/venv`): `bq2_run.py` → `bq2_analyse.py` → `make_report.py`.

**Outputs**
- `results/bq2_measurements.json`, `results/bq2_analysis.json`, `results/bq2_run.log`
- `build/ref/` (git-ignored)
