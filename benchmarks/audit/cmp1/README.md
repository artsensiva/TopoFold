# CMP-1: bounded competitive evaluation of window Cα RMSD retrieval (harness only)

- **Protocol:** [`docs/CMP1_COMPETITIVE_EVALUATION_PROTOCOL.md`](../../../docs/CMP1_COMPETITIVE_EVALUATION_PROTOCOL.md),
  frozen in `results/cmp1_protocol_freeze.json`. Code hashes are in `results/cmp1_pre_measurement_record.json`.
- **Report:** [`REPORT.md`](REPORT.md).
- **Closeout and clarifications:** [`docs/CMP1_CLOSEOUT_2026-10-13.md`](../../../docs/CMP1_CLOSEOUT_2026-10-13.md) (degeneracy-eligibility wording, scoped MDTraj recommendation, memory, CPU accounting).
- **Scope:** harness code only. No production integration, held-out scoring or raw-data download.

| File | Role |
|---|---|
| `requirements_cmp1.in` | Top-level pins of the isolated environment `build/venv` (offline install from a hashed wheelhouse) |
| `cmp1_data.py` | Canonical read-only float64 Å input (real S2 arrays; SYNTHETIC_tiled rows recast through float32) |
| `cmp1_adapters.py` | T-V2b, T-V0 (harness), E-MDT (MDTraj 1.11.2, pre-sliced precentered windows), E-MDA (MDAnalysis 2.10.0 qcprot, prepared float64 windows) |
| `cmp1_accuracy.py` | F64/F32 accuracy tracks, reference, degeneracy suite |
| `cmp1_criteria.py` | Eligibility, E\* per workload and track, outcomes, conclusion |
| `cmp1_devchecks.py` | Adapter-development checks (not deciding) → `results/cmp1_adapter_dev_checks.json` |
| `cmp1_workloads.py` | Frozen workloads and query IDs → `results/cmp1_workloads.json` |
| `cmp1_run.py`, `cmp1_analyse.py`, `make_report.py` | Run, analysis and report |
| `test_cmp1.py` | Decision-logic failure paths and track functions |

- **Environment provenance:** `results/cmp1_environment.json` (install commands, wheel SHA-256 values, MDTraj PyPI
  digest match, inspected source hashes).
- **Git-ignored:** `build/` (venv, wheelhouse, reference arrays).

## Post-hoc diagnostic (not deciding)

`results/cmp1_posthoc_mdtraj_error_diagnostic.json` was computed after the analysis; it took 1.3 s of CPU. On R1_D1,
MDTraj exceeds the 1e-3 Å F32 distance tolerance on **one pair per query out of 210,105**.
- Those pairs have reference RMSDs of about 1.8–2.1 Å, not near zero.
- None of them is in any top-10.
- This is consistent with MDTraj's zero rank non-equivalences on every workload.

It does not change any frozen outcome.
