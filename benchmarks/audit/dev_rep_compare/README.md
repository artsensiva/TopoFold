# DEV-REP-1: replicated descriptive comparison, mdCATH 1ha8A00 (320 K vs 348 K)

Per-window comparison of local Cα geometry (κ, τ and the equal-information Cα(i)–Cα(i+2), Cα(i)–Cα(i+3)
distances) with the replica as the unit: 5 replicas per temperature, frames 43–429 of each (frames 0–42 excluded as
predeclared), histogram Jensen–Shannon distances (base 2), leave-one-replica-out sensitivity and predeclared bin
sensitivity. Data: local mdCATH development domain `1ha8A00` (CC BY 4.0), read through
`python/topofold_inputs/mdcath_bridge.py`; the protected WP-D cohort is not used.

**Evidence level:** descriptive development evidence. Replicas share one starting configuration; pairwise distances
are dependent; categories are not tests. 320 K vs 348 K is a simulation condition, not an experimental endpoint.
Exposure history of the domain is listed in `run_spec.json`.

| File | Content |
|---|---|
| `run_spec.json` | Specification frozen before scoring (sampling, exclusion, bins, budget, replica-level rules, exposure) |
| `results/DEV-REP-1_freeze.json` | Spec and code SHA-256; the runner refuses to run on mismatch |
| `run_rep.py` | `freeze` / `run` |
| `results/DEV-REP-1_results.json`, `results/DEV-REP-1_report.md` | Results and generated report |
| `posthoc_summary.py`, `results/DEV-REP-1_posthoc_summary.json` | **Post-hoc**, descriptive distribution summaries |
| `INTERPRETATION_NOTE_2026-10-23.md` | Closeout: outcome stated as "no window met the complete-separation criterion" (194 channel–window evaluations); how to reproduce DEV-ENS-1 from `644b4cd` |
| `build/` (git-ignored) | Per-replica coordinates and features |

```bash
.venv/bin/python benchmarks/audit/dev_rep_compare/run_rep.py run
.venv/bin/python benchmarks/audit/dev_rep_compare/posthoc_summary.py
```
