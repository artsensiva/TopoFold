# DEV-ENS-1: descriptive per-window comparison of two AdK trajectories (development)

Per-window comparison of local Cα geometry (κ, τ, and the equal-information Cα(i)–Cα(i+2) and Cα(i)–Cα(i+3)
distances) between the AdK open-start and closed-start trajectories, over the matched elapsed interval 0–335.2 ns,
using histogram Jensen–Shannon distances (base 2) and three equal-duration blocks per run.

**Evidence level:** descriptive development comparison. One run per starting condition, so every difference mixes
the starting structure with relaxation from it. Blocks are not replicates. No test, state assignment, population,
causal or biological claim; no claim that one representation is better. Not a reproduction of PENSA.

| File | Content |
|---|---|
| `run_spec.json` | Specification (period, blocks, measure, bins, sensitivity, channels), written before execution |
| `results/DEV-ENS-1_freeze.json` | Spec and code SHA-256 recorded before the run; the runner refuses to run on mismatch |
| `run_compare.py` | `freeze` / `run`: verifies hashes, runs the Cα pipeline on both runs, checks timestamps and mapping, compares |
| `results/DEV-ENS-1_results.json` | Checks, blocks, per-window distances (primary), summaries, rank stability, sensitivity, resources |
| `results/DEV-ENS-1_report.md` | Generated report |
| `posthoc_block_pattern.py`, `results/DEV-ENS-1_posthoc_block_pattern.json` | **Post-hoc**, descriptive: which block pairs drive the block fractions of the top-10 windows |
| `build/` (git-ignored) | Coordinates and features of both runs |

Reproduce (local AdK files from the WP-B G0 record; WP-D audit venv for the adapter stage):

```bash
.venv/bin/python benchmarks/audit/dev_ens_compare/run_compare.py run
.venv/bin/python benchmarks/audit/dev_ens_compare/posthoc_block_pattern.py
```
