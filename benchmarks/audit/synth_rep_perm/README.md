# SYNTHETIC_PERM-1: complete-separation screen vs replica permutation method (SYNTHETIC)

Prospective evaluation of one method, `topofold_inputs.replica_permutation.permutation_test`, against the historical
complete-separation screen.
- Data: the 81 five-versus-five REPCHAR-1 SYNTHETIC cells, with original identifiers, seeds and draws, 200
  repetitions each. The three S6 cells are deferred.
- Specification: `docs/REPLICA_PERMUTATION_EVALUATION_SPEC.md` and `run_spec.json`, frozen with code and environment
  (`results/SYNTHETIC_PERM-1_freeze.json`) before the run.
- Both methods were evaluated on identical draws. The screen reproduced the REPCHAR-1 category counts in all 81
  cells.
- No real data were used.

| File | Content |
|---|---|
| `run_perm.py` | `freeze` / `run`; hard 900 CPU-s cap; checkpoints stay `partial` until all planned cells finish |
| `results/SYNTHETIC_PERM-1_results.json`, `results/SYNTHETIC_PERM-1_report.md` | Per-cell screen categories, permutation rejection frequency, agreement counts, Wilson intervals |

```bash
.venv/bin/python benchmarks/audit/synth_rep_perm/run_perm.py run
```

## What we learned

All figures are from the results file; 200 repetitions per cell, so precision is coarse. The permutation rejection
frequency is shown as a range over cells, at α = 0.05 with a nominal size of 12/252 ≈ 0.048.

| Scenario role | Permutation rejection | Screen calls (C + RD) |
|---|---|---|
| Primary null S0 (6 cells) | 0.035–0.060 (7–12/200); every Wilson interval contains 0.048, including at φ = 0.95 | 0 |
| Heterogeneous exchangeable null S4b, Δ = 0 (6 cells) | 0.050–0.080; 5 of 6 intervals contain 0.048 (the exception, 16/200, has lower bound 0.050) | 0–0.005 |
| Single-replica perturbation S4a Δ = 0, S4c (6 cells) | 0.025–0.055 | 0 |
| Condition-confounded sampling S5 (6 cells) | 0.235–0.870 | 0.015–0.690 |

**Null behaviour.**
- On these generators and replica counts, rejection under the primary null matched the nominal level, unaffected by
  within-replica temporal correlation. That is expected when replicas, not frames, are the exchangeable unit.
- This is per-evaluation behaviour on SYNTHETIC data. It is not general error control and not a multiple-window
  procedure.

**Sensitivity.** In every condition-effect cell the permutation rejection frequency was at least the screen's call
frequency. Examples (permutation rejection vs screen C + RD):
- S1 linear, Δ = 0.25σ, φ = 0: 1.000 vs 0.030 + 0.165;
- Δ = 0.5σ, φ = 0.8: 0.995 vs 0.045 + 0.260;
- Δ = 1σ, φ = 0.95: 0.990 vs 0.075 + 0.245;
- slow population shift S3, w = 0.3, D = 20: 0.455 vs 0 + 0.015.

The screen's calls were nested within the permutation's rejections: no screen call (consistent or replica-dependent)
ever occurred without a permutation rejection on the same draw.

**Not resolved by the method.**
- Condition-confounded sampling (S5) was rejected more often than the screen called it (up to 0.870). The method
  responds to realized distribution differences, whatever their cause; only the design (independent starts per
  condition, and even then not guaranteed) can address attribution.
- A single unusual replica in one condition (S4a Δ = 0, S4c) was rejected at about the nominal level. Neither method
  flags it, and the screen's replica-dependent category remains the only signal of replica influence (in effect
  cells).
- With 5 vs 5 replicas the smallest attainable p is 2/252 ≈ 0.008. That limits any correction for many windows.

## Recommendation

The method **warrants a separately designed real-data evaluation**, not adoption. That design must:
- state why replicas are exchangeable across conditions (mdCATH replicas share one start across temperatures);
- control for many windows (for example FDR with a replica-count-aware minimum p);
- keep the complete-separation screen as a reported companion;
- treat rejections as distribution differences, not causal effects.
