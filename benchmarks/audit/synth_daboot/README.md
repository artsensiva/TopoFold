# SYNTHETIC_DABOOT-1: stationary-bootstrap Kuiper test of equal marginals (SYNTHETIC) — calibration failure

Bounded prototype characterization of one established dependence-aware method. The method is a centered stationary
bootstrap (Politis & Romano 1994), with block length by Politis–White (2004) as corrected by Patton–Politis–White
(2009), applied with a two-sample Kuiper statistic. It is the method selected in
`docs/DEPENDENCE_AWARE_MARGINAL_COMPARISON_SPEC.md` before any outcome.
- **Replica boundaries:** each replica is resampled separately and never joined to another.
- **Cells:** the six NULLSTRESS-1 null cells (identical marginals, unequal temporal correlation; identical draws,
  checked) and four positive controls (Δ = 1, 2; φ 0 vs 0.95). 200 repetitions each; B = 499.
- **Freeze:** spec, code and environment were frozen before the run (`results/SYNTHETIC_DABOOT-1_freeze.json`).
- **Data:** no real data.

```bash
.venv/bin/python benchmarks/audit/synth_daboot/run_daboot.py run
```

## Result (from `results/SYNTHETIC_DABOOT-1_results.json`)

- **Null cells:** rejection 0.095–0.230 at α = 0.05.
  - φ 0 vs 0.8: 0.095 (linear and circular, Wilson lower bound 0.062).
  - φ 0 vs 0.95: 0.230 linear, 0.195 circular.
  - φ 0.8 vs 0.95: 0.180 linear, 0.200 circular.
  - Every null cell meets the predeclared failure condition: above 0.10, or Wilson lower bound above 0.05.
- **Positive controls:** 1.000 in all four cells.
- **Context:** on the same draws, NULLSTRESS-1's permutation rejected 0.085–0.525.
  - The bootstrap reduced over-rejection where correlation differed most from zero (0.23 vs 0.52 at φ 0 vs 0.95).
  - It was higher at φ 0.8 vs 0.95 (0.18–0.20 vs 0.085–0.135).
  - It did not reach calibration in any null cell.
- **Predeclared decision:** calibration failure, so retain descriptive comparisons, record the limitation and stop.
  No method switch and no retuning.

## Recorded limitation

- **The calibration failure:** with 5 replicas × 387 frames, the stationary bootstrap of each condition's empirical
  distribution under-represents the sampling variability of strongly correlated series, so it over-rejects equal
  marginals.
- **Block length (descriptive, not a fix):** the automatic block lengths at φ = 0.95 had a median of about 33
  frames, below the theoretical optimum of about 53 for n = 387. Even an optimal block length leaves only about 7
  blocks per replica.
- **Status:** the prototype (`python/topofold_inputs/dependent_bootstrap.py`) is retained as an audit artifact and is
  not recommended for inference.
