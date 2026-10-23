# SYNTHETIC_NULLSTRESS-1: identical marginals, different temporal correlation (SYNTHETIC)

A small, separately frozen null stress check of the PERM-1 replica permutation method and the complete-separation
screen.
- **Design:** both conditions draw from the same stationary marginal distribution (linear or circular), but their
  replicas have different AR(1) correlation: (φ_A, φ_B) ∈ {(0, 0.8), (0, 0.95), (0.8, 0.95)}. 5 replicas per
  condition, 200 repetitions.
- **Specification:** `docs/NULL_STRESS_CHECK_SPEC.md` and `run_spec.json`, frozen with code and environment
  (`results/SYNTHETIC_NULLSTRESS-1_freeze.json`).
- **Data:** no real data.

```bash
.venv/bin/python benchmarks/audit/synth_null_stress/run_stress.py run
```

## Result

All figures are from `results/SYNTHETIC_NULLSTRESS-1_results.json`; Wilson intervals are in the report.

| Correlation pair | Permutation rejection (nominal 0.048) |
|---|---|
| φ 0 vs 0.8 | 0.125 linear, 0.115 circular |
| φ 0 vs 0.95 | 0.520 linear, 0.525 circular |
| φ 0.8 vs 0.95 | 0.135 linear, 0.085 circular |

- **Why the permutation method rejects:** more correlated replicas give more variable histograms. Median
  within-condition JS distances were roughly 0.10 at φ = 0 against 0.23–0.25 at φ = 0.95, so the replicas are not
  exchangeable even though the marginals are identical.
- **Screen:** no consistent calls and one replica-dependent call in 1,200 draws.

**Implication:** a permutation rejection means the replica-level histograms differ (marginal geometry and/or temporal
structure). It does not by itself mean the geometry differs. The real-data design
(`docs/DEV_PERM_1HA8A00_EVALUATION_DESIGN_DRAFT.md`) states this reading and adds predeclared dynamics diagnostics.
