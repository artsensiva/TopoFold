# SYNTHETIC_DABOOT-1: stationary-bootstrap Kuiper test under unequal temporal correlation (SYNTHETIC)

Generated from `SYNTHETIC_DABOOT-1_results.json`. Status **complete** (10 of 10 cells). α = 0.05, B = 499 bootstrap draws per test, 200 repetitions per cell, which gives coarse Monte Carlo precision. Null draws are identical to SYNTHETIC_NULLSTRESS-1 (checked: True). No real data.

Resources: 260.7 CPU-s (cap 900), wall 260.0 s.

| Cell | Role | Channel | φ_A | φ_B | Δ | Reject [95 % CI] | median block length A / B | NULLSTRESS-1 permutation (same draws) |
|---|---|---|---|---|---|---|---|---|
| 0 | null | linear | 0.0 | 0.8 | 0.0 | 0.095 [0.062, 0.144] | 1.0 / 15.4 | 0.125 |
| 1 | null | linear | 0.0 | 0.95 | 0.0 | 0.230 [0.177, 0.293] | 1.0 / 32.9 | 0.520 |
| 2 | null | linear | 0.8 | 0.95 | 0.0 | 0.180 [0.133, 0.239] | 15.4 / 33.2 | 0.135 |
| 3 | null | circular | 0.0 | 0.8 | 0.0 | 0.095 [0.062, 0.144] | 1.0 / 15.5 | 0.115 |
| 4 | null | circular | 0.0 | 0.95 | 0.0 | 0.195 [0.146, 0.255] | 1.0 / 34.9 | 0.525 |
| 5 | null | circular | 0.8 | 0.95 | 0.0 | 0.200 [0.150, 0.261] | 15.7 / 33.2 | 0.085 |
| 6 | positive_control | linear | 0.0 | 0.95 | 1.0 | 1.000 [0.981, 1.000] | 1.0 / 32.6 |  |
| 7 | positive_control | linear | 0.0 | 0.95 | 2.0 | 1.000 [0.981, 1.000] | 1.0 / 32.4 |  |
| 8 | positive_control | circular | 0.0 | 0.95 | 1.0 | 1.000 [0.981, 1.000] | 1.0 / 33.6 |  |
| 9 | positive_control | circular | 0.0 | 0.95 | 2.0 | 1.000 [0.981, 1.000] | 1.0 / 34.3 |  |

## Predeclared assessment

- Calibration: **failure** (failing null cells: [0, 1, 2, 3, 4, 5]; conservative cells ≤ 0.01: none). Calibration failure = a null cell above 0.10 or with a Wilson lower bound above 0.05. An interval containing 0.05 is not proof of validity.
- Sensitivity (both Δ = 2 controls ≥ 0.80): **useful**.
- Decision by the predeclared rule: **calibration failure: retain descriptive comparisons; record limitation; stop**.
