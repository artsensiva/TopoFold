# SYNTHETIC_NULLSTRESS-1: identical marginals, different temporal correlation (SYNTHETIC)

Generated from `SYNTHETIC_NULLSTRESS-1_results.json`. Status **complete** (6 of 6 cells). α = 0.05; nominal permutation size 12/252 ≈ 0.048 under exact exchangeability, which these cells deliberately violate (the marginals are equal; the temporal correlation is not). 200 repetitions per cell, which gives coarse Monte Carlo precision. No real data.

Resources: 23.9 CPU-s (cap 120), wall 23.2 s.

| Channel | φ_A | φ_B | Permutation reject [95 % CI] | Screen C | Screen RD | Screen I | median within-A / within-B JS distance |
|---|---|---|---|---|---|---|---|
| linear | 0.0 | 0.8 | 0.125 [0.086, 0.178] | 0.000 | 0.000 | 1.000 | 0.094 / 0.135 |
| linear | 0.0 | 0.95 | 0.520 [0.451, 0.588] | 0.000 | 0.005 | 0.995 | 0.095 / 0.241 |
| linear | 0.8 | 0.95 | 0.135 [0.094, 0.189] | 0.000 | 0.000 | 1.000 | 0.135 / 0.231 |
| circular | 0.0 | 0.8 | 0.115 [0.078, 0.167] | 0.000 | 0.000 | 1.000 | 0.102 / 0.143 |
| circular | 0.0 | 0.95 | 0.525 [0.456, 0.593] | 0.000 | 0.000 | 1.000 | 0.103 / 0.246 |
| circular | 0.8 | 0.95 | 0.085 [0.054, 0.132] | 0.000 | 0.000 | 1.000 | 0.145 / 0.250 |
