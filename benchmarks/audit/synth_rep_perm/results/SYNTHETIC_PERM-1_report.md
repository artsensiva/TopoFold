# SYNTHETIC_PERM-1: complete-separation screen vs replica permutation (SYNTHETIC)

Generated from `SYNTHETIC_PERM-1_results.json`. Status **complete** (81 of 81 cells). α = 0.05; nominal permutation size 12/252 ≈ 0.048 under exact exchangeability. 200 repetitions per cell, which gives coarse Monte Carlo precision (Wilson 95 % intervals). Both methods use identical draws. Screen categories: C consistent, RD replica-dependent, I inconclusive. Per-evaluation SYNTHETIC results; not a multiple-window procedure. No real data.

Screen counts reproduce REPCHAR-1 in all cells: True. Resources: 293.8 CPU-s (cap 900), wall 293.2 s.

## Primary identical-distribution null (S0)

| Scenario | Cell | Screen C | Screen RD | Screen I | Permutation reject [95 % CI] | reject & C / RD / I | accept & C / RD |
|---|---|---|---|---|---|---|---|
| S0 | channel=linear, phi=0.0 | 0.000 | 0.000 | 1.000 | 0.045 [0.024, 0.083] | 0 / 0 / 9 | 0 / 0 |
| S0 | channel=linear, phi=0.8 | 0.000 | 0.000 | 1.000 | 0.040 [0.020, 0.077] | 0 / 0 / 8 | 0 / 0 |
| S0 | channel=linear, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.050 [0.027, 0.090] | 0 / 0 / 10 | 0 / 0 |
| S0 | channel=circular, phi=0.0 | 0.000 | 0.000 | 1.000 | 0.050 [0.027, 0.090] | 0 / 0 / 10 | 0 / 0 |
| S0 | channel=circular, phi=0.8 | 0.000 | 0.000 | 1.000 | 0.035 [0.017, 0.070] | 0 / 0 / 7 | 0 / 0 |
| S0 | channel=circular, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.060 [0.035, 0.102] | 0 / 0 / 12 | 0 / 0 |

## Heterogeneous exchangeable null (S4b, Δ = 0)

| Scenario | Cell | Screen C | Screen RD | Screen I | Permutation reject [95 % CI] | reject & C / RD / I | accept & C / RD |
|---|---|---|---|---|---|---|---|
| S4b | channel=linear, delta=0.0, k=0.25, phi=0.0 | 0.000 | 0.005 | 0.995 | 0.050 [0.027, 0.090] | 0 / 1 / 9 | 0 / 0 |
| S4b | channel=linear, delta=0.0, k=0.25, phi=0.8 | 0.000 | 0.000 | 1.000 | 0.080 [0.050, 0.126] | 0 / 0 / 16 | 0 / 0 |
| S4b | channel=linear, delta=0.0, k=0.25, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.050 [0.027, 0.090] | 0 / 0 / 10 | 0 / 0 |
| S4b | channel=linear, delta=0.0, k=0.5, phi=0.0 | 0.000 | 0.000 | 1.000 | 0.050 [0.027, 0.090] | 0 / 0 / 10 | 0 / 0 |
| S4b | channel=linear, delta=0.0, k=0.5, phi=0.8 | 0.000 | 0.005 | 0.995 | 0.065 [0.038, 0.108] | 0 / 1 / 12 | 0 / 0 |
| S4b | channel=linear, delta=0.0, k=0.5, phi=0.95 | 0.000 | 0.005 | 0.995 | 0.050 [0.027, 0.090] | 0 / 1 / 9 | 0 / 0 |

## Single-replica perturbation, not an ordinary null (S4a Δ = 0, S4c)

| Scenario | Cell | Screen C | Screen RD | Screen I | Permutation reject [95 % CI] | reject & C / RD / I | accept & C / RD |
|---|---|---|---|---|---|---|---|
| S4a | channel=linear, delta=0.0, phi=0.0 | 0.000 | 0.000 | 1.000 | 0.045 [0.024, 0.083] | 0 / 0 / 9 | 0 / 0 |
| S4a | channel=linear, delta=0.0, phi=0.8 | 0.000 | 0.000 | 1.000 | 0.050 [0.027, 0.090] | 0 / 0 / 10 | 0 / 0 |
| S4a | channel=linear, delta=0.0, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.055 [0.031, 0.096] | 0 / 0 / 11 | 0 / 0 |
| S4c | channel=linear, D=1 | 0.000 | 0.000 | 1.000 | 0.040 [0.020, 0.077] | 0 / 0 / 8 | 0 / 0 |
| S4c | channel=linear, D=20 | 0.000 | 0.000 | 1.000 | 0.025 [0.011, 0.057] | 0 / 0 / 5 | 0 / 0 |
| S4c | channel=linear, D=100 | 0.000 | 0.000 | 1.000 | 0.055 [0.031, 0.096] | 0 / 0 / 11 | 0 / 0 |

## Condition-confounded sampling, not a null (S5)

| Scenario | Cell | Screen C | Screen RD | Screen I | Permutation reject [95 % CI] | reject & C / RD / I | accept & C / RD |
|---|---|---|---|---|---|---|---|
| S5 | channel=linear, k=0.25, phi=0.0 | 0.295 | 0.110 | 0.595 | 0.735 [0.670, 0.791] | 59 / 22 / 66 | 0 / 0 |
| S5 | channel=linear, k=0.25, phi=0.8 | 0.020 | 0.070 | 0.910 | 0.535 [0.466, 0.603] | 4 / 14 / 89 | 0 / 0 |
| S5 | channel=linear, k=0.25, phi=0.95 | 0.000 | 0.015 | 0.985 | 0.235 [0.182, 0.298] | 0 / 3 / 44 | 0 / 0 |
| S5 | channel=linear, k=0.5, phi=0.0 | 0.615 | 0.075 | 0.310 | 0.870 [0.816, 0.910] | 123 / 15 / 36 | 0 / 0 |
| S5 | channel=linear, k=0.5, phi=0.8 | 0.220 | 0.135 | 0.645 | 0.740 [0.675, 0.796] | 44 / 27 / 77 | 0 / 0 |
| S5 | channel=linear, k=0.5, phi=0.95 | 0.040 | 0.065 | 0.895 | 0.530 [0.461, 0.598] | 8 / 13 / 85 | 0 / 0 |

## Condition effect (S1, S2, S3, S4a Δ = 1, S4b Δ = 1)

| Scenario | Cell | Screen C | Screen RD | Screen I | Permutation reject [95 % CI] | reject & C / RD / I | accept & C / RD |
|---|---|---|---|---|---|---|---|
| S1 | channel=linear, delta=0.1, phi=0.0 | 0.000 | 0.005 | 0.995 | 0.410 [0.344, 0.479] | 0 / 1 / 81 | 0 / 0 |
| S1 | channel=linear, delta=0.1, phi=0.8 | 0.000 | 0.000 | 1.000 | 0.095 [0.062, 0.144] | 0 / 0 / 19 | 0 / 0 |
| S1 | channel=linear, delta=0.1, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.055 [0.031, 0.096] | 0 / 0 / 11 | 0 / 0 |
| S1 | channel=linear, delta=0.25, phi=0.0 | 0.030 | 0.165 | 0.805 | 1.000 [0.981, 1.000] | 6 / 33 / 161 | 0 / 0 |
| S1 | channel=linear, delta=0.25, phi=0.8 | 0.010 | 0.005 | 0.985 | 0.570 [0.501, 0.637] | 2 / 1 / 111 | 0 / 0 |
| S1 | channel=linear, delta=0.25, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.190 [0.142, 0.250] | 0 / 0 / 38 | 0 / 0 |
| S1 | channel=linear, delta=0.5, phi=0.0 | 0.985 | 0.015 | 0.000 | 1.000 [0.981, 1.000] | 197 / 3 / 0 | 0 / 0 |
| S1 | channel=linear, delta=0.5, phi=0.8 | 0.045 | 0.260 | 0.695 | 0.995 [0.972, 0.999] | 9 / 52 / 138 | 0 / 0 |
| S1 | channel=linear, delta=0.5, phi=0.95 | 0.000 | 0.010 | 0.990 | 0.530 [0.461, 0.598] | 0 / 2 / 104 | 0 / 0 |
| S1 | channel=linear, delta=1.0, phi=0.0 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 200 / 0 / 0 | 0 / 0 |
| S1 | channel=linear, delta=1.0, phi=0.8 | 0.770 | 0.210 | 0.020 | 1.000 [0.981, 1.000] | 154 / 42 / 4 | 0 / 0 |
| S1 | channel=linear, delta=1.0, phi=0.95 | 0.075 | 0.245 | 0.680 | 0.990 [0.964, 0.997] | 15 / 49 / 134 | 0 / 0 |
| S1 | channel=linear, delta=2.0, phi=0.0 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 200 / 0 / 0 | 0 / 0 |
| S1 | channel=linear, delta=2.0, phi=0.8 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 200 / 0 / 0 | 0 / 0 |
| S1 | channel=linear, delta=2.0, phi=0.95 | 0.775 | 0.210 | 0.015 | 1.000 [0.981, 1.000] | 155 / 42 / 3 | 0 / 0 |
| S1 | channel=circular, delta=0.1, phi=0.0 | 0.000 | 0.000 | 1.000 | 0.355 [0.292, 0.423] | 0 / 0 / 71 | 0 / 0 |
| S1 | channel=circular, delta=0.1, phi=0.8 | 0.000 | 0.000 | 1.000 | 0.130 [0.090, 0.184] | 0 / 0 / 26 | 0 / 0 |
| S1 | channel=circular, delta=0.1, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.065 [0.038, 0.108] | 0 / 0 / 13 | 0 / 0 |
| S1 | channel=circular, delta=0.25, phi=0.0 | 0.030 | 0.140 | 0.830 | 1.000 [0.981, 1.000] | 6 / 28 / 166 | 0 / 0 |
| S1 | channel=circular, delta=0.25, phi=0.8 | 0.000 | 0.020 | 0.980 | 0.560 [0.491, 0.627] | 0 / 4 / 108 | 0 / 0 |
| S1 | channel=circular, delta=0.25, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.205 [0.155, 0.266] | 0 / 0 / 41 | 0 / 0 |
| S1 | channel=circular, delta=0.5, phi=0.0 | 0.960 | 0.035 | 0.005 | 1.000 [0.981, 1.000] | 192 / 7 / 1 | 0 / 0 |
| S1 | channel=circular, delta=0.5, phi=0.8 | 0.055 | 0.205 | 0.740 | 0.985 [0.957, 0.995] | 11 / 41 / 145 | 0 / 0 |
| S1 | channel=circular, delta=0.5, phi=0.95 | 0.005 | 0.010 | 0.985 | 0.510 [0.441, 0.578] | 1 / 2 / 99 | 0 / 0 |
| S1 | channel=circular, delta=1.0, phi=0.0 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 200 / 0 / 0 | 0 / 0 |
| S1 | channel=circular, delta=1.0, phi=0.8 | 0.855 | 0.125 | 0.020 | 1.000 [0.981, 1.000] | 171 / 25 / 4 | 0 / 0 |
| S1 | channel=circular, delta=1.0, phi=0.95 | 0.060 | 0.215 | 0.725 | 0.980 [0.950, 0.992] | 12 / 43 / 141 | 0 / 0 |
| S1 | channel=circular, delta=2.0, phi=0.0 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 200 / 0 / 0 | 0 / 0 |
| S1 | channel=circular, delta=2.0, phi=0.8 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 200 / 0 / 0 | 0 / 0 |
| S1 | channel=circular, delta=2.0, phi=0.95 | 0.760 | 0.215 | 0.025 | 1.000 [0.981, 1.000] | 152 / 43 / 5 | 0 / 0 |
| S2 | channel=linear, rho=1.25, phi=0.0 | 0.090 | 0.400 | 0.510 | 1.000 [0.981, 1.000] | 18 / 80 / 102 | 0 / 0 |
| S2 | channel=linear, rho=1.25, phi=0.8 | 0.000 | 0.020 | 0.980 | 0.880 [0.828, 0.918] | 0 / 4 / 172 | 0 / 0 |
| S2 | channel=linear, rho=1.25, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.245 [0.191, 0.309] | 0 / 0 / 49 | 0 / 0 |
| S2 | channel=linear, rho=1.5, phi=0.0 | 0.985 | 0.015 | 0.000 | 1.000 [0.981, 1.000] | 197 / 3 / 0 | 0 / 0 |
| S2 | channel=linear, rho=1.5, phi=0.8 | 0.120 | 0.375 | 0.505 | 1.000 [0.981, 1.000] | 24 / 75 / 101 | 0 / 0 |
| S2 | channel=linear, rho=1.5, phi=0.95 | 0.005 | 0.020 | 0.975 | 0.785 [0.723, 0.836] | 1 / 4 / 152 | 0 / 0 |
| S2 | channel=linear, rho=2.0, phi=0.0 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 200 / 0 / 0 | 0 / 0 |
| S2 | channel=linear, rho=2.0, phi=0.8 | 0.930 | 0.070 | 0.000 | 1.000 [0.981, 1.000] | 186 / 14 / 0 | 0 / 0 |
| S2 | channel=linear, rho=2.0, phi=0.95 | 0.020 | 0.160 | 0.820 | 1.000 [0.981, 1.000] | 4 / 32 / 164 | 0 / 0 |
| S3 | channel=linear, w_b=0.4, D=1 | 0.000 | 0.010 | 0.990 | 0.975 [0.943, 0.989] | 0 / 2 / 193 | 0 / 0 |
| S3 | channel=linear, w_b=0.4, D=20 | 0.000 | 0.000 | 1.000 | 0.090 [0.058, 0.138] | 0 / 0 / 18 | 0 / 0 |
| S3 | channel=linear, w_b=0.4, D=100 | 0.000 | 0.000 | 1.000 | 0.060 [0.035, 0.102] | 0 / 0 / 12 | 0 / 0 |
| S3 | channel=linear, w_b=0.3, D=1 | 0.565 | 0.350 | 0.085 | 1.000 [0.981, 1.000] | 113 / 70 / 17 | 0 / 0 |
| S3 | channel=linear, w_b=0.3, D=20 | 0.000 | 0.015 | 0.985 | 0.455 [0.387, 0.524] | 0 / 3 / 88 | 0 / 0 |
| S3 | channel=linear, w_b=0.3, D=100 | 0.000 | 0.010 | 0.990 | 0.160 [0.116, 0.217] | 0 / 2 / 30 | 0 / 0 |
| S3 | channel=linear, w_b=0.1, D=1 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 200 / 0 / 0 | 0 / 0 |
| S3 | channel=linear, w_b=0.1, D=20 | 0.105 | 0.265 | 0.630 | 0.965 [0.930, 0.983] | 21 / 53 / 119 | 0 / 0 |
| S3 | channel=linear, w_b=0.1, D=100 | 0.020 | 0.115 | 0.865 | 0.545 [0.476, 0.613] | 4 / 23 / 82 | 0 / 0 |
| S4a | channel=linear, delta=1.0, phi=0.0 | 0.000 | 1.000 | 0.000 | 1.000 [0.981, 1.000] | 0 / 200 / 0 | 0 / 0 |
| S4a | channel=linear, delta=1.0, phi=0.8 | 0.000 | 0.850 | 0.150 | 1.000 [0.981, 1.000] | 0 / 170 / 30 | 0 / 0 |
| S4a | channel=linear, delta=1.0, phi=0.95 | 0.000 | 0.090 | 0.910 | 0.645 [0.577, 0.708] | 0 / 18 / 111 | 0 / 0 |
| S4b | channel=linear, delta=1.0, k=0.25, phi=0.0 | 0.180 | 0.400 | 0.420 | 1.000 [0.981, 1.000] | 36 / 80 / 84 | 0 / 0 |
| S4b | channel=linear, delta=1.0, k=0.25, phi=0.8 | 0.100 | 0.330 | 0.570 | 1.000 [0.981, 1.000] | 20 / 66 / 114 | 0 / 0 |
| S4b | channel=linear, delta=1.0, k=0.25, phi=0.95 | 0.020 | 0.125 | 0.855 | 0.895 [0.845, 0.930] | 4 / 25 / 150 | 0 / 0 |
| S4b | channel=linear, delta=1.0, k=0.5, phi=0.0 | 0.020 | 0.080 | 0.900 | 0.765 [0.702, 0.818] | 4 / 16 / 133 | 0 / 0 |
| S4b | channel=linear, delta=1.0, k=0.5, phi=0.8 | 0.000 | 0.060 | 0.940 | 0.715 [0.649, 0.773] | 0 / 12 / 131 | 0 / 0 |
| S4b | channel=linear, delta=1.0, k=0.5, phi=0.95 | 0.010 | 0.050 | 0.940 | 0.585 [0.516, 0.651] | 2 / 10 / 105 | 0 / 0 |
