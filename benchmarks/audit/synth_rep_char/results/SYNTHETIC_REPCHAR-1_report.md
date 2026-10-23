# SYNTHETIC_REPCHAR-1: sensitivity of the replica-comparison rule (SYNTHETIC data)

Generated from `SYNTHETIC_REPCHAR-1_results.json`. Status: **complete**. 200 repetitions per cell, which gives coarse Monte Carlo precision (Wilson 95 % intervals shown). Categories: C = consistent, RD = replica-dependent, I = inconclusive; call = C + RD, always with its breakdown. No real data were used; DEV-REP-1 is not re-scored.

Resources: 289.1 CPU-s (cap 900), wall 290.0 s.

## S0 identical-distribution null

| Cell | C | RD | I | call [95 % CI] | median fraction (IQR) | median between/within |
|---|---|---|---|---|---|---|
| channel=linear, phi=0.0 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 1.00 |
| channel=linear, phi=0.8 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 0.99 |
| channel=linear, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 0.98 |
| channel=circular, phi=0.0 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 0.99 |
| channel=circular, phi=0.8 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 0.99 |
| channel=circular, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.04) | 0.98 |

## S5 condition-confounded sampling

| Cell | C | RD | I | call [95 % CI] | median fraction (IQR) | median between/within | median realized offset diff, calls / no calls (σ) |
|---|---|---|---|---|---|---|---|
| channel=linear, k=0.25, phi=0.0 | 0.295 | 0.110 | 0.595 | 0.405 [0.339, 0.474] | 0.50 (0.12–1.00) | 1.39 | 0.45880451559039387 / 0.12241583794470458 |
| channel=linear, k=0.25, phi=0.8 | 0.020 | 0.070 | 0.910 | 0.090 [0.058, 0.138] | 0.16 (0.04–0.44) | 1.18 | 0.5728542436782593 / 0.19996641572905274 |
| channel=linear, k=0.25, phi=0.95 | 0.000 | 0.015 | 0.985 | 0.015 [0.005, 0.043] | 0.08 (0.00–0.16) | 1.05 | 0.6690208190319372 / 0.21411108763556777 |
| channel=linear, k=0.5, phi=0.0 | 0.615 | 0.075 | 0.310 | 0.690 [0.623, 0.750] | 1.00 (0.66–1.00) | 2.29 | 0.686214896753206 / 0.14681337790320115 |
| channel=linear, k=0.5, phi=0.8 | 0.220 | 0.135 | 0.645 | 0.355 [0.292, 0.423] | 0.52 (0.12–0.96) | 1.52 | 0.8242381676674865 / 0.2749882374867254 |
| channel=linear, k=0.5, phi=0.95 | 0.040 | 0.065 | 0.895 | 0.105 [0.070, 0.155] | 0.16 (0.04–0.44) | 1.26 | 1.2382913502352673 / 0.4750229834756646 |

## S1 location shift (units of σ or s)

| Cell | C | RD | I | call [95 % CI] | median fraction (IQR) | median between/within |
|---|---|---|---|---|---|---|
| channel=linear, delta=0.1, phi=0.0 | 0.000 | 0.005 | 0.995 | 0.005 [0.001, 0.028] | 0.08 (0.04–0.16) | 1.08 |
| channel=linear, delta=0.1, phi=0.8 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 1.02 |
| channel=linear, delta=0.1, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 0.99 |
| channel=linear, delta=0.25, phi=0.0 | 0.030 | 0.165 | 0.805 | 0.195 [0.146, 0.255] | 0.64 (0.48–0.76) | 1.46 |
| channel=linear, delta=0.25, phi=0.8 | 0.010 | 0.005 | 0.985 | 0.015 [0.005, 0.043] | 0.16 (0.08–0.32) | 1.16 |
| channel=linear, delta=0.25, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.16) | 1.05 |
| channel=linear, delta=0.5, phi=0.0 | 0.985 | 0.015 | 0.000 | 1.000 [0.981, 1.000] | 1.00 (1.00–1.00) | 2.40 |
| channel=linear, delta=0.5, phi=0.8 | 0.045 | 0.260 | 0.695 | 0.305 [0.245, 0.372] | 0.60 (0.44–0.77) | 1.73 |
| channel=linear, delta=0.5, phi=0.95 | 0.000 | 0.010 | 0.990 | 0.010 [0.003, 0.036] | 0.16 (0.08–0.29) | 1.23 |
| channel=linear, delta=1.0, phi=0.0 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 1.00 (1.00–1.00) | 4.38 |
| channel=linear, delta=1.0, phi=0.8 | 0.770 | 0.210 | 0.020 | 0.980 [0.950, 0.992] | 1.00 (1.00–1.00) | 2.97 |
| channel=linear, delta=1.0, phi=0.95 | 0.075 | 0.245 | 0.680 | 0.320 [0.259, 0.388] | 0.60 (0.44–0.80) | 1.85 |
| channel=linear, delta=2.0, phi=0.0 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 1.00 (1.00–1.00) | 7.39 |
| channel=linear, delta=2.0, phi=0.8 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 1.00 (1.00–1.00) | 5.17 |
| channel=linear, delta=2.0, phi=0.95 | 0.775 | 0.210 | 0.015 | 0.985 [0.957, 0.995] | 1.00 (1.00–1.00) | 2.99 |
| channel=circular, delta=0.1, phi=0.0 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.08 (0.00–0.16) | 1.06 |
| channel=circular, delta=0.1, phi=0.8 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 1.03 |
| channel=circular, delta=0.1, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 1.00 |
| channel=circular, delta=0.25, phi=0.0 | 0.030 | 0.140 | 0.830 | 0.170 [0.124, 0.228] | 0.60 (0.44–0.76) | 1.44 |
| channel=circular, delta=0.25, phi=0.8 | 0.000 | 0.020 | 0.980 | 0.020 [0.008, 0.050] | 0.16 (0.08–0.32) | 1.16 |
| channel=circular, delta=0.25, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.12) | 1.04 |
| channel=circular, delta=0.5, phi=0.0 | 0.960 | 0.035 | 0.005 | 0.995 [0.972, 0.999] | 1.00 (1.00–1.00) | 2.23 |
| channel=circular, delta=0.5, phi=0.8 | 0.055 | 0.205 | 0.740 | 0.260 [0.204, 0.325] | 0.58 (0.40–0.77) | 1.64 |
| channel=circular, delta=0.5, phi=0.95 | 0.005 | 0.010 | 0.985 | 0.015 [0.005, 0.043] | 0.16 (0.08–0.28) | 1.18 |
| channel=circular, delta=1.0, phi=0.0 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 1.00 (1.00–1.00) | 3.99 |
| channel=circular, delta=1.0, phi=0.8 | 0.855 | 0.125 | 0.020 | 0.980 [0.950, 0.992] | 1.00 (1.00–1.00) | 2.89 |
| channel=circular, delta=1.0, phi=0.95 | 0.060 | 0.215 | 0.725 | 0.275 [0.218, 0.341] | 0.60 (0.40–0.77) | 1.81 |
| channel=circular, delta=2.0, phi=0.0 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 1.00 (1.00–1.00) | 6.88 |
| channel=circular, delta=2.0, phi=0.8 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 1.00 (1.00–1.00) | 4.96 |
| channel=circular, delta=2.0, phi=0.95 | 0.760 | 0.215 | 0.025 | 0.975 [0.943, 0.989] | 1.00 (1.00–1.00) | 2.89 |

## S2 scale change

| Cell | C | RD | I | call [95 % CI] | median fraction (IQR) | median between/within |
|---|---|---|---|---|---|---|
| channel=linear, rho=1.25, phi=0.0 | 0.090 | 0.400 | 0.510 | 0.490 [0.422, 0.559] | 0.76 (0.64–0.92) | 1.59 |
| channel=linear, rho=1.25, phi=0.8 | 0.000 | 0.020 | 0.980 | 0.020 [0.008, 0.050] | 0.24 (0.08–0.40) | 1.31 |
| channel=linear, rho=1.25, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.12) | 1.09 |
| channel=linear, rho=1.5, phi=0.0 | 0.985 | 0.015 | 0.000 | 1.000 [0.981, 1.000] | 1.00 (1.00–1.00) | 2.32 |
| channel=linear, rho=1.5, phi=0.8 | 0.120 | 0.375 | 0.505 | 0.495 [0.426, 0.564] | 0.76 (0.55–0.92) | 1.80 |
| channel=linear, rho=1.5, phi=0.95 | 0.005 | 0.020 | 0.975 | 0.025 [0.011, 0.057] | 0.16 (0.08–0.32) | 1.35 |
| channel=linear, rho=2.0, phi=0.0 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 1.00 (1.00–1.00) | 3.29 |
| channel=linear, rho=2.0, phi=0.8 | 0.930 | 0.070 | 0.000 | 1.000 [0.981, 1.000] | 1.00 (1.00–1.00) | 2.52 |
| channel=linear, rho=2.0, phi=0.95 | 0.020 | 0.160 | 0.820 | 0.180 [0.133, 0.239] | 0.48 (0.24–0.68) | 1.66 |

## S3 population shift between fixed modes

| Cell | C | RD | I | call [95 % CI] | median fraction (IQR) | median between/within |
|---|---|---|---|---|---|---|
| channel=linear, w_b=0.4, D=1 | 0.000 | 0.010 | 0.990 | 0.010 [0.003, 0.036] | 0.32 (0.16–0.45) | 1.21 |
| channel=linear, w_b=0.4, D=20 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 0.99 |
| channel=linear, w_b=0.4, D=100 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 0.98 |
| channel=linear, w_b=0.3, D=1 | 0.565 | 0.350 | 0.085 | 0.915 [0.868, 0.946] | 1.00 (0.92–1.00) | 1.74 |
| channel=linear, w_b=0.3, D=20 | 0.000 | 0.015 | 0.985 | 0.015 [0.005, 0.043] | 0.16 (0.08–0.28) | 1.24 |
| channel=linear, w_b=0.3, D=100 | 0.000 | 0.010 | 0.990 | 0.010 [0.003, 0.036] | 0.04 (0.00–0.12) | 1.05 |
| channel=linear, w_b=0.1, D=1 | 1.000 | 0.000 | 0.000 | 1.000 [0.981, 1.000] | 1.00 (1.00–1.00) | 3.31 |
| channel=linear, w_b=0.1, D=20 | 0.105 | 0.265 | 0.630 | 0.370 [0.306, 0.439] | 0.64 (0.40–0.81) | 2.35 |
| channel=linear, w_b=0.1, D=100 | 0.020 | 0.115 | 0.865 | 0.135 [0.094, 0.189] | 0.16 (0.08–0.32) | 1.63 |

## S4a one outlier replica

| Cell | C | RD | I | call [95 % CI] | median fraction (IQR) | median between/within |
|---|---|---|---|---|---|---|
| channel=linear, delta=0.0, phi=0.0 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 0.98 |
| channel=linear, delta=0.0, phi=0.8 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 1.00 |
| channel=linear, delta=0.0, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 1.00 |
| channel=linear, delta=1.0, phi=0.0 | 0.000 | 1.000 | 0.000 | 1.000 [0.981, 1.000] | 0.00 (0.00–0.00) | 4.04 |
| channel=linear, delta=1.0, phi=0.8 | 0.000 | 0.850 | 0.150 | 0.850 [0.794, 0.893] | 0.00 (0.00–0.00) | 2.70 |
| channel=linear, delta=1.0, phi=0.95 | 0.000 | 0.090 | 0.910 | 0.090 [0.058, 0.138] | 0.00 (0.00–0.00) | 1.64 |

## S4b random replica offsets

| Cell | C | RD | I | call [95 % CI] | median fraction (IQR) | median between/within |
|---|---|---|---|---|---|---|
| channel=linear, delta=0.0, k=0.25, phi=0.0 | 0.000 | 0.005 | 0.995 | 0.005 [0.001, 0.028] | 0.04 (0.00–0.08) | 0.96 |
| channel=linear, delta=0.0, k=0.25, phi=0.8 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.04) | 0.99 |
| channel=linear, delta=0.0, k=0.25, phi=0.95 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 0.98 |
| channel=linear, delta=1.0, k=0.25, phi=0.0 | 0.180 | 0.400 | 0.420 | 0.580 [0.511, 0.646] | 0.76 (0.56–0.92) | 2.80 |
| channel=linear, delta=1.0, k=0.25, phi=0.8 | 0.100 | 0.330 | 0.570 | 0.430 [0.363, 0.499] | 0.64 (0.48–0.88) | 2.37 |
| channel=linear, delta=1.0, k=0.25, phi=0.95 | 0.020 | 0.125 | 0.855 | 0.145 [0.103, 0.200] | 0.36 (0.24–0.64) | 1.70 |
| channel=linear, delta=0.0, k=0.5, phi=0.0 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 0.95 |
| channel=linear, delta=0.0, k=0.5, phi=0.8 | 0.000 | 0.005 | 0.995 | 0.005 [0.001, 0.028] | 0.04 (0.00–0.08) | 0.96 |
| channel=linear, delta=0.0, k=0.5, phi=0.95 | 0.000 | 0.005 | 0.995 | 0.005 [0.001, 0.028] | 0.04 (0.00–0.04) | 0.96 |
| channel=linear, delta=1.0, k=0.5, phi=0.0 | 0.020 | 0.080 | 0.900 | 0.100 [0.066, 0.149] | 0.28 (0.16–0.48) | 1.78 |
| channel=linear, delta=1.0, k=0.5, phi=0.8 | 0.000 | 0.060 | 0.940 | 0.060 [0.035, 0.102] | 0.24 (0.16–0.40) | 1.66 |
| channel=linear, delta=1.0, k=0.5, phi=0.95 | 0.010 | 0.050 | 0.940 | 0.060 [0.035, 0.102] | 0.20 (0.08–0.36) | 1.38 |

## S4c rare excursion in one replica

| Cell | C | RD | I | call [95 % CI] | median fraction (IQR) | median between/within |
|---|---|---|---|---|---|---|
| channel=linear, D=1 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 1.00 |
| channel=linear, D=20 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.04 (0.00–0.08) | 1.00 |
| channel=linear, D=100 | 0.000 | 0.000 | 1.000 | 0.000 [0.000, 0.019] | 0.00 (0.00–0.04) | 1.00 |

## S6 replica count (Δ = 1, φ = 0.8)

| Cell | C | RD | I | call [95 % CI] | median fraction (IQR) | median between/within |
|---|---|---|---|---|---|---|
| channel=linear, delta=1.0, n_rep=3, phi=0.8 | 0.930 | 0.070 | 0.000 | 1.000 [0.981, 1.000] | 1.00 (1.00–1.00) | 3.03 |
| channel=linear, delta=1.0, n_rep=5, phi=0.8 | 0.810 | 0.160 | 0.030 | 0.970 [0.936, 0.986] | 1.00 (1.00–1.00) | 3.02 |
| channel=linear, delta=1.0, n_rep=8, phi=0.8 | 0.510 | 0.405 | 0.085 | 0.915 [0.868, 0.946] | 1.00 (0.95–1.00) | 2.98 |

## Reading aid (not a criterion): smallest effect with call frequency ≥ 0.5

- S1 [channel=linear, phi=0.0]: delta=0.5 (C 197, RD 3 of 200)
- S1 [channel=linear, phi=0.8]: delta=1.0 (C 154, RD 42 of 200)
- S1 [channel=linear, phi=0.95]: delta=2.0 (C 155, RD 42 of 200)
- S1 [channel=circular, phi=0.0]: delta=0.5 (C 192, RD 7 of 200)
- S1 [channel=circular, phi=0.8]: delta=1.0 (C 171, RD 25 of 200)
- S1 [channel=circular, phi=0.95]: delta=2.0 (C 152, RD 43 of 200)
- S2 [channel=linear, phi=0.0]: rho=1.5 (C 197, RD 3 of 200)
- S2 [channel=linear, phi=0.8]: rho=2.0 (C 186, RD 14 of 200)
- S2 [channel=linear, phi=0.95]: not reached in the grid
- S3 [channel=linear, D=1]: w_b=0.3 (C 113, RD 70 of 200)
- S3 [channel=linear, D=20]: not reached in the grid
- S3 [channel=linear, D=100]: not reached in the grid
