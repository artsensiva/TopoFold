# DEV-ENS-1: descriptive per-window comparison (AdK open-start vs closed-start)

Generated from `DEV-ENS-1_results.json` (SHA-256 ca203baa104caa0b…). Development evidence only: one run per starting condition, so differences mix the starting structure with relaxation from it. No test, state assignment, population or biological interpretation. Histogram Jensen–Shannon distance (base 2) inspired by established per-feature comparison; not a reproduction of PENSA.

## Observation period and blocks

- Samples: 420 per run, 0–335200 ps every 800 ps; timestamps identical in both runs: True
- Blocks: raw sizes {'open': [140, 140, 140], 'closed': [140, 140, 140]}; equalized to 140 samples (dropped 0); time ranges (ps) [[0.0, 111200.0], [112000.0, 223200.0], [224000.0, 335200.0]]

## Primary configuration (36 bins)

| Channel | Windows | Excluded | Full-run JS distance min / median / max | Windows with all 9 between-run block distances above the largest within-run one | … with none |
|---|---|---|---|---|---|
| kappa | 212 | 2 | 0.020 / 0.138 / 0.767 | 7 | 60 |
| d_i_i2 | 212 | 2 | 0.026 / 0.143 / 0.772 | 6 | 59 |
| tau | 211 | 2 | 0.025 / 0.194 / 0.953 | 6 | 41 |
| d_i_i3 | 211 | 2 | 0.033 / 0.180 / 0.892 | 8 | 42 |

Block fractions are descriptive (blocks are not independent replicates). Full-run and block distances use different sample sizes and are not compared in magnitude.

Out-of-range observations (primary): kappa 0 under / 0 over of 178080, d_i_i2 0 under / 0 over of 178080, d_i_i3 0 under / 0 over of 177240

## Rank stability across time (Spearman between matched block-pair rankings, primary)

- kappa: blocks 0-1: 0.523, blocks 0-2: 0.516, blocks 1-2: 0.611
- d_i_i2: blocks 0-1: 0.468, blocks 0-2: 0.463, blocks 1-2: 0.533
- tau: blocks 0-1: 0.549, blocks 0-2: 0.495, blocks 1-2: 0.659
- d_i_i3: blocks 0-1: 0.580, blocks 0-2: 0.551, blocks 1-2: 0.621

Paired channels on the same windows (Spearman of full-run distances): kappa~d_i_i2 0.944, tau~d_i_i3 0.692

κ and the i→i+2 distance are near-deterministic functions of each other for near-constant bond lengths; the i→i+3 distance cannot distinguish τ from −τ. Larger distances or different rankings do not show that one representation is better.

## Largest full-run distances (primary, descriptive)

- kappa: 120-122 (0.767; 0.67), 134-136 (0.740; 0.67), 41-43 (0.737; 0.67), 143-145 (0.731; 0.67), 9-11 (0.685; 0.67), 124-126 (0.569; 0.67), 152-154 (0.567; 0.89), 122-124 (0.548; 0.78), 135-137 (0.515; 0.44), 119-121 (0.489; 1.00)
- d_i_i2: 120-122 (0.772; 0.67), 134-136 (0.750; 0.67), 41-43 (0.731; 0.67), 143-145 (0.729; 0.67), 9-11 (0.682; 0.67), 152-154 (0.549; 0.56), 124-126 (0.548; 0.67), 122-124 (0.544; 0.89), 119-121 (0.523; 1.00), 135-137 (0.511; 0.44)
- tau: 121-124 (0.953; 1.00), 134-137 (0.830; 0.67), 119-122 (0.821; 0.67), 10-13 (0.788; 0.67), 143-146 (0.773; 0.67), 11-14 (0.769; 0.67), 42-45 (0.754; 0.67), 40-43 (0.724; 0.67), 152-155 (0.669; 0.44), 41-44 (0.665; 0.67)
- d_i_i3: 121-124 (0.892; 1.00), 119-122 (0.825; 0.78), 134-137 (0.777; 0.67), 142-145 (0.770; 0.67), 42-45 (0.750; 0.67), 11-14 (0.745; 0.44), 153-156 (0.715; 1.00), 143-146 (0.697; 0.67), 9-12 (0.679; 0.78), 8-11 (0.673; 0.67)

Each entry: window residues (full-run JS distance; fraction of between-run block distances above the largest within-run block distance).

## Sensitivity (predeclared; not used for selection)

| Configuration | Channel | Spearman vs primary | Top-10 overlap | Full-run median |
|---|---|---|---|---|
| bins18 | kappa | 0.973 | 9/10 | 0.110 |
| bins18 | d_i_i2 | 0.966 | 10/10 | 0.120 |
| bins18 | tau | 0.958 | 9/10 | 0.142 |
| bins18 | d_i_i3 | 0.982 | 10/10 | 0.149 |
| bins72 | kappa | 0.977 | 10/10 | 0.166 |
| bins72 | d_i_i2 | 0.971 | 10/10 | 0.170 |
| bins72 | tau | 0.990 | 10/10 | 0.213 |
| bins72 | d_i_i3 | 0.989 | 10/10 | 0.199 |
| tau_origin_half_bin | kappa | 1.000 | 10/10 | 0.138 |
| tau_origin_half_bin | d_i_i2 | 1.000 | 10/10 | 0.143 |
| tau_origin_half_bin | tau | 0.972 | 10/10 | 0.194 |
| tau_origin_half_bin | d_i_i3 | 1.000 | 10/10 | 0.180 |

## Resources

- Total CPU 28.7 s (cap 1800 s); wall 26.3 s; max RSS per stage (KiB): open 117976.0, closed 117732.0, comparison 60428
