# BQ-1: B2 single-query optimization experiment (harness only, development data)

Protocol `docs/B2_SINGLE_QUERY_OPTIMIZATION_PROTOCOL.md` (SHA-256 `7ffe9d3f8034492c…`, frozen 2026-10-12T19:14:42.762346+00:00), Amendment 1 (`…AMENDMENT_1.md`, `e259c0596dba524f…`), code hashed 2026-10-12T19:20:35.474146+00:00 before any real-data computation. EB-1 versions and outcomes are unchanged. No production code was changed.

CPU used: 0.0729 h of the 1 CPU-h cap (tests, references, gates, calibration, preparation and timing). Stopped early: no.

## Outcome: **R2_B2 not supported by the tested variants on these workloads**

Strict (cold first query, preparation included) qualifying arms: none. Amendment 2 R2_B2 for reference: 25.5 s.

| Arm | status | gate | (a) gate + identity | (b) query ≤ 10 s | (c) prep models | (d) memory ≤ 8 GB | (e) R1 + prep ≤ 96 h | R2 (prepared DB) | R2 (cold) |
|---|---|---|---|---|---|---|---|---|---|
| V1 | timed | True | True | False | False | True | False | **False** | False |
| V2 | rejected_at_gate | False | — | — | — | — | — | **—** | — |
| V2b | timed | True | True | True | False | True | False | **False** | False |

Failed criteria per timed arm (frozen rules; nothing relaxed):

- V1: failed b_query_supported_and_le_10s, c_prep_models_supported, e_R1_B2_plus_prep_le_96h; prep_wall unsupported (fit residual 0.36 > 0.25; worst fit cell n_D·n_win = 5000, measured 0.00296 s, residual 0.36; max fresh confirmation error 0.066); prep_cpu unsupported (fit residual 0.34 > 0.25; worst fit cell n_D·n_win = 5000, measured 0.00296 s, residual 0.34; max fresh confirmation error 0.065); projected prepared-query latency 21.2 s > 10 s.
- V2: rejected_at_gate.
- V2b: failed c_prep_models_supported, e_R1_B2_plus_prep_le_96h; prep_wall unsupported (fit residual 0.28 > 0.25; worst fit cell n_D·n_win = 5000, measured 0.00307 s, residual 0.28; max fresh confirmation error 0.054); prep_cpu unsupported (fit residual 0.28 > 0.25; worst fit cell n_D·n_win = 5000, measured 0.00307 s, residual 0.28; max fresh confirmation error 0.054).

## 1. Checks before timing

- Regression tests (in-run): 17 passed in 1.91s (pass True).
- Reference R_svd vs explicit-rotation Kabsch: max 1.39e-13 Å; vs MDAnalysis qcprot: max 1.64e-12 Å (2,000 sampled pairs per workload). Pass True.
- V0 (harness) gate: True.

## 2. Correctness gates per arm (all workloads, 3 queries each, every window)

| Arm | workload | role | max |Δd| Å | non-equivalent | strict identical / query-windows | near-tie | fallback rows (fraction) |
|---|---|---|---|---|---|---|---|
| V0 | G_nD250_nwin20 | fit_reused | 2.96e-13 | 0 | 60/60 | 0 | 0 (0) |
| V0 | G_nD250_nwin207 | fit_reused | 2.96e-13 | 0 | 621/621 | 1 | 0 (0) |
| V0 | G_nD500_nwin20 | fit_reused | 2.96e-13 | 0 | 60/60 | 0 | 0 (0) |
| V0 | G_nD500_nwin207 | fit_reused | 3.67e-13 | 0 | 621/621 | 2 | 0 (0) |
| V0 | G_nD1015_nwin20 | fit_reused | 2.96e-13 | 0 | 60/60 | 0 | 0 (0) |
| V0 | G_nD1015_nwin207 | fit_reused | 4.6e-13 | 0 | 621/621 | 3 | 0 (0) |
| V0 | CW3 | comparison_reused | 4.6e-13 | 0 | 621/621 | 1 | 0 (0) |
| V0 | CW5 | comparison_reused | 1.47e-13 | 0 | 132/132 | 0 | 0 (0) |
| V0 | FQ1 | confirmation_fresh | 2.86e-13 | 0 | 360/360 | 1 | 0 (0) |
| V0 | FQ2 | confirmation_fresh | 4.02e-13 | 0 | 621/621 | 4 | 0 (0) |
| V0 | FQ3 | confirmation_fresh | 1.29e-13 | 0 | 132/132 | 0 | 0 (0) |
| V0 | FQ4 | confirmation_fresh | 3.77e-13 | 0 | 450/450 | 1 | 0 (0) |
| V1 | G_nD250_nwin20 | fit_reused | 2.96e-13 | 0 | 60/60 | 0 | 0 (0) |
| V1 | G_nD250_nwin207 | fit_reused | 2.96e-13 | 0 | 621/621 | 1 | 0 (0) |
| V1 | G_nD500_nwin20 | fit_reused | 2.96e-13 | 0 | 60/60 | 0 | 0 (0) |
| V1 | G_nD500_nwin207 | fit_reused | 3.67e-13 | 0 | 621/621 | 2 | 0 (0) |
| V1 | G_nD1015_nwin20 | fit_reused | 2.96e-13 | 0 | 60/60 | 0 | 0 (0) |
| V1 | G_nD1015_nwin207 | fit_reused | 4.6e-13 | 0 | 621/621 | 3 | 0 (0) |
| V1 | CW3 | comparison_reused | 4.6e-13 | 0 | 621/621 | 1 | 0 (0) |
| V1 | CW5 | comparison_reused | 1.47e-13 | 0 | 132/132 | 0 | 0 (0) |
| V1 | FQ1 | confirmation_fresh | 2.86e-13 | 0 | 360/360 | 1 | 0 (0) |
| V1 | FQ2 | confirmation_fresh | 4.02e-13 | 0 | 621/621 | 4 | 0 (0) |
| V1 | FQ3 | confirmation_fresh | 1.29e-13 | 0 | 132/132 | 0 | 0 (0) |
| V1 | FQ4 | confirmation_fresh | 3.77e-13 | 0 | 450/450 | 1 | 0 (0) |
| V2 | G_nD250_nwin20 | fit_reused | 4.89e-13 | 0 | 60/60 | 0 | 0 (0) |
| V2 | G_nD250_nwin207 | fit_reused | 8.63e-13 | 0 | 621/621 | 1 | 0 (0) |
| V2 | G_nD500_nwin20 | fit_reused | 4.89e-13 | 0 | 60/60 | 0 | 0 (0) |
| V2 | G_nD500_nwin207 | fit_reused | 9.24e-13 | 0 | 621/621 | 2 | 0 (0) |
| V2 | G_nD1015_nwin20 | fit_reused | 4.89e-13 | 0 | 60/60 | 0 | 0 (0) |
| V2 | G_nD1015_nwin207 | fit_reused | 1.14e-12 | 0 | 621/621 | 3 | 0 (0) |
| V2 | CW3 | comparison_reused | 1.14e-12 | 0 | 621/621 | 1 | 0 (0) |
| V2 | CW5 | comparison_reused | 4.35e-14 | 0 | 132/132 | 0 | 0 (0) |
| V2 | FQ1 | confirmation_fresh | 3.07e-12 | 0 | 360/360 | 1 | 0 (0) |
| V2 | FQ2 | confirmation_fresh | 5.22e-12 | 0 | 621/621 | 4 | 0 (0) |
| V2 | FQ3 | confirmation_fresh | 1.19e-13 | 0 | 132/132 | 0 | 0 (0) |
| V2 | FQ4 | confirmation_fresh | 1.32e-12 | 0 | 450/450 | 1 | 0 (0) |
| V2b | G_nD250_nwin20 | fit_reused | 4.89e-13 | 0 | 60/60 | 0 | 0 (0) |
| V2b | G_nD250_nwin207 | fit_reused | 8.63e-13 | 0 | 621/621 | 1 | 0 (0) |
| V2b | G_nD500_nwin20 | fit_reused | 4.89e-13 | 0 | 60/60 | 0 | 0 (0) |
| V2b | G_nD500_nwin207 | fit_reused | 9.24e-13 | 0 | 621/621 | 2 | 0 (0) |
| V2b | G_nD1015_nwin20 | fit_reused | 4.89e-13 | 0 | 60/60 | 0 | 0 (0) |
| V2b | G_nD1015_nwin207 | fit_reused | 1.14e-12 | 0 | 621/621 | 3 | 0 (0) |
| V2b | CW3 | comparison_reused | 1.14e-12 | 0 | 621/621 | 1 | 0 (0) |
| V2b | CW5 | comparison_reused | 4.35e-14 | 0 | 132/132 | 0 | 0 (0) |
| V2b | FQ1 | confirmation_fresh | 3.07e-12 | 0 | 360/360 | 1 | 0 (0) |
| V2b | FQ2 | confirmation_fresh | 5.22e-12 | 0 | 621/621 | 4 | 0 (0) |
| V2b | FQ3 | confirmation_fresh | 1.19e-13 | 0 | 132/132 | 0 | 0 (0) |
| V2b | FQ4 | confirmation_fresh | 1.32e-12 | 0 | 450/450 | 1 | 0 (0) |

SYNTHETIC adversarial suite (protocol §5.3):

- V0: pass True; degenerate_collinear_planar_mirror True (max 8.43e-08), fallback rows 0; identical_and_exact_duplicates True; permutation True (max 0); rigid_transform True (max 1.53e-14); near_ties True
- V1: pass True; degenerate_collinear_planar_mirror True (max 8.43e-08), fallback rows 0; identical_and_exact_duplicates True; permutation True (max 0); rigid_transform True (max 1.53e-14); near_ties True
- V2: pass False; degenerate_collinear_planar_mirror False (max 0.00031), fallback rows 46; identical_and_exact_duplicates True; permutation True (max 0); rigid_transform True (max 4.75e-14); near_ties True
- V2b: pass True; degenerate_collinear_planar_mirror True (max 8.43e-08), fallback rows 130; identical_and_exact_duplicates True; permutation True (max 0); rigid_transform True (max 4.75e-14); near_ties True

## 3. Direct measurements (same-session V0; medians over rounds; speed-up = per-round paired ratio, median (min–max))

### V1 (timed outputs identical to gate outputs: True)

| Workload | role | n_D × n_win | V0 query s | arm query s | prep s | first query s | speed-up prepared | speed-up cold | prepared MB | peak RSS MB arm / V0 | fallback rows/run |
|---|---|---|---|---|---|---|---|---|---|---|---|
| G_nD250_nwin20 | fit_reused | 250 × 20 | 0.0169 | 0.0134 | 0.00296 | 0.0163 | 1.24 (1.24–1.27) | 1.02 (1.02–1.08) | 1 | 42.5 / 41.6 | 0 |
| G_nD250_nwin207 | fit_reused | 250 × 207 | 0.174 | 0.137 | 0.025 | 0.164 | 1.27 (1.25–1.29) | 1.05 (1.04–1.11) | 10.3 | 51.5 / 41.5 | 0 |
| G_nD500_nwin20 | fit_reused | 500 × 20 | 0.03 | 0.0239 | 0.00542 | 0.0291 | 1.25 (1.22–1.26) | 1.01 (0.993–1.09) | 2 | 44.7 / 42.7 | 0 |
| G_nD500_nwin207 | fit_reused | 500 × 207 | 0.307 | 0.251 | 0.0458 | 0.294 | 1.22 (1.22–1.23) | 1.05 (1.03–1.05) | 20.7 | 62.5 / 42.7 | 0 |
| G_nD1015_nwin20 | fit_reused | 1015 × 20 | 0.0557 | 0.0462 | 0.011 | 0.057 | 1.2 (1.2–1.23) | 0.98 (0.967–0.992) | 4.06 | 49.2 / 45.2 | 0 |
| G_nD1015_nwin207 | fit_reused | 1015 × 207 | 0.574 | 0.474 | 0.0897 | 0.564 | 1.21 (1.21–1.22) | 1.02 (1.02–1.03) | 42 | 85.5 / 45.2 | 0 |
| CW3 | comparison_reused | 2693 × 207 | 1.55 | 1.2 | 0.298 | 1.5 | 1.3 (1.27–1.3) | 1.04 (1.03–1.04) | 111 | 161 / 54.4 | 0 |
| CW5 | comparison_reused | 405 × 44 | 0.0543 | 0.0446 | 0.0098 | 0.0544 | 1.22 (1.21–1.24) | 0.998 (0.989–1.02) | 3.56 | 44.1 / 40.7 | 0 |
| FQ1 | confirmation_fresh | 750 × 120 | 0.251 | 0.205 | 0.038 | 0.241 | 1.23 (1.22–1.23) | 1.04 (1.01–1.04) | 18 | 61.4 / 43.9 | 0 |
| FQ2 | confirmation_fresh | 2000 × 207 | 1.17 | 0.905 | 0.182 | 1.1 | 1.29 (1.29–1.31) | 1.06 (1.05–1.07) | 82.8 | 130 / 50.8 | 0 |
| FQ3 | confirmation_fresh | 405 × 44 | 0.0548 | 0.0443 | 0.00951 | 0.0535 | 1.23 (1.22–1.27) | 1.03 (1.02–1.06) | 3.56 | 44.1 / 40.7 | 0 |
| FQ4 | confirmation_fresh | 2693 × 150 | 1.13 | 0.868 | 0.184 | 1.06 | 1.3 (1.29–1.3) | 1.06 (1.02–1.08) | 80.8 | 131 / 54.4 | 0 |

Query CPU/wall ratio (one thread expected ≈ 1): 0.99–1.

### V2b (timed outputs identical to gate outputs: True)

| Workload | role | n_D × n_win | V0 query s | arm query s | prep s | first query s | speed-up prepared | speed-up cold | prepared MB | peak RSS MB arm / V0 | fallback rows/run |
|---|---|---|---|---|---|---|---|---|---|---|---|
| G_nD250_nwin20 | fit_reused | 250 × 20 | 0.017 | 0.0102 | 0.00307 | 0.0132 | 1.67 (1.57–1.68) | 1.28 (1.25–1.29) | 1 | 41.9 / 41.4 | 0 |
| G_nD250_nwin207 | fit_reused | 250 × 207 | 0.17 | 0.106 | 0.0245 | 0.131 | 1.6 (1.59–1.63) | 1.3 (1.28–1.31) | 10.3 | 50.8 / 41.5 | 0 |
| G_nD500_nwin20 | fit_reused | 500 × 20 | 0.0301 | 0.0138 | 0.00538 | 0.0193 | 2.17 (2.14–2.23) | 1.56 (1.54–1.61) | 2 | 44.1 / 42.7 | 0 |
| G_nD500_nwin207 | fit_reused | 500 × 207 | 0.302 | 0.14 | 0.0451 | 0.188 | 2.17 (2.11–2.19) | 1.64 (1.58–1.66) | 20.7 | 61.8 / 42.7 | 0 |
| G_nD1015_nwin20 | fit_reused | 1015 × 20 | 0.055 | 0.0217 | 0.0108 | 0.0328 | 2.53 (2.5–2.58) | 1.68 (1.66–1.71) | 4.06 | 48.7 / 45.4 | 0 |
| G_nD1015_nwin207 | fit_reused | 1015 × 207 | 0.57 | 0.227 | 0.0888 | 0.316 | 2.54 (2.52–2.6) | 1.83 (1.8–1.88) | 42 | 84.8 / 45.2 | 0 |
| CW3 | comparison_reused | 2693 × 207 | 1.55 | 0.483 | 0.28 | 0.763 | 3.21 (3.19–3.21) | 2.03 (2.02–2.09) | 111 | 160 / 54.4 | 0 |
| CW5 | comparison_reused | 405 × 44 | 0.0543 | 0.0287 | 0.00928 | 0.0383 | 1.89 (1.85–1.92) | 1.43 (1.4–1.44) | 3.56 | 43.4 / 40.7 | 0 |
| FQ1 | confirmation_fresh | 750 × 120 | 0.254 | 0.106 | 0.0391 | 0.145 | 2.4 (2.35–2.43) | 1.75 (1.71–1.78) | 18 | 60.7 / 43.9 | 0 |
| FQ2 | confirmation_fresh | 2000 × 207 | 1.16 | 0.376 | 0.179 | 0.552 | 3.07 (3.05–3.11) | 2.09 (1.98–2.13) | 82.8 | 129 / 50.8 | 0 |
| FQ3 | confirmation_fresh | 405 × 44 | 0.0544 | 0.0293 | 0.00974 | 0.0386 | 1.86 (1.86–1.86) | 1.4 (1.39–1.44) | 3.56 | 43.4 / 40.8 | 0 |
| FQ4 | confirmation_fresh | 2693 × 150 | 1.12 | 0.354 | 0.179 | 0.536 | 3.18 (3.15–3.22) | 2.09 (2.09–2.13) | 80.8 | 131 / 54.5 | 0 |

Query CPU/wall ratio (one thread expected ≈ 1): 0.998–1.

## 4. Models, confirmation and P2 projections (fitted; not measurements)

### V1

| Model | fit residual | supported | confirmation (fresh = deciding): workload rel. error (abs. error) | P2 value | extrapolation × |
|---|---|---|---|---|---|---|
| query_time | 0.013 | True | FQ1 0.014 (0.0029); FQ2 0.004 (0.0036); FQ3 0.0064 (-0.00028); FQ4 0.013 (0.012); CW3 (reused) 0.013 (0.015); CW5 (reused) 0.014 (-0.0006) | 21.2 | 17.8 |
| prep_wall | 0.36 | False | FQ1 0.046 (0.0017); FQ2 0.033 (-0.006); FQ3 0.012 (-0.00011); FQ4 0.066 (-0.012); CW3 (reused) 0.21 (-0.062); CW5 (reused) 0.04 (-0.0004) | 4.18 | 17.8 |
| prep_cpu | 0.34 | False | FQ1 0.036 (0.0014); FQ2 0.041 (-0.0074); FQ3 0.021 (-0.0002); FQ4 0.065 (-0.012); CW3 (reused) 0.22 (-0.064); CW5 (reused) 0.049 (-0.00048) | 4.14 | 17.8 |
| memory | 0.0027 | True | FQ1 0.0028 (-1.8e+05); FQ2 0.002 (-2.8e+05); FQ3 0.0023 (1.1e+05); FQ4 0.0025 (-3.5e+05); CW3 (reused) 0.0027 (-4.5e+05); CW5 (reused) 0.0022 (1e+05) | 2.29e+09 | 17.8 |
| V0_same_session_query_time | 0.015 | True | FQ1 0.0058 (0.0015); FQ2 0.072 (-0.084); FQ3 0.0049 (-0.00027); FQ4 0.067 (-0.076); CW3 (reused) 0.066 (-0.1); CW5 (reused) 0.0053 (0.00029) | 25.2 | 17.8 |
| V0_same_session_memory | 0.0023 | True | FQ1 0.0012 (5.4e+04); FQ2 0.016 (-8.7e+05); FQ3 0.0042 (1.8e+05); FQ4 0.02 (-1.1e+06); CW3 (reused) 0.02 (-1.1e+06); CW5 (reused) 0.0047 (2e+05) | 2.83e+08 | 17.6 |

P2: prepared-query 21.2 s; preparation — s wall, — s CPU; cold first query — s; peak RSS 2.29 GB; R1_B2 + preparation — CPU-h; same-session V0 25.2 s. Amortized per query at N = 1/10/100/1000: — s; break-even N for 10 s: —.

### V2b

| Model | fit residual | supported | confirmation (fresh = deciding): workload rel. error (abs. error) | P2 value | extrapolation × |
|---|---|---|---|---|---|---|
| query_time | 0.023 | True | FQ1 0.0047 (0.0005); FQ2 0.015 (0.0056); FQ3 0.066 (-0.0019); FQ4 0.0062 (0.0022); CW3 (reused) 0.017 (0.0083); CW5 (reused) 0.048 (-0.0014) | 7.92 | 17.8 |
| prep_wall | 0.28 | False | FQ1 0.0037 (0.00014); FQ2 0.027 (-0.0048); FQ3 0.049 (-0.00048); FQ4 0.054 (-0.0096); CW3 (reused) 0.17 (-0.046); CW5 (reused) 0.0019 (-1.7e-05) | 4.13 | 17.8 |
| prep_cpu | 0.28 | False | FQ1 0.0035 (0.00014); FQ2 0.027 (-0.0048); FQ3 0.05 (-0.00048); FQ4 0.054 (-0.0096); CW3 (reused) 0.17 (-0.046); CW5 (reused) 0.0018 (-1.7e-05) | 4.13 | 17.8 |
| memory | 0.0014 | True | FQ1 0.00087 (-5.6e+04); FQ2 0.0033 (-4.5e+05); FQ3 0.0034 (1.5e+05); FQ4 0.0052 (-7.2e+05); CW3 (reused) 0.0038 (-6.5e+05); CW5 (reused) 0.0032 (1.5e+05) | 2.28e+09 | 17.8 |
| V0_same_session_query_time | 0.016 | True | FQ1 0.014 (-0.0036); FQ2 0.064 (-0.074); FQ3 0.0083 (-0.00045); FQ4 0.067 (-0.075); CW3 (reused) 0.069 (-0.11); CW5 (reused) 0.0059 (-0.00032) | 25.2 | 17.8 |
| V0_same_session_memory | 0.0021 | True | FQ1 0.0012 (5.7e+04); FQ2 0.0095 (-5.1e+05); FQ3 0.00088 (3.8e+04); FQ4 0.012 (-6.9e+05); CW3 (reused) 0.011 (-6.4e+05); CW5 (reused) 0.0026 (1.1e+05) | 2.95e+08 | 17.6 |

P2: prepared-query 7.92 s; preparation — s wall, — s CPU; cold first query — s; peak RSS 2.28 GB; R1_B2 + preparation — CPU-h; same-session V0 25.2 s. Amortized per query at N = 1/10/100/1000: — s; break-even N for 10 s: —.

Units: time = wall s; memory = absolute peak RSS bytes. Time models: c0 + c1·(n_D·n_win) + c2·n_win; preparation: p0 + p1·(n_D·n_win); memory: m0 + m1·P.

## 5. Limitations

- Shared desktop, one thread, page caches not dropped; timings are conditional on contention.
- P2 values are extrapolations (factors above), confirmed only on the development workloads.
- Arms are NumPy harness implementations. The result says nothing about compiled or other B2 implementations, and does not generalize to all B2 implementations.
- V2b was added before measurement (Amendment 1) after a SYNTHETIC test exposed frozen V2 on degenerate windows. Its fallback constant was fixed from an error estimate, not tuned.
- Preparation is treated as a database-side stage under the primary verdict. The strict cold verdict counts it in the query latency.

