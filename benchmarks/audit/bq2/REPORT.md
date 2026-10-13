# BQ-2: scaling confirmation of unchanged V2b (harness; development data + SYNTHETIC_tiled performance-only rows)

Protocol `docs/B2_SCALING_CONFIRMATION_PROTOCOL.md` (SHA-256 `6e0915b106a1f426…`, frozen 2026-10-12T19:35:58.995180+00:00); code and analysis hashed 2026-10-12T19:38:48.131002+00:00 before measurement. Informed by BQ-1 (not an independent rediscovery). BQ-1 files and its unsuccessful overall acceptance outcome are unchanged. No production code changed.

CPU: 0.182 h of the 1 CPU-h cap. Stopped early: no. Load average start/end: [['0.80', '0.92', '1.03'], ['0.89', '1.11', '1.36']].

## Conclusion

- **A. Unchanged V2b supports prepared-database operation:** False
- **B. First-query latency still warrants optimization:** False
- **C. Evidence remains insufficient:** True
- Projected violation with supported models: False; measured correctness failures: none; measured resource-limit violations: none.
- Missing evidence: prep_wall model unsupported (fit residual 0.225; max confirmation error 0.411); prep_cpu model unsupported (fit residual 0.209; max confirmation error 0.411).

Scope: a faster implementation of the existing retrieval objective on development and SYNTHETIC performance-only workloads. No biological capability, novelty over existing RMSD software or clinical relevance is implied.

## 1. Correctness gates (before timing)

Tests in-run: 5 passed in 0.22s. Timed outputs identical to gate outputs: True.

| Arm | workload | status | max |Δd| Å | non-equivalent | strict / query-windows | near-tie | fallback rows |
|---|---|---|---|---|---|---|---|
| V2b | F_nD100_nwin5 | ok | 7.37e-13 | 0 | 15/15 | 0 | 0 |
| V2b | F_nD100_nwin20 | ok | 5.77e-13 | 0 | 60/60 | 1 | 0 |
| V2b | F_nD100_nwin60 | ok | 5.47e-13 | 0 | 180/180 | 1 | 0 |
| V2b | F_nD100_nwin207 | ok | 7.81e-13 | 0 | 621/621 | 3 | 0 |
| V2b | F_nD700_nwin5 | ok | 1.05e-12 | 0 | 15/15 | 0 | 0 |
| V2b | F_nD700_nwin20 | ok | 3.02e-12 | 0 | 60/60 | 0 | 0 |
| V2b | F_nD700_nwin60 | ok | 8.87e-13 | 0 | 180/180 | 0 | 0 |
| V2b | F_nD700_nwin207 | ok | 3.02e-12 | 0 | 621/621 | 1 | 0 |
| V2b | F_nD2693_nwin5 | ok | 1.05e-12 | 0 | 15/15 | 0 | 0 |
| V2b | F_nD2693_nwin20 | ok | 3.02e-12 | 0 | 60/60 | 0 | 0 |
| V2b | F_nD2693_nwin60 | ok | 1.35e-12 | 0 | 180/180 | 1 | 0 |
| V2b | F_nD2693_nwin207 | ok | 3.02e-12 | 0 | 621/621 | 2 | 0 |
| V2b | RC1 | ok | 6.84e-12 | 0 | 360/360 | 1 | 0 |
| V2b | RC2 | ok | 1.38e-12 | 0 | 621/621 | 3 | 0 |
| V2b | RC3 | ok | 1.44e-13 | 0 | 132/132 | 1 | 0 |
| V2b | S1 | ok | 6.33e-12 | 0 | 2979/2979 | 13 | 0 |
| V2b | S2 | ok | 6.33e-12 | 0 | 2979/2979 | 17 | 0 |
| V2b | S4_P2slice | ok | 4.62e-07 | 0 | 300/300 | 300 | 0 |
| V2b | S3 | ok | 5.84e-07 | 0 | 2979/2979 | 814 | 0 |
| V0 | RC1 | ok | 3.33e-13 | 0 | 360/360 | 1 | 0 |
| V0 | RC2 | ok | 4.87e-13 | 0 | 621/621 | 3 | 0 |
| V0 | RC3 | ok | 1.3e-13 | 0 | 132/132 | 1 | 0 |
| V0 | S1 | ok | 4.64e-13 | 0 | 2979/2979 | 13 | 0 |
| V0 | S2 | ok | 4.77e-13 | 0 | 2979/2979 | 17 | 0 |
| V0 | S4_P2slice | ok | 1.69e-07 | 0 | 300/300 | 300 | 0 |
| V0 | S3 | ok | 2.38e-07 | 0 | 2979/2979 | 814 | 0 |

## 2. Direct measurements (medians over 3 rounds; speed-up = per-round paired ratio vs same-session V0)

| Workload | kind | role | n_D × n_win | prep s (reps rel. IQR) | V2b query s | first query s | V0 query s | speed-up prepared | speed-up cold | prepared MB | peak RSS MB | fallback rows/run | load avg |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| F_nD100_nwin5 | real | fit | 100 × 5 | 0.00042 (0.11) | 0.00202 | 0.00238 | — | — | — | 0.1 | 40.8 | 0 | 2.6 |
| F_nD100_nwin20 | real | fit | 100 × 20 | 0.00129 (0.11) | 0.00768 | 0.00885 | — | — | — | 0.4 | 41.3 | 0 | 2.3 |
| F_nD100_nwin60 | real | fit | 100 × 60 | 0.00352 (0.082) | 0.0222 | 0.0258 | — | — | — | 1.2 | 42.2 | 0 | 2.6 |
| F_nD100_nwin207 | real | fit | 100 × 207 | 0.0117 (0.089) | 0.0767 | 0.0891 | — | — | — | 4.14 | 44.8 | 0 | 2.3 |
| F_nD700_nwin5 | real | fit | 700 × 5 | 0.00171 (0.12) | 0.00433 | 0.00594 | — | — | — | 0.7 | 44.8 | 0 | 2.4 |
| F_nD700_nwin20 | real | fit | 700 × 20 | 0.00615 (0.032) | 0.0172 | 0.0233 | — | — | — | 2.8 | 46.6 | 0 | 2.5 |
| F_nD700_nwin60 | real | fit | 700 × 60 | 0.0178 (0.07) | 0.0504 | 0.0677 | — | — | — | 8.4 | 51.9 | 0 | 2.6 |
| F_nD700_nwin207 | real | fit | 700 × 207 | 0.0595 (0.058) | 0.176 | 0.232 | — | — | — | 29 | 71.6 | 0 | 2.5 |
| F_nD2693_nwin5 | real | fit | 2693 × 5 | 0.00596 (0.087) | 0.0122 | 0.0176 | — | — | — | 2.69 | 56.9 | 0 | 2.6 |
| F_nD2693_nwin20 | real | fit | 2693 × 20 | 0.0227 (0.1) | 0.047 | 0.0693 | — | — | — | 10.8 | 64.7 | 0 | 2.6 |
| F_nD2693_nwin60 | real | fit | 2693 × 60 | 0.0687 (0.059) | 0.139 | 0.208 | — | — | — | 32.3 | 85.1 | 0 | 2.3 |
| F_nD2693_nwin207 | real | fit | 2693 × 207 | 0.235 (0.017) | 0.472 | 0.71 | — | — | — | 111 | 161 | 0 | 2.6 |
| RC1 | real | confirm_real | 1500 × 120 | 0.0755 (0.027) | 0.175 | 0.25 | 0.514 | 2.94 (2.92–2.96) | 2.06 (2.03–2.08) | 36 | 82.4 | 0 | 2.3 |
| RC2 | real | confirm_real | 2693 × 207 | 0.235 (0.028) | 0.473 | 0.708 | 1.59 | 3.37 (3.33–3.38) | 2.26 (2.23–2.28) | 111 | 161 | 0 | 2.5 |
| RC3 | real | confirm_real | 405 × 44 | 0.00756 (0.072) | 0.0287 | 0.0365 | 0.0547 | 1.88 (1.87–1.93) | 1.52 (1.48–1.52) | 3.56 | 44.3 | 0 | 2.6 |
| S1 | SYNTHETIC_tiled | confirm_stress | 1000 × 993 | 0.407 (0.033) | 1.06 | 1.47 | 2.71 | 2.55 (2.5–2.56) | 1.84 (1.84–1.87) | 199 | 269 | 0 | 2.6 |
| S2 | SYNTHETIC_tiled | confirm_stress | 2500 × 993 | 1.36 (0.095) | 2.17 | 3.54 | 6.48 | 2.99 (2.98–3) | 1.83 (1.75–1.89) | 496 | 576 | 0 | 2.4 |
| S4_P2slice | SYNTHETIC_tiled | confirm_stress | 10000 × 100 | 0.712 (0.08) | 0.803 | 1.52 | 2.58 | 3.21 (3.17–3.29) | 1.73 (1.7–1.75) | 200 | 476 | 0 | 2.5 |
| S3 | SYNTHETIC_tiled | confirm_stress_if_memory | 4000 × 993 | 2.4 (0.027) | 3.28 | 5.69 | 10.3 | 3.14 (3.13–3.15) | 1.81 (1.79–1.82) | 794 | 913 | 0 | 2.4 |
| P2D | SYNTHETIC_tiled | direct_P2_if_memory | not measured | | | | | | | | | | |

**Direct P2 check (P2D, SYNTHETIC_tiled n_D 10,000 × 993):** status skipped_memory_rule; predicted peak 2.6 GB vs MemAvailable 1.73 GB minus 0.6 GB headroom.

## 3. Models (fitted on the real fit grid only) and confirmation

| Model | coefficients | fit residual | confirmation rel. error (abs.) | supported | P2 value | extrapolation × |
|---|---|---|---|---|---|---|
| query_time | c0=0.000592, c1=7.34e-07, c2=0.000309 | 0.24 | RC1 0.031 (-0.0054); RC2 0.002 (0.00097); RC3 0.049 (-0.0014); S1 0.026 (-0.028); S2 0.016 (-0.034); S4_P2slice 0.047 (-0.037); S3 0.017 (-0.055) | **True** | 7.6 | 2.5 |
| prep_wall | c0=0.000282, c1=4.19e-07, c2=4.58e-06 | 0.22 | RC1 0.0094 (0.00071); RC2 0.0017 (-0.00039); RC3 0.051 (0.00039); S1 0.034 (0.014); S2 0.23 (-0.32); S4_P2slice 0.41 (-0.29); S3 0.3 (-0.73) | **False** | 4.16 | 2.5 |
| prep_cpu | c0=0.000274, c1=4.19e-07, c2=4.65e-06 | 0.21 | RC1 0.0094 (0.00071); RC2 0.0013 (-0.0003); RC3 0.05 (0.00038); S1 0.034 (0.014); S2 0.23 (-0.32); S4_P2slice 0.41 (-0.29); S3 0.3 (-0.73) | **False** | 4.16 | 2.5 |
| memory | m0=4.23e+07, m1=0.998 | 0.0041 | RC1 0.001 (8.9e+04); RC2 0.0009 (-1.5e+05); RC3 0.0019 (8.8e+04); S1 0.06 (-1.7e+07); S2 0.0084 (-5.1e+06); S4_P2slice 0.029 (-1.4e+07); S3 0.026 (-2.5e+07) | **True** | 2.29e+09 | 2.53 |

Units: time models in s (c1 per row-window, c2 per window); memory in bytes. BQ-1 cells under the BQ-2 models (diagnostic only; BQ-1 verdict unchanged):

- prep_wall: G_nD250_nwin20 0.2 (-0.00061 s); G_nD250_nwin207 0.066 (-0.0016 s); G_nD500_nwin20 0.15 (-0.00082 s); G_nD500_nwin207 0.012 (-0.00052 s); G_nD1015_nwin20 0.18 (-0.0019 s); G_nD1015_nwin207 0.0051 (0.00045 s); CW3 0.16 (-0.045 s); CW5 0.14 (-0.0013 s); FQ1 0.015 (-0.00059 s); FQ2 0.023 (-0.0041 s); FQ3 0.18 (-0.0018 s); FQ4 0.051 (-0.0092 s).
- query_time: G_nD250_nwin20 0.025 (0.00026 s); G_nD250_nwin207 0.034 (-0.0036 s); G_nD500_nwin20 0.021 (0.00029 s); G_nD500_nwin207 0.0013 (0.00019 s); G_nD1015_nwin20 0.0029 (-6.3e-05 s); G_nD1015_nwin207 0.034 (-0.0076 s); CW3 0.019 (-0.0094 s); CW5 0.05 (-0.0014 s); FQ1 0.019 (-0.002 s); FQ2 0.02 (-0.0074 s); FQ3 0.068 (-0.002 s); FQ4 0.03 (-0.01 s).

## 4. P2 projections (fitted; valid only where the model is supported)

- Prepared-query latency: 7.6 s (slice-scaled from S4_P2slice, assuming linearity in windows only: 7.98 s); limit 10 s.
- Preparation: 4.16 s wall, 4.16 s CPU **(UNSUPPORTED — not a valid projection)**. First query (prep + query): 11.8 s **(UNSUPPORTED — not a valid projection)**.
- Amortized per query at N = 1/10/100/1000: 11.8, 8.02, 7.64, 7.6 s; break-even N for 10 s: 1.73 **(UNSUPPORTED — not a valid projection)**.
- Peak RSS: 2.29 GB (limit 8 GiB). R1_B2 + preparation: 5.79 CPU-h **(UNSUPPORTED — not a valid projection)** (limit 96; not evaluated if unsupported).
- Models supported: {'query_time': True, 'prep_wall': False, 'prep_cpu': False, 'memory': True}.

## 5. Post-hoc observations (not decision-bearing; computed from the direct measurements above)

Preparation wall time per database row-window (s), by workload and row length L:

| Workload | kind | L | n_D × n_win | prep s / (n_D·n_win) |
|---|---|---|---|---|
| F_nD100_nwin5 | real | 214 | 100 × 5 | 8.4e-07 |
| F_nD100_nwin20 | real | 214 | 100 × 20 | 6.46e-07 |
| F_nD100_nwin60 | real | 214 | 100 × 60 | 5.87e-07 |
| F_nD100_nwin207 | real | 214 | 100 × 207 | 5.65e-07 |
| F_nD700_nwin5 | real | 214 | 700 × 5 | 4.87e-07 |
| F_nD700_nwin20 | real | 214 | 700 × 20 | 4.39e-07 |
| F_nD700_nwin60 | real | 214 | 700 × 60 | 4.23e-07 |
| F_nD700_nwin207 | real | 214 | 700 × 207 | 4.11e-07 |
| F_nD2693_nwin5 | real | 214 | 2693 × 5 | 4.42e-07 |
| F_nD2693_nwin20 | real | 214 | 2693 × 20 | 4.22e-07 |
| F_nD2693_nwin60 | real | 214 | 2693 × 60 | 4.25e-07 |
| F_nD2693_nwin207 | real | 214 | 2693 × 207 | 4.22e-07 |
| RC1 | real | 214 | 1500 × 120 | 4.2e-07 |
| RC2 | real | 214 | 2693 × 207 | 4.22e-07 |
| RC3 | real | 51 | 405 × 44 | 4.24e-07 |
| S1 | SYNTHETIC_tiled | 1000 | 1000 × 993 | 4.1e-07 |
| S2 | SYNTHETIC_tiled | 1000 | 2500 × 993 | 5.49e-07 |
| S4_P2slice | SYNTHETIC_tiled | 1000 | 10000 × 100 | 7.12e-07 |
| S3 | SYNTHETIC_tiled | 1000 | 4000 × 993 | 6.03e-07 |

- The per-row-window preparation cost is roughly constant on the real workloads and on S1. It rises on S2, S3 and S4_P2slice, whose resident coordinate arrays (n_D × L) are the largest measured, and the real-data model under-predicts exactly those cells. A plausible structural cause, not tested here, is the strided gather `X[di, s:s+W]` over a resident array larger than the CPU caches. The frozen model has no working-set term, and at P2 the resident array is larger still.
- Slice-scaled preparation at P2 (S4_P2slice × 993/100, linearity in windows assumed; post hoc): 7.07 s; with the slice-scaled query, first query ≈ 15 s. Indicative only.
- The large near-tie counts on S3/S4_P2slice and the ~5e-7 Å maximum distance differences (both arms, within the 1e-6 Å tolerance) come from repeated tiled rows: identical windows give RMSD ≈ 0, where the square root amplifies rounding. Rankings remained equivalent in every query-window.

## 6. Limitations

- SYNTHETIC_tiled rows are performance-only: repeated real frames joined end to end (28 of 993 windows per row are artificial junctions). They test scaling with n_D·n_win, n_win up to 993 and n_D = 10,000. They do not establish performance on real 1000-residue proteins or real 10,000-frame trajectories, nor any retrieval or biological property.
- One thread on a shared desktop; page caches not dropped; load average recorded per cell.
- Unless P2D was measured, P2 values remain projections; the extrapolation factors are listed in §3.

