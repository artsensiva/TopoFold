# EB-1 Amendment 1: corrected accounting and targeted confirmation

Amendment `docs/EB1_AMENDMENT_1.md` (SHA-256 `459aa72c3beaef1a…`, frozen 2026-10-12T17:08:32.557112+00:00, after the original results and before any amended measurement). Original EB-1 protocol, measurements, analysis, report and outcome are unchanged. Development data only.

Cells: 77 of 80 completed. CPU: amendment 0.145 h; aggregate including the original run 0.217 h (cap 2 h).

## Decision under the amended rules: **insufficient_evidence** — missing cells: ["S3 {'d': 'D2', 'n': 1620}: failed", "S3 {'d': 'D2', 'n': 1620}: failed", "S3 {'d': 'D2', 'n': 1620}: failed"]

## 1. Correctness (separate subprocess)

- C1′ closed_start reconstruction vs WP-B derived: max |Δ| 0 Å.
- C2′ features: κ 6.11e-16, τ 4.44e-16 rad.
- C3′ D1+ search: A: 60/60 top-10 agree, near ties 0, max |Δd| 0; B2: 60/60 top-10 agree, near ties 0, max |Δd| 1.74e-13.
- C4′ instrumented search identical to the harness on CW4: True.
- Pass: **True**.

## 2. S4 models and confirmation (W = 8; per method and mode)

| Model | coefficients | fit residual | fit ok | confirmation: workload measured → predicted (rel. err.) | supported | P2 value | P2 extrapolation × |
|---|---|---|---|---|---|---|---|
| S4_A_batch_time | c0=0.0271, c1=1.83e-07, c2=2.96e-06 | 0.33 | False | CW1 21.9 → 22.2 (0.013); CW2 5.18 → 5.71 (0.1); CW4 1.14 → 0.967 (0.15) | **False** | 2.05e+03 | 92.3 |
| S4_A_batch_mem | m0=3.71e+07, m1=1.19 | 0.011 | True | CW1 8.22e+07 → 8.66e+07 (0.054); CW2 1.31e+08 → 1.38e+08 (0.051); CW4 5.29e+07 → 5.2e+07 (0.017) | **True** | 6.69e+08 | 6.26 |
| S4_A_single_time | c0=-7.66e-06, c1=1.92e-07, c2=9.66e-05 | 0.016 | True | CW3 0.128 → 0.127 (0.0045); CW5 0.00742 → 0.00766 (0.033) | **True** | 2 | 17.8 |
| S4_A_single_mem | m0=3.82e+07, m1=1.12 | 0.0045 | True | CW3 5.08e+07 → 4.93e+07 (0.029); CW5 3.85e+07 → 3.87e+07 (0.0046) | **True** | 2.38e+08 | 18 |
| S4_B2_batch_time | c0=0.0309, c1=2.1e-06, c2=2.16e-07 | 0.076 | True | CW1 237 → 226 (0.046); CW2 59.5 → 57.9 (0.026); CW4 9.62 → 9.61 (0.0011) | **True** | 2.08e+04 | 92.3 |
| S4_B2_batch_mem | m0=3.79e+07, m1=1.18 | 0.0076 | True | CW1 1.01e+08 → 1.08e+08 (0.072); CW2 1.91e+08 → 1.91e+08 (0.0012); CW4 5.93e+07 → 6.01e+07 (0.015) | **True** | 9.44e+08 | 5.9 |
| S4_B2_single_time | c0=0.000325, c1=2.55e-06, c2=0.000152 | 0.043 | True | CW3 1.44 → 1.45 (0.0083); CW5 0.0516 → 0.0524 (0.015) | **True** | 25.5 | 17.8 |
| S4_B2_single_mem | m0=3.8e+07, m1=1.61 | 0.0061 | True | CW3 6.45e+07 → 6.34e+07 (0.018); CW5 3.87e+07 → 3.93e+07 (0.015) | **True** | 4.74e+08 | 17.2 |

Units: batch time = CPU s; single time = wall s; memory = absolute peak RSS bytes.

## 3. Other stage fits (S1 reused from the original run with corrected memory predictor; S2/S3 re-measured)

| Fit | supported | c0 | c1 | max rel. residual | P2 projection |
|---|---|---|---|---|---|
| S1_xtc_wall | True | 0.881 | 5.89e-08 | 0.000907 | 104 |
| S1_h5_wall | True | 0.0441 | 9.5e-08 | 0.0062 | 4.22 |
| R4_S1_xtc | True | 1.15e+08 | 0.741 | 0.00256 | 2.13e+08 |
| R4_S1_h5 | True | 4.6e+07 | 3.9 | 0.0164 | 4.32e+09 |
| S2_D1_wall | True | -0.033 | 4.11e-06 | 0.0236 | 45.1 |
| S2_D2_wall | True | 0.00282 | 3.91e-06 | 0.021 | 43 |
| R4_S2 | True | 3.77e+07 | 0.975 | 0.0107 | 8.1e+08 |
| S3_D1_cpu | True | 0.00073 | 1.55e-06 | 0.00301 | 17 |
| S3_D2_cpu | False | — | — | — | — |
| S3_D1_wall_one_frame | True | 0.000715 | 1.55e-06 | 0.00301 | 0.00226 |
| S3_D2_wall_one_frame | False | — | — | — | — |
| R4_S3 | True | 3.77e+07 | 1.01 | 0.00666 | 6.12e+08 |

## 4. Requirements (method-specific composition; thresholds unchanged)

| Requirement | component | projected | limit | unit | ratio |
|---|---|---|---|---|---|
| R1_A | engine | — | 96 | CPU-h | — |
| R2_A | engine | — | 10 | s | — |
| R4_S4_A_batch | engine | 6.53e+05 | 8388608 | kB | 0.0778 |
| R4_S4_A_single | engine | 2.32e+05 | 8388608 | kB | 0.0277 |
| R1_B2 | engine | 5.79 | 96 | CPU-h | 0.0603 |
| R2_B2 | engine | 25.5 | 10 | s | 2.55 |
| R4_S4_B2_batch | engine | 9.21e+05 | 8388608 | kB | 0.11 |
| R4_S4_B2_single | engine | 4.63e+05 | 8388608 | kB | 0.0552 |
| R3_xtc | loader | 150 | 7.2e+03 | s | 0.0208 |
| R3_h5 | loader | 49.3 | 7.2e+03 | s | 0.00685 |
| R4_S1_xtc | loader | 2.08e+05 | 8388608 | kB | 0.0248 |
| R4_S1_h5 | loader | 4.22e+06 | 8388608 | kB | 0.503 |
| R4_S2 | loader | 7.91e+05 | 8388608 | kB | 0.0943 |
| R4_S3_A | engine | 5.98e+05 | 8388608 | kB | 0.0713 |

## 5. S4 time decomposition (share of total; incidental harness slice copies vs required kernel/top-k work)

| Cell | slice copies | kernel | top-k | other loop |
|---|---|---|---|---|
| A/batch/fit/nQ32/nD1015/nwin20 | 0.017 | 0.94 | 0.036 | 0.0019 |
| A/batch/fit/nQ64/nD1015/nwin20 | 0.0086 | 0.96 | 0.031 | 0.0011 |
| A/batch/fit/nQ128/nD1015/nwin20 | 0.0052 | 0.96 | 0.034 | 0.00067 |
| A/batch/fit/nQ192/nD1015/nwin20 | 0.0071 | 0.96 | 0.035 | 0.00057 |
| A/batch/fit/nQ256/nD1015/nwin20 | 0.0053 | 0.96 | 0.034 | 0.00047 |
| A/batch/fit/nQ64/nD1015/nwin5 | 0.0082 | 0.96 | 0.031 | 0.001 |
| A/batch/fit/nQ64/nD1015/nwin80 | 0.0086 | 0.96 | 0.031 | 0.0011 |
| A/single/fit/nQ1/nD1015/nwin20 | 0.2 | 0.64 | 0.13 | 0.027 |
| A/single/fit/nQ1/nD1015/nwin207 | 0.21 | 0.63 | 0.14 | 0.028 |
| A/batch/CW1/nQ512/nD1015/nwin207 | 0.0054 | 0.96 | 0.034 | 0.00032 |
| A/batch/CW2/nQ256/nD2693/nwin40 | 0.0047 | 0.97 | 0.027 | 0.00019 |
| A/single/CW3/nQ1/nD2693/nwin207 | 0.22 | 0.69 | 0.076 | 0.014 |
| A/batch/CW4/nQ256/nD405/nwin44 | 0.0046 | 0.95 | 0.046 | 0.00089 |
| A/single/CW5/nQ1/nD405/nwin44 | 0.14 | 0.61 | 0.21 | 0.046 |
| B2/batch/fit/nQ32/nD1015/nwin20 | 0.0014 | 0.99 | 0.0045 | 0.00029 |
| B2/batch/fit/nQ64/nD1015/nwin20 | 0.0006 | 1 | 0.0037 | 0.00014 |
| B2/batch/fit/nQ128/nD1015/nwin20 | 0.00031 | 1 | 0.0033 | 6.9e-05 |
| B2/batch/fit/nQ192/nD1015/nwin20 | 0.00042 | 1 | 0.0034 | 5.9e-05 |
| B2/batch/fit/nQ256/nD1015/nwin20 | 0.00032 | 1 | 0.0033 | 4.7e-05 |
| B2/batch/fit/nQ64/nD1015/nwin5 | 0.00062 | 1 | 0.0038 | 0.00014 |
| B2/batch/fit/nQ64/nD1015/nwin80 | 0.00062 | 1 | 0.0036 | 0.00015 |
| B2/single/fit/nQ1/nD1015/nwin20 | 0.032 | 0.95 | 0.018 | 0.0039 |
| B2/single/fit/nQ1/nD1015/nwin207 | 0.027 | 0.95 | 0.019 | 0.0038 |
| B2/batch/CW1/nQ512/nD1015/nwin207 | 0.00073 | 0.99 | 0.0084 | 0.00041 |
| B2/batch/CW2/nQ256/nD2693/nwin40 | 0.00029 | 1 | 0.0024 | 1.8e-05 |
| B2/single/CW3/nQ1/nD2693/nwin207 | 0.036 | 0.95 | 0.0089 | 0.0016 |
| B2/batch/CW4/nQ256/nD405/nwin44 | 0.00042 | 0.99 | 0.0056 | 0.00011 |
| B2/single/CW5/nQ1/nD405/nwin44 | 0.015 | 0.94 | 0.038 | 0.0081 |

## 6. What remains valid from EB-1, what was corrected

- **Valid:** original correctness checks C0–C3; S1 ingestion timings and peak RSS (process unaffected; predictor corrected only); the original outcome as a record of the original protocol.
- **Corrected:** S4 memory accounting (effective batch, n_D, method-specific transients, retained outputs, per method and mode); S2/S3 memory (owned inputs; no reference in the S3 process); method-specific paths (B2 without S3); confirmation workloads; time decomposition.

## 7. Limitations

- Shared desktop under memory pressure; page caches not dropped; timings conditional on contention.
- P2 values remain extrapolations beyond the largest confirmed workload (factors in §2–§3); confirmation tests prediction within the measured range plus moderate extrapolation, not at P2.
- Transient-byte coefficients come from a SYNTHETIC tracemalloc calibration at W = 8; RSS includes allocator retention that tracemalloc does not see.
- D1+ adds real closed_start frames for capacity only; it makes no retrieval-quality claim.
- Only existing implementations are measured; no saving is claimed for an implementation that does not exist.

