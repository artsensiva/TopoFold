# EB-1 Amendment 2: remaining defects

Amendment `docs/EB1_AMENDMENT_2.md` (SHA-256 `ab206d3999509069…`, frozen 2026-10-12T18:57:35.278397+00:00, after the read-only preflight and before any Amendment 2 measurement). EB-1 and Amendment 1 files and outcomes are unchanged (Amendment 1: `insufficient_evidence`); Amendment 1 did not select engine work. Development data only.

Cells: 74 of 74 completed. CPU: Amendment 2 0.0325 h; aggregate including EB-1 and Amendment 1 0.249 h (cap 2 h).

## Decision under the frozen rule: **engine**

This is the Amendment 2 outcome. It is not a retrospective outcome of EB-1 or Amendment 1. Under the frozen rule the engine flag comes from the violated engine requirements listed in §1; requirements within limits are projections, not guarantees.

## 1. Validity by group

| Group | valid | status | violated | requirements (projected / limit, ratio) |
|---|---|---|---|---|
| loader | True | within_requirements | — | R3_xtc 150/7.2e+03 s (0.021); R3_h5 49.3/7.2e+03 s (0.0069); R4_S1_xtc 2.08e+05/8388608 kB (0.025); R4_S1_h5 4.22e+06/8388608 kB (0.5); R4_S2 8.07e+05/8388608 kB (0.096) |
| A | True | within_requirements | — | R1_A 0.526/96 CPU-h (0.0055); R2_A 2/10 s (0.2); R4_S4_A_batch 6.54e+05/8388608 kB (0.078); R4_S4_A_single 2.05e+05/8388608 kB (0.024); R4_S3_A 6.02e+05/8388608 kB (0.072) |
| B2 | True | requirements_violated | R2_B2 | R1_B2 5.79/96 CPU-h (0.06); R2_B2 25.5/10 s (2.5); R4_S4_B2_batch 9.25e+05/8388608 kB (0.11); R4_S4_B2_single 2.82e+05/8388608 kB (0.034) |

## 2. Preflight (read-only; not a measurement)

- D1 stored frames 2253 (no burn-in). D2 stored per trajectory [450, 450, 430, 450] = 1780 before burn-in; [405, 405, 387, 405] = 1602 after (array: 1602). D1+ 3706 rows.
- Invalid Amendment 1 requests: S3 {'d': 'D2', 'n': 1620} requested 1620 of 1602; S3 {'d': 'D2', 'n': 1620} requested 1620 of 1602; S3 {'d': 'D2', 'n': 1620} requested 1620 of 1602.
- Amendment 1 memory results affected: 70 of 80; replaced by re-measurement (70); unaffected reused 10.

## 3. Correctness

- New direct-read loaders vs Amendment 1 load paths: 41 of 41 distinct inputs identical (dtype, shape, values); owned: True. FC1 search through `measure_reps` identical to the EB-1 harness: True.
- Every measured cell's inputs owned (no backing array): True.
- Reused: EB-1 correctness True, Amendment 1 correctness True. Overall pass: **True**.

## 4. A batch-time model

Amendment 1 model (with intercept), refitted for reference: coefficients c0=0.0271, c1=1.83e-07, c2=2.96e-06; max fit residual 0.33 (criterion 0.25).

Replacement (frozen, chosen after seeing the Amendment 1 residuals): c1=1.78e-07, c2=1.07e-06, c3=0.00273; terms dropped: none; max fit residual **0.238** (criterion 0.25). Supported: **True**.

| Fit cell n_Q / n_D / n_win | measured CPU s | predicted | abs. error s | rel. error |
|---|---|---|---|---|
| 32 / 1015 / 20 | 0.155 | 0.192 | 0.037 | 0.24 |
| 64 / 1015 / 20 | 0.338 | 0.308 | -0.03 | 0.089 |
| 128 / 1015 / 20 | 0.554 | 0.54 | -0.014 | 0.025 |
| 192 / 1015 / 20 | 0.773 | 0.848 | 0.075 | 0.097 |
| 256 / 1015 / 20 | 1.08 | 1.08 | 0.002 | 0.0019 |
| 128 / 250 / 20 | 0.155 | 0.174 | 0.02 | 0.13 |
| 128 / 500 / 20 | 0.325 | 0.294 | -0.032 | 0.098 |
| 256 / 250 / 20 | 0.31 | 0.348 | 0.038 | 0.12 |
| 256 / 500 / 20 | 0.642 | 0.587 | -0.055 | 0.086 |
| 64 / 1015 / 5 | 0.0797 | 0.077 | -0.0027 | 0.034 |
| 64 / 1015 / 80 | 1.26 | 1.23 | -0.025 | 0.02 |

| Confirmation | role | measured CPU s | predicted | abs. error s | rel. error |
|---|---|---|---|---|---|
| FC1 (D1) | fresh, deciding | 0.131 | 0.125 | -0.0062 | 0.047 |
| FC2 (D1) | fresh, deciding | 0.955 | 0.854 | -0.1 | 0.11 |
| FC3 (D1p) | fresh, deciding | 12.9 | 12.1 | -0.85 | 0.066 |
| FC4 (D2) | fresh, deciding | 0.876 | 0.914 | 0.038 | 0.043 |
| CW1 (D1) | seen, not deciding | 21.9 | 22.3 | 0.44 | 0.02 |
| CW2 (D1p) | seen, not deciding | 5.18 | 5.37 | 0.18 | 0.035 |
| CW4 (D2) | seen, not deciding | 1.14 | 1.09 | -0.05 | 0.043 |

P2 value 1.88e+03 CPU s, extrapolation × 92.3 beyond the largest measured workload.

## 5. All S4 models (W = 8)

| Model | fit residual | confirmation max rel. error | supported | P2 value | extrapolation × |
|---|---|---|---|---|---|
| S4_A_batch_time | 0.24 | 0.11 | True | 1.88e+03 | 92.3 |
| S4_A_batch_mem | 0.011 | 0.054 | True | 6.7e+08 | 6.15 |
| S4_A_single_time | 0.016 | 0.033 | True | 2 | 17.8 |
| S4_A_single_mem | 0.0021 | 0.01 | True | 2.1e+08 | 18 |
| S4_B2_batch_time | 0.076 | 0.046 | True | 2.08e+04 | 92.3 |
| S4_B2_batch_mem | 0.0056 | 0.073 | True | 9.47e+08 | 5.9 |
| S4_B2_single_time | 0.043 | 0.015 | True | 25.5 | 17.8 |
| S4_B2_single_mem | 0.002 | 0.0026 | True | 2.89e+08 | 17.2 |

Units: batch time = CPU s; single time = wall s; memory = absolute peak RSS bytes.

## 6. Re-measured memory and reused times (descriptive)

- S2: peak RSS change A2 − A1 over 18 matched cells: min -712 kB, median 264 kB, max 412 kB.
- S3: peak RSS change A2 − A1 over 15 matched cells: min -212 kB, median 228 kB, max 396 kB.
- S4: peak RSS change A2 − A1 over 34 matched cells: min -1.05e+04 kB, median -60 kB, max 1.06e+03 kB.
- Times in re-measured cells vs reused Amendment 1 times (ratio A2/A1, 67 cells): min 0.951, median 1.05, max 1.41. Re-measured times do not replace reused ones.

## 7. Requirements (thresholds unchanged)

| Requirement | component | projected | limit | unit | ratio |
|---|---|---|---|---|---|
| R1_A | engine | 0.526 | 96 | CPU-h | 0.00548 |
| R2_A | engine | 2 | 10 | s | 0.2 |
| R4_S4_A_batch | engine | 6.54e+05 | 8388608 | kB | 0.078 |
| R4_S4_A_single | engine | 2.05e+05 | 8388608 | kB | 0.0245 |
| R1_B2 | engine | 5.79 | 96 | CPU-h | 0.0603 |
| R2_B2 | engine | 25.5 | 10 | s | 2.55 |
| R4_S4_B2_batch | engine | 9.25e+05 | 8388608 | kB | 0.11 |
| R4_S4_B2_single | engine | 2.82e+05 | 8388608 | kB | 0.0336 |
| R3_xtc | loader | 150 | 7.2e+03 | s | 0.0208 |
| R3_h5 | loader | 49.3 | 7.2e+03 | s | 0.00685 |
| R4_S1_xtc | loader | 2.08e+05 | 8388608 | kB | 0.0248 |
| R4_S1_h5 | loader | 4.22e+06 | 8388608 | kB | 0.503 |
| R4_S2 | loader | 8.07e+05 | 8388608 | kB | 0.0962 |
| R4_S3_A | engine | 6.02e+05 | 8388608 | kB | 0.0718 |

## 8. Limitations

- Shared desktop; page caches not dropped; timings conditional on contention.
- P2 values are extrapolations beyond the largest confirmed workload (factors above). Confirmation tests prediction within and moderately beyond the measured range, not at P2.
- The replacement A batch-time model was chosen after seeing the Amendment 1 residuals; only FC1–FC4 test it prospectively. Its fit residual is close to the 0.25 criterion.
- Memory transients use the SYNTHETIC Amendment 1 calibration; peak RSS includes allocator retention.
- The engine outcome rests on R2_B2 (single-query B2 at P2). No optimization exists or is measured here; the B2 single-query proposal is in `docs/EB1_B2_SINGLE_QUERY_OPTIMIZATION_PROPOSAL.md` (not implemented).

