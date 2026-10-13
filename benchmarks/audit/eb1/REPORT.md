# EB-1 engineering benchmark (development data only)

Protocol `docs/EB1_ENGINEERING_BENCHMARK_PROTOCOL.md` (SHA-256 `666f1591f17c556f…`, frozen before timing). development data only (AdK open_start; mdCATH 1ha8A00 schema sample). Started 2026-10-12T14:55:39.493339+00:00; wall 259.4 s; aggregate CPU 0.0719 h (budget 2 h). Cells: 108 of 108 completed. Engineering only: WP-D retrieval results are not re-estimated.

## Decision (frozen §7): **insufficient_evidence** — unsupported projections: ['R4_S4']

## 1. Correctness at each processing boundary (mandatory)

- C0 stored, D1: bit-identical True (MDAnalysis full-system positions indexed at Cα indices (same decoder); frames [0, 1, 500, 1237, 1238, 2252]). D2: bit-identical True (h5py whole coords dataset read then indexed at Cα indices).
- C1 reconstructed: D1 max |Δ| 0 Å vs WP-B derived open_start_ca_whole.npy (separate run of the same reader and make_whole algorithm); D2 max |Δ| 0 Å vs frozen harness make_whole_r1 along N–CA–C, Cα column.
- C2 features (binding vs NumPy, identical float32-rounded inputs): D1: κ 6.94e-16, τ 4.44e-16 rad; D2: κ 4.44e-16, τ 4.44e-16 rad
- C3 search (100 queries × window per cell; near ties by the frozen rule):

| Cell | top-k agree | near ties | max |Δd| | tolerance | pass |
|---|---|---|---|---|---|
| D1/W5/A | 100/100 | 0 | 0 | 1e-10 | True |
| D1/W5/B2 | 100/100 | 0 | 2.19e-13 | 1e-06 | True |
| D1/W8/A | 100/100 | 0 | 0 | 1e-10 | True |
| D1/W8/B2 | 100/100 | 0 | 4.13e-13 | 1e-06 | True |
| D2/W5/A | 100/100 | 0 | 0 | 1e-10 | True |
| D2/W5/B2 | 100/100 | 0 | 1.03e-13 | 1e-06 | True |
| D2/W8/A | 100/100 | 0 | 0 | 1e-10 | True |
| D2/W8/B2 | 100/100 | 0 | 8.89e-14 | 1e-06 | True |

Overall correctness pass: **True**.

## 2. Measured stage costs (operational path; single process; BLAS threads 1)

| Block | Args | Work x | wall s (median) | CPU s | absolute peak RSS MB | MemAvailable at start MB |
|---|---|---|---|---|---|---|
| M0 | {} | — | 4.02e-07 | 4.02e-07 | 94.5 | 1.36e+03 |
| S1 | {"d": "D1", "n": 563} | 53583525 | 4.04 | 3.99 | 111 | 1.37e+03 |
| S1 | {"d": "D1", "n": 1126} | 107167050 | 7.18 | 7.03 | 112 | 1.47e+03 |
| S1 | {"d": "D1", "n": 2253} | 214429275 | 13.5 | 13.3 | 114 | 1.44e+03 |
| S1 | {"d": "D2", "n": 1} | 91800 | 0.0531 | 0.0498 | 53.1 | 1.4e+03 |
| S1 | {"d": "D2", "n": 2} | 183600 | 0.0612 | 0.0574 | 55.6 | 1.35e+03 |
| S1 | {"d": "D2", "n": 4} | 363120 | 0.0788 | 0.0755 | 56.4 | 1.35e+03 |
| S2 | {"d": "D1", "n": 563} | 120482 | 0.51 | 0.51 | 47.1 | 1.36e+03 |
| S2 | {"d": "D1", "n": 1126} | 240964 | 1.02 | 1.02 | 53.9 | 1.41e+03 |
| S2 | {"d": "D1", "n": 2253} | 482142 | 2.04 | 2.04 | 67.8 | 1.46e+03 |
| S2 | {"d": "D2", "n": 1} | 22950 | 0.0975 | 0.0975 | 36.2 | 1.46e+03 |
| S2 | {"d": "D2", "n": 2} | 45900 | 0.194 | 0.194 | 36.4 | 1.46e+03 |
| S2 | {"d": "D2", "n": 4} | 90780 | 0.394 | 0.39 | 36.9 | 1.44e+03 |
| S3 | {"d": "D1", "n": 563} | 120482 | 0.197 | 0.196 | 74.1 | 1.45e+03 |
| S3 | {"d": "D1", "n": 1126} | 240964 | 0.394 | 0.39 | 101 | 1.47e+03 |
| S3 | {"d": "D1", "n": 2253} | 482142 | 0.784 | 0.784 | 156 | 1.48e+03 |
| S3 | {"d": "D2", "n": 405} | 20655 | 0.0309 | 0.0309 | 41.5 | 1.45e+03 |
| S3 | {"d": "D2", "n": 810} | 41310 | 0.0617 | 0.0617 | 46.4 | 1.43e+03 |
| S3 | {"d": "D2", "n": 1620} | 81702 | 0.123 | 0.123 | 55.2 | 1.41e+03 |
| S4single | {"W": 5, "d": "D1", "frac": 0.25, "m": "A"} | 53130 | 0.026 | 0.026 | 57.1 | 1.41e+03 |
| S4single | {"W": 5, "d": "D1", "frac": 0.5, "m": "A"} | 106470 | 0.0343 | 0.0343 | 57.2 | 1.41e+03 |
| S4single | {"W": 5, "d": "D1", "frac": 1.0, "m": "A"} | 213150 | 0.049 | 0.049 | 57.1 | 1.48e+03 |
| S4batch | {"W": 5, "d": "D1", "m": "A", "nq": 32} | 682080 | 0.0974 | 0.0974 | 56.5 | 1.45e+03 |
| S4batch | {"W": 5, "d": "D1", "m": "A", "nq": 64} | 1364160 | 0.224 | 0.223 | 59.2 | 1.49e+03 |
| S4batch | {"W": 5, "d": "D1", "m": "A", "nq": 128} | 2728320 | 0.423 | 0.422 | 66.7 | 1.46e+03 |
| S4single | {"W": 5, "d": "D1", "frac": 0.25, "m": "B2"} | 53130 | 0.17 | 0.17 | 57.7 | 1.45e+03 |
| S4single | {"W": 5, "d": "D1", "frac": 0.5, "m": "B2"} | 106470 | 0.294 | 0.294 | 57.7 | 1.46e+03 |
| S4single | {"W": 5, "d": "D1", "frac": 1.0, "m": "B2"} | 213150 | 0.542 | 0.541 | 57.7 | 1.44e+03 |
| S4batch | {"W": 5, "d": "D1", "m": "B2", "nq": 32} | 682080 | 1.52 | 1.52 | 65.4 | 1.36e+03 |
| S4batch | {"W": 5, "d": "D1", "m": "B2", "nq": 64} | 1364160 | 3.08 | 3.08 | 76.8 | 1.4e+03 |
| S4batch | {"W": 5, "d": "D1", "m": "B2", "nq": 128} | 2728320 | 6.37 | 6.37 | 99.2 | 1.41e+03 |
| S4single | {"W": 8, "d": "D1", "frac": 0.25, "m": "A"} | 52371 | 0.0286 | 0.0286 | 57.2 | 1.43e+03 |
| S4single | {"W": 8, "d": "D1", "frac": 0.5, "m": "A"} | 104949 | 0.0388 | 0.0388 | 57.2 | 1.43e+03 |
| S4single | {"W": 8, "d": "D1", "frac": 1.0, "m": "A"} | 210105 | 0.0606 | 0.0606 | 57.3 | 1.42e+03 |
| S4batch | {"W": 8, "d": "D1", "m": "A", "nq": 32} | 682080 | 0.19 | 0.19 | 61.4 | 1.4e+03 |
| S4batch | {"W": 8, "d": "D1", "m": "A", "nq": 64} | 1364160 | 0.389 | 0.389 | 69.2 | 1.4e+03 |
| S4batch | {"W": 8, "d": "D1", "m": "A", "nq": 128} | 2728320 | 0.775 | 0.774 | 84.6 | 1.4e+03 |
| S4single | {"W": 8, "d": "D1", "frac": 0.25, "m": "B2"} | 52371 | 0.176 | 0.176 | 57.7 | 1.42e+03 |
| S4single | {"W": 8, "d": "D1", "frac": 0.5, "m": "B2"} | 104949 | 0.311 | 0.311 | 57.7 | 1.41e+03 |
| S4single | {"W": 8, "d": "D1", "frac": 1.0, "m": "B2"} | 210105 | 0.626 | 0.626 | 57.8 | 1.39e+03 |
| S4batch | {"W": 8, "d": "D1", "m": "B2", "nq": 32} | 682080 | 1.6 | 1.6 | 65.6 | 1.38e+03 |
| S4batch | {"W": 8, "d": "D1", "m": "B2", "nq": 64} | 1364160 | 3.14 | 3.14 | 77 | 1.39e+03 |
| S4batch | {"W": 8, "d": "D1", "m": "B2", "nq": 128} | 2728320 | 6.49 | 6.48 | 99.4 | 1.4e+03 |
| S4single | {"W": 5, "d": "D2", "frac": 0.25, "m": "A"} | 4747 | 0.00457 | 0.00457 | 40.4 | 1.4e+03 |
| S4single | {"W": 5, "d": "D2", "frac": 0.5, "m": "A"} | 9494 | 0.00546 | 0.00546 | 40.2 | 1.39e+03 |
| S4single | {"W": 5, "d": "D2", "frac": 1.0, "m": "A"} | 19035 | 0.00679 | 0.00679 | 40.4 | 1.39e+03 |
| S4batch | {"W": 5, "d": "D2", "m": "A", "nq": 32} | 609120 | 0.0908 | 0.09 | 37.8 | 1.39e+03 |
| S4batch | {"W": 5, "d": "D2", "m": "A", "nq": 64} | 1218240 | 0.177 | 0.175 | 39.2 | 1.37e+03 |
| S4batch | {"W": 5, "d": "D2", "m": "A", "nq": 128} | 2436480 | 0.383 | 0.379 | 42.3 | 1.37e+03 |
| S4single | {"W": 5, "d": "D2", "frac": 0.25, "m": "B2"} | 4747 | 0.0218 | 0.0217 | 41 | 1.37e+03 |
| S4single | {"W": 5, "d": "D2", "frac": 0.5, "m": "B2"} | 9494 | 0.0326 | 0.0326 | 40.8 | 1.37e+03 |
| S4single | {"W": 5, "d": "D2", "frac": 1.0, "m": "B2"} | 19035 | 0.0562 | 0.0562 | 41 | 1.35e+03 |
| S4batch | {"W": 5, "d": "D2", "m": "B2", "nq": 32} | 609120 | 1.39 | 1.39 | 41.7 | 1.33e+03 |
| S4batch | {"W": 5, "d": "D2", "m": "B2", "nq": 64} | 1218240 | 2.78 | 2.78 | 46.2 | 1.33e+03 |
| S4batch | {"W": 5, "d": "D2", "m": "B2", "nq": 128} | 2436480 | 5.48 | 5.48 | 55.6 | 1.33e+03 |
| S4single | {"W": 8, "d": "D2", "frac": 0.25, "m": "A"} | 4444 | 0.00491 | 0.00491 | 40.4 | 1.34e+03 |
| S4single | {"W": 8, "d": "D2", "frac": 0.5, "m": "A"} | 8888 | 0.00575 | 0.00575 | 40.3 | 1.32e+03 |
| S4single | {"W": 8, "d": "D2", "frac": 1.0, "m": "A"} | 17820 | 0.00752 | 0.00752 | 40.4 | 1.32e+03 |
| S4batch | {"W": 8, "d": "D2", "m": "A", "nq": 32} | 570240 | 0.154 | 0.154 | 39.6 | 1.32e+03 |
| S4batch | {"W": 8, "d": "D2", "m": "A", "nq": 64} | 1140480 | 0.317 | 0.317 | 43.2 | 1.32e+03 |
| S4batch | {"W": 8, "d": "D2", "m": "A", "nq": 128} | 2280960 | 0.679 | 0.679 | 49.5 | 1.32e+03 |
| S4single | {"W": 8, "d": "D2", "frac": 0.25, "m": "B2"} | 4444 | 0.021 | 0.0211 | 40.9 | 1.32e+03 |
| S4single | {"W": 8, "d": "D2", "frac": 0.5, "m": "B2"} | 8888 | 0.0326 | 0.0322 | 40.9 | 1.32e+03 |
| S4single | {"W": 8, "d": "D2", "frac": 1.0, "m": "B2"} | 17820 | 0.0563 | 0.0555 | 41 | 1.31e+03 |
| S4batch | {"W": 8, "d": "D2", "m": "B2", "nq": 32} | 570240 | 1.36 | 1.36 | 41.8 | 1.32e+03 |
| S4batch | {"W": 8, "d": "D2", "m": "B2", "nq": 64} | 1140480 | 2.71 | 2.71 | 46.2 | 1.31e+03 |
| S4batch | {"W": 8, "d": "D2", "m": "B2", "nq": 128} | 2280960 | 5.34 | 5.34 | 55.8 | 1.34e+03 |
| R | {"part": "targets"} | — | 1.28 | 1.28 | 124 | 1.35e+03 |
| R | {"part": "metrics"} | — | 2.09 | 2.08 | 126 | 1.4e+03 |
| R | {"part": "bootstrap"} | — | 0.413 | 0.412 | 39 | 1.32e+03 |

S1 first access vs repeats (wall s): {'d': 'D1', 'n': 2253} rep0: 13.6; {'d': 'D1', 'n': 2253} rep1: 13.4; {'d': 'D1', 'n': 2253} rep2: 13.5; {'d': 'D2', 'n': 4} rep0: 0.115; {'d': 'D2', 'n': 4} rep1: 0.0788; {'d': 'D2', 'n': 4} rep2: 0.0772

Baseline import-only peak RSS: 94.5 MB (supplementary; decisions use absolute peak RSS).

## 3. Fits and P2 projections (frozen §6; P2 = 1,000 residues, W = 8, 10,000 database frames, 1,000 batch queries)

| Fit | supported | c0 | c1 | max rel. residual | P2 projection | extrapolation × |
|---|---|---|---|---|---|---|
| S1_xtc_wall | True | 0.881 | 5.89e-08 | 0.000907 | 104 | 8.21 |
| S1_h5_wall | True | 0.0441 | 9.5e-08 | 0.0062 | 4.22 | 121 |
| S2_D1_wall | True | -0.00183 | 4.24e-06 | 0.00307 | 46.6 | 22.8 |
| S2_D2_wall | True | -0.00476 | 4.38e-06 | 0.0173 | 48.2 | 121 |
| S3_D1_cpu | True | -0.000151 | 1.63e-06 | 0.00377 | 17.9 | 22.8 |
| S3_D2_cpu | True | -0.000714 | 1.52e-06 | 0.00735 | 16.7 | 135 |
| S3_D1_wall_one_frame | True | 0.00175 | 1.62e-06 | 0.00399 | 0.00337 | 0.00207 |
| S3_D2_wall_one_frame | True | -0.000704 | 1.52e-06 | 0.00731 | 0.000814 | 0.0122 |
| S4single_D1_W8_A_wall | True | 0.0177 | 2.04e-07 | 0.00689 | 2.04 | 47.3 |
| S4single_D2_W8_A_wall | True | 0.00402 | 1.96e-07 | 0.00346 | 1.95 | 557 |
| S4batch_D1_W8_A_cpu | True | -0.00284 | 2.85e-07 | 0.00907 | 2.83e+03 | 3.64e+03 |
| S4batch_D2_W8_A_cpu | True | -0.0268 | 3.08e-07 | 0.0321 | 3.06e+03 | 4.35e+03 |
| S4single_D1_W8_B2_wall | True | 0.0191 | 2.87e-06 | 0.0369 | 28.5 | 47.3 |
| S4single_D2_W8_B2_wall | True | 0.00926 | 2.64e-06 | 0.0036 | 26.2 | 557 |
| S4batch_D1_W8_B2_cpu | True | -0.0748 | 2.4e-06 | 0.0242 | 2.38e+04 | 3.64e+03 |
| S4batch_D2_W8_B2_cpu | True | 0.048 | 2.32e-06 | 0.00636 | 2.31e+04 | 4.35e+03 |
| R4_S1_xtc | True | 1.12e+05 | 0.000743 | 0.00274 | 2.11e+05 | 22.8 |
| R4_S1_h5 | True | 5.4e+04 | 0.00376 | 0.0155 | 5.5e+05 | 121 |
| R4_S2 | True | 3.52e+04 | 0.00123 | 0.108 | 8.46e+05 | 22.8 |
| R4_S3 | True | 3.9e+04 | 0.00637 | 0.0812 | 2.84e+06 | 22.8 |
| R4_S4 | False | — | — | — | — | — |

## 4. Requirements (absolute)

| Requirement | component | projected | limit | unit | ratio | violated |
|---|---|---|---|---|---|---|
| R1_A | engine | 0.855 | 96 | CPU-h | — | None |
| R2_A | engine | 2.04 | 10 | s | — | None |
| R1_B2 | engine | 6.62 | 96 | CPU-h | — | None |
| R2_B2 | engine | 28.5 | 10 | s | — | None |
| R3_xtc | loader | 153 | 7.2e+03 | s | — | None |
| R3_h5 | loader | 52.4 | 7.2e+03 | s | — | None |
| R4_S1_xtc | loader | 2.11e+05 | 8388608 | kB | — | None |
| R4_S1_h5 | loader | 5.5e+05 | 8388608 | kB | — | None |
| R4_S2 | loader | 8.46e+05 | 8388608 | kB | — | None |
| R4_S3 | engine | 2.84e+06 | 8388608 | kB | — | None |
| R4_S4 | engine | — | 8388608 | kB | — | None |

## 5. Research-only path (reported, not deciding)

- targets: wall 1.28 s, CPU 1.28 s, peak RSS 124 MB
- metrics: wall 2.09 s, CPU 2.08 s, peak RSS 126 MB
- bootstrap: wall 0.413 s, CPU 0.412 s, peak RSS 39 MB

M8 parallel efficiency: {'status': 'skipped', 'reason': 'MemAvailable 1337128 kB < 2 GB (protocol §5)'}.

## 6. Post-hoc observations (after results; NOT part of the frozen decision)

- **Why R4_S4 is unsupported:** the S4 memory work variable `resident_bytes` was specified with the fixed batch size 128 rather than each cell's actual query count, so it took only 2 distinct values (one per dataset) and failed the frozen ≥ 3-size rule. Measured absolute peak RSS nevertheless rose with the batch query count (e.g. D1 W8 B2 batch: n_Q 32 → 65.6 MB, n_Q 64 → 77 MB, n_Q 128 → 99.4 MB). This is a specification error in EB-1, not a measurement failure; the frozen outcome stands.
- Largest absolute peak RSS among all operational cells: 156 MB.
- Ratios of the supported projections to their limits (descriptive): R1_A 0.00891; R2_A 0.204; R1_B2 0.0689; R2_B2 2.85; R3_xtc 0.0212; R3_h5 0.00728; R4_S1_xtc 0.0251; R4_S1_h5 0.0656; R4_S2 0.101; R4_S3 0.339.
- Among supported projections, only R2_B2 (single-query, all-window latency for window Cα RMSD) exceeds its limit; no loader projection is close to its limit. Because one required projection is unsupported, this is not a decision.

## 7. Limitations

- Shared desktop under memory pressure (MemAvailable per cell above; swap full); timings are conditional on it and page caches could not be dropped (first access is labelled).
- C0 references use the same decoders (MDAnalysis, h5py): consistency checks, not independent decoders. C1 for D1 is a separate run of the same algorithm.
- P2 projections extrapolate linear fits far beyond the measured sizes (factors above); the linear models are analytically expected for these stages but are not measured at P2.
- Only existing implementations are measured; no saving is claimed for an implementation that does not exist.
- Two development datasets only; other formats and system sizes are untested.

