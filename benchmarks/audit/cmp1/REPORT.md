# CMP-1: bounded competitive evaluation of window Cα RMSD retrieval (harness; development + SYNTHETIC data)

Protocol `docs/CMP1_COMPETITIVE_EVALUATION_PROTOCOL.md` (SHA-256 `3ec41b8b98c3dd42…`, frozen 2026-10-13T14:53:22.046224+00:00); code hashed 2026-10-13T14:53:22.047297+00:00 before deciding measurements. Environment: Python 3.14.4; MDAnalysis==2.10.0, mdtraj==1.11.2, numpy==2.5.3, scipy==1.18.1; isolated venv; wheel hashes in `results/cmp1_environment.json`. CPPTRAJ excluded (protocol §2). This is a bounded comparison of the listed implementations on the listed workloads, not a field-wide best-tool claim.

CPU used: 0.521 h of the 2 CPU-h cap. Stopped early: no. Timed cells completed: 85 of 100.

## Conclusion

- **Primary (common F32 accuracy track): NO_DECISIVE_ADVANTAGE** — per workload: R1_D1 V2b_advantage; R2_D1plus comparator_preferable; R3_D2 comparator_preferable; S1_P2slice V2b_advantage; S2_P2D V2b_advantage.
- Secondary (strict F64 track, float64-accurate implementations only): NO_DECISIVE_ADVANTAGE — R1_D1 V2b_advantage; R2_D1plus V2b_advantage; R3_D2 no_decisive_difference; S1_P2slice V2b_advantage; S2_P2D V2b_advantage. The secondary track cannot by itself justify integration (protocol §7).
- V2b degeneracy suite (F64): True.

## 1. Accuracy (verification outside timing; identical requirements per track)

F64: |d − d_ref| ≤ 1e-6 Å, rank τ = 2e-6 Å. F32: |d − d_ref| ≤ 1e-3 Å, rank τ = 2e-3 Å. Reference: float64 SVD Kabsch on the canonical input. Native ratio = max |d² − d_ref²| / (u·(G_q + G_j)/W), descriptive only.

| Workload | impl | F64 pass | F32 pass | max |Δd| Å | native ratio | rank non-eq. (F32) | strict identical (F32) | near-tie (F32) | input unchanged |
|---|---|---|---|---|---|---|---|---|---|
| R1_D1 | T-V2b | True | True | 2.13e-12 | 1.25e-06 | 0 | 621/621 | 589 | True |
| R1_D1 | T-V0 | True | True | 2.83e-13 | 4.29e-08 | 0 | 621/621 | 589 | True |
| R1_D1 | E-MDT | False | False | 0.00346 | 2.04e+03 | 0 | 608/621 | 589 | True |
| R1_D1 | E-MDA | True | True | 2.53e-12 | 3.63e-06 | 0 | 621/621 | 589 | True |
| R2_D1plus | T-V2b | True | True | 2.4e-12 | 1.17e-06 | 0 | 621/621 | 592 | True |
| R2_D1plus | T-V0 | True | True | 4.86e-13 | 3.94e-08 | 0 | 621/621 | 592 | True |
| R2_D1plus | E-MDT | False | True | 0.000881 | 541 | 0 | 611/621 | 592 | True |
| R2_D1plus | E-MDA | True | True | 1.53e-12 | 7.6e-07 | 0 | 621/621 | 592 | True |
| R3_D2 | T-V2b | True | True | 1.82e-13 | 3.55e-07 | 0 | 132/132 | 101 | True |
| R3_D2 | T-V0 | True | True | 1.05e-13 | 3.32e-08 | 0 | 132/132 | 101 | True |
| R3_D2 | E-MDT | False | True | 7.09e-05 | 136 | 0 | 130/132 | 101 | True |
| R3_D2 | E-MDA | True | True | 1.62e-13 | 3.1e-07 | 0 | 132/132 | 101 | True |
| S1_P2slice | T-V2b | True | True | 3.37e-07 | 1.56e-06 | 0 | 300/300 | 300 | True |
| S1_P2slice | T-V0 | True | True | 2.06e-07 | 4.3e-08 | 0 | 300/300 | 300 | True |
| S1_P2slice | E-MDT | False | False | 0.012 | 1.26e+03 | 0 | 298/300 | 300 | True |
| S1_P2slice | E-MDA | True | True | 5.33e-07 | 1.61e-06 | 0 | 300/300 | 300 | True |
| S2_P2D | T-V2b | True | True | 4.77e-07 | 3.1e-06 | 0 | 2979/2979 | 2979 | True |
| S2_P2D | T-V0 | True | True | 2.92e-07 | 4.64e-08 | 0 | 2979/2979 | 2979 | True |
| S2_P2D | E-MDT | False | False | 0.0109 | 3.88e+03 | 0 | 2969/2979 | 2979 | True |
| S2_P2D | E-MDA | True | True | 7.54e-07 | 2.84e-05 | 0 | 2979/2979 | 2979 | True |

**Competitor- and candidate-specific numerical limitations (listed separately; never sufficient for integration):**

- E-MDT / F64 / R1_D1: max |Δd| 0.00346 Å, rank non-eq. 11
- E-MDT / F32 / R1_D1: max |Δd| 0.00346 Å, rank non-eq. 0
- E-MDT / F64 / R2_D1plus: max |Δd| 0.000881 Å, rank non-eq. 10
- E-MDT / F64 / R3_D2: max |Δd| 7.09e-05 Å, rank non-eq. 1
- E-MDT / F64 / S1_P2slice: max |Δd| 0.012 Å, rank non-eq. 2
- E-MDT / F32 / S1_P2slice: max |Δd| 0.012 Å, rank non-eq. 0
- E-MDT / F64 / S2_P2D: max |Δd| 0.0109 Å, rank non-eq. 8
- E-MDT / F32 / S2_P2D: max |Δd| 0.0109 Å, rank non-eq. 0
- E-MDT / F64 / SYNTHETIC degeneracy suite: failed ['query_identical_to_database_rows', 'collinear_query', 'planar_query', 'mirror_image_query', 'permutation', 'rigid_transform', 'near_ties']; errors query_identical_to_database_rows_max_abs_err_A=0.00173, collinear_query_max_abs_err_A=0.000793, planar_query_max_abs_err_A=0.00598, mirror_image_query_max_abs_err_A=0.00173
- E-MDT / F32 / SYNTHETIC degeneracy suite: failed ['query_identical_to_database_rows', 'planar_query', 'mirror_image_query']; errors query_identical_to_database_rows_max_abs_err_A=0.00173, collinear_query_max_abs_err_A=0.000793, planar_query_max_abs_err_A=0.00598, mirror_image_query_max_abs_err_A=0.00173
- E-MDA / F64 / SYNTHETIC degeneracy suite: failed ['collinear_query']; errors query_identical_to_database_rows_max_abs_err_A=1.19e-07, collinear_query_max_abs_err_A=inf, planar_query_max_abs_err_A=2.92e-07, mirror_image_query_max_abs_err_A=1.19e-07
- E-MDA / F32 / SYNTHETIC degeneracy suite: failed ['collinear_query']; errors query_identical_to_database_rows_max_abs_err_A=1.19e-07, collinear_query_max_abs_err_A=inf, planar_query_max_abs_err_A=2.92e-07, mirror_image_query_max_abs_err_A=1.19e-07

MDTraj semantics (pinned source and adapter-development checks): distances are float32; `precentered` with `atom_indices` is ignored with a warning (True); the `atom_indices=None` path centres the input in place (True); pre-sliced precentered windows agree with subset RMSDs to 1.15e-06 nm. MDTraj unit conversion (Å → nm float32 → Å) representation error on development cases: dev_D1 2.27e-06 Å, dev_D2 1.98e-06 Å, dev_tiled 2.99e-06 Å.

## 2. Performance (fresh process per cell, one thread, 5 seeded rounds; medians, with ranges)

L2 = complete retrieval workflow per query (all windows, shared top-k). L1 = solver-kernel time within the same queries (distance computation only). First-query span = preparation + query 1, excluding imports, warm-up and canonical input loading. Incremental peak RSS = lifetime peak − RSS after imports and warm-up.

| Workload | impl | setup s | prepare s | first-query span s | repeated query L2 s (range) | kernel L1 s | peak RSS MB (incr.) | prepared MB | reps ok | outputs = gate |
|---|---|---|---|---|---|---|---|---|---|---|
| R1_D1 (real) | T-V2b | 0.0223 | 0.087 | 0.312 | 0.224 (0.22–0.294) | 0.216 | 89.1 (49) | 42 | 5/5 | True |
| R1_D1 (real) | T-V0 | 0.0246 | 1.02e-05 | 0.614 | 0.618 (0.557–0.626) | 0.603 | 54.5 (13.3) | — | 5/5 | True |
| R1_D1 (real) | E-MDT | — | — | — | — | — | — (—) | — | 0/5 | None |
| R1_D1 (real) | E-MDA | 0.434 | 0.0873 | 0.591 | 0.489 (0.444–0.504) | 0.466 | 136 (49.8) | 40.3 | 5/5 | True |
| R2_D1plus (real) | T-V2b | 0.0238 | 0.246 | 0.723 | 0.486 (0.471–0.552) | 0.475 | 169 (129) | 111 | 5/5 | True |
| R2_D1plus (real) | T-V0 | 0.0221 | 1.24e-05 | 1.43 | 1.44 (1.43–1.65) | 1.43 | 75.1 (33.8) | — | 5/5 | True |
| R2_D1plus (real) | E-MDT | 0.142 | 0.186 | 0.32 | 0.142 (0.125–0.146) | 0.121 | 135 (78.9) | 55.7 | 5/5 | True |
| R2_D1plus (real) | E-MDA | 0.386 | 0.218 | 1.43 | 1.23 (1.16–1.28) | 1.21 | 216 (130) | 107 | 5/5 | True |
| R3_D2 (real) | T-V2b | 0.0226 | 0.00948 | 0.0386 | 0.0303 (0.0284–0.0336) | 0.0286 | 44.8 (4.87) | 3.56 | 5/5 | True |
| R3_D2 (real) | T-V0 | 0.0226 | 8.84e-06 | 0.0555 | 0.0555 (0.0539–0.0611) | 0.0536 | 43.2 (1.99) | — | 5/5 | True |
| R3_D2 (real) | E-MDT | 0.138 | 0.00785 | 0.0151 | 0.00725 (0.00699–0.00822) | 0.00457 | 59.7 (3.14) | 1.78 | 5/5 | True |
| R3_D2 (real) | E-MDA | 0.393 | 0.00867 | 0.0484 | 0.0402 (0.0391–0.0425) | 0.0379 | 91.9 (5.75) | 3.42 | 5/5 | True |
| S1_P2slice (SYNTHETIC_tiled) | T-V2b | 0.0225 | 0.473 | 1.26 | 0.779 (0.763–0.794) | 0.769 | 637 (597) | 200 | 5/5 | True |
| S1_P2slice (SYNTHETIC_tiled) | T-V0 | 0.0219 | 1.05e-05 | 2.56 | 2.55 (2.52–2.94) | 2.54 | 638 (597) | — | 5/5 | True |
| S1_P2slice (SYNTHETIC_tiled) | E-MDT | — | — | — | — | — | — (—) | — | 0/5 | None |
| S1_P2slice (SYNTHETIC_tiled) | E-MDA | 0.439 | 0.55 | 2.77 | 2.23 (2.07–2.38) | 2.21 | 683 (597) | 192 | 5/5 | True |
| S2_P2D (SYNTHETIC_tiled) | T-V2b | 0.0218 | 6.38 | 14.4 | 7.85 (7.78–9.17) | 7.75 | 2.19e+03 (2.15e+03) | 1.99e+03 | 5/5 | True |
| S2_P2D (SYNTHETIC_tiled) | T-V0 | 0.0235 | 1.1e-05 | 25.5 | 25.5 (25.2–27.8) | 25.4 | 638 (597) | — | 5/5 | True |
| S2_P2D (SYNTHETIC_tiled) | E-MDT | — | — | — | — | — | — (—) | — | 0/5 | None |
| S2_P2D (SYNTHETIC_tiled) | E-MDA | 0.398 | 6.41 | 27.4 | 21.1 (20.9–22.2) | 20.9 | 2.16e+03 (2.07e+03) | 1.91e+03 | 5/5 | True |

## 3. Per-workload comparison (primary F32 track)

| Workload | E* | outcome | V2b/E* median repeated query | paired E*/V2b per round | memory criterion | ratios to V2b: L2, L1, first query, prepare, incr. memory (descriptive) |
|---|---|---|---|---|---|---|
| R1_D1 | E-MDA | V2b_advantage | 0.459 | 1.75, 2.2, 2.22, 2.03, 1.8 | True | T-V0: 2.8, 2.8, 2, 0.00012, 0.27; E-MDT: —, —, —, —, —; E-MDA: 2.2, 2.2, 1.9, 1, 1 |
| R2_D1plus | E-MDT | comparator_preferable | 3.44 | 0.264, 0.261, 0.298, 0.261, 0.274 | False | T-V0: 3, 3, 2, 5e-05, 0.26; E-MDT: 0.29, 0.25, 0.44, 0.76, 0.61; E-MDA: 2.5, 2.5, 2, 0.89, 1 |
| R3_D2 | E-MDT | comparator_preferable | 4.18 | 0.265, 0.242, 0.251, 0.242, 0.226 | False | T-V0: 1.8, 1.9, 1.4, 0.00093, 0.41; E-MDT: 0.24, 0.16, 0.39, 0.83, 0.65; E-MDA: 1.3, 1.3, 1.3, 0.91, 1.2 |
| S1_P2slice | E-MDA | V2b_advantage | 0.349 | 3.05, 2.85, 2.71, 2.67, 2.9 | True | T-V0: 3.3, 3.3, 2, 2.2e-05, 1; E-MDT: —, —, —, —, —; E-MDA: 2.9, 2.9, 2.2, 1.2, 1 |
| S2_P2D | E-MDA | V2b_advantage | 0.372 | 2.4, 2.41, 2.67, 2.69, 2.7 | True | T-V0: 3.3, 3.3, 1.8, 1.7e-06, 0.28; E-MDT: —, —, —, —, —; E-MDA: 2.7, 2.7, 1.9, 1, 0.97 |

Ratios > 1 mean the other implementation takes longer (or uses more memory) than V2b.

## 4. Limitations

- Bounded comparison: four implementations, one machine, one thread, shared desktop. Not a field-wide best-tool claim. CPPTRAJ was not evaluated (protocol §2); multithreaded and GPU paths were not evaluated.
- Workloads are reused development data (AdK, mdCATH 1ha8A00) and SYNTHETIC_tiled performance-only rows. This is engineering evidence, not independent scientific validation.
- The F32 tolerance was re-based on input resolution after adapter-development checks showed the first (a-priori msd-bound) proposal did not model float32 coefficient error. This was disclosed and frozen before any deciding measurement (protocol §4).
- An existing library being preferable for this task would not be a failure of the broader TopoFold mission; TopoFold's scientific questions are separate from which RMSD engine it uses.

