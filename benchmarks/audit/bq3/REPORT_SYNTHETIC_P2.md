# BQ-3: direct P2 qualification of unchanged V2b (SYNTHETIC_tiled, performance-only)

Protocol `docs/B2_DIRECT_P2_QUALIFICATION_PROTOCOL.md` (SHA-256 `4c2b6fecc4decd20…`, frozen 2026-10-13T10:30:19.991107+00:00, starting checkpoint `f55def1`); code hashed 2026-10-13T10:32:40.398214+00:00 before measurement. Harness only. BQ-1, BQ-2 (conclusion C) and EB-1 records are unchanged. Not pushed.

**Workload:** n_D = 10,000, L = 1,000, W = 8, 993 windows, 3 queries. Rows are real D1+ frames joined end to end (SYNTHETIC_tiled, repeated after 3,706 rows). This tests computational scale at P2 dimensions only. It says nothing about real 1,000-residue proteins or 10,000-frame trajectories, retrieval quality or biology.

CPU used: 0.104 h of the 1 CPU-h cap. Stopped early: no.

## Verdicts (frozen rule: every repetition must pass)

- **Prepared-database operation:** qualified
- **First-query operation (preparation + first query):** not_qualified_measured_limit
- Verification: True ; complete: True ; peak RSS within 8 GiB in every repetition: True.

## 1. Verification (outside timing)

- Failure-path tests: 6 passed in 0.03s.
- Reference: status ok, finite True.
- Gate: status ok; candidate finite True; max |Δd| 5.84e-07 Å (tolerance 1e-6); rank non-equivalent 0 of 2979 query-windows; strict identical 2979; near-tie 2979; fallback rows 0; gate peak RSS 2.47 GiB.

## 2. Direct measurements (all repetitions)

| Rep | start (UTC) | P1 s | Q1 s | first-query span s | Q2 s | Q3 s | repeated prep s | peak RSS GiB | fallback rows | load avg start | MemAvailable start GB | SwapFree start MB | outputs = gate |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 2026-10-13T10:35:35.950633+00:00 | 5.34 | 8.05 | 13.4 | 8.03 | 7.8 | 5.2, 5.71 | 2.13 | 0 | 7.65 | 5.78 | 0.24 | True |
| 1 | 2026-10-13T10:36:17.119517+00:00 | 5.07 | 8.01 | 13.1 | 7.97 | 7.98 | 6.06, 5.6 | 2.13 | 0 | 6.30 | 5.86 | 0.04 | True |
| 2 | 2026-10-13T10:36:58.774688+00:00 | 5.13 | 8.01 | 13.1 | 8.1 | 7.76 | 5.03, 5.19 | 2.13 | 0 | 4.27 | 5.81 | 0.04 | True |
| 3 | 2026-10-13T10:37:38.970572+00:00 | 4.91 | 7.92 | 12.8 | 7.78 | 7.73 | 4.85, 4.93 | 2.13 | 0 | 3.05 | 5.87 | 0.068 | True |
| 4 | 2026-10-13T10:38:17.980023+00:00 | 5.1 | 7.71 | 12.8 | 7.73 | 7.87 | 5.28, 5.13 | 2.13 | 0 | 2.57 | 5.88 | 0.012 | True |

- Initial preparation P1: median 5.1 s (min 4.91, max 5.34; n = 5); CPU median 5.07 s (min 4.91, max 5.33; n = 5).
- First query Q1: median 8.01 s (min 7.71, max 8.05; n = 5). First-query span (measured end to end): median 13.1 s (min 12.8, max 13.4; n = 5).
- Subsequent queries Q2–Q3: median 7.83 s (min 7.73, max 8.1; n = 10). All queries: median 7.92 s (min 7.71, max 8.1; n = 15).
- Repeated preparations: median 5.19 s (min 4.85, max 6.06; n = 10); CPU median 5.19 s (min 4.84, max 6.06; n = 10).
- Peak RSS: median 2.13 GiB (min 2.13, max 2.13; n = 5); prepared arrays 1.99 GB, resident rows 0.24 GB.
- Query CPU/wall: median 0.999 (min 0.991, max 1; n = 5). OS threads per process: [1, 1, 1, 1, 1]; threadpools: [{'api': 'openblas', 'num_threads': 1}].

## 3. R1 (not measured here)

R1_B2 remains the EB-1 Amendment 2 batch **projection** (5.79 CPU-h). Adding the directly measured P2 preparation CPU (median 5.07 s) gives 5.79 CPU-h against the 96 CPU-h limit. This mixes a projection with a measurement, is descriptive only, and does not decide.

## 4. Limitations

- SYNTHETIC_tiled rows are performance-only: repeated real frames, with artificial junction windows. Fallback frequency, cache behaviour and timings may differ on real P2-scale data.
- Ranking verification is weak on this workload: 2979 of 2979 query-windows contain near ties, because repeated tiled rows create exact duplicate windows. Rankings were nevertheless strictly identical to the reference in 2979 query-windows. Distance differences up to the reported maximum sit at the square-root rounding floor of RMSD ≈ 0 duplicates.
- One machine, one thread, shared desktop with swap nearly exhausted (recorded per repetition). Page caches were not dropped.
- Five repetitions characterize this machine and session only.
- Disclosure: after the measurement, `make_report.py` gained the ranking-verification limitation above (report text only). The acceptance code (`bq3_lib.py`, `bq3_analyse.py`) is unchanged from its pre-measurement hash.
- A qualification here is an engineering result for a harness implementation. It establishes no biological capability, no novelty over existing RMSD software and no production readiness.

