# EB-1 Amendment 2: Remaining Defects (frozen 2026-10-12)

**Amends:** EB-1 and Amendment 1. Their files, results and outcomes (both **insufficient evidence**) are unchanged.
- Amendment 1 did not select engine work. Its non-deciding observation is not reinterpreted as a decision.
- Recorded after the Amendment 1 results were seen and after a read-only preflight (`results/eb1a2_preflight.json`;
  not a performance measurement), and before any Amendment 2 measurement.
- Thresholds R1–R4, the P2 workload, the confirmation criterion (≤ 25 %), the fit-residual criterion (≤ 0.25) and
  the decision logic are **unchanged**.

## §1. Preflight findings (read-only)

- Frames, from the actual arrays:
  - D1: 2,253 stored frames, no burn-in.
  - D2: 450 / 450 / 430 / 450 stored frames (1,780 total) before burn-in; 405 / 405 / 387 / 405 (1,602) after.
  - D1+: 3,706 rows (1,013 queries + 2,693 database).
- The only invalid Amendment 1 request is **S3 D2 n = 1,620 > 1,602**. All other requested and actual workload sizes
  and window counts agree.
- 70 of the 80 Amendment 1 memory results are affected:
  - S2 D2 by a retained full concatenation;
  - every other S2/S3 cell by a memory-map load transient that coexisted with the copy;
  - S4 cells timed more than once by the previous repetition's outputs staying alive, and S4 cells whose resident
    inputs were comparable to their compute transients by the load transient.
- **Unaffected and reused:** 10 long batch S4 cells (single repetition, compute-dominated) — A: CW1, CW2; B2: fit
  (128/1015/20), (192/1015/20), (256/1015/20), (256/500/20), (64/1015/80), CW1, CW2, CW4.

## §2. Corrections

1. **S3 D2 size:** 1,620 is replaced by **1,600** (≤ 1,602 available). The D2 S3 sizes are {405, 810, 1600}.
2. **Input construction:** every operational cell reads only the rows it uses directly into one preallocated owned
   array (`eb1a2_lib.read_rows`, unbuffered `readinto`; D2 stored frames via `read_d2_stored`). There is no memory map
   and no intermediate full array. Requests beyond the available rows raise; none truncate silently. A constant
   `.npy`-header parsing transient (< 256 KB, size-independent) remains and is disclosed.
3. **Repetitions:** `measure_reps` releases (`del` + `gc.collect`) each output before the next repetition, so at most
   one output is alive. The memory predictor P is unchanged; it counts one retained output.

**What is reused:**
- **Time:** all Amendment 1 time measurements of S2, S3 (except the replaced size) and S4 are reused. The defects
  affect only memory lifetimes, not the timed computation (loading lies outside the timed region). The times taken
  in the re-measured cells are recorded and compared with the reused ones descriptively; they do not replace them.
- **Correctness:** a mandatory, separate subprocess checks that every new-loader input equals the Amendment 1 load
  path (dtype, shape, values) and is owned, and that an FC1 search through `measure_reps` equals the EB-1 harness.
  Earlier correctness results C0–C3 and C1′–C4′ are reused.
- **Memory:** the 10 unaffected cells in §1 are reused. **Every affected memory cell is re-measured** with the
  corrected construction: the S2/S3 cells at the Amendment 1 sizes (and the new D2 S3 size), and the affected S4
  fit/confirmation cells.

## §3. A batch-time model

**Amendment 1 failures** (fit set kept whole; no observation removed):
- n_Q = 32, n_D = 1015, n_win = 20: predicted 0.206 s, measured 0.155 s (absolute error +0.051 s, relative 0.33).
- n_Q = 64, n_D = 1015, n_win = 5: predicted 0.102 s, measured 0.080 s (absolute error +0.022 s, relative 0.27).

**Investigation (from Amendment 1 data only):**
- Per-comparison cost varies between about 0.20 and 0.26 µs across cells.
- The fitted per-database-row coefficient attributes a large share of time to per-call slicing, while the measured
  decomposition puts slicing at about 0.5 % of batch time for A.
- The three-parameter model with intercept therefore misallocates time between collinear work terms, which inflates
  the predictions for the shortest workloads.
- The timed region excludes setup, so a fixed intercept has no physical counterpart.

**Replacement model (frozen):** t = c1·C + c2·(call database rows) + c3·calls, **no intercept**.
- Fitted on the unchanged 11 Amendment 1 fit cells (Amendment 1 CPU times).
- If any coefficient is negative, that term is removed, in the order c3, then c2, and the model is refit.
- The criteria are unchanged: fit residual ≤ 0.25, and every confirmation workload within 25 %.
- Because this model was chosen **after** seeing the Amendment 1 residuals, its decisive confirmation uses **fresh**
  workloads, defined here before they are measured:

| ID | Data | n_Q | n_D | n_win |
|---|---|---|---|---|
| FC1 | D1 | 48 | 1015 | 10 |
| FC2 | D1 | 160 | 750 | 30 |
| FC3 | D1+ | 384 | 2693 | 60 |
| FC4 | D2 | 200 | 405 | 44 |

- The Amendment 1 confirmation errors of the replacement model (CW1, CW2, CW4) are also reported but do not decide.
- FC1–FC4 memory is added to the A batch memory confirmations.
- If support is still inadequate, A remains **insufficient evidence**.

## §4. Reporting

**Validity is reported separately:**
- **loader** (S1, S2; R3, R4-S1/S2);
- **A** (S3, S4_A; R1_A, R2_A, R4-S3, R4-S4_A);
- **B2** (S4_B2; R1_B2, R2_B2, R4-S4_B2).

Each group is valid only if its cells are complete, its correctness passes and its required projections are
supported.

**Global outcome.** The original §7 rule applies: any missing cell or unsupported required projection in any group
means **insufficient evidence**. Per-group status is reported beside it.

## §5. Budget and limits

- **CPU:** aggregate including EB-1 and Amendment 1 ≤ 2 CPU-h. The runner stops at a cumulative 1.8 CPU-h.
- **Memory:** `RLIMIT_AS` = 6 GB; MemAvailable ≥ 500 MB before each cell (wait ≤ 10 min, else skip → missing).
- Fresh subprocess per cell, BLAS threads 1.
- No desktop process is touched.
