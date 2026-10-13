# BQ-2 Scaling Confirmation of Unchanged V2b: Protocol (frozen 2026-10-12)

**Basis and disclosure:** BQ-1 (`docs/B2_SINGLE_QUERY_OPTIMIZATION_PROTOCOL.md`, Amendment 1,
`benchmarks/audit/b2opt/REPORT.md`) is a completed, bounded experiment.
- Its protocol, amendment, code, measurements and unsuccessful overall acceptance outcome are preserved unchanged.
- BQ-2 is **informed by BQ-1** (its timings, models and the preparation-fit failure) and is not an independent
  rediscovery. V2b is the leading candidate; production readiness is not established.
- This is a harness study on development data plus labelled SYNTHETIC stress data.
- Out of scope: loader work, compiled ports, production changes, downloads, held-out scoring, staging, commits and
  tags.

## §1. What BQ-1 established, and what it did not

The three categories are kept distinct:

| Category | BQ-1 status |
|---|---|
| Measured correctness failure | V2 (frozen) on SYNTHETIC degenerate windows. None for V0, V1 or V2b. |
| Measured resource-limit violation | None. No R1, R2, R3 or R4 limit was exceeded by any measurement. |
| Unsupported projection / insufficient evidence | V1 and V2b preparation models (fit residual > 0.25). R1_B2 + preparation therefore **was not evaluated**; it was not violated. |
| Projected violation (supported model) | V1 projected prepared-query latency > 10 s. |

BQ-2 does not change any BQ-1 verdict. BQ-1 cells are used only as disclosed diagnostic comparisons.

## §2. Unchanged candidate

V2b is `benchmarks/audit/b2opt/b2opt_lib.py` at its BQ-1 pre-measurement hash. Unchanged and imported, not copied:
- the solver (QCP Newton with frozen parameters);
- both fallbacks;
- the proper-rotation RMSD objective;
- exhaustive search;
- the `topk_lex` ranking policy and the `rank_check` equivalence rule;
- the tolerances (10⁻⁶ Å distance, 2·10⁻⁶ Å rank).

The baseline V0 is the EB-1 harness kernel, as in BQ-1.

## §3. Preparation scaling (structural)

`prepare()` loops over windows. In each iteration it makes a fancy-index copy of (n_D, W, 3), computes a mean,
subtracts it, writes into the prepared array and computes a norm. The cost therefore has:
- a per-row-per-window term proportional to n_D·n_win;
- a per-window overhead proportional to n_win, from Python and NumPy call setup in each iteration;
- a constant.

BQ-1's model omitted the per-window term. The BQ-1 fit grid had n_win only ∈ {20, 207}, so the per-window term
could not be separated.

**Frozen model:** t_prep = c0 + c1·(n_D·n_win) + c2·n_win.
- This is the same functional form, code and support rules as the single-query time model:
  `eb1a1_analyse.s4_model(kind='single_time')`. c2 is dropped if negative.
- It is fitted separately for preparation wall time and preparation CPU time.
- The fit grid varies n_win at fixed n_D and n_D at fixed n_win (§4).
- No threshold is changed.

## §4. Workloads (`results/bq2_workloads.json`)

**Real query frames:** drawn with seed 20261018 from query positions never used by EB-1, its amendments or BQ-1.

| Set | Data | Configurations | Role |
|---|---|---|---|
| Fit | D1+ (real, L = 214) | n_D ∈ {100, 700, 2693} × n_win ∈ {5, 20, 60, 207}; 3 queries | model fitting only |
| RC1–RC3 | D1+ 1500 × 120; D1+ 2693 × 207; D2 405 × 44 (real, L = 51) | 3 fresh queries each | deciding confirmation, real data |
| S1, S2 | SYNTHETIC_tiled L = 1000: n_D 1000 and 2500 × n_win 993 | 3 tiled queries | deciding confirmation, performance-only |
| S4_P2slice | SYNTHETIC_tiled L = 1000: **n_D = 10,000 (P2) × n_win 100** | 3 tiled queries | deciding confirmation, performance-only; direct per-window cost at the P2 database size |
| S3 | SYNTHETIC_tiled L = 1000: n_D 4000 × n_win 993 | 3 tiled queries | deciding if memory permits (§6); otherwise reported missing |
| P2D | SYNTHETIC_tiled L = 1000: n_D 10,000 × n_win 993 (= P2) | 1 query, V2b only | direct P2 check, only if memory permits |

**SYNTHETIC_tiled (performance-only).** Each row is ⌈1000/214⌉ = 5 real D1+ frames joined end to end, with a 3.8 Å
gap along +x, truncated to 1000 residues. Database row r uses frames (5r + j) mod 3706; query rows use r = 100,000 + q.
- **What it tests:** scaling of preparation, query time and memory with n_D·n_win and n_win at L = 1000 and up to the
  P2 window count (993), and with the P2 database size (10,000).
  - Most windows contain real local Cα geometry. That keeps solver convergence and fallback behaviour close to real
    protein windows; fallback counts are reported.
- **What it cannot establish:**
  - Rows repeat after 3706, and junction windows (28 of 993) are artificial.
  - It cannot establish performance on real L ≈ 1000 proteins or real 10,000-frame trajectories (data
    distribution, cache behaviour, fallback frequency).
  - It cannot establish any retrieval quality or biological property.
  - It is not a substitute for biological validation.

**Diagnostics, not deciding:** the BQ-1 V2b preparation cells, including the smallest (250 × 20), and its query cells
are evaluated under the BQ-2 models.

## §5. Measurements

- **Process:** each (round, workload, arm) runs in a fresh subprocess, one thread (BLAS/OpenMP = 1), `RLIMIT_AS` 6 GB.
- **Inputs:** real inputs are read directly into owned arrays (`eb1a2_lib.read_rows`). Tiled inputs are built
  in-process from the D1+ array before timing.
- **Warm-up:** an identical SYNTHETIC warm-up (as BQ-1).
- **Preparation repetitions:** V2b preparation is repeated r_p times. r_p = 9 if the BQ-1-predicted preparation is
  < 0.05 s, 5 if < 1 s, otherwise 3. Each prepared array is released (`del` + `gc.collect`) before the next
  repetition.
  - Preparation time = median of the repetitions. All repetitions and their dispersion are recorded.
- **Queries:** the last preparation is kept and three queries are timed on it, each output released before the next.
  V0 runs the same three queries with no preparation.
- **Rounds:** 3; P2D has 1. Workloads are shuffled with seed 20261019 and arms are shuffled within each workload. The
  sequence number and UTC start time are recorded.
- **Contention:** at start and end, record `/proc/loadavg`, MemAvailable and query CPU/wall ratios.
- **Memory:** lifetime peak RSS after loading, after preparation and at the end, plus exact bytes of resident,
  prepared and retained arrays.
  - Memory predictor P = resident + prepared + retained + α·n_D, with α = BQ-1 V2b gate calibration (reused,
    recorded).
- **Latencies:**
  - Prepared-query latency = median of the 3 queries.
  - First-query latency = median preparation + query 1.
  - Amortized latency (N) = preparation / N + prepared-query latency.

## §6. Memory headroom and stop rules

- **Start rule:** a cell starts only if MemAvailable − predicted peak ≥ **600 MB**.
  - Predicted peak = 60 MB baseline + 64 MB transient + resident + prepared bytes.
  - Otherwise wait up to 10 min (15 s polls). If the rule is still not met, the cell is **skipped** and recorded as
    missing for memory.
- **No system changes:** no application is closed and no system setting is changed.
- **CPU:** ≤ **1.0 CPU-h** for BQ-2 in total: tests, references, gates and timing, parent and children. Stop before
  any step once the total exceeds 0.9 CPU-h. Steps not run are reported missing.
- **Order:** tests → references → V2b gate → V0 gate → timing rounds → P2D (if allowed) → analysis.
- **Gate failure:** a V2b gate failure (measured correctness failure) stops timing.

## §7. Correctness gate (before timing)

The gate covers every workload except P2D. P2D is checked on its single timed query against a reference computed
only if the start rule allows.
- **Reference:** `b2opt_lib.ref_svd`.
- **Requirements:** for V2b and V0, max |Δd| ≤ 10⁻⁶ Å and zero rank-non-equivalent query-windows (§1.2 of BQ-1).
  Strict identity, near ties and fallback rows are reported.
- **Output identity:** timed outputs must hash-equal the gate outputs.

## §8. Analysis and conclusions (frozen)

**Models (V2b).** All are fitted on the fit set only.
- Query time c0 + c1·C + c2·n_win (C = n_D·n_win) and memory m0 + m1·P use `s4_model`.
- Preparation wall and CPU use the §3 model.
- **Support:** ≥ 3 distinct sizes, leading coefficient > 0, max fit relative residual ≤ 0.25.
- **Confirmation:** every deciding confirmation (RC1–RC3, S1, S2, S4_P2slice, and S3 if measured) within 25 %.
- **Required cells:** RC1–RC3, S1, S2 and S4_P2slice must be measured. If any is missing, the outcome is
  insufficient evidence.

**P2 values** use C = 9.93·10⁶ and n_win = 993.
- They are reported with extrapolation factors relative to the largest measured C.
- Supplementary: a slice-scaled estimate = S4_P2slice prepared-query time × 993/100, which assumes only linearity in
  windows at fixed n_D.
- If P2D is measured, its direct values are reported beside the projections.

**Conclusions.**
- **A — unchanged V2b supports prepared-database operation:** all of the following hold:
  - gates pass and timed outputs are identical to the gate outputs;
  - the query, preparation-wall, preparation-CPU and memory models are supported and confirmed;
  - projected P2 prepared-query latency is ≤ 10 s;
  - projected P2 peak RSS is ≤ 8 GB;
  - R1_B2 (Amendment 2) + projected P2 preparation CPU is ≤ 96 CPU-h;
  - and, if P2D is measured, its direct prepared-query latency is ≤ 10 s.
- **B — first-query latency still warrants optimization:** A holds, and projected (or directly measured) P2
  first-query latency is > 10 s.
- **C — evidence remains insufficient:** any required model is unsupported or unconfirmed, or a required cell is
  missing. The specific missing evidence is listed. No further amendment starts automatically.
- **Other categories, reported if they occur:**
  - measured correctness failure;
  - measured resource-limit violation;
  - projected violation with supported models (A not met although the evidence is sufficient).

**Reporting.** Direct measurements per workload are reported separately from fitted projections. Speed-ups versus
V0 are per-round paired ratios.
