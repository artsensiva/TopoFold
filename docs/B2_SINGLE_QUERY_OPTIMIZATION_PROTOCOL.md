# B2 Single-Query Optimization Experiment (BQ-1): Protocol (frozen 2026-10-12)

**Basis:** EB-1 Amendment 2 is accepted for bounded engineering prioritization. Its only violated requirement is
R2_B2 (`benchmarks/audit/eb1/results/eb1a2_analysis.json`). This protocol refines
`docs/EB1_B2_SINGLE_QUERY_OPTIMIZATION_PROPOSAL.md`; the proposal is unchanged.

**Scope:** EB-1 (all versions), its results and outcomes are preserved. This is a harness-only experiment in
`benchmarks/audit/b2opt/` on development data (D1, D1+, D2 caches).
- Out of scope: production code, downloads, held-out scoring, staging, commits and tags.
- Any production integration needs a separate review.

**Freeze:** this document and `results/b2opt_workloads.json` are frozen before any implementation. Code hashes are
recorded before any timing (`results/b2opt_pre_measurement_record.json`). Thresholds, tolerances, solver parameters
and the workload split do not change after results are seen.

## §1. Objective and fixed semantics

**Objective (unchanged):** minimal Cα RMSD over **proper rotations** between the W = 8 window of one query frame and
the same window of every database frame. Centroids are removed. Coordinates are the S2 rounded float64 arrays. Units
are Å and all computation is float64. Search is **exhaustive**: every database frame and every requested window, with
no pruning or approximation.

### §1.1 Top-k and tie policy (one policy for every arm, the baseline and the reference)

- **Canonical order:** ascending by (d, j), where d is the arm's own float64 distance and j the database position.
- **Top-k:** the first k = 10 entries of that order. Exactly equal distances are broken by the lower position.
- **Implementation:** `topk_lex`, which selects all candidates with d ≤ the k-th partition value and then applies
  `lexsort`.
- **Harness note:** the EB-1 harness `topk` (argpartition plus a stable sort of the selected k) does not fix which
  item is kept on an exact tie at the k-th boundary. The baseline in this experiment is therefore the harness distance
  kernel with `topk_lex`; the cost of `topk_lex` is in every arm.

### §1.2 Ranking equivalence (tested separately from distance agreement)

Let d_ref be the reference distances, s_r the r-th smallest reference distance, and j_r the arm's r-th item.
- **Tolerance-equivalent:** for every rank r = 1…k, |d_ref(j_r) − s_r| ≤ τ_rank, and the k items are distinct.
  τ_rank = 2·10⁻⁶ Å, twice the distance tolerance.
- **Strictly identical:** j_r equals the reference canonical top-k at every rank.
- **Near-tie query-window:** some adjacent reference gap among ranks 1…k+1 is ≤ τ_rank.

The **gate** requires zero non-equivalent query-windows. Strict-identity and near-tie counts are reported.

### §1.3 Reference

**R_svd:** batched Kabsch. SVD of the 3×3 cross-covariance, with the determinant sign correction for a proper
rotation; RMSD² = (G_A + G_B − 2(σ1 + σ2 + sign·σ3)) / W, clamped at 0. This is independent of the Horn eigenvalue
path. It is itself checked against two other references:
- `wpd_lib.kabsch_rmsd`: per pair, explicit rotation and residuals;
- MDAnalysis `lib.qcprot.CalcRMSDRotationalMatrix` (another library; inputs centred first).

The check uses a fixed sample (seed 20261014) of 2,000 pairs per workload plus every SYNTHETIC case. Agreement must
be ≤ 10⁻⁶ Å. If the references disagree, the experiment stops.

## §2. Arms (predeclared; executed in this order)

| Arm | Database preparation | Per-query, per-window computation |
|---|---|---|
| **V0 (baseline)** | none | EB-1 harness `wpd_lib.horn_rmsd`: centres query and database windows, builds the 4×4 key matrix, `eigvalsh` |
| **V1 (simplest)** | per window: centred database coordinates (n_win, n_D, W, 3) and G_B (n_win, n_D), float64 | centre the query window; cross-covariance by `einsum` on the prepared array; key matrix; `eigvalsh` (as V0) |
| **V2** | as V1 | as V1, but the largest key-matrix eigenvalue comes from Newton–Raphson on its characteristic polynomial (QCP), with an `eigvalsh` fallback |

**V2 solver (references verified 2026-10-12):**
- Sources:
  - Horn (1987), J. Opt. Soc. Am. A 4(4):629–642;
  - Theobald (2005), Acta Cryst. A 61:478–480, doi:10.1107/S0108767305015266;
  - Liu, Agrafiotis & Theobald (2010), J. Comput. Chem. 31(7):1561–1563, doi:10.1002/jcc.21439.
- Basis from the sources: the RMSD needs only the largest eigenvalue λ_max of the traceless symmetric key matrix K,
  and Newton–Raphson starts from the upper bound E0 = (G_A + G_B)/2.
- The coefficient expressions were not quoted from the primary text (not accessible here). They are **derived** from
  Newton's identities for a traceless 4×4 matrix: P(λ) = λ⁴ + c2·λ² + c1·λ + c0, with c2 = −½·tr(K²),
  c1 = −⅓·tr(K³) and c0 = det K.
- Regression tests check these coefficients against `numpy.poly(K)` before use. The 4×4 determinant uses the
  2×2-minor expansion, tested against `numpy.linalg.det`.
- **Parameters, fixed a priori and not tuned:** start λ = E0; stop when |Δλ| ≤ 10⁻¹¹·max(λ, 1 Å²); at most 50
  iterations.
- **Fallback to `eigvalsh` for a row** if it does not converge, if λ is not finite, or if λ is outside
  [0, E0·(1 + 10⁻¹²)]. Fallback rows are counted and reported per workload.

**Out of scope:** pruning or lower bounds, float32, compiled ports and multithreading. All runs use one thread.

## §3. Workloads (`results/b2opt_workloads.json`)

**Reused comparisons** (used in earlier EB-1 versions; the query frames are the first query positions 0, 1, 2):
- **Fit grid (model fitting):** D1 single query, n_D ∈ {250, 500, 1015} × n_win ∈ {20, 207}.
- **Reported, not deciding:** CW3 (D1+, n_D 2693, 207 windows) and CW5 (D2, n_D 405, 44 windows).

**Fresh confirmation (reserved):** new query frames drawn with seed 20261013 from query positions not used in any
earlier EB-1, Amendment 1 or Amendment 2 sample, batch prefix or single-query cell. Configurations are new.

| ID | Data | n_D | n_win |
|---|---|---|---|
| FQ1 | D1 | 750 | 120 |
| FQ2 | D1+ | 2000 | 207 |
| FQ3 | D2 | 405 | 44 |
| FQ4 | D1+ | 2693 | 150 |

These are used only for correctness gates and for confirmation. They are not used for tuning (there are no tuned
parameters), model choice or the fit.

## §4. Measurements

- **Process:** each (round, workload, arm) runs in a fresh subprocess.
  - Threads: `OMP/OPENBLAS/MKL/RAYON_NUM_THREADS = 1`; threadpoolctl information is recorded.
  - Limits: `RLIMIT_AS` 6 GB; MemAvailable ≥ 500 MB before each start.
- **Inputs:** read directly into owned arrays (`eb1a2_lib.read_rows`): the query rows plus the first n_D database
  rows.
- **Warm-up:** an identical SYNTHETIC warm-up (one tiny search) in every arm before timing.
- **Timed sequence:**
  - V1, V2: **preparation** (wall and CPU) is timed, then three queries on the prepared database, each timed
    separately. Each output is released before the next query.
  - V0: the same three queries, with no preparation.
- **Latency definitions:**
  - **First-query (cold) latency** = preparation + query 1.
  - **Prepared-database query latency** = median of the three query times.
  - **Amortized latency** over N queries = preparation / N + prepared-query latency.
- **Memory:** lifetime peak RSS (`ru_maxrss`) is recorded after loading, after preparation and at the end. Exact
  bytes are recorded for resident inputs, prepared arrays and retained outputs.
  - The memory predictor is P = resident + prepared + retained + α_arm·n_D.
  - α_arm is the per-window transient bytes per database row. It is calibrated by tracemalloc on SYNTHETIC data in
    the arm's gate, before timing.
- **Rounds:** 3. Within each round, workloads run in seeded-shuffled order and arms in seeded-shuffled order within
  each workload. The sequence number and start time of every subprocess are recorded.
- **Output identity:** each timed run's top-k arrays are hashed and must equal the arm's gate outputs.

## §5. Correctness gates (per arm, before any timing of that arm)

1. **Distances:** on every workload, query and window, max |d_arm − d_ref| ≤ 10⁻⁶ Å.
2. **Ranking:** zero non-equivalent query-windows (§1.2). Strict identity and near ties are reported.
3. **SYNTHETIC adversarial suite:**
   - exact duplicates in the database: bitwise-equal distances, ranked by position;
   - query identical to a database frame: d ≤ 10⁻⁶ Å, and that frame ranked first among the zero-distance group;
   - database permutation: per-item |Δd| ≤ 10⁻⁹ Å, and rankings equivalent after index mapping;
   - rigid transform of the query: |Δd| ≤ 10⁻⁸ Å;
   - collinear, planar and mirror-image windows: within 10⁻⁶ Å of the reference;
   - constructed near ties (perturbations about 10⁻⁹ Å): rankings equivalent.
4. **V2 fallback:** frequency reported. Fallback rows are included in all checks.

V0 must pass the gate. If it does not, the experiment stops: baseline and reference disagree. An arm that fails its
gate is **rejected** and is not timed.

## §6. Analysis and acceptance (thresholds as in EB-1; frozen)

**Models, per arm.** Rules as in EB-1 Amendment 1: ≥ 3 distinct sizes, leading coefficient > 0, max fit relative
residual ≤ 0.25, every fresh confirmation (FQ1–FQ4) within 25 %.
- Query time t = c0 + c1·C + c2·n_win, with C = n_D·n_win. c2 is dropped if negative (`eb1a1_analyse.s4_model`,
  kind single_time). It is fitted on the per-cell medians across rounds of the fit grid.
- Preparation wall time and preparation CPU time each = p0 + p1·(n_D·n_win) (`eb1_analyse.fit`).
- Memory = m0 + m1·P on end-of-process peak RSS (`s4_model`, kind mem).
- **P2 values:** C = 993 × 10,000.
  - Resident = 11,000 frames × L 1,000 × 3 × 8 B.
  - Prepared = n_D·n_win·(W·3 + 1)·8 B.
  - Retained = n_win·k·8 B.

**Connection to R2.** The original R2 is "P2 single-query, all windows: query S3 + S4 wall time ≤ 10 s, single
process". It counts query-side work; database-side stage outputs are resident. For A, the database κ and τ come from
S3 in batch and are not charged to R2.

BQ-1 treats the V1/V2 preparation the same way, as a **database-side stage**, and charges it explicitly:
- **Primary verdict, "R2 supported (prepared database)":** an arm qualifies only if **all** of the following hold:
  - (a) gates pass and the timed outputs match the gate outputs;
  - (b) the query-time model is supported and confirmed, and the projected P2 prepared-query latency is ≤ 10 s;
  - (c) the preparation wall and CPU models are supported and confirmed;
  - (d) the memory model is supported and confirmed, and the projected P2 peak RSS, prepared arrays included, is
    ≤ 8 GB;
  - (e) R1_B2 from Amendment 2 plus the projected P2 preparation CPU is ≤ 96 CPU-h.
- **Strict verdict, "R2 supported cold":** the conditions above, and the projected P2 first-query latency
  (preparation + query) is ≤ 10 s.

**Reported for every arm, not deciding:**
- measured speed-ups versus V0 per workload, as per-round paired ratios (median and min–max) for prepared-query and
  cold latency;
- preparation time and bytes;
- amortized latency at N ∈ {1, 10, 100, 1000};
- break-even N for 10 s;
- confirmation errors;
- fallback frequency;
- CW3/CW5 errors;
- V0 same-session projection next to the Amendment 2 R2_B2 value.

**Measured vs projected.** Direct measurements on the development workloads are reported separately from fitted P2
projections. Every P2 value states its extrapolation factor.

If no arm qualifies, the outcome is "R2_B2 not supported by V1 or V2 on these workloads". It is not generalized to
other B2 implementations.

## §7. Budget and stopping rules

- **CPU:** ≤ 1.0 CPU-h in total for this experiment. That includes tests, reference checks, gates, calibration,
  preparation and timing: the parent process plus all children.
- **Stop:** the runner stops before any step once the cumulative total exceeds 0.9 CPU-h. An arm whose gate or timing
  is not complete is reported as **untested** (or incomplete) and cannot qualify.
- **Order:** tests → reference check → V0 gate → V1 gate → V1 timing (with V0) → V2 gate → V2 timing (with V0) →
  analysis.
- A failed V0 gate or reference check stops everything. A failed arm gate skips that arm's timing.
- No threshold, tolerance, solver parameter or workload is changed after measurement starts.
