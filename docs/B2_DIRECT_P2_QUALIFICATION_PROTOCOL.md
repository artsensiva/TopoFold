# BQ-3 Direct P2 Qualification of Unchanged V2b: Protocol (frozen 2026-10-13)

**Starting checkpoint:** `f55def1a7c70a116e97fd0b6179fd2d249db5446` (`scientific-rebuild`).

**Basis:** the 2026-10-13 preflight stopped for lack of memory
(`docs/B2_DIRECT_P2_QUALIFICATION_PREFLIGHT.md`; record `benchmarks/audit/bq3/results/bq3_memory_preflight.json`,
kept unchanged). A recheck the same day found about 5.85 GB available against the documented 3.16 GB requirement.

**Scope:** this is a **separate, direct measurement**.
- It does not fit any model, change BQ-2's conclusion C, or require BQ-2's unsupported preparation models to become
  supported.
- BQ-1, BQ-2 and EB-1 records stay as they are.
- Out of scope: harness changes to V2b, loader work, compiled ports, production changes, downloads, held-out scoring,
  pushing, tags and releases.

## §1. Candidate and workload (unchanged)

**Candidate:** V2b from `benchmarks/audit/b2opt/b2opt_lib.py`, at the hash in BQ-1's pre-measurement record. That
fixes the solver, both fallbacks, the proper-rotation RMSD objective, exhaustive search and the `topk_lex` /
`rank_check` policy.

**Workload:** the BQ-2 `P2D` definition, `bq2_workloads.json`: n_D = 10,000 database rows, L = 1,000 residues,
W = 8, 993 windows, three query rows.
- The rows are built by the BQ-2 loader (`bq2_run.load`, `bq2_lib.tiled_rows`, hash-checked).
- **SYNTHETIC_tiled, performance-only:** real D1+ frames are joined end to end and repeat after 3,706 rows.
- **Labelling:** every result file name and report title carries SYNTHETIC.
- **Limitations:**
  - These rows test computational scale at P2 dimensions only. They establish nothing about real 1,000-residue
    proteins or real 10,000-frame trajectories (data distribution, fallback frequency, cache behaviour), and nothing
    about retrieval quality or biology.
  - Repeated rows create exact duplicates. RMSD ≈ 0 then sits at the square-root rounding floor and creates many near
    ties (seen in BQ-2).

## §2. Processes (sequential, each a fresh subprocess)

1. **Tests** (`test_bq3.py`): unit tests of the acceptance function's failure paths. All must pass before anything
   else runs.
2. **Reference** (`ref_svd`, outside any timing): distances for all three queries, all windows and all rows. Saved
   under `build/`, and checked to be finite.
3. **Gate** (outside any timing): V2b prepares the database and searches all three queries, keeping distances.
   - Checks: candidate distances finite; max |d − d_ref| ≤ 10⁻⁶ Å; zero rank-non-equivalent query-windows.
   - Reported, not deciding: strict-identity counts, near-tie counts and fallback rows.
   - Records the SHA-256 of each query's top-k array.
4. **Timing repetitions:** R = **5** fresh processes. Each runs, in order:
   1. the identical SYNTHETIC warm-up (BQ-1/BQ-2);
   2. building the P2 rows (not timed as preparation);
   3. **initial preparation (P1)**, then **query 1 (Q1)** immediately after. This is the actual first-query sequence,
      also timed as one span (S = P1 + Q1, measured end to end);
   4. **subsequent queries Q2, Q3** on the same prepared database, each output released before the next;
   5. releasing the prepared arrays and outputs, then **repeated preparations P2, P3**, each released before the
      next.

   Per process, the following are recorded:
   - wall and CPU time of every step;
   - lifetime peak RSS after loading, after P1, after the queries, and at the end;
   - fallback rows;
   - top-k hashes of Q1–Q3;
   - threadpool information (BLAS/OpenMP limited to 1) and the OS thread count;
   - `/proc/loadavg`, MemAvailable and SwapFree at start and end.

**Memory start rule (unchanged from BQ-2):** a process starts only if MemAvailable − predicted peak ≥ 600 MB, with
a wait of at most 10 min. Otherwise it is skipped and recorded as missing. No application is closed and no system
setting is changed.

**CPU cap:** ≤ 1.0 CPU-h for BQ-3 in total: tests, reference, gate and repetitions, parent and children. The runner
stops before any step once the total exceeds 0.9 CPU-h.

## §3. Acceptance (frozen; thresholds unchanged)

**Limits:** R2 = 10 s wall. R4 = 8 GiB peak RSS (8·1024² kB).

**Verification (mandatory for any qualification).** All of the following must hold:
- the reference exists and is finite;
- the gate exists, completed, and passed: candidate finite, distances within 10⁻⁶ Å, zero non-equivalent rankings;
- every timed repetition's Q1, Q2 and Q3 top-k hashes equal the gate hashes for the same queries.

A missing, failed or mismatched check gives **not qualified (verification)**.

**Completeness.** All 5 repetitions must complete. Otherwise the result is **not qualified (insufficient
evidence)**.

**Pass/fail aggregation rule (all repetitions must pass):**
- **Prepared-database operation qualifies** iff verification and completeness hold, every one of the 15 query
  latencies (Q1–Q3 × 5) is ≤ 10 s, and every repetition's peak RSS is ≤ 8 GiB.
- **First-query operation qualifies** iff verification and completeness hold, every repetition's measured first-query
  span S = P1 + Q1 is ≤ 10 s, and every repetition's peak RSS is ≤ 8 GiB.
- The two verdicts are reported separately. Medians and ranges are reported but do not decide.

**R1 is not measured here.** R1_B2 remains the EB-1 Amendment 2 batch projection. The directly measured P2
preparation CPU may be shown beside it, labelled as a mixture of one projection and one direct measurement. It does
not decide.

**Reporting:**
- every repetition and every step;
- direct measurements kept separate from any retained projection;
- the SYNTHETIC/performance-only labelling and its limitations (§1).

**Consequences:**
- If prepared-database operation qualifies, a bounded integration plan may be **proposed** for review. It is not
  implemented.
- If either verdict fails or cannot be reached, the specific limitation is reported. No new experiment starts
  automatically.
- The results are not pushed automatically.
