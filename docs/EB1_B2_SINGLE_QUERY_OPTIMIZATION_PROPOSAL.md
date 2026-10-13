# Proposal: bounded B2 single-query optimization experiment (not implemented)

**Status:** proposal only, 2026-10-12. Nothing here is implemented, measured or frozen. It follows the EB-1 Amendment 2
outcome (`benchmarks/audit/eb1/REPORT_A2.md`). There, the only violated requirement is R2_B2: the projected
single-query B2 wall time at P2 exceeds 10 s. Values are in `results/eb1a2_analysis.json`. The Amendment 1 time
decomposition (`REPORT_A1.md` §5) attributes most B2 time to the RMSD kernel, not to harness slicing or top-k.
This proposal does not claim that any saving exists.

## 1. Question

Can single-query exhaustive B2 search be made to meet R2 (≤ 10 s at P2: L = 1000, W = 8, n_win = 993,
n_D = 10,000), with the same results as the current harness?

## 2. Fixed semantics (must not change)

- **Distance:** optimal-superposition Cα RMSD over the W window atoms, from the Horn quaternion (largest eigenvalue).
  Units are Å and the arithmetic is float64.
- **Ranking:** for every query and window, the top k = 10 database frames by ascending RMSD, using the current
  harness tie rule (`eb1_lib.topk`). Query and database sets, windows and k are unchanged.
- **Equal information:** the same inputs as the harness. There are no precomputed query-specific pruning hints from
  outside the database.

## 3. Candidate variants (one at a time, each a separate arm)

- **V1:** precompute per-window centred coordinates and squared norms of the database once (reused across queries),
  instead of re-slicing and re-centring per call.
- **V2:** batch the 4×4 key-matrix construction across database rows in one vectorised contraction, then compute
  the largest eigenvalue in closed form or by Newton iteration on the characteristic polynomial (as in QCP). Fall back
  to the eigen-solver where convergence fails.
- **V3:** V1 + V2 in a compiled kernel (the existing Rust core or a NumPy-only path; no new dependency without
  approval).

Exact pruning (for example, lower bounds from centroid-free invariants) is **out of scope** for this experiment. It
changes the work performed and needs its own correctness argument.

## 4. Correctness checks (mandatory; any failure stops the arm)

1. **Distances:** for every variant, compared with the harness Horn RMSD and an independent Kabsch (SVD) RMSD on all
   D1, D1+ and D2 confirmation workloads. Requirement: max |Δ| ≤ 1e-6 Å. Edge cases are also tested: identical
   frames (RMSD 0), reflected or collinear windows, and near-degenerate eigenvalues.
2. **Top-k:** identical indices to the harness for every query and window. Differences are allowed only where the
   harness k-th and (k+1)-th distances lie within 1e-9 Å (near ties), and these are counted and reported.
3. **Invariance:** RMSD is unchanged under a random rigid transform of the query and under permutation of the
   database order (after index mapping).

## 5. Representative workloads

- **Fit:** single query (n_Q = 1) on D1 at n_D ∈ {250, 500, 1015} × n_win ∈ {20, 207}, the same grid as EB-1
  Amendment 1, so arms are comparable.
- **Fresh confirmation, frozen before measurement:** D1+ (n_D = 2693, all windows), D2 (n_D = 405, all windows), and one
  new D1 grid point.
- **Per arm:** the harness baseline is re-measured in the same session, the same process type and BLAS threads 1, with
  fresh subprocesses and released outputs (`eb1a2_lib.measure_reps`). Inputs are read directly (`read_rows`).

## 6. Measured acceptance criteria (frozen before measurement)

An arm is accepted only if all of the following hold:
- all checks in §4 pass;
- its single-time model (c0 + c1·C + c2·n_win) is supported under the unchanged EB-1 rules: ≥ 3 sizes, c1 > 0, fit
  residual ≤ 0.25, every fresh confirmation within 25 %;
- the projected R2_B2 at P2 is ≤ 10 s, reported with its extrapolation factor;
- the measured speed-up over the same-session harness is reported on every confirmation workload as the median of
  repetitions with min–max. No speed-up is claimed outside the measured workloads;
- peak RSS stays within R4 under the Amendment 1/2 predictor, including any new precomputed arrays.

If no arm meets every criterion, the outcome is "R2_B2 not met by the tested variants". The result is not to be
described as "B2 is too slow" in general.

## 7. Budget and scope

- At most 1 CPU-h, on development data only (D1, D1+, D2 caches). No downloads, no held-out data and no production
  code change until an arm is accepted and a separate change is approved.
- Each arm is a harness-level implementation in `benchmarks/audit/`. Any port to `crates/` is a later,
  separately reviewed step with its own tests (`cargo test --workspace`).
