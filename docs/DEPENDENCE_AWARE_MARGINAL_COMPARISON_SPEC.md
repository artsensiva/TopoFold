# SYNTHETIC_DABOOT-1: Dependence-Aware Comparison of Marginal Distributions — Method Selection and Specification (2026-10-23)

**Status:** method selected and specification written before any characterization outcome. Frozen with code by
`benchmarks/audit/synth_daboot/run_daboot.py freeze`. Grid: `benchmarks/audit/synth_daboot/run_spec.json`.

**Target:** compare marginal conformational distributions between two conditions while allowing different temporal
correlation. This does not imply equilibrium sampling, identify a causal mechanism, or validate biological function.

## 1. Candidates (primary sources verified)

**(a) Stationary block bootstrap of each condition's empirical distribution.**
- Sources: Politis & Romano, *JASA* 89:1303–1313 (1994); automatic block length by Politis & White, *Econometric
  Reviews* 23:53–70 (2004), corrected by Patton, Politis & White, *Econometric Reviews* 28:372–375 (2009).
- Assumptions: stationarity and weak dependence within each series; block length large relative to the correlation
  length. Each condition is resampled separately, so the two conditions may differ in dependence.
- Fit to 5 × 387 frames per condition: uses all frames. At φ = 0.95 the theoretical optimal block length is about 53
  frames (corrected formula), giving about 7 blocks per replica and 35 per condition. Block bootstraps are known to
  under-represent variance when few blocks are available, so the method may be liberal in exactly the hardest null
  cells.
- Circular features: combined with Kuiper's two-sample statistic (Kuiper, *Proc. KNAW A* 63:38–47, 1960), which is
  invariant to the choice of origin.
- Cost: B resamples × 10 replicas, vectorised (measured below).
- Tunable parameters: expected block length (chosen automatically), number of resamples B, grid.

**(b) Replica-level (cluster) bootstrap.**
- Approach: resample whole replicas within each condition; the only dependence assumption is independence between
  replicas.
- Strength: robust to any within-replica dependence and to non-stationarity within replicas.
- Weakness: with 5 clusters per condition, a replica bootstrap has only 126 distinct resamples per condition and
  shrinks variance by (n − 1)/n. Cameron, Gelbach & Miller (*Rev. Econ. Stat.* 90:414–427, 2008) report that
  conventional cluster-based tests over-reject with few clusters (about 5–30). Their bootstrap-t refinements are
  reported to perform poorly below about 11 clusters in later work (Webb), which I have not verified from the primary
  source and do not rely on.
- Circular features: the same Kuiper statistic applies.
- Cost: negligible. Tunable parameters: B only.

## 2. Selection (before any outcome)

**Selected: (a), the centered stationary bootstrap with the Kuiper statistic.** It is the established method that
directly allows different dependence per condition while using within-replica information. Its five-replica weakness
is milder than (b)'s, whose few-cluster failure is documented for exactly this replica count.

Known risks, stated in advance:
- possible liberal behaviour when the block length is long compared with the replica (φ = 0.95);
- sensitivity to non-stationarity, which the synthetic generators exclude by construction;
- asymptotic, not exact, validity.

Not used: thinning or matching effective sample sizes; any new estimator.

## 3. Prototype (`python/topofold_inputs/dependent_bootstrap.py`)

- **Statistic:** V = max_g D(g) − min_g D(g), D = F_A − F_B. F_c is the ECDF of all frames of condition c (replicas
  pooled for the estimate only) on a fixed grid:
  - linear: 256 equally spaced points on [3, 8];
  - circular: 256 arc ends from −π + 2π/256 to π.
- **Null distribution:** for each of B = 499 draws and each replica, a stationary-bootstrap resample of that replica's
  own 387 frames. Geometric block lengths with mean L_c; a new block starts with probability 1/L_c at a uniform
  position; continuation wraps circularly **within the replica only**. Replicas are never joined.
  - D* = (F*_A − F_A) − (F*_B − F_B), V* = max D* − min D*, p = (1 + #{V* ≥ V − 10⁻¹²}) / (B + 1).
  - Reject at α = 0.05. The resolution of p is 1/500.
- **Block length L_c:** the median over the condition's replicas of the Politis–White stationary-bootstrap
  estimate with the 2009 correction (D̂_SB = 2ĝ(0)²).
  - Constants follow the defaults of the R package `blocklength`: K_N = max(5, ⌈ln n⌉); c = 1.96;
    autocorrelation threshold c·√(log₁₀ n / n); m̂ = smallest m with |ρ̂(m + k)| below threshold for k = 1..K_N;
    M = 2m̂; trapezoidal flat-top window; result clipped to [1, ⌈min(3√n, n/3)⌉] = [1, 60] for n = 387.
  - Circular: the larger of the estimates for cos θ and sin θ.
  - The closed form is tested against the published theoretical values (Patton et al. 2009, Table 1).

## 4. Cells, draws and seeds

- **Null cells 0–5:** the six SYNTHETIC_NULLSTRESS-1 cells (identical marginals, φ pairs (0, 0.8), (0, 0.95),
  (0.8, 0.95), linear and circular), with **identical draws**: generator `SeedSequence([20261024, cell, rep])`,
  series A₀…A₄ then B₀…B₄. The run checks this against the NULLSTRESS-1 histograms for repetition 0 of every null
  cell.
- **Positive controls 6–9:** A at φ = 0 and location μ; B at φ = 0.95 and location μ + Δ·scale (σ = 0.2 Å linear,
  s = 0.3 rad circular), with Δ ∈ {1, 2}, linear and circular. Same seed scheme with the new cell numbers.
- **Bootstrap randomness:** a separate stream, `SeedSequence([20261025, cell, rep])`.
- **Dimensions:** 5 replicas per condition, 387 frames, 200 repetitions per cell.

## 5. Calibration and sensitivity criteria (practical; fixed before the run)

- **Calibration failure:** any null cell with observed rejection > 0.10, or any null cell whose Wilson 95 % lower
  bound exceeds 0.05.
- **Calibration acceptable (practical):** neither condition holds. An interval containing 0.05 is **not** taken as
  proof of validity; it only fails to show miscalibration at 200 repetitions. Strongly conservative cells (rejection
  ≤ 0.01) are reported as such.
- **Sensitivity useful:** rejection ≥ 0.80 in both Δ = 2 positive controls. The Δ = 1 cells are reported.
- **Decision:**
  - calibration failure → report the failure and stop (no method switch, no retuning); recommend retaining
    descriptive comparisons and record the limitation;
  - acceptable calibration with useful sensitivity → recommend a separately frozen real-data development evaluation;
  - acceptable calibration with low sensitivity → record a specific unresolved limitation.
- **Comparison:** NULLSTRESS-1 permutation rejections on the identical null draws are shown alongside for context.
  They are not a criterion.

## 6. Budget and partial results

- **Timing check before freezing** (engineering, synthetic data only, no outcomes recorded): about 0.13 CPU-s per
  test at B = 499, so about 280 CPU-s for 2,000 tests.
- **Cap:** 900 CPU-s, hard `RLIMIT_CPU`. Before each cell, the run stops if used + 1.5 × the mean CPU per repetition
  × 200 would exceed the cap. Cells are all-or-nothing.
- **Checkpoints:** `partial` until all 10 cells finish. A partial result is reported as such, without re-running at
  lower precision.
