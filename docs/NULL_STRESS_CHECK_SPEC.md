# SYNTHETIC_NULLSTRESS-1: Null Stress Check, Different Temporal Correlation Between Conditions (2026-10-23)

**Status:** frozen with code by `benchmarks/audit/synth_null_stress/run_stress.py freeze` before the deciding run. Grid:
`benchmarks/audit/synth_null_stress/run_spec.json`.

**Question:** when both conditions have **identical stationary marginal distributions** but **different temporal
correlation**, how often does the PERM-1 replica permutation method reject, and how often does the complete-separation
screen call?
- The marginal distributions are the same; the conditions are not exchangeable as whole series, because replica
  histograms from more correlated series vary more.
- This informs the real-data design, where temperatures can change dynamics as well as distributions.

**Fixed settings.**
- **Generators:** the frozen `benchmarks/audit/synth_rep_char/synthetic_generators.py`, unchanged.
  - Linear: μ = 5.5 Å, σ = 0.2 Å.
  - Circular: wrapped normal around π, s = 0.3 rad.
  - Latent AR(1) with stationary initialization, so the marginal is identical at every φ.
- **Replicas and frames:** 5 per condition, 387 frames each.
- **Binning and distance:** 36 bins (linear [3, 8] with under/overflow; circular from −π), JS distance base 2.
- **Methods:** `replica_permutation.permutation_test` (all 252 allocations, tie rule as in PERM-1, α = 0.05; nominal
  size 12/252 ≈ 0.048) and `compare.replica_comparison` (categories reported separately).
- **Cells:** channel {linear, circular} × (φ_A, φ_B) ∈ {(0, 0.8), (0, 0.95), (0.8, 0.95)}, giving 6 cells.
- **Repetitions and seeds:** 200 repetitions per cell, seed `SeedSequence([20261024, cell, repetition])`. Draw order:
  A₀…A₄ (φ_A), then B₀…B₄ (φ_B).
- **Reporting:** rejection frequency with Wilson 95 % intervals (coarse at 200 repetitions), the three screen category
  frequencies separately, and the median within-A and within-B replica distance (descriptive).
- **No tuning:** nothing is changed after outcomes are seen.
- **Budget:** 120 CPU-s, hard `RLIMIT_CPU`. Before each cell, the run stops if used + 1.5 × the mean CPU per
  repetition × 200 would exceed the cap. Cells are all-or-nothing. Checkpoints are `partial` until all 6 cells finish.
- **Scope:** SYNTHETIC only; no real data. Historical artifacts are untouched.
