# REPCHAR-1: SYNTHETIC Sensitivity Characterization of the Replica-Comparison Rule — Specification (2026-10-23)

**Status:** frozen with code by `benchmarks/audit/synth_rep_char/run_char.py freeze` before the deciding run.
It supersedes the uncommitted draft of the same date. The machine-readable cell grid is
`benchmarks/audit/synth_rep_char/run_spec.json`.

**Purpose:** on SYNTHETIC data with known truth, measure what the existing rule
`topofold_inputs.compare.replica_comparison` (as used by DEV-REP-1) detects, misses and confuses.

**Boundaries.**
- DEV-REP-1's recorded outcome (no window met the complete-separation criterion) is preserved and not re-scored.
- No real data and no real feature arrays are read.
- Alternative rules (the former P1/P2) are deferred to a separate decision.
- No temperature-contrast changes. Biological validation remains separate; ABL is preferred, other options retained.

## 1. Fixed analysis settings (as DEV-REP-1)

- Replicas: 5 vs 5 (S6 varies this). Frames: 387 per replica.
- Histograms: 36 bins, JS distance base 2, no pseudocounts.
- Linear channel: bins on [3, 8] with under/overflow.
- Circular channel: 36 arcs from origin −π.
- Categories exactly as implemented: *consistent*, *replica-dependent*, *inconclusive*.

## 2. Generators (`synthetic_generators.py`; tested in `tests/test_synthetic_generators.py`)

- **Latent AR(1):** Z₀ ~ N(0, 1) (stationary initialization); Z_t = φ Z_{t−1} + √(1 − φ²) ε_t, ε_t iid N(0, 1). The
  marginal is N(0, 1) at every t for every φ, so correlation never alters the marginal. φ ∈ {0, 0.8, 0.95}.
- **Linear channel:** X_t = μ + σ Z_t, with μ = 5.5 Å and σ = 0.2 Å (marginal N(μ, σ²)).
- **Circular channel:** θ_t = wrap(μ_c + s Z_t) to (−π, π], with μ_c = π and s = 0.3 rad. The marginal is wrapped
  normal with circular SD exactly s; samples straddle ±π.
- **Two-state mixture:** a hidden Markov chain S_t ∈ {0, 1}.
  - Transition probabilities: P(0→1) = w/D and P(1→0) = (1 − w)/D.
  - Stationary weight of state 1 is w. S₀ is drawn from the stationary distribution.
  - State autocorrelation is (1 − 1/D)^k. Mean dwell is D/w in state 0 and D/(1 − w) in state 1.
  - D = 1 gives independent states. D ≥ 1 and 0 < w < 1 keep both probabilities ≤ 1, so every listed (w, D) is
    achievable.
  - Emission: X_t = m_{S_t} + σ_c ε_t, ε_t iid N(0, 1).
- **Effective sample size:** n(1 − φ)/(1 + φ) ≈ 387, 43 and 10 for φ = 0, 0.8, 0.95. This is a **Gaussian
  sample-mean approximation** only. It is not a guarantee for histograms, JS distances or circular features.

## 3. Seeds and draw order

For each cell index c (its position in `run_spec.json`) and repetition r, the generator is
`numpy.random.default_rng(SeedSequence([20261023, c, r]))`. Draws use that generator in this fixed order:
1. condition-level offsets (S5: A then B);
2. replica-level offsets (S4b: A₀…A_{n−1}, then B₀…B_{n−1});
3. the series A₀…A_{n−1}, then B₀…B_{n−1}, each generated completely before the next.

## 4. Scenarios (84 cells × 200 repetitions)

All shifts are in units of σ (linear) or s (circular).

| ID | Scenario | Condition A | Condition B | Grid | Cells |
|---|---|---|---|---|---|
| S0 | Identical-distribution null | μ | μ | channel {linear, circular} × φ | 6 |
| S5 | **Condition-confounded sampling** | μ + s_A | μ + s_B; s_A, s_B iid N(0, (k σ)²), shared by all replicas of a condition | k ∈ {0.25, 0.5} × φ (linear) | 6 |
| S1 | Location shift | μ | μ + Δ | channel × Δ ∈ {0.1, 0.25, 0.5, 1, 2} × φ | 30 |
| S2 | Scale change | σ | ρσ | ρ ∈ {1.25, 1.5, 2} × φ (linear) | 9 |
| S3 | Population shift between fixed modes | w = 0.5 | w ∈ {0.4, 0.3, 0.1} | modes μ ∓ 2σ, σ_c = σ; × D ∈ {1, 20, 100} | 9 |
| S4a | One outlier replica | A₀ at μ + 2, others μ | μ + Δ | Δ ∈ {0, 1} × φ (linear) | 6 |
| S4b | Random replica offsets | μ + o_r | μ + Δ + o_r; o_r iid N(0, (k σ)²) | k ∈ {0.25, 0.5} × Δ ∈ {0, 1} × φ (linear) | 12 |
| S4c | Rare excursion in one replica | A₀: w = 0.1, modes μ and μ + 4σ (σ_c = σ), dwell D; other A iid N(μ, σ²) | iid N(μ, σ²) | D ∈ {1, 20, 100} | 3 |
| S6 | Replica count | μ | μ + 1 | n ∈ {3, 5, 8} per condition, φ = 0.8 (linear) | 3 |

**S5 versus S0.**
- In S0 both conditions have identical generating distributions. Any call is a false call about distributions.
- In S5 both conditions share one generating model, but each condition's replicas share a realized offset. The
  realized distributions **do differ**, by s_A − s_B (recorded per repetition), even though the difference is not
  caused by the condition.
- A call in S5 is therefore a correct description of the sampled data that would be a **confounded** attribution to
  the condition.
- **Analogy to mdCATH is limited.** mdCATH replicas share one starting configuration **across** temperatures. A
  shared component that is common **within** but not **between** conditions arises only if sampling history differs
  by condition, for example condition-specific slow relaxation. S5 measures that situation; it does not model
  mdCATH's actual initialization.

## 5. Outputs

Per cell:
- counts and frequencies of **each** category (*consistent*, *replica-dependent*, *inconclusive*), each with a
  Wilson 95 % interval;
- any combined "call" rate (*consistent* + *replica-dependent*), always shown with its breakdown;
- median and interquartile range of the between-above-within fraction;
- median of the median-between / median-within ratio, and of the pooled distance;
- for S5, the median realized |s_A − s_B| (in σ) among calls and among non-calls.

With 200 repetitions per cell, Monte Carlo precision is coarse (Wilson half-width up to about 0.07). Differences
between nearby cells smaller than that are not interpreted.

**Reading aid (not a criterion):** the report lists, per scenario, the smallest effect at which the combined call
frequency reaches 0.5, together with the breakdown.

## 6. Budget and partial-results procedure

- **Cap:** 900 CPU-s (15 CPU-minutes), enforced by a hard `RLIMIT_CPU`. Estimated cost is about 0.02 CPU-s per
  repetition, so about 6 CPU-minutes in total.
- **Cell order:** the order of `run_spec.json` (controls first: S0, S5, then S1, S2, S3, S4a, S4b, S4c, S6).
- **All-or-nothing cells:** a cell is reported only when all 200 repetitions have completed.
- **Results file:** rewritten after every completed cell.
- **Stopping rule:** before each cell, the runner estimates that cell's cost as 1.5 × the running mean CPU per
  repetition × 200 (× n/5 for S6). If the estimate would exceed the cap, the run stops. The results are then marked
  `partial` and list the cells not run. If the hard limit kills the process, the last written file stands, also
  marked partial.
- **No changes after outcomes:** settings, the grid and the reading aid are not changed after outcomes are seen.
