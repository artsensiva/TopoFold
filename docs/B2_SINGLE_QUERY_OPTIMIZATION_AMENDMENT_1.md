# BQ-1 Amendment 1: V2b Arm (recorded 2026-10-12, before any measurement)

**Amends:** `docs/B2_SINGLE_QUERY_OPTIMIZATION_PROTOCOL.md`. That file is unchanged, and so are its thresholds,
tolerances, workloads, budget and acceptance rules.

**Trigger:** this was found before any real-data computation or timing, during the SYNTHETIC regression tests
written for the protocol. Frozen V2 exceeds the 10⁻⁶ Å distance tolerance on identical **collinear** windows; the
error is about 3·10⁻⁴ Å (`test_b2opt.py`).
- When either window is collinear, the optimal rotation is not unique and λ_max is a multiple root of P(λ).
- Near a double root, the rounding noise in P(λ) (relative ε ≈ 2·10⁻¹⁶) limits the root to about √ε relative
  precision.
- The square root in the RMSD amplifies that error when the RMSD is near 0.
- The frozen V2 fallback conditions (non-convergence, non-finite or out-of-bound λ) do not detect it.

**V2 stays frozen.** It is still run through its gate, which is expected to fail, and is reported as such.

**New arm V2b** = V2 plus one additional row fallback to `eigvalsh`:

  |P′(λ)| ≤ g · max(λ, 1 Å²)³, with g = 10⁻⁴

- **Basis:** P′(λ_max) = Π_{i>1}(λ_max − λ_i) vanishes at a multiple root.
- **Error for an unflagged row:** the root error is about ε·c·λ⁴ / |P′(λ)|. With |P′| > g·λ³, this is about
  c·ε·λ/g ≈ 2·10⁻¹² relative.
- **Choice of g:** fixed from that estimate and not tuned on any workload. It was checked only for the SYNTHETIC
  degenerate case, and for no fallback on generic SYNTHETIC chains.
- **Fallback frequency:** reported per workload, as for V2.

**Execution order:** V1 → V2 → V2b. Each arm has its own gate and its own timing rounds with a same-session V0.
The budget is unchanged at 1 CPU-h in total; arms not reached are reported as untested.
