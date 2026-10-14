# Proposal for Review: Bounded Library Integration of V2b Prepared-Database Search (not implemented)

> **Historical — not approved (status note added 2026-10-14).** This proposal was written after BQ-3 and is kept
> as a record. Following CMP-1 (`docs/CMP1_CLOSEOUT_2026-10-13.md`: NO_DECISIVE_ADVANTAGE; V2b integration
> unsupported), neither option below is approved or pursued. V2b remains an audit artifact. Any reconsideration
> needs a new, reviewed decision.

**Status:** proposal only, 2026-10-13. Nothing is implemented, and no production code has changed.

**Basis:** BQ-3 (`benchmarks/audit/bq3/REPORT_SYNTHETIC_P2.md`), on SYNTHETIC_tiled, performance-only P2 rows.
- **Prepared-database operation qualified:** every one of 15 measured query latencies was ≤ 10 s, peak RSS was within
  8 GiB, and verification passed.
- **First-query operation did not qualify:** the measured preparation + first-query spans exceeded 10 s.

**What this does not mean:** these results describe a NumPy harness on one machine. They establish no biological
capability, no novelty over existing RMSD software, and no production readiness. Integration is an engineering
convenience to be judged on its own merits.

## 1. Scope

- **Proposed function:** one exhaustive, single-query window Cα RMSD top-k search with an explicit prepared database:
  `prepare(database, windows) -> PreparedDB` and `query(PreparedDB, query_frame) -> top-k per window`.
- **Semantics:** unchanged from the harness. Minimal proper-rotation RMSD (Horn key matrix; QCP Newton from
  E0 = (G_A + G_B)/2 with the V2b fallbacks), float64, and the `topk_lex` tie policy: ascending (distance, database
  position).
- **Out of scope:**
  - pruning or approximate search;
  - float32;
  - multithreading;
  - changes to the existing Fréchet index or the bimodality scan;
  - any claim about pocket detection, allostery or biology.

## 2. Options for review (choose one; no default)

- **(a) Python-level module** in `python/topofold/`, ported from the harness. This is the smallest change and keeps
  NumPy performance characteristics.
- **(b) Rust implementation** in `topofold-index` with bindings. This is a larger change; it could address
  first-query latency, but that is unmeasured. It needs its own frozen benchmark against the harness before any
  performance claim.

## 3. Required gates before merge

1. **Parity:** the library output equals the harness `b2opt_lib` output (identical top-k hashes) on the BQ-1/BQ-2
   development workloads and a P2-dimension SYNTHETIC check, under the same verification rule as BQ-3 (§3).
2. **References:** max |Δd| ≤ 10⁻⁶ Å against `ref_svd` and MDAnalysis `qcprot` on samples. Degenerate-geometry,
   duplicate, permutation, rigid-transform and near-tie tests are ported from `test_b2opt.py`.
3. **Memory:** prepared arrays are n_D·n_win·(3W + 1)·8 B, about 2 GB at P2. The API must expose this cost. Peak RSS
   is measured, not assumed.
4. **Documentation:** first-query latency (preparation + query) is documented as **not** meeting the 10 s target at
   P2 on the tested machine, and preparation is a database-side cost.
5. **Test and build discipline:**
   - `cargo test --workspace` and `pytest tests/` pass;
   - no change to existing APIs;
   - the KNOWN_ISSUES MSRV item (m12) is considered.

## 4. Prerequisites outside integration

These remain separate priorities and are not implied by this integration:
- strong equal-information comparisons against established tools (e.g. MDAnalysis, MDTraj) on identical tasks and
  hardware;
- biological validation per application.
