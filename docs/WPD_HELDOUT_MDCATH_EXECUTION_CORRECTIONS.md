# WP-D Held-out (mdCATH): Execution Corrections EC1–EC5

**Status:** execution corrections under the 2026-10-06 user approval of Amendments 2–3. That approval authorized:
- research-harness changes;
- reserve downloads;
- one complete held-out evaluation.

**Scope.** These corrections change no method, weight, endpoint, margin, estimand, resampling scheme or selection
rule. All are recorded **before any held-out retrieval outcome was computed**. Earlier frozen records are preserved.
- The Amendment-3 runner is archived as `benchmarks/audit/wpd_local_retrieval/archive/heldout_runner_a3.py`, and its
  hash equals the A3 record.
- The final implementation hashes are bound in `results/wpd_heldout_scoring_authorization.json` before scoring.

## EC1. Cohort completion (implements original §4.3 replacement rule)

`mdcath_reserves.py` processes the frozen reserves **in order**, at most the 10 listed. For each one it:
1. downloads from the pinned revision;
2. verifies the recorded size and SHA-256 (mandatory; exit 2);
3. applies the unchanged structural eligibility criteria;
4. stops once 30 eligible domains exist.

Rules that are not changed:
- Eligible reserves replace failing primaries in order. The numbering criterion, the frame requirement and the
  selection are unchanged.
- If the queue is exhausted, the original insufficient-evidence rule applies.

Outputs: `results/wpd_mdcath_reserve_ledger.json` and `results/wpd_mdcath_cohort.json`, which records the accepted
domains and the excluded primaries with reasons.

## EC2. Sampled Horn-vs-Kabsch reference check (promised by original §5 via AdK §10 A2; was missing)

- Per window, the first 10 Spearman-sampled (q, d) pairs are recomputed with Kabsch SVD.
- This applies to B2 (Cα RMSD) and, where defined, to the T_pp target.
- A difference above **1e-6 Å** (the frozen tolerance) is a mandatory failure (exit 2).
- The maximum difference is recorded per domain.

## EC3. Executable scoring entry point with authorization validation

`heldout_runner.py score` / `score-heldout` (the official paths) checks, in order:

1. **Authorization record** (`results/wpd_heldout_scoring_authorization.json`). It must exist, parse as JSON and
   contain `authorized_by`, `scope`, `protocol_sha256`, `selection_sha256`, `cohort_sha256`,
   `implementation_sha256` and `cohort_domains`. Every hash must equal the current files:
   - protocol documents: the AdK protocol, the held-out protocol, Amendments 1–3 and this document;
   - the selection;
   - the accepted cohort;
   - the implementation modules: runner, readers, dependence and decision modules, `wpd_lib`, `legacy_lib`,
     `adk_lib`.
   Absent, malformed or mismatched records exit 3.
2. **Binding extension identity:** SHA-256 equal to the WP-A record. Failure exits 2.
3. **Production-source hashes:** the core/python crate sources and Cargo.lock equal the WP-A record. Failure exits 2.
4. **Input integrity:** every cohort file's size and SHA-256. Failure exits 2.
5. **Per-domain scoring:** mandatory checks inside are the reference check, the Horn–Kabsch check and target-oracle
   regret.
6. **Aggregation** with the frozen estimand.

If a mandatory check fails, `<out>.FAILED.json` is written, completed checkpoints are kept and the exit status is 2.

**Check classes are unchanged (Amendment 3 §1):**

| Class | Consequence |
|---|---|
| Mandatory | Exit 2 |
| Authorization | Exit 3 |
| Eligibility | Ledger entry |
| Descriptive | Recorded only |

## EC4. Per-domain checkpoints and safe resumption

- Each domain's result is written atomically to `build/heldout_checkpoints/<domain>.json` with a key: domain,
  domain index, input SHA-256, implementation hashes and protocol hashes.
- A checkpoint is reused **only** if its key equals the current key; otherwise the domain is recomputed.
- Domains are independent and seeds are deterministic (original §6.1), so execution order and parallelism do not
  change results. This is tested: resumed and uninterrupted execution give identical aggregates.

## EC5. Secondary reporting fields (original §7.6; "reported, not deciding")

These are added to each window record:
- the target spread (median T(q, d));
- the random overlap k/n_D.

The final results also contain:
- domain medians of overlap@10, NR@1, raw regret@10 and Spearman;
- target-spread terciles (windows pooled within a cell);
- coverage and undefined-window reasons;
- dependence-class counts.

None of these enters any decision.
