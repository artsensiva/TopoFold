# WP-D Held-out Protocol (mdCATH): Amendment 1

**Amends:** `docs/WPD_HELDOUT_MDCATH_PROTOCOL.md`. That file is unchanged; its SHA-256 `2bc16aa7…` is recorded in
`benchmarks/audit/wpd_local_retrieval/results/wpd_heldout_protocol_freeze.json`.

**Recorded:** after the metadata-only domain draw and before any mdCATH trajectory download or scoring. Bound by
hash in `results/wpd_heldout_amendment1.json`.

**Unchanged:**
- Data source, eligibility and exclusions.
- The draw, including the exact primary and reserve order in `results/wpd_mdcath_selection.json`.
- Trajectory handling, metrics and per-window aggregation.
- Seeds.
- ρ = 0.10 (not re-tuned; no new AdK illustration was computed).
- The AdK evidence.

**Implementation:** `benchmarks/audit/wpd_local_retrieval/heldout_decision.py`.

**Tests:** `test_heldout_decision.py`, on SYNTHETIC fixtures.

Sections below replace the cited parts of the original protocol. Anything not cited stays in force.

---

## §1. Consistent estimands (replaces original §7.1 "Relative loss" and the matching part of §7.2)

For cell c (temperature × W), let NR_m(d) be the per-domain value (original §6.3), and let med_m be the median of
NR_m(d) over the cell's domains.

- **Best method.** b = argmin over {M1, M2, B1, B2, B3} of med_m. Ties within 1e-12 are broken by the order
  B2, B1, B3, M1, M2. This is unchanged.
- **Gain over random.** G_m = 1 − med_m. This is unchanged.
- **Relative loss (new):**

  L_m = (med_m − med_b) / (1 − med_b)

  This is the difference of the **same medians** that select b, divided by the best method's gain.

**Why this replaces the median of paired differences.**
- The original L used median_d[NR_m(d) − NR_b(d)] while b was selected by med. The two are different
  functionals.
- A method whose median ties with the best could therefore show a non-zero loss, and a non-best method could
  show a negative loss. The selection and the margin would then disagree.
- With the amended L, the best method has L = 0 and every other method has L ≥ 0, by construction. L is also
  in the same units as the quantity used to rank methods.

**Cost.** Pairing no longer enters the point estimate. It still enters the uncertainty:
- every resample draws one multiset of domain indices, used for every method (paired resampling);
- **b is re-selected within each resample**, so the interval includes uncertainty about which method is best.

## §2. Invalid denominators (new; applies to the observed sample and to every resample)

The denominator is the best method's gain G_b = 1 − med_b.

- **Invalid denominator:** G_b ≤ ε_G with ε_G = 1e-6. This covers zero, negative (best method no better than
  random) and numerically negligible gain.
  - Per-domain NR values are averages of float64 terms, so their round-off is far below 1e-6.
  - 1e-6 is far below ρ, so this threshold only flags gains that are effectively zero.
  - Gains between ε_G and ρ give valid but large, unstable L values. The intervals show that instability; no
    other cap is applied.
- **Observed sample invalid:** L is **not assessable** for every method in that cell. It is reported as such and
  is never set to 0 or ∞. Gain classes are still computed, because G has no denominator.
- **Resamples:** invalid resamples are **counted, not discarded**. The 95 % interval uses order statistics
  (`numpy.quantile(..., method='inverted_cdf')`), with every invalid resample placed at **−∞ for the lower bound
  and +∞ for the upper bound**.
  - If ≤ 2.5 % of resamples are invalid, both bounds are finite but no narrower than they would be without
    those resamples.
  - If > 2.5 % are invalid, both bounds are infinite, and the cell's loss class becomes **inconclusive**.
  - The invalid fraction is reported for every cell.
- **When the result is inconclusive** for a method's loss in a cell: the observed denominator is invalid
  ("not assessable", which counts as inconclusive in §4), or the conservative interval contains ρ, including
  infinite bounds.

## §3. C1 coverage-matched population (replaces original §6.4 and the C1 parts of §7.1)

**Common population.**
- For window w of domain d, let Q₁(w) be the valid queries with a C1 result: R ≥ 0.1, at least k finite database
  distances, and defined NR.
- For **every** method m, the coverage-matched window value is the mean NR@10 of m over Q₁(w).
- A window enters the matched set only if Q₁(w) is non-empty and every method has a defined value there.
- The domain's matched value NR^cm_m(d) is the median over its matched windows. It is defined only if the domain
  has ≥ 10 matched windows.
- D₁ is the set of domains with defined matched values. The same D₁ is used for every method; the
  implementation rejects inputs where methods differ.

**C1 statistics, on D₁ (with resample multiplicity):**
- G_C1 = 1 − median NR^cm_C1.
- L_C1 = (median NR^cm_C1 − median NR^cm_b) / (1 − median NR^cm_b).
- Numerator and denominator both use the matched population. b is the best method selected on the **primary**
  values (§1) in the same sample or resample.
- The §2 rules apply to the matched denominator.

**Other rules:**
- C1 is never eligible as the best method.
- L_C1 may be negative. This means C1 has a lower matched median than b on its own population. It is reported
  as computed.
- C1 is **not assessable** in a cell if |D₁| < 24. It is reported beside its coverage: D₁ size and the query
  share with a C1 result.

## §4. Labels and interpretation (replaces wording in original §7.3 and §7.5)

### Per cell

| Old label | New label | Rule (unchanged thresholds, ρ = 0.10) |
|---|---|---|
| E "practically equivalent" | **within the declared acceptable-loss margin** | Upper bound of L ≤ ρ |
| I "inferior" | **beyond the declared acceptable-loss margin** | Lower bound of L > ρ |
| U | **inconclusive** | Otherwise; or not assessable (§2, §3) |
| Informative / Uninformative / Inconclusive | unchanged | G interval vs ρ |

### Overall (precedence unchanged)

| Order | Outcome | Rule |
|---|---|---|
| 0 | **Insufficient evidence** | Unchanged |
| 1 | **Meets the declared criterion in all cells** (was "Selected") | Within margin and informative in all four cells |
| 2 | **Fails the declared criterion in all cells** (was "Not selected") | Beyond margin in all four cells, or uninformative in all four |
| 3 | **Cell-dependent** | At least one cell within margin and informative, and at least one beyond margin |
| 4 | **Inconclusive** | Otherwise |

### Interpretation limits (binding on the report)

- "Within the margin" means only that the interval for the loss of gain relative to the best method lies below
  the declared ρ. It is **not** a claim of equivalence, non-inferiority or superiority, and ρ is a stated value
  judgment.
- The intervals are **marginal, unadjusted** 95 % percentile intervals: 6 methods × 4 cells × 2 statistics. They
  have **no family-wise coverage**. Across many classifications, some will be misclassified by chance.
- The overall outcome is a pre-declared summary of these classifications, not a test. "Meets the declared
  criterion" makes a method a **candidate** for engineering consideration. Choosing among candidates uses
  criteria outside this protocol, stated at the time.
- The results describe the drawn superfamily sample, the conditions in original §1 and §5, and cross-replica
  retrieval from a shared starting configuration.

## §5. Frame-count arithmetic (corrects original §4.3)

- The original justification, "numFrames ≥ 200 … so k = 10 is ≤ 5 % of any database after burn-in", is wrong.
- After excluding ⌊0.10 n⌋ frames, a replica with n = 200 keeps 180 frames. n − ⌊0.10 n⌋ is non-decreasing
  in n, so 180 is the minimum for every n ≥ 200.
- So k = 10 is at most 10/180 ≈ 5.6 % of a database. The query set likewise has at least 180 frames.
- The **rule itself (numFrames ≥ 200) is unchanged**; only the stated consequence is corrected.

## §6. Tests (SYNTHETIC; `test_heldout_decision.py`)

| Requirement | Test |
|---|---|
| Tied medians with unequal paired differences | `test_tied_medians_unequal_paired_differences_give_zero_loss` (old estimand would give 0.5) |
| Estimand is the difference of medians, ≥ 0 | `test_loss_is_difference_of_medians_and_nonnegative` |
| Identical methods | `test_identical_methods_zero_loss_in_every_resample` |
| Zero, negative and negligible best gain (observed) | `test_observed_nonpositive_or_negligible_best_gain_not_assessable` (gain 0, −0.2, ε_G/2) |
| Invalid bootstrap denominators | `test_conservative_interval_places_invalid_at_both_tails`, `test_resamples_with_invalid_best_gain_are_counted_not_dropped` |
| C1 coverage matching | `test_c1_uses_common_population_and_is_never_best`, `test_c1_insufficient_coverage_not_assessable_and_mismatched_population_rejected` |
| Frame arithmetic | `test_post_burn_frames_minimum_and_k_fraction` |
| Outcome precedence and labels | `test_overall_outcome_precedence` |

**Sensitivity check (scratch copies only):**
- Replacing the estimand with the median of paired differences fails the first two tests.
- Discarding invalid resamples fails both invalid-denominator tests.
