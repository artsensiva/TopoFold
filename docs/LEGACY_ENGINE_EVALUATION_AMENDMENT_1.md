# LE-1 Amendment 1: Analysis Estimand, H4 Wording, Attribution, Dependence and Verification Gates

**Amends:** `docs/LEGACY_ENGINE_EVALUATION_PROTOCOL.md`. That file is unchanged; SHA-256 `1ee58480…` in
`benchmarks/audit/legacy_engine/results/le_protocol_freeze.json`.

**Recorded:** 2026-10-05, while the LE-1 development run (second attempt) was computing per-window results, and
**before any LE-1 outcome existed or was inspected**.
- The first attempt stopped at 600 of 834 windows, without writing any result; its log contains only progress
  lines.
- The second attempt had produced no result file when this amendment was written.
- The record `results/le_amendment1.json` binds this file by hash.

**Unchanged:**
- The running code (`run_dev.py`, `legacy_lib.py`), methods, targets, splits, endpoints and ρ.
- The original analysis code `analyse_dev.py`. It is preserved unmodified and was never run on AdK data.

**Implementation:** the amended analysis is `analyse_dev_a1.py`; its tests are in `test_analyse_a1.py`
(SYNTHETIC).

---

## §1. Bootstrap estimand equals the reported point estimate (replaces protocol §7 "statistic is recomputed from block-sum matrices")

**Defect.** The point estimate is the **equal-window-weighted** mean: the mean over windows of each window's
mean over its valid queries. The original bootstrap recomputed a **pooled, query-weighted** mean from block
sums. The two differ whenever valid query counts differ across windows, and also where NR is undefined for
some (window, query) units.

**Amended statistic.** For method m in a cell, let x_{w,q} be NR@10 (NaN if undefined) and define

  θ_m = mean over windows w with ≥ 1 defined unit of [ mean over defined q of x_{w,q} ].

**Bootstrap.**
- Per window w and query block b, form S_{w,b} = Σ defined x_{w,q} and C_{w,b} = count of defined units.
- A resample draws query-block multiplicities c_b and window-block multiplicities, which give per-window
  multiplicities c_w.
- The resampled window mean is m*_w = Σ_b c_b S_{w,b} / Σ_b c_b C_{w,b}. It is undefined if the denominator
  is 0.
- The resampled statistic is θ* = Σ c_w m*_w / Σ c_w over windows with defined m*_w.
- With all multiplicities equal to 1, θ* = θ exactly; a regression test checks this.

**Paired contrasts.** For X − Y the statistic is computed on the per-unit difference x^X_{w,q} − x^Y_{w,q}.
It is undefined where either is undefined, and the same equal-window weighting is used. If the two methods'
undefined units coincide, this equals θ_X − θ_Y. Any difference between the two is reported.

**Counts reported.** Valid query counts per window (minimum, maximum, number of distinct values) and undefined
units per method and target. These are reported for each cell before any interval.

## §2. H4 wording (replaces the H4 row of protocol §7)

The original falsification text "Global legacy G_F ≥ local methods" did not state the direction in NR terms.
In NR, lower is better.

| Hypothesis | Contrast | Supported if | Falsified if |
|---|---|---|---|
| H4: local information explains historical advantages (global rankings retrieve worse local conformations than local rankings of the same family) | G_F − F and G_RMSD − B2 | Both intervals entirely above 0 (global worse) | Either interval entirely below 0 (global better than its local counterpart) |

Otherwise H4 is **unresolved**. F − B2 is reported as a separate descriptive contrast and is no longer part of
H4.

## §3. Attribution of F − A (clarifies protocol §7 H2/H3)

**F − A is an overall method comparison.** It combines four differences:
- flexible vs lock-step coupling;
- max vs RMS aggregation;
- omission of the last interior κ;
- per-pair √(Δκ² + Δτ²) vs channel-count-normalized averaging.

It is **not** attributed to residue misregistration alone. The pieces are attributed as follows:

| Contrast | What it isolates |
|---|---|
| **F − F_lock** | Flexible coupling alone: same pairs, same pointwise distance, same max aggregation. This is the only contrast interpreted as the effect of warping, i.e. of permitted misregistration |
| **F_lock − A** | The other three differences together. They are not separated further, because no additional variants are added |

H3 ("Fréchet harms same-sequence retrieval") is assessed with F − F_lock, not with F − A.

## §4. Temporal dependence and block length (replaces the governing-interval rule of protocol §7)

**Both analyses are reported.** Intervals from 100-frame and 250-frame query blocks are reported for every
contrast and method.

The τ_int diagnostic (WP-B `adk_lib.tau_int`, in frames) is computed on the per-query series of window-averaged
differences, in time order. It is reported for each contrast. The statement then follows this rule:

| τ_int | Statement based on |
|---|---|
| ≤ 100 frames | 100-frame interval |
| > 100 and ≤ 250 frames | 250-frame interval |
| > 250 frames, or not estimable | **No interval-based statement.** Reported as "dependence unresolved"; both intervals are shown, but neither governs |

**Limitations:**
- τ_int estimated from about 750–1,000 correlated frames is itself uncertain.
- Blocks longer than the correlation time can still miss slower drift, which WP-B observed in these runs.
- With 250-frame blocks there are only 4–5 query blocks per trajectory, so percentile intervals are coarse.
- The database frames are never resampled.

## §5. Verification gates before any retain/reject conclusion (supplements protocol §8)

An unqualified development outcome (retain, restrict or reject) for A or F is allowed only if **all** of the
following hold:
1. The binding τ/κ reference check passed on every frame. `run_dev.py` stops otherwise, so a missing check
   record counts as a failure.
2. The Rust/NumPy Fréchet difference is ≤ 1e-9 on every window. `run_dev.py` stops otherwise; the recorded
   maximum is re-checked.
3. No non-finite distance occurred for any method in {M1, M2, K, A, F, B1, B2}, in any primary cell.
4. No regret numerical failure occurred.
5. Undefined NR units are below 1 % of units for every method in the best set, in every primary cell.
6. The development classification is computed at **both** block lengths, and the two outcomes agree.
7. In every primary cell, τ_int of the per-query series of NR_X − NR_best (best = observed best) is ≤ 250
   frames. Otherwise the outcome is "dependence unresolved".

If any gate fails, the outcome is reported as **"insufficient evidence (verification gate: <reason>)"**. The
underlying classes are still shown, so that they are not hidden.

## §6. Citation evidence (clarifies protocol §11)

- **Crossref verification** establishes bibliographic metadata only: existence, title, year, DOI.
- It does **not** establish what a paper contains. Statements about the methods in prior work, and hence any
  novelty judgment, are marked as either:
  - **checked against the primary source in this audit** (only the mdCATH paper and code were read here); or
  - **from prior knowledge, not checked against the primary source in this audit**.
- Novelty conclusions that rest on unchecked content are labelled **provisional**.

## §7. Interrupted first attempt

- The first attempt (log preserved as `results/logs/le_run_attempt1.log`) ended at 600 of 834 windows. It had
  no result file and no exit-status line, and the process was not present when checked.
- **The cause was not established.** It coincided with the end of an assistant session, but no
  termination record exists.
- The second attempt runs the identical, hash-checked code and protocol from the start.
- Both attempts are documented in `results/le_run_attempts.json`.
