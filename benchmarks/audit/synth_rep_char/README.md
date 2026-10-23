# SYNTHETIC_REPCHAR-1: sensitivity of the replica-comparison rule (SYNTHETIC)

Characterization of `topofold_inputs.compare.replica_comparison` (the complete-separation rule used by DEV-REP-1) on
SYNTHETIC data with known truth. Specification: `docs/REPLICA_COMPARISON_CHARACTERIZATION_SPEC.md` and `run_spec.json`
(84 cells × 200 repetitions), frozen with the code before the run (`results/SYNTHETIC_REPCHAR-1_freeze.json`). No real
data or real feature arrays were read; DEV-REP-1 is not re-scored. 200 repetitions per cell give coarse Monte Carlo
precision (Wilson half-width up to about 0.07).

| File | Content |
|---|---|
| `synthetic_generators.py` | AR(1), wrapped-normal and Markov-mixture generators with exact stationary marginals (tests: `tests/test_synthetic_generators.py`) |
| `run_char.py` | `freeze` / `run`; hard 900 CPU-s cap and all-or-nothing cells |
| `CLARIFICATION_2026-10-23.md` | Dated clarification: meaning of 0/200 null cells, separate category meanings, limits of independent starts, and the checkpoint-status defect |
| `results/SYNTHETIC_REPCHAR-1_results.json`, `results/SYNTHETIC_REPCHAR-1_report.md` | Per-cell category frequencies (all three, with intervals), statistics, generated report |

```bash
.venv/bin/python benchmarks/audit/synth_rep_char/run_char.py run
```

## What the rule detects, misses and confuses (numbers from the results file)

- **Identical-distribution null (S0):**
  - no consistent or replica-dependent result in any of the 6 cells (0 of 200 each; upper 95 % bound 0.019);
  - the same holds with an outlier replica and no condition effect (S4a, Δ = 0);
  - random replica offsets with no condition effect (S4b, Δ = 0) gave at most 1 call in 200.
- **Detects** location shifts of 0.5σ (linear) or 0.5 s (circular) almost always when frames are independent
  (calls 200/200 and 199/200), scale changes of 1.5× (200/200), and population shifts 50/50 → 70/30 with independent
  states (183/200).
- **Misses, under temporal correlation:** with φ = 0.8 the 0.5σ shift is called in 61/200 repetitions; with φ = 0.95
  even 1σ is called in 64/200. Scale changes of 1.5× at φ = 0.95: 5/200.
- **Misses, slow population shifts:** 50/50 → 70/30 with dwell scale D = 20: 3/200.
- **Misses, rare excursions and heterogeneity:**
  - a rare excursion confined to one replica (S4c): 0/200 at every dwell scale;
  - a 1σ condition shift when replicas carry random offsets of SD 0.5σ: 20/200 at φ = 0.
- **Confuses, condition-confounded sampling (S5):** when each condition's replicas share a realized offset, the rule
  calls a "difference" in 81/200 (offset SD 0.25σ) and 138/200 (0.5σ) repetitions at φ = 0. The calls track the
  realized offset difference, not a condition effect. The rule cannot distinguish the two; only independent starts
  per condition can.
- **Category shifts with replica count (S6, Δ = 1σ, φ = 0.8):**
  - consistent falls from 186/200 (3 replicas) to 102/200 (8 replicas) as replica-dependent rises;
  - total calls stay at 183–200 of 200;
  - more replicas make complete separation harder to reach, because more pairs push the extremes apart.
- **An outlier replica plus a true 1σ shift (S4a)** is classified replica-dependent (200/200 at φ = 0), never
  consistent.

The effective sample sizes quoted for φ (≈ 387, 43, 10) are Gaussian sample-mean approximations, not guarantees for
histogram distances.

## Recommendation

**Retain the rule as a conservative descriptive screen:** it made no calls under the identical-distribution null in
this grid. Report its negatives only as "no window met the complete-separation criterion", alongside this
characterization.

It is not suitable as evidence of absence:
- it has low sensitivity under realistic temporal correlation;
- it is blind to slow population shifts and single-replica rare events.

It is not suitable as evidence of a condition effect without independent starting conditions, because it confounds
condition-shared sampling.

An alternative rule (for example a replica-label permutation or a replica-bootstrap effect size) could be
investigated, but it remains a separate decision. It could improve sensitivity under correlation; no analysis can
remove the confounding in S5.
