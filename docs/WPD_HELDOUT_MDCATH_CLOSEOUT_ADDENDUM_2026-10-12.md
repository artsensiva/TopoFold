# WP-D Held-out (mdCATH): Closeout Addendum, 2026-10-12

**Nature.** These are interpretation corrections only.
- No rerun, no tuning.
- No change to protocols, the scoring implementation, raw results, checkpoints or declared outcomes.
- Values are read by `make_closeout_addendum.py` from the frozen results (SHA-256 `8dbdcf8c813c4a49…`, matching `results/wpd_heldout_run_record.json`).
- Machine-readable record: `benchmarks/audit/wpd_local_retrieval/results/wpd_heldout_closeout_addendum.json`.

The corrected statements were made in the assistant's chat summary of the run (2026-10-06). The generated report `REPORT_HELDOUT.md` contains the classes and intervals and is unchanged.

## C1. B2 at W = 8 on T_φψ: inconclusive, not "beyond the margin"

The summary said that at W = 8 method A beats B2 "with relative-loss intervals that clear the margin". That is wrong for B2. Its relative-loss intervals cross ρ = 0.10, so its declared class there is **inconclusive**:

| Cell | L(B2) | 95 % interval | Contains ρ | Declared class |
|---|---|---|---|---|
| 320/W8/phipsi | 0.114 | [0.075, 0.180] | True | inconclusive |
| 348/W8/phipsi | 0.129 | [0.079, 0.180] | True | inconclusive |

What is supported is narrower:
- A was the observed best in these cells and is within the margin.
- Whether B2 loses more than the declared acceptable fraction of A's gain is **not resolved**.

## C2. W = 5 on T_φψ does not establish equivalence of A and B2

The summary called A and B2 at W = 5 "statistically indistinguishable". The classes say only that both relative-loss intervals lie below ρ, i.e. both are within the declared acceptable-loss margin (Amendment 1 §4).
- That is neither an equivalence nor a non-inferiority claim.
- The intervals are marginal and unadjusted.
- The best method switched between A and B2 across resamples:

| Cell | Observed best | Best counts over resamples | L(A) [95 %] | L(B2) [95 %] |
|---|---|---|---|---|
| 320/W5/phipsi | A | {'B2': 650, 'B1': 5, 'B3': 0, 'M1': 0, 'M2': 0, 'A': 1345} | 0.000 [0.000, 0.022] | 0.002 [0.000, 0.025] |
| 348/W5/phipsi | B2 | {'B2': 1023, 'B1': 0, 'B3': 0, 'M1': 0, 'M2': 0, 'A': 977} | 0.001 [0.000, 0.030] | 0.000 [0.000, 0.020] |

## C3. C1 failing the competitive criterion is not absence of gain over random

The summary said that "every method except C1 is clearly better than random". That is wrong.
- C1's declared outcome *fails the criterion in all cells* refers to its relative loss against the best method (loss class *beyond the margin* everywhere).
- Its gain over random retrieval, G = 1 − median NR (coverage-matched), has a 95 % interval **above 0 in every cell**.
- On T_φψ its gain class is inconclusive or uninformative relative to ρ, not "no gain":

| Cell | G(C1) | 95 % interval | Lower bound > 0 | Gain class vs ρ | Loss class |
|---|---|---|---|---|---|
| 320/W5/phipsi | 0.097 | [0.084, 0.112] | True | inconclusive | beyond_margin |
| 320/W8/phipsi | 0.087 | [0.069, 0.096] | True | uninformative | beyond_margin |
| 348/W5/phipsi | 0.117 | [0.081, 0.140] | True | inconclusive | beyond_margin |
| 348/W8/phipsi | 0.088 | [0.076, 0.114] | True | inconclusive | beyond_margin |
| 320/W5/pp | 0.175 | [0.153, 0.228] | True | informative | beyond_margin |
| 320/W8/pp | 0.144 | [0.134, 0.165] | True | informative | beyond_margin |
| 348/W5/pp | 0.204 | [0.161, 0.233] | True | informative | beyond_margin |
| 348/W8/pp | 0.159 | [0.135, 0.176] | True | informative | beyond_margin |

"Uninformative" in Amendment 1 means the upper bound of G is at or below ρ. That is a gain smaller than the declared threshold, not a gain of zero.

## C4. Standing interpretation

**Candidates for their respective geometric tasks:**
- **A** (residue-aligned κ–τ): Cα-only retrieval judged by an angular φ/ψ target. Declared T_φψ outcome: `meets_criterion_in_all_cells`.
- **B2** (window Cα RMSD): Cα-only retrieval judged by a Cartesian peptide-plane target. Declared T_pp outcome: `meets_criterion_in_all_cells`.

**Scope.** These are within-dataset generalization results:
- mdCATH, CHARMM22*, 320/348 K;
- 30 superfamilies;
- cross-replica retrieval from a shared starting configuration.

**Not evidence of:**
- the legacy κ–τ Fréchet index or its writhe cascade (neither was part of this evaluation);
- any detector (bimodality scan, pocket, allostery);
- biological function;
- clinical or cancer-treatment efficacy.

