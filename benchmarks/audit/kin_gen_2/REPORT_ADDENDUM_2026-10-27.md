# KIN-GEN-2 correction addendum (2026-10-27)

The outcome in `REPORT.md` (**insufficient evidence**) and every frozen artifact are unchanged. This addendum qualifies how the study may be read; it does not upgrade or downgrade the outcome. Prompted by an external review of the code and records, each point below was checked against the repository.

## 1. Correction of the error-count wording

`REPORT.md` says the differences between the representations "correspond to one or two entries". That is wrong as a description of the errors. From `results/eval_results.json` (rows with prediction different from label):

| Representation | Weighted BA | Errors |
|---|---|---|
| T (κ/τ) | 0.964 | 5 false positives (in predicted out), all FGFR1 |
| B1 (window distances) | 0.979 | 3 false positives, all FGFR1 |
| B2 (local RMSD) | 0.958 | 1 false negative (out predicted in), MET |
| B3 (all 36 distances) | 0.993 | 1 false positive, FGFR1 |

The error types differ (T, B1, B3 err on in entries of one kinase; B2 misses one of eight out entries), and the similar weighted scores partly reflect the kinase × class cell weighting (a single FGFR1 or MET entry carries the weight of its cell). Eight kinases were processed, but PLK1 gave no usable entry, only seven contributed predictions, and only JAK2 and MET contributed the eight out entries. These are descriptive counts, not a comparison.

## 2. Methodological qualifications

| Issue (verified in code or protocol) | Consequence | Needed in a future study |
|---|---|---|
| `Model` scales with all development kinases, and `loko` reuses that model, so each held-out kinase influenced the scaling used to select k and the primary baseline. This is inside development cross-validation; no evaluation entry influenced scaling, k, boundaries or baseline choice. | Development leave-one-kinase-out figures (T 0.975 etc.) are optimistic and must not be read as performance. | Fit preprocessing inside each training fold. |
| Development selection of k and the primary baseline uses pooled balanced accuracy over entries (`balanced_accuracy`), whereas the primary evaluation metric is the kinase × class weighted one. The protocol says "leave-one-kinase-out balanced accuracy" without defining the pooling. | Heavily sampled kinases (AURKA, PRKACA) influence selection more than the evaluation metric would. | State the selection metric in the protocol and make it match the evaluation weighting. |
| The αC Glu anchor comes from contacts observed across the candidate kinase's own entries. | Eligibility depends partly on the conformation under study; CDK2 (development) and PLK1 (evaluation) gave no entry. | Establish the αC correspondence from sequence or structure without requiring the conformation being classified. |
| Labels are D1/D2 clusters fitted on development data. | A surrogate target; boundary behaviour and biological meaning are not independently established. | A frozen, independently validated reference rule. |
| "Comparable" was defined as a point-estimate difference of at most 0.05 in weighted BA. | A small observed difference is not equivalence or non-inferiority. | A paired comparison with enough precision and a justified margin. |
| Evaluation downloaded 96.5 MB; the protocol stopping rule capped it at 90 MB (overall limit 100 MB). | The 100 MB limit held; the 90 MB stopping rule was not enforced, because it was checked before each kinase. | Check the remaining allowance before each download. |

No direct tuning on evaluation scores and no use of evaluation structures as training examples was found.

## 3. Reading the next step

Collecting more DFG-out structures would address only the count gate and leave the points above in place. Selecting DFG-out-rich kinases defines a different, enriched population: it would be a new, separately frozen study, not an amendment that upgrades this outcome. This corrects my earlier statement that more out entries were "the only thing missing".

## 4. Test coverage

The clustering test needs SciPy and was skipped where SciPy is absent. CI now installs SciPy in each Python matrix job and runs `tests/test_kin_gen_2.py` with unexpected skips treated as failures (only the test needing git-ignored public data may skip; the matrix checkout has full history for the freeze test).
