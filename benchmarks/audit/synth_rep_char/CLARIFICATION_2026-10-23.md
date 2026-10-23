# SYNTHETIC_REPCHAR-1 Clarification (2026-10-23)

**Scope:** interpretation and one documented defect. The frozen specification, runner, results and report are
unchanged. Nothing was rerun.

## Interpretation

- **Zero calls in a null cell** means 0 of 200 repetitions observed, with an upper 95 % Wilson bound of about 1.9 %
  (exactly 3.8415 / (200 + 3.8415) = 0.0188). It describes those generators and settings only. It does **not**
  establish general false-positive control for other distributions, correlation structures, replica counts or
  multiple-window use.
- **Consistent and replica-dependent outcomes mean different things.**
  - *Consistent*: complete separation of the two conditions with all replicas included.
  - *Replica-dependent*: separation appears only after one specific replica is left out, so the description hinges
    on that replica (for example an outlier).
  - The two are reported separately and should not be pooled into one "detection" without their breakdown.
- **Independent starting conditions** per condition would remove the shared-sampling confound that S5 models. That
  strengthens a design but does not by itself guarantee causal attribution: other condition-correlated differences in
  setup, protocol or sampling can remain.

## Defect: intermediate files reported "complete"

- **What:** `run_char.py` initializes `status = "complete"` and rewrites the results file after every finished cell.
  Every intermediate checkpoint therefore carried `"status": "complete"` before all cells had finished. Only a
  stop by the budget rule would have set `"partial"`. A process killed by the hard CPU limit would have left a file
  wrongly marked complete.
- **Effect on the result:** none. The run finished all 84 cells within the cap (289 CPU-s). The final file lists 84
  cells, and its "complete" status is correct.
- **Handling:** the frozen runner is not edited. Later runners mark checkpoints `partial` until every planned cell has
  finished (first applied in SYNTHETIC_PERM-1).
