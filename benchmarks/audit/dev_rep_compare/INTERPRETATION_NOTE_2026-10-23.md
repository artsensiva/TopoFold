# DEV-REP-1 Interpretation Note and Closeout (2026-10-23)

**Scope:** wording and reproducibility only. DEV-REP-1 was not rerun; its specification, freeze, results, report and
post-hoc summary are unchanged. The original outcome stands as recorded.

## Corrected interpretation

- **Recorded outcome:** all 194 channel–window evaluations were classified *inconclusive* under the frozen rule
  (κ 49, Cα(i)–Cα(i+2) 49, τ 48, Cα(i)–Cα(i+3) 48). The four channels share windows: 49 three-residue windows serve κ
  and its distance baseline, and 48 four-residue windows serve τ and its distance baseline. The 194 evaluations are
  therefore not 194 independent findings.
- **Correction:** the closeout summary of 2026-10-23 described this as "no difference beyond replica variation". The
  correct statement is: **no window met the complete-separation criterion**. That criterion requires every one of the
  25 between-temperature replica-pair distances to exceed every one of the 20 within-temperature distances.
- **What it does not mean:** it is not evidence that 320 K and 348 K give the same local geometry, nor an estimate
  of effect size. How sensitive the criterion is has not been characterised; see the proposed characterization
  (`docs/REPLICA_COMPARISON_CHARACTERIZATION_SPEC.md`).
  *Later the same day:* the criterion was characterised on SYNTHETIC data (SYNTHETIC_REPCHAR-1,
  `benchmarks/audit/synth_rep_char/`). That study describes the rule's behaviour; it does not re-score DEV-REP-1.
- **Unchanged limits:** the five replicas share one starting configuration; pairwise distances are dependent; the
  comparison is of per-window marginal distributions; 320 K vs 348 K is a simulation condition with no experimental
  endpoint. Development evidence only.

## Reproducing DEV-ENS-1 after the shared-code change

DEV-REP-1 appended `replica_comparison` to `python/topofold_inputs/compare.py` and one test to
`tests/test_ensemble_compare.py`. The change is append-only (+54 and +23 lines, no line removed), but those two files
are in the DEV-ENS-1 freeze, so `run_compare.py run` on the current tree stops with a hash mismatch. This is intended:
DEV-ENS-1 must be reproduced from its own commit. All eight files in `DEV-ENS-1_freeze.json` match their recorded
SHA-256 at commit `644b4cd` (checked 2026-10-23).

```bash
git worktree add ../TopoFold-dev-ens-1 644b4cd
cd ../TopoFold-dev-ens-1
# git-ignored inputs and the adapter environment are not in the worktree; link the existing ones
ln -s "$OLDPWD/benchmarks/audit/adk_feasibility/data" benchmarks/audit/adk_feasibility/data
mkdir -p benchmarks/audit/wpd_local_retrieval/build
ln -s "$OLDPWD/benchmarks/audit/wpd_local_retrieval/build/venv" benchmarks/audit/wpd_local_retrieval/build/venv
"$OLDPWD/.venv/bin/python" benchmarks/audit/dev_ens_compare/run_compare.py run
```

- The runner reads the committed `results/DEV-ENS-1_freeze.json` in the worktree and checks the spec and code
  against it. It overwrites the worktree's copies of the results, not the main tree's.
- The `topofold` extension comes from the main `.venv`. Its native build may differ from the one used originally;
  the features provenance records the wrapper and native-extension hashes separately for comparison.
- Remove the worktree afterwards with `git worktree remove ../TopoFold-dev-ens-1`.
