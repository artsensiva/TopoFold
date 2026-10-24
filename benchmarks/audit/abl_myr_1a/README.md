# ABL-MYR-1a: post-hoc amendment of ABL-MYR-1 (initialization audit; verdict V4, no comparison run)

**This is not the original frozen test.** ABL-MYR-1 (`../abl_myr_1/`: specification, runner, code freeze, blocked G1 record,
member audit) is preserved unchanged; its hashes are re-verified inside this amendment's freeze record. Everything here is
labelled amended. Engineering audit of starting files only: **no trajectory was opened, no descriptor computed, no condition
difference scored.** Nothing here is evidence about Ayaz WT/M472I, HDX, drug response or biological function.

## Why an amendment

The published model archives contain no `Set2/system.gro`, which the original specification required (original run: G0
passed, G1 failed; `../abl_myr_1/results/ABL-MYR-1_g1_member_audit.json`). This amendment does not substitute
`system_em.gro` for `system.gro` in any original rule. It first defines, freezes and runs a bounded audit of the files that
do exist, and fixes in advance what each can and cannot establish.

## What was frozen, and when

| Item | Record |
|---|---|
| Audit specification `audit_spec.json`, audit code `run_init_audit.py`, tests `tests/test_abl_myr_1a_audit.py` | `results/ABL-MYR-1a_audit_freeze.json` (2026-10-24T21:18:00Z) |
| State at freeze | no member of the model archives extracted; only zip member names and the dataset README had been seen |
| Fixed in advance | file capabilities and limits, measures, categories, thresholds, stop conditions S1–S4, ordered verdict rules, and the consequence of each verdict |
| Original artifacts | 7 files hashed in the freeze record (and unchanged at run time) |

## Result (generated: `results/ABL-MYR-1a_audit_results.json`, `results/ABL-MYR-1a_audit_report.md`)

**Verdict (by the pre-declared ordered rules): `V4_set1_comparability_failed`. No comparison of any kind was run.**

| Question | Observed (see the generated report for the numbers) |
|---|---|
| Set1 `system.gro`, Myr vs noMyr (Q1) | **identical** protein coordinates (all protein atoms, max difference 0.0) |
| Set2 `system_em.gro` vs Set1 `system_em.gro`, per condition (Q4) | **identical** for Myr and for noMyr: the Set2 post-minimization protein coordinates are copies of the Set1 ones |
| Set2 start, derived (`set2_start_direct_evidence`) | **true** (Q4 identical in both conditions and Q1 identical) |
| `production_md.mdp`, Myr vs noMyr (Q5), both sets | **not identical**: the only differing lines are `compressed-x-grps` and `tc_grps`, which name the `MYRP` group in the Myr system |
| `topol.top` after removing Myr lines (Q6), both sets | **not identical**: the water count differs (`SOL` 29949 vs 29938) |
| Set1 vs Set2 parameters (Q5) | differ only in `tcoupl` (V-rescale vs Nosé–Hoover), as the paper states |
| `md.tpr` run inputs (Q7) | all 40 readable with MDAnalysis; one protein-parameter digest per set and condition, equal between conditions; MDAnalysis also returned positions (recorded, not used by the frozen rules) |
| Initialization stated in the mdp (Q8) | `gen_vel = no`, `continuation = yes` in all four files: velocities were carried over from equilibration; seeds and replicate-specific velocities are not stated and not verified |
| Q2/Q3: `system_em.gro` Myr vs noMyr | category `differs_large` in both sets with identical numbers, as expected from Q4. The frozen measure has **no superposition**, so box placement or periodic wrapping alone can produce this; the audit did not test that. It does not enter the verdict when Q4 is identical |

**Why the verdict is V4:** the rule requires all three components of the original G4 for Set1 (coordinates, mdp, topology).
Coordinates pass; the mdp and topology components fail. These two failures are not defects of the Set1 data: a Myr system
must name its Myr group in the index-group lines, and must contain a different number of waters from a system without Myr.
The original G4 text ("production_md.mdp identical", "every non-Myr line of topol.top identical, including water counts")
cannot be satisfied by any Myr vs noMyr pair, in either set. That is a design flaw of the original specification, found only
by reading the start files. The frozen rule was applied as written and was **not relaxed after seeing the data**.

## What is and is not established

- **Established (direct, audit):**
  - the Set1 conditions start from identical protein coordinates;
  - the Set2 post-minimization structures are copies of Set1's;
  - the protein force-field parameters of the production inputs match between conditions;
  - the parameter differences are confined to the Myr group names and the water count.
- **Not established:** that these differences are immaterial for the dynamics; the pre-minimization Set2 start (unpublished,
  an inference from the copy evidence that a different one is unlikely is untested); replicate-specific velocities or seeds; anything about trajectories.

## Owner decision

Any continuation is a **second dated post-hoc amendment** (ABL-MYR-1b), designed after seeing these start files, and before
any trajectory is opened. None was started. Options:

1. **Amend the parameter rule for a Myr vs noMyr contrast:** accept mdp differences confined to index-group names that
  mention `MYRP`, and topology differences confined to the Myr molecule and the water count, while keeping the coordinate
  and protein-parameter checks. With the Set2 start now established directly, a two-set comparison could then be specified
  and frozen before scoring. The amended specification would have to disclose that its rule was written after the start
  files were read.
2. **Stop ABL-MYR-1.**
3. **Additionally audit replicate starts** from the readable `md.tpr` coordinates (not frozen or done here).

## Command and resources

```bash
PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1a/run_init_audit.py run \
    --adapter-python benchmarks/audit/wpd_local_retrieval/build/venv/bin/python
PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1a/run_init_audit.py report
```

Resources: 65.7 CPU-s, 64 s wall, 164 MB extracted (54 members), nothing downloaded; budget 0 B download, 1 GB extracted,
600 CPU-s. The extracted audit files are in the git-ignored `data_dryad/extracted_audit/`.
