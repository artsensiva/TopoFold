# ABL-MYR-1b: post-hoc comparability check of the production run inputs (outcome: not_comparable, K1; no protocol written)

**Post-hoc engineering check of run inputs and parameter files. No trajectory was opened, no descriptor computed, no
condition difference in dynamics scored.** ABL-MYR-1 (stopped at G1) and ABL-MYR-1a (verdict V4) are preserved unchanged
and remain valid: 12 of their files are hash-recorded in the freeze records here and verified at every run. Nothing here is
evidence about Ayaz WT/M472I, HDX, drug response or biological function, and no causal Myr claim is made.

## Sequence

| Step | Record | Result |
|---|---|---|
| Check frozen: `check_spec.json`, `run_comp_check.py`, `tests/test_abl_myr_1b_check.py` | `results/ABL-MYR-1b_check_freeze.json` (2026-10-24T21:44:11Z) | Questions, rules, constants, criteria K1–K5, stop conditions and budget fixed before any md.tpr coordinate, `nvt.mdp`, `npt.mdp` or `index.ndx` content was read |
| Attempt 1 (frozen code) | `results/ABL-MYR-1b_check_results_attempt1.json`, `…_report_attempt1.md`, `…_production_start_ca_attempt1.npz` | **Stopped at S3**: `duplicate index group 'MYRP'`. No md.tpr was loaded |
| Amendment a1 (dated): `check_amendment_a1.json`, `run_comp_check_a1.py`, `tests/test_abl_myr_1b_check_a1.py` | `results/ABL-MYR-1b_check_a1_freeze.json` (2026-10-24T21:46:49Z) | Implementation defect: GROMACS index files may repeat group names (here the `make_ndx` default residue-name groups). a1 accepts repeated names, treats a tc group as ambiguous if its occurrences differ, and otherwise imports the frozen check byte-identical |
| a1 run | `results/ABL-MYR-1b_check_a1_results.json`, `results/ABL-MYR-1b_check_a1_report.md` (generated) | **`not_comparable`: K1 failed** (0 of 40 run inputs validated); K2, K3, K4 and K5 passed |

## Why K1 failed (observed cause)

- MDAnalysis 2.10.0, the only tpr reader available, returns finite positions, velocities and the full topology for all 40
  md.tpr files (94,607 / 94,530 atoms and 64,527 / 64,462 bonds, Myr / noMyr), but **no simulation box** (`dimensions`
  is None in 40 of 40).
- The tpr headers declare a stored box in 40 of 40 runs (`bBox = 1`; also `bX = 1`, `bV = 1`, `ngtc = 2`; GROMACS 2020.1;
  `results/ABL-MYR-1b_tpr_header_diagnostic.json`, from the post-run script `tpr_header_diagnostic.py`). The box is in the
  files; the frozen reading path does not expose it.
- The frozen K1 requires a box for every run, and the bond-length check and the whole-chain Cα trace use it. Neither
  could therefore be applied, and the Cα RMSD to `system_em.gro` was not computed for any run. The atom-order validation
  did pass for 40 of 40 runs (full-system atom names identical to the text-parsed `system_em.gro`).
- GROMACS is not installed, so `gmx dump`/`gmx check` could not provide the independent reading.
- **Not done:** a second amendment of the validation path after this failure. The task excluded repeated changes after a
  failed criterion; a1 corrected a parser defect before any criterion was evaluated.

## What the check established and what it did not

| Question | Status |
|---|---|
| R1 usable production-start coordinates | **Not established** (K1). Positions are present but unvalidated |
| R2 production-start geometry; within- vs between-condition variation | **Not computed**: there are no validated whole Cα traces, so no pairs (all pair counts are 0). The reported "distinct starts = 10" is an artefact of zero pairs |
| R3 tc-grps content | **Verified from `index.ndx` against the tpr atom order**, both sets: solute group = 280 protein residues + ATP + MG (+ MYRP in Myr); solvent group = SOL + NA + CL; the two groups partition the system |
| R4 thermostat and coupling | **Comparable**: the only differences between conditions in `nvt.mdp`, `npt.mdp` and `production_md.mdp` are group names that add MYRP. ref_t 300 K for both groups everywhere. NVT: V-rescale, tau_t 0.5 ps, `gen_vel = yes`, `gen_seed = -1`, position restraints (`define = -DPOSRES`). NPT: V-rescale, Berendsen barostat, restraints. Production: V-rescale (Set1) or Nosé–Hoover (Set2), tau_t 1.0 ps, Parrinello–Rahman, `continuation = yes` |
| R5 composition | **Inherent differences only**: MYRP present only with Myr; SOL 29949 vs 29938 (11 waters, 0.037 %); NA 100, CL 87, ATP 1, MG 1 in both; net charge 4e-06 e in both |
| R6 classification | inherent to the Myr model: the group names, MYRP, the 11 waters; potential confounds: **none found** in parameters, groups or composition; production-start geometry: not separable from Myr presence during equilibration (and not measured here) |
| K5 replicates not copies | Passed **only through the velocity route**: 10 distinct velocity sets per set and condition, read by the same unvalidated reader. Consistent with `gen_vel = yes`, `gen_seed = -1` in each replicate's NVT step; not independently verified |

**What the evidence supports about shared starts:**
- **From ABL-MYR-1a:** both conditions derive from identical pre-minimization protein coordinates (Set1), and the Set2
  structures are copies of Set1's.
- **From this check:**
  - the replicates were equilibrated with the same, comparable thermostat and restraint settings, with random initial
    velocities stated;
  - the coupling groups and composition differ only in ways inherent to the Myr model.
- **Not supported yet:** what the validated production-start coordinates are, and how much they vary within and between
  conditions.

## Owner decision

1. **Install GROMACS** (for example the distribution package) so that `gmx dump`/`gmx check` read the box, coordinates and
   velocities independently. The frozen K1–K5 could then be re-evaluated as a dated amendment, a2, without changing any
   threshold. This is the remedy the check was designed for.
2. **Authorise amendment a2 without GROMACS:** read the box with the tpr-header routine of the same library
   (`extract_box_info`). It is less independent; the amendment would be written after this failure and must say so.
3. **Stop** ABL-MYR-1 here.

Any later ABL-MYR-1b protocol would still have to be written and frozen separately, before any trajectory is opened.

## Commands and resources

```bash
PY=benchmarks/audit/wpd_local_retrieval/build/venv/bin/python   # MDAnalysis interpreter
PYTHONPATH=python $PY benchmarks/audit/abl_myr_1b/run_comp_check.py run      # frozen check (attempt 1)
PYTHONPATH=python $PY benchmarks/audit/abl_myr_1b/run_comp_check_a1.py run   # amendment a1
PYTHONPATH=python $PY benchmarks/audit/abl_myr_1b/run_comp_check_a1.py report
```

Resources: attempt 1, 2.7 CPU-s; a1, 90.6 CPU-s and 89.7 s wall, 177 MB extracted (64 members), nothing downloaded.
Budgets were 1 GB extracted and 1800 CPU-s per run. The extractions are in the git-ignored `data_dryad/extracted_1b*/`.
