# ABL-MYR-1b: post-hoc comparability check, trajectory-comparison protocol, and its bounded comparison

**Closed 2026-10-25: the comparison result is final for this protocol** (no rerun, no threshold change, no frozen file
edited; provenance inventory and recovery steps in [`CLOSEOUT_2026-10-25.md`](CLOSEOUT_2026-10-25.md), D66).

**Two parts.** The comparability check (attempt 1, a1, a2) read run inputs and parameter files only: no trajectory was opened
and no dynamics scored in it. The post-hoc trajectory comparison was then run once (section "ABL-MYR-1b comparison"); it is
descriptive evidence on simulated ensembles. ABL-MYR-1 (stopped at G1) and ABL-MYR-1a (verdict V4) are preserved unchanged and
remain valid; their files are hash-recorded in the freeze records here and verified at every run. Nothing here is evidence
about Ayaz WT/M472I, HDX, drug response or biological function, and no causal Myr claim is made.

## Clarification (2026-10-25) of the amendment-a1 outcome

The a1 K1 outcome means that **production-start comparability was unverified, because MDAnalysis did not expose the stored
simulation box** (the run-input headers declare one) and GROMACS was not available to read it independently. It does **not**
demonstrate that the production starts differ, and it does not demonstrate that they match. Its label `not_comparable` is the
frozen name for "criteria not met". For the parameters, thermostat settings, coupling groups and composition that were
inspected, the differences detected between conditions are those expected from the presence of Myr (group names, the MYRP
molecule, 11 more waters), and **no unexpected difference was detected**. That is not proof that no confound exists: only the
listed items were inspected. Wording in frozen files that reads more strongly (for example "potential_confound: []" in the
generated reports, which means that no item matched the predeclared rules) is read with this qualification; see
[`CLARIFICATIONS_2026-10-25.md`](CLARIFICATIONS_2026-10-25.md). The a1 files are unchanged. Also, the a1 loader labelled
MDAnalysis's run-input coordinates as Å, while the reader returns nm; no a1 result used them numerically, and a2 compares both
readers in nm.

## Status

- **a2 (one bounded attempt, GROMACS 2024.5; timestamps in the records are UTC): criteria K1–K5 passed (`comparable_for_1b`).** This means the
  frozen criteria were met, not that the production starts are identical or that no confound exists.
- **The post-hoc ABL-MYR-1b trajectory-comparison protocol** ([`PROTOCOL_1B.md`](PROTOCOL_1B.md), `protocol_1b.json`, frozen
  2026-10-24T22:17:50Z UTC) was reviewed by an **independent AI scientific review** (not human peer review;
  [`REVIEW_2026-10-25_independent_AI.md`](REVIEW_2026-10-25_independent_AI.md)). A dated amendment settled one execution
  ambiguity (a bilateral G5 stop), and the wrapper was frozen before `PCA.zip` was opened.
- **The bounded comparison was then run once; all gates G0–G6 passed.** Frozen-rule outcome (post-hoc, descriptive): **baselines
  cover TopoFold**; no TopoFold-only residue survives the frozen descriptive rule. See "ABL-MYR-1b comparison" below.

## Sequence

| Step | Record | Result |
|---|---|---|
| Check frozen: `check_spec.json`, `run_comp_check.py`, `tests/test_abl_myr_1b_check.py` | `results/ABL-MYR-1b_check_freeze.json` (2026-10-24T21:44:11Z) | Questions, rules, constants, criteria K1–K5, stop conditions and budget fixed before any md.tpr coordinate, `nvt.mdp`, `npt.mdp` or `index.ndx` content was read |
| Attempt 1 (frozen code) | `results/ABL-MYR-1b_check_results_attempt1.json`, `…_report_attempt1.md`, `…_production_start_ca_attempt1.npz` | **Stopped at S3**: `duplicate index group 'MYRP'`. No md.tpr was loaded |
| Amendment a1 (dated): `check_amendment_a1.json`, `run_comp_check_a1.py`, `tests/test_abl_myr_1b_check_a1.py` | `results/ABL-MYR-1b_check_a1_freeze.json` (2026-10-24T21:46:49Z) | Implementation defect: GROMACS index files may repeat group names (here the `make_ndx` default residue-name groups). a1 accepts repeated names, treats a tc group as ambiguous if its occurrences differ, and otherwise imports the frozen check byte-identical |
| a1 run | `results/ABL-MYR-1b_check_a1_results.json`, `results/ABL-MYR-1b_check_a1_report.md` (generated) | **`not_comparable`: K1 failed** (0 of 40 run inputs validated); K2, K3, K4 and K5 passed |
| Amendment a2 (dated, one attempt): `check_amendment_a2.json`, `run_comp_check_a2.py`, `tests/test_abl_myr_1b_check_a2.py` | `results/ABL-MYR-1b_check_a2_freeze.json` (2026-10-24T22:12:26Z; GROMACS 2024.5 isolated install, acquisition in `results/ABL-MYR-1b_gromacs_acquisition.json`) | Changes only how the box, positions and velocities are read: `gmx dump` against the earlier reader and MDAnalysis's `extract_box_info`; K1–K5 and their thresholds unchanged; adds validation (d), reader agreement |
| a2 run | `results/ABL-MYR-1b_check_a2_results.json`, `results/ABL-MYR-1b_check_a2_report.md` (generated), `results/ABL-MYR-1b_production_start_ca_a2.npz` | **`comparable_for_1b`: K1–K5 all passed** |
| Protocol | `protocol_1b.json`, `PROTOCOL_1B.md` (generated by `write_protocol_1b.py` from the a2 results), `results/ABL-MYR-1b_protocol_freeze.json` (2026-10-24T22:17:50Z) | Frozen for independent review; no trajectory opened |
| Independent AI review | `REVIEW_2026-10-25_independent_AI.md` | Recorded with date and attribution, findings verbatim; not human peer review |
| Protocol amendment a1 (dated) | `protocol_amendment_a1.json`, erratum `protocol_amendment_a1_ERRATUM.md`, `results/ABL-MYR-1b_protocol_amendment_a1_freeze.json` (2026-10-24T22:34:10Z UTC) | Bilateral G5 stop; leave-one-run-out counts every retained run; wording per the review; αI rule. The erratum corrects one inaccurate phrase (the stop comes after extraction, before any frame is analysed) |
| Wrapper | `run_abl_myr_1b.py`, `tests/test_abl_myr_1b_wrapper.py`, `results/ABL-MYR-1b_wrapper_freeze.json` (2026-10-24T22:40:06Z UTC) | Imports the frozen runner unchanged; replaces only G4 (by G4b), the Set2 member list, G5 (bilateral) and the wording; frozen before `PCA.zip` was opened |
| Comparison run | `results/ABL-MYR-1b_results.json`, `…_report.md` (generated), `…_windows_Set{1,2}.csv`, `…_residue_summary.json` | All gates passed; outcome below |

## Why K1 failed in a1 (observed cause)

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
- **History:** at that point no further amendment of the validation path was made; a1 had corrected a parser defect before
  any criterion was evaluated. The owner then authorised one bounded amendment, a2 (below).

## Amendment a2: what was measured

- **Independent reader:** GROMACS 2024.5 (conda-forge, CPU, installed in user space; source, hashes and commands in the acquisition
  record) read the box, positions and velocities of all 40 production `md.tpr`. The earlier reader and GROMACS agree for 40 of
  40 files: max |Δx| 9e-07 nm, max |Δv| 2e-07 nm/ps (tolerances 0.001); GROMACS box vs MDAnalysis `extract_box_info` 4.9976e-05 nm
  (tolerance 1e-4). Velocities are present in both readers for 40 of 40 runs.
- **Validation:** atom order equal to `system_em.gro` in 40 of 40; protein bond lengths 0.9475–1.9623 Å (limit 0.9–2.3);
  Cα RMSD of each production start to its `system_em.gro` 0.265–0.309 Å (limit 1.0).
- **Production-start geometry** (Cα RMSD after optimal superposition, with the box from GROMACS and chains made whole):
  - Within Myr, within noMyr and between conditions, in both sets, every pair of runs differs by 0.296–0.372 Å (medians
    0.316–0.341 Å; heavy-atom medians 0.467–0.502 Å).
  - **No two starts are closer than 0.01 Å**: the 10 runs of each set and condition start from distinct equilibrated
    structures.
  - The within-condition and between-condition distributions overlap. The contrast statistic E of the analysis, applied to
    the starts (negative control): +0.0099 Å in Set1 and −0.0001 Å in Set2. Full table: `PROTOCOL_1B.md` and the a2 report.
- **Run-input settings** (read from the run inputs themselves, descriptive): ref_t 300 K for both groups, V-rescale (Set1) or
  Nosé–Hoover (Set2) in production; the box volumes, seeds and molecule blocks are in the a2 JSON. Each run has its own box,
  as expected from NPT equilibration, and its own random seed.
- **Resources:** 162.6 CPU-s, 160 s wall, 177 MB extracted (64 members), nothing downloaded. GROMACS acquisition: about
  80 MB downloaded, 14 s, 204 MB installed, no sudo.

## What the a1 check established and what it did not (kept as recorded; superseded for K1 by a2)

| Question | Status |
|---|---|
| R1 usable production-start coordinates | **Not established** (K1). Positions are present but unvalidated |
| R2 production-start geometry; within- vs between-condition variation | **Not computed**: there are no validated whole Cα traces, so no pairs (all pair counts are 0). The reported "distinct starts = 10" is an artefact of zero pairs |
| R3 tc-grps content | **Verified from `index.ndx` against the tpr atom order**, both sets: solute group = 280 protein residues + ATP + MG (+ MYRP in Myr); solvent group = SOL + NA + CL; the two groups partition the system |
| R4 thermostat and coupling | **No unexpected difference detected** in the inspected settings: the only differences between conditions in `nvt.mdp`, `npt.mdp` and `production_md.mdp` are group names that add MYRP. ref_t 300 K for both groups everywhere. NVT: V-rescale, tau_t 0.5 ps, `gen_vel = yes`, `gen_seed = -1`, position restraints (`define = -DPOSRES`). NPT: V-rescale, Berendsen barostat, restraints. Production: V-rescale (Set1) or Nosé–Hoover (Set2), tau_t 1.0 ps, Parrinello–Rahman, `continuation = yes` |
| R5 composition | **Only differences expected from the Myr model detected**: MYRP present only with Myr; SOL 29949 vs 29938 (11 waters, 0.037 %); NA 100, CL 87, ATP 1, MG 1 in both; net charge 4e-06 e in both |
| R6 classification | expected from the Myr model: the group names, MYRP, the 11 waters; **no unexpected difference detected** in the inspected parameters, groups and composition (not proof that no confound exists); production-start geometry: not separable from Myr presence during equilibration (and not measured here) |
| K5 replicates not copies | Passed **only through the velocity route**: 10 distinct velocity sets per set and condition, read by the same unvalidated reader. Consistent with `gen_vel = yes`, `gen_seed = -1` in each replicate's NVT step; not independently verified |

**What the evidence supports about shared starts:**
- **From ABL-MYR-1a:** both conditions derive from identical pre-minimization protein coordinates (Set1), and the Set2
  structures are copies of Set1's.
- **From this check:**
  - the replicates were equilibrated with the same, comparable thermostat and restraint settings, with random initial
    velocities stated;
  - the inspected coupling groups and composition show no unexpected difference (the differences found follow from the presence of Myr); this is not proof of no confound.
- **Not supported yet:** what the validated production-start coordinates are, and how much they vary within and between
  conditions.

## ABL-MYR-1b comparison (run once, 2026-10-24 22:40 UTC start; all numbers read from `results/ABL-MYR-1b_results.json`)

Post-hoc, descriptive evidence from simulated Abl kinase-domain ensembles (de Buhr and Gräter, Dryad). It is separate from Ayaz
WT/M472I, HDX, drug response and biological function, and no causal, clinical or novelty claim is made.

- **Gates:** G0 archives and all 17 lineage hashes verified; G1 107 members; G2 topology `KD_noIa_bb_ref.pdb`, 753 atoms
  (N, CA, C of 251 residues), φ/ψ comparator available; G3 identity 1.0000 to P00519-2, reference span 254–504, no unobserved
  range; G4b evidence verified (it does not establish identical production starts); G5 40 of 40 runs retained, 501 analysis
  frames each (files cover 0–2000 ns, window from 1000 ns), both sets analysed; G6 pbc `none`, worst excluded-window fraction
  0.004 (terminal φ/ψ).
- **αI:** residues 518–531 are **not represented** (the topology ends at residue 504). The analysis covers the represented
  kinase-domain residues only; no claim is made about the αI helix or the Myr-binding region.
- **Outcome (frozen rule, restricted label): "baselines cover TopoFold (post-hoc, descriptive)", step 4.** The windows that are
  robust in both thermostat sets (E > 0 in the full set and in every leave-one-run-out recomputation; 20 recomputations per set)
  number κ 13, τ 6, d(i,i+2) 12, d(i,i+3) 13, RMSF 3 and, for the stronger-information comparator, φ 10, ψ 12.
  - TopoFold residues (R_T): 34 residues in 6 segments (isoform 1b numbering: 277–279, 346–348, 358–365, 445–447, 452–454,
    478–491).
  - Equal-information Cα baseline residues (R_B): 40 residues in 9 segments.
  - **TopoFold residues outside the baselines' ±2-residue neighbourhoods: none (0).** No TopoFold-only residue survives the
    frozen descriptive rule. The baselines flag one residue with no TopoFold window nearby (317).
  - φ/ψ comparator residues, reported separately: 15.
- **Per set:** the per-set counts of robust windows, run-to-run variability, temporal diagnostics (lag-1 autocorrelation and
  decorrelation within 50 frames) and the production-start descriptors of the 20 retained runs per set are in the report.
  The median E of every channel is negative in both sets: typical windows show no cross-condition separation beyond
  run-to-run variability.
- **Resources:** about 16.5 CPU-min (wall 15 min 46 s), 293 MB extracted, nothing downloaded; budget 2 CPU-h.

**Limits.** E > 0, leave-one-run-out stability and the three-residue rule are descriptive filters with no error-rate control or
significance. The thermostat sets share their post-minimization structures, so agreement between them is a consistency check,
not independent replication. Backbone frames every 2 ns support a geometry description, not kinetics; temporal dependence is
shown in the report. The production-start Cα RMSD diagnostic is not a calibrated control for E. One modelled kinase domain
(residues 254–504), one force field. The comparability criteria passed, which does not rule out unexamined confounds.

## Owner decisions

1. **Decided 2026-10-25 by the owner: the result is final for this protocol** and is preserved as a negative outcome for
   TopoFold-specific signal (D66). Any later human review or new analysis is a separate, newly frozen study that does not
   replace this record.
2. **No further amendment is open.** The wrapper did not need a workaround; no gate failed.
3. The separate Ayaz M472I route (DESRES follow-up, HDX request) is unaffected by this study and remains pending.

## Commands and resources

```bash
PY=benchmarks/audit/wpd_local_retrieval/build/venv/bin/python   # MDAnalysis interpreter
GMX=<isolated GROMACS 2024.5>/env/bin/gmx                        # see results/ABL-MYR-1b_gromacs_acquisition.json
PYTHONPATH=python $PY benchmarks/audit/abl_myr_1b/run_comp_check.py run          # frozen check (attempt 1)
PYTHONPATH=python $PY benchmarks/audit/abl_myr_1b/run_comp_check_a1.py run       # amendment a1
TOPOFOLD_GMX=$GMX PYTHONPATH=python $PY benchmarks/audit/abl_myr_1b/run_comp_check_a2.py run      # amendment a2
python3 benchmarks/audit/abl_myr_1b/write_protocol_1b.py write                   # protocol from the a2 results
PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1b/run_abl_myr_1b.py run --adapter-python $PY   # comparison (once)
PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1b/run_abl_myr_1b.py report
python3 benchmarks/audit/abl_myr_1b/summarize_results_1b.py                      # residue numbering summary
```

Resources of the check: attempt 1, 2.7 CPU-s; a1, 90.6 CPU-s; a2, 162.6 CPU-s; each extracted 177 MB; GROMACS acquisition about 80 MB
and 14 s, no sudo. The comparison used about 16.5 CPU-min and 293 MB; extractions are in the git-ignored `data_dryad/extracted*/`.

Closeout inventory (hashes only; nothing rerun): `PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1b/closeout_inventory.py --copy-dir ~/Downloads`,
then `sha256sum -c benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_closeout_SHA256SUMS.txt`.
