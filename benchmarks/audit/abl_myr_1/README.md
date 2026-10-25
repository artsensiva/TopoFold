# ABL-MYR-1: Myr vs noMyr kinase-domain ensembles from public replicated MD (spec and runner frozen; original run BLOCKED at G1)

**Status (2026-10-24): specification and runner frozen; the original run passed G0 and stopped at G1; nothing was
computed.** All three required archives are present and verified, but both model archives lack the required member
`Set2/system.gro`, as published (details below). Descriptive real-MD evidence only once run, never biological
validation.

Later (kept separate from this original run): the post-hoc ABL-MYR-1a audit and ABL-MYR-1b check and comparison
(`../abl_myr_1a/`, `../abl_myr_1b/`); ABL-MYR-1b is closed (2026-10-25, D66).

## Original run (2026-10-24T20:59:43Z, frozen runner, unmodified)

| Gate | Result |
|---|---|
| G0 | **Passed.** Verified twice: with `sha256sum` against the values in `run_spec.json` (sizes and SHA-256 of `PCA.zip`, `kinase_kinked_Myr.zip`, `kinase_kinked_noMyr.zip` all match), and by the runner |
| G1 | **Failed.** `member 'Set2/system.gro': 0 matches in the archive (exactly one required)`. Record: `results/ABL-MYR-1_blocked.json` (frozen runner and spec hashes inside) |
| G2–G6 | Not reached; nothing was extracted, no frame read, no statistic computed |

- **Observed cause, in the data:** `results/ABL-MYR-1_g1_member_audit.json` (from `g1_member_audit.py`, which reads zip
  directory listings only): `PCA.zip` has 43 of 43 required members; each model archive has 31 of 32, missing only
  `Set2/system.gro`. The dataset README (SHA-256 equals Dryad's listed digest) lists `system.gro`, `system_box.gro`,
  `system_em.gro`, `system_solv.gro` and `system_solv_ions.gro` for `Set1`, but only `system_em.gro` and `topol.top`
  for `Set2`, in both models. So the member is absent in the published dataset, not through a bad download.
- **Not an implementation defect:** the runner behaved as frozen ("a missing member stops the run, no substitute"). No code
  was changed and no amendment exists. `kinase_straight_noMyr.zip` and the other archives were not used.
- **My error in the specification:** the frozen spec requires `system.gro` in both sets. The README tree available at intake
  showed that `Set2` lacks it, and I did not check the member list file by file before freezing.
- **Why `system_em.gro` is not a free substitute:** it is the energy-minimised structure of each system, so protein
  coordinates can differ between the Myr and noMyr systems after minimisation. The frozen 0.001 nm identity rule was
  defined for `system.gro`; using another file would change the test and needs a dated amendment.

## Owner decision needed

Any continuation is a **dated post-hoc amendment** in a separate directory (the original attempt stays as recorded),
decided before any archive member is extracted. Options, none started:

1. **Amend G4 for Set2:** check `topol.top` and `production_md.mdp` identity (Set2 has both) and state plainly that
   starting-coordinate identity is not verified for Set2. Set1 keeps the frozen test. Weaker evidence of a shared start for
   Set2.
2. **Analyse Set1 only:** the frozen "reproduced" rule needs both sets, so the outcome could not exceed Set1 evidence.
3. **Stop ABL-MYR-1.**

- **Data:** de Buhr & Gräter, Dryad [10.5061/dryad.9cnp5hqnx](https://doi.org/10.5061/dryad.9cnp5hqnx) (CC0), study
  eLife [10.7554/eLife.85216](https://doi.org/10.7554/eLife.85216). Intake record:
  [`../debuhr_abl_intake/intake_record.json`](../debuhr_abl_intake/intake_record.json) (conditional-go).
- **Separate question.** This is the myristoyl/allostery question in a modelled isoform-1b kinase domain with
  ATP/Mg²⁺. It is unrelated to the Ayaz M472I/imatinib-resistance question and its apo HDX endpoint, and no result
  here is evidence for the ABL WT vs M472I preregistration.
- **Freeze:** `run_spec.json` was frozen in `results/ABL-MYR-1_freeze.json` (2026-10-24T19:54:07Z, SHA-256
  `f8825840…`) before any trajectory, archive or topology was downloaded (`input_files_present_at_freeze: []`).
- **Code freeze:** the runner, its tests and the library modules it uses are hash-recorded in
  `results/ABL-MYR-1_code_freeze.json` (2026-10-24T20:32:08Z), before any study archive was opened. The runner refuses to
  run if any of these hashes, the spec or the configuration differs. Files that had already appeared in `data_dryad/` are
  listed there by name, size and time only; none was hashed, read or opened.

## Why this subset

| Item | Finding (from the Dryad file list, the rendered README and the eLife Methods) |
|---|---|
| Smallest matched, replicated comparison | `kinase_kinked_Myr` vs `kinase_kinked_noMyr`: the same kinase-only model with a kinked αI helix, with or without Myr |
| Replicates | 10 runs per condition and thermostat set: Set1 (v-rescale) replicates 11–20; Set2 (Nosé–Hoover) 21–30. 2000 ns each, of which the paper analyses 1000–2000 ns |
| Conditions | GROMACS 2018.5/2020.3, CHARMM36m (March 2019), TIP3P, 150 mM NaCl, 300 K; ATP and Mg²⁺ bound; models built from 2FO0 (D355N reverted, inhibitor replaced by ATP) |
| Trajectories TopoFold can read | `PCA.zip` (875 MB) holds `traj_{model}_{rep}.xtc` for all 7 models × 20 runs (140 files): kinase-domain backbone, one frame per 2 ns. Topology candidate: `KD_noIa_bb_ref.pdb` (backbone reference; whether αI is included is unverified) |
| Model archives (`kinase_*.zip`, 123–133 MB) | **Simulation inputs and analysis outputs only:** `.mdp`, `.top`/`.itp`, starting `.gro`, per-replicate `.tpr`, `.xvg`/`.dat`. **No trajectory files.** They are needed only to verify the shared start (gate G4) |
| Larger alternatives | `FDA.zip` (26.3 GB): protein-only XTC every 200 ps, first 1000 ns removed; outside the budget |
| Numbering | isoform 1b (paper); verified 1b = 1a + 19 on UniProt P00519-2 (K290, E305, H380, D400, M491) |
| Reuse terms | CC0-1.0 (Dryad record); paper CC BY |
| Size needed | 1,139,923,203 bytes (PCA.zip + both kinked model archives) |

## Frozen design (`run_spec.json`)

- **Gates:**
  - **G0:** access and SHA-256 against the Dryad digests;
  - **G1:** archive members;
  - **G2:** topology and atom count;
  - **G3:** mapping to P00519-2 with `kinase_roles` under declared isoform-1b numbering;
  - **G4:** identical starting protein coordinates and run parameters for Myr and noMyr;
  - **G5:** analysis frames (t ≥ 1000 ns, 2-ns interval, ≥ 450 frames per run);
  - **G6:** breaks and exclusions (STOP above 5 % of windows).
- **Method:** the existing `topofold_inputs.ensemble` analysis. Whole runs are the replicates. Set1 is primary, Set2 an
  independent replication; the sets are never pooled.
- **Channels:**
  - TopoFold κ, τ;
  - equal-information Cα baselines: d(i,i+2), d(i,i+3) and per-residue Cα RMSF;
  - φ/ψ as a separately labelled stronger-information comparator.
- **Statistic:** per window, E = cross-condition run-pair distance minus within-condition run-pair distance (medians).
  A window is reproduced if E > 0 in the full set and in all 20 leave-one-run-out subsets, in both sets.
- **Outcomes**, ordered: blocked; inconclusive; TopoFold adds reproduced residues (≥ 3 outside the baseline residues
  ± 2); baselines cover TopoFold; partial. Temporal diagnostics are reported. There is no inference (D33).
- **Budget:** ≤ 1.2 GB download, ≤ 4 GB storage, ≤ 2 CPU-h.

## Engineering preparation (SYNTHETIC tests only; no scientific result)

- **Runner:** `run_abl_myr.py` has the subcommands `freeze-code`, `run` and `report`. `run` applies gates G0–G6 exactly as
  frozen, then the comparison. A failed gate stops it and writes `ABL-MYR-1_blocked.json` with the gate and the observed
  cause. The orchestrator needs the TopoFold core; a `_worker` subcommand runs in the MDAnalysis interpreter, because no one
  environment has both.
- **Interpretations fixed before any study file was opened** (the spec text is ambiguous on them; listed in the runner
  docstring): members are found by path suffix with exactly one match; "backbone" atom names; "differ only by the Myr
  molecule" means every non-Myr line of `topol.top`, including water and ion counts, is identical; a bad time stamp or
  interval stops the run while only a run with fewer than 450 analysis frames is excluded; the spec's pbc `make_whole` is
  the adapter's `make_whole_sequential`.
- **Tests:** `tests/test_abl_myr_runner.py` (16 pure tests that run everywhere, plus 4 end-to-end tests on SYNTHETIC archives
  that need the core and `TOPOFOLD_ADAPTER_PYTHON`) and `tests/test_abl_myr_runner_mda.py` (6 worker tests in the CI adapter
  job). They cover archive-member selection, traversal-safe extraction, topology-to-trajectory atom counts, frame time
  stamps, residue mapping, shared-start checks, exclusion reporting, φ/ψ against an independent formula, RMSF against an
  independent superposition, and the G4 and G6 stops.
- **Issues the synthetic runs found before the freeze:** the adapter's make-whole policy has a different name from the spec's;
  a `.gro` difference of exactly the tolerance failed through floating-point noise; and the fixtures themselves needed
  intact chain geometry (noise large enough to break Cα–Cα distances correctly drives the run into G6).

## Present in `data_dryad/`

At the code freeze only metadata was listed (`PCA.zip` and `kinase_kinked_noMyr.zip` were absent). Later all three required
archives arrived and passed G0 (see the original run above). Also present but not part of the specification and not used:
`README.md` (documentation, read for the G1 audit), `kinase_active.zip`, `kinase_straight_Myr.zip`,
`kinase_straight_noMyr.zip` and `SH3-SH2-KD_Myr.zip`.

## Command of the original run

The command that produced the original attempt (from the repository root; it would stop at the same G1 gate again):

```bash
PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1/run_abl_myr.py run \
    --adapter-python benchmarks/audit/wpd_local_retrieval/build/venv/bin/python
```

`.venv` has the TopoFold core; the second interpreter has MDAnalysis. G0 verifies the three archives against the frozen
SHA-256 values before anything is extracted. The `report` command only renders a completed `ABL-MYR-1_results.json`; there is none.

## Blocker and remedy

- **Earlier blocker (G0, observed 2026-10-24; resolved when the owner downloaded the files):** the Dryad file API returned
  HTTP 401 ("must have current bearer token") and the web route served a proof-of-work bot challenge. Neither was bypassed.
- **Remedy (owner), either:**
  - download `PCA.zip`, `kinase_kinked_Myr.zip` and `kinase_kinked_noMyr.zip` from the dataset page in a browser into
    `benchmarks/audit/abl_myr_1/data_dryad/`; or
  - provide a Dryad API token through an environment variable, not in the repository.
- **After that:** the runner is written, code-frozen and executed against the frozen gates.

## Figshare single-run set (loader check only)

[10.6084/m9.figshare.4496795](https://doi.org/10.6084/m9.figshare.4496795) (CC BY 4.0): one 1-µs run of ABL1 242–493
(AMBER99SB-ILDN, 300 K, 4000 frames at 250 ps) in MDTraj HDF5 (141,949,663 B).
- **Not downloaded:** reading it needs an HDF5/MDTraj path that the Dryad XTC route does not use, so it would not
  reduce a risk of this path.
- **Never replicated evidence:** it is a single run.
