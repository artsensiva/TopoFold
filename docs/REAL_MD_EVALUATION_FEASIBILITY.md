# Real-MD Evaluation Feasibility

**Status:** DRAFT for independent protocol review. No analysis has been run and no trajectory
has been downloaded.
**Baseline:** `295d739a783cc1d178131c4c61293042a0c56f7d` (branch `scientific-rebuild`)
**Checked:** 2026-10-04

Only metadata and small files were retrieved:
- repository APIs;
- the two single-frame `.gro` files of candidate A (13 MB in total), into a scratch directory
  outside the repository;
- the first 4 MB of each candidate A trajectory, to read frame headers.

Each fact below states how it was verified.

## 1. Task definition

**Task: backbone-state localization, not pocket detection.** Given an MD ensemble of one protein,
each method ranks residues (or residue windows) by evidence that their local backbone geometry
changes between conformational states. No method is told where the change is.

Pocket detection would need cavity analysis and pocket ground truth, and it is out of scope. A
residue that changes backbone state is not thereby part of a pocket.

## 2. Candidates (two, both *E. coli* adenylate kinase)

### Reference structures (label source for both candidates)

Verified via the RCSB REST API:

| Entry | Title | Organism | Length | Chains | Resolution | Citation |
|-------|-------|----------|--------|--------|------------|----------|
| 4AKE | "ADENYLATE KINASE" | *E. coli* | 214 | 2 | 2.2 Å | Structure 1996, doi:10.1016/S0969-2126(96)00018-4 |
| 1AKE | Complex with the inhibitor Ap5A | *E. coli* | 214 | 2 | 2.0 Å | J. Mol. Biol. 1992, doi:10.1016/0022-2836(92)90582-5 |

The two canonical sequences are identical. 4AKE is used as the open reference and 1AKE as the
closed reference, as in both candidate sources.

### Candidate A: Orädd & Andersson, unrestrained MD

**Source.** Zenodo [10.5281/zenodo.5583119](https://doi.org/10.5281/zenodo.5583119), CC BY 4.0,
published 2021-10-20. Paper: Orädd et al., *Sci. Adv.* 7, eabi5514 (2021),
doi:10.1126/sciadv.abi5514. Methods read via
[PMC8597995](https://pmc.ncbi.nlm.nih.gov/articles/PMC8597995/).

| Property | Value | How verified |
|----------|-------|--------------|
| Simulations | Two unrestrained trajectories: 450 ns starting from 4AKE (open) and 335 ns starting from 1AKE (closed); heavy-atom restraints only during equilibration | Paper Methods |
| Engine and conditions | CHARMM36, GROMACS 2016.4, 303.15 K, Nose–Hoover thermostat, Parrinello–Rahman **semi-isotropic** barostat (unusual for a soluble protein; rationale not found) | Paper Methods |
| Files | `ORADD_open_state_traj.xtc` 780,301,896 B; `ORADD_closed_state_traj.xtc` 580,049,352 B; `ORADD_{open,closed}_state.gro` ≈ 6.57 MB each; ≈ 1.37 GB in total, with md5 checksums listed | Zenodo API |
| System contents | One chain of 214 residues (numbered 1–214), TIP3 water, K⁺ (`POT`). **No Ap5A or any other ligand in either system**, so the closed-start simulation is apo. 95,175 / 95,196 atoms; initial cubic box 10 nm | Downloaded `.gro` files; md5 matches Zenodo |
| Sequence | Identical to 4AKE/1AKE, no mutations | `.gro` Cα residues vs RCSB sequence |
| Saved data | All atoms. Frame interval 200 ps (frames at 0, 200, 400 ps in both files). xtc precision 1000/nm, i.e. 0.01 Å quantization | xtc frame headers from the first 4 MB |
| Frame count | ≈ 2,250 (open start) and ≈ 1,675 (closed start). Estimated from file size / compressed frame size and consistent with 450 / 335 ns at 200 ps; **not yet verified** | Estimate |

### Candidate B: Seyler et al., transition path ensembles

**Source.** figshare [10.6084/m9.figshare.7165306.v2](https://doi.org/10.6084/m9.figshare.7165306.v2),
CC BY 4.0. Data from Seyler, Kumar, Thorpe & Beckstein, *PLoS Comput. Biol.* 11(1) (2015).
Verified via the figshare API and the
[MDAnalysisData page](https://www.mdanalysis.org/MDAnalysisData/adk_transitions.html).

| Property | Value |
|----------|-------|
| Contents | 200 DIMS MD paths (`DIMS.tar.gz`, 731,492,423 B; 3,341 atoms; `adk4ake.psf`) and 200 FRODA paths (`FRODA.tar.gz`, 521,237,730 B; 1,656 atoms, no hydrogens) |
| Transition | Apo *E. coli* AdK, closed (1AKE) → open (4AKE) |
| DIMS conditions (primary source) | CHARMM22/CMAP, **implicit solvent (ACE in CHARMM)**, CHARMM c36b2, 300 K. Progress variable: heavy-atom RMSD to the target. Start 1AKE chain A, target 4AKE chain A. 200 paths. Source: Seyler et al. 2015, Methods, [doi:10.1371/journal.pcbi.1004568](https://doi.org/10.1371/journal.pcbi.1004568) |
| FRODA conditions (primary source) | Geometric targeting algorithm on heavy atoms, with no explicit solvent or force field. Progress variable: Cα RMSD to the target. Same start/target |
| Nature | **Biased path sampling**, not equilibrium MD |
| Licence discrepancy | figshare record: CC BY 4.0. The paper's data statement points to the PSAnalysisTutorial repository (Zenodo 10.5281/zenodo.31457) under GNU GPL v3. Which terms govern which copy is unresolved |

## 3. Assessment

**Candidate B is not suitable for evaluation.** Its paths are driven toward the same crystal
structure that defines the labels, so a method recovering the 1AKE→4AKE differences would partly
recover the bias. It has no equilibrium populations or kinetics. At most it could serve as a
plumbing test (reading files, running every method end to end), with no evaluative claim.

**Candidate A is conditionally suitable for a single-system development pilot.**

Strengths:
- unbiased production MD;
- open licence (CC BY 4.0 permits commercial use with attribution);
- all atoms, so φ/ψ baselines are possible;
- the label source (crystal structures) is independent of both TopoFold and the trajectories.

Weaknesses:
- one trajectory per starting state, so no replicate variance;
- 450 / 335 ns may not contain open↔closed transitions;
- the closed-start run lacks the ligand that stabilizes 1AKE;
- a single protein gives no held-out evaluation.

**Recommendation.** Use candidate A only as a development dataset, and only if the feasibility
gate G1 (§6) passes. No evaluation claim can rest on it. A held-out system is not identified;
it must be chosen and frozen before any development results are examined.

Excluded: D. E. Shaw Research trajectories (including BPTI). According to secondary sources,
access is by request and for non-commercial use; this was not verified with DESRES. Per
instructions, BPTI/Met52 is not assumed to be ground truth.

## 4. Observable and labels (to be frozen before any MD analysis)

### Observables

- **Per-trajectory:** residue/window scores computed separately on each trajectory.
- **Pooled:** scores computed on frames from both trajectories, with no origin label given to any
  method. Any pooled signal partly reflects the different starting structures. This is a
  construction effect and must be reported as such.

Neither design can establish kinetics. Pooling two trajectories launched from different crystal
states mixes them by construction, so it cannot show that transitions occur spontaneously or
estimate their rates. A single trajectory per start can at most record whether transitions
occurred in that run; it cannot estimate rates or equilibrium populations.

### Independent residue labels from the crystal structures

No residue range is assumed. Two label definitions are pre-specified; both are reported.

- **L1 (backbone dihedrals):** residue i is *changed* if the circular difference in φ or ψ
  between states exceeds the within-state chain-A vs chain-B difference by a pre-registered
  margin.
- **L2 (local Cα geometry):** residue i is *changed* if the window-superposed Cα RMSD between
  states exceeds the within-state chain-to-chain value by a pre-registered margin.

The margin and the window must be fixed before MD scores are computed. Each label definition
favours baselines that share its descriptor; this bias is reported.

**Unlabelled residues are not verified negatives.** A residue without a crystal-state difference
may still move in solution. Results are reported as positive-unlabelled:
- recall of labelled positives;
- *labelled-positive* precision at k, i.e. the fraction of top-k residues that are labelled
  positives.

The second quantity is not the true precision: unlabelled residues among the top k may be true
positives, so the true precision is unknown. No specificity, false-positive rate or ROC-based
specificity claim is made.

Domain annotations from the literature (e.g. LID/NMP) may be cited only after they are verified
against the primary structure papers. They are never substituted for L1/L2.

### Frame-level state (needed only for gate G1)

Cα RMSD to 4AKE and to 1AKE after superposition on a pre-registered core. This label is
independent of TopoFold. It favours Cartesian baselines, which is disclosed.

## 5. Equal-information baselines and statistics

**Equal information**
- All methods get the same frames, the same 214 residues, and no prior region.
- Window-based methods use the same window sizes, from a pre-registered grid.
- Primary comparison is Cα-only for every method. Full-backbone baselines (φ/ψ) are a secondary
  comparison, labelled as having more information.

**Baselines (strongest reasonable)**
- per-residue φ/ψ as (sin, cos), scored by the same multimodality/variability statistic as
  TopoFold;
- windowed Cα–Cα distance features;
- window-superposed Cα RMSD to the ensemble mean;
- RMSF after global superposition;
- TICA on Cα distances, across a pre-registered lag grid (lag dependence reported, not tuned).

**Statistics**
- Frames 200 ps apart are autocorrelated.
- Uncertainty comes from a block bootstrap over contiguous time blocks, with a pre-registered
  block-length grid (e.g. 1, 5, 10, 25 ns), plus an effective-sample-size estimate.
- With one trajectory per start, between-trajectory variance cannot be estimated, and no
  significance claim is made across trajectories.

**Coordinate precision.** The xtc files store coordinates rounded to 0.01 Å (±0.005 Å per
coordinate). This adds an input error to every descriptor, but the resulting angular error
depends on geometry: for torsions it grows as 1/sin κ near collinear triples, and it is
correlated across descriptors that share atoms. 0.005·√3 / 3.8 ≈ 2e-3 rad is an order of
magnitude for typical, non-collinear geometry only. **It is not a universal bound.** The actual
effect must be measured on these trajectories, for example by recomputing descriptors after
re-quantizing the coordinates, and reported with any small-variance result.

M5 (Fréchet distance) is not used on this path and is not a prerequisite.

## 6. Development vs held-out, and feasibility gates

| Gate | Action | Condition |
|------|--------|-----------|
| **G0 (data)** | Download candidate A (≈ 1.37 GB, needs authorization). Verify md5 against Zenodo and record frame counts and times | Integrity |
| **G1 (states)** | Compute per-frame Cα RMSD to 4AKE and 1AKE, independent of TopoFold. Record whether either trajectory visits the other state. Determine which crystal chain was used as the start by comparing frame 0 with each chain | If neither trajectory changes state, per-trajectory localization has no expected positive signal and only the pooled construction remains. This limits what AdK can test. It is not evidence of engine failure. Report it rather than proceed silently |
| **G2 (labels)** | Freeze L1/L2 from the crystal structures, before any TopoFold or baseline score is computed on MD | — |
| **Held-out** | A second, independent system with an open licence, unbiased MD and an experimental two-state label source must be selected and frozen before AdK results are seen. **Not yet identified** | — |

Splitting an AdK trajectory in time is not a held-out evaluation.

### WP-B execution specification (after review and download authorization)

**Data requirements**
- Exactly the four candidate A files, about 1.37 GB:
  - `ORADD_open_state_traj.xtc` (780,301,896 B);
  - `ORADD_closed_state_traj.xtc` (580,049,352 B);
  - two `.gro` files (6,567,118 and 6,568,567 B).
- 4AKE and 1AKE coordinate files from RCSB.

**Provenance checks**
- md5 against the Zenodo record.
- SHA-256 recorded locally.
- Frame count, first and last time, and time step read from the files, compared with the paper's
  450 / 335 ns.
- All of this is recorded in a provenance results file, e.g. `results/adk_data_provenance.json`.
  It is not labelled SYNTHETIC, because these are third-party MD data; it states the source DOI.

**Preprocessing**
- TopoFold has no XTC reader. Use one pinned external reader (MDAnalysis or MDTraj; version
  recorded) to extract Cα and backbone (N, Cα, C, O) coordinates to float64 arrays, with residue
  numbering 1–214.
- Hash the arrays.
- Check for chain breaks (none expected: 214 contiguous residues; verify).
- Waters and ions are dropped.
- Raw files are not modified.

**Independent state assessment (G1)**
- Per-frame Cα RMSD to 4AKE and 1AKE (each chain) after superposition on a pre-registered core
  definition.
- Time series and histograms are stored. No TopoFold or baseline scores are computed in WP-B.

**Proceed / stop**
- **Proceed** to a separately reviewed scoring protocol only if G0 passes and G1 is documented.
  A no-transition outcome restricts the design to the pooled construction.
- **Stop** and report if the md5 checks fail, the frame count disagrees with the paper beyond
  one frame per trajectory, or the trajectories cannot be read reproducibly.

### Held-out selection without looking at performance

1. Before any AdK scoring, write `held_out_selection.md` with the eligibility criteria and hash
   it. Criteria:
   - open licence;
   - unbiased MD of total length stated;
   - experimental structures in at least two states for the same sequence;
   - no overlap with AdK;
   - no prior TopoFold analysis.
2. Enumerate candidates from the literature and databases, recording each source.
3. A reviewer checks eligibility from metadata only.
4. If several candidates are eligible, choose by a recorded seeded random draw.
5. No TopoFold or baseline computation is run on any candidate before the choice is frozen.

## 7. Unresolved

1. Exact frame counts and file integrity (G0).
2. Whether open↔closed transitions occur (G1).
3. The reason for the semi-isotropic barostat.
4. Which crystal chain started each run.
5. The effect of simulating 1AKE without Ap5A.
6. No replicates.
7. No held-out system.
8. Candidate B licence terms (CC BY 4.0 on figshare vs GPL v3 in the paper's data statement).
9. DESRES access terms not verified first-hand.
10. Candidate A water model as stated by the authors (the `.gro` files show residue name `TIP3`).
11. Exact frame counts remain estimates until G0.

## 8. Inputs and decisions still needed to start WP-B gates G0–G1

WP-B is independent of WP-A and of the ordered-representation candidate (WP-C). It runs
whether or not that candidate is retained.

1. **Authorization** to download candidate A, ≈ 1.37 GB (§6).
2. **XTC reader.** Choose MDAnalysis or MDTraj, pin its version, and approve adding it to an
   isolated WP-B environment. TopoFold has no XTC reader, and adding one to the project
   environment is a dependency change.
3. **Storage location** for raw and derived arrays, with a narrow ignore rule like the WP-A
   `build/` rule. Raw files are never modified.
4. **RMSD core definition for G1,** pre-registered before any frame is inspected: the residue
   set used for superposition, and which 4AKE/1AKE chains serve as references.
5. **Recording format for G1 outputs** (per-frame RMSD to each reference, time series). Real MD
   data, so not labelled SYNTHETIC; the source DOI and CC BY 4.0 attribution are recorded.
6. **For G2, not needed for G0–G1:** L1/L2 label margins and window sizes, fixed before any
   representation is scored.
7. **Held-out eligibility file** (§6), written before any development scoring.

No TopoFold or baseline scoring takes place in G0–G1.

## 9. Execution of G0–G1

G0–G1 are executed in `benchmarks/audit/adk_feasibility/` under a specification frozen before
any trajectory inspection (`SPEC.md`; hashes in `results/adk_spec_freeze.json`). That
specification supersedes one rule in §6. A frame count that differs from the paper's stated
duration is now a **recorded discrepancy**, not a stop condition: md5 establishes file
integrity, and the stated durations are source statements. Results: `REPORT.md` in that
directory.
