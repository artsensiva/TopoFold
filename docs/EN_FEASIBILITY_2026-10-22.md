# EN Feasibility: Real-Data Ensemble-Comparison Candidates (2026-10-22)

**Status:** feasibility only (roadmap §12.3, workstream EN). No data were downloaded and nothing was run, implemented or
requested. C1 remains blocked and resumable; this note does not replace it. ABL remains conditional on provider
permission and data. *Revised 2026-10-22: permissions separated by provider (§3); GPCRmd frame counts and the
permutation-test limit corrected; access-request draft (`ABL_DESRES_DATA_REQUEST_DRAFT.md`) and evaluation checklist (§5) added.*

**Sources** (read 2026-10-22, metadata and text only):
- Ayaz P, Lyczek A, Paung Y, … Shan Y, Shaw DE. *Structural mechanism of a drug-binding process involving a large
  conformational change of the protein target.* Nat Commun 14, 1885 (2023), doi:10.1038/s41467-023-36956-5,
  PMC10076256 (article CC BY).
- Saleh T, Rossi P, Kalodimos CG. *Atomic view of the energy landscape in the allosteric regulation of Abl kinase.*
  Nat Struct Mol Biol 24, 893–901 (2017), doi:10.1038/nsmb.3470, PMC5745040. This is an NMR study and supplies **no MD
  trajectories**; the analyst assessment's pairing of it with MD data is not supported.
- Vögele M, Thomson NJ, Truong ST, McAvity J, Zachariae U, Dror RO. *Systematic analysis of biomolecular
  conformational ensembles with PENSA.* J Chem Phys 162, 014101 (2025), arXiv:2212.02714. Code: github.com/drorlab/pensa
  (MIT licence).
- GPCRmd simulation reports for dynamics IDs 151 and 169, the GPCRmd legal notice and data-download documentation,
  and the PENSA tutorial (`pensa.readthedocs.io`, "tut-1-intro").

Labels as in roadmap §12: [O] observed in a source, [D] deduction, [A] assumption, [P] proposal.

## 1. Candidates

### ABL — imatinib-resistance mutants vs wild type (Ayaz et al. 2023)

- **Question [P]:** do Abl kinase-domain ensembles of the resistance mutants M472I and F486S differ from wild type in
  C-lobe and activation-loop geometry, and do the residues that differ coincide with the experimentally measured HDX-MS
  differences (M472I vs WT)?
- **Permitted claim:** a level-B statement that a method locates replicated ensemble differences, plus a level-C
  *association* with independent HDX-MS data. Not that a method explains resistance, predicts drug response or has any
  clinical value.
- **Data [O]:** unbiased all-atom MD on Anton 2 (a99SB-disp protein, GAFF imatinib, 310 K). The data availability
  statement says the trajectories "are available for non-commercial use through contacting
  trajectories@deshawresearch.com". Relevant sets, as stated in the paper: from PDB 1OPJ (DFG-out, A-loop closed): WT,
  M472I and F486S, apo and imatinib-bound; from 3K5V: imatinib + GNF-2 for WT, M472I and F486S, 10 × 10 µs each; from
  2F4J: WT, M472I, F486S, 3 × 20 µs each.
- **Inconsistencies [O]:** run counts and lengths differ between Results and Methods for several sets (e.g. 1OPJ: 5 ×
  20 µs WT and 10 × 10 µs mutants in Methods vs 10 × 10 µs in Results). The simulated residue range is not stated (the
  NMR construct is residues 251–533). These must be resolved with the data providers before any protocol is frozen.
- **Size [A]:** unknown; frame interval not stated. At an assumed 1 ns interval and protein heavy atoms only
  (~2,300 atoms, float32), one 10 × 10 µs set is ~100,000 frames, ~2.8 GB. Full solvated frames would be far larger.
- **Atoms required:** Cα for TopoFold; backbone and side-chain heavy atoms for the PENSA baseline and for DFG
  side-chain state.
- **Replicates and dependence [O/D]:** within a set, runs share a starting structure (independent velocities assumed),
  so they are replicates of one initial condition, not independent samples of the equilibrium ensemble. 10-µs runs
  are long compared with local backbone relaxation but may be short compared with DFG or A-loop transitions; the
  inferential unit is the run, not the frame.
- **Confounds [D]:** different starting structures across sets; protonation choices (Asp381, imatinib N34 varied by
  stage); run length and count may differ between WT and mutants; the authors' own interpretation (C-lobe instability)
  is published, so the expected answer is known. Evaluation outcomes must be judged against the HDX-MS data, not
  against the paper's narrative.
- **External experimental evidence [O]:** HDX-MS (M472I vs WT, apo; also A344V), thermal denaturation (8–13 K Tm
  reduction for M343T, M351T, M472I, F486S), stopped-flow imatinib binding kinetics, NMR, cell IC50. These come from
  the same paper; they are experimental and independent of the simulations, but not of the authors.
- **Strongest baseline:** PENSA on backbone and side-chain torsions and distances with the same frames and residues;
  local Cα RMSD and Cα–Cα distances as equal-information baselines; per-residue RMSF difference as a simple reference.
- **Stated terms [O]:** the paper offers the trajectories "for non-commercial use" on request. Whether internal method
  evaluation by a project with future commercial ambitions, publication of derived results, or later commercial use is
  covered is for D. E. Shaw Research to answer (§3).

### GPCRmd μ-opioid receptor — apo vs BU72-bound (PENSA reference workflow)

- **Why this alternative [D]:** it is the dataset PENSA's own tutorial uses. That makes it the cleanest way to compare
  a TopoFold ensemble-comparison method with the strongest established baseline on identical, openly downloadable
  inputs. It has no oncology relevance; it would test method behaviour only.
- **Question [P]:** with replicate-level uncertainty, does a Cα-based comparison locate the same receptor regions that
  PENSA flags as differing between apo and BU72-bound simulations, and does it miss or add regions?
- **Permitted claim:** level B, method agreement or disagreement with PENSA on one system. No pharmacological or
  oncology claim; PENSA is a reference method here, not ground truth.
- **Data [O]:** GPCRmd dynamics 151 (apo; 70,430 atoms; 3 replicates; 1.1 µs accumulated) and 169 (BU72; 70,461
  atoms; 3 replicates; 1.25 µs accumulated). Both start from PDB 5C1M (active state), CHARMM36m, 4 fs timestep, frame
  interval 0.2 ns, unbiased NVT, run with ACEMD. Files: `.xtc` trajectories with `.psf` and `.pdb` topologies; per-file
  download from each report page (IDs listed in the PENSA tutorial).
- **Terms [O]:** the GPCRmd legal notice states the dataset "is generally open for personal, academic, and
  non-commercial use"; commercial use needs a separate licence; "Any download of large datasets must be requested in
  advance"; acknowledgement is required. No Creative Commons licence.
- **Frames [D]:** at 0.2 ns per saved frame, 1.1 µs accumulated is ~5,500 frames for the apo condition and 1.25 µs is
  ~6,250 frames for BU72, ~11,750 in total. Per replicate, ~1,800 (apo) and ~2,100 (BU72) if the three runs are of
  equal length, which the report does not state.
- **Size [A]:** not measured. Full-system frames of ~70,000 atoms in XTC are likely to total several GB; the GPCRmd
  advance-request rule for large downloads may apply.
- **Replicates and dependence [O/D]:** 3 vs 3 replicates from the same starting structure, each roughly 0.4 µs if of equal length
  (accumulated time ÷ 3). For an **exact two-sample permutation test that treats the six runs as exchangeable units**
  (one summary statistic per run, condition labels permuted among runs), there are C(6,3) = 20 assignments, so the
  smallest attainable p-value is 1/20 = 0.05 one-sided and 2/20 = 0.1 two-sided. This limit applies to that test only;
  other tests (parametric or hierarchical models across runs) can give smaller p-values only by adding distributional
  assumptions that 3 runs per condition cannot check. Frame-level tests that ignore run structure would overstate
  certainty.
- **Confounds [D]:** both conditions start from the agonist-bound active structure, so apo runs may
  be relaxing away from it; differences can reflect starting-state relaxation rather than ligand effect. Different
  simulation lengths across conditions.
- **External experimental evidence [O/D]:** crystal structures (active 5C1M; an inactive structure exists and would
  need verification before use) and receptor pharmacology. No residue-level experimental measurement of the apo-vs-BU72
  ensemble difference was identified, so no level-C test is possible.
- **Strongest baseline:** PENSA (its published workflow on these exact files).

## 2. Comparison

| | ABL mutants vs WT | GPCRmd μOR apo vs BU72 |
|---|---|---|
| Scientific value | Oncology-relevant; association with independent HDX-MS possible (B + C association) | Method calibration only (B); no oncology link |
| Access | Email request to DESRES (third-party message; owner action) | Direct per-file download; large downloads need advance request |
| Stated terms | "Non-commercial use" on request; scope for this project unanswered | Academic and non-commercial use; commercial licence on request |
| Size | Unknown; several GB per set if protein-only (assumption) | ~11,750 frames in total; several GB (assumption) |
| Replicates | 3–10 runs per condition, shared starts; counts inconsistent in the paper | 3 vs 3; exact run-level permutation test: min p 0.05 (one-sided) / 0.1 (two-sided) |
| Main confounds | Different starts and protonation across sets; published expected answer | Start-state relaxation; unequal lengths |
| Strongest baseline | PENSA plus local RMSD and Cα distances | PENSA (reference workflow) |
| Smallest informative evaluation | One pair (WT vs M472I, one starting structure): per-residue difference localization with run-level uncertainty, compared with HDX-MS regions and PENSA | Reproduce the PENSA tutorial result, then a Cα-based comparison on the same frames; region-level agreement with run-level uncertainty |
| Prerequisites | Owner authorization to request; DESRES permission covering the intended uses; resolved run inventory; XTC/DCD adapter via MDAnalysis; frozen protocol | Owner authorization; GPCRmd terms covering the intended uses (commercial licence if needed); isolated PENSA environment; frozen protocol |
| Estimated cost [A] | 5–10 days after data arrive; < 20 CPU-h; 5–20 GB | 3–6 days; < 5 CPU-h; < 5 GB |
| Ready now? | No: access and provider permission | No: provider terms for this project |

## 3. Permissions, by provider

Owner authorization to pursue a dataset is a project decision. It is **not** permission from the data provider, and
neither substitutes for the other.

| Provider | Dataset | What the source states [O] | Unanswered for this project | Who answers |
|---|---|---|---|---|
| Rascore maintainers | C1 SW1 labels | Code Apache-2.0; article CC BY-NC-ND 4.0 with unstated coverage of its datasets; no data-specific statement in the repository | All label-data terms. They are **not** established as non-commercial-only; they are unresolved (`C1_FEASIBILITY_CHECKS_2026-10-13.md` §1) | Rascore maintainers (request drafted, `C1_RASCORE_TERMS_REQUEST_DRAFT.md`) |
| D. E. Shaw Research | ABL trajectories (Ayaz et al. 2023) | Available "for non-commercial use" on request | Internal method evaluation by a prototype with commercial ambitions; publication of derived results; redistribution; later commercial use | DESRES (request drafted, `ABL_DESRES_DATA_REQUEST_DRAFT.md`) |
| GPCRmd | μOR dynamics 151, 169 | Open for "personal, academic, and non-commercial use"; commercial licence available on request; large downloads requested in advance | Whether this project's evaluation counts as academic non-commercial use | GPCRmd (no request drafted) |

## 4. Recommendation

**Neither candidate is ready to execute.** Each needs its own provider's answer (§3); ABL also needs the data request
to be sent and answered.

**Recommended candidate: ABL (WT vs M472I).** It is the only candidate that connects an ensemble difference to
independent experimental measurements in an oncology target. GPCRmd μOR is retained as a cheaper method-calibration
alternative, useful mainly for checking agreement with PENSA, not for the mission question.

**Next gate:** the owner reviews and, if appropriate, sends `ABL_DESRES_DATA_REQUEST_DRAFT.md`. ABL proceeds only if
DESRES's reply covers the intended uses and the run inventory supports a matched apo WT vs M472I comparison from one
starting structure. Otherwise ABL stays conditional and nothing is substituted without its own gates.

**What a positive ABL result would and would not show.** It could show that a method locates replicated structural
ensemble differences that coincide with HDX-MS differences. It would not show that those differences cause
resistance, predict drug response or have therapeutic value.

## 5. Evaluation checklist (ABL, before any protocol is frozen)

1. **Experimental reference:** peptide-level HDX-MS uptake for apo WT and M472I, with replicate uncertainty and the
   peptide map, obtained before any TopoFold result on these trajectories is seen. Regions are compared at peptide
   resolution, not residue resolution.
2. **No tuning to the published answer:** representation, window, distance, aggregation and thresholds are fixed in
   the frozen protocol before scoring. The published C-lobe interpretation and HDX figures are not used to choose
   them; any choice informed by them is disclosed and the result reported as exploratory.
3. **Replication:** the run is the statistical unit; WT and M472I are compared over matched analysis intervals
   (same starting structure, same discarded equilibration, same length per run, same frame interval).
4. **Equal-information baselines:** Cα-only local RMSD, Cα–Cα distances and per-residue Cα RMSF difference, on the
   same residues, runs and frames as TopoFold.
5. **Practical comparator, reported separately:** full-atom PENSA (backbone and side-chain torsions), labelled as
   having more information than the Cα-only methods.
