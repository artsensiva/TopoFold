# ABL WT vs M472I evaluation: preregistration DRAFT (2026-10-24)

**Status: DRAFT. Not frozen, not executed. No ABL trajectory or uptake value has been obtained, inspected or scored.**
Revised 2026-10-24 after a scientific audit; see §12. **Corrected 2026-10-24 (D55): the premise "no apo M472I
simulation is described" was wrong; see §0 and §13.**
It is frozen only after the inputs arrive and pass `benchmarks/audit/abl_prereg/validate_inputs.py`. Every value marked
*[fix at freeze]* is chosen then, before any comparison is computed.

**Source:** Ayaz et al., "Structural mechanism of a drug-binding process involving a large conformational change of the
protein target", *Nat Commun* 14:1885 (2023), doi:10.1038/s41467-023-36956-5. Everything below is from its text, not its
figures.

## 0. Simulated systems described in the paper (corrected 2026-10-24, D55)

**Correction.** An earlier version of this section stated: "**No apo M472I simulation is described.** … An apo-vs-apo
MD/HDX comparison is therefore possible only if unpublished apo M472I runs exist." That was wrong. It read only the
Methods subsection "Simulations starting from the DFG-out/AL-closed conformation" and missed the next subsection,
"Simulations starting from the DFG-in/AL-open conformation". Re-read 2026-10-24 from PMC10076256 (HTML SHA-256
`d1fcf3c7…`, fetched that day). The superseded statement is kept here, and in D47, as the record of the error.

The Methods describe, with their words in quotation marks:

| Methods subsection and start | System | Runs | Ligand as written in that subsection |
|---|---|---|---|
| DFG-out/AL-closed, 1OPJ | WT + imatinib | 5 × 20 µs, plus 5 × 10 µs | imatinib |
| same | M472I + imatinib; F486S + imatinib | 10 × 10 µs each | imatinib |
| same | "apo WT in the DFG-out/AL-closed conformation" | 5 × 20 µs | none ("apo") |
| same, 3K5V pose | WT, M472I, F486S "each simultaneously bound with imatinib and GNF-2" | 10 × 10 µs each | imatinib + GNF-2 |
| DFG-in/AL-open, 2F4J | "(1) WT, (2) M472I, and (3) F486S Abl kinase" | "Three independent 20-µs simulations … for each system" | **none specified** |

The 2F4J systems are "~80 Å per side, containing ~65,000 atoms total" and were run "following the same methods used for
the G0 simulations". The G0 Methods state: a99SB-disp force field; Anton 2; NVT at 310 K with a Nosé–Hoover
thermostat; 2.5-fs step. For G0 they also state that the DFG Asp was protonated, and that three imatinib molecules were
placed in the solvent. **Whether those two G0 preparation steps apply to the 2F4J-start systems is not stated.**

**Established, and not established:**
- **Stated in the paper:**
  - WT and M472I systems prepared from 2F4J, three runs of 20 µs each;
  - no ligand is named in that subsection;
  - the main text does not discuss these runs.
- **Not stated, and not verified by us:**
  - that the systems are ligand-free (no ligand specified ≠ a verified ligand-free supplied system);
  - VX-680 removal;
  - reversal of 2F4J's engineered H396P mutation to the WT His396;
  - treatment of the M228 cloning artefact and the construct range;
  - the DFG Asp381 protonation;
  - how the three runs were initialised;
  - the frame interval and file sizes;
  - whether the WT and M472I systems share every preparation step except residue 472.
- **The start structure itself:** PDB 2F4J is "the kinase domain of an imatinib-resistant Abl mutant in complex with the
  Aurora kinase inhibitor VX-680" (Young et al. 2006). In the local copy (`benchmarks/audit/kinase_descriptors/
  data_pdb/2F4J.pdb`, SHA-256 `50fd239b…`), chain A is residues 227–513, `SEQADV` lists PRO 396 as an engineered
  mutation of HIS and MET 228 as a cloning artefact, and VX-680 (VX6) is bound. Its conformation is DFG-in/AL-open
  (active), so the start is not the solution ensemble that apo HDX measures.
- **No Ayaz trajectory or derived value is held, inspected or scored.** The DESRES request (sent) and the HDX request
  (sent) are pending. The sent DESRES text contains the superseded premise; an optional follow-up naming the 2F4J sets
  is drafted, not sent (`docs/ABL_DESRES_FOLLOWUP_DRAFT.md`).

## 1. Question

Does a per-peptide summary of TopoFold's local backbone-geometry differences (κ/τ distributions) between M472I and WT
ensembles rank the peptides with increased apo HDX (M472I − WT) above the other peptides? Is that ranking better than
equal-information baselines?

The question is about **localization agreement** with an independent experimental endpoint. It is not about mechanism,
causality, drug resistance or clinical relevance.

**Branches**, chosen at freeze from the inventory and the verification gates (§2), before any trajectory frame is
analysed:
- **A (preferred, corrected 2026-10-24):** WT vs M472I started from 2F4J (3 × 20 µs each, per the Methods), matched to
  apo HDX. **Only** if the supplied topology and preparation records verify gates V1–V5 (§2): ligand-free, sequence and
  mutation identity, numbering, and shared starting conditions.
  - Limits that stay with the result: three runs per variant; one active-state start (DFG-in/AL-open), so the ensembles
    depend on that start over 20 µs; the HDX apo state is a solution ensemble of unknown composition.
  - The 1OPJ-start apo WT runs have no M472I counterpart from that start; they are **not** paired with 2F4J-start M472I.
- **B (exploratory cross-state evidence only):** imatinib-bound WT vs imatinib-bound M472I (1OPJ start), compared with
  apo HDX.
  - The MD ligand state differs from the HDX state, so a branch-B result is **exploratory cross-state evidence**. It is
    never reported as an evaluation outcome of the primary question.
  - It uses the separate exploratory outcome labels in §9, not the primary categories.
  - The paper makes this cross-state comparison qualitatively. Branch B cannot show agreement between MD and HDX in the
    same state.
- **C (sensitivity only):** the imatinib + GNF-2 sets (3K5V start).

## 2. Required trajectory inventory

Schema: `benchmarks/audit/abl_prereg/schemas/trajectory_inventory.schema.json`.

**Verification gates for branch A (from the supplied topology and preparation records, before any frame is
analysed):**
- **V1 ligand-free:** the topology of every WT and M472I system contains no residue other than protein, water and ions
  (in particular no VX-680 and no imatinib). The provider's word alone is recorded but is not sufficient.
- **V2 sequence and mutation:** the protein sequence of each system is aligned to P00519 with
  `topofold_inputs.kinase_roles` under an explicitly declared numbering. WT must carry His396 (2F4J's H396P reversed)
  and Met472; M472I must differ from WT only at 472 (Ile). Any other difference from P00519, including residue 228,
  must be identical in both systems and is recorded.
- **V3 numbering:** the construct range and numbering are confirmed by that alignment, not by author numbers.
- **V4 shared starting conditions:** same 2F4J-derived start coordinates, force field, temperature, ensemble, ion
  concentration, box size within the stated range, and Asp381 protonation, documented for both systems.
- **V5 replicate identity:** three runs per variant, each with a run identifier, length and initialisation record
  (e.g. new velocities from the same start).

**Minimum to proceed:**
- written permission for internal evaluation **and** for publication of derived results;
- for the chosen branch, at least 3 runs per variant from the identical start, with the same force field, temperature
  and protonation;
- protein heavy atoms, which the protection-factor baseline needs; Cα only is not enough;
- the construct range and numbering, with an offset to c-Abl isoform 1a if a different numbering is used;
- file SHA-256 values.

**Recorded, not required:**
- frame interval and per-run length;
- whether runs share a start: shared-start runs are conditionally independent, not independent starting states;
- any discrepancy with the paper's counts (Methods vs Results).

## 3. Residue-number mapping

- **Reference:** c-Abl isoform 1a, UniProt P00519 canonical (`reference/UniProt_P00519.fasta`, SHA-256 `f28f5df8…`).
  Landmarks verified: K271, E286, T315, A344, D381-F382-G383, M472.
- **HDX peptides:** `benchmarks/audit/abl_prereg/hdx_peptide_mapping.json`, made by `map_peptides.py` from sequences
  only.
  - c-Abl = construct + 228 (76 of 80 peptides match exactly and uniquely).
  - The 4 non-matching peptides are the A344V and M472I mutant peptides at exactly those positions.
  - HDX coverage: c-Abl 238–512, 252 of 275 positions.
- **MD:** the construct numbering is taken from the inventory and confirmed by exact sequence matching of the topology to
  P00519. A mismatch stops the analysis (§10).

## 4. Exclusions (fixed in advance)

- **Peptides spanning residue 472:** their sequences differ between WT and M472I, so their uptake is not directly
  comparable. They are excluded from the primary analysis and reported separately.
- **Peptides outside the MD construct:** a peptide with any residue outside the MD construct, or without a complete κ
  window inside it, is excluded. The causes are recorded with the existing exclusion machinery.
- **Insufficient replication:** peptides without at least two replicate values (or a provider-given spread) per state and
  time are excluded.
- **Equilibration:** the first 10 % of each run is excluded *[fix at freeze; no tuning]*.
- **Imatinib branch:** no ligand atoms are used; the protein only.

## 5. Independent endpoint (HDX)

Schema for one row of the processed table: `benchmarks/audit/abl_prereg/schemas/hdx_uptake.schema.json`. Table-level
rules are enforced by `validate_inputs.py`:
- unique keys;
- at least two replicates;
- identical peptide and time coverage across the compared states;
- sequences consistent with P00519 at offset 228, allowing only the stated mutation.

**Classification:**
- **Primary:** a peptide is "increased" if the HDX authors' own significance call, if supplied, marks M472I > WT at any
  labelling time.
- **Otherwise:** the hybrid significance test of Hageman & Weis, *Anal Chem* 91:8008–8016 (2019),
  doi:10.1021/acs.analchem.9b01325, with its published settings and α = 0.05, computed from the replicate values. This is
  an established test belonging to the endpoint, not a TopoFold statistic.
- **Fallback, only if processed values never arrive:** a coarse text-level endpoint. The regions the paper names as
  increased in M472I are Leu452–Glu470, Tyr353–Leu370 and Ile402–Ser410. It would be used as a region-level label with
  lower status, and reported as such.
- **Not used:** values read from figures (Fig. 4E, Supplementary Fig. 4B).

## 6. Methods compared

**6a. Equal-information Cα baselines (the primary comparison).** These methods receive exactly what TopoFold receives:
the same Cα frames, residues and peptide list, and **no HDX information** (AGENTS rule 6).
- **TopoFold:** per-window pooled Jensen–Shannon distance between M472I and WT for κ and τ (`topofold_inputs.ensemble`,
  existing binning).
- **Window Cα distances:** the same windows, the same JS.
- **Local Cα RMSD:** the distribution difference of window RMSD to a common reference *[reference fixed at freeze: the
  start structure]*.
- **Per-residue Cα RMSF difference**, between the M472I and WT ensembles.
- **φ/ψ distributions:** these use backbone N and C as well as Cα. They are reported in a separately labelled
  backbone row, only if backbone extraction (capability "B") is available at freeze; otherwise listed as not run.

**6b. Stronger-information comparator (reported separately, not part of the primary comparison).**
- **Backbone-amide protection-factor model:** Best & Vendruscolo, *Structure* 14:97–106 (2006),
  doi:10.1016/j.str.2005.09.012, with ln P = βc·Nc + βh·Nh.
  - The parameters are taken from the primary source at freeze, with no fitting to ABL data.
  - It uses heavy atoms and amide hydrogens (built if absent, by a method fixed at freeze), which is more information
    than TopoFold's Cα input. It is therefore not an equal-information baseline.
  - It is the most directly HDX-related physical model. It is reported as a reference point, and TopoFold's result
    against it is stated separately (§9).

**Peptide score:** each method's per-residue or per-window values, in 6a and 6b alike, are aggregated identically. The mean is taken over the
windows fully inside the peptide *[fix at freeze]*.

## 7. Metrics

- **Primary:** AUROC of the peptide score for "increased" versus other peptides.
- **Secondary:**
  - Spearman correlation between the peptide score and the maximum ΔD over time;
  - precision among the top-k peptides, with k = the number of increased peptides.
- **Overlapping peptides:** their scores are not independent. Counts of peptides and of distinct covered residues are
  both reported.

## 8. Uncertainty

Calibration-based inference is closed under D33, so there are **no p-values or confidence intervals** for the
TopoFold-vs-baseline comparison.

**Reported instead:**
- Leave-one-run-out recomputation per variant: the min–max of each metric over the leave-one-out sets.
- Temporal diagnostics per run, from the existing tool.

A peptide-level bootstrap is not used: overlapping peptides and spatial correlation violate its independence
assumption. An inferential analysis could be added only as a separately frozen protocol, with a method calibrated in
advance on SYNTHETIC data matched to the delivered inventory (D33 conditions).

## 9. Outcomes

Thresholds are provisional, recorded here before any data are seen, and confirmed or changed only at freeze, never after
looking. "Best Cα baseline" means the highest AUROC among the equal-information baselines of §6a.

**Primary outcome (branch A).** The steps are applied in order and the first one that matches is the outcome, so the
categories are mutually exclusive and cover every case.

| Step | Outcome | Rule |
|---|---|---|
| 1 | **Inconclusive** | Fewer than 5 "increased" or fewer than 5 other peptides remain after exclusions, or the leave-one-run-out range of TopoFold's AUROC includes 0.5 |
| 2 | **Supportive (descriptive)** | All three hold: TopoFold AUROC ≥ 0.70; it exceeds the best Cα baseline by ≥ 0.05; its leave-one-run-out minimum is ≥ 0.60 |
| 3 | **Not supportive** | TopoFold AUROC < 0.60, or it is more than 0.05 below the best Cα baseline |
| 4 | **Comparable** | Every remaining case. For example, TopoFold within 0.05 of the best Cα baseline, or above it with AUROC between 0.60 and 0.70 |

**Separate statement against the stronger-information comparator (§6b).** This is reported for every primary outcome and
does not change it. It takes exactly one of three values: TopoFold AUROC higher than the protection-factor model by
≥ 0.05; within 0.05 of it; or lower by more than 0.05.

**Branch B (exploratory cross-state evidence).** It uses the same ordered steps but different labels:
- **inconclusive (exploratory)**;
- **exploratory agreement**, in place of "supportive";
- **exploratory disagreement**, in place of "not supportive";
- **exploratory, comparable to baselines**.

These are never reported as an evaluation of the primary question, and they cannot be converted into the primary
categories later.

**Reporting:** every outcome is reported, including "not supportive" and "exploratory disagreement". In any branch the
result is development evidence of localization agreement, not biological validation of a mechanism.

## 10. Stop conditions (no analysis, or analysis halted)

1. Permission for internal evaluation or publication is refused or unknown.
2. Fewer than 3 runs per variant from a matched start for every branch.
3. The inventory contradicts the paper's description and the provider cannot resolve it.
4. Sequence mapping of the MD topology to P00519 fails, or the HDX table fails validation.
5. HDX lacks apo M472I or apo WT, or lacks replicate information.
5a. **Branch A only (added 2026-10-24, D55):** any of gates V1–V5 (§2) fails or cannot be checked from supplied files.
    In particular, "no ligand specified" in the Methods is not accepted as ligand-free. Branch A then stops; it is not
    silently replaced by another branch.
5b. The supplied 2F4J-start runs are fewer than three per variant, or differ in length or start between WT and M472I
    without a documented reason.
6. Branch B is the only option and the owner has not approved running it as exploratory cross-state evidence.
7. Any change to this protocol after data receipt that is not recorded as a dated amendment made before scoring.

## 11. Freeze procedure

1. Validate the deliveries:
   `validate_inputs.py --inventory … --hdx …`.
2. Record the branch decision and fill every *[fix at freeze]* item.
3. Write the frozen spec and code hashes, in the same pattern as the KRAS and DABOOT freezes.
4. Only then extract coordinates and compute.

## 12. Audit record (2026-10-24, `/scientific-audit`, draft not frozen)

| Section | Result | Evidence or action |
|---|---|---|
| 1 Status | pass | Labelled "DRAFT. Not frozen, not executed" |
| 2 Freeze | not applicable | No freeze yet. The freeze procedure (§11) follows the hash pattern |
| 3 Evidence class | **fixed** | Branch B relabelled as exploratory cross-state evidence with its own outcome labels |
| 4 Preregistration | **fixed** | The outcome categories overlapped (for example, a case could be both "comparable" and "not supportive"); they are now an ordered, mutually exclusive procedure. The protection-factor model was listed as an "equal-information baseline" although it uses heavy atoms and hydrogens; it is now a separate stronger-information comparator (§6b). Thresholds remain provisional *[fix at freeze]* |
| 5 Provenance | partial | Inputs that exist are hashed: peptide map `2fd6b624…`, UniProt FASTA `f28f5df8…`. MD and HDX inputs are pending delivery. Output manifests are to be added in the freeze |
| 6 Corruption tests | open | To be predeclared at freeze, for example a peptide-label permutation that preserves sequence contiguity, and a residue-offset shift of the mapping. Each should give AUROC near 0.5 |
| 7 Frozen artifacts | not applicable | Nothing frozen |
| 8 Negative outcomes | pass | "Not supportive", "inconclusive" and the exploratory negatives are reported |
| 9 Claim boundaries | pass | Development evidence only; no mechanism, causation or clinical claim |

## 13. Correction record (2026-10-24, D55)

| Item | Before | After |
|---|---|---|
| §0 premise | "No apo M472I simulation is described"; apo-vs-apo possible only with unpublished runs | The Methods describe WT, M472I and F486S systems from 2F4J, 3 × 20 µs each, with no ligand specified; ligand-free status is unverified |
| Branch A | apo WT vs apo M472I, existence unknown | WT vs M472I from 2F4J, conditional on gates V1–V5 |
| Stop conditions | 1–7 | adds 5a (gates) and 5b (replicate match) |
| Unchanged | | §§1, 3–9, 11 and 12; branches B and C; the outcome rules and thresholds. The draft remains unfrozen and nothing was scored |

Cause of the error: one Methods subsection was read and the adjacent one was not. Lesson recorded: inventory claims
of the form "X is not described" are checked against every Methods subsection and the data-availability statement
before they are used.

