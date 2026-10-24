# ABL WT vs M472I evaluation: preregistration DRAFT (2026-10-24)

**Status: DRAFT. Not frozen, not executed. No ABL trajectory or uptake value has been obtained, inspected or scored.**
Revised 2026-10-24 after a scientific audit; see §12.
It is frozen only after the inputs arrive and pass `benchmarks/audit/abl_prereg/validate_inputs.py`. Every value marked
*[fix at freeze]* is chosen then, before any comparison is computed.

**Source:** Ayaz et al., "Structural mechanism of a drug-binding process involving a large conformational change of the
protein target", *Nat Commun* 14:1885 (2023), doi:10.1038/s41467-023-36956-5. Everything below is from its text, not its
figures.

## 0. Finding that shapes this design (public-data observation)

The Methods list these Anton 2 systems with WT/M472I relevance: "(1) WT bound with imatinib; (2) and (3) M472I and F486S
mutants bound with imatinib; and (4) apo WT in the DFG-out/AL-closed conformation", all from PDB 1OPJ.

| System (start 1OPJ) | Runs |
|---|---|
| WT + imatinib | 5 × 20 µs, plus 5 × 10 µs |
| M472I + imatinib | 10 × 10 µs |
| apo WT | 5 × 20 µs |
| WT, M472I and F486S, each with imatinib and GNF-2 (start 3K5V) | 10 × 10 µs each |

**No apo M472I simulation is described.** The HDX-MS comparison, by contrast, is apo: "M472I increases HDX in the C-lobe
of the kinase in the apo form".
- An apo-vs-apo MD/HDX comparison is therefore possible only if unpublished apo M472I runs exist.
- The matched sets that are described are imatinib-bound.
- The owner has sent the trajectory request to D. E. Shaw Research, which asks about apo M472I directly, and the
  processed-HDX request. Drafts: `docs/ABL_DESRES_DATA_REQUEST_DRAFT.md` and `docs/ABL_HDX_DATA_REQUEST_DRAFT.md`.
- **Whether apo M472I trajectories exist remains unknown until D. E. Shaw Research replies.**

## 1. Question

Does a per-peptide summary of TopoFold's local backbone-geometry differences (κ/τ distributions) between M472I and WT
ensembles rank the peptides with increased apo HDX (M472I − WT) above the other peptides? Is that ranking better than
equal-information baselines?

The question is about **localization agreement** with an independent experimental endpoint. It is not about mechanism,
causality, drug resistance or clinical relevance.

**Branches**, chosen at freeze from the inventory alone, before any data are opened:
- **A (preferred):** apo WT vs apo M472I from the same start, matched to apo HDX.
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

