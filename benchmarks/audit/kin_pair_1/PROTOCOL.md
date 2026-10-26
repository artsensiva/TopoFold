# KIN-PAIR-1: structpair on previously unused human kinases (protocol, frozen 2026-10-26)

**Status: frozen before panel selection and before any structpair execution** (freeze record
`results/KIN-PAIR-1_freeze_1.json`, which hashes this file, `protocol.json`, the selection and helper code and the two input files
`results/candidate_targets_verified.json` and `results/exposure_history.json`). Machine-readable parameters: `protocol.json`.
A change after the freeze is a dated amendment at the end of this file, made before any re-scoring.

## 1. Question and permitted conclusions

**Question.** Can the existing `structpair` workflow (v0.3.0, commit `2af92d5`) correctly map, analyze and explain structural
differences between public experimental structures of human kinases not used during TopoFold development, with transparent
exclusions and agreement with independent numerical references?

**This is a structural-analysis and engineering evaluation.** It does not test drug efficacy, resistance, causal allostery,
function, kinetics or clinical benefit, and the structures are static.

**Permitted conclusions** (each limited to the panel actually executed):
1. For each pair: the residue correspondence is verified against the UniProt sequence and the deposited identifiers, or the pair
   is an *input or mapping failure*.
2. For each executed pair: whether structpair's window values, exclusion set, global RMSD, displacement and invariance behaviour
   agree with independent references within the stated tolerance, or the pair is an *implementation failure*.
3. Descriptive: whether the TopoFold local geometry (d_kappa, d_tau) is *redundant with* or *complementary to* the
   equal-information baselines (Cα-distance difference, local RMSD), by the rule in section 8, and how the displacement after
   alignment differs from local shape change.
4. Which exclusions occurred and their recorded causes.

**Not permitted:** any claim about biological state, ligand effect, mutation effect, resistance, allostery or mechanism; novelty
or biological superiority of any metric; user time savings; generality to kinases outside the panel; independence of pairs that
share a target, a ligand class, a laboratory or a crystal form. Chains of one crystal and several pairs of one target are not
independent replicates.

## 2. Exposure history (recorded, not hidden)

`results/exposure_history.json` lists every PDB-ID-like token in the repository's development records (docs, audit READMEs, specs,
results, tracked structure files) that RCSB confirms, with the UniProt accessions of their polymer entities (171 entries, 26
accessions, which include false-positive tokens that happen to be valid IDs). Targets with an exposed accession and entries with
an exposed ID are excluded. Kinases used in development include ABL1, CDK2, MAPK14, AURKA, BRAF and PRKACA. None of the 23
declared candidate accessions is among the exposed ones. The tool itself was developed and verified on KRAS and Abl entries, so
the *workflow* is not blind; the *targets and entries* are previously unused. The panel is therefore prospective only in that
sense. The selecting agent (an AI) knows these kinases from general knowledge and computed nothing on them before the freeze.

## 3. Candidate targets (declared before any structure query)

23 oncology-relevant human protein kinases, declared in `results/candidate_targets_verified.json` with gene, accession and the
UniProt verification (gene, organism 9606, name, length; 23 of 23 verified before the freeze): EGFR, ERBB2, MET, ALK, KIT, FLT3,
FGFR1, KDR, RET, SRC, LCK, JAK2, BTK, MAP2K1, CHEK1, PLK1, AKT1, CDK6, MAPK1, SYK, IRAK4, PIM1, PDPK1. The list is not chosen by
expected TopoFold performance. Lipid kinases and non-human kinases are not candidates.

## 4. Inclusion and exclusion criteria

*Entry level (metadata, RCSB search and data API):* X-ray diffraction; deposited resolution at most 2.5 Å; a protein polymer
entity whose reference sequence identifier is the target accession (RCSB/SIFTS); a PDB-format file exists at RCSB; the entry ID is not
in the exposure list.

*Entry level (input quality, from the downloaded file; fixed-column parser, no TopoFold score):* the target chain (chain rule,
below) maps to the UniProt canonical sequence by local alignment (the existing `kinase_roles` alignment) with identity at least
0.98 over aligned residues and at least 0.80 of its residues aligned; a DBREF record of the chain names the target accession and
every aligned residue inside a DBREF segment agrees with the DBREF numbering; every residue-name difference from UniProt has a
SEQADV record; at least 200 mapped residues are read by the adapter (Cα present, adapter protein residue name).

*Pair level:* different primary ligands; at least 200 shared mapped positions.

**Not eligible:** NMR, cryo-EM, entries over 2.5 Å, entries without a PDB-format file (mmCIF-only entries need a reader that
does not exist in structpair, and no adapter is built for this study), apo-only pairs (a primary ligand is required for both
entries, only to make a metadata-based contrast), and anything that fails a quality rule. Residues the adapter does not read
(for example phosphorylated TPO, SEP, PTR) are *not read*: they keep their place in the sequence mapping but are not in the
correspondence, and windows across them are excluded as discontinuities. This is a stated limitation of the workflow.

## 5. Deterministic enumeration and selection (`select_panel.py`)

1. For each non-exposed candidate target, one RCSB search (accession, X-ray, resolution at most 2.5 Å), sorted by resolution
   ascending then PDB ID; entries with an exposed ID are removed; the count is the metadata-eligible count.
2. Targets are ranked by that count (descending), ties by accession. A target needs at least 4 eligible entries.
3. Walking the ranking, a target is examined (at most 8 targets); the first 3 that yield 2 usable disjoint pairs form the panel.
   Maximum 3 targets and 6 pairs.
4. Pairs: entries ordered by (resolution, ID); E1 is the first usable entry not in an earlier pair; its partner is the first later
   usable entry with a primary ligand different from E1's that shares at least 200 positions. Primary ligand = the non-polymer
   component, outside a fixed additive list (`protocol.json`), with the largest RCSB formula weight (ties: component ID). Entry
   examination stops at 25 entries per target and 40 downloaded files in total.
5. **Chain rule.** Among the chains of the target entity, the one with the most adapter-read residues; ties alphabetical.
6. Every candidate target and every examined entry is recorded with its status and reason in
   `results/KIN-PAIR-1_selection.json`. If fewer than 3 targets qualify, the panel is whatever qualified; if none, that is reported
   once as *not executable* and the rules are not relaxed.

## 6. Mapping and structpair policies

- **Sequence correspondence:** each chain is aligned to the UniProt canonical sequence; the correspondence is the set of residues read
  in both chains that map to the same UniProt index, ordered by index, written as explicit `pairs` (author number to author number).
  Numbering is never assumed: author numbers are checked against DBREF, and offsets are reported.
- **Missing residues:** unobserved residues are absent from the correspondence; windows across them are excluded by structpair
  (`numbering_discontinuity` or `missing_CA`) and `on_missing` is `exclude`. The mapped-but-unread span positions are counted and listed.
- **Altlocs:** `blank_then_A` for both entries; `altloc_sensitivity: highest_occupancy`.
- **Residue names:** `resname_mismatch: flag`.
- **Second chain:** per `protocol.json`; chain-to-chain variation is descriptive only (not a noise floor, not replicates).
- **No advisory thresholds** (`min_difference` empty), so no separation flags; numbers only.
- **Recorded per pair, never used to adjust anything:** construct span (DBREF), engineered mutations (SEQADV and RCSB mutation count),
  primary and other ligands, space group, resolution, release dates, number of chains of the target entity.

## 7. Numerical references, tolerances and failure conditions

All references are independent of the workflow under test: fixed-column parsing of the same file text in float64 and numpy only
(`kinpair.py`), with the existing alignment helper used only for sequence mapping.

| Check | Reference | Tolerance | Failure |
|---|---|---|---|
| R1 mapping | UniProt sequence, DBREF, SEQADV (section 4) | exact agreement | *input/mapping failure* of the pair |
| R2 window set | independent window validity (numbering adjacency, Cα–Cα at most 4.2 Å, runs of at least 4 positions) | identical set per window size | any unexplained difference: *implementation failure* |
| R3 κ/τ-based values | d_kappa, d_i_i2, d_tau, d_i_i3 recomputed from raw Cα (turning angle, osculating-plane dihedral, distances) | 1e-4 (rad or Å) | any value beyond tolerance |
| R4 local RMSD | rmsd3, rmsd4 by Kabsch with proper rotation (SVD) | 1e-4 Å | any value beyond tolerance |
| R5 global RMSD | Kabsch over all positions read in both | 1e-4 Å | beyond tolerance |
| R6 displacement | trimmed Kabsch (20 % trim, ties lower position first) vs structpair `disp` where a multi-chain run provides it | 1e-4 Å | beyond tolerance; not applicable (stated) for single-chain pairs |
| R7 rigid-motion invariance | structpair re-run with A and B moved by exact proper rotations (coordinate permutations and sign changes, exact in the file) and exact translations | 1e-4 on all numeric window columns and the global RMSD | beyond tolerance |
| R8 corruption control | the mapping shifted by one UniProt position | identity below 0.5 | control cannot discriminate |
| R9 exclusions | structpair's excluded windows against R2 | identical | unexplained difference |

The tolerance of 1e-4 is declared because the adapter holds Cα positions as float32 (cause identified for the K1 comparison,
2026-10-25); the reference uses float64 from the file text. A crash of structpair on inputs that passed R1 is an implementation failure.
Equal-information baselines: for d_kappa the Cα distance across the same 3-residue window and rmsd3; for d_tau the Cα distance across
the same 4-residue window and rmsd4. Same chains, residues and windows for every metric.

## 8. Classification and study outcome

*Per pair and window size* (secondary, descriptive): let O be the number of the 10 largest windows of the geometry metric (ties:
lower position first) that are among the 10 largest windows of either baseline. O of 8 or more: **redundant** with the baselines;
O of 7 or fewer: **complementary descriptive information**; fewer than 50 valid windows: **insufficient evidence**. Spearman correlations
are reported. This is not novelty, not biological superiority and not a test.

*Study outcome* (first match applies): **NOT EXECUTABLE** (no pair selected); **IMPLEMENTATION FAILURE** (any executed pair fails R2 to R9,
or structpair crashes on a mapped pair); **SUPPORTED ON THIS PANEL** (at least 4 executed pairs from at least 2 targets, all passing
R2 to R9); **PARTIAL** (all executed pairs pass, fewer than 4 pairs or 2 targets). Mapping failures reduce the executed pairs and are
reported separately; they are not implementation failures.

## 9. Outputs and resources

Under `results/`: `KIN-PAIR-1_selection.json` (all candidates and reasons), `KIN-PAIR-1_results.json` (all numbers; text and figures read
from it), compact structpair outputs per pair, `report.md`, and at most two figures with captions in `figures/`. Downloads: at most 100 MB
(expected under 15 MB); execution at most 2 CPU-hours, installation excluded and not counted; resources are measured and reported.
Raw PDB files stay in the git-ignored `data_pdb/` with SHA-256 values in the results.

## 10. Stop rules

Stop and report if the selection cannot fit the download limit, if no pair is selected, or if a frozen hash differs. The protocol is not
changed after outcomes are seen, and the study is not repeated with different rules to obtain a favourable result.

## Amendments

- **Amendment 1 (2026-10-26), made after the first selection run and before any structpair execution or any score.** The first run of
  `select_panel.py run` (record kept as `results/KIN-PAIR-1_selection_run1_defective.json`, log `results/select_run1.log`) stopped after
  two targets (EGFR, PIM1; 4 pairs) because the code compared the number of *pairs* with `targets_needed` (3) instead of the number
  of selected *targets*. This contradicts section 5 (the first 3 targets that yield 2 usable disjoint pairs). The code was corrected
  (`n_selected` counts targets) and a regression test added (`test_selection_stops_after_three_targets_not_three_pairs`). Nothing else
  changed: the ranking, criteria, quality rules and tolerances are as frozen. The corrected run uses the same deterministic ranking, so
  EGFR and PIM1 are selected again by the same rules and the third target is the next one that qualifies. The selection saw only RCSB
  metadata, downloaded PDB files and sequence checks; no TopoFold score, window metric or plot had been computed. New freeze record:
  `results/KIN-PAIR-1_freeze_1a.json`.

