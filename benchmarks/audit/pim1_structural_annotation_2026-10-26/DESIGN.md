# Verified structural-role annotation of the existing PIM1 comparisons: rules fixed before the annotation is run (2026-10-26)

**Interpretation task on existing results; no score is recomputed, no threshold changed, no structure selected.** KIN-PAIR-1, the phospho regression and all earlier records are unchanged.
Question: where do the seven newly retained windows around residue 261 (per pair) lie in the verified PIM1 sequence and structural architecture, and what descriptive interpretation, if any, does that add?

## Inputs
- Structures: 8R0Q (PIM1-1 A), 6NO9 (PIM1-1 B), 6VRV (PIM1-2 A), 6L17 (PIM1-2 B); chain A of each; local files, hash-verified against `kin_pair_1/results/KIN-PAIR-1_selection.json`.
- Reference: UniProt **P11309, canonical isoform P11309-1 (PIM1_HUMAN), 313 residues**, the sequence recorded with hash in the selection record (`4de04d20…`); numbering convention **reference index = UniProt canonical position** (`Numbering.reference_index()`); author numbering is preserved separately and its offset to the reference is reported, never applied.
- Comparison outputs: `structpair_phospho_support_2026-10-26/results/runs/PIM1-{1,2}/{new,old_with_numeric_record}/` (windows identified by complete residue identifiers; no recomputation of metrics).
- Primary-source annotation (small retrieval, provenance recorded): the UniProt P11309 JSON entry (features: protein kinase domain 38–290 with PROSITE rule PRU00159, active site 167 proton acceptor, ATP binding site 67, modified residue 261 phosphoserine with its cited evidence).

## Anchors for the existing `kinase_roles` helper (declared here, not tuned afterwards)
| Anchor | Declared as | Basis |
|---|---|---|
| β3 Lys | reference index 67 | UniProt binding site 67 (ATP; cited literature in the entry); the helper checks that the reference has K there |
| HRD | HRD_H at reference index 165 | UniProt active site 167 (proton acceptor): the helper checks H, R, D at 165–167 |
| DFG | `MotifWindow(168, 290)` | from after the UniProt active site to the end of the UniProt kinase domain; the helper requires exactly one DFG in the window and refuses otherwise |
| αC Glu | the unique Glu of reference window 77–97 (β3 Lys + 10 to + 30) whose carboxylate O (OE1 or OE2) is within 4.0 Å of the Lys 67 NZ in **every** one of the four structures | UniProt has no αC annotation; the rule is a coordinate contact rule fixed now. If not exactly one residue qualifies, αC-Glu stays unverified and the helper (which needs all four anchors) is not run |
The helper is run with its defaults (`min_identity` 0.90, `min_aligned` 30, role mutations refused). If it refuses, the exact reason is reported and nothing is changed to make it pass; a retry with the documented `allow_role_mutations` flag is only allowed if the refusal is solely a role/motif mutation, and both outcomes are then reported.
Successful anchor mapping verifies the four anchors, the ambiguity checks and the mapped roles; it does **not** verify any helix or loop boundary.

## Region labels
Only: the UniProt protein kinase domain (38–290), the four anchors, and **deposited** secondary-structure records (HELIX and SHEET lines of each local file, labelled as deposited annotation, not independently verified). Conventional kinase element names (αC, activation loop, αH, αI, N-lobe, C-lobe boundaries and so on)
are used only where one of those sources supports them: a helix is called "the deposited helix containing the αC-Glu anchor" and nothing more. Anything else stays "unverified". No label is inferred from residue-number ranges.

## Annotation of the existing outputs
The seven new windows per pair (3 three-residue and 4 four-residue, identified by complete residue identifiers as present in the new numeric record and absent from the old) are reported with: residue identities (deposited number, reference index, deposited name, modification),
the verified location or the remaining uncertainty, their existing values and their rank among all windows of the same metric in the new run (descending; ties by position), and whether they enter the previously reported top-10 lists.
The existing KIN-PAIR-1 descriptive classification (redundant, complementary, insufficient: the recorded rule) is re-applied to the existing numeric values only to see whether the previously reported descriptive conclusion changes; the maximum trimmed-fit displacement and its location are not recomputed.
Modification annotation (SEP, MODRES, UniProt "modified residue") is kept separate from any functional interpretation; these pairs do not isolate phosphorylation, ligand or mutation effects (the two entries of a pair differ in construct, ligand and crystal form).

## Output
`results/pim1_annotation.json`, a generated `REPORT.md` with one captioned schematic of verified landmarks and the added windows (uncertain boundaries marked), and updates to the roadmap and decision log. No new structure or trajectory download.
