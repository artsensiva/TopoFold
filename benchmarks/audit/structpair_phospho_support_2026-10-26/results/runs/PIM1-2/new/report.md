# Two-structure comparison: KIN-PAIR-1 PIM1-2: 6VRV vs 6L17 (PIM1) (derived regression spec: position 261 added, modified_residues phospho_pdb)

**Descriptive, static.** Two static structures: the values are direct geometric differences between these coordinate sets. They can reflect crystal form and packing, crystallization conditions, resolution, refinement, alternate conformations and coordinate error as well as any biological difference. No ensemble, population, kinetic, functional-state or treatment-response claim follows. Agreement or disagreement between metrics does not show that one representation is better. Crystal copies of one entry are not independent replicates.

## Inputs

| | Label | File (SHA-256, first 12) | Chain | Altloc policy |
|---|---|---|---|---|
| A | 6VRV chain A | `6VRV.pdb` (5fb4eae4ba16) | A | blank_then_A |
| B | 6L17 chain A | `6L17.pdb` (3ad04425727d) | A | blank_then_A |

## Within-entry crystal-copy variation

**Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of them. The within-entry crystal-copy variation is a descriptive range of differences among the chains the user selected. It is not a calibrated noise floor, not a statistical uncertainty, not a nucleotide or ligand effect and not an HDX replicate estimate. An advisory flag appears only for a channel whose minimum difference the user configured; it is an advisory rule, not a gate, not an adequacy verdict, not significance and not biological separation.**

**Unavailable: one chain per entry.** Name every chain of an entry in `structures.X.chains` (with `primary_chain`) to obtain it.


## Mapping

- Correspondence: pairs (271 positions; 271 present in both). Missing-residue policy: `exclude`.
- Residue-name mismatches (flag): A:261 SER→A:261 SEP.
- Alternate locations resolved: none.
- Breaks (adapter): A none, B none (row index after which the chain breaks).

## Global Cα RMSD (context only)

- 0.450 Å over 271 Cα after optimal superposition (MDAnalysis QCP); independent numpy SVD/Kabsch check 0.450 Å (|difference| 8.4e-15 Å).
- Published RMSD values often use other atom sets, residue ranges or outlier rejection; compare only with a stated definition.

## Largest local differences (top 10; descriptive)

| Metric | Unit | Windows (A numbering; value; * altloc, † residue-name mismatch in window) |
|---|---|---|
| d_kappa | rad | A:100..A:102 (0.748), A:98..A:100 (0.622), A:97..A:99 (0.570), A:96..A:98 (0.433), A:75..A:77 (0.333), A:99..A:101 (0.277), A:81..A:83 (0.251), A:101..A:103 (0.156), A:95..A:97 (0.117), A:50..A:52 (0.112) |
| d_i_i2 | Å | A:98..A:100 (1.433), A:100..A:102 (1.336), A:96..A:98 (0.885), A:97..A:99 (0.828), A:99..A:101 (0.633), A:81..A:83 (0.600), A:75..A:77 (0.577), A:101..A:103 (0.465), A:95..A:97 (0.293), A:73..A:75 (0.226) |
| rmsd3 | Å | A:100..A:102 (0.804), A:98..A:100 (0.733), A:97..A:99 (0.581), A:96..A:98 (0.492), A:75..A:77 (0.356), A:99..A:101 (0.325), A:81..A:83 (0.303), A:101..A:103 (0.218), A:95..A:97 (0.146), A:50..A:52 (0.112) |
| d_tau | rad | A:100..A:103 (3.041), A:98..A:101 (2.561), A:99..A:102 (0.984), A:302..A:305 (0.850), A:101..A:104 (0.828), A:97..A:100 (0.789), A:95..A:98 (0.762), A:74..A:77 (0.667), A:73..A:76 (0.566), A:75..A:78 (0.432) |
| d_i_i3 | Å | A:95..A:98 (1.562), A:96..A:99 (1.500), A:99..A:102 (1.125), A:302..A:305 (1.125), A:100..A:103 (0.901), A:101..A:104 (0.746), A:72..A:75 (0.485), A:80..A:83 (0.447), A:98..A:101 (0.429), A:136..A:139 (0.369) |
| rmsd4 | Å | A:100..A:103 (1.561), A:98..A:101 (1.410), A:99..A:102 (1.076), A:95..A:98 (0.902), A:97..A:100 (0.751), A:96..A:99 (0.723), A:302..A:305 (0.674), A:101..A:104 (0.610), A:73..A:76 (0.524), A:75..A:78 (0.472) |

## Agreement between equal-information metrics (Spearman over windows; descriptive)

| Pair | Spearman | Windows |
|---|---|---|
| d_kappa~d_i_i2 | 0.435 | 269 |
| d_kappa~rmsd3 | 0.748 | 269 |
| d_i_i2~rmsd3 | 0.716 | 269 |
| d_tau~d_i_i3 | 0.331 | 268 |
| d_tau~rmsd4 | 0.732 | 268 |
| d_i_i3~rmsd4 | 0.566 | 268 |

Three-residue metrics share one triangle, so close agreement is expected; d(i−1, i+2) is blind to the sign of τ.

## Exclusions

No window inside the mapped range was excluded.

## Modified residues

Policy `phospho_pdb`: SEP, TPO and PTR in PDB-format files with MODRES records; no other modified residue or format. The deposited residue name is kept; the parent is documentation for the correspondence and is not a rename (SEP versus SER stays a visible name difference). A modification annotation is not evidence of functional state or of any effect.

| Structure | Residue | Deposited name | Parent (source) | Annotation | Record |
|---|---|---|---|---|---|
| B | A:261 | SEP | SER (MODRES record for A:261 (parent SER) agrees with the policy table SEP->SER (wwPDB chemical components)) | MODIFIED RESIDUE | HETATM |

## Alternate-conformation sensitivity

not applicable: no alternate locations at the mapped Cα atoms

## Figure

![structure pair](figure.png)

Files: `windows.csv` (all metrics with units in column names), `residues.csv`, `exclusions.csv`, `summary.json`, `manifest.json` (hashes and software), `numeric_record.json` (computed values before display rounding). structpair 0.5.0.
