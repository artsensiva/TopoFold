# Two-structure comparison: KIN-PAIR-1 PIM1-1: 8R0Q vs 6NO9 (PIM1) (derived regression spec: position 261 added, modified_residues phospho_pdb)

**Descriptive, static.** Two static structures: the values are direct geometric differences between these coordinate sets. They can reflect crystal form and packing, crystallization conditions, resolution, refinement, alternate conformations and coordinate error as well as any biological difference. No ensemble, population, kinetic, functional-state or treatment-response claim follows. Agreement or disagreement between metrics does not show that one representation is better. Crystal copies of one entry are not independent replicates.

## Inputs

| | Label | File (SHA-256, first 12) | Chain | Altloc policy |
|---|---|---|---|---|
| A | 8R0Q chain A | `8R0Q.pdb` (6c9bd62287f5) | A | blank_then_A |
| B | 6NO9 chain A | `6NO9.pdb` (0115a741b19f) | A | blank_then_A |

## Within-entry crystal-copy variation

**Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of them. The within-entry crystal-copy variation is a descriptive range of differences among the chains the user selected. It is not a calibrated noise floor, not a statistical uncertainty, not a nucleotide or ligand effect and not an HDX replicate estimate. An advisory flag appears only for a channel whose minimum difference the user configured; it is an advisory rule, not a gate, not an adequacy verdict, not significance and not biological separation.**

**Unavailable: one chain per entry.** Name every chain of an entry in `structures.X.chains` (with `primary_chain`) to obtain it.


## Mapping

- Correspondence: pairs (267 positions; 267 present in both). Missing-residue policy: `exclude`.
- Residue-name mismatches (flag): A:250 GLY→A:250 ARG, A:261 SEP→A:261 SER.
- Alternate locations resolved: A:291 (A: A).
- Breaks (adapter): A [46], B [46] (row index after which the chain breaks).

## Global Cα RMSD (context only)

- 0.555 Å over 267 Cα after optimal superposition (MDAnalysis QCP); independent numpy SVD/Kabsch check 0.555 Å (|difference| 1.0e-13 Å).
- Published RMSD values often use other atom sets, residue ranges or outlier rejection; compare only with a stated definition.

## Largest local differences (top 10; descriptive)

| Metric | Unit | Windows (A numbering; value; * altloc, † residue-name mismatch in window) |
|---|---|---|
| d_kappa | rad | A:97..A:99 (0.806), A:98..A:100 (0.360), A:96..A:98 (0.359), A:33..A:35 (0.326), A:251..A:253 (0.264), A:49..A:51 (0.258), A:50..A:52 (0.185), A:103..A:105 (0.182), A:104..A:106 (0.160), A:59..A:61 (0.157) |
| d_i_i2 | Å | A:97..A:99 (1.557), A:59..A:61 (1.027), A:98..A:100 (0.881), A:96..A:98 (0.827), A:60..A:62 (0.613), A:251..A:253 (0.508), A:103..A:105 (0.407), A:49..A:51 (0.346), A:34..A:36 (0.320), A:99..A:101 (0.309) |
| rmsd3 | Å | A:97..A:99 (0.886), A:59..A:61 (0.502), A:33..A:35 (0.448), A:98..A:100 (0.436), A:96..A:98 (0.425), A:60..A:62 (0.403), A:251..A:253 (0.291), A:49..A:51 (0.259), A:103..A:105 (0.212), A:50..A:52 (0.190) |
| d_tau | rad | A:44..A:47 (0.582), A:248..A:251 (0.569†), A:95..A:98 (0.560), A:250..A:253 (0.482†), A:49..A:52 (0.449), A:43..A:46 (0.436), A:46..A:49 (0.418), A:45..A:48 (0.381), A:104..A:107 (0.380), A:249..A:252 (0.361†) |
| d_i_i3 | Å | A:96..A:99 (1.685), A:248..A:251 (1.617†), A:95..A:98 (1.494), A:49..A:52 (1.378), A:98..A:101 (0.877), A:60..A:63 (0.620), A:99..A:102 (0.576), A:43..A:46 (0.437), A:59..A:62 (0.429), A:45..A:48 (0.424) |
| rmsd4 | Å | A:96..A:99 (0.887), A:97..A:100 (0.788), A:95..A:98 (0.748), A:248..A:251 (0.674†), A:49..A:52 (0.624), A:98..A:101 (0.553), A:58..A:61 (0.469), A:59..A:62 (0.459), A:250..A:253 (0.449†), A:33..A:36 (0.423) |

## Agreement between equal-information metrics (Spearman over windows; descriptive)

| Pair | Spearman | Windows |
|---|---|---|
| d_kappa~d_i_i2 | 0.920 | 263 |
| d_kappa~rmsd3 | 0.977 | 263 |
| d_i_i2~rmsd3 | 0.967 | 263 |
| d_tau~d_i_i3 | 0.437 | 261 |
| d_tau~rmsd4 | 0.792 | 261 |
| d_i_i3~rmsd4 | 0.681 | 261 |

Three-residue metrics share one triangle, so close agreement is expected; d(i−1, i+2) is blind to the sign of τ.

## Exclusions

| Window (A) | Size | Structure | Cause | Location |
|---|---|---|---|---|
| A:78..A:86 | 3 | A | numbering_discontinuity | A:79→A:86 |
| A:78..A:86 | 3 | B | numbering_discontinuity | A:79→A:86 |
| A:79..A:87 | 3 | A | numbering_discontinuity | A:79→A:86 |
| A:79..A:87 | 3 | B | numbering_discontinuity | A:79→A:86 |
| A:77..A:86 | 4 | A | numbering_discontinuity | A:79→A:86 |
| A:77..A:86 | 4 | B | numbering_discontinuity | A:79→A:86 |
| A:78..A:87 | 4 | A | numbering_discontinuity | A:79→A:86 |
| A:78..A:87 | 4 | B | numbering_discontinuity | A:79→A:86 |
| A:79..A:88 | 4 | A | numbering_discontinuity | A:79→A:86 |
| A:79..A:88 | 4 | B | numbering_discontinuity | A:79→A:86 |

## Modified residues

Policy `phospho_pdb`: SEP, TPO and PTR in PDB-format files with MODRES records; no other modified residue or format. The deposited residue name is kept; the parent is documentation for the correspondence and is not a rename (SEP versus SER stays a visible name difference). A modification annotation is not evidence of functional state or of any effect.

| Structure | Residue | Deposited name | Parent (source) | Annotation | Record |
|---|---|---|---|---|---|
| A | A:261 | SEP | SER (MODRES record for A:261 (parent SER) agrees with the policy table SEP->SER (wwPDB chemical components)) | MODIFIED RESIDUE | HETATM |

## Alternate-conformation sensitivity

Policy `highest_occupancy`: largest change per metric d_kappa 0.0000, d_i_i2 0.0000, rmsd3 0.0000, d_tau 0.0000, d_i_i3 0.0000, rmsd4 0.0000; global RMSD 0.555 Å.

## Figure

![structure pair](figure.png)

Files: `windows.csv` (all metrics with units in column names), `residues.csv`, `exclusions.csv`, `summary.json`, `manifest.json` (hashes and software), `numeric_record.json` (computed values before display rounding). structpair 0.5.0.
