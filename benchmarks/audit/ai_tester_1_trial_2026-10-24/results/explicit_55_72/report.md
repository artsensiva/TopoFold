# Two-structure comparison: AI tester 1 technical missing-residue case: explicit 55–72 identity mapping

**Descriptive, static.** Two static structures: the values are direct geometric differences between these coordinate sets. They can reflect crystal form and packing, crystallization conditions, resolution, refinement, alternate conformations and coordinate error as well as any biological difference. No ensemble, population, kinetic, functional-state or treatment-response claim follows. Agreement or disagreement between metrics does not show that one representation is better. Crystal copies of one entry are not independent replicates.

## Inputs

| | Label | File (SHA-256, first 12) | Chain | Altloc policy |
|---|---|---|---|---|
| A | 5US4 chain A | `5US4.pdb` (a1d0c9105887) | A | reject |
| B | 5USJ chain A | `5USJ.pdb` (dab6b79531b0) | A | reject |

## Mapping

- Correspondence: pairs (18 positions; 13 present in both). Missing-residue policy: `exclude`.
- Residue-name mismatches (reject): none.
- Alternate locations resolved: none.
- Breaks (adapter): A [4, 5, 6, 7, 8, 9], B [5, 6, 7, 8] (row index after which the chain breaks).

## Global Cα RMSD (context only)

- 1.919 Å over 13 Cα after optimal superposition (MDAnalysis QCP); independent numpy SVD/Kabsch check 1.919 Å (|difference| 1.8e-15 Å).
- Published RMSD values often use other atom sets, residue ranges or outlier rejection; compare only with a stated definition.

## Largest local differences (top 10; descriptive)

| Metric | Unit | Windows (A numbering; value; * altloc, † residue-name mismatch in window) |
|---|---|---|
| d_kappa | rad | A:57..A:59 (0.270), A:56..A:58 (0.162), A:65..A:67 (0.092), A:69..A:71 (0.078), A:67..A:69 (0.043), A:66..A:68 (0.029), A:68..A:70 (0.020), A:70..A:72 (0.017), A:55..A:57 (0.012) |
| d_i_i2 | Å | A:57..A:59 (0.460), A:56..A:58 (0.431), A:65..A:67 (0.220), A:69..A:71 (0.189), A:67..A:69 (0.154), A:66..A:68 (0.056), A:55..A:57 (0.034), A:68..A:70 (0.029), A:70..A:72 (0.028) |
| rmsd3 | Å | A:57..A:59 (0.288), A:56..A:58 (0.205), A:65..A:67 (0.111), A:69..A:71 (0.094), A:67..A:69 (0.066), A:66..A:68 (0.033), A:68..A:70 (0.026), A:70..A:72 (0.023), A:55..A:57 (0.019) |
| d_tau | rad | A:56..A:59 (0.186), A:69..A:72 (0.134), A:65..A:68 (0.097), A:68..A:71 (0.078), A:66..A:69 (0.029), A:67..A:70 (0.011), A:55..A:58 (0.008) |
| d_i_i3 | Å | A:56..A:59 (0.637), A:65..A:68 (0.558), A:55..A:58 (0.342), A:68..A:71 (0.335), A:66..A:69 (0.133), A:69..A:72 (0.095), A:67..A:70 (0.075) |
| rmsd4 | Å | A:56..A:59 (0.353), A:65..A:68 (0.207), A:55..A:58 (0.201), A:69..A:72 (0.137), A:68..A:71 (0.136), A:66..A:69 (0.071), A:67..A:70 (0.059) |

## Agreement between equal-information metrics (Spearman over windows; descriptive)

| Pair | Spearman | Windows |
|---|---|---|
| d_kappa~d_i_i2 | 0.950 | 9 |
| d_kappa~rmsd3 | 1.000 | 9 |
| d_i_i2~rmsd3 | 0.950 | 9 |
| d_tau~d_i_i3 | 0.393 | 7 |
| d_tau~rmsd4 | 0.571 | 7 |
| d_i_i3~rmsd4 | 0.893 | 7 |

Three-residue metrics share one triangle, so close agreement is expected; d(i−1, i+2) is blind to the sign of τ.

## Exclusions

| Window (A) | Size | Structure | Cause | Location |
|---|---|---|---|---|
| A:58..A:60 | 3 | A | missing_CA | A:59→A:60 |
| A:59..A:61 | 3 | A | missing_coordinates | A:60 |
| A:59..A:61 | 3 | A | missing_coordinates | A:61 |
| A:59..A:61 | 3 | B | missing_CA | A:60→A:61 |
| A:60..A:62 | 3 | A | missing_coordinates | A:60 |
| A:60..A:62 | 3 | A | missing_coordinates | A:61 |
| A:60..A:62 | 3 | A | missing_coordinates | A:62 |
| A:60..A:62 | 3 | B | missing_coordinates | A:61 |
| A:60..A:62 | 3 | B | missing_coordinates | A:62 |
| A:61..A:63 | 3 | A | missing_coordinates | A:61 |
| A:61..A:63 | 3 | A | missing_coordinates | A:62 |
| A:61..A:63 | 3 | A | missing_coordinates | A:63 |
| A:61..A:63 | 3 | B | missing_coordinates | A:61 |
| A:61..A:63 | 3 | B | missing_coordinates | A:62 |
| A:61..A:63 | 3 | B | missing_coordinates | A:63 |
| A:62..A:64 | 3 | A | missing_coordinates | A:62 |
| A:62..A:64 | 3 | A | missing_coordinates | A:63 |
| A:62..A:64 | 3 | A | missing_coordinates | A:64 |
| A:62..A:64 | 3 | B | missing_coordinates | A:62 |
| A:62..A:64 | 3 | B | missing_coordinates | A:63 |
| A:63..A:65 | 3 | A | missing_coordinates | A:63 |
| A:63..A:65 | 3 | A | missing_coordinates | A:64 |
| A:63..A:65 | 3 | B | missing_CA | A:63→A:64 |
| A:64..A:66 | 3 | A | missing_CA | A:64→A:65 |
| A:57..A:60 | 4 | A | missing_coordinates | A:60 |
| A:58..A:61 | 4 | A | missing_coordinates | A:60 |
| A:58..A:61 | 4 | A | missing_coordinates | A:61 |
| A:58..A:61 | 4 | B | missing_coordinates | A:61 |
| A:59..A:62 | 4 | A | missing_coordinates | A:60 |
| A:59..A:62 | 4 | A | missing_coordinates | A:61 |
| A:59..A:62 | 4 | A | missing_coordinates | A:62 |
| A:59..A:62 | 4 | B | missing_coordinates | A:61 |
| A:59..A:62 | 4 | B | missing_coordinates | A:62 |
| A:60..A:63 | 4 | A | missing_coordinates | A:60 |
| A:60..A:63 | 4 | A | missing_coordinates | A:61 |
| A:60..A:63 | 4 | A | missing_coordinates | A:62 |
| A:60..A:63 | 4 | A | missing_coordinates | A:63 |
| A:60..A:63 | 4 | B | missing_coordinates | A:61 |
| A:60..A:63 | 4 | B | missing_coordinates | A:62 |
| A:60..A:63 | 4 | B | missing_coordinates | A:63 |
| A:61..A:64 | 4 | A | missing_coordinates | A:61 |
| A:61..A:64 | 4 | A | missing_coordinates | A:62 |
| A:61..A:64 | 4 | A | missing_coordinates | A:63 |
| A:61..A:64 | 4 | A | missing_coordinates | A:64 |
| A:61..A:64 | 4 | B | missing_coordinates | A:61 |
| A:61..A:64 | 4 | B | missing_coordinates | A:62 |
| A:61..A:64 | 4 | B | missing_coordinates | A:63 |
| A:62..A:65 | 4 | A | missing_coordinates | A:62 |
| A:62..A:65 | 4 | A | missing_coordinates | A:63 |
| A:62..A:65 | 4 | A | missing_coordinates | A:64 |
| A:62..A:65 | 4 | B | missing_coordinates | A:62 |
| A:62..A:65 | 4 | B | missing_coordinates | A:63 |
| A:63..A:66 | 4 | A | missing_coordinates | A:63 |
| A:63..A:66 | 4 | A | missing_coordinates | A:64 |
| A:63..A:66 | 4 | B | missing_coordinates | A:63 |
| A:64..A:67 | 4 | A | missing_coordinates | A:64 |

## Alternate-conformation sensitivity

not applicable: no alternate locations at the mapped Cα atoms

## Figure

![structure pair](figure.png)

Files: `windows.csv` (all metrics with units in column names), `residues.csv`, `exclusions.csv`, `summary.json`, `manifest.json` (hashes and software). structpair 0.1.0.
