# Two-structure comparison: Input-handling demonstration (values not interpreted): 9O0N chain A (WT) vs 5US4 chain B (G12D), residue-name mismatch flagged

**Descriptive, static.** Two static structures: the values are direct geometric differences between these coordinate sets. They can reflect crystal form and packing, crystallization conditions, resolution, refinement, alternate conformations and coordinate error as well as any biological difference. No ensemble, population, kinetic, functional-state or treatment-response claim follows. Agreement or disagreement between metrics does not show that one representation is better. Crystal copies of one entry are not independent replicates.

## Inputs

| | Label | File (SHA-256, first 12) | Chain | Altloc policy |
|---|---|---|---|---|
| A | 9O0N A WT GDP | `9O0N.pdb` (b8ab08fedf0c) | A | blank_then_A |
| B | 5US4 B G12D GDP | `5US4.pdb` (a1d0c9105887) | B | reject |

## Mapping

- Correspondence: same_numbering (166 positions; 159 present in both). Missing-residue policy: `exclude`.
- Residue-name mismatches (flag): A:12 GLY→B:12 ASP.
- Alternate locations resolved: A:44 (A: A), A:79 (A: A), A:87 (A: A), A:106 (A: A), A:126 (A: A), A:142 (A: A), A:154 (A: A).
- Breaks (adapter): A none, B [58, 59, 60, 61, 62, 63, 64, 65] (row index after which the chain breaks).

## Global Cα RMSD (context only)

- 0.932 Å over 159 Cα after optimal superposition (MDAnalysis QCP); independent numpy SVD/Kabsch check 0.932 Å (|difference| 2.0e-14 Å).
- Published RMSD values often use other atom sets, residue ranges or outlier rejection; compare only with a stated definition.

## Largest local differences (top 5; descriptive)

| Metric | Unit | Windows (A numbering; value; * altloc, † residue-name mismatch in window) |
|---|---|---|
| d_kappa | rad | A:120..A:122 (0.631), A:50..A:52 (0.202), A:51..A:53 (0.184), A:31..A:33 (0.174), A:73..A:75 (0.159) |
| d_i_i2 | Å | A:120..A:122 (1.355), A:73..A:75 (0.417), A:51..A:53 (0.383), A:72..A:74 (0.332), A:31..A:33 (0.306) |
| rmsd3 | Å | A:120..A:122 (0.723), A:51..A:53 (0.208), A:73..A:75 (0.200), A:50..A:52 (0.196), A:31..A:33 (0.186) |
| d_tau | rad | A:120..A:123 (1.214), A:122..A:125 (0.754), A:121..A:124 (0.676), A:119..A:122 (0.628), A:72..A:75 (0.594) |
| d_i_i3 | Å | A:31..A:34 (0.880), A:73..A:76 (0.725), A:119..A:122 (0.669), A:71..A:74 (0.562), A:34..A:37 (0.504) |
| rmsd4 | Å | A:120..A:123 (0.963), A:119..A:122 (0.732), A:122..A:125 (0.660), A:121..A:124 (0.580), A:72..A:75 (0.525) |

## Agreement between equal-information metrics (Spearman over windows; descriptive)

| Pair | Spearman | Windows |
|---|---|---|
| d_kappa~d_i_i2 | 0.799 | 155 |
| d_kappa~rmsd3 | 0.931 | 155 |
| d_i_i2~rmsd3 | 0.882 | 155 |
| d_tau~d_i_i3 | 0.451 | 153 |
| d_tau~rmsd4 | 0.737 | 153 |
| d_i_i3~rmsd4 | 0.644 | 153 |

Three-residue metrics share one triangle, so close agreement is expected; d(i−1, i+2) is blind to the sign of τ.

## Exclusions

| Window (A) | Size | Structure | Cause | Location |
|---|---|---|---|---|
| A:58..A:60 | 3 | B | missing_CA | B:59→B:60 |
| A:59..A:61 | 3 | B | missing_coordinates | B:60 |
| A:59..A:61 | 3 | B | missing_coordinates | B:61 |
| A:60..A:62 | 3 | B | missing_coordinates | B:60 |
| A:60..A:62 | 3 | B | missing_coordinates | B:61 |
| A:60..A:62 | 3 | B | missing_coordinates | B:62 |
| A:61..A:63 | 3 | B | missing_coordinates | B:61 |
| A:61..A:63 | 3 | B | missing_coordinates | B:62 |
| A:61..A:63 | 3 | B | missing_coordinates | B:63 |
| A:62..A:64 | 3 | B | missing_coordinates | B:62 |
| A:62..A:64 | 3 | B | missing_coordinates | B:63 |
| A:62..A:64 | 3 | B | missing_coordinates | B:64 |
| A:63..A:65 | 3 | B | missing_coordinates | B:63 |
| A:63..A:65 | 3 | B | missing_coordinates | B:64 |
| A:63..A:65 | 3 | B | missing_coordinates | B:65 |
| A:64..A:66 | 3 | B | missing_coordinates | B:64 |
| A:64..A:66 | 3 | B | missing_coordinates | B:65 |
| A:64..A:66 | 3 | B | missing_coordinates | B:66 |
| A:65..A:67 | 3 | B | missing_coordinates | B:65 |
| A:65..A:67 | 3 | B | missing_coordinates | B:66 |
| A:66..A:68 | 3 | B | missing_CA | B:66→B:67 |
| A:57..A:60 | 4 | B | missing_coordinates | B:60 |
| A:58..A:61 | 4 | B | missing_coordinates | B:60 |
| A:58..A:61 | 4 | B | missing_coordinates | B:61 |
| A:59..A:62 | 4 | B | missing_coordinates | B:60 |
| A:59..A:62 | 4 | B | missing_coordinates | B:61 |
| A:59..A:62 | 4 | B | missing_coordinates | B:62 |
| A:60..A:63 | 4 | B | missing_coordinates | B:60 |
| A:60..A:63 | 4 | B | missing_coordinates | B:61 |
| A:60..A:63 | 4 | B | missing_coordinates | B:62 |
| A:60..A:63 | 4 | B | missing_coordinates | B:63 |
| A:61..A:64 | 4 | B | missing_coordinates | B:61 |
| A:61..A:64 | 4 | B | missing_coordinates | B:62 |
| A:61..A:64 | 4 | B | missing_coordinates | B:63 |
| A:61..A:64 | 4 | B | missing_coordinates | B:64 |
| A:62..A:65 | 4 | B | missing_coordinates | B:62 |
| A:62..A:65 | 4 | B | missing_coordinates | B:63 |
| A:62..A:65 | 4 | B | missing_coordinates | B:64 |
| A:62..A:65 | 4 | B | missing_coordinates | B:65 |
| A:63..A:66 | 4 | B | missing_coordinates | B:63 |
| A:63..A:66 | 4 | B | missing_coordinates | B:64 |
| A:63..A:66 | 4 | B | missing_coordinates | B:65 |
| A:63..A:66 | 4 | B | missing_coordinates | B:66 |
| A:64..A:67 | 4 | B | missing_coordinates | B:64 |
| A:64..A:67 | 4 | B | missing_coordinates | B:65 |
| A:64..A:67 | 4 | B | missing_coordinates | B:66 |
| A:65..A:68 | 4 | B | missing_coordinates | B:65 |
| A:65..A:68 | 4 | B | missing_coordinates | B:66 |
| A:66..A:69 | 4 | B | missing_coordinates | B:66 |

## Alternate-conformation sensitivity

Not requested.

## Figure

![structure pair](figure.png)

Files: `windows.csv` (all metrics with units in column names), `residues.csv`, `exclusions.csv`, `summary.json`, `manifest.json` (hashes and software). structpair 0.1.0.
