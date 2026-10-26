# Two-structure comparison: KIN-PAIR-1 EGFR-2: 5UG9 vs 5HG8 (EGFR)

**Descriptive, static.** Two static structures: the values are direct geometric differences between these coordinate sets. They can reflect crystal form and packing, crystallization conditions, resolution, refinement, alternate conformations and coordinate error as well as any biological difference. No ensemble, population, kinetic, functional-state or treatment-response claim follows. Agreement or disagreement between metrics does not show that one representation is better. Crystal copies of one entry are not independent replicates.

## Inputs

| | Label | File (SHA-256, first 12) | Chain | Altloc policy |
|---|---|---|---|---|
| A | 5UG9 chain A | `5UG9.pdb` (24013b98bdeb) | A | blank_then_A |
| B | 5HG8 chain A | `5HG8.pdb` (8556a280ab40) | A | blank_then_A |

## Within-entry crystal-copy variation

**Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of them. The within-entry crystal-copy variation is a descriptive range of differences among the chains the user selected. It is not a calibrated noise floor, not a statistical uncertainty, not a nucleotide or ligand effect and not an HDX replicate estimate. An advisory flag appears only for a channel whose minimum difference the user configured; it is an advisory rule, not a gate, not an adequacy verdict, not significance and not biological separation.**

**Unavailable: one chain per entry.** Name every chain of an entry in `structures.X.chains` (with `primary_chain`) to obtain it.


## Mapping

- Correspondence: pairs (280 positions; 280 present in both). Missing-residue policy: `exclude`.
- Residue-name mismatches (flag): none.
- Alternate locations resolved: A:766 (B: A), A:767 (A: A), A:790 (A: A), A:802 (B: A), A:804 (A: A), A:834 (B: A), A:854 (B: A), A:881 (A: A), A:925 (A: A), A:952 (B: A), A:953 (A: A).
- Breaks (adapter): A [46], B [46] (row index after which the chain breaks).

## Global Cα RMSD (context only)

- 0.590 Å over 280 Cα after optimal superposition (MDAnalysis QCP); independent numpy SVD/Kabsch check 0.590 Å (|difference| 9.4e-14 Å).
- Published RMSD values often use other atom sets, residue ranges or outlier rejection; compare only with a stated definition.

## Largest local differences (top 10; descriptive)

| Metric | Unit | Windows (A numbering; value; * altloc, † residue-name mismatch in window) |
|---|---|---|
| d_kappa | rad | A:872..A:874 (0.471), A:873..A:875 (0.422), A:737..A:739 (0.189), A:713..A:715 (0.165), A:735..A:737 (0.132), A:852..A:854 (0.109*), A:703..A:705 (0.105), A:710..A:712 (0.102), A:935..A:937 (0.100), A:939..A:941 (0.099) |
| d_i_i2 | Å | A:872..A:874 (0.949), A:983..A:985 (0.787), A:873..A:875 (0.619), A:713..A:715 (0.362), A:737..A:739 (0.339), A:735..A:737 (0.311), A:710..A:712 (0.243), A:703..A:705 (0.222), A:939..A:941 (0.220), A:935..A:937 (0.220) |
| rmsd3 | Å | A:872..A:874 (0.525), A:983..A:985 (0.460), A:873..A:875 (0.429), A:737..A:739 (0.204), A:713..A:715 (0.191), A:735..A:737 (0.157), A:710..A:712 (0.122), A:852..A:854 (0.120*), A:703..A:705 (0.119), A:935..A:937 (0.116) |
| d_tau | rad | A:873..A:876 (1.063), A:872..A:875 (0.858), A:874..A:877 (0.531), A:875..A:878 (0.378), A:866..A:869 (0.364), A:753..A:756 (0.349), A:867..A:870 (0.347), A:870..A:873 (0.290), A:871..A:874 (0.242), A:702..A:705 (0.231) |
| d_i_i3 | Å | A:871..A:874 (1.694), A:982..A:985 (0.831), A:872..A:875 (0.734), A:761..A:764 (0.494), A:710..A:713 (0.479), A:806..A:809 (0.397), A:866..A:869 (0.390), A:912..A:915 (0.379), A:745..A:748 (0.325), A:737..A:740 (0.318) |
| rmsd4 | Å | A:871..A:874 (0.720), A:872..A:875 (0.634), A:873..A:876 (0.597), A:982..A:985 (0.427), A:874..A:877 (0.382), A:875..A:878 (0.257), A:753..A:756 (0.256), A:866..A:869 (0.222), A:867..A:870 (0.214), A:713..A:716 (0.205) |

## Agreement between equal-information metrics (Spearman over windows; descriptive)

| Pair | Spearman | Windows |
|---|---|---|
| d_kappa~d_i_i2 | 0.944 | 276 |
| d_kappa~rmsd3 | 0.989 | 276 |
| d_i_i2~rmsd3 | 0.970 | 276 |
| d_tau~d_i_i3 | 0.406 | 274 |
| d_tau~rmsd4 | 0.778 | 274 |
| d_i_i3~rmsd4 | 0.677 | 274 |

Three-residue metrics share one triangle, so close agreement is expected; d(i−1, i+2) is blind to the sign of τ.

## Exclusions

| Window (A) | Size | Structure | Cause | Location |
|---|---|---|---|---|
| A:747..A:753 | 3 | A | numbering_discontinuity | A:748→A:753 |
| A:747..A:753 | 3 | B | numbering_discontinuity | A:748→A:753 |
| A:748..A:754 | 3 | A | numbering_discontinuity | A:748→A:753 |
| A:748..A:754 | 3 | B | numbering_discontinuity | A:748→A:753 |
| A:746..A:753 | 4 | A | numbering_discontinuity | A:748→A:753 |
| A:746..A:753 | 4 | B | numbering_discontinuity | A:748→A:753 |
| A:747..A:754 | 4 | A | numbering_discontinuity | A:748→A:753 |
| A:747..A:754 | 4 | B | numbering_discontinuity | A:748→A:753 |
| A:748..A:755 | 4 | A | numbering_discontinuity | A:748→A:753 |
| A:748..A:755 | 4 | B | numbering_discontinuity | A:748→A:753 |

## Alternate-conformation sensitivity

Not requested.

## Figure

![structure pair](figure.png)

Files: `windows.csv` (all metrics with units in column names), `residues.csv`, `exclusions.csv`, `summary.json`, `manifest.json` (hashes and software). structpair 0.3.0.
