# Two-structure comparison: KIN-PAIR-1 EGFR-1: 8A27 vs 8A2D (EGFR)

**Descriptive, static.** Two static structures: the values are direct geometric differences between these coordinate sets. They can reflect crystal form and packing, crystallization conditions, resolution, refinement, alternate conformations and coordinate error as well as any biological difference. No ensemble, population, kinetic, functional-state or treatment-response claim follows. Agreement or disagreement between metrics does not show that one representation is better. Crystal copies of one entry are not independent replicates.

## Inputs

| | Label | File (SHA-256, first 12) | Chain | Altloc policy |
|---|---|---|---|---|
| A | 8A27 chain A | `8A27.pdb` (9fba519399f7) | A | blank_then_A |
| B | 8A2D chain A | `8A2D.pdb` (3fb105982609) | A | blank_then_A |

## Within-entry crystal-copy variation

**Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of them. The within-entry crystal-copy variation is a descriptive range of differences among the chains the user selected. It is not a calibrated noise floor, not a statistical uncertainty, not a nucleotide or ligand effect and not an HDX replicate estimate. An advisory flag appears only for a channel whose minimum difference the user configured; it is an advisory rule, not a gate, not an adequacy verdict, not significance and not biological separation.**

**Unavailable: one chain per entry.** Name every chain of an entry in `structures.X.chains` (with `primary_chain`) to obtain it.


## Mapping

- Correspondence: pairs (269 positions; 269 present in both). Missing-residue policy: `exclude`.
- Residue-name mismatches (flag): A:858 LEU→A:858 ARG, A:948 VAL→A:948 ARG.
- Alternate locations resolved: A:707 (A: A), A:713 (A: A), A:730 (A: A), A:745 (A: A), A:762 (B: A), A:775 (A: A), A:785 (B: A), A:804 (B: A), A:806 (A: A), A:812 (A: A), A:816 (A: A), A:832 (A: A), A:836 (A: A), A:841 (A: A), A:853 (A: A), A:908 (B: A), A:925 (A: A), A:939 (B: A), A:945 (B: A), A:953 (A: A), A:957 (A: A), A:960 (A: A), A:962 (A: A), A:963 (A: A), A:967 (A: A), A:971 (A: A), A:973 (A: A), A:974 (A: A), A:977 (A: A).
- Breaks (adapter): A [158], B [158] (row index after which the chain breaks).

## Global Cα RMSD (context only)

- 0.606 Å over 269 Cα after optimal superposition (MDAnalysis QCP); independent numpy SVD/Kabsch check 0.606 Å (|difference| 3.0e-13 Å).
- Published RMSD values often use other atom sets, residue ranges or outlier rejection; compare only with a stated definition.

## Largest local differences (top 10; descriptive)

| Metric | Unit | Windows (A numbering; value; * altloc, † residue-name mismatch in window) |
|---|---|---|
| d_kappa | rad | A:855..A:857 (0.504), A:856..A:858 (0.439†), A:854..A:856 (0.216), A:701..A:703 (0.199), A:915..A:917 (0.177), A:748..A:750 (0.134), A:768..A:770 (0.122), A:918..A:920 (0.101), A:877..A:879 (0.079), A:719..A:721 (0.073) |
| d_i_i2 | Å | A:855..A:857 (1.089), A:856..A:858 (0.906†), A:854..A:856 (0.484), A:915..A:917 (0.374), A:701..A:703 (0.371), A:768..A:770 (0.305), A:916..A:918 (0.232), A:857..A:859 (0.214†), A:877..A:879 (0.205), A:922..A:924 (0.158) |
| rmsd3 | Å | A:855..A:857 (0.582), A:856..A:858 (0.499†), A:854..A:856 (0.253), A:701..A:703 (0.219), A:915..A:917 (0.203), A:768..A:770 (0.151), A:748..A:750 (0.131), A:918..A:920 (0.106), A:916..A:918 (0.100), A:877..A:879 (0.099) |
| d_tau | rad | A:855..A:858 (2.939†), A:856..A:859 (1.625†), A:854..A:857 (1.623), A:722..A:725 (0.548), A:723..A:726 (0.501), A:915..A:918 (0.397), A:876..A:879 (0.342), A:914..A:917 (0.325), A:747..A:750 (0.320), A:916..A:919 (0.278) |
| d_i_i3 | Å | A:856..A:859 (2.409†), A:855..A:858 (2.172†), A:854..A:857 (1.599), A:853..A:856 (0.597*), A:767..A:770 (0.354), A:723..A:726 (0.324), A:913..A:916 (0.323), A:701..A:704 (0.315), A:766..A:769 (0.307), A:876..A:879 (0.275) |
| rmsd4 | Å | A:855..A:858 (1.888†), A:856..A:859 (1.658†), A:854..A:857 (1.610), A:914..A:917 (0.369), A:876..A:879 (0.347), A:915..A:918 (0.341), A:853..A:856 (0.269*), A:701..A:704 (0.267), A:723..A:726 (0.241), A:916..A:919 (0.240) |

## Agreement between equal-information metrics (Spearman over windows; descriptive)

| Pair | Spearman | Windows |
|---|---|---|
| d_kappa~d_i_i2 | 0.545 | 265 |
| d_kappa~rmsd3 | 0.814 | 265 |
| d_i_i2~rmsd3 | 0.771 | 265 |
| d_tau~d_i_i3 | 0.338 | 263 |
| d_tau~rmsd4 | 0.667 | 263 |
| d_i_i3~rmsd4 | 0.624 | 263 |

Three-residue metrics share one triangle, so close agreement is expected; d(i−1, i+2) is blind to the sign of τ.

## Exclusions

| Window (A) | Size | Structure | Cause | Location |
|---|---|---|---|---|
| A:858..A:876 | 3 | A | numbering_discontinuity | A:859→A:876 |
| A:858..A:876 | 3 | B | numbering_discontinuity | A:859→A:876 |
| A:859..A:877 | 3 | A | numbering_discontinuity | A:859→A:876 |
| A:859..A:877 | 3 | B | numbering_discontinuity | A:859→A:876 |
| A:857..A:876 | 4 | A | numbering_discontinuity | A:859→A:876 |
| A:857..A:876 | 4 | B | numbering_discontinuity | A:859→A:876 |
| A:858..A:877 | 4 | A | numbering_discontinuity | A:859→A:876 |
| A:858..A:877 | 4 | B | numbering_discontinuity | A:859→A:876 |
| A:859..A:878 | 4 | A | numbering_discontinuity | A:859→A:876 |
| A:859..A:878 | 4 | B | numbering_discontinuity | A:859→A:876 |

## Alternate-conformation sensitivity

Not requested.

## Figure

![structure pair](figure.png)

Files: `windows.csv` (all metrics with units in column names), `residues.csv`, `exclusions.csv`, `summary.json`, `manifest.json` (hashes and software). structpair 0.3.0.
