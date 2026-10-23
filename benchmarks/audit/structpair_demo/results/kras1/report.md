# Two-structure comparison: KRAS-STRUCT-1 reproduction: 9O0N (GDP) vs 9O0O (GMPPNP), wild-type KRAS, both with MRTX1133

**Descriptive, static.** Two static structures: the values are direct geometric differences between these coordinate sets. They can reflect crystal form and packing, crystallization conditions, resolution, refinement, alternate conformations and coordinate error as well as any biological difference. No ensemble, population, kinetic, functional-state or treatment-response claim follows. Agreement or disagreement between metrics does not show that one representation is better. Crystal copies of one entry are not independent replicates.

## Inputs

| | Label | File (SHA-256, first 12) | Chain | Altloc policy |
|---|---|---|---|---|
| A | 9O0N GDP | `9O0N.pdb` (b8ab08fedf0c) | A | blank_then_A |
| B | 9O0O GMPPNP | `9O0O.pdb` (8ac5d4be7827) | A | blank_then_A |

## Mapping

- Correspondence: same_numbering (168 positions; 168 present in both). Missing-residue policy: `reject`.
- Residue-name mismatches (reject): none.
- Alternate locations resolved: A:44 (A: A), A:79 (A: A), A:87 (A: A), A:106 (A: A), A:118 (B: A), A:126 (A: A), A:142 (A: A), A:154 (A: A).
- Breaks (adapter): A none, B none (row index after which the chain breaks).

## Global Cα RMSD (context only)

- 0.718 Å over 168 Cα after optimal superposition (MDAnalysis QCP); independent numpy SVD/Kabsch check 0.718 Å (|difference| 1.1e-13 Å).
- Published RMSD values often use other atom sets, residue ranges or outlier rejection; compare only with a stated definition.

## Largest local differences (top 10; descriptive)

| Metric | Unit | Windows (A numbering; value; * altloc, † residue-name mismatch in window) |
|---|---|---|
| d_kappa | rad | A:105..A:107 (0.528*), A:28..A:30 (0.468), A:31..A:33 (0.422), A:64..A:66 (0.395), A:30..A:32 (0.359), A:103..A:105 (0.347), A:63..A:65 (0.321), A:36..A:38 (0.302), A:35..A:37 (0.285), A:37..A:39 (0.216) |
| d_i_i2 | Å | A:28..A:30 (0.959), A:105..A:107 (0.866*), A:103..A:105 (0.741), A:36..A:38 (0.703), A:63..A:65 (0.683), A:31..A:33 (0.656), A:64..A:66 (0.633), A:30..A:32 (0.608), A:35..A:37 (0.529), A:34..A:36 (0.522) |
| rmsd3 | Å | A:105..A:107 (0.553*), A:28..A:30 (0.527), A:31..A:33 (0.435), A:64..A:66 (0.410), A:103..A:105 (0.399), A:30..A:32 (0.380), A:63..A:65 (0.365), A:36..A:38 (0.358), A:35..A:37 (0.311), A:34..A:36 (0.253) |
| d_tau | rad | A:36..A:39 (0.898), A:28..A:31 (0.894), A:31..A:34 (0.869), A:107..A:110 (0.720), A:59..A:62 (0.695), A:37..A:40 (0.654), A:32..A:35 (0.555), A:26..A:29 (0.547), A:30..A:33 (0.524), A:34..A:37 (0.521) |
| d_i_i3 | Å | A:102..A:105 (1.915), A:27..A:30 (0.960), A:103..A:106 (0.938*), A:30..A:33 (0.867), A:26..A:29 (0.821), A:63..A:66 (0.783), A:60..A:63 (0.745), A:32..A:35 (0.703), A:29..A:32 (0.638), A:62..A:65 (0.613) |
| rmsd4 | Å | A:28..A:31 (0.888), A:102..A:105 (0.764), A:35..A:38 (0.712), A:107..A:110 (0.666), A:105..A:108 (0.597*), A:30..A:33 (0.590), A:36..A:39 (0.588), A:103..A:106 (0.574*), A:32..A:35 (0.564), A:27..A:30 (0.535) |

## Agreement between equal-information metrics (Spearman over windows; descriptive)

| Pair | Spearman | Windows |
|---|---|---|
| d_kappa~d_i_i2 | 0.908 | 166 |
| d_kappa~rmsd3 | 0.973 | 166 |
| d_i_i2~rmsd3 | 0.950 | 166 |
| d_tau~d_i_i3 | 0.478 | 165 |
| d_tau~rmsd4 | 0.786 | 165 |
| d_i_i3~rmsd4 | 0.753 | 165 |

Three-residue metrics share one triangle, so close agreement is expected; d(i−1, i+2) is blind to the sign of τ.

## Exclusions

No window inside the mapped range was excluded.

## Alternate-conformation sensitivity

Policy `highest_occupancy`: largest change per metric d_kappa 0.0234, d_i_i2 0.0773, rmsd3 0.0337, d_tau 0.0450, d_i_i3 0.0975, rmsd4 0.0364; global RMSD 0.718 Å.

## Figure

![structure pair](figure.png)

Files: `windows.csv` (all metrics with units in column names), `residues.csv`, `exclusions.csv`, `summary.json`, `manifest.json` (hashes and software). structpair 0.1.0.
