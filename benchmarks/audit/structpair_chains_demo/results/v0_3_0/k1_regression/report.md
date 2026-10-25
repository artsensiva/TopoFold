# Two-structure comparison: Regression of the recorded K1 crystal calculation: 5US4 A,B vs 5USJ A-F, residues 4-59 and 68-164, K1 settings

**Descriptive, static.** Two static structures: the values are direct geometric differences between these coordinate sets. They can reflect crystal form and packing, crystallization conditions, resolution, refinement, alternate conformations and coordinate error as well as any biological difference. No ensemble, population, kinetic, functional-state or treatment-response claim follows. Agreement or disagreement between metrics does not show that one representation is better. Crystal copies of one entry are not independent replicates.

## Inputs

| | Label | File (SHA-256, first 12) | Chain | Altloc policy |
|---|---|---|---|---|
| A | 5US4 GDP | `5US4.pdb` (a1d0c9105887) | A (primary of A, B) | reject |
| B | 5USJ GMPPNP | `5USJ.pdb` (dab6b79531b0) | A (primary of A, B, C, D, E, F) | reject |

## Within-entry crystal-copy variation

**Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of them. The within-entry crystal-copy variation is a descriptive range of differences among the chains the user selected. It is not a calibrated noise floor, not a statistical uncertainty, not a nucleotide or ligand effect and not an HDX replicate estimate. An advisory flag appears only for a channel whose minimum difference the user configured; it is an advisory rule, not a gate, not an adequacy verdict, not significance and not biological separation.**

- Status: **available**. Chains A ['A', 'B'] (primary A), B ['A', 'B', 'C', 'D', 'E', 'F'] (primary A). Pairs: 12 cross, 1 within A, 15 within B. The legacy tables above are the primary pair.
- The numeric between-entry and within-entry values are always given in `chain_variation.csv`; an empty cell is unavailable, never zero, and the `*_n` columns count the chain pairs that contributed.
- Advisory rule (not significance, not biological separation), configured by the user for disp (minimum difference 1 Å): the between-entry **mean** is at least max(minimum difference, 2 × the within-entry maximum); equality meets the rule and a value within 1e-9 of the threshold is marked as at the boundary.
- Displacement trim fraction 0.2.

| Channel | Unit | Positions with a value | Advisory rule | Rule met | Rule not met | Not assessed |
|---|---|---|---|---|---|---|
| disp | Å | 153 | configured | 8 | 145 | 0 |
| d_kappa | rad | 149 | not configured | 0 | 0 | 149 |
| d_i_i2 | Å | 149 | not configured | 0 | 0 | 149 |
| rmsd3 | Å | 149 | not configured | 0 | 0 | 149 |
| d_tau | rad | 147 | not configured | 0 | 0 | 147 |
| d_i_i3 | Å | 147 | not configured | 0 | 0 | 147 |
| rmsd4 | Å | 147 | not configured | 0 | 0 | 147 |

Largest between-entry displacement (median over cross pairings) with the within-entry maximum (top 10):

| A residue | Between median [min–max] (Å) | Within-entry max (Å) | Scope | Advisory rule |
|---|---|---|---|---|
| A:68 | 4.860 [3.681–5.281] | 4.130 | both | not met |
| A:32 | 4.509 [4.157–4.955] | 1.506 | both | met |
| A:69 | 3.739 [3.174–4.017] | 2.778 | both | not met |
| A:36 | 3.663 [3.461–3.890] | 0.456 | both | met |
| A:70 | 3.308 [3.163–3.558] | 2.256 | both | not met |
| A:35 | 3.051 [2.908–3.107] | 0.273 | both | met |
| A:71 | 2.514 [2.207–2.943] | 2.274 | both | not met |
| A:122 | 2.156 [1.349–2.302] | 1.354 | both | not met |
| A:74 | 1.946 [1.717–2.360] | 0.882 | both | met |
| A:33 | 1.594 [1.490–1.718] | 0.420 | both | met |

Files: `chain_variation.csv` (long format; every statistic with its count), `chains.csv`, `chain_variation.png`.


## Mapping

- Correspondence: same_numbering (153 positions; 153 present in both). Missing-residue policy: `reject`.
- Residue-name mismatches (reject): none.
- Alternate locations resolved: none.
- Breaks (adapter): A [55], B [55] (row index after which the chain breaks).

## Global Cα RMSD (context only)

- 0.973 Å over 153 Cα after optimal superposition (MDAnalysis QCP); independent numpy SVD/Kabsch check 0.973 Å (|difference| 5.6e-14 Å).
- Published RMSD values often use other atom sets, residue ranges or outlier rejection; compare only with a stated definition.

## Largest local differences (top 10; descriptive)

| Metric | Unit | Windows (A numbering; value; * altloc, † residue-name mismatch in window) |
|---|---|---|
| d_kappa | rad | A:34..A:36 (0.537), A:33..A:35 (0.414), A:51..A:53 (0.281), A:57..A:59 (0.270), A:35..A:37 (0.266), A:50..A:52 (0.264), A:36..A:38 (0.238), A:37..A:39 (0.230), A:73..A:75 (0.212), A:74..A:76 (0.203) |
| d_i_i2 | Å | A:34..A:36 (1.067), A:33..A:35 (0.924), A:35..A:37 (0.651), A:51..A:53 (0.623), A:73..A:75 (0.539), A:57..A:59 (0.460), A:36..A:38 (0.460), A:56..A:58 (0.431), A:25..A:27 (0.426), A:72..A:74 (0.411) |
| rmsd3 | Å | A:34..A:36 (0.599), A:33..A:35 (0.484), A:51..A:53 (0.326), A:35..A:37 (0.323), A:57..A:59 (0.288), A:50..A:52 (0.271), A:36..A:38 (0.264), A:73..A:75 (0.263), A:37..A:39 (0.219), A:74..A:76 (0.212) |
| d_tau | rad | A:33..A:36 (2.894), A:35..A:38 (2.497), A:31..A:34 (1.812), A:34..A:37 (1.277), A:29..A:32 (1.182), A:30..A:33 (1.009), A:37..A:40 (0.772), A:50..A:53 (0.712), A:72..A:75 (0.696), A:49..A:52 (0.645) |
| d_i_i3 | Å | A:33..A:36 (4.082), A:35..A:38 (2.715), A:31..A:34 (1.802), A:29..A:32 (1.563), A:34..A:37 (1.208), A:73..A:76 (0.775), A:56..A:59 (0.637), A:32..A:35 (0.602), A:24..A:27 (0.516), A:36..A:39 (0.492) |
| rmsd4 | Å | A:33..A:36 (2.334), A:35..A:38 (2.004), A:31..A:34 (1.641), A:34..A:37 (1.289), A:29..A:32 (0.906), A:30..A:33 (0.870), A:72..A:75 (0.629), A:36..A:39 (0.606), A:37..A:40 (0.528), A:73..A:76 (0.490) |

## Agreement between equal-information metrics (Spearman over windows; descriptive)

| Pair | Spearman | Windows |
|---|---|---|
| d_kappa~d_i_i2 | 0.789 | 149 |
| d_kappa~rmsd3 | 0.934 | 149 |
| d_i_i2~rmsd3 | 0.910 | 149 |
| d_tau~d_i_i3 | 0.378 | 147 |
| d_tau~rmsd4 | 0.817 | 147 |
| d_i_i3~rmsd4 | 0.622 | 147 |

Three-residue metrics share one triangle, so close agreement is expected; d(i−1, i+2) is blind to the sign of τ.

## Exclusions

| Window (A) | Size | Structure | Cause | Location |
|---|---|---|---|---|
| A:58..A:68 | 3 | A | numbering_discontinuity | A:59→A:68 |
| A:58..A:68 | 3 | B | numbering_discontinuity | A:59→A:68 |
| A:59..A:69 | 3 | A | numbering_discontinuity | A:59→A:68 |
| A:59..A:69 | 3 | B | numbering_discontinuity | A:59→A:68 |
| A:57..A:68 | 4 | A | numbering_discontinuity | A:59→A:68 |
| A:57..A:68 | 4 | B | numbering_discontinuity | A:59→A:68 |
| A:58..A:69 | 4 | A | numbering_discontinuity | A:59→A:68 |
| A:58..A:69 | 4 | B | numbering_discontinuity | A:59→A:68 |
| A:59..A:70 | 4 | A | numbering_discontinuity | A:59→A:68 |
| A:59..A:70 | 4 | B | numbering_discontinuity | A:59→A:68 |

## Alternate-conformation sensitivity

Not requested.

## Figure

![structure pair](figure.png)

Files: `windows.csv` (all metrics with units in column names), `residues.csv`, `exclusions.csv`, `summary.json`, `manifest.json` (hashes and software). structpair 0.3.0.
