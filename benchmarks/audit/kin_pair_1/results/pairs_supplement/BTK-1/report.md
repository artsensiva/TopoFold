# Two-structure comparison: KIN-PAIR-1 BTK-1: 6E4F vs 4RFZ (BTK)

**Descriptive, static.** Two static structures: the values are direct geometric differences between these coordinate sets. They can reflect crystal form and packing, crystallization conditions, resolution, refinement, alternate conformations and coordinate error as well as any biological difference. No ensemble, population, kinetic, functional-state or treatment-response claim follows. Agreement or disagreement between metrics does not show that one representation is better. Crystal copies of one entry are not independent replicates.

## Inputs

| | Label | File (SHA-256, first 12) | Chain | Altloc policy |
|---|---|---|---|---|
| A | 6E4F chain A | `6E4F.pdb` (ad770d384948) | A | blank_then_A |
| B | 4RFZ chain A | `4RFZ.pdb` (bc8dc4307281) | A | blank_then_A |

## Within-entry crystal-copy variation

**Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of them. The within-entry crystal-copy variation is a descriptive range of differences among the chains the user selected. It is not a calibrated noise floor, not a statistical uncertainty, not a nucleotide or ligand effect and not an HDX replicate estimate. An advisory flag appears only for a channel whose minimum difference the user configured; it is an advisory rule, not a gate, not an adequacy verdict, not significance and not biological separation.**

**Unavailable: one chain per entry.** Name every chain of an entry in `structures.X.chains` (with `primary_chain`) to obtain it.


## Mapping

- Correspondence: pairs (263 positions; 263 present in both). Missing-residue policy: `exclude`.
- Residue-name mismatches (flag): A:489 MET→A:489 ALA, A:492 ARG→A:492 ALA, A:624 GLU→A:624 ALA, A:625 LYS→A:625 ALA.
- Alternate locations resolved: A:453 (A: A), A:484 (A: A), A:571 (A: A), A:596 (A: A).
- Breaks (adapter): A none, B none (row index after which the chain breaks).

## Global Cα RMSD (context only)

- 0.706 Å over 263 Cα after optimal superposition (MDAnalysis QCP); independent numpy SVD/Kabsch check 0.706 Å (|difference| 2.1e-13 Å).
- Published RMSD values often use other atom sets, residue ranges or outlier rejection; compare only with a stated definition.

## Largest local differences (top 10; descriptive)

| Metric | Unit | Windows (A numbering; value; * altloc, † residue-name mismatch in window) |
|---|---|---|
| d_kappa | rad | A:552..A:554 (0.403), A:413..A:415 (0.271), A:616..A:618 (0.206), A:551..A:553 (0.161), A:465..A:467 (0.158), A:410..A:412 (0.151), A:617..A:619 (0.140), A:451..A:453 (0.132*), A:553..A:555 (0.131), A:424..A:426 (0.131) |
| d_i_i2 | Å | A:552..A:554 (0.539), A:465..A:467 (0.446), A:571..A:573 (0.325*), A:413..A:415 (0.284), A:451..A:453 (0.273*), A:551..A:553 (0.271), A:616..A:618 (0.271), A:611..A:613 (0.269), A:487..A:489 (0.268†), A:617..A:619 (0.263) |
| rmsd3 | Å | A:552..A:554 (0.406), A:413..A:415 (0.259), A:616..A:618 (0.206), A:465..A:467 (0.206), A:551..A:553 (0.172), A:571..A:573 (0.159*), A:617..A:619 (0.153), A:451..A:453 (0.151*), A:553..A:555 (0.148), A:422..A:424 (0.143) |
| d_tau | rad | A:408..A:411 (1.271), A:413..A:416 (1.182), A:407..A:410 (1.007), A:551..A:554 (0.889), A:552..A:555 (0.608), A:409..A:412 (0.604), A:550..A:553 (0.472), A:616..A:619 (0.446), A:412..A:415 (0.441), A:537..A:540 (0.410) |
| d_i_i3 | Å | A:550..A:553 (1.079), A:413..A:416 (0.970), A:408..A:411 (0.828), A:553..A:556 (0.762), A:444..A:447 (0.712), A:465..A:468 (0.492), A:464..A:467 (0.462), A:450..A:453 (0.404*), A:430..A:433 (0.392), A:410..A:413 (0.377) |
| rmsd4 | Å | A:408..A:411 (0.727), A:413..A:416 (0.626), A:551..A:554 (0.572), A:552..A:555 (0.523), A:550..A:553 (0.500), A:407..A:410 (0.425), A:409..A:412 (0.362), A:615..A:618 (0.327), A:412..A:415 (0.320), A:444..A:447 (0.297) |

## Agreement between equal-information metrics (Spearman over windows; descriptive)

| Pair | Spearman | Windows |
|---|---|---|
| d_kappa~d_i_i2 | 0.706 | 261 |
| d_kappa~rmsd3 | 0.883 | 261 |
| d_i_i2~rmsd3 | 0.861 | 261 |
| d_tau~d_i_i3 | 0.445 | 260 |
| d_tau~rmsd4 | 0.760 | 260 |
| d_i_i3~rmsd4 | 0.729 | 260 |

Three-residue metrics share one triangle, so close agreement is expected; d(i−1, i+2) is blind to the sign of τ.

## Exclusions

No window inside the mapped range was excluded.

## Alternate-conformation sensitivity

Not requested.

## Figure

![structure pair](figure.png)

Files: `windows.csv` (all metrics with units in column names), `residues.csv`, `exclusions.csv`, `summary.json`, `manifest.json` (hashes and software). structpair 0.3.0.
