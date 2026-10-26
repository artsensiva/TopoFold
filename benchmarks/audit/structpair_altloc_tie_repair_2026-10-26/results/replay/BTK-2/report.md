# Two-structure comparison: KIN-PAIR-1 BTK-2: 7L5O vs 6J6M (BTK)

**Descriptive, static.** Two static structures: the values are direct geometric differences between these coordinate sets. They can reflect crystal form and packing, crystallization conditions, resolution, refinement, alternate conformations and coordinate error as well as any biological difference. No ensemble, population, kinetic, functional-state or treatment-response claim follows. Agreement or disagreement between metrics does not show that one representation is better. Crystal copies of one entry are not independent replicates.

## Inputs

| | Label | File (SHA-256, first 12) | Chain | Altloc policy |
|---|---|---|---|---|
| A | 7L5O chain A | `7L5O.pdb` (f3affba350c4) | A | blank_then_A |
| B | 6J6M chain A | `6J6M.pdb` (f8e8fe158f09) | A | blank_then_A |

## Within-entry crystal-copy variation

**Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of them. The within-entry crystal-copy variation is a descriptive range of differences among the chains the user selected. It is not a calibrated noise floor, not a statistical uncertainty, not a nucleotide or ligand effect and not an HDX replicate estimate. An advisory flag appears only for a channel whose minimum difference the user configured; it is an advisory rule, not a gate, not an adequacy verdict, not significance and not biological separation.**

**Unavailable: one chain per entry.** Name every chain of an entry in `structures.X.chains` (with `primary_chain`) to obtain it.


## Mapping

- Correspondence: pairs (263 positions; 263 present in both). Missing-residue policy: `exclude`.
- Residue-name mismatches (flag): none.
- Alternate locations resolved: A:400 (A: A), A:431 (A: A), A:437 (A: A), A:444 (A: A), A:450 (A: A), A:453 (A: A), A:455 (A: A), A:459 (A: A), A:460 (A: A), A:463 (A: A), A:470 (A: A), A:481 (A: A), A:484 (A: A), A:488 (A: A), A:492 (A: A), A:496 (A: A), A:500 (A: A), A:510 (A: A), A:513 (A: A), A:516 (A: A), A:571 (A: A), A:586 (A: A), A:595 (A: A), A:596 (A: A), A:600 (A: A), A:604 (A: A), A:608 (A: A), A:624 (A: A), A:635 (A: A), A:639 (A: A), A:640 (A: A), A:649 (A: A), A:655 (A: A), A:658 (A: A).
- Breaks (adapter): A [17], B [17] (row index after which the chain breaks).

## Global Cα RMSD (context only)

- 0.275 Å over 263 Cα after optimal superposition (MDAnalysis QCP); independent numpy SVD/Kabsch check 0.275 Å (|difference| 1.5e-13 Å).
- Published RMSD values often use other atom sets, residue ranges or outlier rejection; compare only with a stated definition.

## Largest local differences (top 10; descriptive)

| Metric | Unit | Windows (A numbering; value; * altloc, † residue-name mismatch in window) |
|---|---|---|
| d_kappa | rad | A:414..A:416 (0.290), A:466..A:468 (0.261), A:408..A:410 (0.128), A:432..A:434 (0.117), A:467..A:469 (0.108), A:407..A:409 (0.090), A:429..A:431 (0.086*), A:468..A:470 (0.079*), A:430..A:432 (0.072*), A:451..A:453 (0.067*) |
| d_i_i2 | Å | A:414..A:416 (0.464), A:466..A:468 (0.270), A:407..A:409 (0.197), A:432..A:434 (0.186), A:408..A:410 (0.184), A:468..A:470 (0.169*), A:465..A:467 (0.165), A:463..A:465 (0.164*), A:430..A:432 (0.158*), A:429..A:431 (0.157*) |
| rmsd3 | Å | A:414..A:416 (0.303), A:466..A:468 (0.250), A:408..A:410 (0.129), A:432..A:434 (0.123), A:407..A:409 (0.104), A:467..A:469 (0.098), A:429..A:431 (0.094*), A:430..A:432 (0.084*), A:468..A:470 (0.083*), A:463..A:465 (0.080*) |
| d_tau | rad | A:393..A:396 (0.346), A:407..A:410 (0.290), A:431..A:434 (0.272*), A:479..A:482 (0.248*), A:466..A:469 (0.242), A:432..A:435 (0.220), A:467..A:470 (0.202*), A:480..A:483 (0.154*), A:430..A:433 (0.141*), A:478..A:481 (0.132*) |
| d_i_i3 | Å | A:393..A:396 (0.561), A:465..A:468 (0.413), A:431..A:434 (0.350*), A:433..A:436 (0.319), A:554..A:557 (0.272), A:481..A:484 (0.263*), A:467..A:470 (0.258*), A:414..A:417 (0.257), A:546..A:549 (0.244), A:548..A:551 (0.227) |
| rmsd4 | Å | A:414..A:417 (0.353), A:393..A:396 (0.334), A:431..A:434 (0.287*), A:465..A:468 (0.244), A:466..A:469 (0.221), A:407..A:410 (0.182), A:432..A:435 (0.166), A:430..A:433 (0.166*), A:467..A:470 (0.159*), A:429..A:432 (0.141*) |

## Agreement between equal-information metrics (Spearman over windows; descriptive)

| Pair | Spearman | Windows |
|---|---|---|
| d_kappa~d_i_i2 | 0.726 | 259 |
| d_kappa~rmsd3 | 0.890 | 259 |
| d_i_i2~rmsd3 | 0.859 | 259 |
| d_tau~d_i_i3 | 0.390 | 257 |
| d_tau~rmsd4 | 0.766 | 257 |
| d_i_i3~rmsd4 | 0.659 | 257 |

Three-residue metrics share one triangle, so close agreement is expected; d(i−1, i+2) is blind to the sign of τ.

## Exclusions

| Window (A) | Size | Structure | Cause | Location |
|---|---|---|---|---|
| A:409..A:414 | 3 | A | numbering_discontinuity | A:410→A:414 |
| A:409..A:414 | 3 | B | numbering_discontinuity | A:410→A:414 |
| A:410..A:415 | 3 | A | numbering_discontinuity | A:410→A:414 |
| A:410..A:415 | 3 | B | numbering_discontinuity | A:410→A:414 |
| A:408..A:414 | 4 | A | numbering_discontinuity | A:410→A:414 |
| A:408..A:414 | 4 | B | numbering_discontinuity | A:410→A:414 |
| A:409..A:415 | 4 | A | numbering_discontinuity | A:410→A:414 |
| A:409..A:415 | 4 | B | numbering_discontinuity | A:410→A:414 |
| A:410..A:416 | 4 | A | numbering_discontinuity | A:410→A:414 |
| A:410..A:416 | 4 | B | numbering_discontinuity | A:410→A:414 |

## Alternate-conformation sensitivity

**UNAVAILABLE, not zero.** The requested sensitivity (policy `highest_occupancy`) was not computed: the `highest_occupancy` policy cannot choose between alternate locations of equal occupancy at 31 residue(s); the primary comparison (primary altloc policy) is complete and unaffected. A partial sensitivity would need a new rule for the tied residues (choose one or skip it), which changes the windows around them and the global fit, so the whole optional sensitivity is unavailable. This does not mean the structures are insensitive to the alternate-location choice. No value is reported.

| Structure | Residue | In the correspondence | Alternate locations | Occupancies |
|---|---|---|---|---|
| A | A:400 | yes | A, B | 0.5, 0.5 |
| A | A:431 | yes | A, B | 0.5, 0.5 |
| A | A:437 | yes | A, B | 0.5, 0.5 |
| A | A:453 | yes | A, B | 0.5, 0.5 |
| A | A:455 | yes | A, B | 0.5, 0.5 |
| A | A:459 | yes | A, B | 0.5, 0.5 |
| A | A:460 | yes | A, B | 0.5, 0.5 |
| A | A:463 | yes | A, B | 0.5, 0.5 |
| A | A:470 | yes | A, B | 0.5, 0.5 |
| A | A:484 | yes | A, B | 0.5, 0.5 |
| A | A:488 | yes | A, B | 0.5, 0.5 |
| A | A:492 | yes | A, B | 0.5, 0.5 |
| A | A:496 | yes | A, B | 0.5, 0.5 |
| A | A:500 | yes | A, B | 0.5, 0.5 |
| A | A:510 | yes | A, B | 0.5, 0.5 |
| A | A:513 | yes | A, B | 0.5, 0.5 |
| A | A:516 | yes | A, B | 0.5, 0.5 |
| A | A:571 | yes | A, B | 0.5, 0.5 |
| A | A:586 | yes | A, B | 0.5, 0.5 |
| A | A:595 | yes | A, B | 0.5, 0.5 |
| A | A:596 | yes | A, B | 0.5, 0.5 |
| A | A:600 | yes | A, B | 0.5, 0.5 |
| A | A:604 | yes | A, B | 0.5, 0.5 |
| A | A:608 | yes | A, B | 0.5, 0.5 |
| A | A:624 | yes | A, B | 0.5, 0.5 |
| A | A:635 | yes | A, B | 0.5, 0.5 |
| A | A:639 | yes | A, B | 0.5, 0.5 |
| A | A:640 | yes | A, B | 0.5, 0.5 |
| A | A:649 | yes | A, B | 0.5, 0.5 |
| A | A:655 | yes | A, B | 0.5, 0.5 |
| A | A:658 | yes | A, B | 0.5, 0.5 |

## Figure

![structure pair](figure.png)

Files: `windows.csv` (all metrics with units in column names), `residues.csv`, `exclusions.csv`, `summary.json`, `manifest.json` (hashes and software). structpair 0.3.1.
