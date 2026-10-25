# Two-structure comparison: Descriptive demonstration: 1OPJ chains A,B vs 1OPK chain A (mouse Abl1, different inhibitor, construct, engineered D363N, crystal form)

**Descriptive, static.** Two static structures: the values are direct geometric differences between these coordinate sets. They can reflect crystal form and packing, crystallization conditions, resolution, refinement, alternate conformations and coordinate error as well as any biological difference. No ensemble, population, kinetic, functional-state or treatment-response claim follows. Agreement or disagreement between metrics does not show that one representation is better. Crystal copies of one entry are not independent replicates.

## Inputs

| | Label | File (SHA-256, first 12) | Chain | Altloc policy |
|---|---|---|---|---|
| A | 1OPJ (imatinib; kinase domain) | `1OPJ.pdb` (6d3444a1f081) | A (primary of A, B) | reject |
| B | 1OPK (PD166326; SH3-SH2-kinase) | `1OPK.pdb` (867ed2b25b80) | A | reject |

## Within-entry crystal-copy variation

**Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of them. The within-entry crystal-copy variation is a descriptive range of differences among the chains the user selected. It is not a calibrated noise floor, not a statistical uncertainty, not a nucleotide or ligand effect and not an HDX replicate estimate. An advisory flag appears only for a channel whose minimum difference the user configured; it is an advisory rule, not a gate, not an adequacy verdict, not significance and not biological separation.**

- Status: **partial**. Chains A ['A', 'B'] (primary A), B ['A'] (primary A). Pairs: 2 cross, 1 within A, 0 within B. The legacy tables above are the primary pair.
- The numeric between-entry and within-entry values are always given in `chain_variation.csv`; an empty cell is unavailable, never zero, and the `*_n` columns count the chain pairs that contributed.
- **No advisory rule is configured, so no flag is assessed.** Set `chain_variation.advisory.min_difference` to a finite positive number for a channel (in its unit) to obtain a flag for that channel.
- Displacement trim fraction 0.2.

| Channel | Unit | Positions with a value | Advisory rule | Rule met | Rule not met | Not assessed |
|---|---|---|---|---|---|---|
| disp | Å | 282 | not configured | 0 | 0 | 282 |
| d_kappa | rad | 280 | not configured | 0 | 0 | 280 |
| d_i_i2 | Å | 280 | not configured | 0 | 0 | 280 |
| rmsd3 | Å | 280 | not configured | 0 | 0 | 280 |
| d_tau | rad | 279 | not configured | 0 | 0 | 279 |
| d_i_i3 | Å | 279 | not configured | 0 | 0 | 279 |
| rmsd4 | Å | 279 | not configured | 0 | 0 | 279 |

Largest between-entry displacement (median over cross pairings) with the within-entry maximum (top 8):

| A residue | Between median [min–max] (Å) | Within-entry max (Å) | Scope | Advisory rule |
|---|---|---|---|---|
| A:409 | 33.054 [32.853–33.255] | 1.084 | A_only | not assessed |
| A:408 | 31.720 [31.664–31.776] | 1.387 | A_only | not assessed |
| A:410 | 30.137 [29.987–30.287] | 0.646 | A_only | not assessed |
| A:407 | 29.042 [28.894–29.190] | 0.462 | A_only | not assessed |
| A:406 | 24.619 [24.538–24.699] | 0.247 | A_only | not assessed |
| A:411 | 24.128 [24.102–24.154] | 0.426 | A_only | not assessed |
| A:405 | 20.997 [20.964–21.029] | 0.353 | A_only | not assessed |
| A:412 | 19.351 [19.274–19.428] | 0.277 | A_only | not assessed |

Files: `chain_variation.csv` (long format; every statistic with its count), `chains.csv`, `chain_variation.png`.


## Mapping

- Correspondence: same_numbering (282 positions; 282 present in both). Missing-residue policy: `reject`.
- Residue-name mismatches (flag): A:382 ASP→A:382 ASN.
- Alternate locations resolved: none.
- Breaks (adapter): A none, B none (row index after which the chain breaks).

## Global Cα RMSD (context only)

- 5.121 Å over 282 Cα after optimal superposition (MDAnalysis QCP); independent numpy SVD/Kabsch check 5.121 Å (|difference| 2.7e-14 Å).
- Published RMSD values often use other atom sets, residue ranges or outlier rejection; compare only with a stated definition.

## Largest local differences (top 8; descriptive)

| Metric | Unit | Windows (A numbering; value; * altloc, † residue-name mismatch in window) |
|---|---|---|
| d_kappa | rad | A:410..A:412 (0.693), A:249..A:251 (0.645), A:403..A:405 (0.613), A:405..A:407 (0.590), A:417..A:419 (0.535), A:398..A:400 (0.510), A:408..A:410 (0.509), A:409..A:411 (0.478) |
| d_i_i2 | Å | A:410..A:412 (1.467), A:249..A:251 (1.339), A:403..A:405 (1.271), A:405..A:407 (1.244), A:408..A:410 (1.087), A:409..A:411 (0.999), A:398..A:400 (0.968), A:417..A:419 (0.967) |
| rmsd3 | Å | A:410..A:412 (0.789), A:249..A:251 (0.728), A:403..A:405 (0.693), A:405..A:407 (0.671), A:408..A:410 (0.581), A:417..A:419 (0.577), A:398..A:400 (0.558), A:409..A:411 (0.542) |
| d_tau | rad | A:410..A:413 (3.030), A:296..A:299 (2.917), A:405..A:408 (2.718), A:250..A:253 (2.683), A:400..A:403 (2.577), A:408..A:411 (2.555), A:409..A:412 (2.465), A:403..A:406 (2.388) |
| d_i_i3 | Å | A:403..A:406 (4.555), A:410..A:413 (4.001), A:404..A:407 (3.736), A:250..A:253 (3.432), A:407..A:410 (3.032), A:248..A:251 (2.945), A:252..A:255 (2.833), A:409..A:412 (2.649) |
| rmsd4 | Å | A:410..A:413 (2.298), A:250..A:253 (2.187), A:403..A:406 (2.166), A:409..A:412 (1.985), A:296..A:299 (1.860), A:404..A:407 (1.852), A:402..A:405 (1.810), A:297..A:300 (1.792) |

## Agreement between equal-information metrics (Spearman over windows; descriptive)

| Pair | Spearman | Windows |
|---|---|---|
| d_kappa~d_i_i2 | 0.937 | 280 |
| d_kappa~rmsd3 | 0.981 | 280 |
| d_i_i2~rmsd3 | 0.973 | 280 |
| d_tau~d_i_i3 | 0.527 | 279 |
| d_tau~rmsd4 | 0.804 | 279 |
| d_i_i3~rmsd4 | 0.757 | 279 |

Three-residue metrics share one triangle, so close agreement is expected; d(i−1, i+2) is blind to the sign of τ.

## Exclusions

No window inside the mapped range was excluded.

## Alternate-conformation sensitivity

Not requested.

## Figure

![structure pair](figure.png)

Files: `windows.csv` (all metrics with units in column names), `residues.csv`, `exclusions.csv`, `summary.json`, `manifest.json` (hashes and software). structpair 0.3.0.
