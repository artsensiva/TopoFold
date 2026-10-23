# KRAS-STRUCT-1: GDP-bound (9O0N) vs GMPPNP-bound (9O0O) wild-type KRAS, both with MRTX1133

Generated from `KRAS-STRUCT-1_results.json` and `KRAS-STRUCT-1_table.csv`. **Static structures: direct geometric differences only.** Confounded by crystal form (P 21 21 21 vs P 2 3) and crystallization conditions, resolution (1.40 vs 1.90 Å), alternate conformations and coordinate error; both contain MRTX1133. No ensemble, population, kinetic or functional-state claim; nucleotide identity does not establish functional state.

- Mapping: chain A, residues 1–168 (168 residues), names identical, all present; breaks: GDP none, GMPPNP none.
- Global Cα RMSD after optimal superposition (MDAnalysis QCP, 168 Cα): 0.718 Å (the authors report 0.3 Å; atom selection may differ).
- Alternate locations: primary policy blank-then-A; sensitivity (highest occupancy): ok; max change per metric {"d_kappa": 0.0234, "d_i_i2": 0.0773, "rmsd3": 0.0337, "d_tau": 0.045, "d_i_i3": 0.0975, "rmsd4": 0.0364}

## Largest local differences (descriptive)

| Metric | Unit | Ten largest windows (value; * = window contains an altloc residue) |
|---|---|---|
| d_kappa | rad | 105-107 (0.528*), 28-30 (0.468), 31-33 (0.422), 64-66 (0.395), 30-32 (0.359), 103-105 (0.347), 63-65 (0.321), 36-38 (0.302), 35-37 (0.285), 37-39 (0.216) |
| d_i_i2 | Å | 28-30 (0.959), 105-107 (0.866*), 103-105 (0.741), 36-38 (0.703), 63-65 (0.683), 31-33 (0.656), 64-66 (0.633), 30-32 (0.608), 35-37 (0.529), 34-36 (0.522) |
| rmsd3 | Å | 105-107 (0.553*), 28-30 (0.527), 31-33 (0.435), 64-66 (0.410), 103-105 (0.399), 30-32 (0.380), 63-65 (0.365), 36-38 (0.358), 35-37 (0.311), 34-36 (0.253) |
| d_tau | rad | 36-39 (0.898), 28-31 (0.894), 31-34 (0.869), 107-110 (0.720), 59-62 (0.695), 37-40 (0.654), 32-35 (0.555), 26-29 (0.547), 30-33 (0.524), 34-37 (0.521) |
| d_i_i3 | Å | 102-105 (1.915), 27-30 (0.960), 103-106 (0.938*), 30-33 (0.867), 26-29 (0.821), 63-66 (0.783), 60-63 (0.745), 32-35 (0.703), 29-32 (0.638), 62-65 (0.613) |
| rmsd4 | Å | 28-31 (0.888), 102-105 (0.764), 35-38 (0.712), 107-110 (0.666), 105-108 (0.597*), 30-33 (0.590), 36-39 (0.588), 103-106 (0.574*), 32-35 (0.564), 27-30 (0.535) |

## Agreement between equal-information metrics (Spearman over windows; descriptive)

| Pair | Spearman | Windows |
|---|---|---|
| d_kappa~d_i_i2 | 0.908 | 166 |
| d_kappa~rmsd3 | 0.973 | 166 |
| d_i_i2~rmsd3 | 0.950 | 166 |
| d_tau~d_i_i3 | 0.478 | 165 |
| d_tau~rmsd4 | 0.786 | 165 |
| d_i_i3~rmsd4 | 0.753 | 165 |

Three-point windows: κ, d(i−1, i+1) and the 3-residue RMSD all depend on the same triangle, so close agreement is expected. Four-point windows add torsion; d(i−1, i+2) is blind to the sign of τ. Agreement or disagreement does not show that one representation is better.

Resources: 12.8 CPU-s, 4.6 s wall.
