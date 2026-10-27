# Draft feature-export specification for the MET inhibitor DMS (unfrozen; not implemented)

Purpose: the smallest exporter whose output can be used for the exploratory comparison in `SCREEN.md`. No model, no feature selection, no outcome reading.

## Rows and files
- **Unit:** (position p in 1..287, substituted residue, inhibitor). Inhibitors used for the structure-dependent blocks: crizotinib (2WGJ), savolitinib (6SDE), merestinib (4EEV), AMG-458 (5T3Q), tepotinib (4R1V), NVP-BVU972 (3QTI). Tivantinib (MET-independent activity per the article) and the four docked inhibitors (no poses available) are not exported with structure blocks; this rule uses article statements only.
- **Outputs:** `features.csv` (inputs only) with a `features_manifest.json` (units, definitions, input SHA-256, template construction, missingness counts). Outcome tables are produced by a separate script and are never read by the exporter.

## Mapping rules
- DMS position p = UniProt P08581 residue p + 1058 (wild-type identity verified at all 287 positions; the exporter asserts it from the file hash).
- Template residues are matched by PDB residue number (equal to UniProt numbering in the entries checked); a template residue whose name differs from P08581 is flagged (3QTI 1194 and 1234, 3DKC 1194, 1234, 1235) and its structural features set missing.
- Chain: the chain with the most resolved residues in 1059-1345, ties to the first chain; altloc blank or A; hydrogens and waters ignored; unresolved residue -> missing value plus a missingness flag; no imputation in the exporter.

## Feature blocks (fixed in advance, units explicit)
- **M, mutation properties:** wild-type and substituted residue; volume change (Å^3), hydropathy change (Kyte-Doolittle), charge change, BLOSUM62 score; proline and glycine involvement flags.
- **L, template site context (per inhibitor):** minimum heavy-atom distance between the residue and the ligand (Å); count of residue atoms within 4.5 Å of the ligand; neighbour count within 10 Å of the Cα.
- **T, TopoFold local geometry:** curvature, torsion (sin and cos) of the 9-residue Cα window centred on p from the template (existing `compute_ca_features`); window with any unresolved residue -> missing.
- **C, equal-information Cα baseline:** the 36 pairwise Cα distances of the same window (Å), same residues and frames as T.
- **D, drug baselines:** inhibitor identity and class; no outcome statistic is exported (training means are computed inside folds by the analysis, not here).

## Tests (written before any real-data run)
1. Numbering: offset 1058 reproduces wild-type identity on the real file hash; a shifted offset fails.
2. Template identity: the flagged residues above are exactly those found; no other mismatch.
3. Distances: ligand-site and window distances equal an independent computation (scipy cdist) on synthetic coordinates, and one real template spot-check.
4. Missing data: a synthetic window with a gap yields missing values and a flag, never an imputed value.
5. Outcome independence: output bytes are identical with the outcome file absent, shuffled or altered; the exporter module imports no outcome reader.
6. Invariance: T is unchanged under rigid motion; a mirror image changes torsion sign only (existing core tests reused).
7. Determinism: two runs give identical hashes; no column named like an outcome or a published prediction.

Stopping condition: stop and report if numbering or template-identity checks fail beyond the listed engineered residues, or if the owner does not confirm the licence basis.
