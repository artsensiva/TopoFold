# KRAS-STRUCT-1: GDP- vs GMPPNP-bound wild-type KRAS crystal structures (descriptive, static)

A small experimental-structure comparison. Pair selection, mapping, metrics and permitted interpretation were frozen
in `run_spec.json` (`results/KRAS-STRUCT-1_freeze.json`) before any comparison score was computed.
- **Inputs:** PDB 9O0N (GDP) and 9O0O (GMPPNP). Both are wild-type KRAS(1–169), same construct, both with the
  inhibitor MRTX1133. Bonsor et al., *Nat Commun* 17:1614 (2026).
- **Licence:** PDB data CC0.
- **Local copies:** the structures are git-ignored in `data_pdb/` (1.6 MB mmCIF + 0.8 MB PDB downloaded; SHA-256 in
  `run_spec.json`).

```bash
.venv/bin/python benchmarks/audit/kras_struct_cmp/run_kras.py run
```

## Selection record

- **Screen:** RCSB, KRAS (P01116), X-ray, GDP- or GNP-bound: 431 entries. Matched-sequence GNP/GDP pairs from one
  citation: 49 combinations.
- **A (9O0N/9O0O): selected.** Matched construct, wild type, one chain each.
  - Its first filter used a summary field that omitted MRTX1133. Full ligand lists then showed MRTX1133 in both
    structures, so the question was narrowed to "both with MRTX1133" before the freeze.
- **B (9U5T/9USB): not usable here.** It is mmCIF-only, and no mmCIF reader is available in the existing
  environments.
- **5US4/5USJ: identified but not examined.** This G12D pair is the only matched pair without a pharmacological
  ligand. It was left unexamined under the two-pair limit.

## Results (from `results/KRAS-STRUCT-1_results.json` and `results/KRAS-STRUCT-1_table.csv`)

**Mapping.**
- Chain A, residues 1–168, identical residue names, no breaks.
- Alternate locations: 7 residues in 9O0N and 1 in 9O0O; primary policy blank-then-A.
- The highest-occupancy sensitivity changes each metric by at most 0.023–0.097 (radians or Å).

**Global Cα RMSD.**
- 0.718 Å after optimal superposition (MDAnalysis QCP, 168 Cα). The authors report 0.3 Å.
- The difference is not resolved here; atom selection or the superposition subset may differ.

**Where local differences are largest (descriptive).** All six metrics place their largest values in residues
26–40, 59–66 and 102–110. Examples:
- κ: 105–107 (0.528 rad), 28–30 (0.468 rad);
- 3-residue RMSD: 105–107 (0.553 Å), 28–30 (0.527 Å);
- τ: 36–39 (0.898 rad), 28–31 (0.894 rad);
- d(i−1, i+2): 102–105 (1.915 Å);
- 4-residue RMSD: 28–31 (0.888 Å).

Windows touching the 9O0N alternate location at residue 106 are flagged.

**Equal-information agreement** (Spearman over windows):
- κ ~ d(i−1, i+1): 0.908;
- κ ~ 3-residue RMSD: 0.973;
- τ ~ d(i−1, i+2): 0.478;
- τ ~ 4-residue RMSD: 0.786.

The three-point metrics are close to redundant, which is expected because they share one triangle. The four-point
metrics agree less; d(i−1, i+2) is blind to the sign of τ. Neither pattern shows that one representation is better.

## Permitted interpretation

- These are direct geometric differences between two crystal structures. They are confounded by:
  - crystal form (P 2₁2₁2₁ vs P 2 3) and crystallization conditions;
  - resolution (1.40 vs 1.90 Å);
  - alternate conformations and coordinate error.
- Both structures contain MRTX1133, so this is not apo KRAS.
- No ensemble, population, kinetic or functional-state claim follows; nucleotide identity does not establish
  functional state.
- The concentration in residues 26–40 and 59–66 agrees in location with the authors' qualitative statement of minor
  Switch-I/II differences. That is a description, not a validation.
- Unused here: Rascore labels, MD, pocket detection, and inference.
