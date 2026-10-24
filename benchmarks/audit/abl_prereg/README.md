# ABL preregistration preparation (WP1, 2026-10-24)

Supports `docs/ABL_WT_M472I_PREREGISTRATION_DRAFT.md`. **Nothing here scores ABL data.** No trajectory or uptake value
has been obtained.

| File | Purpose |
|---|---|
| `schemas/trajectory_inventory.schema.json` | JSON Schema for the trajectory inventory requested from D. E. Shaw Research: systems, runs, files, hashes, permissions |
| `schemas/hdx_uptake.schema.json` | JSON Schema for one row of the processed HDX table requested from the HDX authors |
| `validate_inputs.py` | Validates both against the schemas, plus table-level rules. It reports comparable WT/M472I pairs and blocking permissions |
| `map_peptides.py` → `hdx_peptide_mapping.json` | Maps the PRIDE PXD034008 peptide list (construct numbering) to c-Abl 1a by exact sequence matching |
| `reference/UniProt_P00519.fasta` | Reference sequence (UniProt P00519 canonical, CC BY 4.0), SHA-256 `f28f5df8…` |

```bash
python3 benchmarks/audit/abl_prereg/map_peptides.py           # needs the git-ignored peptide map (abl_access/)
.venv/bin/python benchmarks/audit/abl_prereg/validate_inputs.py --inventory INV.json --hdx UPTAKE.csv
.venv/bin/python -m pytest tests/test_abl_prereg_inputs.py     # SYNTHETIC documents only
```

**Mapping result (sequence-based):**
- c-Abl = construct + 228; 76 of 80 peptides match P00519 exactly and uniquely.
- The 4 that do not are the A344V (2) and M472I (2) mutant peptides, differing exactly at 344 and 472.
- Coverage: c-Abl 238–512, 252 of 275 positions.
