# KRAS K1 (proposed; DRAFT; not frozen)

The study design is in [`PROTOCOL_DRAFT.md`](PROTOCOL_DRAFT.md). Nothing was scored, no MD coordinate was read and no
trajectory was downloaded; the endpoint has not been computed.

| File | Content |
|---|---|
| `PROTOCOL_DRAFT.md` | Question, evidence checks, design, ordered decision rule, stop conditions, budget, interpretation limits |
| `protocol_params_DRAFT.json` | Design parameters (single source for the protocol and the schematic); not calibrated |
| `check_k1_provenance.py` | Parses S1 Table of Pantsar et al. 2018, compares it with the Zenodo file names, reads XTC header fields (no coordinates) of the 15 G12D files and the residue labels of the authors' PyInteraph files |
| `results/k1_provenance_checks.json` | Generated record of those checks |
| `plot_design_schematic.py` | Draws the design schematic into `docs/figures/kras_k1/` from the parameter file and the records |

```bash
python3 -I benchmarks/audit/kras_k1/check_k1_provenance.py --fetch-dir <new empty directory outside the repo>   # needs pdftotext
PYTHONPATH=python .venv/bin/python benchmarks/audit/kras_k1/plot_design_schematic.py
```

Sources and terms: Zenodo record 10.5281/zenodo.1346073 (CC BY 4.0); Pantsar et al., PLoS Comput Biol 14:e1006458,
doi:10.1371/journal.pcbi.1006458 (CC BY); crystal entries 5US4 and 5USJ (PDB, CC0). Nothing was sent externally.
