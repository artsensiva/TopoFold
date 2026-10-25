# KRAS K1 (protocol v1 frozen 2026-10-25; crystal endpoint step S1 only)

The study design is in [`PROTOCOL.md`](PROTOCOL.md) (frozen at S0 by `results/KRAS-K1_freeze.json`). No MD trajectory has
been downloaded, opened or scored. The 5US4/5USJ crystal contrast is an external, already-exposed static structural proxy,
not held out and not functional ground truth.

| File | Content |
|---|---|
| `PROTOCOL.md` | Question, evidence checks, design, ordered decision rule C1–C4, stop conditions, budget, interpretation limits |
| `protocol_params.json` | Frozen design parameters (single source for the protocol, the endpoint code and the schematic); not calibrated |
| `k1_endpoint.py` | Frozen step S1: `freeze` (hashes, input identities) and `run` (the crystal-only endpoint, once; refuses a changed file or an existing record) |
| `check_k1_provenance.py`, `results/k1_provenance_checks.json` | S0 evidence: S1 Table of Pantsar et al. 2018 against the Zenodo names, XTC header fields of the 15 G12D files, PyInteraph labels |
| `plot_design_schematic.py` | Draws the design schematic into `docs/figures/kras_k1/` |
| `results/KRAS-K1_freeze.json` | Hashes of the protocol, parameters, endpoint code, tests and inputs, recorded before the endpoint was computed |
| `results/KRAS-K1_endpoint.json`, `results/KRAS-K1_endpoint_report.md`, `figures/` | The endpoint record, its short report and its figure (generated after the freeze commit) |

```bash
.venv/bin/python benchmarks/audit/kras_k1/k1_endpoint.py freeze      # once, before the endpoint
.venv/bin/python benchmarks/audit/kras_k1/k1_endpoint.py run         # once
PYTHONPATH=python .venv/bin/python benchmarks/audit/kras_k1/plot_design_schematic.py
python3 -I benchmarks/audit/kras_k1/check_k1_provenance.py --fetch-dir <new empty directory outside the repo>   # needs pdftotext
```

Sources and terms: Zenodo record 10.5281/zenodo.1346073 (CC BY 4.0); Pantsar et al., PLoS Comput Biol 14:e1006458,
doi:10.1371/journal.pcbi.1006458 (CC BY); crystal entries 5US4 and 5USJ (PDB, CC0). Nothing was sent externally.
