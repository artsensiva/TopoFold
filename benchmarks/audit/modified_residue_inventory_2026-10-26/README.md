# Modified-residue coverage inventory of the KIN-PAIR-1 inputs (2026-10-26)

Bounded, read-only inventory of non-standard amino acids in the 12 original KIN-PAIR-1 input files (hash-verified against the selection record), the
installed adapter's behaviour on them (the real adapter command run as a diagnostic) and the 3- and 4-residue windows involved. No production code was
changed, nothing was downloaded beyond re-reading the three UniProt reference sequences (hash-checked against the selection record), and no comparison score was computed.
Frequencies apply to this local, selected collection only. The generated [`REPORT.md`](REPORT.md) holds the captioned coverage tables and the recommendation.

| File | Content |
|---|---|
| `inventory.py`, `results/modified_residue_inventory.json` | Per-residue records (identity, MODRES, coordinates, neighbours, adapter and study inclusion, category), gap composition, window counts, denominators, recommendation |
| `make_report.py`, `REPORT.md` | Report generated from the JSON |

Tests: `tests/test_modified_residue_inventory.py` (SYNTHETIC; includes the real adapter on SEP, TPO, PTR and MSE). Reproduce: `.venv/bin/python benchmarks/audit/modified_residue_inventory_2026-10-26/inventory.py`
then `make_report.py` (needs the git-ignored KIN-PAIR-1 `data_pdb/`).
