# structpair numeric record (KNOWN_ISSUES m13b), 2026-10-26

Opt-in `numeric_record.json` (spec `"numeric_record": true`, structpair 0.4.0): the computed float64 window values before the six-significant-figure display rounding of the
CSV files, with stable position and residue identifiers, units and the status of every missing value; listed with its hash in the manifest and checked by `structpair_verify`.
Existing outputs were inspected first: `summary.json` already keeps full precision but only for aggregates, and `work/` intermediates are not hash-verified, so window-level values existed only in the rounded CSV.

| File | Content |
|---|---|
| `btk2_regression.py`, `results/btk2_numeric_record_regression.json` | Separate regression check on the original and transformed BTK-2 inputs (read-only use of the KIN-PAIR-1 supplement specs) |
| `make_report.py`, `REPORT.md` | Generated table: computed value → rounded CSV → precision-preserving record |

Tests: `tests/test_structpair_numeric_record.py` (SYNTHETIC). Reproduce: `.venv/bin/python benchmarks/audit/structpair_numeric_record_2026-10-26/btk2_regression.py` then `make_report.py`
(needs the git-ignored KIN-PAIR-1 `data_pdb/` and `work_supplement/` directories; no download).
