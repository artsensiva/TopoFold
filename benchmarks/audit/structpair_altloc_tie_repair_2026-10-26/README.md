# structpair optional altloc-sensitivity tie repair (KNOWN_ISSUES m13), 2026-10-26

Engineering repair verified with synthetic cases and a regression replay of the six KIN-PAIR-1 specifications. **Not a new evaluation**; KIN-PAIR-1 keeps its frozen-rule outcome and all its files.
The generated [`REPORT.md`](REPORT.md) holds the before/after table and the BTK-2 precision trace (numbers read from `results/*.json`).

| File | Content |
|---|---|
| `DESIGN.md` | The behaviour, written and committed before the implementation |
| `record_repair.py`, `results/repair_before_after.json` | `tests/test_structpair_altloc_ties.py` on the commit before the repair and on the repaired code |
| `replay_kin_pair_1.py`, `results/replay_summary.json`, `results/replay/` | The six original KIN-PAIR-1 specifications (optional request kept) run with structpair 0.3.1; compact outputs |
| `diagnose_btk2_precision.py`, `results/btk2_precision_diagnosis.json` | BTK-2 invariance deviation traced through input, computation, formatting and comparison, with a SYNTHETIC control |
| `make_report.py`, `REPORT.md` | Report generated from the JSON records |

Reproduce: `.venv/bin/python benchmarks/audit/structpair_altloc_tie_repair_2026-10-26/replay_kin_pair_1.py`, then `diagnose_btk2_precision.py`, `record_repair.py`, `make_report.py` (needs the git-ignored KIN-PAIR-1 `data_pdb/` and `work/` directories, no download).
