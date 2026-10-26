# KIN-PAIR-1: structpair on previously unused human kinases (2026-10-26)

**Structural-analysis and engineering evaluation; not a test of drug efficacy, resistance, causal allostery or clinical benefit.**
Question: can structpair map, analyze and explain structural differences on previously unused human kinases, with transparent exclusions and
agreement with independent numerical references? The results are in the generated [`report.md`](report.md) (all numbers read from `results/*.json`).

| File | Content |
|---|---|
| `PROTOCOL.md`, `protocol.json` | The protocol, frozen before selection and execution; dated amendments 1 and 2 at the end |
| `results/KIN-PAIR-1_freeze_1.json`, `_1a.json`, `_2.json` | Hashes of the protocol, code and inputs at each freeze (1: selection; 1a: after amendment 1; 2: analysis code and panel before execution) |
| `results/candidate_targets_verified.json`, `results/exposure_history.json` | The 23 declared candidate kinases (verified against UniProt) and the exposure history that is excluded |
| `select_panel.py`, `results/KIN-PAIR-1_selection.json` | Deterministic selection and every candidate with its reason; `..._selection_run1_defective.json` and `select_run1.log` keep the first run (amendment 1) |
| `kinpair.py` | Independent PDB parsing, UniProt mapping checks, numpy references |
| `run_study.py`, `results/KIN-PAIR-1_results.json`, `results/pairs/` | Primary execution (frozen code) and compact structpair outputs |
| `supplement_altloc_ties.py`, `results/KIN-PAIR-1_supplement.json`, `results/pairs_supplement/` | Post-hoc supplement (amendment 2) for the 4 pairs refused on altloc occupancy ties |
| `make_report.py`, `report.md`, `figures/` | Report and two figures generated from the JSON files |

Reproduce (network for RCSB/UniProt, about 8 MB): `.venv/bin/python benchmarks/audit/kin_pair_1/select_panel.py run` then `run_study.py run`, `supplement_altloc_ties.py`,
`make_report.py`. Raw PDB files (`data_pdb/`) and work directories are git-ignored; their SHA-256 values are in the results. Tests: `tests/test_kin_pair_1.py` (SYNTHETIC fixtures).
