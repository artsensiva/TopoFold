# Verified structural-role annotation of the existing PIM1 comparisons (2026-10-26)

Interpretation of existing outputs for 8R0Q, 6NO9, 6VRV and 6L17: where do the seven new windows per pair around residue 261 lie in the verified PIM1 sequence and structure?
No score was recomputed, no threshold changed, no structure selected. Rules (anchors, region-label policy) were fixed in `DESIGN.md` before the run.

| File | Content |
|---|---|
| `DESIGN.md` | Anchors, region-label rules and plan, committed before the run |
| `annotate_pim1.py`, `run_roles.py` | The annotation (UniProt record, αC contact rule, existing `kinase_roles` helper via the adapter interpreter, windows by complete residue identifiers) |
| `results/pim1_annotation.json`, `results/kinase_roles_reports.json` | The machine-readable annotation and the helper's full reports |
| `make_report.py`, `REPORT.md`, `figures/pim1_windows_schematic.png` | Generated report and the captioned schematic (verified landmarks and the added windows only) |

Tests: `tests/test_pim1_annotation.py` (consistency of the record). Reproduce: `.venv/bin/python benchmarks/audit/pim1_structural_annotation_2026-10-26/annotate_pim1.py`, then `make_report.py`
(needs the git-ignored KIN-PAIR-1 `data_pdb/`; one small UniProt retrieval, no structure download).
