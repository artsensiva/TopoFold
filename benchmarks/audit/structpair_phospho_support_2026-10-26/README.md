# Phosphorylated-residue support (SEP, TPO, PTR) in the PDB-to-structpair workflow (KNOWN_ISSUES m14), 2026-10-26

Opt-in policy `modified_residues: "phospho_pdb"` (structpair 0.5.0, adapter 0.2.0): PDB-format files only; SEP, TPO, PTR with parents SER, THR, TYR; deposited identity, parent, source and annotation preserved; excluded cases recorded.
**Engineering regression, not a new evaluation:** KIN-PAIR-1 and its frozen outcome are unchanged; the PIM1 pairs are replayed with new derived specifications and the originals are preserved.

| File | Content |
|---|---|
| `DESIGN.md` | The behaviour, written and committed before the implementation |
| `results/component_parents_verified.json` | wwPDB chemical-component parents of SEP, TPO, PTR (metadata) |
| `replay_pim1.py`, `results/phospho_replay.json` | The replay: derived specs (`results/derived_specs/`), expected membership written before the comparison (`results/expected/`), checks, compact outputs (`results/runs/`) |
| `make_report.py`, `REPORT.md` | Report with the captioned before/after residue-and-window coverage tables, generated from the JSON |

Tests: `tests/test_phospho_residues.py` (SYNTHETIC; correct SEP, TPO, PTR; missing Cα and Cα-only; gaps; missing and contradictory annotation; insertion codes; altlocs; a free component and calcium; non-PDB format; SEP versus SER; unchanged defaults; independent references and rigid motion; viewer).
Reproduce: `.venv/bin/python benchmarks/audit/structpair_phospho_support_2026-10-26/replay_pim1.py` then `make_report.py` (needs the git-ignored KIN-PAIR-1 `data_pdb/` and `work/`; no download).
