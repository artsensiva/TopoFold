"""Writes README.md from results/manifest.json (no number typed here). Run with the repository .venv."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
m = json.loads((HERE / "results/manifest.json").read_text())
cov = m["coverage"]
t = ["| drug | PDB | chain | ligand (component, residue) | engineered residues | L ok | window ok (T and C) | altloc residues | incomplete side chains |", "|---|---|---|---|---|---|---|---|---|"]
for code, e in m["entries"].items():
    t.append(f"| {code} | {e['pdb']} | {e['chain']} | {e['ligand']['component']} {e['ligand']['resseq']} | {', '.join(f'{q} {n}' for q, n in e['engineered_detail'].items()) or 'none'} | {cov[code]['L']['ok']}/{cov[code]['L']['positions']} | {cov[code]['T']['ok']}/{cov[code]['T']['positions']} | {len(e['altloc_residues'])} | {len(e['incomplete_side_chain_residues'])} |")
ru = m["row_universe"]
md = f"""# MET wild-type template feature export (frozen exporter, generated from `results/manifest.json`)

**What it is:** a deterministic export of structural-site-context features for the MET kinase domain, from six experimental inhibitor complexes and the pinned UniProt P08581 sequence. **What it is not:** the features describe the wild-type template; they do not measure mutation-induced structural change. Masking an engineered template residue does not remove possible wider effects of that residue on the template. No assay table was read, no model trained.

- Contract: `CONTRACT.json`; exporter `met_export.py`; freeze `freeze_exporter.json` (committed before the first export; `run` refuses if frozen files differ). Row universe: {ru['positions']} positions × {ru['substitutions_per_position']} substitutions × {ru['drugs']} inhibitors = {ru['rows']} rows, generated from the reference sequence alone.
- Tables: `features_structure.csv` ({m['tables']['features_structure.csv']['rows']} position×inhibitor rows; identical across the 19 substitutions by construction), `features_mutation.csv` ({m['tables']['features_mutation.csv']['rows']} wild-type/substituted pairs), `coverage.csv`, `manifest.json`; `row_universe.csv` is regenerated and hash-checked (git-ignored). Join: `met_export.join_flat`.
- Blocks: mutation properties (volume, Kyte-Doolittle hydropathy, charge, BLOSUM62; pinned and cited), ligand/site context (4 columns), seven κ and six τ (sin, cos) from each nine-Cα window, all 36 pairwise Cα distances from exactly the same nine coordinates; explicit reason for every missing value.

{chr(10).join(t)}

![coverage](figures/met_feature_coverage.png)

*Figure 1. Purpose: show where each feature block is available, per inhibitor template, along the scanned span. Shows the primary reason per position (blue available, light blue glycine without side chain, orange missing residue, dark red engineered template residue, purple incomplete side chain, black missing Cα). Limitation: window masking propagates to every position whose nine-residue window contains an affected residue; unresolved loops differ by template, so missingness is not random. Source: `results/features_structure.csv`.*

Terms: PDB entries CC0 (wwPDB policy); UniProt sequence attributed (CC BY 4.0 per UniProt; its licence page was not re-read); literature scales cited; BLOSUM62 from NCBI. The exporter does not consume the assay tables, whose terms stay unresolved; no commercial data rights are inferred from any approval.

Reproduce: `.venv/bin/python benchmarks/audit/met_feature_export_2026-10-27/met_export.py run` (about 3 CPU-seconds), `make_coverage_figure.py`, `make_plan.py`, `make_readme.py`. The next evaluation is specified in `EVALUATION_PLAN.md`.
"""
(HERE / "README.md").write_text(md)
print("readme written")
