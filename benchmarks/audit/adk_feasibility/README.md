# WP-B G0–G1: AdK real-MD data acquisition and feasibility

Third-party MD data: Orädd et al., *Sci. Adv.* 2021, Zenodo
[10.5281/zenodo.5583119](https://doi.org/10.5281/zenodo.5583119), CC BY 4.0. Reference structures
4AKE and 1AKE are from the RCSB PDB. Not SYNTHETIC. No TopoFold or candidate-feature scoring takes
place here.

- Specification, frozen before trajectory inspection: [`SPEC.md`](SPEC.md) and `wpb_spec.json`,
  with hashes in `results/adk_spec_freeze.json`.
- Amendments 1–2 ([`AMENDMENTS.md`](AMENDMENTS.md), `wpb_spec_amendment{1,2}.json`) were recorded
  before G0/G1 ran, with freeze records `results/adk_spec_amendment{1,2}_freeze.json`. Where the
  amendments differ from `SPEC.md`, they govern.
- Findings: [`REPORT.md`](REPORT.md).

## Reproduce (from this directory)

```bash
# isolated environment (CPython 3.14.4)
python3 -m venv build/venv
build/venv/bin/pip install "numpy==2.5.3" "MDAnalysis==2.10.0" "pytest==9.1.1"
build/venv/bin/pip freeze > build/pip_freeze.txt

# authorized downloads (~1.37 GB + two PDB files); raw files are never modified
mkdir -p data/raw && cd data/raw
for f in ORADD_open_state.gro ORADD_closed_state.gro ORADD_open_state_traj.xtc ORADD_closed_state_traj.xtc; do
  curl -sS -L --fail --retry 3 -o "$f" "https://zenodo.org/records/5583119/files/$f?download=1"; done
for id in 4AKE 1AKE; do curl -sS -L --fail --retry 3 -o "$id.pdb" "https://files.rcsb.org/download/$id.pdb"; done
cd ../..

build/venv/bin/python -B -m pytest -p no:cacheprovider test_g1.py   # alignment, PBC, core, classification, precedence
build/venv/bin/python -B g0.py      # integrity, enforced mapping, frames, timestamps, cadence, boxes
build/venv/bin/python -B g1.py      # binds to G0; converged cores; per-frame RMSDs; variants; outcome
build/venv/bin/python -B make_report.py
```

RCSB may remediate PDB files over time. The SHA-256 of the files used is recorded in
`results/adk_g0_provenance.json`; a later download may differ and must be compared with it.

## Outputs

| File | Kept? | Content |
|------|-------|---------|
| `results/adk_spec_freeze.json` | yes | Spec hashes and freeze time |
| `results/adk_g0_provenance.json` | yes | Checksums, environment, topology/trajectory facts, stops and discrepancies |
| `results/adk_g1_summary.json` | yes | Core, d_ref, observable summaries, classification, sensitivity, outcome |
| `results/adk_g1_frames_{open,closed}.csv` | yes | Per-frame continuous observables |
| `data/raw/`, `data/derived/` | no (git-ignored) | Raw downloads; whole Cα arrays (hashes in the summary) |
| `build/` | no (git-ignored) | Isolated environment |
