# DEV-CA-1: Cα pipeline development exercise (AdK open-start, local data)

Engineering development evidence for `python/topofold_inputs` (adapter → κ/τ features → reports) on one already-local,
already-authorized third-party trajectory (Zenodo 10.5281/zenodo.5583119, CC BY 4.0; WP-B G0 record). It shows that the
pipeline runs on real MD input with consistent mapping and complete, finite features. It makes **no** biological,
state, novelty or superiority claim and is not a comparison with other methods.

| File | Content |
|---|---|
| `run_spec.json` | Specification written before execution (dataset hashes, frame subset, policies, expected checks, cap) |
| `run_dev.py` | Runner: verifies source hashes, runs the pipeline under GNU time with a CPU cap, checks expectations, cross-checks κ/τ against an independent NumPy reference |
| `results/DEV-CA-1_summary.json` | Checks, resources, reference agreement, output hashes |
| `results/DEV-CA-1.json`, `results/DEV-CA-1_features.json`, `results/DEV-CA-1_report.md` | Adapter mapping/provenance report, feature counts/provenance, generated inspection report |
| `PROVENANCE_CORRECTION_2026-10-23.md` | Correction: the recorded `core.binary_sha256` identifies the Python wrapper; the historical native binary is unknown |
| `build/` (git-ignored) | Coordinate and feature arrays; regenerate with the command below |

Reproduce (needs the local AdK files listed in `run_spec.json` and the WP-D audit venv for the adapter stage):

```bash
.venv/bin/python benchmarks/audit/dev_ca_pipeline/run_dev.py
```
