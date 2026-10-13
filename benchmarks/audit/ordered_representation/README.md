# WP-C ordered angular representation: reproduction

SYNTHETIC evaluation of residue-ordered (cos τ, sin τ) features. Findings: [`REPORT.md`](REPORT.md).
Protocol: [`docs/ORDERED_ANGULAR_REPRESENTATION_PROTOCOL.md`](../../../docs/ORDERED_ANGULAR_REPRESENTATION_PROTOCOL.md).

## Inputs and outputs

| Path | Kept? | Content |
|------|-------|---------|
| `*.py` | yes | `features.py`, `retrieval.py`, `inputs.py`, `controls.py`, `run.py`, `make_report.py`, `test_ordered.py` |
| `fixtures/SYNTHETIC_wpc_controls.json` | yes | N1/N2 parameters, draw order and SHA-256 |
| `results/SYNTHETIC_ordered_representation_results.json` | yes | All checks, retrieval scores, verdicts, input hashes |
| `results/SYNTHETIC_ordered_representation_provenance.json` | yes | HEAD, versions, binding identity, commands and exit codes |
| `REPORT.md` | yes | Rendered from the results |
| `build/` | **no** (git-ignored) | N1/N2 binaries and the comparator outputs (hashes in the results) |

**WP-A inputs are read-only.**
- WP-A fixtures come from `../angular_characterization/build/fixtures/`, checked against the
  frozen WP-A manifest.
- If those files are absent, the cases are regenerated **in memory** with the WP-A generator
  and must reproduce the manifest hashes. Only the comparator's scratch copies are written, and
  they go under this directory's `build/`.
- The WP-A pipeline is never rerun, and nothing is written into the WP-A directory.

## Environment

This uses the WP-A isolated interpreter, which holds the binding built from the current source
(see `../angular_characterization/README.md` for how to create it). `run.py` records the
binding's identity and checks that the production sources are unchanged since that build.
A Rust toolchain ≥ 1.88 is needed for the comparator (see KNOWN_ISSUES m12).

## Commands (from the repository root)

```bash
PY=$PWD/benchmarks/audit/angular_characterization/build/venv/bin/python
$PY -B benchmarks/audit/ordered_representation/run.py            # full run; writes results/
$PY -B benchmarks/audit/ordered_representation/make_report.py    # renders REPORT.md
$PY -B benchmarks/audit/ordered_representation/run.py --verify   # recompute, compare, write nothing
cd benchmarks/audit/ordered_representation && $PY -B -m pytest -p no:cacheprovider test_ordered.py
```

`--verify` uses the comparator outputs in `build/`. It refuses to proceed if their hashes differ
from those recorded in the results. Without `build/`, run `run.py` first; the outputs are
deterministic for the same environment.
