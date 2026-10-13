# WP-A angular characterization: reproduction

SYNTHETIC characterization of the scalar circular-mean/BC torsion scan. Findings:
[`REPORT.md`](REPORT.md). Protocol:
[`docs/ANGULAR_REPRESENTATION_CHARACTERIZATION_PROTOCOL.md`](../../../docs/ANGULAR_REPRESENTATION_CHARACTERIZATION_PROTOCOL.md).

## What is kept and what is regenerated

| Path | Kept? | Content |
|------|-------|---------|
| `*.py` | yes | Harness: generator, trace builder, NumPy reference, P1/REF path, comparison, report, orchestrator, tests |
| `../../../crates/topofold-core/examples/angular_characterization.rs` | yes | Rust paths R1/R2 (a Cargo example, not production code) |
| `fixtures/SYNTHETIC_fixtures.json` | yes | Every case's parameters, seed, RNG draw order and SHA-256 of its binary inputs |
| `results/SYNTHETIC_angular_characterization_{rust,python,summary,provenance}.json` | yes | Original run outputs and provenance |
| `results/SYNTHETIC_angular_characterization_editorial_log.json` | yes | Post-run editorial changes to the report (separate from run provenance) |
| `REPORT.md` | yes | Rendered from the results JSON by `make_report.py` |
| `build/` | **no** (git-ignored) | Fixture binaries, per-frame arrays, wheel, isolated venv. Regenerable, ≈ 140 MB |

**The saved JSON results reference binary arrays under `build/`.** The per-frame means,
resultant lengths and deviations sit in `build/arrays/*.f64`, stored as `[offset, length]`
references. Summary values (BC, moments, recovery errors, check outcomes) are in the JSON
itself. Any verification that needs the per-frame arrays, e.g. rerunning `compare.py`, requires
one of:
- regenerating them with the commands below. Fixture regeneration is deterministic and is
  checked against the hashes in `fixtures/SYNTHETIC_fixtures.json`; the arrays file hash
  `arrays_sha256` is in the Python results; or
- a separately retained evidence bundle containing `build/fixtures/` and `build/arrays/`, with its
  hashes checked against the JSON.

The original `build/` from the recorded run is not version-controlled.

## Environment (as recorded in the provenance file)

- Rust toolchain: rustc/cargo 1.98.1. Rust ≥ 1.88 is needed regardless of the declared
  `rust-version`; see KNOWN_ISSUES m12.
- maturin 1.15.0. The recorded run used the repository's `.venv/bin/maturin`; installing
  `maturin==1.15.0` into any environment works the same way, provided `run_wpa.py`'s `MATURIN`
  path points to it.
- Isolated Python environment: CPython 3.14.4, NumPy 2.5.3, pytest 9.1.1.

Create the isolated environment (from the repository root):

```bash
python3 -m venv benchmarks/audit/angular_characterization/build/venv
benchmarks/audit/angular_characterization/build/venv/bin/pip install "numpy==2.5.3" "pytest==9.1.1"
```

## Reproduce

```bash
# full pipeline: fixtures -> Rust R1/R2 -> wheel build + install into build/venv ->
# P1/REF -> compare -> report; writes results/*.json and REPORT.md, records provenance
benchmarks/audit/angular_characterization/build/venv/bin/python -B \
    benchmarks/audit/angular_characterization/run_wpa.py

# harness tests (isolated copies for all corruption tests)
cd benchmarks/audit/angular_characterization
build/venv/bin/python -B -m pytest -p no:cacheprovider test_harness.py
```

`run_wpa.py` exits nonzero if any command fails. It also fails if the extension loaded by
Python is not byte-identical to the `.so` inside the wheel it has just built.

**A full rerun overwrites the results files.** To preserve the recorded run, copy `results/`
first, or rerun in a separate checkout.

To verify the recorded results without overwriting them:
1. Run only `generate_fixtures.py` and confirm the regenerated binaries match the manifest
   hashes. This rewrites the manifest with identical content; compare its hash with the
   provenance file.
2. Run the Rust example with temporary output paths. It takes the output paths as arguments.
3. `run_python_path.py` always writes to `results/`. Run it in a separate checkout, or copy
   `results/` first.
4. Run `compare.py --rust … --python … --out …` with temporary outputs, then compare the result
   with the recorded summary.
