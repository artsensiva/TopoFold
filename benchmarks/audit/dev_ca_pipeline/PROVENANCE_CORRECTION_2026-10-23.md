# DEV-CA-1 Provenance Correction (2026-10-23)

**Scope:** one provenance field. No DEV-CA-1 measurement, array, check or original record has been changed. The
pipeline was not rerun and the native extension was not rebuilt for this correction.

## Wrong field

- **Record:** `results/DEV-CA-1_features.json` (SHA-256 `30e458c6fa5e5d345507a78c3ff3e775f2741034884a2b4e33fb479aaeefa956`,
  committed in `6a9a64a`; the identical copy is in the git-ignored `build/`).
- **Field:** `provenance.core.binary_sha256` = `49780093095f27c89bdf234be15e1f5a9bac957d5a9dbf5ae0f93c856ee1590d`.
- **What it identifies:** the Python wrapper `topofold/__init__.py`, not a binary. `features.py` 0.1.0 hashed
  `topofold.__file__`, which in this maturin layout is the 115-byte wrapper (`from .topofold import *`). That value
  equals the SHA-256 of the installed wrapper and the wrapper entry in the installed `RECORD`.
  `compute_invariants` is defined by the native module `topofold.topofold`
  (`topofold.cpython-314-x86_64-linux-gnu.so`).

## Can the native binary used by DEV-CA-1 be established?

**No. Historical native-binary identity: unknown.**

| Evidence | Observation |
|---|---|
| Run time | `created_utc` in the features record: 2026-10-22T22:21:02Z |
| Installed `.so` | mtime 2026-10-22T22:23:05Z (00:23:05 +0200), from the `maturin develop --release` run during later verification; SHA-256 `d90e7ed4…` matches the current `RECORD`. It replaced the build that DEV-CA-1 loaded |
| Cargo outputs | `target/release/libtopofold.so` and `target/release/deps/libtopofold.so` have the same 00:23:05 mtime (overwritten) |
| Wheels | `maturin develop` wheels are written to temporary directories; every remaining `/tmp/.tmp*` directory is empty, and no `topofold*.whl` exists in `/tmp` or `~/.cache`. `target/wheels/…cp314…whl` (2026-10-03) and `dist/topofold-0.1.0…whl` predate the MI fix and the build in question |
| Source | The DEV-CA-1 build was compiled from the working tree after the MI fix. Later, before the rebuild, only Rust doc comments changed. Equal source does not establish an equal binary, so the current binary is **not** substituted as the historical one |

## What remains valid

- The features source hash recorded by DEV-CA-1 (`c88dee13…`) equals `python/topofold_inputs/features.py` at `b721eaf`.
- The run's numerical checks are separate evidence and are unaffected: all pre-specified checks passed, and κ/τ agree
  with an independent NumPy reference to 5.6e-16 rad (`results/DEV-CA-1_summary.json`).

## Fix for future exports

`features.py` 0.1.1 (schema `topofold_inputs.features/2`) replaces `binary_sha256` with `core.wrapper` and
`core.native_extension`, each with its own path and SHA-256. The native module is taken from
`compute_invariants.__module__` and accepted only if Python loaded it with an extension-module loader from a file with
an extension suffix. Otherwise it is recorded as `unresolved`, with a reason and no hash.
