# WP-A Angular Characterization Report (SYNTHETIC)

Characterization of **one component**, the scalar circular-mean → centered-deviation → Sarle BC torsion scan in `compute_bimodality_profile`, on SYNTHETIC descriptor-level fixtures. It is not validation of the engine or of any biological application, and it changes no production behaviour. Protocol: `docs/ANGULAR_REPRESENTATION_CHARACTERIZATION_PROTOCOL.md` (amendments §9). Every number below is read from `results/SYNTHETIC_angular_characterization_*.json` by `make_report.py`.

## Provenance

- Git HEAD `295d739a783cc1d178131c4c61293042a0c56f7d`. The working tree had 13 uncommitted entries at run time (documentation drafts and this harness); they are listed in the provenance file. Production sources are identified by hash, not by HEAD alone.
- Toolchain: rustc 1.98.1 (48a229cea 2026-09-01); cargo 1.98.1 (797e8a9bc 2026-08-05); maturin 1.15.0; isolated venv Python 3.14.4, NumPy 2.5.3, pytest 9.1.1.
- Python extension: wheel `topofold-0.8.3-cp314-cp314-manylinux_2_34_x86_64.whl` built by this run from the current source. Loaded extension `/home/artem/TopoFold/benchmarks/audit/angular_characterization/build/venv/lib/python3.14/site-packages/topofold/topofold.cpython-314-x86_64-linux-gnu.so` has SHA-256 `732909cc99ec88835f0b8a1f60152b7f134ae4a353d152be7bd9f048e232ba01`; equal to the `.so` inside the freshly built wheel: **True**. Loaded from the isolated venv: **True**.
- Fixtures: 422 cases (105 with coordinates); manifest SHA-256 `17d5bd2827688e917eff84f7fed9404929214d6549a7fadec03ad20fab08c928`. Seed rule: default_rng(20261004 + k), k = base-fixture index; rotations reuse base draws; F6 uses default_rng(1) replay.
- Run status: `completed`. Commands and exit codes:

| Exit | Command |
|---:|---|
| 0 | `git rev-parse HEAD` |
| 0 | `git status --porcelain` |
| 0 | `rustc --version` |
| 0 | `cargo --version` |
| 0 | `.venv/bin/maturin --version` |
| 0 | `benchmarks/audit/angular_characterization/build/venv/bin/python -c import sys;print(sys.version)` |
| 0 | `benchmarks/audit/angular_characterization/build/venv/bin/python -c import numpy;print(numpy.__version__)` |
| 0 | `benchmarks/audit/angular_characterization/build/venv/bin/python -c import pytest;print(pytest.__version__)` |
| 0 | `benchmarks/audit/angular_characterization/build/venv/bin/python -B benchmarks/audit/angular_characterization/generate_fixtures.py` |
| 0 | `cargo run --release -q -p topofold-core --example angular_characterization -- benchmarks/audit/angular_characterization/build/fixtures/SYNTHETIC_fixtu` |
| 0 | `.venv/bin/maturin build --release -i benchmarks/audit/angular_characterization/build/venv/bin/python --out benchmarks/audit/angular_characterization/b` |
| 0 | `benchmarks/audit/angular_characterization/build/venv/bin/python -m pip install -q --force-reinstall --no-deps benchmarks/audit/angular_characterizatio` |
| 0 | `benchmarks/audit/angular_characterization/build/venv/bin/python -B -c import topofold, topofold.topofold as e, hashlib, pathlib;p=pathlib.Path(e.__fil` |
| 0 | `benchmarks/audit/angular_characterization/build/venv/bin/python -B benchmarks/audit/angular_characterization/run_python_path.py` |
| 0 | `benchmarks/audit/angular_characterization/build/venv/bin/python -B benchmarks/audit/angular_characterization/compare.py` |
| 0 | `benchmarks/audit/angular_characterization/build/venv/bin/python -B benchmarks/audit/angular_characterization/make_report.py` |

Source hashes:

| File | SHA-256 |
|---|---|
| `crates/topofold-core/src/circular.rs` | `94244f9dd563470784b5abe3d2a65211de0b4b261f377519a23ca60fdd440ae3` |
| `crates/topofold-core/src/bimodality.rs` | `0c093efd9cf607093e03920e5e10cc5a5c23c5494867c8a325b0ccc5edc6aad0` |
| `crates/topofold-core/src/geometry.rs` | `1076afcb996893919e14027e354ab8bf922185fb9e489f04cb581ad4525ea1ba` |
| `crates/topofold-core/src/types.rs` | `ae9c9ea98bb7246374c43518c7e605a36118b0e210ac4d6088fda7cd8a74839f` |
| `crates/topofold-core/src/lib.rs` | `6024f14ff67b8e17941a6f671899baff89fa9b0c6bb25d5c2dfaf6f8b3dd5ca2` |
| `crates/topofold-python/src/lib.rs` | `4d102f1d60126a13d5961ac2592b73566b029594409e1eacf49471ab083d3757` |
| `Cargo.lock` | `ed9396593b5d5e19eb9fa9440b84569ea738665e44b586c377fca51384a44a61` |
| `crates/topofold-core/examples/angular_characterization.rs` | `aae324beb91c74631e6f39e59bfcdbb22ef55545493d508486abec43cf47825f` |
| `docs/ANGULAR_REPRESENTATION_CHARACTERIZATION_PROTOCOL.md` | `0f0ff86d4ec68f796d7af2d1ca88cb3cad811a4ef179eacb59102a7f598712ae` |

## 1. Implementation correctness (harness)

| Check | Result |
|---|---|
| ref arrays hash matches python results | PASS |
| coverage: rust results cover all fixtures | PASS |
| coverage: python results cover all fixtures | PASS |
| coverage: R2 present exactly for coordinate fixtures | PASS |
| coverage: P1 present exactly for coordinate fixtures | PASS |
| fixture files match manifest hashes | PASS |
| R2 construction recovery <= 1e-09 rad | PASS |
| R2 replication bitwise equal to compute_bimodality_profile | PASS |
| R1 sarles_bimodality_coefficient equals StreamingMoments BC bitwise | PASS |
| R1/R2 vs REF agree wherever conditioned (protocol §4) | PASS |
| analytical controls (F1, F3, F4 noiseless) reproduced by REF and R1 | PASS |
| F6 replay consistent with values documented in KNOWN_ISSUES M3 | PASS |
| F6b coordinate path equals scalar path within conditioned tolerance (no wraps) | PASS |

- Worst construction recovery: R2 (f64) 4.441e-16 rad; P1 (float32) 3.769e-07 rad; P1 cases above the 1e-5 flag: 0.
- Comparisons with REF: 449 conditioned (compared, all within the pre-declared bound); 69 not conditioned (recorded only, because a mean was numerically undefined, the deviation set was too narrow, or the guard was ambiguous).
- The R2 bitwise check applies only to the exact replication of production operations on identical f64 inputs. REF and float32 inputs are never compared bitwise.
- Coverage: R1 and REF ran on all 422 cases; R2 and P1 ran only on the 105 cases with coordinates. In particular, F4 rotations other than α = 0 and the F6a scalar inputs have R1/REF results only.

## 2. Representation information loss

**F1: coordinate-wise states with equal circular mean.**

| Case | BC R1 | BC R2 | BC P1 | guard (R2) | scalar state-mean distance (rad) | embedding state-mean distance |
|---|---:|---:|---:|---|---:|---:|
| F1a_s0 | 0 | 0 | 0 | True | 0 | 1.414 |
| F1a_s2deg | 0.3357 | 0.3357 | 0.3357 | False | 0.002281 | 1.417 |
| F1b_s0 | 0 | 0 | 0 | True | 0 | 1.414 |
| F1b_s2deg | 0.3254 | 0.3254 | 0.3254 | False | 0.001894 | 1.414 |
| F1c_s0 | 0.9911 | 0.9911 | 0.9911 | False | 1.047 | 1.414 |
| F1c_s2deg | 0.9824 | 0.9824 | 0.9824 | False | 1.046 | 1.412 |

In F1a/F1b, two states that the coordinate-wise embedding separates (distance 1.414) have identical per-frame summaries (scalar distance 0). Noiseless, the scan returns BC 0 via the variance guard. With 2° noise it returns 0.3357, a value of the size produced by unimodal noise (F6a Gaussian control: 0.3348). The positive control F1c, where the states differ in circular mean, gives 0.9911. The two-state structure of F1a/F1b is **not represented** in the scan output. This is shown for these constructions only.

**F6: BC non-specificity carried through the angular path.**

| Input | BC (REF) | BC (R1) |
|---|---:|---:|
| F6a_rare_p0.5 | 0.9581 | 0.9581 |
| F6a_rare_p0.1 | 0.9622 | 0.9622 |
| F6a_rare_p0.02 | 0.8392 | 0.8392 |
| F6a_rare_p0.005 | 0.7436 | 0.7436 |
| F6a_exponential | 0.5769 | 0.5769 |
| F6a_lognormal_s0.8 | 0.5122 | 0.5122 |
| F6a_gamma_k1 | 0.5698 | 0.5698 |
| F6a_chi2_df1 | 0.5898 | 0.5898 |
| F6a_gaussian | 0.3348 | 0.3348 |

| F6b case | BC R2 (coordinates) | BC P1 | R2 − scalar REF |
|---|---:|---:|---:|
| F6b_rare_p0.5 | 0.9581 | 0.9581 | -5.44e-15 |
| F6b_rare_p0.1 | 0.9622 | 0.9622 | 5.55e-16 |
| F6b_rare_p0.02 | 0.8392 | 0.8392 | 2.11e-15 |
| F6b_rare_p0.005 | 0.7436 | 0.7436 | 1.89e-15 |

Rare excursions and skewed unimodal inputs give BC values above those of the Gaussian control (KNOWN_ISSUES M3). The angular path reproduces the scalar values. The replay of `bc_check.py` draws matches the values documented in KNOWN_ISSUES M3 to three decimals; this is a consistency check, not proof of identical samples.

## 3. Numerical conditioning and sensitivity (sampled)

**F2: paired rotations** (24 cases). Max |ΔBC| vs α = 0: R1 9.44e-16, R2 9.44e-16, P1 7.31e-08. Max non-wrap deviation change (R1): 2.66e-15 rad. Wrap events: 0. Deviations within 1e-6 of ±π: 0.
F2c was intended to stress the ±π wrap but produced no deviation near ±π, so wrap behaviour under rotation is **not exercised by F2**. F4, below, does place deviations near ±π; its BC range across the 8 rotations was at most 1.78e-15.

**F3: within-frame near-zero resultant** (R1 frame-0 mean shift under a perturbation δ of one torsion; BC per path).

| ε | σ_n | per-frame R | shift δ=+1e-9 | shift δ=+1e-6 | BC R1 | BC R2 | BC P1 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| +1e-01 | 0e+00 | 5.00e-02 | 5.00e-10 | 5.00e-07 | 0.0000 | 0.0000 | 0.0000 |
| +1e-01 | 1e-06 | 5.00e-02 | 5.00e-10 | 5.00e-07 | 0.3545 | 0.3545 | 0.3563 |
| +1e-01 | 2e-02 | 5.00e-02 | 5.00e-10 | 5.00e-07 | 0.3712 | 0.3712 | 0.3712 |
| +1e-03 | 0e+00 | 5.00e-04 | 5.00e-10 | 5.00e-07 | 0.0000 | 0.0000 | 0.0000 |
| +1e-03 | 1e-06 | 5.00e-04 | 5.00e-10 | 5.00e-07 | 0.3293 | 0.3293 | 0.3301 |
| +1e-03 | 2e-02 | 5.00e-04 | 5.00e-10 | 5.00e-07 | 0.3626 | 0.3626 | 0.3626 |
| +1e-06 | 0e+00 | 5.00e-07 | 5.41e-10 | 2.99e+00 | 0.0000 | 0.0000 | 0.0000 |
| +1e-06 | 1e-06 | 5.00e-07 | 4.99e-10 | 2.71e+00 | 0.2976 | 0.2976 | 0.2964 |
| +1e-06 | 2e-02 | 5.00e-07 | 5.27e-10 | 6.13e-03 | 0.3261 | 0.3261 | 0.3261 |
| +1e-09 | 0e+00 | 5.00e-10 | 3.14e-02 | 3.14e+00 | 0.0000 | 0.0000 | 0.0000 |
| +1e-09 | 1e-06 | 5.00e-10 | 3.14e-02 | 3.14e+00 | 0.3487 | 0.3475 | 0.3127 |
| +1e-09 | 2e-02 | 5.00e-10 | 2.73e+00 | 3.14e+00 | 0.3758 | 0.3758 | 0.2571 |
| +0e+00 | 0e+00 | 1.49e-16 | 3.11e+00 | 3.11e+00 | 0.0000 | 0.0000 | 0.0000 |
| +0e+00 | 1e-06 | 5.55e-17 | 3.11e+00 | 3.11e+00 | 0.8916 | 0.3807 | 0.2233 |
| +0e+00 | 2e-02 | 5.55e-17 | 4.43e-01 | 4.43e-01 | 0.7514 | 0.2155 | 0.2887 |
| -1e-09 | 0e+00 | 5.00e-10 | 4.74e-08 | 5.83e-07 | 0.0000 | 0.0000 | 0.0000 |
| -1e-09 | 1e-06 | 5.00e-10 | 3.80e-08 | 5.70e-07 | 0.3301 | 0.3325 | 0.2569 |
| -1e-09 | 2e-02 | 5.00e-10 | 5.35e-08 | 4.40e-07 | 0.3517 | 0.3517 | 0.2170 |
| -1e-06 | 0e+00 | 5.00e-07 | 4.80e-10 | 5.00e-07 | 0.0000 | 0.0000 | 0.0000 |
| -1e-06 | 1e-06 | 5.00e-07 | 5.50e-10 | 5.00e-07 | 0.3464 | 0.3464 | 0.3473 |
| -1e-06 | 2e-02 | 5.00e-07 | 5.00e-10 | 5.00e-07 | 0.3237 | 0.3237 | 0.3237 |
| -1e-03 | 0e+00 | 5.00e-04 | 5.00e-10 | 5.00e-07 | 0.0000 | 0.0000 | 0.0000 |
| -1e-03 | 1e-06 | 5.00e-04 | 5.00e-10 | 5.00e-07 | 0.3945 | 0.3945 | 0.3958 |
| -1e-03 | 2e-02 | 5.00e-04 | 5.00e-10 | 5.00e-07 | 0.3967 | 0.3967 | 0.3967 |
| -1e-01 | 0e+00 | 5.00e-02 | 5.00e-10 | 5.00e-07 | 0.0000 | 0.0000 | 0.0000 |
| -1e-01 | 1e-06 | 5.00e-02 | 5.00e-10 | 5.00e-07 | 0.4181 | 0.4181 | 0.4180 |
| -1e-01 | 2e-02 | 5.00e-02 | 5.00e-10 | 5.00e-07 | 0.3422 | 0.3422 | 0.3422 |

Sampled observations:
- A perturbation of 1e-6 rad moved the frame-0 mean by more than 1 rad at ε ∈ [0.0, 1e-09, 1e-06].
- A perturbation of 1e-9 rad moved it by more than 1e-3 rad at ε ∈ [0.0, 1e-09].
- These flips are expected mathematically: the mean jumps by about π when a perturbation of about R reverses the sign of the resultant. They are not round-off.
- Max |BC P1 − BC R2| per sampled ε (float32 coordinates vs f64): -1e-01: 0.000114, -1e-03: 0.00133, -1e-06: 0.000842, -1e-09: 0.135, +0e+00: 0.157, +1e-09: 0.119, +1e-06: 0.00125, +1e-03: 0.000734, +1e-01: 0.00173.
- The scan reports no resultant length, so these cases cannot be recognized from its output. Behaviour between sampled ε values is not established.

**F4: across-frame cancellation and antipodal modes** (α = 0; BC per path).

| Case | across R (R1) | BC R1 | BC R2 | BC P1 |
|---|---:|---:|---:|---:|
| F4_p0.50_e0_n0_a+000 | 6.123e-17 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.50_e0_n1deg_a+000 | 2.604e-04 | 0.9906 | 0.9906 | 0.9906 |
| F4_p0.50_e+1e-06_n0_a+000 | 5.000e-07 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.50_e+1e-06_n1deg_a+000 | 1.284e-04 | 0.9906 | 0.9906 | 0.9906 |
| F4_p0.50_e-1e-06_n0_a+000 | 5.000e-07 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.50_e-1e-06_n1deg_a+000 | 4.507e-04 | 0.9906 | 0.9906 | 0.9906 |
| F4_p0.50_e+1e-03_n0_a+000 | 5.000e-04 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.50_e+1e-03_n1deg_a+000 | 2.429e-05 | 0.9906 | 0.9906 | 0.9906 |
| F4_p0.50_e-1e-03_n0_a+000 | 5.000e-04 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.50_e-1e-03_n1deg_a+000 | 6.680e-04 | 0.9906 | 0.9906 | 0.9906 |
| F4_p0.50_e+1e-01_n0_a+000 | 4.998e-02 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.50_e+1e-01_n1deg_a+000 | 5.046e-02 | 0.9906 | 0.9906 | 0.9906 |
| F4_p0.50_e-1e-01_n0_a+000 | 4.998e-02 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.50_e-1e-01_n1deg_a+000 | 5.024e-02 | 0.9905 | 0.9905 | 0.9905 |
| F4_p0.49_e0_n0_a+000 | 2.000e-02 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.49_e0_n1deg_a+000 | 2.001e-02 | 0.4875 | 0.4875 | 0.4875 |
| F4_p0.49_e+1e-06_n0_a+000 | 2.000e-02 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.49_e+1e-06_n1deg_a+000 | 1.999e-02 | 0.4878 | 0.4878 | 0.4878 |
| F4_p0.49_e-1e-06_n0_a+000 | 2.000e-02 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.49_e-1e-06_n1deg_a+000 | 2.000e-02 | 0.9352 | 0.9352 | 0.9352 |
| F4_p0.49_e+1e-03_n0_a+000 | 2.001e-02 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.49_e+1e-03_n1deg_a+000 | 2.001e-02 | 0.5173 | 0.5173 | 0.5173 |
| F4_p0.49_e-1e-03_n0_a+000 | 2.001e-02 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.49_e-1e-03_n1deg_a+000 | 2.003e-02 | 0.9906 | 0.9906 | 0.9906 |
| F4_p0.49_e+1e-01_n0_a+000 | 5.382e-02 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.49_e+1e-01_n1deg_a+000 | 5.345e-02 | 0.9905 | 0.9905 | 0.9905 |
| F4_p0.49_e-1e-01_n0_a+000 | 5.382e-02 | 0.9911 | 0.9911 | 0.9911 |
| F4_p0.49_e-1e-01_n1deg_a+000 | 5.419e-02 | 0.9905 | 0.9905 | 0.9905 |
| F4_p0.45_e0_n0_a+000 | 1.000e-01 | 0.9914 | 0.9914 | 0.9914 |
| F4_p0.45_e0_n1deg_a+000 | 9.998e-02 | 0.4462 | 0.4462 | 0.4462 |
| F4_p0.45_e+1e-06_n0_a+000 | 1.000e-01 | 0.9914 | 0.9914 | 0.9914 |
| F4_p0.45_e+1e-06_n1deg_a+000 | 9.998e-02 | 0.4470 | 0.4470 | 0.4470 |
| F4_p0.45_e-1e-06_n0_a+000 | 1.000e-01 | 0.9914 | 0.9914 | 0.9914 |
| F4_p0.45_e-1e-06_n1deg_a+000 | 9.999e-02 | 0.4472 | 0.4472 | 0.4472 |
| F4_p0.45_e+1e-03_n0_a+000 | 1.000e-01 | 0.9914 | 0.9914 | 0.9914 |
| F4_p0.45_e+1e-03_n1deg_a+000 | 9.999e-02 | 0.4451 | 0.4451 | 0.4451 |
| F4_p0.45_e-1e-03_n0_a+000 | 1.000e-01 | 0.9914 | 0.9914 | 0.9914 |
| F4_p0.45_e-1e-03_n1deg_a+000 | 9.999e-02 | 0.4467 | 0.4467 | 0.4467 |
| F4_p0.45_e+1e-01_n0_a+000 | 1.117e-01 | 0.9914 | 0.9914 | 0.9914 |
| F4_p0.45_e+1e-01_n1deg_a+000 | 1.120e-01 | 0.9909 | 0.9909 | 0.9909 |
| F4_p0.45_e-1e-01_n0_a+000 | 1.117e-01 | 0.9914 | 0.9914 | 0.9914 |
| F4_p0.45_e-1e-01_n1deg_a+000 | 1.118e-01 | 0.9909 | 0.9909 | 0.9909 |
| F4_control120_n0_a+000 | 5.000e-01 | 0.9911 | 0.9911 | 0.9911 |
| F4_control120_n1deg_a+000 | 4.997e-01 | 0.9901 | 0.9901 | 0.9901 |

Sampled observations:
- At p = 0.5, ε = 0 the across-frame resultant is 6.1e-17, so the across-frame mean is undefined. The deviations are [[-1.5707963267948966, 500], [1.5707963267948966, 500]]. BC is 0.9911 and did not depend on the arbitrary centering: a two-point set separated by π gives the same BC wherever the center falls.
- With 1° per-frame noise, these cases changed BC by more than 0.1 relative to the noiseless construction: p=0.45, ε=-1e-03: 0.991 → 0.447; p=0.45, ε=-1e-06: 0.991 → 0.447; p=0.45, ε=+0e+00: 0.991 → 0.446; p=0.45, ε=+1e-06: 0.991 → 0.447; p=0.45, ε=+1e-03: 0.991 → 0.445; p=0.49, ε=+0e+00: 0.991 → 0.487; p=0.49, ε=+1e-06: 0.991 → 0.488; p=0.49, ε=+1e-03: 0.991 → 0.517.
- **Draws are not paired across F4 configurations.** The 44 F4 configurations (p, ε, σ) used 44 distinct seeds. The minority-mode frame assignment and the noise draws therefore differ between configurations, including between +ε and −ε, and between noisy and noiseless cases. Differences between ±ε rows **cannot be attributed to the sign of ε**: draw-to-draw variation is not separated from ε. The noiseless BC values do not depend on the assignment (two-point sets with fixed counts), but every noisy value does.
- Proposed mechanism, consistent with the recorded deviations but not isolated by these fixtures: when the across-frame mean lies near the majority mode, the antipodal minority mode lies near ±π, and per-frame noise can place its members on both sides of the wrap. How often this happens for a given p, ε and σ is not established. A paired ±ε control (shared assignment and noise) is specified in `docs/ORDERED_ANGULAR_REPRESENTATION_PROTOCOL.md`.

**κ channel.** Non-zero `bc_kappa`: R2 0/105, P1 0/105. Cases where the score came from κ (P1): 0. With κ = 1.2 rad and these coordinate scales, float32 quantization did not lift the κ variance above the absolute guard. Other geometries were not tested.

**Float32 path vs f64 path.** |BC P1 − BC R2| over 105 cases: median 2.41e-11, max 0.157; exactly equal in 14 cases. The large differences occur in the F3 near-cancellation cases above.

## 4. Implications for this component only, and decision record

**Harness correctness:** established for the tested paths. Every check in §1 passed: bitwise replication of production, recovery, REF agreement where conditioned, analytical controls, and the F6 replay.

**Information loss (scientific limitation of the representation):** demonstrated for the F1 constructions. Coordinate-wise state changes that preserve the window circular mean are invisible to the scan. BC non-specificity (F6) is carried unchanged into the angular path. Absence of loss for other constructions is not shown.

**Numerical sensitivity:** in the sampled cases, the scan output changed substantially under two conditions:
- per-frame resultant comparable to input perturbations: F3, including float32 quantization;
- an antipodal minority mode placed at the wrap cut under 1° noise: F4.
Neither condition is reported by the scan. Rotation invariance held at round-off level in every sampled case.

**Decision: CHANGE** (protocol §8). The losses and sensitivities above affect configurations that cannot be detected from the scan output. Candidate changes, each needing its own protocol and validation (none adopted here):
1. a per-coordinate (sin, cos) representation instead of the per-frame scalar mean;
2. carrying and reporting within- and across-frame resultant lengths;
3. a centering or summary rule that does not place a minority mode at the wrap cut.

Until a change is validated, a low BC from this scan should not be read as absence of multimodality, and a high BC should not be read as evidence of two states.
This decision concerns this component only. Geometry descriptors and retrieval are not affected.

**Limitations:**
- SYNTHETIC constructions with constant κ = 1.2 rad and W ∈ {5, 8};
- a closed, finite fixture matrix;
- sensitivity reported only at sampled ε, p and noise values;
- F2 did not exercise ±π wraps under rotation;
- no real-MD data were used.

**Next concrete decision required:** choose which candidate change, if any, to specify as its own protocol (with equal-information baselines). The alternative is to restrict the current scan to an explicitly documented, non-inferential screening role.

## 5. Post-execution notes (recorded, checks unchanged)

- Post-run editorial changes to this report (wording only; no rerun, no change to fixtures or results files) are logged separately in `results/SYNTHETIC_angular_characterization_editorial_log.json`. The run provenance file is the original and is not modified by editorial changes.
- The protocol (§5, F3) stated that |ε| = 1e-9 was "expected to be numerically undefined". By the pre-declared diagnostic, with exact f64 inputs, it is not: R = 5.00e-10, b = 6.7e-07 rad. The observed sensitivity there comes from input perturbations of about R, not from round-off. That expectation statement was wrong; no check depended on it.
- F2c did not reach its intended wrap stress (see §3).
- Artifact layout: regenerable binaries (fixtures, arrays, wheel, isolated venv) are under `build/`. Hashes of the fixtures are in `fixtures/SYNTHETIC_fixtures.json`; hashes of the arrays are in the results files.
- During development, before any fixture was generated, the harness trace builder and the independent NumPy dihedral each had a sign error. The recovery check detected both, and both were fixed (`test_independent_dihedral_detects_sign_error` guards against regression).
