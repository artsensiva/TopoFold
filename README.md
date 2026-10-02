# TopoFold

**SE(3)-invariant local geometry of biomolecular backbones, in Rust with Python bindings.**

[![Rust](https://img.shields.io/badge/rust-1.80%2B-orange.svg)](https://www.rust-lang.org)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org)
[![unsafe forbidden](https://img.shields.io/badge/unsafe-forbidden-success.svg)](https://github.com/rust-secure-code/safety-dance)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](#license)
[![Status](https://img.shields.io/badge/status-research%20prototype-yellow.svg)](docs/KNOWN_ISSUES.md)

TopoFold computes discrete differential-geometry descriptors of Cα (and RNA phosphorus) traces — virtual bond angles (curvature κ), virtual torsions (τ), local writhe and a Cβ orientation angle — and uses them for alignment-free comparison, indexing and screening of conformational ensembles. The descriptors depend only on internal geometry, so no structural superposition is needed, and the descriptors of a residue window depend only on that window and its immediate neighbours.

> **Project status: research prototype under scientific re-validation.**
> The software runs and its core geometry is tested. **No current claims about biomolecular dynamics, thermodynamics, kinetics, pocket detection, or comparative method superiority have yet been validated on independently generated trajectory data.** XCL1 remains a real static structural example, but its current distance-based interpretation requires correction. See the correction notice below and [`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md).

---

## Correction notice (October 2026)

Versions up to and including **v0.8.3** of this README, `docs/MANUSCRIPT_DRAFT.md` and `docs/whitepaper.md` made claims that are **withdrawn**:

- Benchmarks were described as "authentic explicit-solvent MD", "real BPTI MD" or "real crystal data". In fact, nearly all ensembles used for the reported numbers (BPTI, Abl1 DFG flip, SARS-CoV-2 Mpro, α-synuclein, adenine riboswitch, PROTAC complexes, rotamer test) were **generated or interpolated by the benchmark scripts** — by interpolation between experimental structures, rigid rotations and added noise, or fully synthetically — with the states and labels built in. They are controlled synthetic tests, not MD. The XCL1 benchmark uses real experimental NMR structures compared statically, but overstates the evidence (the distance metric used is not SE(3)-invariant). See [`benchmarks/README.md`](benchmarks/README.md).
- The reported "activation free-energy barriers" (3.43–6.17 kBT) are artefacts of the estimator and of the data generators, not physical quantities.
- The PROTAC comparison used a wrong PDB entry (5T3E is a non-ribosomal peptide synthetase domain, not a PROTAC ternary complex); the "non-productive" complex was built from 5T35 coordinates and its coupling was hard-coded. This benchmark is withdrawn.
- The Cβ orientation angle θβ cannot detect side-chain rotamer changes (χ1 rotation does not move Cβ). The "rotamer gating" claim is withdrawn.
- Claims of superiority over PCA, dPCA, TICA, Foldseek, PocketMiner and AlphaFold were not supported by fair comparisons and are withdrawn.
- The writhe implementation (`crates/topofold-core/src/writhe.rs`) in v0.8.3 divided by 4π and used a reversed sign convention, returning −0.5 times the standard Gauss-integral writhe. **Corrected on the scientific-rebuild branch:** sign and normalization validated against an independently implemented polygonal reference and converged numerical Gauss quadrature. The public v0.8.3 baseline contained the −0.5 defect.

All figures in `assets/` belong to v0.8.3 and carry the old labels; they are kept only for the record and will be regenerated after the fixes listed in [`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md).

---

## What TopoFold computes

| Descriptor | Definition | Status |
| --- | --- | --- |
| Curvature κᵢ | Turning angle between consecutive Cα–Cα tangents, κᵢ ∈ [0, π] (virtual bond angle supplement) | Implemented, tested (SE(3) invariance, reflection) |
| Torsion τᵢ | Signed dihedral between consecutive osculating planes, τᵢ ∈ (−π, π] (Cα virtual torsion) | Implemented, tested |
| Writhe (total, local window) | Discrete Gauss integral over segment pairs via solid angles (van Oosterom–Strackee) | Implemented; **normalisation issue** (audit reproduction yields −0.500 ratio vs numerical reference; normalization/sign convention validation pending Phase 3) |
| Cβ orientation θβ | Angle of the Cα→Cβ vector in the local Frenet frame | Implemented; a backbone descriptor, **not** a side-chain rotamer descriptor |
| RNA invariants | κ, τ on P atoms; base vector C1′→N9/N1 | Implemented; base vector lies along the glycosidic bond and cannot measure χ (syn/anti), fix pending |
| Window bimodality scan | Sarle's bimodality coefficient of window-averaged descriptors over an ensemble | Implemented as a **heuristic**; no circular statistics, no calibrated null model yet |
| Subcurve search | Full-frame `query()`: approximate VP-tree cascade over writhe spectra + (κ, τ) Fréchet. Subcurve `query_subcurve()`: exhaustive search using invariant Fréchet distance | Implemented; full-frame cascade is **approximate** (recall not yet measured); subcurve is exhaustive |
| Correlation network | Gaussian (covariance-based) mutual information between per-residue (κ, τ, θβ) | Implemented; no significance testing yet |

## What TopoFold does not (yet) do

- It does not detect cryptic pockets: a bimodal backbone window is not a pocket. Cavity analysis and an external benchmark are required.
- It does not estimate free energies or kinetics.
- It has not been compared fairly with local-feature baselines (local RMSD, Cα–Cα distances, φ/ψ, TICA/VAMP on the same residues). On the project's own synthetic benchmark, simple local Cα–Cα distances separate the states as well as TopoFold's (κ, τ) features (see `benchmarks/audit/`).

---

## Installation

Prerequisites: Rust 1.80+, Python 3.9+, [maturin](https://www.maturin.rs/).

```bash
git clone https://github.com/artsensiva/TopoFold.git
cd TopoFold
pip install maturin
maturin develop --release
```

Run the tests:

```bash
cargo test --workspace
pytest tests/
```

Note: some Python tests are skipped when benchmark data files are absent.

## Quick start

```python
import numpy as np
import topofold as tf

# Cα coordinates of one structure, shape (N, 3), float32
ca = tf.read_pdb("structure.pdb", chain="A")

# Per-residue curvature, torsion and local writhe spectrum
kappa, tau, writhe_spectrum = tf.compute_invariants(ca)

# Ensemble of frames, shape (F, N, 3), float32
traj = tf.read_dcd("trajectory.dcd")

# Window-wise bimodality profile (screening heuristic, see KNOWN_ISSUES)
scores, bc_tau, bc_kappa = tf.compute_bimodality_profile(traj, window_size=8)

# Exhaustive subcurve nearest-neighbour search across frames
# (query_subcurve uses invariant Fréchet distance; full-frame query() uses VP-tree cascade)
index = tf.ConformationalIndex(window_size=8, coarse_radius=0.5)
index.fit(traj)
hits = index.query_subcurve(start_res=9, end_res=17, query_coords=ca, k=10)
```

Treat the bimodality scan as a way to rank windows for inspection, not as a statistical test.

---

## Repository layout

```
crates/
  topofold-core/    κ, τ, writhe, θβ, RNA invariants, bimodality, correlation network
  topofold-index/   VP-tree and discrete Fréchet search
  topofold-io/      DCD and multi-model PDB readers
  topofold-python/  PyO3 bindings
benchmarks/         controlled synthetic tests (see benchmarks/README.md) and audit scripts
docs/               working manuscript draft, known issues, architecture decision records
apps/               Streamlit viewer
assets/             v0.8.3 figures (superseded, kept for the record)
```

## Roadmap

1. Fix the issues in [`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md), each with a regression test.
2. Re-run the controlled synthetic tests with fair, equal-information baselines.
3. Validate on real, independently generated MD trajectories with ground truth that does not come from TopoFold itself.
4. Measure index recall against brute force and performance against vectorised MDTraj/MDAnalysis.
5. Release with a Zenodo DOI and a preprint.

## Development note

TopoFold is developed by Artem Galukhin with substantial help from AI coding assistants (Gemini in VS Code for versions up to 0.8.3). The v0.8.3 benchmark scripts and claims were produced in that workflow without independent verification, which led to the errors listed above. Rules for AI-assisted development of this repository are in [`AGENTS.md`](AGENTS.md).

## Citation

If you use TopoFold, please cite the software (see [`CITATION.cff`](CITATION.cff)):

```bibtex
@software{galukhin_topofold,
  author  = {Galukhin, Artem},
  title   = {TopoFold: SE(3)-invariant local geometry of biomolecular backbones},
  year    = {2026},
  url     = {https://github.com/artsensiva/TopoFold},
  version = {0.8.3}
}
```

## License

Licensed under either of [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at your option, as declared in `Cargo.toml` and `pyproject.toml`. (The v0.8.3 README mentioned AGPL-3.0 with a commercial option; no AGPL license file was ever included, and that statement is withdrawn.)

Contact: artem.galukhin@gmail.com
