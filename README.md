# TopoFold

**SE(3)-invariant local geometry of biomolecular backbones, in Rust with Python bindings.**

[![Rust](https://img.shields.io/badge/rust-1.80%2B-orange.svg)](https://www.rust-lang.org)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org)
[![unsafe forbidden](https://img.shields.io/badge/unsafe-forbidden-success.svg)](https://github.com/rust-secure-code/safety-dance)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](#license)
[![Status](https://img.shields.io/badge/status-research%20prototype-yellow.svg)](docs/KNOWN_ISSUES.md)

TopoFold computes discrete differential-geometry descriptors of Cα (and RNA phosphorus) traces — virtual bond angles (curvature κ), virtual torsions (τ), local writhe and a Cβ orientation angle — and uses them for alignment-free comparison, indexing and screening of conformational ensembles. The descriptors depend only on internal geometry, so no structural superposition is needed, and the descriptors of a residue window depend only on that window and its immediate neighbours.

> **Project status: research prototype under scientific re-validation.**
> The software runs and its core geometry is tested. **No current claims about biomolecular dynamics, thermodynamics, kinetics, pocket detection, allostery or biological function have been validated.** The one completed held-out evaluation (October 2026, see *Status at a glance* below) concerns geometric local-conformation retrieval only. XCL1 remains a real static structural example, but its current distance-based interpretation requires correction. See the correction notice below and [`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md).

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
| Writhe (total, local window) | Discrete Gauss integral over segment pairs via solid angles (van Oosterom–Strackee) | Implemented; normalisation and sign corrected on `scientific-rebuild` (M1, unreleased), regression-tested against independent polygonal and numerical Gauss references. Values differ from v0.8.3 by a factor of −2 |
| Cβ orientation θβ | Angle of the Cα→Cβ vector in the local Frenet frame | Implemented; a backbone descriptor, **not** a side-chain rotamer descriptor |
| RNA invariants | κ, τ on P atoms; base vector C1′→N9/N1 | Implemented; base vector lies along the glycosidic bond and cannot measure χ (syn/anti), fix pending |
| Window bimodality scan | Sarle's bimodality coefficient of window-averaged descriptors over an ensemble | Implemented as a **heuristic**. Circular means for τ and θβ on `scientific-rebuild` (M2, unreleased). BC is not a statistical test: no calibrated null model or multiple-testing control (M3). Characterized on SYNTHETIC fixtures in WP-A (2026-10-04, decision CHANGE): coordinate-wise changes that keep the window circular mean are invisible, and the output is numerically sensitive near cancellation, so a low BC does not mean one state and a high BC does not mean two. The later replica-level studies (REPCHAR-1, PERM-1) concern JS-distance comparisons and do not validate BC |
| Subcurve search | Full-frame `query()`: approximate VP-tree cascade over writhe spectra + (κ, τ) Fréchet. Subcurve `query_subcurve()`: exhaustive search using invariant Fréchet distance | Implemented; full-frame cascade is **approximate** (recall measured only on AdK development data, LE-1); subcurve is exhaustive; the legacy Fréchet coupling was deprioritized for same-sequence windows (LE-1) |
| Correlation network | Gaussian (covariance-based) mutual information between per-residue (κ, τ, θβ) | Implemented; no significance testing yet |

Numerical tests (invariance, locality, reference comparisons) are not scientific validation. No capability above has yet been validated on independent data; see the capability matrix in [`docs/ENGINE_MISSION_AND_ROADMAP.md`](docs/ENGINE_MISSION_AND_ROADMAP.md).

Implementation notes:
- **Locality.** Window descriptors depend only on the atoms of that window (tested to 1e-12), so motions elsewhere, including flexible termini, do not change them. This does not extend to global writhe, full-frame search or any pocket claim.
- **Invariance precision.** SE(3) invariance is tested at 1e-6 rad in f64. The Python bindings take float32 coordinates and convert them to f64, so input precision is float32.
- **Memory.** The Rust DCD reader works frame by frame, but `read_dcd`, the bimodality profile and the index hold all frames in memory. The bindings copy input arrays; they are not zero-copy.
- **Parameters.** Window size, BC threshold, coarse radius, writhe tolerance, gap tolerance and MI regularizer are user-set parameters with defaults.

## Status at a glance (reconciled 2026-10-13)

| Category | What it covers |
| --- | --- |
| **Long-term mission and commercial objectives** | A reproducible engine for locating, comparing and retrieving conformational changes, with oncology-relevant applications as the long-term goal. Commercial targets (blind pocket identification, automated reports, throughput, no per-dataset tuning) are **hypotheses**, not capabilities. See [`docs/ENGINE_MISSION_AND_ROADMAP.md`](docs/ENGINE_MISSION_AND_ROADMAP.md) §1 and §11, and [`docs/ONCOLOGY_VALIDATION_REQUIREMENTS.md`](docs/ONCOLOGY_VALIDATION_REQUIREMENTS.md). |
| **Experimentally supported** | One pre-registered held-out evaluation of **geometric** local-conformation retrieval (below). Nothing biological, clinical or pocket-related. |
| **Harness-only implementations** | Engineering benchmarks and the V2b single-query RMSD search variant live in `benchmarks/audit/` as research harness code. They are **not** in the library or the Python API and are not production-ready. |
| **Unresolved limitations** | Everything in [`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md), the capability matrix in the roadmap, and the open engineering items below. |

### Experimentally supported: held-out geometric retrieval (October 2026)

One pre-registered, held-out evaluation has been completed: **Cα-only local-conformation retrieval** on 30 CATH domains from mdCATH. The data are third-party MD (CHARMM22*, 320/348 K) spanning 30 superfamilies. The task is cross-replica retrieval from a shared starting configuration, judged against two geometric targets reported separately. Results, declared outcomes and corrections: [`benchmarks/audit/wpd_local_retrieval/REPORT_HELDOUT.md`](benchmarks/audit/wpd_local_retrieval/REPORT_HELDOUT.md) and [`docs/WPD_HELDOUT_MDCATH_CLOSEOUT_ADDENDUM_2026-10-12.md`](docs/WPD_HELDOUT_MDCATH_CLOSEOUT_ADDENDUM_2026-10-12.md).

- **A**, a residue-aligned (κ, τ) window distance, met the declared criterion on the angular φ/ψ target.
- **B2**, window Cα RMSD, met it on the Cartesian peptide-plane target.
- Each is a candidate for its own geometric task; neither is preferred across both. B2 at W = 8 on φ/ψ is inconclusive, and W = 5 does not establish equivalence (closeout addendum).
- These are within-dataset generalization results for a geometric benchmark. They are **not** validation of the legacy index, any detector, biological function, or any clinical or cancer-treatment use.

### Harness-only engineering results (development data; October 2026)

| Study | Outcome under its frozen rules | Report |
| --- | --- | --- |
| EB-1 engineering benchmark, Amendments 1–2 | EB-1 and Amendment 1: **insufficient evidence**. Amendment 2: decision **engine**, from the projected single-query B2 latency at the P2 target (an extrapolation, not a P2 measurement). The loader and A were within requirements as projected. | [`benchmarks/audit/eb1/`](benchmarks/audit/eb1/README.md) |
| BQ-1 B2 single-query optimization | Correctness preserved by V0, V1 and V2b; frozen V2 rejected on degenerate geometry. **V2b measured faster than the harness baseline on every tested workload.** Overall acceptance **not met**: the preparation projection was unsupported. | [`benchmarks/audit/b2opt/REPORT.md`](benchmarks/audit/b2opt/REPORT.md) |
| BQ-2 scaling confirmation of unchanged V2b | Query-time and memory models supported and confirmed, including on SYNTHETIC performance-only stress rows. Preparation model unconfirmed on the largest stress cells. **Conclusion C: insufficient evidence.** | [`benchmarks/audit/bq2/REPORT.md`](benchmarks/audit/bq2/REPORT.md) |
| BQ-3 direct-P2 qualification | First attempt **stopped at memory preflight**. After a memory recheck, the direct run on SYNTHETIC P2-dimension rows (performance-only, repeated frames) gave:<br>• **Prepared-database operation: qualified.** All query latencies ≤ 10 s, peak memory within 8 GiB, verification passed.<br>• **Preparation + first query: did not meet 10 s.** That span excludes input loading and row construction.<br>BQ-2's conclusion C is unchanged. | [`docs/B2_DIRECT_P2_QUALIFICATION_CLOSEOUT_2026-10-13.md`](docs/B2_DIRECT_P2_QUALIFICATION_CLOSEOUT_2026-10-13.md), [`benchmarks/audit/bq3/REPORT_SYNTHETIC_P2.md`](benchmarks/audit/bq3/REPORT_SYNTHETIC_P2.md) |
| CMP-1 bounded competitive evaluation (V2b vs MDTraj 1.11.2 and MDAnalysis 2.10.0, one thread) | **No decisive advantage** under the frozen rules (common F32 track):<br>• **MDTraj** passed the common accuracy track on two real workloads and was about 3–4× faster than V2b there.<br>• **MDTraj** exceeded that track's 1e-3 Å distance tolerance on rare pairs in the other three workloads. Its rankings stayed equivalent.<br>• **V2b** was about 2–3× faster than the remaining eligible comparator, MDAnalysis per-pair QCP.<br>• **Integration of V2b is not supported by this evidence.**<br>Bounded comparison on development and SYNTHETIC data, not a field-wide best-tool claim; CPPTRAJ was not evaluated. | [`benchmarks/audit/cmp1/REPORT.md`](benchmarks/audit/cmp1/REPORT.md) |

What these results do and do not show:
- V2b computes the **same** exhaustive proper-rotation Cα RMSD top-k as the baseline, faster on the tested development and SYNTHETIC workloads. Per-workload speed-ups are in the BQ-1 and BQ-2 reports.
- They show no new biological capability, no technological novelty over existing RMSD software, and no production readiness.
- A successful engineering qualification would not establish any of those either.
- CMP-1 is the only equal-information comparison against external tools (MDTraj, MDAnalysis) so far. It is bounded to harness single-query Cα RMSD search on development and SYNTHETIC workloads; the library has not been compared.

### Next tasks

The development programme, evidence levels (A–E) and the current queue are in
[`docs/ENGINE_MISSION_AND_ROADMAP.md`](docs/ENGINE_MISSION_AND_ROADMAP.md) §12; recorded decisions are in
[`docs/DECISION_LOG.md`](docs/DECISION_LOG.md). In short:
- **KRAS C1** (RAS switch-I structural label recovery) is blocked and unfrozen pending the Rascore label-data terms
  ([`docs/C1_DECISION_NOTE_2026-10-13.md`](docs/C1_DECISION_NOTE_2026-10-13.md),
  [`docs/C1_EVALUATION_DESIGN_DRAFT.md`](docs/C1_EVALUATION_DESIGN_DRAFT.md)). A positive result would be structural
  (level B), not functional.
- **V2b** stays an audit artifact; CMP-1 does not support integration
  ([`docs/CMP1_CLOSEOUT_2026-10-13.md`](docs/CMP1_CLOSEOUT_2026-10-13.md)).
- The MI network's double regularization is fixed (2026-10-22, KNOWN_ISSUES M7); its outputs can differ from earlier
  versions. This is a numerical correction and does not validate allostery.
- An optional Cα input pipeline (`python/topofold_inputs`: MDAnalysis adapter → core κ/τ with explicit validity →
  NPZ/JSON/inspection report) is engineering-ready, with a development run on local AdK data
  ([`benchmarks/audit/dev_ca_pipeline/`](benchmarks/audit/dev_ca_pipeline/README.md)). It is not a scientific
  validation.
- Ensemble comparison, kept in three separate evidence classes:
  - **Development evidence (real MD, descriptive):** DEV-ENS-1 compared one AdK open-start and one closed-start run
    (single runs, so the starting structure is confounded); DEV-REP-1 compared 5 vs 5 mdCATH `1ha8A00` replicas at
    320 vs 348 K, and no window met the complete-separation criterion
    ([`benchmarks/audit/dev_rep_compare/`](benchmarks/audit/dev_rep_compare/INTERPRETATION_NOTE_2026-10-23.md)).
  - **SYNTHETIC calibration:** REPCHAR-1 (complete-separation screen: no null calls, low sensitivity under temporal
    correlation), PERM-1 (exact replica permutation: near-nominal under identical distributions, more sensitive) and
    NULLSTRESS-1 (the permutation also rejects when only temporal correlation differs), in
    [`benchmarks/audit/`](benchmarks/audit/README.md).
  - **Biological validation:** none yet. ABL remains the preferred route.
  - The real-data permutation design for `1ha8A00` is **deferred, not scored**
    ([`docs/DEV_PERM_1HA8A00_EVALUATION_DESIGN_DRAFT.md`](docs/DEV_PERM_1HA8A00_EVALUATION_DESIGN_DRAFT.md)). The
    permutation method responds to histogram differences caused by temporal correlation alone, so it is not validated
    for marginal-geometry comparison.

## What TopoFold does not (yet) do

- It does not detect cryptic pockets: a bimodal backbone window is not a pocket. Cavity analysis and an external benchmark are required.
- It does not estimate free energies or kinetics.
- Outside the single geometric retrieval task below, it has not been compared fairly with local-feature baselines (local RMSD, Cα–Cα distances, φ/ψ, TICA/VAMP on the same residues). On the project's own synthetic benchmark, simple local Cα–Cα distances separate the states as well as TopoFold's (κ, τ) features (see `benchmarks/audit/`).
- The legacy (κ, τ) Fréchet index, its writhe cascade, the bimodality scan and every detector remain unvalidated; the held-out evaluation did not test them.

---

## Installation

Prerequisites: a recent stable Rust toolchain, Python 3.9+ and [maturin](https://www.maturin.rs/). `Cargo.toml` declares Rust 1.80, but the code uses APIs stabilized in Rust 1.87–1.88, so the declared minimum is not honoured (KNOWN_ISSUES m12).

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
benchmarks/         controlled synthetic tests (see benchmarks/README.md); audit harnesses, protocols' code and results in benchmarks/audit/ (index: benchmarks/audit/README.md)
docs/               working manuscript draft, known issues, architecture decision records
apps/               Streamlit viewer
assets/             v0.8.3 figures (superseded, kept for the record)
```

## Roadmap

1. Fix the issues in [`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md), each with a regression test.
2. Engine-choice review after CMP-1 (see *Next tasks*). Direct P2 qualification (BQ-3) and a bounded competitive evaluation (CMP-1) were completed with scoped outcomes.
3. Validate on real, independently generated MD trajectories, with ground truth that does not come from TopoFold itself, for each intended application.
4. Measure full-frame index recall against brute force, and performance against vectorised MDTraj/MDAnalysis on identical tasks.
5. Release with a Zenodo DOI and a preprint.

Detailed status and sequencing: [`docs/ENGINE_MISSION_AND_ROADMAP.md`](docs/ENGINE_MISSION_AND_ROADMAP.md).

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
