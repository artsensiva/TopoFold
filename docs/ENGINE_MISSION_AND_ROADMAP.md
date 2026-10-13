# TopoFold Engine: Mission, Capabilities and Roadmap

**Status:** DRAFT for review (2026-10-04; WP-D status reconciled 2026-10-12; engineering status reconciled 2026-10-13)
**Baseline:** `295d739` (branch `scientific-rebuild`)

This document states what TopoFold is for, what it can do today and the evidence for that,
and what must happen next. Where it conflicts with older material (ADRs, whitepaper, v0.8.3
benchmark descriptions), this document and [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) take precedence.
Rules for contributors stay in [`../AGENTS.md`](../AGENTS.md).

---

## 1. Mission

> Build TopoFold into a reproducible geometry and search engine for biomolecular conformational
> ensembles: helping researchers locate, compare, retrieve, and investigate structural changes
> across protein and RNA trajectories, with traceable evidence and application-specific
> validation.

Three consequences follow:

- **The engine is the product; analyses are components.** Geometry, comparison/search and
  provenance are shared. Each analysis (e.g. the window bimodality scan) is one component with
  its own validation status. When one component fails, only that component's status changes;
  unrelated geometry or retrieval capabilities are not invalidated.
- **Validation is per application.** A capability validated for one question (e.g. locating
  local backbone changes) carries no claim for another (e.g. finding pockets).
- **Evidence must be traceable.** Every reported scientific number is produced by code and
  stored in a machine-readable results file (AGENTS.md rule 3).

## 2. Intended users and tasks

**Users**
- Computational biophysicists and structural biologists analysing MD ensembles or multi-model
  structure sets.
- Method developers comparing representations.
- Later, and only after application-specific validation: discovery teams investigating
  conformational change in specific targets.

**Concrete tasks the engine should support**

| ID | Task | Example |
|----|------|---------|
| U1 | Compute local geometric descriptors per residue, window and frame | κ, τ, local writhe along a trajectory |
| U2 | Locate regions whose local conformation changes across an ensemble, ranked for inspection | Which windows switch between two backbone states? |
| U3 | Compare two ensembles or states region by region | Apo vs ligand-bound simulation |
| U4 | Retrieve frames whose selected region resembles a query conformation | Frames where residues 30–45 resemble a reference structure |
| U5 | Compare static structures (e.g. NMR models, crystal forms) with the same descriptors | Metamorphic protein folds |
| U6 | Export traceable outputs | Inputs hashed, parameters recorded, results as JSON |

## 3. Shared infrastructure: current state

| Layer | Present | Missing or limited |
|-------|---------|--------------------|
| Input and provenance | PDB (Cα, optional Cβ), multi-model PDB, DCD reader/writer, RNA PDB (`topofold-io`, `topofold-core::pdb`) | No XTC/TRR reader (most public GROMACS data is XTC). Chain breaks are treated as bonds (m1). No input hashing or residue-mapping manifest in the library. The Python readers load whole trajectories into memory |
| Representations | Cα κ, τ, local and total writhe, Cβ angle θβ (backbone-only), RNA P-trace κ, τ, RNA base vector | No φ/ψ, χ, inter-domain geometry, contacts or cavity descriptors. RNA base vector cannot measure χ (M4) |
| Comparison and search | (κ, τ) discrete Fréchet and DTW, writhe-spectrum distance, VP-tree, exhaustive subcurve search | Full-frame cascade is approximate; its recall is measured only on AdK development data (LE-1) (M6). No residue-aligned κ–τ distance (A) or window RMSD in the library (the V2b window-RMSD search exists only as harness code in `benchmarks/audit/b2opt/`). The writhe-spectrum distance is not a metric for unequal lengths (M6). The Python `discrete_frechet_distance` is Cartesian and not SE(3)-invariant (M5). Fréchet weights are fixed at 1.0 in queries |
| Analysis | Window BC scan with circular means (M2), Gaussian MI network, IDP density, RNA hinge scan | BC is not a test (M3). MI estimator undocumented and untested for significance (M7). "Spectral" density has no ablation (m4) |
| Outputs | Results JSON convention (AGENTS.md rule 3), audit tooling for M3 | No CLI, no report generator, no library-level provenance record |

## 4. Capability matrix

Definitions:
- **Implemented:** code exists in the library.
- **Numerically tested:** regression tests against independent references, invariances or edge
  cases.
- **Scientifically validated:** shown on independent data, against independent labels, with
  equal-information baselines, for a stated question. *(Status 2026-10-12: no **library** capability is in this
  column. The WP-D held-out result concerns two harness methods for geometric local-conformation retrieval within
  mdCATH only: A, a residue-aligned distance on library κ/τ, and B2, window Cα RMSD, which is not a library
  function.)*
- **Proposed:** not implemented.

| Capability | Impl. | Numerically tested | Validated | Notes |
|------------|:-----:|--------------------|:---------:|-------|
| Cα κ, τ | ✔ | SE(3) invariance (random rigid motions, tolerance 1e-6 rad); subcurve locality (1e-12); golden-dataset tests | ✘ | τ near ±π needs circular handling |
| Writhe (total, local) | ✔ | 8 regression tests vs independent polygonal and numerical Gauss references; reflection and reversal (M1, unreleased) | ✘ | Absolute values changed ×(−2) vs v0.8.3; thresholds need recalibration |
| θβ (Cβ orientation) | ✔ | Unit tests | ✘ | Backbone descriptor only; cannot see χ1 (C4). Glycine pseudo-Cβ inconsistent |
| RNA P-trace κ, τ / base vector | ✔ | Unit tests | ✘ | Base vector blind to χ (M4) |
| Circular mean / centered deviations | ✔ | Unit tests incl. branch cut (M2, unreleased) | — | Characterization pending (WP-A) |
| Window BC scan | ✔ | Formula unit tests | ✘ | Heuristic; no null model or FDR (M3) |
| Product-vM mixture ΔBIC detector | Python audit prototype only | M3-A2 audit: provisional O2 | ✘ | Not in the engine |
| Subcurve search (`query_subcurve`) | ✔ | One integration test; LE-1 Rust-vs-NumPy agreement on AdK | ✘ | Exhaustive linear scan, so exact for its distance by construction. *2026-10-12:* the legacy flexible (Fréchet) coupling was deprioritized for same-sequence windows (LE-1, AdK development; engineering decision; formal outcome inconclusive) |
| Full-frame search (`query`) | ✔ | One integration test | ✘ | Approximate cascade. *2026-10-12:* Recall@10 vs exhaustive measured only on AdK development data at the default radius (LE-1), with no speed-up observed there; unmeasured elsewhere (M6) |
| MI coupling network | ✔ | 2 unit tests | ✘ | Estimator inconsistent with docs; no significance; autocorrelation unhandled (M7) |
| IDP topological density | ✔ | 3 unit tests | ✘ | Naming and ablation (m4) |
| DCD reader/writer | ✔ | Python I/O tests (skipped when data absent) | — | Frame-by-frame reader in Rust; Python `read_dcd` collects all frames |
| Throughput | — | EB-1 (2026-10-12): harness implementations on development data, per stage, with correctness checks. EB-1 and Amendment 1: insufficient evidence. Amendment 2: decision engine (B2 single-query R2 projection). BQ-1: V2b harness variant correct and faster than the harness baseline on tested workloads; acceptance not met. BQ-2: conclusion C (insufficient evidence). BQ-3: stopped at memory preflight, no P2 measurement | — | Harness-only; not in the library. Earlier µs/frame figures are unverified (m9) |
| XTC reader, φ/ψ/χ, cavity features, CLI, reports | Proposed | — | — | See §6 and §8 |

## 5. Representation requirements

### 5.1 Phenomena and what can observe them

| Phenomenon | What changes | Observed by current representations? | Complementary features needed |
|------------|--------------|--------------------------------------|-------------------------------|
| **Local backbone deformation** | Virtual bond/torsion angles of a few consecutive residues | κ, τ, local writhe: yes, if the change alters the Cα trace. Changes that leave the Cα trace nearly unchanged (e.g. peptide-plane flips) may be weak or invisible | φ/ψ as (sin, cos) for full-backbone data |
| **Relative motion of internally rigid domains** | Orientation/position of domains; local geometry changes only at hinges | Local descriptors see the hinge residues, not the amplitude or direction of the motion. By construction they are blind to anything outside the window. Global writhe changes, but non-specifically | Inter-domain distances/orientations, domain decomposition |
| **Side-chain / base rearrangements** | χ angles, base orientation, pairing, stacking | Cα-only: no. θβ: no (C4). RNA base vector: not χ (M4) | χ angles (sin, cos), side-chain centroid vectors; RNA χ, sugar pucker, base-pair/stacking descriptors |
| **Cavity formation** | Solvent-accessible volume between atoms | No representation in the engine. Backbone change is neither necessary nor sufficient for a pocket | Cavity detection on all heavy atoms, pocket-lining residues, independent pocket ground truth |
| **Statistical coupling** | Co-variation of descriptors within an ensemble | MI network, with the M7 caveats | Correct estimator, significance and autocorrelation handling |
| **Causal allostery** | Effect of a perturbation at one site on another | Not observable from one unperturbed ensemble | Perturbation data (mutations, ligands) or perturbation simulations |

### 5.2 Testable hypotheses and decision criteria

None of these prescribes an architecture. Each will be decided by evidence.

| ID | Hypothesis | Test | Decision |
|----|------------|------|----------|
| H1 | The per-frame scalar circular mean of a window discards coordinate-wise state differences | WP-A fixture F1 | If confirmed: restrict the BC scan to documented use, or evaluate a per-coordinate embedding as a separate candidate with its own validation |
| H2 | Window-level outputs are numerically unreliable near zero resultant or antipodal cancellation | WP-A F3–F5 | If confirmed: carry R or flag affected windows; quantify the affected region |
| H3 | Local window descriptors are unaffected by flexible regions outside the window | Already true by construction; subcurve locality tested to 1e-12. Full-frame search with disordered termini is untested | Claim locality only for window descriptors; test full-frame search separately before any claim |
| H4 | Cα-level descriptors localize real backbone-state changes at least as well as equal-information baselines | WP-B, then a held-out system | Margin and metric fixed before results. No claim without held-out evidence |
| H5 | Side-chain-only changes are invisible to Cα descriptors | Constructed fixtures (future) | If confirmed, as expected from geometry: tracks 1 and 4 need χ features before claims |
| H6 | Full-frame cascade recall is acceptable for retrieval | Recall@k vs brute force on the same index | Below a pre-declared target: label it approximate or remove the coarse filter |

## 6. Application tracks

None of these tracks is solved, and no earlier result for them stands: the v0.8.3 benchmarks
were synthetic or interpolated (C1).

| Track | Required observables | Independent evidence (to be verified per AGENTS.md rule 9) | Baselines | Blockers that actually affect it |
|-------|----------------------|------------------------------------------------------------|-----------|----------------------------------|
| **1. Abl kinase conformational change** (DFG flip, αC position, A-loop) | Local backbone (DFG, A-loop), χ of Phe (DFG), αC rigid-body position | Experimental structures in defined states; published unbiased MD with state definitions | φ/ψ, χ, inter-residue distances (e.g. the αC salt bridge), local RMSD, TICA with lag grid | No χ (C4/H5), no inter-element distances, no real MD. Earlier PDB IDs (2GQG, 1IEP) must be re-verified before reuse |
| **2. Allosteric coupling networks** | Ensemble co-variation (statistical); perturbation response (causal) | Mutational/perturbation data for the same system | Linear and rank correlation, MI with a validated estimator, dihedral-based networks | M7; autocorrelation handling; no perturbation data |
| **3. α-synuclein / IDP ensembles** | Distributions of local geometry, transient secondary structure, compaction | Ensemble-averaged experiments (e.g. NMR chemical shifts, SAXS), with entries to be identified | Secondary-structure propensities, Rg, contact maps | m4; writhe scale changed (M1); no real IDP trajectories in the repo |
| **4. KRAS and Mpro conformational/pocket investigations** | Switch-region backbone change (local); pocket opening (cavity) | Experimental apo/holo structures; curated cryptic-pocket sets (e.g. Meller et al. 2023, to be verified) | Local RMSD/φψ for switches; fpocket/MDpocket and PocketMiner for pockets | No cavity representation; C4/H5; m8. Backbone localization can proceed before pocket claims |
| **5. RNA riboswitch conformational change** | P-trace geometry, χ, pucker, pairing/stacking | Experimental apo/holo structures, with entries to be identified; chemical probing | η/θ pseudo-torsions, χ, base-pair annotations | M4; no real RNA MD |

**Sequence.** Shared protein local-conformation analysis and retrieval comes first (§7). It
underpins tracks 1 and 4 (backbone parts) and partly track 3. Track 5 needs M4 first; track 2
needs M7 first; pocket claims in track 4 need cavity features. These can be planned in
parallel; none blocks the first deliverable.

## 7. First deliverable: validated protein local-conformation analysis and retrieval

**Scope**
- U1 (κ, τ, local writhe per window and frame);
- U4 restricted to exhaustive subcurve retrieval;
- U2/U3 restricted to localizing regions of local backbone-state change in proteins.

**Evidence required before calling it validated**
1. Numerical: existing invariance and locality tests; WP-A characterization of the angular
   components; a decision on the same-sequence retrieval distance (M5); and, if the full-frame
   cascade is included, Recall@k vs brute force (M6).
2. Development: WP-B on a development dataset (AdK is the candidate), with labels from
   experimental structures and equal-information baselines.
3. Held-out: at least one further protein, selected by the procedure in §9 before any TopoFold
   result on it is seen, with metrics and margins fixed in advance.

Pocket detection, kinetics, free energies and allostery are outside this deliverable.

## 8. AdK: conditional development role

*E. coli* adenylate kinase (open 4AKE / closed 1AKE; MD from Orädd & Andersson 2021, Zenodo
10.5281/zenodo.5583119) is a candidate **development** dataset for the first deliverable, not a
mission target. Its use is conditional on the gates in
[`REAL_MD_EVALUATION_FEASIBILITY.md`](REAL_MD_EVALUATION_FEASIBILITY.md):

- If neither trajectory changes state, that limits what AdK can test; it says nothing about
  the engine.
- AdK results alone cannot validate anything. *(Status 2026-10-12: a held-out system set has since been evaluated; see §9 WP-D status.)*

## 9. Execution plan

### WP-A: Angular characterization harness

Fully specified in
[`ANGULAR_REPRESENTATION_CHARACTERIZATION_PROTOCOL.md`](ANGULAR_REPRESENTATION_CHARACTERIZATION_PROTOCOL.md)
§7–§8.

- **Files:** fixture generator, trace builder, Rust example harness, Python/NumPy path,
  comparison script, harness tests, report. No change to `crates/*/src`.
- **Interfaces exercised:** `circular_mean`, `mean_resultant_length`, `centered_deviations`,
  `StreamingMoments`, `sarles_bimodality_coefficient`, `compute_bimodality_profile`,
  `extract_curve_invariants`, and the Python bindings `compute_bimodality_profile` and
  `compute_invariants`.
- **Completion:** all fixtures run on all paths; recovery and bitwise-replication checks pass;
  results JSON written; and a retain / restrict-use / change / insufficient-evidence decision
  recorded per protocol §8 (amendments: protocol §9).
- **Dependency:** protocol review only. Needs no data download.

**WP-A status (2026-10-04): executed.** Decision **CHANGE** for the scalar circular-mean/BC
scan component only. See `benchmarks/audit/angular_characterization/REPORT.md` and its editorial
log. Production behaviour is unchanged.

### WP-C: Ordered angular representation (candidate)

Protocol: [`ORDERED_ANGULAR_REPRESENTATION_PROTOCOL.md`](ORDERED_ANGULAR_REPRESENTATION_PROTOCOL.md).
It evaluates a per-coordinate, residue-ordered (cos τ, sin τ) representation with transparent
same-sequence distances on SYNTHETIC fixtures, reusing the immutable WP-A fixtures plus two
controls. Fidelity, retrieval and detection are kept separate; detection is out of scope.
- **Outcome:** retain / restrict / reject / insufficient evidence.
- **Dependency:** protocol review. No production change.
- WP-B does not depend on it.

**WP-B G0–G1 status (2026-10-04): executed.** The outcome is *unresolved pending specified
evidence* for global-state tasks (`benchmarks/audit/adk_feasibility/REPORT.md`). Next-experiment
decision memo: [`EVALUATION_DESIGN_DECISION_MEMO.md`](EVALUATION_DESIGN_DECISION_MEMO.md), which
recommends WP-D (local backbone retrieval against an independent φ/ψ target).

### WP-B: Real-MD data feasibility (AdK)

Fully specified in [`REAL_MD_EVALUATION_FEASIBILITY.md`](REAL_MD_EVALUATION_FEASIBILITY.md) §6–§7.

- Download ≈ 1.37 GB (4 files, md5-verified).
- Convert XTC to Cα/backbone arrays with a pinned external reader, recording versions and
  hashes.
- Assess states independently (RMSD to 4AKE/1AKE).
- Proceed/stop gates G0–G2. No TopoFold or baseline scoring in WP-B.
- **Dependency:** protocol review and explicit download authorization. Independent of WP-A.

**WP-D status (2026-10-12): held-out evaluation completed.** Protocol chain:
`WPD_LOCAL_RETRIEVAL_PROTOCOL.md` → `WPD_HELDOUT_MDCATH_PROTOCOL.md` → Amendments 1–3 → `WPD_HELDOUT_MDCATH_EXECUTION_CORRECTIONS.md`;
results `benchmarks/audit/wpd_local_retrieval/REPORT_HELDOUT.md`; interpretation corrections
`WPD_HELDOUT_MDCATH_CLOSEOUT_ADDENDUM_2026-10-12.md`.
- Within mdCATH (30 superfamilies, 320/348 K), A (residue-aligned κ–τ) met the declared criterion for the angular
  φ/ψ target and B2 (window Cα RMSD) for the Cartesian peptide-plane target. Both are retained as candidates for their
  respective geometric tasks. This addresses §7 evidence item 3 for **geometric local-conformation retrieval only**.
- B2 at W = 8 on φ/ψ is inconclusive against the margin; W = 5 results do not establish equivalence; C1 failing the
  competitive criterion does not mean no gain over random (addendum C1–C3).
- Not covered: the legacy Fréchet index and writhe cascade (deprioritized for this task, LE-1), the bimodality scan and
  every detector, U2/U3 localization of state changes, biological function, and any clinical or oncology use.
- EB-1 executed 2026-10-12 (`EB1_ENGINEERING_BENCHMARK_PROTOCOL.md`; `benchmarks/audit/eb1/REPORT.md`): correctness passed;
  frozen decision **insufficient evidence** (one memory projection unsupported because of a specification error in
  EB-1). Amendment 1 (`EB1_AMENDMENT_1.md`; `REPORT_A1.md`): accounting corrected, B2 models confirmed; still
  **insufficient evidence** (3 cells failed by a frame-count design error; A batch-time fit residual above limit).
  Amendment 2 (`EB1_AMENDMENT_2.md`; `REPORT_A2.md`): valid sizes, direct-read loaders, released repeated outputs, a
  replacement A batch-time model confirmed on fresh workloads; frozen decision **engine**, driven only by the projected
  B2 single-query time at P2 (R2_B2). The loader and A are within requirements as projected. This is an extrapolated
  development-data result; an optimization experiment is proposed, not implemented
  (`EB1_B2_SINGLE_QUERY_OPTIMIZATION_PROPOSAL.md`). BQ-1 (harness experiment, frozen protocol plus Amendment 1):
  R2_B2 not supported by the tested NumPy variants; see `benchmarks/audit/b2opt/REPORT.md`. BQ-2
  (`B2_SCALING_CONFIRMATION_PROTOCOL.md`; `benchmarks/audit/bq2/REPORT.md`): unchanged V2b; query-time and memory models
  confirmed; preparation model unconfirmed on the largest SYNTHETIC stress cells; **conclusion C (insufficient
  evidence)**. BQ-3 (`B2_DIRECT_P2_QUALIFICATION_PREFLIGHT.md`): stopped at its memory preflight; **no direct P2
  measurement**. These are harness-only engineering results; they imply no biological capability, novelty or production
  readiness. Oncology evidence requirements: `ONCOLOGY_VALIDATION_REQUIREMENTS.md` (no target chosen).

**Held-out selection, specified before any AdK scoring**
1. Write the eligibility criteria into a frozen file: open licence, unbiased MD, experimental
   two-state structures for the same sequence, no TopoFold output used.
2. Enumerate candidates from literature and databases by those criteria, recording sources.
3. Choose among eligible candidates by a recorded rule (e.g. a seeded random draw). No TopoFold
   or baseline run happens on any candidate before the choice is frozen and hashed.

### What can proceed now

*Status 2026-10-13 (replaces the 2026-10-12 table; that table's content is preserved in the git history):*

| Done | Next (bounded) | Subsequent priorities | Needs more |
|------|----------------|-----------------------|------------|
| WP-A, WP-C, WP-B (G0–G1), WP-D AdK tier, LE-1, WP-D held-out evaluation and closeout; EB-1 with Amendments 1–2; BQ-1 (acceptance not met; V2b correct and faster on tested workloads); BQ-2 (conclusion C); BQ-3 (stopped at memory preflight) | Direct P2 qualification of unchanged V2b, on a host with ≥ about 3.2 GB available memory; protocol and mandatory-verification acceptance rule frozen before execution | Strong equal-information competitive comparisons (external tools and strongest reasonable baselines); biological validation per application. Library integration of any harness variant only after a separate review | Any oncology task needs the evidence in `ONCOLOGY_VALIDATION_REQUIREMENTS.md`; track-specific work needs the blockers in §6 |

Historical issues block only the capabilities mapped to them in §6 and §10. They do not block
everything.

## 10. Blocker map

| Issue | Capability affected | Tracks affected |
|-------|---------------------|-----------------|
| C1 | Credibility of all v0.8.3 benchmark claims; no code capability | All (prior evidence void) |
| C2 | Evaluation design | All evaluations |
| C3 | PMF estimator, removed from scope | None |
| C4 | θβ as side-chain descriptor | 1, 4 |
| C5 | PROTAC benchmark, withdrawn | None |
| M1 (resolved, unreleased) | Writhe values; thresholds | 3, full-frame search |
| M2 (resolved, unreleased) | Circular means in BC scan | Screening (U2) |
| M3 | BC scan as a statistical test | Screening (U2), not retrieval |
| M4 | RNA χ | 5 |
| M5 | Same-sequence distance; Cartesian Fréchet API | Retrieval (U4, U5) |
| M6 | Full-frame search exactness | Full-frame retrieval only |
| M7 | MI network | 2 |
| M8 | Reporting discipline | All reports |
| m1 | Chain breaks in real inputs | All real data with gaps |
| m2, m3, m9 | Precision/speed/zero-copy wording | Documentation and performance claims |
| m4 | IDP density | 3 |

## 11. Commercial hypotheses and targets (not established capabilities)

The project has commercial ambitions. They are recorded here as hypotheses and targets, each
with the measurement that would test it. None is currently established. They are not
contributor rules.

| Target / hypothesis | Current state | What would test it |
|---------------------|---------------|--------------------|
| Blind identification of experimentally known pockets from third-party data, delivered as an automated report | No cavity representation; no report generator; no CLI | Cavity features plus a held-out pocket benchmark with independent ground truth; an implemented report pipeline that reads only results JSON |
| Robustness to disordered termini (local analysis unaffected by tail motion) | True by construction for window descriptors; untested for full-frame search and for any pocket claim | WP-A/WP-B plus tail-perturbation tests of full-frame search |
| Throughput of 50,000–100,000 frames/s on a standard server | Unmeasured for the library. Python path copies float32 to f64; trajectories are loaded fully into memory. EB-1/BQ-1/BQ-2 measured harness-only single-query search, not library frame throughput | Benchmark against vectorized MDTraj/MDAnalysis on identical tasks and hardware, with the configuration recorded |
| Analyses usable without per-dataset tuning | Not currently true: window size, BC threshold, coarse radius, writhe tolerance, gap and regularizer are parameters | Documented defaults plus sensitivity analysis across them |
| Automated pipeline removing manual transcription of results (M8) | Results-JSON convention exists for audits only | Report generator that reads results files and has tests |

Competitive positioning against MDTraj, PyEMMA, MDAnalysis or commercial tools may be claimed
only after equal-information comparisons (AGENTS.md rule 6).
