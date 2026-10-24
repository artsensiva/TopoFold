# TopoFold Engine: Mission, Capabilities and Roadmap

**Status:** DRAFT for review (2026-10-04; WP-D status reconciled 2026-10-12; engineering status reconciled 2026-10-13;
development programme added 2026-10-22, §12)
**Baseline:** `295d739` (branch `scientific-rebuild`); programme reconciled against `7ff4b90`

This document states what TopoFold is for, what it can do today and the evidence for that,
and what must happen next. Where it conflicts with older material (ADRs, whitepaper, v0.8.3
benchmark descriptions), this document and [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) take precedence.
Rules for contributors stay in [`../AGENTS.md`](../AGENTS.md).
Recorded decisions and superseded claims are listed in [`DECISION_LOG.md`](DECISION_LOG.md). This is the only
programme document; other files link here instead of restating status.

---

## 1. Mission

> Build TopoFold into a reproducible geometry and search engine for biomolecular conformational
> ensembles: helping researchers locate, compare, retrieve, and investigate structural changes
> across protein and RNA trajectories, with traceable evidence and application-specific
> validation.

**Long-term purpose.** Contribute to cancer research and, ultimately, treatment, with appropriate commercial returns
and scientific recognition. This purpose sets direction; it is not a capability or a near-term claim.

**Near-term mission (2026-10-22).**

> Develop reproducible methods for locating, comparing and retrieving molecular conformational changes, and establish
> whether they improve oncology-relevant research decisions beyond established tools.

**Evidence levels.** Every result states the level it supports. Passing one level does not establish the next.

| Level | Meaning | Results at this level so far (observation) |
|---|---|---|
| A | Mathematical and numerical correctness | Writhe and circular-mean fixes with reference tests (M1, M2); WP-A characterization; V2b numerical correctness against SVD/QCP references (BQ-1, CMP-1) |
| B | Structural/geometric task performance | WP-D held-out, within mdCATH only: A on the φ/ψ target, B2 on the peptide-plane target (scoped geometric retrieval). LE-1: AdK development tier, inconclusive. C1, if run, would be level B (expert structural labels) |
| C | Biological relevance | None |
| D | Practical value to researchers | None. EB-1, BQ-1–3 and CMP-1 are harness engineering measurements, not user value |
| E | Therapeutic or clinical benefit | None, and outside the near-term scope (`ONCOLOGY_VALIDATION_REQUIREMENTS.md` §2) |

Not validated at any level as detectors: the scalar BC scan (M3) and the Product-vM diagnostic (M3-A2, provisional O2).

Kept outside levels B–E (2026-10-23):
- **Development evidence on real MD, descriptive:** DEV-CA-1, DEV-ENS-1, DEV-REP-1.
- **SYNTHETIC calibration of replica-level comparison rules:** REPCHAR-1, PERM-1, NULLSTRESS-1.

None of these is biological validation, and none validates the BC scan.

Distinctions kept in every report: harness implementation vs library capability; synthetic control vs real MD vs
experimental measurement; structural label vs independent functional evidence; association vs causation; engineering
improvement vs methodological novelty.

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
| Outputs | Results JSON convention (AGENTS.md rule 3); Python CLIs for the Cα adapter, κ/τ pipeline and descriptive group comparison (`python -m topofold_inputs`, `.pipeline`, `.ensemble`, `.structpair`), which write tables, hashed manifests and generated reports | No CLI for the Rust core; no CI (m11); provenance records exist for the Python workflows, not for the core library |

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
| Circular mean / centered deviations | ✔ | Unit tests incl. branch cut (M2, unreleased) | — | Characterized in WP-A (SYNTHETIC, 2026-10-04): decision CHANGE for the scalar circular-mean → BC path (information loss and numerical sensitivity shown); no change adopted |
| Window BC scan | ✔ | Formula unit tests | ✘ | Heuristic; no null model or FDR (M3). WP-A: low BC is not absence of multimodality, high BC is not two states. The replica permutation work (PERM-1) does not apply to or validate BC |
| Product-vM mixture ΔBIC detector | Python audit prototype only | M3-A2 audit: provisional O2 | ✘ | Not in the engine |
| Subcurve search (`query_subcurve`) | ✔ | One integration test; LE-1 Rust-vs-NumPy agreement on AdK | ✘ | Exhaustive linear scan, so exact for its distance by construction. *2026-10-12:* the legacy flexible (Fréchet) coupling was deprioritized for same-sequence windows (LE-1, AdK development; engineering decision; formal outcome inconclusive) |
| Full-frame search (`query`) | ✔ | One integration test | ✘ | Approximate cascade. *2026-10-12:* Recall@10 vs exhaustive measured only on AdK development data at the default radius (LE-1), with no speed-up observed there; unmeasured elsewhere (M6) |
| MI coupling network | ✔ | 2 unit tests | ✘ | Estimator inconsistent with docs; no significance; autocorrelation unhandled (M7). *2026-10-22:* double regularization fixed (ε now added once; regularized I_ε, analytical regression tests); module description corrected (no copula transform). Remaining M7 items open |
| IDP topological density | ✔ | 3 unit tests | ✘ | Naming and ablation (m4) |
| DCD reader/writer | ✔ | Python I/O tests (skipped when data absent) | — | Frame-by-frame reader in Rust; Python `read_dcd` collects all frames |
| Throughput | — | EB-1 (2026-10-12): harness implementations on development data, per stage, with correctness checks. EB-1 and Amendment 1: insufficient evidence. Amendment 2: decision engine (B2 single-query R2 projection). BQ-1: V2b harness variant correct and faster than the harness baseline on tested workloads; acceptance not met. BQ-2: conclusion C (insufficient evidence). BQ-3: first attempt stopped at memory preflight; the direct run on SYNTHETIC P2-dimension rows qualified prepared-database operation, while preparation + first query (excluding input loading) did not meet 10 s. CMP-1 (bounded, one thread): no decisive advantage for V2b over MDTraj 1.11.2 / MDAnalysis 2.10.0 under the frozen criteria | — | Harness-only; not in the library. Earlier µs/frame figures are unverified (m9) |
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

Current status, next task and trigger for each track, plus PROTAC/molecular glue: §12.3.

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
  evidence)**. BQ-3: the first attempt stopped at its memory preflight (`B2_DIRECT_P2_QUALIFICATION_PREFLIGHT.md`). After
  a memory recheck, the direct run (`B2_DIRECT_P2_QUALIFICATION_PROTOCOL.md`; closeout
  `B2_DIRECT_P2_QUALIFICATION_CLOSEOUT_2026-10-13.md`) on SYNTHETIC, performance-only P2-dimension rows **qualified
  prepared-database operation**. **Preparation plus the first query did not meet 10 s**; that span excludes input
  loading and row construction. BQ-2's conclusion C stands. CMP-1 (`CMP1_COMPETITIVE_EVALUATION_PROTOCOL.md`; `benchmarks/audit/cmp1/REPORT.md`): bounded
  comparison with MDTraj 1.11.2 and MDAnalysis 2.10.0 on development and SYNTHETIC workloads. Frozen conclusion: **no
  decisive advantage**, so V2b integration is not supported. MDTraj was faster where it met the common F32 accuracy
  track; V2b was faster than MDAnalysis per-pair QCP elsewhere. These are harness-only engineering results; they imply no biological capability, novelty or production
  readiness. Oncology evidence requirements: `ONCOLOGY_VALIDATION_REQUIREMENTS.md` (no target chosen).

**Held-out selection, specified before any AdK scoring**
1. Write the eligibility criteria into a frozen file: open licence, unbiased MD, experimental
   two-state structures for the same sequence, no TopoFold output used.
2. Enumerate candidates from literature and databases by those criteria, recording sources.
3. Choose among eligible candidates by a recorded rule (e.g. a seeded random draw). No TopoFold
   or baseline run happens on any candidate before the choice is frozen and hashed.

### What can proceed now

*Status 2026-10-22:* the queue is in §12.7. Earlier versions of this table (2026-10-12, 2026-10-13) are preserved in
the git history.

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
| Blind identification of experimentally known pockets from third-party data, delivered as an automated report | No cavity representation; no pocket report (the existing generated reports cover descriptive ensemble comparison only); no CLI | Cavity features plus a held-out pocket benchmark with independent ground truth; an implemented report pipeline that reads only results JSON |
| Robustness to disordered termini (local analysis unaffected by tail motion) | True by construction for window descriptors; untested for full-frame search and for any pocket claim | WP-A/WP-B plus tail-perturbation tests of full-frame search |
| Throughput of 50,000–100,000 frames/s on a standard server | Unmeasured for the library. Python path copies float32 to f64; trajectories are loaded fully into memory. EB-1/BQ-1/BQ-2 measured harness-only single-query search, not library frame throughput | Benchmark against vectorized MDTraj/MDAnalysis on identical tasks and hardware, with the configuration recorded |
| Analyses usable without per-dataset tuning | Not currently true: window size, BC threshold, coarse radius, writhe tolerance, gap and regularizer are parameters | Documented defaults plus sensitivity analysis across them |
| Automated pipeline removing manual transcription of results (M8) | Results-JSON convention for audits; the descriptive ensemble command generates its report from its tables | Report generator that reads results files and has tests |

Competitive positioning against MDTraj, PyEMMA, MDAnalysis or commercial tools may be claimed
only after equal-information comparisons (AGENTS.md rule 6).

## 12. Development programme (2026-10-22)

Planning and evidence reconciliation only: nothing in this section is implemented, run or authorized beyond what it
says. Statement labels: **[O]** observation (in this repository), **[D]** deduction, **[L]** literature fact,
**[A]** assumption, **[P]** proposal. **[L?]** marks a literature statement taken from the external analyst assessment
of 2026-10-15 (an untracked file at the repository root) that has not yet been verified under AGENTS.md rule 9; it
must not be cited as fact until verified. That assessment was treated as evidence to assess, not as instructions.

### 12.1 Position

- [O] Level B, scoped: WP-D supports geometric local-conformation retrieval within mdCATH only (A on φ/ψ, B2 on the
  peptide-plane target). Nothing reaches level C, D or E.
- [O] Engineering chain closed: BQ-3 qualified prepared-database operation on SYNTHETIC P2-dimension rows; preparation
  plus first query did not meet 10 s. CMP-1 found no decisive advantage; V2b integration is unsupported and V2b stays
  an audit artifact.
- [O] KRAS C1 is unfrozen and paused. The Rascore terms request is drafted
  (`C1_RASCORE_TERMS_REQUEST_DRAFT.md`); the repository records neither sending nor a reply.
- [O] The MI estimator regularized the joint covariance twice; with covariances of order ε this returned 0 for
  correlated features. Fixed 2026-10-22 (D16). It is an engineering correction; it does not validate allostery.
- [D] The scalar BC scan and Product-vM diagnostics are not validated state detectors (M3).
- [L?] PENSA is an established tool for comparing ensembles and is the expected strongest baseline for ensemble work;
  it was not considered in earlier planning.

### 12.2 Phases

Each phase passes only on its own verification. "Inconclusive" never counts as passing.

| Phase | Entry | Deliverables | Verification and exit | Failure / inconclusive | Parallel / deferred |
|---|---|---|---|---|---|
| **0 Consolidate** | — | This section; `DECISION_LOG.md`; stale status wording fixed | Links resolve; frozen artifacts unchanged; no duplicated status. Exit: owner accepts | Correct the documents; no other work waits on it | Parallel: none. Deferred: all experiments |
| **1 Minimum foundation for the next task** | Phase 0 accepted; one defect authorized at a time | Fixes mapped to the capability they block (table below), each with a regression test against an independent reference | `cargo test --workspace`, `pytest tests/`; analytical reference passes. Exit: the defects blocking the next scheduled evaluation are closed | Fix fails its reference: revert, document, keep the capability unvalidated | Parallel: one lightweight feasibility activity. Deferred: defects that block no scheduled evaluation |
| **2 First bounded application (KRAS C1)** | Written permission covering internal evaluation and publication; design frozen with the §12.5 checklist; adapter-only foundation | Frozen protocol, smallest adapter, results JSON, report | Outcome from the frozen table (`C1_EVALUATION_DESIGN_DRAFT.md` §5). Exit: any of its four outcomes, reported at level B | Permission denied or unresolved: C1 stays **blocked**, not closed, and can resume when permission is clarified. An alternative evaluation runs only after passing its own data-access, feasibility, protocol and budget gates (§12.5). Materially inferior or inconclusive: reported as such; no re-scoring on the same data | Parallel: Phase 1 engineering, Phase 3 feasibility. Deferred: other RAS loops, nucleotide analysis |
| **3 Reproducible ensemble changes** | One justified system chosen by a recorded rule, with verified data access and licence, a feasibility note, a frozen protocol and a budget; PENSA-class baseline available | Feasibility note, then a frozen protocol comparing two ensembles of one system against established workflows | Replicate-level uncertainty with temporal dependence handled; baselines with equal information. Exit: a level-B result on ensemble differences, or a recorded negative | No accessible data with usable terms: stop and record; do not substitute generated ensembles as MD | Any justified system qualifies; ABL is one candidate, assessed without assuming data access or commercial usability. Deferred: kinetics, adaptive sampling |
| **4 One research workflow** | A level-B result from Phase 2 or 3; at least a few user conversations recorded | One workflow (API/CLI, report generated from results JSON) and a measurement of end-to-end user value against the user's current tools | Pre-declared task-time or decision-quality endpoint (level D). Exit: measured value or a recorded negative | Value not shown: keep the method as a library function; no product claim | Deferred: multi-target platform, commercial packaging |
| **5 Selective expansion** | A level-D result, or a level-B result with a recorded reason to extend | One new application or shared path per decision, from §12.3 | As for the corresponding phase | As for the corresponding phase | Everything not selected stays deferred |

**Phase 1 defect map** (the fix is scheduled only when the capability it blocks is needed):

| Defect | Capability it blocks | Needed by |
|---|---|---|
| M7 MI estimator (double regularization fixed 2026-10-22); no significance; Gaussian and linear-angle assumptions | MI coupling network | Allostery path only |
| m1 chain breaks treated as bonds | Real structures with gaps | Phase 2 excludes chains with missing Cα, so not blocking; Phase 3 needs it |
| No real-structure adapter (mmCIF, residue mapping, altloc) | Any real-data evaluation | Phase 2 (C1 adapter, after freeze), Phase 3 (via MDTraj/MDAnalysis readers) |
| A and B2 not in the library | Library-level retrieval | Phase 4; Phase 2 can use the harness A with reference tests and an established RMSD library |
| M5 Cartesian Fréchet naming and units; M6 approximate cascade | Full-frame and Fréchet retrieval | Not needed by Phases 2–3 |
| M3 BC scan not a test | State detection | Conditional state-detection path only |
| M4 RNA χ | RNA path | Phase 5 only |
| m12 MSRV | Reproducible builds | Before any release; small |

**Dependencies**

```
Phase 0 ──► Phase 1 (only defects needed next)
              │
              ├─► Phase 2 KRAS C1 ◄── permission + frozen design + adapter
              │        │
              └─► Phase 3 ensemble changes ◄── system choice + data access + PENSA-class baseline
                       │
        Phase 2 or 3 level-B result ──► Phase 4 workflow ◄── user research
                                              │
                                              ▼
                                        Phase 5 selective expansion
```

Phases 2 and 3 do not depend on each other. Phase 3 does not wait for C1, and starting it does not close C1.

### 12.3 Portfolio

Status: **active**, **ready**, **blocked**, **deferred**, **closed**. Budgets are provisional ranges assuming one
contributor with the current workstation, excluding waiting on third parties; they are not commitments.

| ID and question | Capability and evidence | Limitations and dependencies | Data and permissions | Strongest baseline | Smallest next task | Effort / compute / storage | Accept / negative / inconclusive | Status | Trigger |
|---|---|---|---|---|---|---|---|---|---|
| **KR** Does SW1 Cα geometry recover Rascore ON/OFF labels? | Harness A; level B only if run | Cα sees part of the label rule; confounds (`C1_DECISION_NOTE_2026-10-13.md`) | PDB CC0; Rascore label terms unresolved | B2 via established library; B1; Y32–G12 proxy | Freeze design after permission; adapter | 3–8 days / < 2 CPU-h / < 1 GB | Frozen outcome table | blocked | Written permission |
| **AB** Do kinase state changes (DFG, αC, A-loop) appear in ensemble comparisons? (one candidate system for EN) | KIN-DESC-1 (D48, four public Abl chains, static): side-chain descriptors (D1/D2, A-loop H-bond) reproduce the stated labels; Cα κ/τ and Cα distances separate the A-loop contrast in many windows but DFG only at a few four-point windows; αC not assessable (all αC-in) | Needs χ and inter-element distances (C4, H5) | [L?] published ABL trajectories reported available on request for non-commercial use; unverified | PENSA-class comparison; φ/ψ, χ, distances, local RMSD | Metadata-only feasibility: data access, terms, state definitions | 1–3 days / none / none | Accessible data with usable terms / none / unclear terms | ready (feasibility only) | Lightweight slot free |
| **MP-b** Mpro backbone conformational change | Same as other local backbone retrieval (level B scoped, mdCATH only) | m8 (chain description); no Mpro data selected | To be identified | Local RMSD, φ/ψ, PENSA-class comparison | None | — | — | deferred | Selected as a Phase 3 system |
| **MP-p** Mpro pocket opening | None | Requires cavity observables, which the engine lacks (§5.1) | To be identified | fpocket/MDpocket; [L?] PocketMiner | None until cavity features are scheduled | — | — | deferred | Cavity features scheduled |
| **AL** Statistical coupling, then perturbation response | MI network implemented, not validated | M7; association ≠ causation; needs perturbation data for causal claims | None needed for the fix | Linear/rank correlation; validated MI estimator | None scheduled; regularization fix done (D16). Significance and estimator choice remain open | — | — | deferred | Allostery selected as an application |
| **ID** IDP / α-synuclein ensemble distributions | IDP density implemented, not validated | m4; writhe scale (M1) | Real ensembles and ensemble-averaged experiments to be identified | Secondary-structure propensities, Rg, contact maps | None | — | — | deferred | Phase 3 method exists |
| **RN** RNA riboswitch changes | P-trace κ, τ | M4 χ blindness | To be identified | η/θ, χ, base-pair annotations | None | — | — | deferred | M4 fixed and a data source found |
| **PG** PROTAC / molecular glue | None; earlier benchmark withdrawn (C5) | Ternary-complex geometry outside current representations | To be identified | To be defined | None | — | — | deferred | A defined geometric question with data |
| **RP** Representation and feature adequacy | κ, τ, writhe at level A | Cα-only identifiability limit; [D] writhe adds information only relative to an incomplete representation | — | φ/ψ, χ as (sin, cos), distances | Physics statement per feature when an application needs it (AGENTS.md rule 7) | — | — | deferred | Phase 3 system chosen |
| **RT** Exhaustive retrieval; approximate indexing only if needed | A, B2 at level B (scoped); V2b harness | CMP-1: no decisive advantage | — | MDTraj, MDAnalysis | None; reuse established RMSD | — | — | closed (optimization); ANN deferred | Measured bottleneck in a real workflow |
| **IO** Adapters; streaming only if needed | Cα pipeline (`python/topofold_inputs`): MDAnalysis adapter → core κ/τ with validity reasons → NPZ/JSON/report; engineering-ready (SYNTHETIC tests; DEV-CA-1 development run on local AdK). Not scientifically validated | Cα only; frames held in memory; frame identity unverifiable for coordinate-only formats; no mmCIF/SIFTS mapping | Per dataset | MDTraj/MDAnalysis readers | Backbone (N, Cα, C) only when a selected comparison needs it (package README lists activation triggers) | — | — | ready (engineering) | Selected comparison needing backbone atoms |
| **EN** Calibrated ensemble-distribution comparison on any justified system | Descriptive per-window JS-distance comparison (`python/topofold_inputs/compare.py`); DEV-ENS-1 development run on local AdK (one run per condition); DEV-REP-1 replicated comparison (mdCATH 1ha8A00; no window met the complete-separation criterion); replica-level methods characterised on SYNTHETIC data (REPCHAR-1, PERM-1: exact replica permutation near nominal under the null, more sensitive, does not resolve confounding; NULLSTRESS-1: it also rejects when only temporal correlation differs). Real-data permutation design on `1ha8A00` deferred, not scored (permutation not validated for isolating marginal geometry). Dependence-aware stationary-bootstrap Kuiper prototype failed SYNTHETIC calibration (DABOOT-1, null rejection 0.095–0.230); descriptive comparisons only | Temporal dependence; replicates | Per system | [L?] PENSA; [L?] MMD with block or wild bootstrap | Done 2026-10-22 (`EN_FEASIBILITY_2026-10-22.md`); next: owner sends the ABL request (`ABL_DESRES_DATA_REQUEST_DRAFT.md`) | — | — | blocked (provider permission) | DESRES reply covering the intended uses and a matched run inventory |
| **FS** Frame selection for expensive downstream analysis | None | [A] gain depends on cost ratio of screening to analysis and the fraction kept | — | Random and clustering-based selection | None | — | — | deferred | Phase 4 workflow needs it |
| **SD** Conditional state detection | BC scan, Product-vM: not validated | M3; no calibrated null | — | Dip test or mixture ΔBIC with calibrated null | None | — | — | deferred | A selected application needs discrete states; ensemble differences alone are not sufficient |
| **KI** Adaptive sampling and kinetics | None | [D] structural retrieval does not determine kinetics | — | MSM workflows | None | — | — | deferred (long-term) | Level-B ensemble result and a kinetic question |
| **UX** API/CLI, reporting, user research | Python CLIs and generated reports (adapter, pipeline, `ensemble` v0.3.0, `structpair` v0.1.0); local read-only viewer for structpair results (`apps/structpair_viewer`); pilot handoff; separate-agent review | No Rust-core CLI; m11 no CI; no independent human trial yet (one external AI technical trial and one AI confirmation completed, D45–D46) | — | Users' current scripts and tools | A few structured user conversations (no claims) | 1–3 days / none / none | Recorded needs / no interest / mixed | ready (lightweight) | Lightweight slot free |

Short-lived Git branches are used for each authorized task.

### 12.4 Concurrency and cost

- At most one active scientific evaluation, one bounded engineering task and one lightweight feasibility or
  user-research activity.
- **Essential:** Phase 0; the C1 gates; the defect blocking the next evaluation. **Optional:** UX
  conversations, m12. **Deferred:** everything marked deferred in §12.3.
- Effort is given as provisional ranges with assumptions; no monetary figures, dates or success probabilities.

### 12.5 Verification rules

**Pre-measurement checklist** (frozen and hashed before any scoring): question; estimand; endpoint; permitted inputs;
provenance; leakage controls; development/evaluation separation; baselines with equal information; correctness
references and tolerances; statistical unit, weighting and uncertainty; temporal dependence; multiple testing; margin
justified separately from sample size; precision assumptions; resource limits and stop rules; permitted conclusions.

**Records:** software versions and file hashes; reports generated from recorded results only.

**Outcome categories** (exactly one per endpoint): implementation failure; data or design failure; materially
inferior; useful non-inferiority; superior; insufficient evidence. Non-inferiority to a weak method does not establish
adequacy. Ensemble work includes PENSA and established workflows as baselines; RMSD uses established implementations.

### 12.6 Adaptation rules

At each decision gate:
1. Compare the result with the frozen criteria.
2. Classify it with the outcome categories above.
3. List the smallest next actions.
4. Compare their information value, cost, dependencies and risk.
5. Recommend proceed, restrict, switch, defer or close, and record it in `DECISION_LOG.md`.

Prohibited: changing thresholds after results; relabelling failures as success; erasing negative results; changing
population or endpoint silently; reusing exposed evaluation data as fresh held-out evidence; repeated amendments
without a named uncertainty. A follow-up needs a concrete question and a bounded budget. Changes to the claim, data,
acceptance criteria, major dependencies or resource ceiling need a new recorded decision. Routine choices inside an
authorized package are made without escalation.

Milestone updates state: where we are; what happened; what it means; what comes next; any decision needed. Decision
tables use Option | Expected benefit | Cost/effort | Main risk | Reversibility | Recommendation, with at most three
alternatives.

### 12.7 Queue

| Slot | Task | Status | Condition |
|---|---|---|---|
| Scientific evaluation | **KR** (C1): freeze the design, then the smallest adapter, then scoring | blocked | Written permission from the Rascore maintainers |
| Bounded engineering (parallel) | Done: **AL** MI fix (D16); **A** Cα pipeline (D21); descriptive group-comparison command `python -m topofold_inputs.ensemble` (D34), corrected to v0.2.0 (D36). Calibration phase closed (D33): inference reopens only under the recorded conditions. Next: practical user pilot of the descriptive tool (`docs/PILOT_HANDOFF_ENSEMBLE_TOOL.md`; internal walkthrough done, independent user run pending; D37); **B** backbone atoms only when a selected comparison requires them | — | Selected comparison or user need |
| Biological route (parallel) | **ABL** preferred (D35): run inventory and permissions need the owner to send the D. E. Shaw Research request; processed HDX uptake was not found in the inspected resources and the authors' processing pipeline is unavailable locally (raw PRIDE data are public); requests to D. E. Shaw Research and the HDX authors drafted, not sent. Other routes resumable (C1, GPCRmd, mdCATH method work) | blocked (owner actions) | Data and permissions received |

**Application paths: next actions (2026-10-23).** Each path waits on its own dependency; at most one runs as the
active scientific evaluation at a time.

| Path | Next action | Depends on (capability / data) | Status |
|---|---|---|---|
| ABL | Preregistration drafted, not frozen (`docs/ABL_WT_M472I_PREREGISTRATION_DRAFT.md`, D47): no apo M472I MD is described in the paper, so the branch (apo or imatinib-bound) is decided at freeze from the inventory; input schemas, validator and peptide-to-c-Abl mapping ready | Trajectories, run inventory and permissions (D. E. Shaw Research request); processed HDX uptake (request to the HDX authors) | waiting for replies: both requests sent by the owner; apo M472I availability unknown |
| KRAS / Mpro | Static-structure milestone closed (D42; KRAS-STRUCT-1/-2). Reusable `structpair` command delivered (D43). No new KRAS case now; possible later: a GDP/GMPPNP pair in one crystal form, only if a specific question needs it | Existing Cα geometry and RMSD components (available); mmCIF-only entries need a reader (not installed); pocket claims need cavity observables (not available) | ready (KRAS backbone); Mpro: metadata search done (D49): F@h COVID Moonshot set is the one conditionally suitable candidate (CC0; ligand-site and IC50 endpoints, not conformational); pocket work deferred |
| Allostery | None scheduled; MI estimator correction done (D16) | Perturbation data for causal claims; validated significance (closed, D33) | deferred |
| IDP / α-synuclein | None scheduled | Real IDP ensembles with ensemble-averaged experiments; m4 | deferred |
| RNA riboswitches | None scheduled | RNA χ descriptor (M4); RNA trajectories | deferred |
| Outside-user trial | AI technical trial and confirmation closed (D45–D46); viewer work stopped. Human trial kept as a future product-usability milestone: send the invitation and feedback form; record it with the form's §0 (unaided vs assisted), with consent before publishing or identifying the tester | Owner nominates a computational scientist | deferred (future milestone; owner action) |
