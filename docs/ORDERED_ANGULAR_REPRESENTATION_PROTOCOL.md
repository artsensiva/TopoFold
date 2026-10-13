# Ordered Angular Representation Protocol (WP-C, candidate)

**Status:** Amended (A–E, §11) before execution; executed as WP-C. See
`benchmarks/audit/ordered_representation/REPORT.md`. No production change.
**Baseline:** `295d739` plus the WP-A package (`benchmarks/audit/angular_characterization/`).
**Scope:** SYNTHETIC evaluation of one candidate representation. It does not choose a
multimodality detector, does not set production thresholds, and makes no claim about modes,
states or biology.

## 1. Question

Does keeping each torsion at its residue position remove the representation collisions and
branch sensitivities demonstrated in WP-A, while giving useful distances for comparing local
conformations of the same sequence?

Three questions are answered separately and must not be merged:

| Question | What it asks | Evidence |
|----------|--------------|----------|
| **A. Representation fidelity** | Does the representation distinguish what WP-A showed the scalar summary confuses, and does it vary continuously across the ±π branch? | Analytical identities and WP-A fixtures |
| **B. Retrieval quality** | Do same-sequence distances retrieve frames of the same constructed state, compared with equal-information baselines? | Exhaustive nearest neighbours on labelled synthetic fixtures |
| **C. Statistical detection** | Can modes or states be detected? | **Out of scope.** No detector, mode count or BC is computed on the embedding |

## 2. Representation

**Window.** Residues (Cα indices) i … i+W−1, as in the existing scan.

**Torsions.** τ_j for j = i … i+W−4 (m = W−3 values). τ_j is the dihedral of Cα atoms
j, j+1, j+2, j+3, in the convention of production `compute_torsions` (the IUPAC sign; WP-A
recovery checks established agreement to 4.4e-16 rad). τ_j is attributed to the virtual bond
(j+1, j+2). The feature order is ascending j, so position within the window is explicit.

**Ordered torsion features (dimensionless):**

  x = (cos τ_i, sin τ_i, cos τ_{i+1}, sin τ_{i+1}, …, cos τ_{i+m−1}, sin τ_{i+m−1}) ∈ ℝ^{2m}

**Curvature block (radians), kept separate:** k = (κ_i, …, κ_{i+W−3}), with W−2 values; κ_j is
the turning angle at Cα j+1. κ ∈ [0, π] is not periodic and needs no embedding.

**No combined torsion-plus-curvature distance is defined.** Combining them needs a relative
weight, which would be an unvalidated parameter. Each block is reported on its own.

### Same-sequence distances (residue to residue, no alignment or reparameterization)

For two frames a, b of the same window, let d_j = |wrap(τ_j^a − τ_j^b)| ∈ [0, π] be the circular
distance of torsion j.

- **Chord distance (primary):** D_ch(a, b) = √( (1/m) Σ_j ‖u(τ_j^a) − u(τ_j^b)‖² ) = ‖x_a − x_b‖ / √m,
  with u(θ) = (cos θ, sin θ). Since ‖u(α) − u(β)‖ = 2|sin(d/2)|,
  D_ch = √( (1/m) Σ_j 4 sin²(d_j/2) ) ∈ [0, 2].
- **Circular RMS distance (transparent comparator):** D_circ(a, b) = √( (1/m) Σ_j d_j² ) ∈ [0, π].
- **Relation:** for each torsion, (2/π)·d ≤ 2 sin(d/2) ≤ d, so (2/π)·D_circ ≤ D_ch ≤ D_circ.
  For small differences, 2 sin(d/2) = d − d³/24 + …, so the two coincide to first order. D_ch is a
  Euclidean distance in ℝ^{2m}; D_circ is the L2 combination of circle geodesics. Both are
  metrics.
- **Normalization across window sizes:** the 1/m factor puts both distances on a common
  per-torsion RMS scale. It does **not** make different window sizes equally difficult, and it
  does not make a distance value or threshold transferable between window sizes.
- **Curvature block distance:** D_κ = √( (1/(W−2)) Σ (κ^a − κ^b)² ), in radians.

This same-sequence correspondence differs from the Fréchet distance used by the index (M5):
there is no reparameterization, so residues cannot be misregistered.

### What is preserved

- Each torsion modulo 2π, together with its position in the window.
- Chirality: the sign of τ, through the sine components.
- Continuity: each feature is 1-Lipschitz in its torsion (|Δx_j| ≤ |Δτ_j|), so there is no
  branch cut and no dependence on a centering direction.
- Invariance under proper rigid motions, because τ is.

### What is still omitted

- Virtual bond lengths, so chain breaks are invisible (m1).
- κ, unless the separate block is used.
- Everything outside the window, including relations between windows and domain motions.
- Side chains and bases; cavities.
- Torsions at near-collinear triples are undefined; production raises an error there. Such
  windows are flagged, not imputed.

### Expected differences from other representations

- **Common shift.** The representation is not invariant to adding a common α to every torsion,
  because that is a different conformation. Only paired invariance is expected:
  D(τ_a + α, τ_b + α) = D(τ_a, τ_b).
- **Reflection.** Reflection maps τ → −τ, which the representation detects. The Cα-distance
  baseline does not.

### Comparator and diagnostics

- **Historical comparator:** the existing per-frame scalar circular mean, with frame distance
  |wrap(m_a − m_b)|, and the WP-A scan outputs. These are reported unchanged.
- **Diagnostics only:** the within-window resultant length per frame is recorded. It is not part
  of any distance and cannot restore positional information.

## 3. Equal-information baseline

A local Cα-distance vector computed from the **same window atoms**: all Cα–Cα distances (p, q)
with q − p ≥ 2, inside i … i+W−1. Pairs with q − p = 1 are bonds and are constant here. The
baseline distance is the RMS difference over pairs, in Å.

Information differences, stated in advance:
- The baseline is reflection-invariant.
- Its (p, p+2) distances encode κ.
- Its (p, p+3) distances encode cos τ, so they are blind to the sign of τ.
- The candidate and the baseline use identical frames, residues and coordinates.

## 4. Fixtures

### Reused, immutable WP-A fixtures

Inputs are read through the frozen WP-A manifest `fixtures/SYNTHETIC_fixtures.json`, and
every file is hash-verified before use. If WP-A's `build/` files are absent, the inputs are
regenerated **in memory** with the WP-A generator functions and must reproduce the manifest
hashes. Nothing is ever written into the WP-A directory, and the WP-A pipeline is not rerun.

| Family | Cases | Used for |
|--------|------:|----------|
| F1 | 6 | A: collisions; B: retrieval with stored state labels |
| F2 | 24 | A: paired-shift invariance (pairs of frames under the same α) |
| F3 | 27 | A: continuity at near-zero within-window resultant |
| F4, α = 0 only | 44 | B: retrieval. Labels are the φ₂ assignment, reproduced from the recorded seed and draw order and verified by reproducing the torsion hash. A: behaviour near ±π |
| F6b | 4 | A: continuity sanity with equal torsions |

### New controls (only these two, each justified)

| ID | Construction | Justification |
|----|--------------|---------------|
| **N1: paired ±ε** | W = 5, F = 1,000, p = 0.49, σ = 1°, ε ∈ {±1e-6, ±1e-3, ±1e-1}: 6 cases. **One** permutation (minority assignment) and **one** per-frame noise draw, from `default_rng(20261104)`, shared by all six. Torsions as in WP-A F4 | WP-A F4 drew independently per configuration, so effects of ε sign could not be separated from draw variation. N1 isolates ε for both the candidate and the historical comparator |
| **N2: mirror** | F1a_s0 with every τ negated: 1 case, deterministic | Distinguishes torsion features (chirality-sensitive) from the Cα-distance baseline (reflection-invariant) |

Coordinates for N1/N2 are built with the unchanged WP-A `build_traces.py`.

### Torsion sources (no new production or Rust code)

| Source | What it is |
|--------|------------|
| T0 | Prescribed fixture torsions |
| T1 | Production torsions through the freshly built binding `compute_invariants`, with float32 input (the user-facing path) |
| T2 | Independent NumPy dihedral of the f64 coordinates |

T1 vs T0 measures the effect of the float32 path; T2 vs T0 is the recovery check.

**Historical scalar comparator.** For every retrieval fixture, the per-frame circular means
and scan BCs come from the unchanged WP-A Rust example (production functions, R1/R2). It is run
on a WP-C manifest that points at the hash-verified WP-A fixture files and at the N1/N2 files,
and it writes only into the WP-C `build/`. The WP-A REF formula is evaluated alongside as a
check.

## 5. Analytical checks (independent of the implementation under test)

### 5a. Feature-contract checks (convention)

Distances are preserved by any consistent permutation of positions, by swapping cos and sin, and
by reversing every sine sign. Distance checks therefore cannot detect those errors, so the
feature values are checked directly:

- **FC1, values.** The asymmetric known torsions τ = (0.1, 0.7, −2.0, 3.0, −0.4) rad (W = 8,
  window start i = 0) must give exactly
  x = (cos 0.1, sin 0.1, cos 0.7, sin 0.7, cos −2.0, sin −2.0, cos 3.0, sin 3.0, cos −0.4, sin −0.4).
  Expected entries are computed independently with `math.cos` and `math.sin`; tolerance 1e-15.
- **FC2, window offset.** From a 12-residue chain with torsions τ_j = 0.1·(j+1) for j = 0…8, the
  window starting at residue i = 3 with W = 6 uses τ_3, τ_4, τ_5. Its features are
  (cos 0.4, sin 0.4, cos 0.5, sin 0.5, cos 0.6, sin 0.6).
- **FC3, residue mapping.** The same window reports its torsions on the virtual bonds
  (4, 5), (5, 6), (6, 7).
- **FC4, chirality.** Feature sine entries carry the sign of τ, so for τ = −0.7 the sine entry is
  −sin 0.7.

### 5b. Distance-geometry checks

1. **Chord identity** ‖u(α) − u(β)‖ = 2|sin(d/2)| on random pairs, including d near 0 and π
   (1e-12).
2. **Bounds** (2/π)·D_circ ≤ D_ch ≤ D_circ on every computed pair (slack 1e-12).
3. **F1 noiseless, state A vs state B:**

   | Case | D_ch | D_circ | Scalar comparator |
   |------|------|--------|-------------------|
   | F1a_s0 | 1 | π/3 | 0 |
   | F1b_s0 | √(2/5) | (π/3)·√(2/5) | 0 |
   | F1c_s0 | 1 | π/3 | π/3 |

   Tolerance 1e-12.
4. **Paired shift.** D(τ_a + α, τ_b + α) = D(τ_a, τ_b) on F2 frame pairs (1e-12).
5. **Continuity on F3.** Perturbing one torsion by δ changes x by at most |δ| (slack 1e-15),
   reported next to the WP-A scalar flips. For noiseless F3 frames, the frames at +ε and −ε
   differ only in τ₂ (by 2ε), so D_ch = √2·|sin ε|. The WP-A scalar means instead differ by
   about π.
6. **N2.** D_ch(A, mirror A) = √((1/m) Σ 4 sin² τ_j). The baseline distance between them is 0
   (≤ 1e-9 Å, f64 coordinates).
7. **Branch.** τ = π − 1e-12 vs τ = −π + 1e-12 gives D_ch ≤ 3e-12.
8. **F6b.** With equal torsions x in both positions, D_ch between two frames is 2|sin(Δx/2)|.

### 5c. Precision rules

| Inputs | Rule |
|--------|------|
| T0 (prescribed, f64) | Exact analytical tolerances as stated: 1e-12 for distances, 1e-15 for feature entries |
| T2 (NumPy recovery from f64 coordinates) and T1 (production, float32 input) | Recovery-aware. Let δ be the measured maximum circular torsion-recovery error of the case. Each feature entry must lie within δ + 1e-15 of its T0 value (1-Lipschitz), and each distance within 2δ + 1e-12 of its T0 value. T1 retrieval scores are reported separately from T0, never pooled |

## 6. Retrieval (question B)

- For every query frame, the k nearest other frames of the same fixture are found by
  **exhaustive** search, with k ∈ {1, 10}.
- Scored separately for D_ch, D_circ, the historical scalar distance, the Cα-distance baseline
  and D_κ.
- The score is the fraction of neighbours sharing the query's construction label. Labels are
  independent of every representation.
- **Ties, label-independent.**
  - The query itself is excluded by index; other frames at distance 0 remain candidates.
  - Let d_k be the k-th smallest distance. Candidates within tol = 1e-12·(1 + d_k) of d_k form
    the boundary tie set S₌; candidates below d_k − tol form S₍<₎. The tolerance is fixed:
    rounding of O(1) distances is ~1e-15.
  - The score is the **expected** purity under uniform selection of the k − |S₍<₎| remaining
    neighbours from S₌:
    purity = ( same(S₍<₎) + (k − |S₍<₎|)·same(S₌)/|S₌| ) / k.
  - Labels are never used to break ties; frame index order is never used.
  - Tie frequency is reported: the fraction of queries where |S₌| > k − |S₍<₎|.
- **Reference:** the chance purity of each fixture (same-label share among the other frames).
- **Descriptive only.** Scores are reported per fixture and per representation, for these fixed
  fixtures. No confidence intervals, significance claims or population inference. No
  replication study.
- Retrieval fixtures:
  - F1 (6);
  - F4 α = 0 (44);
  - N1 (6);
  - N2 pooled with F1a_s0, with 4 labels (state × chirality).

## 7. Corruption-test set (bounded)

Each corruption must violate an observable contract and be detected. Tests use isolated
temporary copies or substituted functions.

| # | Corruption | Contract violated | Detected by |
|---|------------|-------------------|-------------|
| 1 | Torsion positions permuted in the feature builder | Positional order | FC1 (FC2) |
| 2 | cos/sin swapped | Entry convention | FC1 |
| 3 | Sine sign reversed | Chirality convention | FC1 and FC4 |
| 4 | Window offset shifted by one residue | Residue mapping | FC2, FC3 |
| 5 | Tie handling breaks ties by frame index or label | Label-independent ties | Tie unit test with a constructed fully tied set (expected purity must equal chance) |
| 6 | Fixture byte altered | Input hash | Hash verification |
| 7 | A stored result value altered | Reproducibility | `run.py --verify` (recomputes and compares) |
| 8 | A case missing from the results | Coverage | `run.py --verify` |

## 8. Execution scope and completion criteria

**Files**

| Path | Content |
|------|---------|
| `benchmarks/audit/ordered_representation/` | `features.py` (candidate, distances, feature contracts), `inputs.py` (WP-A input verification, label reproduction), `controls.py` (N1/N2), `retrieval.py` (tie-aware purity), `run.py` (execution, results, `--verify`), `make_report.py`, `test_ordered.py`, `README.md`, `REPORT.md`, `results/SYNTHETIC_ordered_representation_*.json` |
| its `build/` | N1/N2 binaries and comparator outputs. Covered by its own narrow ignore rule |

No change to `crates/*/src`, the WP-A package or its results.

**Complete when**
- every case in §4 has run on T0/T1/T2;
- all §5 checks and §7 tests pass;
- retrieval scores and tie frequencies are written (descriptive, no intervals);
- the report quotes only numbers from the results JSON;
- `cargo test --workspace` and `pytest tests/` still pass;
- the decision record below is filled in for review.

## 9. Decision record (filled in after execution, then reviewed)

Three separate verdicts are recorded, then one decision follows from them by the ordered rules
below.

**V1. Harness integrity.** Pass when coverage, input hashes, recovery, tie-handling tests and the
corruption set all behave as specified.

**V2. Fidelity.** Pass when every §5a feature-contract check and every §5b distance check passes
on T0, and every §5c recovery-aware check passes on T2 and T1.

**V3. Fixed-fixture retrieval** (T0, k = 10, descriptive). For each retrieval fixture, compare
D_ch purity with:
- the historical scalar comparator;
- the Cα-distance baseline.

A **material deficit** is a D_ch purity more than 0.05 below a comparator's. The 0.05 is a
pre-declared descriptive margin for these fixtures; it is not a production threshold. Ties and
small differences are recorded as such. They are neither failures nor evidence of superiority.

**Decision**, applied in order; the first rule that matches decides:

| Order | Decision | Rule |
|---|----------|------|
| 1 | **Insufficient evidence** | V1 fails |
| 2 | **Reject** | V2 fails |
| 3 | **Restrict** | V3 shows a material deficit against the historical scalar comparator on any fixture. A replacement must not lose what it replaces. Name the affected fixtures |
| 4 | **Retain: eligible for a real-MD pilot** | Otherwise |

Deficits against the Cα-distance baseline do not block pilot eligibility: a working
representation need not beat every baseline. They are reported, and the pilot must include
that baseline. No superiority or generalization is claimed from these fixtures.

The decision concerns representation and retrieval on synthetic fixtures only. It says nothing
about detecting modes or biological states.

## 10. Connection to a later real-MD pilot (not part of this task)

If retained, the candidate can enter a separately reviewed pilot on the WP-B development data,
after WP-B's gates pass. Equal information applies throughout:
- the same Cα coordinates, frames and windows (W from a pre-registered grid);
- representations: ordered features (D_ch, D_circ), the κ block, the local Cα-distance vector,
  and the historical scan as comparator.

| Task | What it measures | Precondition |
|------|------------------|--------------|
| Retrieval | Query windows taken from the experimental reference structures; retrieve MD frames; score against the independent G1 frame states (Cα RMSD to 4AKE/1AKE) | Meaningful only if G1 shows both states |
| Localization | One per-window between-ensemble separation statistic, defined once and applied identically to every representation, scored against crystal-derived residue labels (L1/L2) | Positive-unlabelled metrics and a time-block bootstrap |

**WP-B does not depend on this protocol.** Its data-integrity and state-assessment gates proceed
whether the candidate is retained or not.

## 11. Amendment record (before execution)

| ID | Amendment | Reason |
|----|-----------|--------|
| A | Direct feature-contract checks FC1–FC4, kept separate from distance-geometry checks; each corruption mapped to the contract it violates | Distances are invariant under consistent permutation, cos/sin swap and sine reversal, so distance checks could not detect those errors |
| B | Label-independent tie handling with expected purity at the k-th boundary; query excluded by index; tie frequency reported | Avoids index-order or label leakage at ties |
| C | Descriptive retrieval scores only; bootstrap removed | Leave-one-out query scores are not independent observations, so an ordinary bootstrap over queries is not justified; no population inference from fixed fixtures |
| D | Separate verdicts V1–V3 and ordered, non-conflicting decision rules; deficits against the Cα baseline do not block a pilot | Interval-based rules conflicted and treated ties as failure |
| E | Exact tolerances for T0; recovery-aware bounds for T2/T1; RMS normalization described as a common scale only | Precision differs by input path; normalization does not imply transferability |
| — | Historical comparator computed by the unchanged WP-A Rust example writing into the WP-C `build/` | Production comparator values without touching WP-A outputs |
