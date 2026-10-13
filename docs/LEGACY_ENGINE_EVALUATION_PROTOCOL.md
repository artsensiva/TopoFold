# Legacy κ–τ Engine Evaluation: Development Protocol (LE-1)

**Written:** 2026-10-05, before any new AdK outcome in this protocol was computed. Frozen by SHA-256 in
`benchmarks/audit/legacy_engine/results/le_protocol_freeze.json`.

**Scope:** bounded development evaluation on already-authorized AdK data. **Not** held-out validation.

**Earlier protocols preserved unchanged:**
- WP-A: `ANGULAR_REPRESENTATION_CHARACTERIZATION_PROTOCOL.md`;
- WP-C: `ORDERED_ANGULAR_REPRESENTATION_PROTOCOL.md`;
- WP-D AdK: `WPD_LOCAL_RETRIEVAL_PROTOCOL.md`;
- WP-D held-out: `WPD_HELDOUT_MDCATH_PROTOCOL.md` and `WPD_HELDOUT_MDCATH_AMENDMENT_1.md`.

This protocol changes none of them, and it does **not** change the held-out method set or domain selection.

## 1. Questions (kept separate)

| Question | What this protocol can address | What it cannot |
|---|---|---|
| Representation fidelity | Which Cα geometric elements each feature uses or omits (§3, synthetic tests) | Information content for arbitrary proteins |
| Distance suitability | Residue-aligned vs flexible (Fréchet) coupling for same-sequence windows | Cross-sequence or indel matching |
| Retrieval utility | Agreement of each method's neighbours with a φ/ψ target and a peptide-plane target, on AdK | Generalization to other proteins (held-out tier) |
| Computational utility | Measured runtime and memory of the implementations used here | Optimized-implementation comparisons beyond those listed |
| Biological relevance | None directly. Both targets are geometric reconstructions from the same trajectory | Function, states, kinetics, allostery, pockets |
| Novelty | Literature and software comparison (§9) | — |

The task is a **geometric local-conformation retrieval benchmark**. No inference about multimodality,
metastability, kinetics, allostery or pockets is made from it.

## 2. Dataset roles and contamination history

| Data | Role | Prior influence |
|---|---|---|
| SYNTHETIC fixtures (constructed in `test_legacy.py`; WP-A fixtures by hash) | Mechanism and implementation tests only | WP-A/WP-C designs |
| AdK, Orädd et al. 2021 (Zenodo 10.5281/zenodo.5583119): `open_start`, `closed_start` | Development | WP-B G0/G1 (inspected); WP-D AdK tier: M1, M2, B1, B2, B3 (= K here) and C1 were **already scored** against the φ/ψ target (`REPORT_ADK.md`). Known pattern: B2 best in every cell; M1/M2 behind B1/B2. **New in this protocol:** methods A, F, F_lock, G_RMSD, G_F, and the peptide-plane target |
| mdCATH (frozen selection, not downloaded) | Held-out | Metadata only. Not used here |

Because the AdK φ/ψ results for M1, M2, K, B1 and B2 are known, any contrast involving only those methods is
**not** a new test. Contrasts involving A or F are new for both targets.

## 3. Geometric elements (Cα trace r_0 … r_{N−1}, 0-based)

**Elements:**
- Bond b_j = r_{j+1} − r_j.
- Turning angle (production "curvature") κ_v = angle between b_{v−1} and b_v, at vertex v = 1…N−2. Range
  [0, π]. Production array index v−1.
- Torsion τ_j = dihedral(r_j, r_{j+1}, r_{j+2}, r_{j+3}), about bond b_{j+1}. Range (−π, π]. Production array
  index j = 0…N−4.

**A window of W residues (start s = i − 1)** contains:
- interior vertices s+1 … s+W−2, giving **W − 2 κ**;
- dihedrals s … s+W−4, giving **W − 3 τ**.

**Legacy pairing (verified in code: `distance.rs`, `discrete_frechet_invariants`):**
- Element j of a window couples (κ at vertex s+1+j, τ about bond b_{s+1+j}), for j = 0…W−4. The sequence
  length is min(#κ, #τ) = W − 3.
- **The last interior κ (vertex s+W−2) is omitted.**
- Each τ_j depends on both vertices s+1+j and s+2+j but is paired with only the first.
- Bond lengths are not used.

## 4. Methods

**Common input.** Every method receives the same whole-made Cα coordinates: float32 values held in float64,
loaded as in WP-D (`run_adk.load_traj`), with the same windows, query frames, candidate frames and exclusions.

**Notation:**
- d_κ = |κ^a − κ^b|;
- d_τ = |wrap(τ^a − τ^b)| ∈ [0, π];
- m = W − 3 and M = W − 2.

**Methods:**

| ID | Distance | Status |
|---|---|---|
| M1 | τ chord-RMS: √((1/m) Σ_j 4 sin²(d_τ,j / 2)) | Existing (WP-C/WP-D) |
| M2 | τ circular-RMS: √((1/m) Σ_j d_τ,j²) | Existing |
| K | κ-only: √((1/M) Σ_v d_κ,v²) | Existing (= WP-D B3) |
| F | **Legacy equal-weight κ–τ discrete Fréchet**, `weight_kappa = weight_tau = 1`, over the W − 3 legacy pairs. Pointwise distance √(d_κ² + d_τ²); minimax coupling. Computed by the **actual Rust path**: `ConformationalIndex(window_size=8, coarse_radius=0.5).fit(database frames)` then `query_subcurve(s, s+W−1, query_frame, k = n_D)` (exhaustive) | Legacy |
| F_lock | Lock-step version of F: max_j √(d_κ,j² + d_τ,j²) over the same W − 3 legacy pairs (diagonal coupling only) | Attribution diagnostic |
| A | **Residue-aligned equal-weight κ–τ**: √(½·(1/M) Σ_v d_κ,v² + ½·(1/m) Σ_j d_τ,j²), over **all** W − 2 κ and W − 3 τ | New combination of existing channels |
| B1 | Local Cα–Cα distances (q − p ≥ 2), RMS difference (Å) | Existing |
| B2 | Window Cα RMSD, Horn quaternion (proper rotations) | Existing |
| G_RMSD | Global whole-chain Cα RMSD (214 atoms, Horn); the same ranking for every window | Global control |
| G_F | Legacy full-chain κ–τ Fréchet, exhaustive: `query_subcurve(0, 213, …)`; the same ranking for every window | Global legacy control |

**Channel normalization in A.**
- Each channel is a mean over its own element count, then the two channels are averaged with weight ½. A
  window's τ and κ therefore contribute equally whatever their counts.
- **Equal weighting in rad² is a modeling choice**, not an absence of hyperparameters. κ varies over a
  narrower range than τ in proteins, so equal weights do not make the channels equally influential.
- **No weighting search** is performed. A further variant is allowed only if §6 tests reveal a specific defect.
  It would be declared as an amendment before computation, with every variant reported.

**Excluded:**
- C1 (scalar circular mean): already characterized in WP-A and WP-D.
- Legacy DTW: same coupling family as F; not part of the engine's query path.
- Writhe-spectrum distance as a retrieval method: approximate filter only.

**Writhe cascade (computational check only).** The full-frame `ConformationalIndex.query(q, k=10)` (default
radius 0.5) is compared with exhaustive G_F. Reported as recall@10 against G_F's tie-aware top 10, using the
same final distance. It is not used as a retrieval method.

## 5. Targets, queries, candidates and endpoints

**Splits and windows** (WP-D §5, unchanged):
- query = first ⌊0.45 n⌋ frames; gap ⌊0.10 n⌋ excluded; database = the rest;
- trajectories never pooled; W ∈ {5, 8}, all start positions.
- The gap reduces temporal adjacency; it does not make query and database frames independent.

**Primary target T_φψ** (WP-D §3, unchanged): circular RMS over interior φ and ψ (rad).

**Complementary target T_pp** (prespecified):
- Atoms: C_r, O_r and N_{r+1} for r = i … i+W−2, i.e. the **peptide-plane atoms of the window's W − 1 peptide
  bonds, excluding every Cα** (3(W − 1) atoms).
- Distance: their optimal-superposition RMSD in Å (Horn).
- O is made whole relative to its own C by minimum image; changes are counted.
- **Bias statement:** T_φψ is angular and may favour angular representations. T_pp is Cartesian and may
  favour Cartesian superposition (B2). Neither target replaces the other, and both are reported for every
  method.

**Endpoints:**
- **Primary:** normalized regret NR@10 (WP-D §6: regret / (random − oracle); 0 = oracle, 1 = random). It is
  computed per (window, query) with the WP-D tie rule (tol = 1e-12·(1 + d_k)). Primary target T_φψ.
- **Secondary:**
  - NR@10 under T_pp;
  - NR@1 and overlap@10 against the random reference k/n_D, both targets;
  - raw regret@10 (rad or Å).

## 6. Verification before real-data scoring (SYNTHETIC)

**`test_legacy.py` must pass.** It covers:
- independent rigid motion of each compared structure;
- torsion wrap at ±π;
- κ-only and τ-only perturbations;
- a scalar-mean collision;
- reflection (τ flips sign; κ and B1 unchanged; B2 distinguishes);
- boundary κ omission and unequal κ/τ lengths;
- repeated, collinear and near-collinear points;
- float32 vs float64;
- Fréchet-induced residue misregistration;
- empty and invalid inputs.

**Legacy Rust path vs the research implementation:**
- On AdK, the vectorized NumPy Fréchet (`legacy_lib.frechet_legacy`) must match the Rust `query_subcurve`
  distances for **every** window and query within 1e-9. Both are float64 from identical values;
  differences can only come from libm rounding (≈1e-15).
- A larger difference is reported and the **Rust values remain the primary F**.
- A Rust probe (`rust_probe/`, separate cargo workspace, path dependency on `topofold-core`, offline) records
  the engine's behaviour on empty, NaN and unequal-length inputs. Defects are documented, not repaired.

## 7. Statistical analysis

**Cell** = trajectory × W × target. Trajectories and W are never pooled.

**Summaries:**
- Per (window, query): NR@10 for each method.
- Cell statistic: the mean over windows of the per-window mean over queries (primary). The median across
  windows of window means (the WP-D convention) is also reported.

**Paired contrasts.** ΔNR(X − Y) = difference of cell means on identical (window, query) units. Pre-declared:

| Hypothesis | Contrast(s) | Supported if | Falsified if |
|---|---|---|---|
| H1: κ adds to τ | A − M2, A − K | Both intervals below 0 | Either interval above 0 |
| H2: residue alignment avoids harmful misregistration | F − A > 0; mechanism split F − F_lock (warping) and F_lock − A (max aggregation + omitted κ) | F − A interval above 0 | F − A interval below 0 |
| H3: Fréchet harms same-sequence retrieval | as H2 (same-sequence case only). "Helps other tasks" is **not testable** here and stays untested | | |
| H4: historical advantages reflect local vs global information | G_F − F, G_RMSD − B2, and F vs B2 | Global methods worse than their local counterparts, and F not better than B2 | Global legacy G_F ≥ local methods |
| H5: SE(3) invariance is shared | B1 and B2 are invariant by construction (tested in §6). Empirical: B1 or B2 − A ≤ 0 | — | — |
| H6: the target favours angular methods | Method ranks under T_φψ vs T_pp | — | — |

An interval that contains 0 is reported as **unresolved**. These are descriptive interval statements on one
protein, not significance tests; no multiplicity adjustment is made.

**Uncertainty: two-way block bootstrap.**
- Query frames are resampled in contiguous blocks of 100 frames (20 ns), and window start positions in
  contiguous blocks of 20. The two dimensions are resampled independently with replacement.
- The statistic is recomputed from block-sum matrices. The **same resample is used for every method**, so
  pairing is preserved.
- B = 2,000 resamples, seed 20261012, 95 % percentile intervals.
- **Sensitivity:** query blocks of 250 frames.
- **Diagnostic:** integrated autocorrelation time (WP-B `adk_lib.tau_int`) of the per-query series of window-
  averaged ΔNR, for the H1 and H2 contrasts.
- **What it estimates:** variability due to which query frames and which residue positions of *this* AdK
  trajectory are sampled.
- **What it cannot estimate:** database-frame variability (the candidate pool is fixed), between-trajectory,
  between-protein or force-field variability.
- If the τ_int diagnostic exceeds the 100-frame block length, the 250-frame intervals govern the H1/H2
  statements. That rule is fixed in advance.

**Practical usefulness.** There is **no task-based minimum useful gain** for AdK geometric retrieval, so no
binary superiority claim is made. Effects are reported as continuous paired differences with intervals, plus
absolute NR and the gap to random. The held-out acceptable-loss margin ρ = 0.10 (relative loss L, Amendment 1
estimand) is used only for the descriptive development classification below. It is not a usefulness
threshold.

**Low-signal and invalid cases:**
- Queries with an undefined NR denominator (≤ 1e-12) are excluded and counted, per method and target.
- A cell with best-method gain 1 − NR_best ≤ 1e-6 makes L not assessable (Amendment 1 §2).
- Non-finite distances are counted and ranked last (+∞). Any non-finite legacy F on AdK is reported as a
  defect.
- No window or query is dropped for being unfavourable.

## 8. Development classification (for A and F only; governs nothing in the held-out protocol)

**Loss statistic.** Per primary cell (trajectory × W, T_φψ):

L_X = (NR_X − NR_best) / (1 − NR_best)

- Cell means are used. best = argmin over {M1, M2, K, A, F, B1, B2}, re-selected per resample.
- Intervals use the conservative invalid-resample rule of Amendment 1 §2.
- Each cell is classed **within margin** (upper bound ≤ ρ), **beyond margin** (lower bound > ρ) or
  **inconclusive**.

**Outcomes:**

| Outcome | Rule |
|---|---|
| **Insufficient evidence** | §6 verification fails, a legacy-path failure, or any cell not assessable |
| **Retain as held-out candidate** | Within margin in all four primary cells. Adding the method to the held-out set still requires a reviewed amendment **before** held-out scoring |
| **Restrict** | Within margin in some cells and beyond in others; the role is limited to the within-margin cells |
| **Reject for the local same-sequence task** | Beyond margin in all four primary cells |
| **Inconclusive** | Otherwise |

T_pp results are reported beside the classification and do not change it.

## 9. Computation (separate from retrieval)

**Timing setup:**
- Measured on a fixed subset (open_start, W = 8, windows i = 1…20, all queries), single process.
- Steady state is timed separately from setup: Rust wheel build (from the WP-A provenance), index `fit`, and
  feature extraction.

**Implementations timed:**
- NumPy vectorized: M1, M2, K, A, F_lock, B1, B2, and NumPy Fréchet;
- the Rust legacy path, per `query_subcurve` call, which includes per-call query invariant extraction and
  per-frame subcurve allocation;
- the writhe cascade `query`.

**Memory:** peak resident memory per implementation, measured in separate subprocesses.

**Attribution:** implementation differences (Rust vs NumPy) are **not** attributed to the representation.

## 10. Failure handling and budget

**Stop rules:**
- Any WP-D input-hash or binding-identity failure (reused `run_adk.verify_inputs`).
- A legacy `fit` error.
- A reference-check violation.

Each stops the run and is reported. No substitute data is synthesized.

**Budget:** ≤ 2 h CPU on 4 workers, ≤ 8 GB RAM, ≤ 2 GB disk under the audit `build/` directory (git-ignored).
No downloads except the literature metadata checks in §11.

## 11. Prior work

Citations are checked against Crossref by DOI or title and recorded with the returned metadata. Nothing is
cited unverified.

## 12. Amendments

Any change after the freeze is recorded in a separate amendment file, stating whether relevant results had
already been observed.
