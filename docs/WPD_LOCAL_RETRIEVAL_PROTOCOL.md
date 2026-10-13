# WP-D Protocol: Cα-only Local-Conformation Retrieval against a Full-Backbone Reference Target

**Status:** AdK development tier amended (§10) and frozen before computation; hashes in
`benchmarks/audit/wpd_local_retrieval/results/wpd_protocol_freeze.json`. The mdCATH held-out arm
(§7–§8) is **provisional, not frozen and not executed**. Original status line: implementation-ready
draft. Freezing
(hashes recorded) must precede any computation of targets, features or scores.
**Baseline:** HEAD `295d739` plus the WP-A, WP-C and WP-B packages, all preserved unchanged.
**Eligibility (done):** `benchmarks/audit/wpd_local_retrieval/results/wpd_adk_backbone_eligibility.json`.
It establishes atom availability (N, Cα, C once per residue 1–214), finite coordinates, and backbone
bond continuity (≤ 2.0 Å) in raw coordinates. It does **not** establish numerical validity of any
dihedral; the harness checks that (§10, A1).

## 1. Claim scope

**Claimed (if supported):** how well each Cα-only window representation retrieves, from a
separate time range of the same trajectory (AdK) or another replica of the same domain (mdCATH),
frames whose **full-backbone local conformation** is close to the query's.

**Not claimed:**
- the target is not experimental ground truth;
- no pockets, no state detection, no biological discovery, no kinetics;
- no statistical non-inferiority.

**Nature of the target.** It is a structural reference computed from backbone N, Cα and C atoms.
No compared method receives N or C. The target is nevertheless **geometrically coupled** to the
Cα trace, which those atoms determine. The experiment measures agreement with a richer structural
description, not independent truth.

## 2. Inputs common to all methods

- **Coordinates:** Cα coordinates only, W consecutive residues per window, read from the
  trajectory in Å. AdK xtc data are quantized to 0.01 Å; float64 arithmetic does not restore
  precision.
- **Windows:** start residue i with residues i…i+W−1 (resid numbering 1–214 for AdK),
  W ∈ {5, 8}, all positions with 1 ≤ i ≤ 215 − W (210 windows for W = 5, 207 for W = 8). The two
  window sizes are analysed and reported separately.
- **Periodic boundaries:** Cα windows are made whole by sequential minimum image
  (`adk_lib.make_whole`, WP-B). For AdK the eligibility check found no backbone bond longer than
  2.0 Å in raw coordinates, so the reconstruction is expected to be the identity. That is
  verified per frame at run time; any change is recorded.

## 3. Target (full-backbone reference)

- **Dihedrals per residue r:**
  - φ_r = dihedral(C_{r−1}, N_r, Cα_r, C_r), defined for r = 2…214;
  - ψ_r = dihedral(N_r, Cα_r, C_r, N_{r+1}), defined for r = 1…213.
  - IUPAC convention (`build_traces.dihedral`, independently checked in WP-A), radians in
    (−π, π].
- **Target components:** all residues r = i+1 … i+W−2, the window's interior residues (M = W−2),
  with **both** φ_r and ψ_r.
  - The set is fixed by window geometry alone: every interior residue has both dihedrals
    defined from atoms inside the window's residue span.
  - It is not chosen for agreement.
  - ω and side chains are excluded, because the target is the standard Ramachandran description
    of local backbone conformation. ω is recorded as a diagnostic: count of |ω| < π/2,
    i.e. cis peptides.
- **Distance (periodic):**
  D_T(a, b) = √( (1/(2M)) Σ_r [ d(φ_r^a, φ_r^b)² + d(ψ_r^a, ψ_r^b)² ] ), with
  d(x, y) = |wrap(x − y)| ∈ [0, π].
- **Terminal exclusions:** windows lie within 1…214, so interior residues lie within 2…213.
  Every interior residue therefore has both dihedrals; no further terminal exclusion is needed.
- **Missing atoms:** any missing N, Cα or C makes every window containing that residue
  ineligible for all methods. AdK: none (eligibility check).
- **Chain breaks:** any N–Cα, Cα–C or C–N₊₁ bond > 2.0 Å after periodic reconstruction makes
  the windows spanning it ineligible in that frame, for all methods. AdK: none in raw coordinates.
- **Periodic boundaries for the target:** backbone atoms are made whole along the
  N–Cα–C–N₊₁ order by sequential minimum image before dihedrals are computed.

## 4. Methods (identical Cα windows)

All distances are computed per window between one query frame and one database frame. Each
method is used only to rank its own distances; no distances are mixed across methods.

| ID | Distance | Definition | Units |
|----|----------|------------|-------|
| M1 | Ordered chord D_ch | √((1/m) Σ_j ‖u(τ_j^a) − u(τ_j^b)‖²), u = (cos, sin), m = W − 3 Cα torsions τ_j (Cα j…j+3) | dimensionless, [0, 2] |
| M2 | Ordered circular D_circ | √((1/m) Σ_j d(τ_j^a, τ_j^b)²) | rad |
| B1 | Local Cα–Cα distances | RMS difference over all pairs (p, q), q − p ≥ 2, inside the window | Å |
| B2 | Window Cα RMSD | Optimal superposition (Kabsch, proper rotation) of the W Cα | Å |
| B3 | Curvature block D_κ | RMS difference of the W − 2 virtual bond angles κ | rad |
| C1 | Historical scalar-mean comparator | \|wrap(m_a − m_b)\|, m = per-frame circular mean of the window's τ_j | rad |

- C1 is a **frame-to-frame distance** derived from the scan's per-frame summary. It is **not**
  the across-frame BC scan: the BC is one number per window over an ensemble and cannot rank
  frames, so it is not used here.
- τ and κ come from the production binding `topofold.compute_invariants`; identity and reference
  check are in §10, A2.

## 5. AdK splits (per trajectory; exact, nonoverlapping)

Rule: the query is the first ⌊0.45 n⌋ frames, then a gap of ⌊0.10 n⌋ frames is excluded, and the
database is the remaining frames. Frame indices are 0-based; the cadence is 200 ps (WP-B G0).

| Trajectory | n | Query | Excluded gap | Database |
|------------|---|-------|--------------|----------|
| open_start | 2,253 | 0–1012 (1,013 frames; 0.0–202.4 ns) | 1013–1237 (225 frames; 202.6–247.4 ns) | 1238–2252 (1,015 frames; 247.6–450.4 ns) |
| closed_start | 1,678 | 0–754 (755 frames; 0.0–150.8 ns) | 755–921 (167 frames; 151.0–184.2 ns) | 922–1677 (756 frames; 184.4–335.4 ns) |

Trajectories are never pooled.

The gap reduces temporal proximity between query and database frames. It does **not** make them
independent: WP-B τ_int estimates for global observables reach about 300 frames, and the runs
drift (non-stationary). Query and database frames may therefore come from different conformational
regimes. That is a property of the test, reported as such, not corrected by burn-in removal.

## 6. Metrics

For window w, query q and database set D (|D| = n_D), each method gives a distance R(q, ·) and the
target gives T(q, ·).

**Tie-aware membership.** For distance X and rank k, let d_k be the k-th smallest value and
tol = 1e-12·(1 + d_k):
- p_X(d) = 1 if X(q, d) < d_k − tol;
- p_X(d) = (k − n_<)/n_= if |X(q, d) − d_k| ≤ tol;
- p_X(d) = 0 otherwise.

This is uniform selection from the boundary tie set, as in WP-C, and uses no labels.

| Metric | Definition |
|--------|------------|
| **Overlap@k** (k = 1, 10) | Σ_d p_R(d)·p_T(d) / k: expected shared membership under independent tie draws |
| **Retrieval regret@k** | Σ_d p_R(d)·T(q,d)/k − Σ_d p_T(d)·T(q,d)/k: mean target distance of the method's neighbours minus the target-optimal neighbours' (≥ 0; rad) |
| **Normalized regret@k** | regret / (mean_d T(q,d) − Σ_d p_T(d)·T(q,d)/k). 0 = target-optimal, 1 = random retrieval in expectation. Undefined if the denominator ≤ 1e-12, which is recorded as a degenerate query |
| **Random reference** | Expected overlap@k = k/n_D; expected normalized regret = 1 |
| **Rank correlation** | Spearman ρ between R and T over 20,000 (q, d) pairs per window, drawn uniformly from Q × D with a fixed seed (recorded). Pairs are shared by all methods |

**Aggregation, in fixed order:**
1. Mean over queries within a window.
2. Per trajectory, window size and method: median and interquartile range across windows.

Windows are also reported in terciles of target spread (median T(q, d)) to show where the target
barely varies. No window is excluded post hoc. **AdK trajectories are not pooled.** For mdCATH
(§8): per domain, temperature, window size and method, take the median across windows; then
across domains, the median with a domain-level bootstrap.

## 7. Decisions

**AdK (development tier): descriptive only.** Results are reported per trajectory and window size,
with no retain/reject decision. AdK is not held-out validation.

**mdCATH (held-out tier): PROVISIONAL rules.** Not frozen, not executed. They must be reviewed
and frozen before any mdCATH data are scored.
- Primary metric: overlap@10.
- Comparison: M1 vs the best of B1–B3 (best = highest median), and M1 vs C1.
- Analysed separately per temperature (320 K, 348 K) and per window size.
- Margin 0.05: a **pre-declared descriptive margin**, not a statistical non-inferiority test.
  Bootstrap intervals are reported but do not decide.

| Order | Decision | Rule |
|---|----------|------|
| 1 | **Insufficient evidence** | Fewer than 24 of 30 planned domains eligible and analysed; or harness checks fail; or the target is nearly constant: median T(q, d) < 0.05 rad in more than half of all windows at both temperatures. The 0.05 rad floor is a provisional, uncalibrated choice; its value must be justified before the held-out arm is frozen |
| 2 | **Reject** | M1 median overlap@10 is below the best baseline by more than 0.05 in every temperature × window-size cell |
| 3 | **Restrict** | M1 is below the best baseline by more than 0.05 in at least one cell but not all, or M1 is not above C1 in at least one cell. The role is restricted accordingly, e.g. used alongside distance features |
| 4 | **Retain for the first deliverable** | Otherwise: within 0.05 of, or above, the best baseline, and above C1, in every cell |

## 8. Held-out design (mdCATH; no download yet)

**Facts from the primary publication** (Mirarchi et al., *Sci. Data* 2024):
- per-domain HDF5 with `coords` F × N × 3 in Å over the N protein atoms (`pdbProteinAtoms`);
  `resid`, `resname`, `chain` and `element` per atom; per-residue, per-frame 8-state DSSP;
  `box` in nm;
- one frame per ns; most trajectories 500 ns (mean 464, SD 76 ns);
- CHARMM22*, TIP3P with Na⁺/Cl⁻, ACEMD, 4 fs steps with hydrogen-mass repartitioning;
  constant-volume production;
- 5 temperatures × 5 replicas, **all started from the same equilibrated configuration**;
- CC BY 4.0, total > 3 TB.

Not stated there: whether N includes hydrogens, and per-domain file sizes. Both are to be read
from file metadata.

**Sample and budget**
- 30 domains, one per CATH superfamily (30 distinct superfamilies), so each domain is its own
  cluster.
- Domain length 50–300 residues; contiguous chain segments only. Windows spanning internal breaks
  are excluded per §3.
- Budget: per-domain size × 30. The ~0.6 GB average implied by 3 TB / 5,398 domains is a crude
  estimate, to be replaced by measured sizes before authorization.

**Superfamily separation**
- Exclude every superfamily that CATH assigns to *E. coli* adenylate kinase (to be looked up and
  recorded at selection time), and any domain with > 30 % sequence identity to AdK.
- Selection: eligibility file written and hashed, then a seeded draw. No scoring before the draw
  is frozen.

**Replicas and temperatures**
- Query = replica 0, database = replica 1, at 320 K and at 348 K, analysed **separately**.
  379–450 K are excluded (the paper reports a shift toward coil at 413 and 450 K).
- **Limitation:** replicas share their starting configuration, so they are not independent
  initial conditions.

**Uncertainty:** percentile bootstrap over domains (2,000 resamples, recorded seed) for the
median difference M1 − best baseline. Reported, not used as a significance test.

## 9. Execution scope (after freeze)

**Code:** `benchmarks/audit/wpd_local_retrieval/`, harness only:
- target module (dihedrals, D_T);
- methods module, reusing WP-C `features.py` and `retrieval.py` tie rules;
- splits;
- metrics;
- runner with input-hash binding to WP-B G0 and to the eligibility artifact;
- tests: dihedral references, wrap and terminal cases, tie-aware overlap and regret on
  constructed cases, split boundaries;
- report.

No production change.

**Cost (AdK):** per window and method, ≈ 1.0e6 (open) or 5.7e5 (closed) query × database
distances. 417 windows × 6 methods + target per trajectory: minutes in vectorized NumPy.

**Stop rule:** any eligibility, hash or recovery-check failure stops the run and is reported.

## 10. AdK development-tier amendments (A1–A5, recorded before any feature or score was computed)

### A1. Degenerate dihedrals and ill-conditioned C1 means

**Degenerate dihedral.** A dihedral over points p0..p3 is degenerate in a frame if either bond
angle (p0, p1, p2) or (p1, p2, p3) has sin < 1e-3, i.e. a nearly collinear triple. This applies to:
- Cα torsions τ: the Cα virtual bond angles at the window's vertices;
- target φ/ψ: the backbone angles C–N–Cα, N–Cα–C, Cα–C–N.

A (window, frame) pair with any degenerate dihedral in its Cα torsions **or** in its target
dihedrals is excluded for **every** method and for the target alike, because the criterion is
geometric and fixed in advance. Exclusions are counted and reported as coverage. If the production
binding raises an error for a frame (collinear segments), all windows of that frame are excluded,
recorded as numerical failures, and counted.

**C1 conditioning.** Each per-frame circular mean m has a resultant length R (mean of the window's
unit torsion vectors).
- The mean is **ill-conditioned** if R < 0.1. With a declared input-uncertainty scale δ = 0.01 rad
  for xtc-quantized torsions (an order of magnitude, not a calibrated value), δ/R > 0.1 rad then;
  this is the WP-A diagnostic with δ in place of round-off.
- A query with an ill-conditioned mean has **no C1 result**: it is excluded from C1 aggregates and
  counted.
- A database frame with an ill-conditioned mean is given C1 distance +∞, i.e. ranked last and
  never treated as a valid angle; the count is reported.
- If fewer than k database frames are finite, that query's C1 result is undefined and counted.
- The distribution of R (min, 5 %, median) and the ill-conditioned fraction are reported per
  trajectory and window size.

### A2. Identical Cα inputs; separate target; binding identity

- Every Cα method receives the same Cα coordinates: the float32 values stored by MDAnalysis from
  the xtc, made whole by sequential minimum image and held in float64. The binding receives the
  same values as float32 (exactly representable, so no information differs between methods).
- The target uses backbone N, Cα, C and is computed separately with NumPy (IUPAC dihedral, WP-A
  `build_traces.dihedral`).
- **Binding identity:** the WP-A wheel `topofold-0.8.3-cp314-cp314-manylinux_2_34_x86_64.whl`,
  installed into an isolated WP-D environment. The loaded `.so` SHA-256 must equal the WP-A
  recorded value, and the production sources must equal the WP-A recorded hashes; otherwise stop.
- **Reference check:** binding τ and κ vs a NumPy recomputation from the same float64 coordinates,
  circular difference ≤ 1e-9 rad, since both compute in float64 from identical values. Applied to
  every frame. Violations are numerical failures and stop the run.
- **B2 (window Cα RMSD)** uses Horn's quaternion method (proper rotations), checked against Kabsch
  SVD on sampled pairs (≤ 1e-6 Å).

### A3. Tolerances and degenerate metric cases

- **Ties:** tol = 1e-12·(1 + d_k), as in WP-C.
- **Regret:** values in [−1e-10, 0) are rounding and are set to 0 (counted). Values < −1e-10 are
  numerical failures (counted and reported; they stop the run if they occur in the target-oracle
  computation itself).
- **Normalized regret:** undefined when (random − oracle) ≤ 1e-12, e.g. a constant target. The
  query is excluded from the normalized aggregate and counted.
- **Required tests** (SYNTHETIC fixtures): identical rankings; complete ties in a method;
  boundary ties (hand-computed); random retrieval; constant target; deliberately reversed
  retrieval. Plus dihedral, Horn-vs-Kabsch, Spearman-vs-SciPy and split-boundary checks.

### A4. Eligibility wording

See the corrected statement under the title. Numerical dihedral validity is checked by the
harness (A1), not by the eligibility artifact.

### A5. mdCATH arm

§7–§8 are provisional, and the unsupported quantization-based justification of the 0.05 rad floor
has been removed. Nothing in the mdCATH arm is frozen or executed in this task.

### Spearman specification (unchanged intent, made exact)

Per window: 20,000 (q, d) pairs drawn uniformly with replacement from the valid Q × D using
`numpy.random.default_rng(seed)`, with seed = 20261005 + 1000·W + i. Spearman ρ is the Pearson
correlation of average ranks. Pairs undefined for C1 are dropped for C1 only and counted.
