# Angular Representation Characterization Protocol

**Status:** Amended after review (amendments A–F, §9), recorded **before** any fixture was generated.
Executed as WP-A; see `benchmarks/audit/angular_characterization/REPORT.md`.
**Baseline:** `295d739a783cc1d178131c4c61293042a0c56f7d` (branch `scientific-rebuild`)
**Scope:** characterization of a fixed implementation. This is not detector validation. It sets
no accuracy targets, no monotonicity requirements and no production thresholds.

**Bounds**
- The fixture matrix in §5 is closed: adding a fixture requires protocol review.
- The task ends with the decision record in §8.
- Agreement between implementations (Rust vs NumPy) is a correctness check of the code against
  its stated formulas. It is not scientific validation of the representation.

All fixtures are SYNTHETIC descriptor-level constructions. A Cα trace built from prescribed
(κ, τ) values is a geometric test object, not a protein and not physically validated.
Units: radians in code and results; degrees only in prose.

## 1. Question

What information does the current torsion pipeline preserve or lose, and where is its numerical
behaviour sensitive? The pipeline is: per-frame circular mean over a window, then centered
deviations across frames, then Sarle's BC.

## 2. Implementation under test (as of the baseline)

`compute_bimodality_profile` (`crates/topofold-core/src/bimodality.rs`) does the following for
window `i` of size `W` (residues `i … i+W−1`):

1. `extract_curve_invariants` (`geometry.rs`) is computed per frame: tangents, then κ
   (`compute_curvatures`, length N−2), then τ (`compute_torsions` via `compute_binormals`,
   length N−3). τ_j uses atoms j…j+3.
2. Per frame:
   - `tau_mean = circular_mean(torsions[i .. i+W−3])`, i.e. W−3 torsions reduced to one angle;
   - `kappa_mean` is the arithmetic mean of `curvatures[i .. i+W−2]`;
   - θβ uses `circular_mean` only when Cβ coordinates are supplied.
3. Across frames:
   - `centered_deviations(tau_means)` (`circular.rs`) wraps each value to (−π, π] about the
     across-frame circular mean;
   - `StreamingMoments` turns the deviations into population skewness and kurtosis, plus the
     n-dependent term in the denominator.
4. `bimodality_coefficient` has fixed guards. It returns 0 if n < 4, if the variance is
   ≤ 1e-14 (an absolute threshold, in rad² for τ) or if the denominator is ≤ 1e-12. Otherwise it
   clamps the result to [0, 1].
5. `score = max(bc_tau, bc_kappa, bc_theta)`.

Two properties are visible from the code and need characterizing:

- **The scan does not compute or use resultant magnitude to assess reliability.** `circular_mean`
  forms the mean sine and cosine and returns their `atan2`. The magnitude of that vector
  (the resultant length R) is never computed in the scan, so a window whose torsions nearly
  cancel yields a direction set by round-off, used like any other. `mean_resultant_length` exists
  in `circular.rs` but is not called by the scan. `circular_mean` returns 0.0 for empty input.
- **The across-frame mean is undefined when the frame means cancel.** The centering then depends
  on round-off.

The Python binding `topofold.compute_bimodality_profile` converts float32 input to f64 before
calling the same function (`array3_to_traces`, `crates/topofold-python/src/lib.rs`). It does
not expose per-frame means.

### Production functions exercised

| Function | Path |
|----------|------|
| `circular_mean`, `mean_resultant_length`, `centered_deviations` | `topofold_core::circular` |
| `StreamingMoments`, `sarles_bimodality_coefficient`, `compute_bimodality_profile` | `topofold_core::bimodality` |
| `extract_curve_invariants` → `compute_tangents`, `compute_curvatures`, `compute_binormals`, `compute_torsions` | `topofold_core::geometry` |
| `topofold.compute_bimodality_profile`, `topofold.compute_invariants` (float32 input) | `crates/topofold-python` |

Not exercised: `detect_bistable_segments*` (thresholds and clustering), the θβ channel, RNA
functions and the index.

## 3. Execution paths

| Path | What it tests | Precision |
|------|---------------|-----------|
| **R1** Rust function level | Calls `circular_mean`, `mean_resultant_length`, `centered_deviations` and `StreamingMoments` directly on fixture angle arrays | f64 |
| **R2** Rust coordinate level | Builds Cα traces from fixture (κ, τ) and calls `compute_bimodality_profile`. It also recomputes the per-frame means by calling `extract_curve_invariants` and `circular_mean` on the same slices. The recomputed BC must equal the profile's `bc_tau` **bitwise** (justification below the table) | f64 |
| **P1** Python binding | Same traces as R2, passed as float32 through the binding | float32 coordinates |
| **REF** independent NumPy reference | Circular mean `np.angle(mean(exp(iθ)))`, explicit wrap, two-pass moments, and the BC formula transcribed from the documented definition | f64 |

REF agreeing with R1/R2 shows that the Rust code computes the stated formulas. It does not show
the formulas are appropriate. REF is not a test of the Rust implementation beyond that
equivalence.

**Why bitwise for R2.** The harness's replication calls the same production functions, on the same f64
inputs, with the same per-frame slices and the same sequential accumulation order. In
`compute_bimodality_profile` the frame invariants are computed in parallel, but each frame is
independent and the moment accumulation is sequential. The computation is therefore
deterministic, and identical operation order gives identical bits. A mismatch means the
harness is not reproducing production (wrong slice, order or function); it is not a numerical
tolerance question. The run is then invalid, not "close". Bitwise equality is required only
for this replication. Comparisons with REF use the tolerances below, because REF uses a
different summation order.

**Trace construction (harness code, not production).** Natural-extension placement with a
constant virtual bond length of 3.8 Å, prescribed κ (turning angle) and τ. Before any fixture is
used, the construction must reproduce the prescribed values through production
`extract_curve_invariants`, and through an independent NumPy dihedral, within the R2 tolerance
below. The sign convention of τ is established by this check, not assumed. All fixtures use
constant κ = 1.2 rad, so no segments are collinear.

## 4. Tolerances

| Comparison | Tolerance | Justification |
|------------|-----------|---------------|
| Constructed vs recovered κ, τ (R2, f64) | 1e-9 rad | Coordinates ≤ ~30 Å; f64 relative ε ≈ 2.2e-16 gives position error ~1e-14 Å. Angle error scales as error / (3.8 Å · sin κ) with sin 1.2 ≈ 0.93, so 1e-9 leaves margin of several orders |
| Constructed vs recovered κ, τ (P1, float32) | 1e-5 rad (flag) | Applies **only to these fixtures** (κ = 1.2 rad, coordinates ≤ ~30 Å). float32 relative ε ≈ 6e-8 gives ~2e-6 Å, i.e. ~1e-6 rad there. τ error grows as 1/sin κ near collinearity, so this is not a general bound for real structures. Values above 1e-5 are flagged, not failed |
| R1/R2 vs REF: angles (per-frame and across-frame means) | Circular distance ≤ 4·b + 1e-12 rad, where b is the conditioning bound below | See *Conditioning diagnostic*. Compared only where b ≤ 0.1 rad; otherwise the mean is classed as numerically undefined, and differences are recorded as characterization, not compared |
| R1/R2 vs REF: BC | \|ΔBC\| ≤ 10·z³·r + 1e-12, with r = e_dev / s and z = max \|standardized deviation\| | See *Conditioning diagnostic*. Compared only where the deviation set is conditioned; otherwise recorded only |
| Paired-rotation invariance (R1) | Same bound as R1 vs REF, applied to the α vs α = 0 difference | Deviations are rotation-invariant in exact arithmetic. Differences above the bound are reported as sensitivity, with the deviations nearest ±π |
| Paired-rotation invariance (P1) | Recorded only | The float32 quantization of rotated coordinates differs between rotations |

Every tolerance is a numerical-agreement criterion. None is a criterion for detector quality.

All angle comparisons use the circular distance d(a, b) = |wrap(a − b)|, with wrap to (−π, π].

**Conditioning diagnostic (harness only; fixed before results were examined).**

*Mean of n angles.* With u = 2⁻⁵³, the computed mean-resultant components each carry absolute
error ≤ (n + 1)·u: recursive summation of n terms bounded by 1, plus one rounding per sin/cos.
Perturbing every input angle by at most δ_in moves the mean vector by at most δ_in. The
direction error is therefore bounded, to first order, by

  b = ((n + 1)·u + δ_in) / R,

where R is the resultant length of that mean (b = ∞ if R = 0). δ_in = 0 for R1, which has
identical inputs to REF. For R2, δ_in is the measured maximum construction-recovery error of
the case.

- Within-frame and across-frame resultants are both recorded.
- A mean with b > 0.1 rad is classed **numerically undefined**: R is too small for its direction
  to be determined at this input precision. Implementation differences there characterize that
  undefined mean; they are not compared.

*BC.* Deviation errors are bounded by e_dev = max over frames of b_frame + b_across.
The deviation set is **conditioned** when all three hold:
1. every b ≤ 0.1 rad;
2. no deviation lies within e_dev of ±π, where wrapping is ambiguous;
3. r = e_dev / s ≤ 1e-6, with s the deviation standard deviation.

For conditioned sets, perturbing standardized values by r changes skewness and kurtosis by at
most ~4·z³·r (z = max |standardized deviation|), hence the BC tolerance 10·z³·r + 1e-12. When
the REF and Rust variances fall on opposite sides of the production guard (1e-14), the case is
recorded as *guard-ambiguous* and not compared.

None of these numbers is a production threshold. They decide only which harness comparisons are
meaningful.

## 5. Fixture matrix

Seeds: NumPy `default_rng(20261004 + k)` for fixture k. Wrapped normal noise is
`wrap(θ + N(0, σ²))`. F is the number of frames. Each fixture is a single window spanning the
whole trace (N = W), except where stated.

### F1. Coordinate-wise state differences hidden by averaging

**Construction**
- F1a: W = 5, so the window has 2 torsions. Frames alternate between state A, τ = (30°, 90°), and
  state B, τ = (90°, 30°). F = 1,000 (500 per state, interleaved).
- Noise σ ∈ {0, 2°}, iid per torsion.
- F1b: as F1a with W = 8 (5 torsions): A = (30°, 90°, 60°, 60°, 60°) and
  B = (90°, 30°, 60°, 60°, 60°).
- F1c (positive control): W = 5, A = (30°, 30°), B = (90°, 90°), σ ∈ {0, 2°}.

**Analytical expectations (σ = 0)**
- F1a: both states have per-frame circular mean exactly 60° and R = cos 30° ≈ 0.866.
- F1a: the across-frame deviations are zero up to round-off, so the variance guard may trigger.
  Whether it triggers is recorded, not predicted.
- F1b: per-frame means are equal for both states (60°).
- In the coordinate-wise embedding (cos τ₁, sin τ₁, cos τ₂, sin τ₂), the state means are separated
  by chord 2 sin 30° = 1 in each of coordinates 1 and 2, i.e. Euclidean distance √2.
- F1c: per-frame means are 30° and 90°, so the scalar summary separates the states. This shows
  the harness can register a difference.

**Record**
- Per-frame τ, `tau_mean`, R and `kappa_mean`.
- Across-frame mean, deviations, variance, skewness, excess kurtosis, BC and whether a guard
  triggered.
- Per-state means of the scalar summary and of the embedding.

**Decision informed.** Can the scalar summary distinguish coordinate-wise states that the
embedding separates? If it cannot, the per-frame circular mean is a lossy representation for
this class of change. The choice is then between keeping it as a documented narrow heuristic and
moving to a per-coordinate (sin, cos) representation. This fixture alone does not decide which.

### F2. Paired angular rotation of identical wrapped observations

**Construction.** Base samples, with rotation applied identically to every torsion of the
**same** draws: τ′ = wrap(τ + α), α ∈ {0°, 45°, 90°, 135°, 170°, 179°, 180°, −179°}.
- F2a: the F1a σ = 2° draws.
- F2b: W = 5, each frame has both torsions equal to one draw from a two-state mixture with
  centers ±40° about 0°, σ = 10°, equal weights, F = 1,000.
- F2c (wrap stress): as F2b with centers ±100° and σ = 30°, so that some across-frame deviations
  approach ±π.

**Analytical expectations.** The circular mean is rotation-equivariant and the centered
deviations are rotation-invariant, except where a deviation lies within round-off of ±π and can
wrap to the opposite sign. BC is therefore invariant up to round-off, except through such wraps.
(F2c also changes the geometry of the constructed traces under rotation; κ is unchanged.)

**Record**
- For each α: BC, max |Δ deviation| vs α = 0, and the count and values of deviations with
  |d| > π − 1e-6.
- Run on R1, R2 and P1.

**Decision informed.**
- Invariance violated beyond tolerance away from ±π → defect in centering or wrapping; fix
  before any further use.
- Violations confined to ±π wraps → a characterized sensitivity, to be documented.

### F3. Within-frame near-zero resultant

**Construction.** W = 5, so each frame has two torsions. Frame f has
τ_f = (θ₀ + η_f, θ₀ + η_f + π − ε), both wrapped to (−π, π], with:
- θ₀ = 20°;
- ε ∈ {1e-1, 1e-3, 1e-6, 1e-9, 0, −1e-9, −1e-6, −1e-3, −1e-1} rad;
- F = 200;
- η_f ~ N(0, σ_n²), σ_n ∈ {0, 1e-6 rad, 1°}. One draw per frame, added to **both** torsions, so
  the within-frame separation stays π − ε up to the rounding of the wrap.

**Analytical expectations (exact arithmetic)**
- The per-frame resultant vector is sin(ε/2)·exp(i[θ₀ + η_f + (π − ε)/2]), so R = |sin(ε/2)|.
- For ε > 0 the per-frame mean is θ₀ + η_f + (π − ε)/2.
- For ε < 0, sin(ε/2) < 0 and the mean is θ₀ + η_f + (π − ε)/2 + π (an extra π branch).
- At ε = 0 the mean is undefined.
- For |ε| = 1e-9, R ≈ 5e-10, so b ≫ 0.1 rad: the computed mean is expected to be numerically
  undefined even though it is mathematically defined.

**Record**
- Per-frame mean, R, the across-frame quantities as in F1, and BC.
- Sensitivity: change in per-frame mean when one torsion is perturbed by ±1e-9 and ±1e-6 rad.
- Run on R1, R2 and P1 (on P1 the float32 quantization acts as an additional perturbation).

**Decision informed.** If outputs are dominated by round-off or quantization at low R, the
pipeline needs to carry R (masking, weighting or reporting) or change representation. No R
threshold is chosen here.

### F4. Across-frame cancellation (kept separate from F3)

**Construction.** W = 5, both torsions equal within a frame, so per-frame R = 1. Frame values
are either φ₁ = 0° or φ₂ = π − ε, plus a per-frame noise draw ν_f ~ N(0, σ²) added to both
torsions (σ ∈ {0, 1°}). Parameters:
- fraction p ∈ {0.5, 0.49, 0.45} of frames at φ₂, assigned as the first round(p·F) frames of a
  seeded permutation;
- ε ∈ {0, ±1e-6, ±1e-3, ±1e-1} rad;
- F = 1,000.

Non-antipodal control: φ₂ = 120°, p = 0.5.

Each F4 case is also evaluated on R1 and REF at every F2 rotation α, applied to the same draws.

*Matrix amendment C, review-authorized:* the negative-ε cases are added. R2/P1 run at α = 0 only.

**Analytical expectations (noiseless constructions only)**
- Across-frame resultant: |(1 − p) + p·e^{i(π−ε)}|. At p = 0.5 this is |sin(ε/2)|, which is 0 for
  ε = 0. There the across-frame mean, and every deviation, is determined by round-off; BC is not
  predicted.
- For noisy cases no exact resultant is claimed; the finite-sample resultant is recorded.

**Record**
- Across-frame R and mean.
- The deviation set: for σ = 0, the distinct values and their counts.
- BC, and its change under the F2 rotations.

**Decision informed.** Does the centering step make BC depend on an arbitrary reference direction
when frame means cancel? If yes, either the centering rule must change or this configuration must
be detected and reported.

### F5. Small perturbations around antipodal configurations

F3 (within frame) and F4 (across frames) already include ε grids on both sides of exact
antipodality. F5 is the joint analysis of those runs; it adds no new data.

**Record:** outputs as a function of the sampled ε values, including where they change abruptly.

**Decision informed.** Only sampled sensitivity bounds are reported, e.g. "at |ε| ≤ 1e-6 the
mean was undefined; at |ε| ≥ 1e-3 it agreed with the analytical value". The width of the
unreliable region is not established between sampled values.

### F6. Existing controls from `benchmarks/audit/bc_check.py`

**Construction**
- F6a (scalar BC, R1 vs REF): the inputs of `bc_check.py` sections 2–3, reproduced without
  executing that script by replaying its exact RNG draw order on `numpy.random.default_rng(1)`:
  1. Section 1 is consumed and discarded: five calls `normal(mu, sd, size=(2500, 5))` with
     (mu, sd) = (−170, 15), (−175, 15), (180, 10), (−120, 15), (50, 10).
  2. Section 2, for p in (0.5, 0.1, 0.02, 0.005): `x = normal(0, 0.05, 2500)`, then
     `random(2500) < p` selects the excursions, which get +1.0.
  3. Section 3, in list order: `exponential(1, 2500)`, `lognormal(0, 0.8, 2500)`,
     `gamma(1, 1, 2500)`, `chisquare(1, 2500)`, `normal(size=2500)`.

  The claim of identical historical samples rests on this replay and on NumPy's stable
  `PCG64`/`Generator` streams. As a consistency check, not proof, the resulting BC values must
  round to the values documented in KNOWN_ISSUES M3 (0.744, 0.839, 0.577, 0.590). If they do not,
  the controls are relabelled as newly generated and no historical identity is claimed.
  - Gaussian(0, 0.05) + p excursions of +1.0, p ∈ {0.5, 0.1, 0.02, 0.005}.
  - Exponential(1), lognormal(s = 0.8), gamma(k = 1), χ²₁, Gaussian.
- F6b (angular path, R2/P1): rare-excursion inputs only, used as per-frame τ values. Both torsions
  equal x; W = 5; F = 2,500.
  - The values stay inside (−π, π]; this is checked and recorded.
  - Skewed controls are excluded from the angular path because exponential and χ² values can
    exceed π and would wrap, which changes the construction.

**Analytical expectations**
- R1 and REF agree within the §4 tolerance on identical inputs.
- On F6b, BC of the coordinate path equals BC of the scalar path applied to the same values,
  within tolerance. Two facts give this: the per-frame circular mean of two equal torsions is the
  value itself, and moments are shift-invariant. It holds only if no centered deviation wraps
  past ±π; this is checked and recorded.
- No BC values are predicted.

**Record:** BC, skewness, excess kurtosis and n per input on every path.

**Decision informed.** Confirms whether the documented non-specificity of BC (KNOWN_ISSUES M3)
carries through the angular path unchanged. This bears on retain vs replace; it does not decide
it.

### Cross-cutting record: the κ channel

All fixtures use constant κ. Record `bc_kappa` and the κ variance on R2 and P1:
- on R2 the variance is expected ≤ 1e-14 (guard → 0);
- on P1, float32 quantization may leave a variance above the guard. The 1e-14 guard is
  absolute.

Any non-zero `bc_kappa` there enters `score = max(…)`. This is recorded as a finding about the
guard and the max composition, not as a property of τ.

## 6. Interpretation rules

- Report three things separately:
  1. **harness correctness:** recovery, replication, reference agreement where conditioned,
     analytical controls;
  2. **observed numerical sensitivity:** behaviour where means are ill-conditioned or undefined;
  3. **scientific limitations of the representation:** information loss.
- Observed information loss is demonstrated **for the tested constructions only**. Absence of
  loss elsewhere is not shown and must not be inferred.
- Sensitivity is reported as sampled bounds, never as an established region width.
- Each outcome describes this implementation on these fixtures only. Do not generalize from an
  isolated fixture.
- Do not infer the number of modes from BC, from resultant length, or from mixture component
  counts.
- Do not convert any recorded value into a threshold.
- Possible overall conclusions remain open:
  - keep the scalar pipeline as a narrow, documented heuristic;
  - change the representation (e.g. per-coordinate embedding, carrying R);
  - replace the detector.

## 7. Execution and artifact scope for the subsequent implementation task (WP-A)

| Artifact | Purpose |
|----------|---------|
| `benchmarks/audit/angular_characterization/generate_fixtures.py` | Writes `fixtures/SYNTHETIC_*.f64` (raw little-endian f64) and `fixtures/SYNTHETIC_fixtures.json` (shapes, seeds, parameters, SHA-256) |
| `benchmarks/audit/angular_characterization/build_traces.py` | Harness trace construction plus the recovery check (§3) |
| `crates/topofold-core/examples/angular_characterization.rs` | R1/R2 harness. Reads the fixtures and calls only public production functions. Writes `results/SYNTHETIC_angular_characterization_rust.json`. A Cargo example: no change to `src/` |
| `benchmarks/audit/angular_characterization/run_python_path.py` | P1 via the installed bindings, and REF. Writes `results/SYNTHETIC_angular_characterization_python.json` |
| `benchmarks/audit/angular_characterization/compare.py` | Applies the §4 tolerances and writes `results/SYNTHETIC_angular_characterization_summary.json` |
| `benchmarks/audit/angular_characterization/test_harness.py` | Tests of the harness itself: construction recovery, fixture hashes, bitwise R2 replication |
| `benchmarks/audit/angular_characterization/REPORT.md` | Reads every number from the results JSON |

**Tests in `test_harness.py`**
- The trace construction recovers prescribed κ, τ within §4 tolerances, through production and
  through the independent NumPy dihedral.
- Fixture files match the hashes recorded in `SYNTHETIC_fixtures.json`.
- The R2 replication equals `compute_bimodality_profile` bitwise for one fixture per family.
- The REF implementation reproduces analytical values: the F1 per-frame mean of 60°; R = cos 30°;
  R = |sin(ε/2)| in F3; the F4 across-frame resultant.
- The comparison script fails (nonzero exit) on a deliberately altered copy of a results file in
  a temporary directory.

**Completion criteria**
- Every fixture in §5 has run on every applicable path, with outputs in the three results JSON
  files.
- All harness tests pass. `cargo test --workspace` and `pytest tests/` still pass.
- `REPORT.md` quotes only numbers present in the results JSON.
- The §8 decision record is completed for review.

Constraints for that task:
- no change to `crates/*/src`;
- no thresholds;
- no real data;
- `cargo test --workspace` and `pytest tests/` must still pass;
- all result files carry SYNTHETIC in their names.

## 8. Decision record (completed at the end of WP-A, then reviewed)

The task ends with one recorded decision for the scalar circular-mean/BC component: retain,
restrict use, change, or **insufficient evidence**. The decision applies to this component
only; geometry and retrieval capabilities are unaffected.

| Decision | Evidence that supports it |
|----------|---------------------------|
| **Retain** as a documented screening heuristic | F1 loss is confined to the documented class; F2 invariance holds up to round-off except ±π wraps; F3–F5 unreliable regions are narrow and detectable; F6 paths agree |
| **Restrict use** | Losses or sensitivities are real but bounded and detectable. Specify the restriction, e.g. report R or flag antipodal windows, exclude the τ channel for affected window sizes, or report `bc_kappa` separately from `max(…)` |
| **Change** | The losses or instabilities affect configurations that cannot be detected from the outputs, or invariance fails away from ±π. Name the candidate change (e.g. per-coordinate (sin, cos) representation, carrying R, a different centering rule). The candidate needs its own protocol and validation; it is not adopted by this task |
| **Insufficient evidence** | The fixtures do not discriminate between the options above, or harness correctness could not be established |

The record must cite the results-file entries behind the decision, list the fixtures that did
not inform it, and state that the outcome is a property of this implementation on these
fixtures. It is not evidence about real proteins.

## 9. Amendment record (review round 2, before fixture generation)

| ID | Amendment | Reason |
|----|-----------|--------|
| A | F3 analytical mean corrected: resultant sin(ε/2)·exp(i[θ₀+(π−ε)/2]); extra π branch for ε < 0; undefined at ε = 0. Per-frame noise is a common shift applied to both torsions | Previous text omitted the branch and left the noise construction ambiguous |
| B | Angle comparisons use circular distance. Both resultants are recorded. A conditioning diagnostic, a pre-declared error bound, replaces the variance-only criterion | Variance alone does not determine whether a mean is defined |
| C | F4 includes negative ε; exact resultants claimed for noiseless cases only; F4 rotations specified (R1/REF only) | Symmetry of the antipodal neighbourhood; noisy resultants are finite-sample quantities |
| D | Production description: the scan does not compute or use resultant magnitude | "Computed but discarded" was inaccurate |
| E | F6 RNG draw order specified, with replay of section 1 draws and a consistency check | Identity with historical samples must be shown, not assumed |
| F | Interpretation separates correctness, sensitivity and limitations; sampled bounds only; "insufficient evidence" allowed | Prevents over-generalization |

The R2 bitwise requirement is unchanged and applies only to exact replication of the same
production operations on identical f64 inputs. It is not applied to REF or to float32 inputs.

### Post-execution notes (no check was changed)

- **F3 expectation corrected.** The statement that |ε| = 1e-9 is "expected to be numerically
  undefined" was wrong under the pre-declared diagnostic: with exact f64 inputs, b is far below
  0.1 rad. The observed sensitivity comes from input perturbations comparable to R. See
  `benchmarks/audit/angular_characterization/REPORT.md` §5.
- **F2c wrap stress not reached.** F2c produced no deviation within 1e-6 of ±π. Wrap behaviour
  was exercised by F4 instead.
- **Artifact layout.** Regenerable binaries are under
  `benchmarks/audit/angular_characterization/build/`. Fixture hashes are in
  `fixtures/SYNTHETIC_fixtures.json`; array hashes are in the results files.
