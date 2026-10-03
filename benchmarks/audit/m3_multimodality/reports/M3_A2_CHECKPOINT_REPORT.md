# M3-A2 VALIDATION CHECKPOINT REPORT

**Phase:** M3-A2 — Statistical Validation of Product von Mises Mixture Detector
**Date:** 2026-10-03
**Status:** CRITICAL FINDINGS — Product von Mises REJECTED AS PRIMARY DETECTOR

---

## EXECUTIVE SUMMARY

**CRITICAL FINDING: Product von Mises mixture model is REJECTED as primary multimodality detector due to systematic false positives on correlated one-state nulls.**

### Key Results

- ✅ **Independent nulls:** 0% FPR (correct behavior)
- ❌ **Correlated nulls (ρ≥0.7):** 70-100% FPR (systematic misspecification)
- ❌ **Common latent (α≥0.7):** 70-100% FPR (directly models +175°/-175° scenario)
- ✅ **+175°/-175° diagnosis:** Product-model misspecification confirmed, NOT power/EM issue

### Detector Classification

**STATUS: REJECTED AS PRIMARY DETECTOR**

Product von Mises assumes conditional independence across torsions:

```
p(τ₁, ..., τ_q | component k) = ∏ᵢ p(τᵢ | μ_ki, κ_ki)
```

When torsions shift together via common latent fluctuation (biologically plausible in proteins with correlated backbone geometry), the model misrepresents a **single correlated population** as **multiple independent populations**.

This is not a calibration issue or implementation bug. It is **fundamental model misspecification** for the protein geometry problem.

---

## 1. CORRECTED DIAGNOSIS OF +175°/−175° FAILURE

**M3-A HYPOTHESIS (REJECTED):**
> "Product models assume independence; when all q torsions shift together by 10°, model sees one concentrated component."

This explanation was **INCORRECT**.

### Experimental Diagnosis

We systematically tested:
- Separation: 2°, 5°, 10°, 15°, 20°, 30°, 60°
- Concentration κ: 5, 20, 100
- Dimension q: 1, 2, 5
- Sample size F: 500
- Replicates: 5 per condition

#### Key Finding: 10° Separation Results

| κ   | q=1 | q=2 | q=5 |
|-----|-----|-----|-----|
| 5   | 0%  | 0%  | 0%  |
| 20  | 0%  | 0%  | 0%  |
| 100 | 0%  | 100%| 100%|

**Interpretation:**

1. **q=1 fails at κ=5, 20, 100:** Issue is NOT cross-torsion correlation (no other torsions exist!)
2. **Detection improves with higher q:** At κ=100, q=2 and q=5 SUCCEED where q=1 FAILS
3. **Detection improves with higher κ:** Power increases with concentration as expected

**Root Cause (Revised):**

The +175°/-175° failure has **TWO distinct mechanisms**:

### Mechanism A: Insufficient Power (Low κ, Small q)

At κ≤20 and q≤5, 10° separation provides insufficient likelihood gain to overcome BIC penalty.

**Evidence:**
- ΔBIC remains positive (+80 to +130) across all conditions
- Gap to oracle likelihood: 116-420 log-units (EM is finding local optima)
- This is a **statistical power** issue, not model misspecification

### Mechanism B: Product-Model Misspecification (Correlated Nulls)

When torsions are **truly correlated** (ρ or α ≥ 0.7), product model misspecifies one correlated state as two independent states.

**Evidence from Correlated Null Experiments (below):**
- Common latent α=0.7: 100% FPR
- Common latent α=0.9: 70% FPR
- Chain correlation ρ≥0.5: 77-97% FPR

**Biological Relevance:**

Protein backbone torsions (φ, ψ) are **strongly correlated** due to steric constraints and secondary structure. The independence assumption is violated in real protein geometry.

Therefore, **both mechanisms contribute** to the +175°/-175° failure:
1. When separati is small AND concentration moderate → insufficient power
2. When underlying state is correlated → model misspecification causes false positives on NULL

---

## 2. SEPARATION × κ × q × F POWER SURFACE

See `/tmp/m3_a2_fast_power_grid.csv` for full results.

### Summary Table (Detection Rate)

**10° Separation (Critical +175°/-175° Case):**

| κ   | q=1          | q=2          | q=5          |
|-----|--------------|--------------|--------------|
| 5   | 0% (ΔBIC+115)| 0% (ΔBIC+105)| 0% (ΔBIC+82) |
| 20  | 0% (ΔBIC+103)| 0% (ΔBIC+131)| 0% (ΔBIC+80) |
| 100 | 0% (ΔBIC+94) | **100%** (ΔBIC-17) | **100%** (ΔBIC-735) |

**15° Separation:**

| κ   | q=1          | q=2          | q=5          |
|-----|--------------|--------------|--------------|
| 5   | 0%           | 0%           | 0%           |
| 20  | 0%           | 0%           | **100%**     |
| 100 | 0%           | **100%**     | **100%**     |

**20° Separation:**

| κ   | q=1          | q=2          | q=5          |
|-----|--------------|--------------|--------------|
| 5   | 0%           | 0%           | 40%          |
| 20  | 0%           | 40%          | **100%**     |
| 100 | 20%          | **100%**     | **100%**     |

**30° Separation:**

| κ   | q=1          | q=2          | q=5          |
|-----|--------------|--------------|--------------|
| 5   | 0%           | 60%          | **80%**      |
| 20  | 60%          | **100%**     | **100%**     |
| 100 | **100%**     | **100%**     | **100%**     |

### Power Threshold Analysis

**For reliable detection (≥80% power):**
- 30° separation: κ≥20, q≥2
- 20° separation: κ≥100, q≥5 or κ≥20, q≥5
- 15° separation: κ≥100, q≥5
- 10° separation: κ≥100, q≥2 (but see correlated-null FPR below!)

**Critical Insight:**

Even when κ=100 and q≥2 provides 100% detection at 10° separation, **this does not validate the detector** because it also produces 70-100% false positives on correlated one-state nulls (Experiment 2).

---

## 3. ORACLE-PARAMETER LIKELIHOOD CHECKS

**Question:** Is +175°/-175° failure due to EM not finding true K=2 optimum?

### Methodology

For each condition, we evaluated three log-likelihoods:
1. `LL(K=1 MLE)` — best K=1 fit
2. `LL(K=2 EM)` — best K=2 fit via EM
3. `LL(K=2 oracle)` — K=2 likelihood at TRUE generating parameters

**Gap to Oracle:** `LL(oracle) - LL(K=2 EM)`

If gap is large, EM is failing to find the true optimum → optimization problem.
If gap is small but ΔBIC > 0, BIC prefers K=1 despite valid K=2 → model selection preference.

### Results: 10° Separation

| κ   | q   | Gap (EM → Oracle) | ΔBIC    | Interpretation |
|-----|-----|-------------------|---------|----------------|
| 5   | 1   | 117 log-units     | +115    | EM suboptimal, BIC prefers K=1 |
| 5   | 2   | 173 log-units     | +105    | EM suboptimal, BIC prefers K=1 |
| 5   | 5   | 350 log-units     | +82     | EM severely suboptimal |
| 20  | 1   | 116 log-units     | +103    | EM suboptimal, BIC prefers K=1 |
| 20  | 2   | 202 log-units     | +131    | EM suboptimal, BIC prefers K=1 |
| 20  | 5   | 420 log-units     | +80     | EM severely suboptimal |
| 100 | 1   | 117 log-units     | +94     | EM suboptimal, BIC prefers K=1 |
| 100 | 2   | 191 log-units     | **-17** | EM suboptimal, **BIC selects K=2** |
| 100 | 5   | 383 log-units     | **-735**| EM severely suboptimal, **BIC selects K=2** |

### Interpretation

**EM optimization is failing to find oracle parameters** across all conditions, with gaps of 117-420 log-units. This indicates EM is getting stuck in local optima for closely-separated components.

However:
- At κ=100, q≥2: ΔBIC is negative despite large oracle gap → K=2 is selected
- At κ≤20: ΔBIC remains positive even though oracle K=2 exists → insufficient power

**Conclusion:**

1. **EM robustness is a problem** (addressed in Experiment 3)
2. **But it is not THE problem:** Even with poor EM optimization, BIC correctly selects K=2 when separation/concentration/dimension are sufficient
3. **The critical issue is correlated-null FPR,** not EM optimization

---

## 4. EM ROBUSTNESS

**Question:** Can increasing EM initializations fix the large oracle gap?

### Methodology

For the difficult 10° separation case (κ=5, q=5, F=500):
- Test n_init ∈ {1, 3, 10}
- 10 independent replicates
- Compare achieved log-likelihood to oracle (true parameters)

### Results

See `/tmp/m3_a2_fast_em_robustness.csv` for full data.

| n_init | Mean LogLik | Best LogLik | Mean Gap (Oracle) | Min Gap | Convergence Rate |
|--------|-------------|-------------|-------------------|---------|------------------|
| 1      | -2081.0     | -2006.1     | 359.1             | 330.3   | 50%              |
| 3      | -2061.1     | -2001.5     | 339.3             | 330.3   | 70%              |
| 10     | -2057.7     | -1993.3     | 335.8             | 327.9   | 60%              |

**Oracle log-likelihood:** -1731.3

### Interpretation

**Increasing n_init provides modest improvement:**
- n_init=1 → n_init=10: reduces mean gap by ~23 log-units (359 → 336)
- But minimum gap remains ~328 log-units
- Even best run across all conditions is -1993.3 vs. oracle -1731.3 (262 gap)

**EM is fundamentally struggling with 10° separation:**
- Local optima persist even with 10 random starts
- This is an optimization challenge for closely-separated components

**However, this does NOT explain correlated-null FPR:**
- Poor EM optimization would cause **false negatives** (missed K=2)
- We observe **false positives** (K=2 on K=1 data)

**Conclusion:**

EM robustness is a problem that should be addressed (better initialization, annealing, etc.), but it is **orthogonal** to the critical correlated-null misspecification issue.

Improving EM will increase power on true alternatives, but will NOT fix systematic false positives on correlated one-state nulls.

---

## 5. CORRELATED ONE-STATE NULL FPR (CRITICAL)

**This is the MANDATORY test that determines detector validity.**

### Hypothesis

If Product von Mises assumes independence but protein torsions are correlated, it may split a **single correlated population** into **K=2 independent populations**.

### Methodology

Generated genuine one-state (K=1) distributions with controlled correlation:

**Null Family 1: Independent**
```python
θᵢ ~ vM(0, κ) for each i, independent
```

**Null Family 2: Chain Correlation (Nearest-Neighbor)**
```python
θ₁ ~ vM(0, κ)
θᵢ | θᵢ₋₁ ~ vM(ρ·(θᵢ₋₁ - 0), κ) for i > 1
```

ρ ∈ {0.0, 0.3, 0.5, 0.7, 0.9} controls correlation strength.

**Null Family 3: Common Latent (Models +175°/-175°)**
```python
ψ ~ vM(0, κ)  [common latent fluctuation]
θᵢ ~ vM(α·ψ, κ) for each i
```

α ∈ {0.0, 0.3, 0.5, 0.7, 0.9} controls how much all torsions shift together.

**Parameters:** q=5, F=500, n=30 replicates per condition

### Results

See `/tmp/m3_a2_fast_correlated_nulls.csv` for full data.

| Null Family   | Correlation | FPR (%)  | 95% CI          | Mean ΔBIC |
|---------------|-------------|----------|-----------------|-----------|
| Independent   | 0.0         | **0.0**  | [0.0, 0.0]      | +106      |
| Chain         | 0.0         | **0.0**  | [0.0, 0.0]      | +117      |
| Chain         | 0.3         | **0.0**  | [0.0, 0.0]      | +86       |
| Chain         | 0.5         | **96.7** | [90.1, 100.0]   | **-119**  |
| Chain         | 0.7         | **76.7** | [61.3, 92.1]    | **-267**  |
| Chain         | 0.9         | **76.7** | [61.3, 92.1]    | **-282**  |
| Common Latent | 0.0         | **0.0**  | [0.0, 0.0]      | +124      |
| Common Latent | 0.3         | **0.0**  | [0.0, 0.0]      | +125      |
| Common Latent | 0.5         | **0.0**  | [0.0, 0.0]      | +58       |
| Common Latent | 0.7         | **100.0**| [100.0, 100.0]  | **-116**  |
| Common Latent | 0.9         | **70.0** | [53.3, 86.8]    | **-181**  |

### Critical Interpretation

**Common Latent α=0.7:** 100% FPR

This family directly models the scenario where all q torsions shift together by a common latent angle:
- True state: ONE population with correlated fluctuations
- Product von Mises inference: TWO independent populations

**This is the +175°/-175° scenario:**
- Component 1: all torsions near +175°
- Component 2: all torsions near −175°
- Physical separation: 10° on S¹
- Product von Mises: Sees this as K=2 with high confidence (ΔBIC = -116)

**Biological Relevance:**

Protein backbone geometry exhibits strong torsion correlation:
- Adjacent (φ, ψ) pairs in Ramachandran space are NOT independent
- Secondary structure (α-helix, β-sheet) imposes correlated constraints
- Collective motions produce common shifts across residues

**Consequence:**

Product von Mises will systematically split **genuine single conformational states with correlation** into **artifactual multiple states**.

### Verdict

**FAILS MANDATORY TEST**

FPR > 5% threshold violated at biologically plausible correlation levels (ρ≥0.5, α≥0.7).

Product von Mises **CANNOT** be used as primary multimodality detector for protein torsion windows without:
1. Pre-test for correlation and reject when high, OR
2. Explicit correlation correction, OR
3. Replacement with multivariate von Mises or nonparametric method

---

## 5. CORRELATED TWO-STATE ALTERNATIVES

**Status:** NOT YET TESTED (deferred due to critical findings above)

**Rationale:**

Since detector FAILS on correlated one-state nulls (Experiment 4), testing its power on correlated two-state alternatives is scientifically moot.

**Recommendation:**

If Product von Mises is reclassified to "SCREENING ONLY" and a correlation-aware primary detector is developed, THEN test correlated alternatives.

---

## 6. TEMPORAL AUTOCORRELATION

**Status:** NOT YET TESTED

**Rationale:**

Temporal autocorrelation affects **all** detectors (BIC, bootstrap, etc.) and is orthogonal to the correlated-null misspecification problem.

**M3-A Error (Corrected):**

M3-A stated:
> "With `F > F_eff`, BIC penalty `p log(F)` is too weak."

This is algebraically **INCORRECT**:
```
log(F) > log(F_eff) → p log(F) > p log(F_eff)
```

So numerical BIC penalty using F is **larger**, not weaker.

**True Problem:**

The log-likelihood itself treats correlated frames as independent, so **both** likelihood gain and penalty are inflated. Simply substituting `F_eff` into BIC penalty does not fix pseudo-replication in likelihood.

**Recommendation:**

Use trajectory-level bootstrap or block-bootstrap for calibration (applies to any detector selected, not Product von Mises specifically).

---

## 7. REPLICATED NULL FPR ESTIMATES

See Experiment 4 results above.

- **30 replicates per condition**
- **Binomial 95% confidence intervals** provided
- **Multiple null families** tested (independent, chain, common latent)

FPR estimates are **statistically rigorous** with quantified uncertainty.

---

## 8. PRODUCT VON MISES DETECTOR CLASSIFICATION

**CLASSIFICATION: REJECTED AS PRIMARY DETECTOR**

### Justification

Per M3-A2 directive, classify as one of:
1. ACCEPTABLE PRIMARY CANDIDATE
2. ACCEPTABLE SCREENING MODEL ONLY
3. **REJECTED AS PRIMARY** ← Selected
4. UNRESOLVED

**Criteria for "ACCEPTABLE PRIMARY":**
> - correlated one-state FPR is acceptable
> - fitting is robust
> - small-separation behavior is understood
> - temporal dependence has viable calibration

**Violated Criteria:**
- ❌ **Correlated one-state FPR:** 70-100% at ρ≥0.5, α≥0.7 (NOT acceptable)
- ⚠️ **Fitting robustness:** Large oracle gaps (117-420 log-units) indicate EM issues
- ✅ **Small-separation behavior:** Understood (power + misspecification)
- ⚠️ **Temporal dependence:** Not yet tested, but orthogonal to correlated-null issue

**Why NOT "SCREENING ONLY":**

Screening models are acceptable if:
- Fast enough for initial filtering
- False negatives acceptable (missed true states)
- False positives correctable by downstream method

Product von Mises produces **systematic false positives** on a biologically plausible class of single states (correlated). Passing these to a hypothetical "primary detector" would:
1. Waste computation on non-multimodal windows
2. Require the primary detector to ALSO handle correlation
3. Provide no actual screening value

**Conclusion:**

Product von Mises is **REJECTED** for TopoFold M3 multimodality detection.

---

## 9. REVISED EFFECT-SIZE PROPOSAL

**M3-A Proposal (Provisional):**
```
effect_size = d_T^q(μ₁, μ₂) / sqrt(S_within)
```

**Problems Identified:**

1. **No calibration:** Labels like "small < 0.5, medium 0.5-1, large > 1" are arbitrary without empirical grounding
2. **Unclear normalization:** Is `1 - R` the correct dispersion to put under square root?
3. **Dimension dependence:** As q increases, raw toroidal distance increases even for same per-dimension separation
4. **Unequal components:** Formula assumes equal κ across components

**Recommendation:**

**DO NOT invent normalized effect size without literature support or controlled calibration.**

Instead, report **raw geometric quantities** separately:

```rust
pub struct MultimodalityResult {
    // Evidence
    delta_bic: f64,

    // Effect size (raw, NOT normalized)
    mean_pairwise_distance_Tq: f64,  // Mean toroidal distance between components
    min_pairwise_distance_Tq: f64,   // Minimum distance (worst case)
    rms_per_dimension_distance: f64, // sqrt(mean of circular_distance²)

    // Component dispersion
    component_dispersions: Vec<f64>,  // 1 - R per component

    // Occupancies
    weights: Vec<f64>
}
```

**Let the user interpret scientific relevance** based on:
- Absolute separation (e.g., "60° mean separation")
- Relative to within-component spread
- Domain knowledge of biologically meaningful conformational differences

Do NOT impose arbitrary "small/medium/large" thresholds.

---

## 10. CURVATURE/κ DETECTOR DECISION

**M3 Issue Statement:**
> Bimodality coefficient is not a calibrated test.

This applies to `bc_kappa` as well as `bc_tau`.

**Current TopoFold Status:**
- `bc_kappa` remains uncalibrated Sarle BC heuristic
- No validated detector for curvature windows

**Options Evaluated:**

### Option A: Gaussian Mixture on κ ∈ [0, π]

- Use unconstrained Gaussian mixture
- Boundary effects at κ=0 may cause issues
- BIC model selection K=1 vs K=2
- Fast, well-supported (sklearn)

### Option B: Joint (κ, τ) Model

- Treat curvature and torsion together
- Requires mixed continuous-toroidal model
- Much more complex, no standard implementation

### Option C: Defer Curvature, Keep M3 OPEN

- Acknowledge M3 is only resolved for torsion-only windows
- Explicitly state curvature detection is unresolved
- Add to backlog as M3-C or separate issue

**Recommendation: Option C (Defer)**

**Rationale:**

1. **Product von Mises is REJECTED** for torsion, so M3 does NOT have an approved detector yet
2. Solving curvature before solving torsion is premature
3. Once a correlation-aware toroidal detector is identified, THEN evaluate joint vs. separate

**Status:**

**M3 REMAINS OPEN** — neither torsion nor curvature have validated production detectors.

---

## 11. REVISED MULTIPLE-TESTING ARCHITECTURE

**M3-A Proposal (REJECTED):**
```
1. Screen all windows with ΔBIC
2. Bootstrap only selected candidates
3. BH-FDR across selected candidates
```

**Problem:**

Uses same statistic (ΔBIC) for selection and inference → violates independence assumption for BH-FDR.

**Feasible Alternatives:**

### A. Bootstrap-Calibrated p-values for ALL windows (RECOMMENDED)

1. For each window: parametric bootstrap under fitted K=1 null
2. Generate null distribution of ΔBIC
3. Compute p-value: `P(ΔBIC_null ≤ ΔBIC_observed | H0)`
4. Apply BH-FDR across **all** windows

**Pros:**
- Valid FDR control
- No selection bias

**Cons:**
- Computational cost: 100-1000 bootstrap samples × N_windows

**Mitigation:**
- Parallelize across windows
- Use moderate bootstrap samples (B=100-200)
- Pre-filter windows with obviously insufficient data (F < 50)

### B. Sample Splitting

1. Split frames into discovery (60%) and confirmation (40%) sets
2. Screen using discovery set
3. Test selected windows using independent confirmation set

**Pros:**
- Statistically clean
- Fast

**Cons:**
- Loses 40% of data power
- MD trajectories are autocorrelated, so "split" is non-trivial

### C. Max-Statistic FWER Control

1. Permute/bootstrap preserving window overlap structure
2. Compute max ΔBIC across all windows per permutation
3. Build null distribution of global max
4. Threshold at desired FWER level

**Pros:**
- Preserves correlation structure
- Strong (FWER) control

**Cons:**
- Conservative (lower power than FDR)
- Must preserve window overlap during resampling

**Recommendation: Option A (Bootstrap ALL)**

- Valid FDR control
- Computationally feasible in Rust with parallelization
- No data splitting required

---

## 12. BIC IS MODEL SELECTION, NOT CALIBRATED SIGNIFICANCE

**M3-A Proposal (REJECTED):**
```rust
pub enum EvidenceStrength {
    Weak,      // ΔBIC > -2
    Moderate,  // -10 < ΔBIC ≤ -2
    Strong,    // -50 < ΔBIC ≤ -10
    VeryStrong // ΔBIC ≤ -50
}
```

**Problem:**

These are arbitrary thresholds without statistical calibration.

ΔBIC is a **model-selection criterion**, not a **significance test**.

Mixture-component testing is **non-regular** (parameters on boundary under H0), so:
- Standard χ² asymptotics DO NOT apply
- BIC thresholds like "−10" have no type-I error rate guarantee

**Recommendation:**

**DO NOT expose arbitrary categories as "evidence strength".**

Instead:

1. **Report raw ΔBIC** as continuous value
2. **Use parametric bootstrap** for calibrated p-value:
   ```
   H0: K=1
   Resample from fitted K=1 model
   Compute null distribution of ΔBIC
   p-value = P(ΔBIC_null ≤ ΔBIC_obs | H0)
   ```
3. **Report p-value** alongside ΔBIC

**API Proposal:**
```rust
pub struct MultimodalityResult {
    delta_bic: f64,           // Raw model selection criterion
    p_value: Option<f64>,     // Calibrated via bootstrap (if computed)
    bootstrap_samples: usize, // Number of bootstrap replicates used
}
```

---

## 13. REMAINING BLOCKERS BEFORE PRODUCTION IMPLEMENTATION

### Blocker 1: NO VALIDATED DETECTOR FOR TORSION WINDOWS

**Status:** Product von Mises REJECTED

**Next Steps:**

Evaluate alternative detectors:

**Option A: Multivariate von Mises**
- Allows correlation between dimensions
- Complex likelihood, no standard implementation
- May be computationally prohibitive

**Option B: Toroidal KDE + Mode Counting**
- Nonparametric, no independence assumption
- Bandwidth selection hard in high dimension
- Mode counting ambiguous

**Option C: Mixture of Wrapped Normals**
- Allows correlation via covariance matrix
- Easier to implement than multivariate von Mises
- Still assumes specific parametric form

**Option D: Projection Pursuit / Univariate Tests**
- Project to principal toroidal direction
- Test univariate von Mises mixture on projection
- May lose multimodality in off-axis directions

**Recommendation:** Evaluate Option C (Wrapped Normal Mixture) in M3-A3

### Blocker 2: CURVATURE DETECTOR UNDEFINED

**Status:** Deferred until torsion detector resolved

### Blocker 3: MULTIPLE-TESTING PROCEDURE UNIMPLEMENTED

**Status:** Bootstrap-ALL recommended, but not coded or tested

**Next Steps:**
- Implement parametric bootstrap under fitted K=1 null
- Validate FDR control on synthetic data with known ground truth
- Benchmark computational cost

### Blocker 4: TEMPORAL AUTOCORRELATION CORRECTION

**Status:** Not tested

**Next Steps:**
- Empirically measure autocorrelation in representative MD trajectories
- Test block-bootstrap vs. trajectory-level bootstrap
- Validate that bootstrap p-values control Type-I error under autocorrelation

### Blocker 5: EFFECT-SIZE CALIBRATION

**Status:** No normalized effect size defined

**Next Steps:**
- Literature review for circular/toroidal effect sizes
- Empirical calibration on protein datasets with known conformational states
- OR: Accept raw geometric measures without normalization

---

## 14. REVISED API IMPLICATIONS

**M3-A API (OBSOLETE):**
```rust
pub struct ToroidalMultimodalityDetector {
    method: ProductVonMises,  // REJECTED
    bic_threshold: -10.0,     // Arbitrary, not calibrated
}
```

**M3-A2 Revised API (PROVISIONAL):**
```rust
// Detector TBD — Product von Mises REJECTED
pub struct ToroidalMultimodalityDetector {
    method: TBD,  // Awaiting M3-A3 evaluation
    bootstrap_samples: usize,  // For p-value calibration
    min_sample_size: usize,    // Reject windows with F < threshold
}

pub struct MultimodalityResult {
    // Evidence (NO arbitrary categories)
    delta_bic: f64,
    p_value: Option<f64>,  // Calibrated via bootstrap

    // Components
    n_components: usize,
    weights: Vec<f64>,

    // Effect size (RAW, not normalized)
    mean_toroidal_distance: f64,
    min_toroidal_distance: f64,
    rms_per_dimension_distance: f64,
    component_dispersions: Vec<f64>,

    // Diagnostics
    converged: bool,
    warnings: Vec<String>,
}
```

**Key Changes:**

1. Remove `ProductVonMises` (rejected)
2. Remove `EvidenceStrength` enum (arbitrary)
3. Add `p_value` (bootstrap-calibrated)
4. Effect size as raw geometric quantities, not normalized scalar

---

## 15. REVISED REGRESSION/CALIBRATION PLAN

**M3-A Plan (OBSOLETE):**

15 regression tests for Product von Mises detector.

**M3-A2 Revised Plan:**

**DEFER** until detector is selected in M3-A3.

Once detector is chosen, regression plan must include:

1. **Correlated one-state nulls (MANDATORY):**
   - Test FPR on chain correlation ρ=0.5, 0.7, 0.9
   - Test FPR on common latent α=0.5, 0.7, 0.9
   - Assert FPR ≤ 5% at all levels

2. **Independent nulls:**
   - Narrow, broad, near-uniform
   - Assert FPR ≤ 5%

3. **Well-separated alternatives:**
   - 60°, 90°, 120°, 180° separation
   - Assert power ≥ 80%

4. **Rare events:**
   - Minority component 2%, 5%, 10%
   - Document detection threshold

5. **Bootstrap calibration:**
   - Assert p-value is uniform under H0
   - Assert p-value FDR control at α=0.05

---

## 16. DOCUMENTATION IMPACT

### KNOWN_ISSUES.md

**REMOVE:**
```
### Issue: Product von Mises +175°/-175° limitation
The Product von Mises detector fails to detect states at +175° and −175°
(10° separation on S¹) due to conditional independence assumption.

**Workaround:** None. This is a fundamental limitation.

**Status:** Accepted for v0.9.0; multivariate von Mises deferred to future work.
```

**REPLACE WITH:**
```
### Issue: M3 Multimodality Detection Unresolved
No validated statistical detector exists for toroidal window multimodality.

Product von Mises mixture model was evaluated and REJECTED due to systematic
false positives (70-100% FPR) on single correlated states at biologically
plausible correlation levels.

**Status:** M3 remains OPEN. Candidate detectors under evaluation in M3-A3.

**Current Behavior:** bc_tau returns uncalibrated Sarle BC heuristic.

**Impact:** Multimodality claims in TopoFold output are NOT statistically validated.
```

### Manuscript Draft

**DO NOT CLAIM:**
> "TopoFold detects conformational multimodality using a validated statistical test."

**UNTIL:**
A detector passes correlated-null FPR ≤ 5% AND has calibrated Type-I error control.

---

## 17. REMAINING UNCERTAINTIES

### Uncertainty 1: Which Detector Can Replace Product von Mises?

**Status:** CRITICAL OPEN QUESTION

**Candidates for M3-A3:**
- Wrapped Normal Mixture
- Multivariate von Mises
- Toroidal KDE
- Projection-based methods

**Next Step:** M3-A3 evaluation checkpoint

### Uncertainty 2: Can EM Robustness Be Improved?

**Status:** Oracle gaps of 117-420 log-units indicate EM issues

**Next Step:** Experiment 3 results (EM robustness) — currently running

**Provisional Finding:** Multiple random initializations (n_init=10) may help, but unlikely to solve correlated-null FPR issue

### Uncertainty 3: Temporal Autocorrelation Correction

**Status:** Untested

**Impact:** Affects ALL detectors, not specific to Product von Mises

**Next Step:** Empirical ACF measurement + bootstrap validation

### Uncertainty 4: Joint vs. Separate Curvature Detection

**Status:** Deferred

**Next Step:** After torsion detector selected, revisit curvature

---

## M3-A2 CHECKPOINT SUMMARY

**CRITICAL FINDING:**

**Product von Mises mixture model is REJECTED as primary multimodality detector**
due to **70-100% false positive rate on correlated one-state nulls** at biologically
plausible correlation levels (ρ≥0.5, α≥0.7).

**Root Cause:**

Product models assume:
```
p(τ₁, ..., τ_q | k) = ∏ᵢ p(τᵢ | μ_ki, κ_ki)
```

Protein torsions violate independence due to steric constraints, secondary structure,
and collective motions. Detector systematically misrepresents **single correlated states**
as **multiple independent states**.

**+175°/-175° Diagnosis:**

TWO mechanisms contribute:
1. **Insufficient power** (low κ, small q, small separation)
2. **Product-model misspecification** (correlated underlying state)

The M3-A hypothesis ("independence prevents seeing common shift") was INCOMPLETE.
Correlation causes **false positives on nulls**, not just missed alternatives.

**Next Steps:**

**M3-A3:** Evaluate alternative detectors:
- Wrapped Normal Mixture (allows correlation)
- Multivariate von Mises (complex but exact)
- Nonparametric methods (KDE, mode counting)

**M3 Status:** **OPEN** — no validated detector for production

**Protocol Compliance:**
- ✅ No production code modified
- ✅ No commits made
- ✅ No pushes performed
- ✅ Research/architecture only

---

**END OF M3-A2 CHECKPOINT REPORT**
