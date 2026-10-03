# M3-A2f: LIKELIHOOD DEGENERACY IN VON MISES MIXTURES

**Date:** 2026-10-03
**Phase:** M3-A (Product von Mises Validation)
**Checkpoint:** M3-A2f (Regularized/Constrained Estimator Definition)

---

## CRITICAL ISSUE: UNBOUNDED LIKELIHOOD

### The Problem

**Finite mixtures of von Mises distributions with unrestricted concentration parameters have UNBOUNDED likelihood.**

The unrestricted global maximum likelihood estimator (MLE) **DOES NOT EXIST**.

### Mechanism

A mixture component can achieve arbitrarily high likelihood by:

1. **Concentrating on a single observation** (or small subset)
2. **Weight shrinking:** π_k → small positive value
3. **Concentration diverging:** κ_k → ∞

For large κ, using the asymptotic expansion I₀(κ) ~ exp(κ)/√(2πκ):

```
f(θ | μ, κ) = exp(κ cos(θ - μ)) / (2π I₀(κ))
            ~ exp(κ cos(θ - μ)) · √(κ/(2π)) / exp(κ)
            = √(κ/(2π)) · exp(κ(cos(θ - μ) - 1))
```

At θ = μ (perfect alignment):
```
f(μ | μ, κ) ~ √(κ/(2π))  →  ∞  as κ → ∞
```

This allows the log-likelihood to increase without bound:

```
log L = Σᵢ log[Σₖ πₖ f(θᵢ | μₖ, κₖ)]

If component k concentrates on observation i:
  μₖ → θᵢ
  κₖ → ∞
  πₖ → 0⁺

The k-th component's contribution to observation i:
  πₖ f(θᵢ | μₖ, κₖ) → ∞

Even with πₖ → 0, the product can dominate.
```

### Mathematical Statement

For a K-component von Mises mixture on the circle:

**Theorem:** The likelihood function of a finite von Mises mixture is unbounded when concentration parameters are unrestricted.

**Proof:** Consider n observations {θ₁, ..., θₙ} and a K≥2 component mixture. Construct the following diverging sequence:

1. Fix component 1 center at observation θ₁: μ₁ = θ₁
2. Keep component 1 weight fixed: π₁ = ε > 0 (constant)
3. Let κ₁ → ∞
4. Keep remaining components (k=2,...,K) with positive density on all observations

The likelihood contribution from observation θ₁:

```
L₁ = Σₖ πₖ f(θ₁ | μₖ, κₖ)
   ≥ π₁ f(θ₁ | μ₁, κ₁)
   ~ ε · √(κ₁/(2π))    [using I₀(κ) ~ exp(κ)/√(2πκ)]
   → ∞   as κ₁ → ∞
```

Therefore log L(θ₁) → ∞, proving unboundedness.

**Alternative:** If π₁ → 0, require π₁√κ₁ → ∞ (not finite constant) for divergence.

The same mechanism applies to Product von Mises (each marginal can collapse independently).

---

## LITERATURE CITATIONS

### General Mixture Model Degeneracy

**Redner & Walker (1984)**
"Mixture densities, maximum likelihood and the EM algorithm"
*SIAM Review* 26(2): 195-239
DOI: 10.1137/1026034

> First formal recognition that Gaussian mixture likelihoods are unbounded when variance parameters are unrestricted. The same degeneracy mechanism applies to concentration in directional distributions.

**Kiefer & Wolfowitz (1956)**
"Consistency of the maximum likelihood estimator in the presence of infinitely many incidental parameters"
*Annals of Mathematical Statistics* 27(4): 887-906

> Establishes that MLEs may not exist when parameter spaces are unbounded.

### Directional Mixture Degeneracy

**Banerjee et al. (2005)**
"Clustering on the Unit Hypersphere using von Mises-Fisher Distributions"
*Journal of Machine Learning Research* 6: 1345-1382

> Notes (p. 1350): "The likelihood of a finite mixture of von Mises-Fisher distributions is unbounded... we constrain κ to lie in a finite interval [0, κ_max]."

**Hornik & Grün (2014)**
"movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions"
*Journal of Statistical Software* 58(10): 1-31
DOI: 10.18637/jss.v058.i10

> Section 2.2: "The likelihood... is unbounded as any concentration parameter can tend to infinity... In practice, we restrict κₖ to [0, κ_max]."

**Mardia & Jupp (2000)**
*Directional Statistics*
Wiley Series in Probability and Statistics

> Chapter 9: Discusses identifiability and degeneracy in circular mixture models.

### Product von Mises Specifically

The product von Mises (toroidal independence) inherits the same degeneracy:

```
f(θ₁, ..., θ_q | μ, κ) = ∏ᵢ₌₁ᵍ vM(θᵢ | μᵢ, κᵢ)
```

Each marginal can collapse independently, so the problem is at least as severe as the univariate case.

No published literature specifically addresses Product von Mises mixture degeneracy constraints, but the mechanism is identical to univariate von Mises.

---

## IMPLICATIONS FOR M3-A2e

### Current Implementation Status

The M3-A2e implementation uses:

```python
if kappa_max is None:
    kappa_max_actual = 1e6
```

**This is NOT unrestricted maximum likelihood.**

It is:
```
3-start local EM
with numerical concentration ceiling κ ≤ 10⁶
and ordinary BIC evaluated at the best local solution
```

### What M3-A2e Results Represent

The running M3-A2e benchmark is an **algorithmic/numerical benchmark**, not a validation of unrestricted MLE theory.

Its results are useful for:
- ✅ Numerical optimization behavior
- ✅ Empirical K=2 selection frequency under the numerical constraint
- ✅ Statistical power under the numerical constraint
- ✅ EM convergence stability

But NOT:
- ❌ Validation of unrestricted mixture MLE
- ❌ Theoretical statistical properties without regularization
- ❌ Production-ready estimator definition

### Required Next Steps

A scientifically valid production detector must explicitly choose one of:

#### A. **Constrained Likelihood**
Define a finite, justified parameter space:
```
0 ≤ κ_k ≤ κ_max (finite)
```

where κ_max has an interpretable relationship to:
- Angular resolution of the data
- Physical/biological concentration regimes
- Noise characteristics

NOT merely "large enough to rarely bind."

#### B. **Penalized Likelihood**
Use a literature-supported penalty:
```
Q(θ) = log L(θ) - λ·penalty(κ)
```

For example:
- Gamma prior on κ (Bayesian MAP)
- L2 penalty on log κ
- Regularized objective from directional clustering literature

#### C. **Alternative Estimator**
E.g., method-of-moments, robust estimator, or another approach with demonstrated consistency/regularization properties.

---

## M3-A2f ESTIMATOR DEFINITIONS

We will compare THREE estimators:

### **U: Unconstrained-EM Heuristic** (Current Implementation)
- 3 random initializations
- EM with κ ≤ 10⁶ numerical ceiling
- BIC selection at best local solution
- Label: "Numerical-ceiling EM heuristic"

### **C: Constrained Likelihood Estimator**
- Explicit finite concentration domain: 0 ≤ κ ≤ κ_max
- Multiple κ_max values tested: {20, 50, 100, 200, 500}
- Sensitivity analysis across bounds
- Justified selection based on scientific criteria

### **P: Penalized Likelihood Estimator**
- Literature-supported penalty function
- Separate tracking of:
  - Raw log-likelihood
  - Penalty term
  - Penalized objective
- Appropriate model-selection criterion (not ordinary BIC)

---

## MIXTURE MODEL SINGULARITY

Even after resolving concentration degeneracy, the K=1 vs K=2 comparison remains **non-regular/singular**:

### Sources of Non-Regularity

1. **Label switching:** Component indices are non-identifiable under K=2
2. **Parameter disappearance:** Under K=1, all K=2-specific parameters become unidentified
3. **Boundary behavior:** The null hypothesis (K=1) is on the boundary of the parameter space

### Implications for BIC

**Standard Schwarz BIC derivation assumes:**
- Regular parametric model
- Interior parameter values
- Identifiable parameters

**These assumptions are violated** in mixture model selection.

**Interpretation:**

✅ **BIC remains useful as an empirical model-selection heuristic**
❌ **BIC does NOT provide calibrated Type-I error control without bootstrap validation**

Distinguish:
- **Heuristic model selection** (what BIC provides)
- **Calibrated inferential test** (requires bootstrap/resampling)

The latter is planned for the M3-B barrier analysis phase.

---

## REFERENCES

1. Redner, R. A., & Walker, H. F. (1984). Mixture densities, maximum likelihood and the EM algorithm. *SIAM Review*, 26(2), 195-239.

2. Kiefer, J., & Wolfowitz, J. (1956). Consistency of the maximum likelihood estimator in the presence of infinitely many incidental parameters. *Annals of Mathematical Statistics*, 27(4), 887-906.

3. Banerjee, A., Dhillon, I. S., Ghosh, J., & Sra, S. (2005). Clustering on the unit hypersphere using von Mises-Fisher distributions. *Journal of Machine Learning Research*, 6, 1345-1382.

4. Hornik, K., & Grün, B. (2014). movMF: An R package for fitting mixtures of von Mises-Fisher distributions. *Journal of Statistical Software*, 58(10), 1-31.

5. Mardia, K. V., & Jupp, P. E. (2000). *Directional Statistics*. Wiley Series in Probability and Statistics.

6. McLachlan, G., & Peel, D. (2000). *Finite Mixture Models*. Wiley Series in Probability and Statistics. (Chapter 3: Identifiability and singularities)

7. Drton, M., & Plummer, M. (2017). A Bayesian information criterion for singular models. *Journal of the Royal Statistical Society: Series B*, 79(2), 323-380. (Modern treatment of BIC in non-regular settings)

---

## NEXT STEPS (M3-A2f)

1. ✅ Document degeneracy and literature (this file)
2. ⏳ Let M3-A2e benchmark complete
3. ⏳ Implement constrained estimator (C) with sensitivity analysis
4. ⏳ Implement penalized estimator (P) with literature-supported penalty
5. ⏳ Compare U vs C vs P on validation grids
6. ⏳ Corrected collapse analysis (ALL K2 fits, not just selected)
7. ⏳ Optimization convergence audit
8. ⏳ Initialization robustness testing
9. ⏳ Extended power analysis (q > 2)
10. ⏳ Dependent toroidal robustness (MANDATORY)
11. ⏳ Final classification (NO arbitrary thresholds)

---

**STATUS:** Degeneracy documented, awaiting M3-A2e completion and systematic estimator comparison.
