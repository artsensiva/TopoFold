# M3: Circular/Toroidal Multimodality Detection

Scientific audit artifacts for M3 (multimodality detection on circular and toroidal geometry).

**Current Status:** M3-A2f CLOSED — [M3_STATUS.md](M3_STATUS.md)

**Current Classification:**
```
VALID ALGORITHMIC HEURISTIC /
OPTIMIZATION ROBUSTNESS UNRESOLVED
```

---

## Purpose

M3 validates statistical methods for detecting multimodal distributions on circular (S¹) and toroidal (S¹×...×S¹) geometries, targeting protein backbone dihedral angle distributions.

**Scientific Question:**
> Given a sample of points θ ∈ (-π, π]^q from a toroidal distribution, does the data support a unimodal (K=1) or multimodal (K≥2) hypothesis?

---

## Current Status Summary

### What Has Been Validated ✅

1. **Core Likelihood Implementation**
   - Fast Newton-refined κ MLE for von Mises distribution
   - Numerically stable log(I₀) using exponentially scaled Bessel functions
   - EM algorithm for Product von Mises K=2 mixtures
   - Canonical ΔBIC convention documented and verified

2. **Exact-Family Algorithmic Behavior**
   - Tested: 2250 Product-vM K=1 null samples
   - Result: 3/2250 = **0.13%** K=2 selections
   - Interpretation: Excellent algorithmic specificity under exact model

3. **K=2 Clear-Separation Recovery (q=2 only)**
   - Tested: Bimodal Product-vM with separations ≥60°
   - Result: **1620/1620 = 100%** recovery
   - Limitation: Only tested for q=2 (bivariate toroidal)

### Primary Blocker ⚠️

**Optimizer Convergence Robustness**
- Wrapped-normal σ≈1.0: 97% EM convergence failure
- Most K=2 selections are non-converged local optima, not valid detector decisions
- Cannot deploy until convergence rates characterized and initialization strategies validated

### Other Unresolved Requirements

- Unbounded likelihood requires explicit regularization strategy (U/C/P)
- Initialization robustness
- Extended K=2 power (q=1, q=5)
- **Dependent-toroidal robustness (MANDATORY for production)**
- Temporal correlation sensitivity
- Calibrated significance testing
- Overlapping-window multiple testing correction

**Product-vM Status:**
- NOT production-ready
- NOT validated screening detector
- NOT primary detector candidate

---

## Directory Layout

```
benchmarks/audit/m3_multimodality/
├── M3_STATUS.md              Current status & classification
├── MANIFEST.md               Provenance, checksums, regenerability
├── README.md                 This file
│
├── reports/                  Checkpoint reports & documentation
│   ├── M3_A2f_CHECKPOINT_REPORT.md         [CURRENT] Final validated state
│   ├── M3_A2f_CLOSURE.md                   Checkpoint closure document
│   ├── M3_A2f_LIKELIHOOD_DEGENERACY.md     Unbounded likelihood proof
│   ├── M3_A2e_CHECKPOINT_REPORT.md         [Superseded] Full benchmark
│   └── ... (historical checkpoints)
│
├── scripts/                  Diagnostic & reproducibility code
│   ├── product_vm_fast_accurate.py         [Definitive] Corrected fitter (U estimator)
│   ├── product_vm_three_estimators.py      [Definitive] U/C/P framework
│   ├── m3_a2f_analysis_framework.py        Corrected analysis tools
│   ├── investigate_wrapped_normal_anomaly.py  Convergence investigation
│   └── ... (benchmark and historical scripts)
│
└── results/                  Validated benchmark outputs
    ├── m3_a2e_checkpoint6_exact_family.csv     0.13% exact-family K2
    ├── m3_a2e_checkpoint9_k2_recovery.csv      100% K2 power (≥60°)
    └── m3_a2e_checkpoint10_wrapped_normal.csv  σ=1.0 anomaly
```

---

## Reproducing Key Results

### M3-A2e Exact-Family Validation (0.13%)

```bash
# From repository root
python benchmarks/audit/m3_multimodality/scripts/m3_a2e_full_benchmark.py
```

**Expected output:**
- `m3_a2e_checkpoint6_exact_family.csv` — 3/2250 K=2 selections
- **Runtime:** ~10-30 minutes

### M3-A2f Wrapped-Normal Convergence Investigation

```bash
python benchmarks/audit/m3_multimodality/scripts/investigate_wrapped_normal_anomaly.py
```

**Expected findings:**
- σ≈1.0: 33-50% raw K=2 selections
- Convergence rate: ~3% (97% failures)
- **Runtime:** ~5-10 minutes

---

## Toroidal Independence Warning

**Product von Mises assumes:**
```
p(θ₁, θ₂, ..., θ_q) = ∏ᵢ₌₁^q p(θᵢ)
```

This is a **strong simplifying assumption**. Real protein dihedral distributions may exhibit angular correlations.

**Validation under dependent-toroidal distributions is MANDATORY** before production deployment.

---

## Relationship to TopoFold Production Code

**IMPORTANT:** None of these audit artifacts have been promoted to production TopoFold code.

**Audit scripts are diagnostic code only.**

Product von Mises detector is:
- NOT in production pipeline
- NOT available as public API
- NOT recommended for use

---

**Last Updated:** 2026-10-03
**Audit Branch:** scientific-rebuild
**Promotion HEAD:** 95244fd486be1a96e29cc0cb37bb5cb244fa9407
