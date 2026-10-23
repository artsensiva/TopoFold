# DEV-PERM-1 Real-Data Evaluation Design — DRAFT, not frozen, not scored (2026-10-23)

> **DEFERRED, not scored (notice dated 2026-10-23).**
> - The replica permutation method (PERM-1) responds to differences in finite-trajectory histogram distributions,
>   including differences caused solely by temporal correlation (NULLSTRESS-1). It is not validated for isolating
>   differences in marginal geometry, so this design is deferred.
> - The descriptive comparison pipeline remains available.
> - Two statements below were corrected under the same notice: the meaning of a shared starting structure
>   (assumption 2) and the definition of exchangeability (assumption 1).

**Purpose:** one bounded development evaluation of the PERM-1 replica permutation method on real MD data, with a
single fixed global statistic, synchronized whole-replica permutations across windows, and descriptive residue maps.
- **Evidence level:** development evidence on an exposed development domain. It is not confirmatory, not a
  multiple-window significance procedure, and not biological validation.
- **Gating:** nothing is scored until this draft is reviewed, revised and frozen.

## Data and exposure

- **Domain:** mdCATH `1ha8A00` (local; CC BY 4.0; SHA-256 `e06d8cdb…`; read through
  `python/topofold_inputs/mdcath_bridge.py`).
- **Conditions:** 320 K vs 348 K, replicas 0–4 at each.
- **Frames:** 43–429, the same frames as DEV-REP-1. The initial-frame exclusion rule is inherited, not re-chosen.
- **Exposure:** WP-D schema check and runner validation, EB-1, CMP-1, and DEV-REP-1. DEV-REP-1 already computed
  per-window replica JS distances on exactly these frames; no window met the complete-separation criterion. The
  outcome of this evaluation is therefore read as development evidence about the method's behaviour on real data,
  not as an independent test.
- **Protected cohort:** the WP-D held-out cohort is not touched.

## Fixed analysis

- **Channels:** κ (49 three-residue windows) and τ (48 four-residue windows) are primary. Cα(i)–Cα(i+2) and
  Cα(i)–Cα(i+3) distances are the equal-information channels, reported separately.
- **Histograms and distances:** identical to DEV-REP-1 (36 bins; κ linear on [0, π]; distances linear with
  under/overflow; τ circular from −π; JS distance base 2).
- **Per-window statistic:** T_w = mean between-temperature replica JS distance − ½(within-320 mean + within-348 mean),
  as in PERM-1.
- **Single global statistic (deciding):** G = mean of T_w over the 97 κ and τ windows. The same quantity over the
  97 distance-channel windows, G_dist, is reported as a companion and does not decide.
- **Synchronized permutations:** each of the 252 labelled allocations of the ten whole replicas into two groups of
  five is applied identically to every window and channel. This keeps the dependence between windows within a
  replica.
  - The exact upper-tail probability of G uses the PERM-1 tie rule, α = 0.05, minimum attainable p = 2/252.
- **Residue maps (descriptive only):** per window, T_w,obs and the fraction of synchronized allocations with
  T_w ≥ T_w,obs. Both are plotted along the sequence for κ, τ and both distance channels. No per-window significance
  or FDR statement is made, and the map is not used to choose windows.
- **Companion:** DEV-REP-1's complete-separation categories, already recorded and not recomputed for decisions.

## Exchangeability assumptions (explicit)

1. **Exchangeability (definition):** the joint distribution of the quantities entering the statistic, here the ten
   replica-level histogram sets, is invariant under every permutation of replica labels across the two temperatures.
   That invariance is what the exact permutation probability requires.
   - Identical marginal geometry and identical dynamics in independent, identically distributed replicas is one
     sufficient condition. It is not a universally necessary one.
   - NULLSTRESS-1 shows that different temporal correlation alone can break it.
2. **Shared starting configuration:** all ten replicas start from one equilibrated configuration.
   - With independent stochastic seeds, replicas are conditionally independent given that start. A shared fixed start
     does not by itself make them dependent.
   - The limitation is **limited starting-state coverage**: every replica, and every conclusion, is conditional on
     that one start, and slow components may stay correlated with it within 430 ns.
3. Frame identity cannot be verified for the HDF5 coordinates, and time is inferred (1 ns per frame). Equal frame
   spacing at both temperatures is assumed.
4. Replica heterogeneity is assumed to be common to both temperatures. A replica unusual at one temperature only
   breaks exchangeability, as in PERM-1 S4a and S4c.

## How SYNTHETIC_NULLSTRESS-1 changes this design

With identical marginals and different temporal correlation, the permutation method rejected at 0.115–0.525
(φ 0 vs 0.8 / 0.95) and 0.085–0.135 (0.8 vs 0.95), well above the nominal 0.048. The screen made no consistent
call and one replica-dependent call in 1,200 draws. Different dynamics alone produce rejections, and 320 K and
348 K can be expected to differ in dynamics. Therefore:

- **Reading a rejection:** "replica-level histogram distributions differ between temperatures (marginal geometry
  and/or temporal structure)", never "local geometry differs".
- **Reading no rejection:** not evidence of equal geometry.
- **Predeclared diagnostics (reported, never used for selection or for re-running):** per temperature and channel,
  - the median within-temperature replica JS distance;
  - the median lag-1 autocorrelation of each window's feature series (κ and distances linear; τ via the circular
    analogue on sin/cos).

  If these differ between temperatures, the report states that dynamics could account for a rejection.
- **Not adopted here:** thinning frames to reduce correlation, or matching effective sample sizes. Either would need
  its own SYNTHETIC check first and is a separate decision.

## Outcomes and reporting

- G and its exact p, with the decision at α = 0.05 under the reading above.
- G_dist and its p, as a companion.
- The diagnostics.
- Residue maps for all four channels.
- Exposure, exchangeability and the limits above, stated in the report.

## Budget

- 10 replica feature sets already reproducible through the existing bridge and features pipeline.
- About 97 × 45 JS distances per channel group, then 252 vectorized allocations: seconds.
- Proposed cap: 300 CPU-s, with checkpoints `partial` until complete.

## Not in scope

- Scoring before the design is frozen.
- Other domains.
- The protected cohort.
- Biological interpretation.

ABL remains the preferred biological-validation route; other replicated datasets remain options.
