# Evidence Requirements for One Future Oncology Validation Task (outline, 2026-10-12)

**Status:** requirements outline only. No target, dataset or protein has been chosen. Nothing here is a claim of
oncology relevance, biological function or treatment efficacy.

**Revision 2026-10-12:**
- **Scope.** Clinical validation is **outside the current work package**. It is not excluded from the project's
  long-term mission; reaching it would need the clinical-tier evidence in §2, obtained by appropriate clinical
  research.
- **Holdouts must match the claim** (§3):
  - within-protein claims (new variants of a protein already studied) need held-out variants of that protein,
    simulated independently;
  - across-protein claims need held-out proteins, or sequence/structure families, never seen in development.
- **Baselines** (§4) are of two kinds, reported separately:
  - *equal-information baselines* receive exactly TopoFold's inputs (AGENTS.md rule 6) and test whether the method
    adds anything;
  - *stronger-information practical comparators* may use more (sequence conservation, structures, experimental
    annotations, all-atom features) and test whether the result is useful in practice.

**Current evidence base.** The WP-D held-out result establishes geometric local-conformation retrieval within
mdCATH. Everything below would need new, independent evidence.

## 1. Choosing the task without favouring TopoFold

1. **Start from the question, not the method.** State the oncology question and its intended use first, for
   example: "does a variant change the population of a local backbone state that is experimentally linked to
   activity?". Record who would use the answer and what decision it would inform.
2. **Eligibility criteria, frozen before candidates are enumerated:**
   - open, licensed data;
   - experimentally measured ground truth that is independent of any TopoFold output;
   - at least the minimum number of independent units (variants or proteins) set in §4;
   - MD that exists, or can be generated, under a documented protocol;
   - no prior TopoFold analysis of the candidate.
3. **Enumerate eligible candidates** from literature and databases, recording sources and identifiers. Verify every
   identifier (AGENTS.md rule 9).
4. **Select by a recorded rule that does not use expected TopoFold performance:** a seeded draw among eligible
   candidates, or a ranking on pre-declared criteria such as clinical relevance, label quality and data
   completeness. Freeze and hash the choice before any TopoFold or baseline run.
5. **Disclose** every candidate considered and why each was excluded.

## 2. Ground truth and evidence tiers

| Tier | Example evidence | What a positive result could support |
|---|---|---|
| Biophysical | Orthogonal experimental measurements of the relevant conformational state or equilibrium (structural, spectroscopic or biochemical), independent of the simulations | "TopoFold features track an experimentally measured conformational property" |
| Functional | Cellular or biochemical activity assays for the same variants | "…and that property relates to measured function" |
| Clinical | Patient outcomes from appropriate clinical studies | Outside this project's scope. **No treatment-efficacy claim can follow from tiers 1–2** |

- Labels must not be derived from TopoFold or from the simulations being evaluated (AGENTS.md rule 5).
- Label noise and assay limits must be stated.

## 3. Simulation and data requirements

- MD provenance:
  - force field, water model, ensemble and temperature;
  - multiple replicas per system, with **independent starting conditions** (unlike mdCATH);
  - equilibration and burn-in rules;
  - convergence diagnostics.
- Leakage control:
  - split by protein (and by sequence-similarity clusters), not by frame;
  - the held-out units are never used for development.
- Every synthetic or interpolated ensemble is labelled SYNTHETIC and excluded from validation evidence.

## 4. Comparators and analysis

**Baselines with equal information** (AGENTS.md rule 6), at minimum:
- the strongest geometric local baselines from WP-D (B2, B1, A);
- φ/ψ-based features;
- established MD analyses (e.g., dihedral PCA, Markov-state features);
- a non-dynamic baseline, e.g., a sequence- or structure-based variant-effect predictor.

The question is whether dynamics-derived features add value **beyond** these.

**Statistical plan, pre-registered:**
- primary endpoint and minimum useful effect, justified by the intended use;
- resampling unit = protein or variant;
- multiplicity control;
- calibration;
- a stated insufficient-evidence outcome.

**Sample size:** the minimum number of independent units is set from a power or precision calculation before data
access, not after.

## 5. Claims discipline

- Results are reported at the tier actually tested, with all methods and negative findings.
- Commercial or clinical language is excluded unless supported by evidence of that tier.
- The roadmap's commercial hypotheses (`ENGINE_MISSION_AND_ROADMAP.md` §11) remain hypotheses until such evidence
  exists.
