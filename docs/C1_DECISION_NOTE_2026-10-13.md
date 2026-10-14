# Decision Note: C1 (RAS SW1 Structural Label Recovery) — 2026-10-13

**Decision: C1 remains BLOCKED.** Nothing has been scored, implemented or downloaded beyond the metadata already
recorded. No other task is substituted.

**Why it is blocked.** The terms for using the Rascore label data are unresolved. This is a description of the gap,
not a legal conclusion:
- The code is Apache-2.0 and PDB coordinates are CC0.
- The label tables carry no data-specific statement in the repository.
- The article that hosts the same data is CC BY-NC-ND 4.0, with unstated coverage of its supplementary datasets.
- The live database shows no terms.

Internal evaluation, publication of derived labels or results, and commercial use are all unanswered. A request to the
maintainers is drafted (`C1_RASCORE_TERMS_REQUEST_DRAFT.md`) and **not sent**.

**What the metadata shows (observed):**
- SW1 ON vs OFF: 272 / 190 X-ray chains in 86 crystal forms.
- Within KRAS: ON 133 chains in 21 crystal forms, one of which holds 67 chains; OFF 166 chains in 43 crystal forms.
- HRAS has only 18 OFF chains (5 crystal forms). SW2 is too fragmented for a primary task.

**What remains assumed (not shown):**
- that crystal forms are independent units;
- that leave-one-crystal-form-out prevents leakage. It does not by itself: constructs, ligand series, laboratories
  and whole-set labelling can link units.

**Confounding.** ON vs OFF is confounded with isoform overall. Within KRAS, it is still confounded with ligand site,
engineered constructs, dominant mutants, and the concentration of ON chains in one crystal form. A within-KRAS
analysis would reduce one confound and is a candidate, not adopted.

**Scope of any future claim.** Recovery of an expert **structural** classification from Cα geometry. The labels
depend on a Y32 side-chain distance rule, nucleotide state and φ/ψ/ω clusters, which Cα-only input observes only
partly. Such a result would not show activity, signalling, drug response, pockets or clinical relevance.
Nucleotide-state analysis is deferred.

**Precision.** Under stated assumptions (86 independent units; method disagreement 10–30 %), the 95 % half-width of
a paired accuracy difference would be 0.067–0.116. This is not demonstrated power, and the effective number of units
may be smaller. The scientific margin (proposed: 0.05, for review) is set separately and is not tuned to these
numbers.

**Before any scoring, all of the following are required:**
1. written clarification of the label-data terms covering internal evaluation, publication of derived results or
   labels, and commercial context;
2. a frozen preregistration resolving the ten choices in `C1_FEASIBILITY_CHECKS_2026-10-13.md` §6, including the
   population (all isoforms vs KRAS-only), grouping and leakage rules, baselines, margin, CI method and the
   insufficient-evidence rule;
3. a precision or simulation analysis under explicit dependence assumptions, frozen alongside the margin.

**Owner actions:**
- review and, if appropriate, send the draft request;
- choose the margin and the population;
- decide whether a within-KRAS analysis should be primary or a sensitivity analysis.

**Engineering consequence (unchanged).** No adapter, loader or kernel work until the block is cleared and the
preregistration is frozen. Established RMSD libraries would be reused, and V2b remains an audit artifact.
