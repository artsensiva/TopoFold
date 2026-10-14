# C1 Feasibility Checks: RAS Switch-Region Retrieval (2026-10-13, revision 2)

**Status: BLOCKED.**

**Scope:** the three checks required by `docs/BIOLOGICAL_VALIDATION_DECISION_MEMO_2026-10-13.md` §5. Nothing was
scored, no adapter was implemented, no coordinates were downloaded, and nothing was staged or committed.

**Revision 2 (same day):**
- observed counts are separated from sampling-independence assumptions;
- structural label recovery is separated from biological function;
- assumption-dependent precision is separated from demonstrated power, and "attainable" and "predictable outcome"
  wording is removed;
- within-isoform counts are added;
- the limits of leave-one-crystal-form-out are stated;
- the choices still requiring preregistration are listed;
- the licensing scope is described without a legal conclusion;
- the request to the maintainers is drafted, not sent (`docs/C1_RASCORE_TERMS_REQUEST_DRAFT.md`).

**Inputs:** only the metadata already retrieved. Sources and SHA-256 values are in
`benchmarks/audit/c1_feasibility/results/c1_sources.json`, re-verified unchanged for this revision. Counts come from
code:
- `c1_counts.py` → `results/c1_counts.json`;
- `c1_groups.py` → `results/c1_groups.json`;
- `c1_isoform.py` → `results/c1_isoform.json`.

## 1. Label terms: what is established and what is unresolved

| Item | What the primary source states |
|---|---|
| PDB coordinates | wwPDB usage policy: archive data files "are available under the CC0 1.0 Universal (CC0 1.0) Public Domain Dedication" |
| Rascore code | Repository `mitch-parker/rascore` at commit `8f77dc2` (2024-03-19), the repository the paper cites: `LICENSE` is Apache License 2.0, and the README states "License: Apache License 2.0" |
| Rascore labels, repository copy (`entry.tsv`, cluster tables) | The files sit inside that repository. No statement in the repository addresses the data files specifically |
| Rascore labels, article copy | Parker et al. 2022: "All data are available as static tables in Supplementary Datasets S1-S4". The article's permissions statement says it is distributed under CC BY-NC-ND 4.0. The statement does not say whether it extends to the supplementary datasets |
| Rascore live database | No terms visible in the static page; it is a JavaScript application |

**Unresolved scope.** This is a description of what is not answered, not a legal conclusion.
- Whether the repository's licence is intended to cover the label data files.
- Whether the article's licence applies to Supplementary Datasets S1–S4.
- Which terms govern if the two copies differ.
- What terms apply to the live database.
- Whether publishing per-chain labels, splits or results derived from them is covered.
- Whether use within a project that has commercial objectives (roadmap §11) is covered.

These questions are addressed to the maintainers in the unsent draft request. The project should not assume an answer
to any of them.

## 2. Observed counts (facts about the metadata table)

**Table:** `entry.tsv` at commit `8f77dc2`.
- 767 chains in 442 entries;
- deposited 1990-04-30 to 2022-07-11;
- KRAS 472, HRAS 284, NRAS 11;
- X-ray 738, NMR 26, EM 2, neutron 1.

"Outlier" (complete loop, unclassified) and "Disordered" (incomplete loop) carry no class label.

| Loop | Class | Chains (KRAS / HRAS / NRAS) | Entries | Crystal forms |
|---|---|---|---|---|
| SW1 | Y32in.3P-ON | 277 (137 / 138 / 2) | 183 | 41 |
| SW1 | Y32out.2P-OFF | 193 (169 / 18 / 6) | 118 | 51 |
| SW1 | Y32out.0P-GEF | 23 (0 / 23 / 0) | 23 | 1 |
| SW1 | Outlier / Disordered | 156 / 118 | — | — |
| SW2 | 9 classes, 313 labelled chains | 12–105 per class | 5–91 | 1–26 (7 of 9 classes ≤ 7) |

**ON vs OFF, X-ray only, within each isoform** (`results/c1_isoform.json`):

| Isoform | ON chains / entries / CFs (largest CF) | OFF chains / entries / CFs (largest CF) | CFs with both classes |
|---|---|---|---|
| KRAS | 133 / 52 / 21 (**67**) | 166 / 100 / 43 (19) | 3 |
| HRAS | 137 / 125 / 20 (42) | 18 / 11 / 5 (8) | 1 |
| NRAS | 2 / 1 / 1 (2) | 6 / 4 / 4 (2) | 0 |

Across isoforms: 86 X-ray crystal forms, 4 of which contain both classes and 4 of which span isoforms.

**Composition within KRAS (observed):**
- **ON:** 57 of 133 chains are Q61H and 32 are G12D. 49 carry SP12-site ligands.
- **OFF:** 84 of 166 chains carry SP2-site ligands. 59 are one engineered construct (G12C, C51S, C80L, C118S).
- **Concentration:** one crystal form holds 67 of the 133 ON chains.

## 3. Sampling-independence assumptions (not facts)

Treating crystal forms, or any other group, as independent units is an **assumption**:
- Crystal forms can share laboratories, constructs, ligand series and refinement practices.
- The labels themselves were produced by one clustering procedure over the whole set.
- No metadata field here establishes independence.

**Leave-one-crystal-form-out (LOCFO)** removes same-lattice neighbours from the reference set. **It does not by itself
establish independence or prevent all leakage.** Remaining routes include:
- near-identical constructs, ligands or mutants in different crystal forms;
- the same laboratory or deposition series;
- NMR or other structures of the same complex;
- label assignment that used information from the whole set (DBSCAN plus nearest-neighbour assignment).

**Within-isoform comparison.** Restricting to KRAS removes the isoform imbalance: ON is about half HRAS overall,
while OFF is about 87 % KRAS. It would **reduce one confound** but leave others visible in §2:
- ligand occupancy (SP2 vs SP12);
- engineered constructs;
- dominant mutants;
- one crystal form holding about half of the ON chains.

HRAS alone has too few OFF chains (5 crystal forms). Within-isoform analysis is therefore a **candidate primary or
sensitivity design, not adopted here**; the choice is listed in §6.

## 4. Structural label recovery vs biological function

The labels are expert **structural** classes: a Y32 side-chain distance rule, nucleotide state, and SW1 backbone
dihedral clusters (residues 25–40, φ/ψ/ω, EDIA > 0.4).

| Label component | Cα-only observable? |
|---|---|
| Y32(OH)–G12(CA) < 10.5 Å | No (side-chain atom) |
| Nucleotide state (clustering done per state) | No |
| φ/ψ/ω cluster of residues 25–40 | Partly (backbone path only) |
| Loop completeness | Yes (missing Cα) |
| DBSCAN noise / outlier status | No |

**What a positive Cα-only result would support:** "SW1 Cα geometry recovers this expert structural classification
for structures outside the reference crystal forms".

**What it would not support:**
- activity, signalling state or effector binding;
- drug response;
- pocket presence;
- clinical relevance.

Rascore's association of classes with nucleotide or partner states is the authors' finding, not something this task
would test. Nucleotide-state analysis remains deferred.

## 5. Precision: assumption-dependent estimates, not demonstrated power

**Assumptions:**
- 86 crystal forms treated as independent units, one evaluation per unit;
- a paired method-disagreement rate p_d of 0.1, 0.2 or 0.3 (assumed, not measured);
- a normal approximation, 1.96·√(p_d/n).

**Result under those assumptions:** a 95 % half-width of 0.067, 0.095 or 0.116 for the difference in accuracy between
two methods.

**Not demonstrated:** power, which depends on unmeasured quantities:
- the true disagreement rate;
- within- and between-group correlation;
- class balance and weighting;
- the effective number of independent units, which may be well below 86 (§3).

These numbers indicate scale only. Whether a given margin δ could be resolved must be decided with explicitly stated
assumptions, and preferably with a simulation under several dependence scenarios, before any scoring.

**Margin δ (scientific, set separately from precision).** Proposed for review: δ = 0.05 balanced-accuracy points
against the strongest Cα-only baseline, justified from the triage use case. This is a judgement. It is not to be
derived from, or adjusted to, the precision estimates.

## 6. Choices still requiring preregistration (before any scoring)

1. **Label snapshot:** repository commit `8f77dc2`, article S1–S4 or a live export; chain/residue mapping (renumbered
   vs author numbering).
2. **Population:** which isoforms (all vs KRAS-only); X-ray only; resolution or quality filters; inclusion of
   engineered constructs, mutants and ligand-bound chains; treatment of the GEF class (excluded: one crystal form).
3. **Grouping unit and leakage rules:** crystal form plus entry; whether to add ligand series, construct or deposition
   group; handling of mixed-class and cross-isoform crystal forms.
4. **Input definition:** SW1 residues 25–40 Cα; missing-atom and altloc policy; model selection.
5. **Methods:** TopoFold distance(s) and settings; label-transfer rule (k, ties); reference-set composition.
6. **Baselines:** SW1 Cα RMSD (established library), Cα–Cα distances, and the Y32 CA–G12 CA distance proxy. The
   label rule itself (Y32 OH, nucleotide) as a separately reported stronger-information comparator.
7. **Estimand and metric:** balanced accuracy; chain vs crystal-form weighting.
8. **Inference:** non-inferiority margin δ and its justification; CI method (group-level bootstrap or similar) and
   dependence assumptions; minimum units; insufficient-evidence rule.
9. **Confound reporting:** stratification by isoform, ligand site, construct, mutation and deposition era.
10. **Secondary analyses:** temporal sensitivity and within-KRAS analysis, and the rule that they cannot override the
    primary result.

## 7. Feasibility summary

| Check | Observed | Interpretation |
|---|---|---|
| 1. Label terms | Code Apache-2.0; PDB CC0; label-data terms not stated for the repository copy; article CC BY-NC-ND 4.0 with unstated coverage of the datasets; live-database terms not visible | **Unresolved scope** (§1). C1 stays blocked |
| 2. Counts | SW1 ON 272 / OFF 190 X-ray chains in 86 CFs; KRAS ON 133 (21 CFs, one of 67 chains) / OFF 166 (43 CFs); HRAS OFF 18 (5 CFs); SW2 fragmented | Enough labelled chains exist for an SW1 two-class design. Independence of units is an assumption (§3) |
| 3. Label definitions | Side-chain OH rule, nucleotide state, φ/ψ/ω clusters with EDIA filtering | Cα-only sees only the backbone path. Claim limited to structural label recovery (§4) |
| Precision | 0.067–0.116 half-width under stated assumptions | Not a power demonstration. δ to be set scientifically and frozen first (§5) |

**Recommendation: BLOCKED** until the label-data terms are clarified by the maintainers (§1). No other task is
substituted.

Even if unblocked, scoring requires a frozen preregistration resolving the choices in §6.
