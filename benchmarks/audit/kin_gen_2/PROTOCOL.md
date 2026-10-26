# KIN-GEN-2: do Cα-only TopoFold features recover a defined kinase conformation better than equal-information baselines on previously unseen kinases? (protocol, written 2026-10-26)

**Replaces the unfrozen draft `benchmarks/audit/kin_gen_1/PROTOCOL_DRAFT.md` (preserved; it is superseded, not edited in substance).** Machine-readable parameters: `protocol.json`.
**Evidence scope: structural-label recovery only.** The labels are geometric groups of static public structures. Nothing here establishes function, drug response or biological validation.

## Why the draft was replaced
The draft depended on KinCoRe label files whose terms are unresolved, scored "AUROC of a representation" without a scoring procedure, selected best windows with a permutation control whose exchangeability was not justified, and mixed full-atom implementation checks with Cα-only comparisons.

## 1. Reference route (selected before any structure is scored): B, independent labels
`results/reference_source_check.json` (official pages re-read 2026-10-26) records: the KinCoRe website data (live resource, "Last updated 2026-05-25") states no licence or terms; the GitHub standalone software repositories have no licence file; the article is CC BY-NC-ND 4.0 and
that licence is not established for the website data or the supplementary datasets. **Route A (official labels) is therefore not used and nothing is assumed about its terms (unresolved).** The published definition names the distances D1 and D2 and the method (average-linkage clustering into three groups) but publishes no numeric boundaries.
Route B: D1 and D2 are computed here from public PDB coordinates with an independently written implementation (full-atom, uses the DFG-Phe Cζ), the three groups are fitted by average-linkage clustering **on development data only**, and evaluation entries take the nearest development group. The groups are called "D1D2 groups" (in, inter, out), never KinCoRe labels, never biological measurements.
No official label, table or code is copied or read. The full-atom rule defines the labels, so its agreement with them is 1 by construction; it is reported only as the reference and never as a performance claim.

## 2. Question and the one contrast
Primary task: **DFGin versus DFGout** (inter reported, excluded from scoring). Question: does the Cα-only TopoFold representation recover the D1D2 group of an unseen kinase's structure better than, or as well as, equal-information Cα distance and local-RMSD baselines? Spatial DFG classification is chosen because its definition is precise (two named distances) and the DFG window is local and mappable.

## 3. Inputs, eligibility and deterministic selection
- Development kinases (previously exposed, so they cannot count as unseen): ABL1, CDK2, MAPK14, AURKA, BRAF, PRKACA, the first 40 X-ray entries (≤ 2.5 Å, resolution then ID) per kinase plus three controls (1IEP, 1KV2, 1FIN).
- Evaluation kinases: the 23 declared oncology-relevant candidates of KIN-PAIR-1 minus every exposed accession, the KIN-PAIR-1 panel targets (EGFR, PIM1, BTK) and the development kinases, ranked by RCSB entry count, processed in rank order, first 30 entries (≤ 2.5 Å, resolution then ID) per kinase, until at least 15 entries labelled out from at least 3 kinases are collected, or 8 kinases or 90 MB are reached.
- One chain per entry (the alphabetically first chain of the target entity that passes every rule). Mapping by the existing `kinase_roles` helper against the UniProt canonical sequence (never by residue number), defaults, role mutations refused.
- Window: 9 residues, DFG-Asp − 3 to DFG-Asp + 5 (reference indices), all observed with Cα, consecutive in reference and deposited numbering, Cα–Cα ≤ 4.2 Å, standard names or MSE; D1 and D2 atoms observed (altloc blank, else A). Any other non-standard residue name in the window excludes the entry (recorded). Missing or excluded entries are listed with reasons.
- Anchors per kinase are fixed by rule (`protocol.json`): β3 Lys and HRD from UniProt features, DFG by the unique DFG in a UniProt-bounded window, αC Glu by a coordinate contact rule over that kinase's own mapped entries; a kinase without a unique anchor set is excluded.

## 4. Representations and the one scoring procedure
Representations from the same window: **T** (TopoFold core κ and τ, 19 values), **B1** (window-matched Cα distances, 13), **B2** (local RMSD), **B3** (all 36 Cα distances). Every representation is scored by the same procedure: k-nearest neighbours among development entries (distance as defined), k chosen on development
by leave-one-kinase-out balanced accuracy, decision "out if at least half of the neighbours are out". Scaling, k, group boundaries and the primary baseline (the best of B1, B2, B3 on development) are fixed from development data only. **No window or feature is selected, no direction or threshold is tuned on evaluation data.**

## 5. Metric, outcomes and inference
Primary: weighted balanced accuracy (each kinase × class cell equal weight within its class); secondary: pooled balanced accuracy, per-kinase counts, AUROC of the neighbour score, and the subset of kinases below the family-flag identity (kinase-domain identity to any development kinase < 0.40; related kinases are flagged, not removed).
Outcomes are those of `protocol.json` (insufficient evidence; comparable; advantage; baseline advantage; unclear), plus a separate Cα-only limitation statement. **No p-value, interval or permutation test**: the independent groups (kinases) are too few and exchangeability is not justified; results are descriptive.

## 6. Leakage controls and correctness gates (before evaluation scoring)
Disjoint kinases and entries (G7); development-only scaling, k, boundaries (by construction and asserted); gates G1–G8 in `protocol.json` (core κ/τ against numpy, rigid-motion invariance, Kabsch against MDAnalysis QCP, role identities and a shifted-mapping control, the three control structures, group sizes, disjointness, D1/D2 against an MDAnalysis recomputation). Any failed gate stops the study with that gate named.

## 7. Resources and freezes
Evaluation: ≤ 100 MB of scientific downloads (cap 90 MB in the stopping rule) and ≤ 2 CPU-hours; preparation (development downloads and fitting) is recorded separately. Freeze 1: this protocol and the source check (before any structure is downloaded). Freeze 2: the analysis code, the development model (groups, scaling, k, primary baseline) and the evaluation-selection code, before any evaluation structure is selected or scored.
A negative or insufficient-evidence outcome completes the work package. Earlier studies are not reopened.

## Amendments

- **Amendment 1 (2026-10-26), made during development, after the first development preparation run and before any evaluation structure was selected, downloaded or scored** (`amendment_1.json`; protocol.json is unchanged and stays hashed by freeze 1).
  The first run (log kept as `results/dev_prep_run1.log`) failed gate G5 and G2 and excluded MAPK14, AURKA and CDK2 for anchor reasons: UniProt active-site annotations are not reliable across kinases (MAPK14's lists the DFG Asp, with ILD in the reference at the HRD position; AURKA has two ATP-binding Lys),
  and the CDK2 contact rule fails because almost all CDK2 entries are αC-out. Changes: (1) G2 allows 1e-9 for the float64 representations and 1e-4 for T (the core receives float32); (2) anchors: the β3 Lys is the unique ATP-binding-site Lys at UniProt kinase-domain start + 20 to + 40, HRD and DFG are unique motif hits in windows
  measured from the β3 Lys (+70 to +125 and +95 to +140), the UniProt active site is no longer used; the windows come from the development kinases' spacing and are applied unchanged to evaluation kinases; (3) controls: out = 1IEP, 1KV2, 1UWH; in = at least 90 % of the PRKACA development entries. Nothing about evaluation kinases, scoring, metric or outcomes changes.
- **Amendment 2 (2026-10-26), same stage (second development preparation run, before any evaluation structure was selected, downloaded or scored)** (`amendment_2.json`). After amendment 1 the preparation gave 116 development entries from five kinases (CDK2 yielded none: almost all its entries are αC-out, so the contact rule finds no unique αC Glu;
  its 1FIN control is therefore absent) and gate G5 failed only because the control 1KV2 (p38) is excluded by the window rule (a chain break in its DFG region). Change: a control that fails the eligibility rules is reported as "not testable" and does not fail G5; G5 passes if every testable out-control (1IEP, 1KV2, 1UWH) is assigned out, at least two are testable, and at least 90 % of the PRKACA development entries are assigned in.
  Nothing about evaluation kinases, scoring, metric or outcomes changes. Run-1 and run-2 development records are kept (`results/dev_model_run1.json`, `results/dev_prep_run1.log`).

