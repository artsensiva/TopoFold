# TopoFold Scientific Audit Artifacts

This directory contains scientific validation audit trails for TopoFold core methods.

**Purpose:** Preserve the evidence chain for methodological decisions, including failed hypotheses, implementation bugs, and superseded conclusions.

---

## Workstream index (reconciled 2026-10-13)

Each directory holds harness code, frozen-protocol references, results files and a report. **Harness code is not production code** (see below). Frozen protocols live in [`../../docs/`](../../docs/).

| Directory | Workstream | Data | Entry point | Status (short) |
|---|---|---|---|---|
| [`m3_multimodality/`](m3_multimodality/) | M3 Product-vM multimodality detector audit | SYNTHETIC | [`M3_STATUS.md`](m3_multimodality/M3_STATUS.md) | Authoritative status in `M3_STATUS.md`; not production-ready |
| [`angular_characterization/`](angular_characterization/) | WP-A angular characterization of the BC torsion scan | SYNTHETIC | [`REPORT.md`](angular_characterization/REPORT.md) | Characterization of one component |
| [`ordered_representation/`](ordered_representation/) | WP-C ordered angular representation | SYNTHETIC | [`REPORT.md`](ordered_representation/REPORT.md) | Candidate representation, fixtures only |
| [`adk_feasibility/`](adk_feasibility/) | WP-B AdK real-MD feasibility (G0–G1) | Third-party MD (CC BY 4.0) | [`REPORT.md`](adk_feasibility/REPORT.md) | Data integrity and reference assessment; development data |
| [`wpd_local_retrieval/`](wpd_local_retrieval/) | WP-D local retrieval: AdK development tier and mdCATH held-out evaluation | Third-party MD derivatives (CC BY 4.0) | [`REPORT_HELDOUT.md`](wpd_local_retrieval/REPORT_HELDOUT.md), [`REPORT_ADK.md`](wpd_local_retrieval/REPORT_ADK.md) | Held-out evaluation completed; see closeout addendum |
| [`legacy_engine/`](legacy_engine/) | LE-1 legacy κ–τ engine evaluation | AdK derivative + SYNTHETIC | [`REPORT.md`](legacy_engine/REPORT.md) | Legacy Fréchet coupling deprioritized (engineering decision) |
| [`eb1/`](eb1/) | EB-1 engineering benchmark, Amendments 1–2 | Development data | [`README.md`](eb1/README.md) | EB-1/A1 insufficient evidence; A2 decision engine (projection) |
| [`b2opt/`](b2opt/) | BQ-1 B2 single-query optimization (V0/V1/V2/V2b) | Development data | [`REPORT.md`](b2opt/REPORT.md) | Acceptance not met; V2b correct and faster on tested workloads |
| [`bq2/`](bq2/) | BQ-2 scaling confirmation of unchanged V2b | Development data + SYNTHETIC_tiled (performance-only) | [`REPORT.md`](bq2/REPORT.md) | Conclusion C (insufficient evidence) |
| [`c1_feasibility/`](c1_feasibility/) | C1 feasibility checks (RAS switch labels): terms, counts, label definitions | Public metadata only (no coordinates) | [`docs/C1_FEASIBILITY_CHECKS_2026-10-13.md`](../../docs/C1_FEASIBILITY_CHECKS_2026-10-13.md) | Blocked on label-data terms; SW1 two-class feasible on counts |
| [`cmp1/`](cmp1/) | CMP-1 bounded competitive evaluation (V2b vs MDTraj 1.11.2, MDAnalysis 2.10.0) | Development + SYNTHETIC_tiled | [`REPORT.md`](cmp1/REPORT.md) | No decisive advantage under frozen criteria; V2b integration not supported |
| [`bq3/`](bq3/) | BQ-3 direct-P2 qualification | SYNTHETIC_tiled (performance-only) | [`REPORT_SYNTHETIC_P2.md`](bq3/REPORT_SYNTHETIC_P2.md), [`README.md`](bq3/README.md) | Preflight stop (first attempt); direct run: prepared-database operation qualified, first query did not meet 10 s (excludes input loading) |
| [`dev_ca_pipeline/`](dev_ca_pipeline/) | DEV-CA-1 Cα pipeline development exercise (adapter → κ/τ → reports) | Third-party MD (CC BY 4.0), local AdK | [`README.md`](dev_ca_pipeline/README.md) | Development evidence: all 15 pre-specified checks pass; engineering only |
| [`dev_ens_compare/`](dev_ens_compare/) | DEV-ENS-1 descriptive per-window comparison, AdK open- vs closed-start (κ, τ, Cα distances; JS distance; time blocks) | Third-party MD (CC BY 4.0), local AdK | [`README.md`](dev_ens_compare/README.md) | Descriptive development comparison; one run per condition; no test or biological claim |
| [`dev_rep_compare/`](dev_rep_compare/) | DEV-REP-1 replicated descriptive comparison, mdCATH 1ha8A00 320 K vs 348 K (5 replicas each; κ, τ, Cα distances) | Third-party MD (CC BY 4.0), local development domain | [`INTERPRETATION_NOTE_2026-10-23.md`](dev_rep_compare/INTERPRETATION_NOTE_2026-10-23.md) | Closed: no window met the complete-separation criterion (194 channel–window evaluations); criterion sensitivity uncharacterised |
| [`synth_rep_char/`](synth_rep_char/) | SYNTHETIC_REPCHAR-1 sensitivity characterization of the replica-comparison rule (84 cells × 200 repetitions) | SYNTHETIC | [`README.md`](synth_rep_char/README.md) | Complete; conservative under the identical-distribution null, low sensitivity under temporal correlation, confounded by condition-shared sampling |
| [`synth_rep_perm/`](synth_rep_perm/) | SYNTHETIC_PERM-1 prospective replica permutation method vs complete-separation screen (81 cells x 200) | SYNTHETIC | [`README.md`](synth_rep_perm/README.md) | Complete; near-nominal null rejection, higher sensitivity; confounding unresolved; real-data evaluation needs its own design |
| [`synth_null_stress/`](synth_null_stress/) | SYNTHETIC_NULLSTRESS-1: identical marginals, different temporal correlation between conditions (6 cells x 200) | SYNTHETIC | [`README.md`](synth_null_stress/README.md) | Complete; permutation rejects 0.085–0.525 (nominal 0.048); screen almost never calls |
| [`synth_daboot/`](synth_daboot/) | SYNTHETIC_DABOOT-1: stationary-bootstrap Kuiper test of equal marginals under unequal temporal correlation (10 cells x 200) | SYNTHETIC | [`README.md`](synth_daboot/README.md), [`CLOSEOUT_2026-10-23.md`](synth_daboot/CLOSEOUT_2026-10-23.md) | Calibration failure (null rejection 0.095–0.230); positive controls 1.000; stopped as predeclared |
| [`dev_ensemble_report/`](dev_ensemble_report/) | Descriptive ensemble-comparison command demonstrated on mdCATH 1ha8A00 and AdK | Third-party MD (CC BY 4.0), local | [`README.md`](dev_ensemble_report/README.md) | Descriptive only; reproducible; agrees with DEV-REP-1/DEV-ENS-1; v0.3.0 exact exclusion causes (`results_v3/`) |
| [`pilot_review_2026-10-23/`](pilot_review_2026-10-23/) | Separate-agent usability review of the descriptive tool (archived evidence) | Development data (local mdCATH) | [`README.md`](pilot_review_2026-10-23/README.md) | Qualified pass; internal pilot closed; not independent human validation |
| [`abl_access/`](abl_access/) | ABL biological route: access status (trajectories, permissions, HDX data) | Public metadata; CC0 peptide map (git-ignored) | [`ACCESS_STATUS_2026-10-23.md`](abl_access/ACCESS_STATUS_2026-10-23.md) | Trajectories and processed HDX uptake not obtained; owner actions listed |
| [`kras_struct_cmp/`](kras_struct_cmp/) | KRAS-STRUCT-1: GDP vs GMPPNP wild-type KRAS crystal structures (both with MRTX1133), local Cα geometry vs RMSD | PDB (CC0), git-ignored local copies | [`README.md`](kras_struct_cmp/README.md) | Descriptive static comparison; differences at residues 26–40, 59–66, 102–110; confounded by crystal form; global RMSD independently verified (closeout) |
| [`ai_tester_1_trial_2026-10-24/`](ai_tester_1_trial_2026-10-24/) | AI tester 1: external AI technical trial of checkpoint 6ce7b49 (install, structpair, viewer) and bounded confirmation of 049a47f (`confirmation/`) | Tester's own trial directory (copies verified identical) | [`README.md`](ai_tester_1_trial_2026-10-24/README.md), [`REPORT.md`](ai_tester_1_trial_2026-10-24/REPORT.md) | Technical execution evidence only; not human or biological validation |
| [`abl_prereg/`](abl_prereg/) | ABL WT vs M472I preregistration support: input schemas, validator, HDX peptide → c-Abl mapping | PRIDE PXD034008 peptide map (CC0, git-ignored); UniProt P00519 (CC BY 4.0) | [`README.md`](abl_prereg/README.md) | Preparation only; no ABL data scored |
| [`kin_rob_1/`](kin_rob_1/) | KIN-ROB-1 (+ amendment `kin_rob_1a/`): robustness of named-atom extraction and kinase descriptors on 9 public kinase chains | PDB (CC0), UniProt (CC BY 4.0), git-ignored local copies | [`README.md`](kin_rob_1/README.md) | Engineering only; 7 pass, 1 pass with exclusions, 1 fail (spec numbering error caught by C1; amendment passes) |
| [`debuhr_abl_intake/`](debuhr_abl_intake/) | Intake: de Buhr & Gräter Abl MD (Dryad, CC0) | Dryad metadata, README, eLife text; archives later downloaded by the owner (local, git-ignored; D58, D66) | [`README.md`](debuhr_abl_intake/README.md) | conditional-go (access met; gates per study) |
| [`abl_myr_1/`](abl_myr_1/) | ABL-MYR-1: Myr vs noMyr kinase-domain ensembles, frozen spec | PCA.zip backbone XTC + two model archives (obtained 2026-10-24, local and git-ignored, hashes verified) | [`README.md`](abl_myr_1/README.md) | G0 passed; original run **stopped at G1** (`Set2/system.gro` absent in both model archives as published); nothing computed |
| [`abl_myr_1a/`](abl_myr_1a/) | ABL-MYR-1a: post-hoc amendment of ABL-MYR-1 — frozen initialization audit of the published Set1/Set2 start files (not the original frozen test) | model archives' start files only; no trajectory | [`README.md`](abl_myr_1a/README.md) | Verdict **V4** (Set1 original mdp/topology checks fail: Myr group names and water count differ); Set2 structures are copies of Set1's; no comparison run |
| [`abl_myr_1b/`](abl_myr_1b/) | ABL-MYR-1b: post-hoc comparability check of the 40 production run inputs (attempt 1 stopped at S3; amendments a1, a2), then the post-hoc trajectory comparison (run once; **closed 2026-10-25**, [`CLOSEOUT_2026-10-25.md`](abl_myr_1b/CLOSEOUT_2026-10-25.md); later lineage note: [`ADDENDUM_2026-10-25_closeout_lineage.md`](abl_myr_1b/ADDENDUM_2026-10-25_closeout_lineage.md)) | check: run inputs and mdp/ndx files; comparison: PCA.zip backbone XTC (40 runs) | [`README.md`](abl_myr_1b/README.md) | a1: K1 not established (MDAnalysis did not expose the stored box; not proof of incomparability); **a2 (GROMACS 2024.5): K1–K5 passed**; post-hoc ABL-MYR-1b protocol, independent AI review, dated amendment and wrapper frozen before `PCA.zip` was opened; **comparison run once, all gates passed: baselines cover TopoFold (post-hoc, descriptive), no TopoFold-only residue; αI not represented** |
| [`next_app_feasibility/`](next_app_feasibility/) | Next-application feasibility (Mpro prioritised; KRAS alternative): metadata only | Zenodo/RCSB/UniProt/S3 metadata; one 38 kB topology; XTC and HDF5 headers via range requests (no coordinates) | [`README.md`](next_app_feasibility/README.md) | Recommended conditional M1 (Moonshot apo/holo vs crystal fragment contacts); alternative K1 (KRAS G12D GDP/GTP); nothing scored or frozen (D67) |
| [`kras_k1/`](kras_k1/) | KRAS K1 (protocol v1 frozen 2026-10-25): G12D GDP vs GTP simulations against the external, already-exposed 5US4/5USJ crystal contrast, TopoFold vs Cα baselines at equal budget | Zenodo 10.5281/zenodo.1346073 (CC BY 4.0), S1 Table, XTC header fields; no coordinate read | [`PROTOCOL.md`](kras_k1/PROTOCOL.md) | Crystal endpoint computed once: gate failed, **K1_ENDPOINT_INADEQUATE, application closed**; no method was tested; no MD downloaded ([report](kras_k1/results/KRAS-K1_endpoint_report.md); D68–D70) |
| [`oncology_search_2026-10-25/`](oncology_search_2026-10-25/) | One bounded oncology evaluation search: Abl WT vs M472I (simulations + HDX-MS), KRAS mutant panel (simulations + biochemistry), KRAS Rascore labels; metadata and protocol feasibility only | Crossref, PRIDE and Europe PMC metadata and article text; no data files | [`README.md`](oncology_search_2026-10-25/README.md), [`docs/ONCOLOGY_EVALUATION_SEARCH_2026-10-25.md`](../../docs/ONCOLOGY_EVALUATION_SEARCH_2026-10-25.md) | **No candidate passes a feasibility gate today**; Abl route pending access; recommended next task is a product improvement (D72) |
| [`abl_hdx_lead_2026-10-25/`](abl_hdx_lead_2026-10-25/) | Inventory of the small public Abl HDX-MS files (PRIDE PXD034008 DynamX export; article Source Data and Supplementary Information) | PRIDE (CC0) and article supplementary files; no raw mass-spectrometry data | [`README.md`](abl_hdx_lead_2026-10-25/README.md), [`docs/ABL_HDX_PUBLIC_DATA_INVENTORY_2026-10-25.md`](../../docs/ABL_HDX_PUBLIC_DATA_INVENTORY_2026-10-25.md) | **Insufficient for an endpoint**: a six-column peptide list with no state, time, uptake or replicate (D73) |
| [`structpair_phospho_support_2026-10-26/`](structpair_phospho_support_2026-10-26/) | Opt-in SEP, TPO, PTR support (m14): design, PIM1 regression replay with derived specs, before/after coverage table | Local KIN-PAIR-1 inputs (no download), SYNTHETIC fixtures | [`REPORT.md`](structpair_phospho_support_2026-10-26/REPORT.md) | Engineering regression only; KIN-PAIR-1 unchanged (D85) |
| [`modified_residue_inventory_2026-10-26/`](modified_residue_inventory_2026-10-26/) | Modified-residue coverage inventory of the KIN-PAIR-1 inputs: adapter and correspondence omissions, affected windows | Local KIN-PAIR-1 files (hash-verified; no download), UniProt reference re-read and hash-checked | [`REPORT.md`](modified_residue_inventory_2026-10-26/REPORT.md) | Local collection only; no scoring; recommends a narrow support task (D84, m14) |
| [`structpair_numeric_record_2026-10-26/`](structpair_numeric_record_2026-10-26/) | Opt-in full-precision numeric record (m13b): BTK-2 original/transformed regression, computed → CSV → record table | Local KIN-PAIR-1 inputs (no download), SYNTHETIC fixtures | [`REPORT.md`](structpair_numeric_record_2026-10-26/REPORT.md) | Computational precision only; KIN-PAIR-1 not reclassified (D83) |
| [`structpair_altloc_tie_repair_2026-10-26/`](structpair_altloc_tie_repair_2026-10-26/) | Repair of KNOWN_ISSUES m13 (optional altloc-sensitivity ties): design, tests failing before and passing after, replay of the six KIN-PAIR-1 specifications, BTK-2 precision diagnosis | Local KIN-PAIR-1 inputs (no download), SYNTHETIC fixtures | [`REPORT.md`](structpair_altloc_tie_repair_2026-10-26/REPORT.md) | Engineering regression only; KIN-PAIR-1 outcome unchanged (D82) |
| [`kin_pair_1/`](kin_pair_1/) | KIN-PAIR-1: structpair on previously unused human kinases (EGFR, PIM1, BTK; 6 pairs); frozen protocol, independent references | RCSB (CC0), UniProt; about 8 MB | [`README.md`](kin_pair_1/README.md), [`report.md`](kin_pair_1/report.md) | Frozen-rule outcome IMPLEMENTATION FAILURE (4 of 6 pairs refused on altloc ties); engineering evidence only (D80, D81) |
| [`ai_trials_structpair_0_3_0_2026-10-25/`](ai_trials_structpair_0_3_0_2026-10-25/) | Two independent AI technical trials of structpair 0.3.0 (19b5268): evidence matrix, findings, hash inventory | Trial records (outside the repo; compact copy) | [`README.md`](ai_trials_structpair_0_3_0_2026-10-25/README.md) | AI trials, not human; 195 original files hashed and unchanged (D77) |
| [`structpair_viewer_folder_switch_2026-10-25/`](structpair_viewer_folder_switch_2026-10-25/) | Viewer focus after a Results-folder change and back: reproduction on 19b5268, fix, verification | Local result copies | [`README.md`](structpair_viewer_folder_switch_2026-10-25/README.md) | Reproduced 3 of 3 on the old code; fixed (D77) |
| [`structpair_advisory_correction_2026-10-25/`](structpair_advisory_correction_2026-10-25/) | Advisory output correction (structpair 0.3.0, schema /3): regression cases failing on 0.2.0, passing on 0.3.0 | Synthetic test inputs; local K1 records | [`README.md`](structpair_advisory_correction_2026-10-25/README.md) | 13 of 19 cases fail on 0.2.0 and pass on 0.3.0; 6 pass on both (D75) |
| [`structpair_abl_demo_2026-10-25/`](structpair_abl_demo_2026-10-25/) | structpair 0.3.0 on Abl: 1OPJ chains A, B vs 1OPK chain A (descriptive; partial within-entry status) | Local 1OPJ, 1OPK (PDB, CC0); RCSB entry metadata | [`README.md`](structpair_abl_demo_2026-10-25/README.md) | Descriptive only; no attribution to ligand, nucleotide, mutation or resistance (D76) |
| [`structpair_chains_demo/`](structpair_chains_demo/) | STRUCTPAIR-CHAINS-VERIFY-1: structpair 0.2.0 with several named chains, checked against the frozen KRAS-STRUCT-2 and K1 records and the four 0.1.0 demonstrations | Local KRAS structures (PDB, CC0) | [`README.md`](structpair_chains_demo/README.md) | Implementation consistency only: 0 mismatches in 4,740 + 1,224 values, 4 of 4 single-chain runs byte-identical (D74) |
| [`kin_map_1/`](kin_map_1/) | KIN-MAP-1: kinase residue-role mapping helper (`kinase_roles.py`) on the KIN-ROB-1 panel | KIN-ROB-1 local copies (no download); committed derived fixture | [`README.md`](kin_map_1/README.md) | Engineering only, not blind; 21/21 cases as frozen (19 accepted, 2 refused incl. the original PKA spec) |
| [`kinase_descriptors/`](kinase_descriptors/) | KIN-DESC-1: kinase-state descriptors vs Cα κ/τ, distances and local RMSD on four public Abl chains | PDB (CC0), git-ignored local copies | [`README.md`](kinase_descriptors/README.md) | Representation capability only; frozen before computing |
| [`mpro_dataset_search/`](mpro_dataset_search/) | Mpro MD dataset search: S3 metadata listing of the F@h COVID Moonshot set | AWS Open Data (CC0), metadata only | [`README.md`](mpro_dataset_search/README.md) | No trajectory downloaded |
| [`structpair_demo/`](structpair_demo/) | STRUCTPAIR-VERIFY-1: the two-structure command checked against KRAS-STRUCT-1/-2; fresh-clone smoke test | PDB (CC0), git-ignored local copies | [`README.md`](structpair_demo/README.md) | Engineering verification: exact reproduction; internal smoke test only |
| [`kras_struct_g12d/`](kras_struct_g12d/) | KRAS-STRUCT-2: GDP vs GMPPNP KRAS G12D crystal structures without a modelled pharmacological inhibitor (5US4/5USJ), all chain pairings | PDB (CC0), git-ignored local copies | [`README.md`](kras_struct_g12d/README.md) | Descriptive static comparison; residues ~29–39 differ in all 12 chain pairings and recur from STRUCT-1; 100–125 varies between crystal copies |

Large or third-party inputs that are not in Git, and how to retrieve or regenerate them: [`../../docs/CHECKPOINT_2026-10-13_PROVENANCE.md`](../../docs/CHECKPOINT_2026-10-13_PROVENANCE.md).

---

## Audit Principles

### Why Version Scientific Audit Artifacts?

**Transparency and Reproducibility:**
- Document how conclusions were reached, including errors discovered
- Preserve evidence for methodological decisions
- Enable independent verification of validation claims
- Support publication reproducibility requirements

**NOT a Clean Success Narrative:**
- Failed hypotheses are scientifically valuable (document what doesn't work)
- Implementation bugs and their discovery are part of the evidence
- Superseded conclusions remain in the historical record with supersession markers
- The messy scientific process is intentionally visible

### Distinction: Audit vs. Production Code

| Audit Artifacts | Production Code |
|----------------|-----------------|
| **Location:** `benchmarks/audit/` | **Location:** `src/`, `python/topofold/` |
| **Purpose:** Diagnostic, validation, evidence | **Purpose:** Production implementation |
| **Quality:** Research/experimental quality | **Quality:** Production quality |
| **Promotion:** Explicit approval required | **Promotion:** Not automatic from audit |

**IMPORTANT:** Audit artifact inclusion does NOT imply production readiness.

---

## Current Workstreams

### M3: Toroidal Multimodality Detection

> **Superseded summary (marker added 2026-10-13).** The status lines below date from M3-A2f (2026-10-03). Later checkpoints (M3-A2g) changed the classification. The authoritative current status is [`m3_multimodality/M3_STATUS.md`](m3_multimodality/M3_STATUS.md). The text below is kept for the record.

**Location:** [m3_multimodality/](m3_multimodality/)

**Current Status:** **M3-A2f CLOSED**

**Classification:**
```
VALID ALGORITHMIC HEURISTIC /
OPTIMIZATION ROBUSTNESS UNRESOLVED
```

**Product von Mises detector:**
- Exact-family validation: 0.13% (3/2250)
- K2 power (q=2, ≥60°): 100% (1620/1620)
- **PRIMARY BLOCKER:** 97% EM convergence failure
- **NOT production-ready**

See [m3_multimodality/M3_STATUS.md](m3_multimodality/M3_STATUS.md)

---

## Historical Checkpoints and Supersession

### What is a "Superseded" Checkpoint?

A checkpoint whose **conclusions were invalidated** by later work, but **retained for provenance**.

**Supersession Header Format:**
```markdown
SUPERSEDED SCIENTIFIC CONCLUSION

This report is retained for provenance.

Invalidated by: [checkpoint]
Reason: [one-line explanation]
Current status: see [STATUS.md]
```

### Why Preserve Failed Hypotheses?

1. **Document negative results**
2. **Prevent repeated mistakes**
3. **Methodological transparency**
4. **Complete evidence chain**

---

## Promotion to Production

**Audit artifacts are NOT automatically production code.**

Before promoting any audit implementation to `src/` or `python/topofold/`:

1. **Scientific validation complete** (all MANDATORY requirements met)
2. **External review and approval** (documented in ADR)
3. **Code quality standards** (production-quality implementation, tests, documentation)
4. **API design** (public interface, backwards compatibility)
5. **Integration testing** (works with existing TopoFold pipeline)

---

**Repository:** TopoFold
**Branch:** scientific-rebuild
**Last Updated:** 2026-10-03
