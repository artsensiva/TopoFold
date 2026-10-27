"""Writes results/abl_response_corrections_2026-10-27.json, CORRECTIONS_2026-10-27.md and EVALUATION_PROPOSAL_UNFROZEN_2026-10-27.md
from results/*.json (no inspected number is typed here; quotes and identifiers are literature or notice text).

    .venv/bin/python -I benchmarks/audit/abl_response_intake_2026-10-27/make_corrections.py
"""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
R = HERE / "results"
coh = json.loads((R / "endpoint_cohort_2026-10-27.json").read_text())
ins = json.loads((R / "inputs_inspection.json").read_text())
old = json.loads((R / "abl_response_intake_record.json").read_text())
P, S = coh["primary_candidate_after_exclusion"], coh["primary_sensitivity_subset_without_flagged"]
C = coh["cohorts"]
ab, sv = coh["abstract_vs_sheet"], coh["si_vs_github_sheet"]
t = ins["endpoint_table"]
dwt = coh["davis_wt_in_hauser_reproduction_nM"]
extra = [HERE / "data_raw/corrections_2026-10-27/hauser2018_SI_MOESM1.pdf", HERE / "data_raw/corrections_2026-10-27/crossref_schrock2013.json"]
extra_bytes = sum(p.stat().st_size for p in extra)
sha = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in extra}
prior = old["download_ledger"]["scientific_download_total_bytes_upper_bound"]
heavy = max(P["records_per_position"].items(), key=lambda kv: kv[1])
notices = [
    {"item": "Hauser et al. 2018 article", "notice": "'Open Access This article is licensed under a Creative Commons Attribution 4.0 International License ... The images or other third party material in this article are included in the article's Creative Commons license, unless indicated otherwise in a credit line to the material.'", "establishes": "CC BY 4.0 for the article and its included material", "does_not_establish": "terms for the GitHub files, the supplementary PDF or the third-party source values"},
    {"item": "Hauser Data availability statement", "notice": "'All relevant data are publicly available: compiled experimental datasets, input files for Prime and FEP+, and computational results that support our findings can be found at GitHub ...'", "establishes": "the authors made these files publicly available", "does_not_establish": "any licence or reuse terms"},
    {"item": "GitHub repository kehauser/Predicting-resistance-... (commit d5c6c508)", "notice": "licence field null in the API; no licence file in the 6-file tree; README: 'houses supplementary data files associated with the publication'", "establishes": "public availability", "does_not_establish": "reuse, redistribution or commercial terms"},
    {"item": "Official supplementary PDF (Springer ESM MOESM1)", "notice": "text searched for licence, copyright, creative, permission, reuse: no match", "establishes": "no notice in the file", "does_not_establish": "its terms (the article's licence statement does not mention it)"},
    {"item": "Schrock et al. 2013 abstract (Crossref metadata)", "notice": "no licence element in the Crossref record; the publisher page returned a challenge page to a non-browser client", "establishes": "a metadata abstract was read", "does_not_establish": "terms"},
    {"item": "Pemovska et al. 2015 and Davis et al. 2011", "notice": "not accessed (publisher supplementary URLs returned HTTP 403; articles not open access per Europe PMC)", "establishes": "nothing", "does_not_establish": "assay details, tables or terms"},
]
import math
RT = 0.0019872041 * 300.0
recs = {r['id']: r for r in coh['records']}
erl_wt, gef_wt = recs['erlotinib:Q252H']['wt_reference_nM'], recs['gefitinib:Q252H']['wt_reference_nM']
g = recs['ponatinib:G250E']
g_ic, g_ddg, g_wt = g['value'], g['ddg_sheet_kcal_mol'], coh['si_wt_ic50_nM']['ponatinib']
g_nm = g_wt * math.exp(g_ddg / RT)
t315 = ', '.join(f"{i.split(':')[0]} {recs[i]['value']:.0f}" for i in coh['primary_flagged_censored_or_bound_records'] if i.endswith('T315I'))
corrections = [
    {"id": "C1", "original": "144 variant-drug ΔΔG records presented as one endpoint with a binding-free-energy reading in places", "correction": f"Three separate cohorts by assay and source: {C['primary_candidate_bafs_viability']['records']} Ba/F3 cell-viability IC50 records for five inhibitors (Schrock 2013), {C['secondary_axitinib_cell_ic50']['records']} axitinib IC50 records (Pemovska 2015; assay not verified), {C['secondary_binding_davis_kd']['records']} erlotinib/gefitinib binding-constant records (Davis 2011; phosphorylated Abl; constructs not verified). Their transformed values are not pooled. Even within the cell-viability cohort ΔΔG is a log-potency shift expressed in kcal/mol at 300 K under the authors' assumptions, not a binding free energy."},
    {"id": "C2", "original": "wild-type IC50 'implied' from the sheet; 'the sheet gives no wild-type IC50'", "correction": f"Two different things are separated. (i) Reproducing the published transformation: the official SI states the wild-type IC50 for six inhibitors ({', '.join(f'{k} {v:.0f}' for k, v in coh['si_wt_ic50_nM'].items())} nM) and these agree with the values implied by the sheet (e.g. {ab['bosutinib']['wt_in_SI_nM']:.0f} vs {ab['bosutinib']['implied_wt_from_sheet_nM']} nM); SI and GitHub sheet agree on {sv['six_inhibitors_ic50_equal']}/{sv['six_inhibitors_records']} IC50 values and {sv['table1_ddg_equal_within_0.006']}/{sv['all_records']} ΔΔG values. For erlotinib and gefitinib the denominators (inferred {erl_wt:.1f} and {gef_wt:.1f} nM) equal the phosphorylated-Abl wild-type values in the Hauser reproduction of the Davis data ({dwt['erlotinib']['phosphorylated']:.0f} and {dwt['gefitinib']['phosphorylated']:.0f} nM), not the Davis paper itself. (ii) Independently verifying the original wild-type and mutant measurements: not done; only the abstract's stated ranges and the compilations were seen."},
    {"id": "C3", "original": "ponatinib G250E recorded as an unresolved discrepancy within the sheet", "correction": f"The same IC50 ({g_ic} nM) and ΔΔG ({g_ddg} kcal/mol) appear in the official SI, and {g_ic} nM lies outside the ponatinib range {ab['ponatinib']['abstract_range'][0]}-{ab['ponatinib']['abstract_range'][1]} nM stated in the Schrock abstract, while {g_ddg} kcal/mol corresponds to about {g_nm:.1f} nM given the wild-type value {g_wt:.0f} nM. The record is excluded from the primary cohort; its origin remains unresolved."},
    {"id": "C4", "original": f"'Only a large improvement could be detected' from a site-level AUROC standard error of {t['site_level_auroc_standard_error_if_true_auroc_0.8_hanley_mcneil']:.2f} at a true AUROC of 0.8", "correction": "Withdrawn. The calculation assumed a true AUROC of 0.8, a binary label per position and independent positions; none of these was the proposed prediction unit, endpoint or metric, so it does not establish what improvement is detectable, and no sample-count threshold is treated as evidence of adequate power."},
    {"id": "C5", "original": "'Internal analysis: conditional'", "correction": "Replaced by: unresolved. No source notice read states reuse terms for the data files; see the notice table. No legal conclusion is drawn and no permission is assumed."},
    {"id": "C6", "original": "two censored values; five inhibitors trace to a conference abstract that was not read", "correction": f"The abstract was read (Crossref metadata): assay 'engineered Ba/F3 cells, cell viability at 72 hours', 21 single mutants, five inhibitors, IC50 ranges that match the sheet (table below), and T315I reported only as '>4000' for the non-ponatinib inhibitors while the compilation gives numbers ({t315}). These four, plus imatinib L248R stored at 10,000 nM, are flagged as censored or bound-like and are treated in a sensitivity subset; whether the numbers are exact values or ceilings is unresolved."},
]
rec = {"schema": "topofold.abl-response-corrections/1", "date": "2026-10-27", "supersedes_statements_in": ["REPORT.md (intake of 2026-10-27, preserved unchanged)", "results/abl_response_intake_record.json (preserved unchanged)", "docs/DECISION_LOG.md D89"],
       "corrections": corrections, "permission_notices": notices,
       "permissions_status": {"internal_evaluation": "unresolved", "publication_of_derived_results": "unresolved for the data files (the article text is CC BY 4.0)", "redistribution_of_derived_files": "unresolved", "later_commercial_use": "unresolved"},
       "downloads": {"new_files": {p.name: {"bytes": p.stat().st_size, "sha256": sha[p.name]} for p in extra}, "new_bytes": extra_bytes, "cumulative_bytes_upper_bound": prior + extra_bytes, "correction_ceiling_bytes": 20_000_000,
                     "blocked_not_bypassed": ["Springer ESM for Pemovska 2015 and Davis 2011 (HTTP 403)", "ASH page for the Schrock abstract (challenge page)", "Zenodo Data S3 (832 MB, not downloaded)"]},
       "primary_candidate": {"definition": "Ba/F3 cell-viability IC50 (72 h) of the five Schrock inhibitors, 21 substitutions, minus the excluded record", "after_exclusion": P, "sensitivity_subset": S, "flagged": coh["primary_flagged_censored_or_bound_records"], "excluded": coh["excluded"]},
       "implementation_decision": {"data_dependent_feature_export": "not run", "reasons": ["reuse terms for the structure and endpoint files are unresolved", "original per-record values and wild-type measurements are not independently seen", "one record is an unresolved discrepancy and five are censored or bound-like"],
                                   "synthetic_only_tests": "none added: no builder exists yet; feature definitions to be fixed after the terms and cohort decisions"},
       "missing_evidence": ["a statement of reuse terms for the GitHub files, the supplementary PDF and the third-party values (owner decision; no one was contacted)", "the Schrock figure data (per-record Ba/F3 IC50 and the wild-type measurement) from the full source", "the origin of the ponatinib G250E value", "whether the T315I values are exact IC50 values or ceilings", "the Pemovska and Davis tables for the secondary cohorts"]}
(R / "abl_response_corrections_2026-10-27.json").write_text(json.dumps(rec, indent=1))


def table(rows, head):
    return "\n".join(["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(str(x) for x in r) + " |" for r in rows])


t_c = table([[k, v["records"], v["substitutions"], v["positions"], v["drugs"], v["source_groups"], v["resistant_records"], len(v["positions_with_resistant_record"])] for k, v in C.items()] +
            [["primary candidate after exclusion", P["records"], P["substitutions"], P["positions"], P["drugs"], P["source_groups"], P["resistant_records"], len(P["positions_with_resistant_record"])],
             ["primary sensitivity subset (flagged removed)", S["records"], S["substitutions"], S["positions"], S["drugs"], S["source_groups"], S["resistant_records"], len(S["positions_with_resistant_record"])]],
            ["cohort", "records", "substitutions", "positions", "drugs", "source groups", "resistant records (ΔΔG ≥ 1.36)", "positions with a resistant record"])
t_ab = table([[k, v["abstract_range"], v["sheet_min_max_excluding_T315I"], v["wt_in_SI_nM"], v["implied_wt_from_sheet_nM"]] for k, v in ab.items()], ["inhibitor", "IC50 range stated in the abstract (nM, excluding T315I)", "sheet min-max excluding T315I (nM)", "wild-type IC50 in SI (nM)", "wild-type implied by sheet (nM)"])
t_n = table([[n["item"], n["notice"], n["establishes"], n["does_not_establish"]] for n in notices], ["source", "notice read", "establishes", "does not establish"])
t_rec = table([[k, ", ".join(v)] for k, v in coh["records_with_unresolved_flags"].items()], ["record", "flags"])
t_pos = ", ".join(f"{k}: {v}" for k, v in P["records_per_position"].items())
md = f"""# Correction addendum to the Abl response intake (dated 2026-10-27)

Generated by `make_corrections.py` from `results/endpoint_cohort_2026-10-27.json` and the earlier inspection. The original `REPORT.md`, `results/abl_response_intake_record.json` and decision-log entry D89 are preserved unchanged; the statements corrected below are superseded by this addendum (D90). **No model was trained or scored; nothing was sent externally.** Evidence scope: assay-record structure and feasibility only.

## Corrections

{(chr(10) * 2).join(f"**{c['id']}.** Original: {c['original']}." + chr(10) + f"Corrected: {c['correction']}" for c in corrections)}
![source to endpoint](figures/abl_response_source_to_endpoint.png)

*Figure 1 (unchanged from the intake; reading corrected by this addendum). Purpose: trace records to sources. Read the three left boxes as three separate cohorts (cell viability, axitinib cell IC50 of unverified assay type, binding constants), and the central box as a compilation, not one homogeneous endpoint; the Cα-trace statement on the right is unchanged. Limitation: assay and source attributions for the axitinib and binding cohorts are unverified.*

## Cohorts and the proposed primary cohort

{t_c}

*Table 1. Purpose: count the meaningful units per cohort after exclusions. Source: `results/endpoint_cohort_2026-10-27.json`. Limitation: counts say nothing about adequacy.*

**Proposed primary cohort:** the Ba/F3 cell-viability (72 h) IC50 records of the five Schrock inhibitors (bosutinib, dasatinib, imatinib, nilotinib, ponatinib; templates 3UE4, 4XEY, 1OPJ, 3CS9, 3OXZ): {P['records']} records, {P['substitutions']} substitutions at {P['positions']} positions, {P['source_groups']} source group, after excluding ponatinib G250E. Per-record source, assay, wild-type reference, units, censoring flags and discrepancies are in the `records` list of the JSON. The axitinib and binding cohorts stay separate and are not analysed with it.

{t_ab}

*Table 2. The sheet agrees with the abstract's ranges, except ponatinib (the G250E value 0.021 nM) and the imatinib lower end, which the abstract's 201 matches to the wild-type value, not a mutant.*

Records with unresolved flags:

{t_rec}

## Is an exploratory position-grouped comparison interpretable?

Only as a pilot. The {P['records']} records form {P['positions']} position groups (records per position: {t_pos}); the largest group ({heavy[0]}, {heavy[1]} records) holds {100 * heavy[1] / P['records']:.0f}% of the records, and {len(P['positions_with_resistant_record'])} positions carry all {P['resistant_records']} resistant records ({S['resistant_records']} in the sensitivity subset). Records of one position share the residue and its template context across drugs, so they are strongly dependent; there is one assay, one laboratory, and the drug is perfectly confounded with its template. An exploratory result could show a gross failure (no consistent direction across positions) or a consistent direction that deserves confirmation; it could not support a statement of improvement, equivalence or generalization beyond this assay. No sample-count threshold is used as evidence of power.

## Terms: what the notices establish

{t_n}

The notices establish public availability and CC BY 4.0 for the article text. They do not state reuse terms for the data files. No legal conclusion is drawn; no permission is assumed. Statuses: internal evaluation, publication of derived results (data files), redistribution and commercial use are all recorded as **unresolved**.

## Decision on implementation

**The data-dependent feature export is not run.** Missing evidence, named exactly: {'; '.join(rec['missing_evidence'])}. No builder was written: the feature definitions should be fixed after the terms and cohort decisions, and a builder tested only on synthetic data would commit to a design before those decisions. The unfrozen evaluation proposal is in `EVALUATION_PROPOSAL_UNFROZEN_2026-10-27.md`.

Downloads in this correction: {extra_bytes / 1e3:.0f} kB (official supplementary PDF and the Crossref record of the abstract; ceiling 20 MB); cumulative {(prior + extra_bytes) / 1e6:.1f} MB. Publisher pages that returned HTTP 403 or a challenge page were not circumvented; the 832 MB Zenodo archive was not downloaded.
"""
(HERE / "CORRECTIONS_2026-10-27.md").write_text(md)
prop = f"""# Unfrozen evaluation proposal (2026-10-27): not a protocol, not frozen, not run

Generated by `make_corrections.py`. Everything below is a proposal for discussion. **No claim of independent validation is made or possible:** all {t['rows_with_full_values']} records, the {t['mutations_with_values']} substitutions and the eight templates were inspected during the intake, and the cohort definition itself (the exclusion of one record, the censoring flags) was made after seeing them.

## Question, unit, endpoint, metric
- **Generalization question:** for a residue position not used in training, does a structure-derived site descriptor of the wild-type complex add information about the cell-viability potency shift of a substitution, beyond the drug and the mutation properties?
- **Prediction unit:** a variant-drug record (substitution × inhibitor); the **grouping unit** is the residue position (all substitutions and drugs at a position are held out together).
- **Endpoint:** ΔΔG_cell = RT ln(IC50_mut/IC50_WT) per drug within the Ba/F3 cohort ({P['records']} records; sensitivity subset {S['records']} without the {len(coh['primary_flagged_censored_or_bound_records'])} censored or bound-like records), in kcal/mol at 300 K, labelled a log-potency shift and not a binding free energy. The wild-type denominator is the SI value per drug.
- **Primary metric (to be confirmed by the owner):** for each drug, the Spearman correlation between out-of-fold predictions and ΔΔG_cell across its records, averaged over the {P['drugs']} drugs; secondary: out-of-fold RMSE; resistant/susceptible classification at 1.36 kcal/mol is descriptive only ({len(P['positions_with_resistant_record'])} positions carry resistant records).

## Design
- **Outer evaluation:** leave-one-position-out over {P['positions']} positions; the model sees no record of the held-out position.
- **Training-only steps:** scaling, any hyperparameter choice and any imputation inside each training fold; tuning by inner position-grouped cross-validation on the training positions only; one fixed regularized linear model family, no feature selection, no use of the held-out fold for any choice.
- **Equal-information comparisons (paired on identical folds):** B0 drug mean from training positions; B1 B0 + mutation properties (volume, hydrophobicity, charge change, BLOSUM62); B2 B1 + wild-type ligand-site context (site-ligand distance, contact counts) from the template; B3a B2 + Cα window distances; B3b B2 + TopoFold κ/τ window of the same Cα window (same residues, same frames). The primary contrast is B3b against B2 and against B3a; B2 against B1 shows what site context adds.
- **Excluded inputs:** published Prime/FEP+ predictions, resistance annotations, hotspot lists, any feature chosen by looking at outcomes.
- **Uncertainty:** report per-position paired differences and a leave-one-position-out influence (jackknife) of the primary contrast, plus a paired bootstrap over the {P['positions']} positions as a descriptive interval whose coverage is not assumed. A conclusion needs a direction that is consistent across positions, not a small p-value or interval. The minimum useful improvement is to be set by the owner from decision utility, before any scoring.

## Confounding and limits
- **Template/drug confounding:** each drug has one template, so template-derived features are constant within a drug and differ across drugs; within the cohort they cannot be separated from drug chemistry, and the wild-type reference per drug removes only the drug's baseline potency.
- **Single source:** one assay and laboratory, so there is no between-source variation and no test of transfer to other assays; the axitinib and binding cohorts are not substitutes.
- **Dependence:** records at one position are dependent ({P['positions']} groups; the largest has {heavy[1]} records); the pooled metric must be position-weighted.
- **Reuse of exposed data:** the same records cannot serve as an untouched test. A claim beyond a pilot needs measurements not seen here.
- **Not testable:** no mutant-specific backbone effect (no mutant coordinates exist), no cross-kinase generalization.

## Preconditions before any implementation
{chr(10).join('- ' + m for m in rec['missing_evidence'])}
"""
(HERE / "EVALUATION_PROPOSAL_UNFROZEN_2026-10-27.md").write_text(prop)
print("wrote corrections JSON, CORRECTIONS and proposal")
