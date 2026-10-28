"""Writes eval_plan.json and EVALUATION_PLAN.md (the executable, unscored evaluation plan) from results/manifest.json,
results/features_structure.csv and the candidate-screen record. No assay table is read here: eligibility counts come from the screen record.
Run with the repository .venv."""
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
man = json.loads((HERE / "results/manifest.json").read_text())
scr = json.loads((HERE.parent / "dms_candidate_screen_2026-10-27/results/met_screen.json").read_text())
rows = list(csv.DictReader(open(HERE / "results/features_structure.csv")))
DR = [d["code"] for d in json.loads((HERE / "CONTRACT.json").read_text())["drugs"]]
ok = {d: {int(r["uniprot_position"]) for r in rows if r["drug"] == d and r["T_reason"] == "ok" and r["L_reason"] in ("ok", "no_side_chain")} for d in DR}
common = set.intersection(*ok.values())
tol = {d: scr["per_condition"][d]["paired_and_not_DMSO_LOF"] for d in DR}
plan = {
    "schema": "topofold.met-eval-plan/1", "status": "unfrozen plan; no outcome table has been read for it; scoring not started",
    "inputs": {"features": {k: v["sha256"] for k, v in man["tables"].items()}, "exporter_freeze": man["freeze"],
               "assay_table": {"file": "WT_rosace_effect_all.tsv", "repository": "fraser-lab/MET_kinase_Inhibitor_DMS", "commit": "c2052b97e12859e71cebd4ccd947c02e1084bc51", "sha256": scr["files"]["WT_rosace_effect_all.tsv"],
                               "terms": "article CC BY; repository MIT licence file; whether the MIT text governs the data tables is not spelled out; owner confirmation of this basis is required before the run; no commercial data rights are inferred"},
               "uncertainty_table": {"file": "met_scores_filtered.tsv", "sha256": "1a866db3cbc93f7d1d649d085c5e3a59e294c9cf464802583a0236a499660b86"}},
    "construct": "Ba/F3 cells expressing the wild-type (exon 14 retained) TPR-MET intracellular domain fusion, IL-3 withdrawal; the exon 14-skipped library (EX_ files) is not used",
    "endpoint": {"name": "gamma", "definition": "ROSACE_effects(inhibitor) - ROSACE_effects(DMSO) for the same variant (the article's gamma = beta_inh - beta_DMSO; log2 fitness scale; wild-type synonymous = 0; each inhibitor dosed at its own IC50, so drugs are modelled separately)",
                 "drugs": DR, "excluded_drugs": ["tivantinib (MET-independent activity per the article)", "capmatinib, glumetinib, cabozantinib, glesatinib analog (no experimental complex)"]},
    "eligibility": {"rule": "missense variant with numeric effects in both the inhibitor and the DMSO condition and ROSACE.LOF == FALSE in DMSO (the file's own flag); position inside 1059-1345; no outcome-based filtering otherwise", "expected_eligible_variants_per_drug_from_screen": tol},
    "missing_data": {"policy": "exporter masks are final; no imputation. Primary comparison on the common complete-case rows of each drug (L, T and C all available); the same rows and folds for every feature set. The only analysis-side fill: L_min_sidechain_dist_A at glycine sites takes L_min_heavy_dist_A, with an indicator column (documented, training-independent). Missingness (unresolved loops, engineered template residues) is not random: results apply to resolved regions only; coverage by region is reported",
                     "positions_with_all_blocks_per_drug": {d: len(v) for d, v in ok.items()}, "positions_with_all_blocks_in_every_drug": len(common)},
    "uncertainty_handling": {"posterior_sd": "approximate sd_gamma = sqrt(sd_inh^2 + sd_DMSO^2) from the filtered table, independence assumed although controls are shared; used only for a sensitivity restricted to variants with sd_gamma at or below the median, and as context for the MSE scale. It does not establish a smallest detectable improvement."},
    "feature_sets": {"F0": "mutation properties (14 numeric M columns) + drug baseline (per-drug model intercept)", "F1": "F0 + ligand/site context L (4 columns + glycine indicator)", "F2": "F1 + 36 Cα distances (C)", "F3": "F1 + TopoFold angular features T (7 kappa, 6 sin tau, 6 cos tau)",
                     "F4_optional_external": "F1 + ESM-1b log-likelihood ratio taken from the authors' pre-binning feature file (inputs only; their predictions and model outputs excluded); requires one pinned download of that file (about 7 MB) and verification of the column before use"},
    "comparisons": [{"id": "C1", "test": "F1 vs F0", "meaning": "what ordinary ligand/site context adds"}, {"id": "C2", "test": "F2 vs F1", "meaning": "Cα distances beyond context"},
                    {"id": "C3", "test": "F3 vs F1", "meaning": "TopoFold angular features beyond context (primary)"}, {"id": "C4", "test": "F3 vs F2", "meaning": "equal-coordinate head-to-head of the two geometry encodings (primary)"},
                    {"id": "C5_optional", "test": "F4+T vs F4", "meaning": "angular features beyond a language-model baseline"}],
    "folds": {"outer": "10 contiguous blocks of positions (29 or 28 each), identical for all feature sets, drugs and models", "buffer": "training excludes positions within 8 residues of the held-out block on both sides (windows span 4 residues each side, so windows overlap up to separation 8)",
              "inner": "5 contiguous blocks over the training positions with the same buffer, for ridge alpha only", "grouping": "all substitutions and all drugs at a position are held out together",
              "limits": "blocking is in sequence; residues close in space but distant in sequence stay correlated; shared DMSO control and libraries make drugs correlated, so the six drugs are not six independent replicates"},
    "models": {"ridge": "per-drug ridge, features standardized on training rows only, alpha in {0.1, 1, 10, 100, 1000} chosen by the inner blocks", "gbm": "per-drug HistGradientBoostingRegressor(max_depth=3, learning_rate=0.05, max_iter=200, min_samples_leaf=20, l2_regularization=1.0, random_state=0), no tuning",
               "no_other_models": True},
    "primary_metric": {"name": "mean squared error of gamma", "aggregation": "within a drug, squared errors averaged per position then over held-out positions with equal position weight; then equal-weight mean over the six drugs", "difference": "Delta = MSE(a) - MSE(b); negative favours a",
                       "secondary": ["relative MSE reduction", "per-drug Pearson and Spearman correlation (comparable in kind to the correlations the authors report, not on their splits)", "per-block and per-drug Delta"]},
    "decision_rule": {"consistent_improvement": "Delta < 0 AND the 95% percentile interval of a moving-block bootstrap (blocks of 20 consecutive positions, all drugs resampled together, 2000 resamples) excludes 0 AND Delta < 0 in at least 8 of 10 outer blocks AND in at least 5 of 6 drugs",
                      "no_improvement": "Delta >= 0, or Delta < 0 in at most 5 of 10 outer blocks", "inconclusive": "otherwise",
                      "claim_requires": "the same classification for both ridge and GBM; five pre-specified comparisons, no multiplicity correction claimed; intervals are descriptive (coverage not assumed)",
                      "practical_usefulness": "separate from the statistical statement: the owner sets a practical threshold on relative MSE reduction before scoring (proposed default 5%); a statistically consistent improvement below it is reported as not practically useful"},
    "external_benchmark": {"authors_model": "XGBoost on 14 interpretable features including ESM-1b LLR and structural features computed from Umol-predicted inhibitor-bound mutant complexes; trained with 8192 feature combinations; test = all variants of one type I and one type II inhibitor, 20% of positions and two amino acids",
                           "fair_and_cheap": ["arm F4 (ESM LLR as an input feature) against F1 under our protocol", "descriptive side-by-side of per-inhibitor correlations with the published values, labelled as different splits", "no recreation of their model or of Umol structures"],
                           "not_permitted": "any claim of superiority over the authors' model; their strongest competitor uses mutant-complex structure predictions that our wild-type-template features do not have"},
    "budget": {"cpu_minutes_max": 20, "ram_gb_max": 2, "downloads": "none required; optional one file of about 7 MB for arm F4", "estimate": "per-drug ridge and shallow GBM on a few thousand rows: about 2 CPU-minutes"},
    "acceptance_conditions": ["plan, evaluation code and tests hash-frozen before the first scoring run", "all input hashes match this plan; exporter freeze verified", "identical row sets and folds across all feature sets (asserted)", "buffer and grouping tests pass on synthetic data (no training row within 8 residues of a held-out block; all substitutions and drugs of a position in the same fold)",
                              "no outcome-dependent feature or row selection; scaling and tuning inside training folds only (tested with a leakage probe on synthetic data)", "results JSON contains every comparison, per-block and per-drug values, coverage and exclusions; report generated from it",
                              "stop and report if any gate fails; a negative or inconclusive outcome completes the task"],
    "commands_for_the_next_task": ["ulimit -t 1800   # 30 CPU-minute hard cap", ".venv/bin/python benchmarks/audit/met_eval_2026-10-28/evaluate.py freeze", ".venv/bin/python benchmarks/audit/met_eval_2026-10-28/evaluate.py run", ".venv/bin/python -m pytest tests/test_met_eval.py -q"],
}
(HERE / "eval_plan.json").write_text(json.dumps(plan, indent=1))
sha = hashlib.sha256((HERE / "eval_plan.json").read_bytes()).hexdigest()
md = f"""# MET exploratory evaluation: executable plan (unfrozen, unscored; generated by `make_plan.py` from `eval_plan.json`, sha256 `{sha[:16]}…`)

**Evidence scope:** a position-grouped, exploratory comparison of feature sets on public DMS data. It cannot support claims about mutation-induced structural change, transfer to other assays or kinases, or superiority over the authors' published model. No outcome table was read to write this plan.

## Construct, endpoint, eligibility
- {plan['construct']}.
- Endpoint γ = {plan['endpoint']['definition']}. Drugs: {', '.join(DR)}.
- Eligibility: {plan['eligibility']['rule']}. Expected eligible variants per drug (from the screen): {', '.join(f'{d} {n}' for d, n in tol.items())}.
- Uncertainty: {plan['uncertainty_handling']['posterior_sd']}
- Missing data: {plan['missing_data']['policy']}. Positions with all blocks available: {', '.join(f'{d} {n}' for d, n in plan['missing_data']['positions_with_all_blocks_per_drug'].items())} of 287; in every drug: {len(common)}.

## Feature sets and paired comparisons
{chr(10).join(f'- **{k}:** {v}' for k, v in plan['feature_sets'].items())}

{chr(10).join(f"- **{c['id']}** {c['test']}: {c['meaning']}" for c in plan['comparisons'])}

## Folds, models, metric, decision rule
- Outer: {plan['folds']['outer']}; buffer: {plan['folds']['buffer']}; inner: {plan['folds']['inner']}; grouping: {plan['folds']['grouping']}. Limits: {plan['folds']['limits']}.
- Models: ridge ({plan['models']['ridge']}); GBM ({plan['models']['gbm']}).
- Primary metric: {plan['primary_metric']['name']}; {plan['primary_metric']['aggregation']}; {plan['primary_metric']['difference']}. Secondary: {'; '.join(plan['primary_metric']['secondary'])}.
- Consistent improvement: {plan['decision_rule']['consistent_improvement']}. No improvement: {plan['decision_rule']['no_improvement']}. Otherwise inconclusive. {plan['decision_rule']['claim_requires']}.
- Practical usefulness: {plan['decision_rule']['practical_usefulness']}.

## External benchmark
Authors' model: {plan['external_benchmark']['authors_model']}. Fair and cheap: {'; '.join(plan['external_benchmark']['fair_and_cheap'])}. Not permitted: {plan['external_benchmark']['not_permitted']}.

## Budget, acceptance conditions, commands
Budget: {plan['budget']['cpu_minutes_max']} CPU-minutes, {plan['budget']['ram_gb_max']} GB RAM; downloads: {plan['budget']['downloads']}.
{chr(10).join('- ' + a for a in plan['acceptance_conditions'])}

Commands (the evaluation script is the next task):
```
{chr(10).join(plan['commands_for_the_next_task'])}
```
"""
(HERE / "EVALUATION_PLAN.md").write_text(md)
print("plan written", {d: len(v) for d, v in ok.items()}, len(common))
