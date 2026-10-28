"""Writes plan.json and PLAN.md (corrected, frozen evaluation plan, v2) from evaluate.DESIGN and the pinned inputs. Run with the repository .venv."""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ev", HERE / "evaluate.py")
ev = importlib.util.module_from_spec(spec)
sys.modules["ev"] = ev
spec.loader.exec_module(ev)
me = ev.load_export()
ref = me.read_fasta(ev.EXPORT / "data/P08581.fasta")
esm = ev.read_esm(HERE / "data_raw/all_features_all_data_pre_binning.csv")
bl = me.read_blosum62(ev.EXPORT / "data/BLOSUM62.txt")
pairs = [(esm[k], bl[(ref[k[0] + 1058 - 1], k[1])]) for k in esm if ref[k[0] + 1058 - 1] != k[1]]
rho = float(spearmanr([p[0] for p in pairs], [p[1] for p in pairs])[0])
data = ev.assemble()
D = ev.DESIGN
plan = {
    "schema": "topofold.met-eval-plan/2", "date": "2026-10-28", "supersedes": "benchmarks/audit/met_feature_export_2026-10-27/EVALUATION_PLAN.md (v1, unfrozen; preserved)",
    "authorization": "The owner authorized, on 2026-10-28, a limited internal evaluation using the pinned public assay tables on the documented article (CC BY) and repository (MIT licence file) basis. The authorization does not establish commercial or redistribution rights; the source assay tables are not published here; nobody is contacted.",
    "target_population": "variants passing the measured no-drug eligibility rule (missense, numeric effects in both conditions, DMSO ROSACE.LOF == FALSE). The evaluation does not demonstrate prediction for variants whose baseline function is unknown.",
    "construct_and_endpoint": {"construct": "Ba/F3 TPR-MET wild-type (exon 14 retained) intracellular domain, IL-3 withdrawal; WT_rosace_effect_all.tsv", "gamma": "ROSACE_effects(inhibitor) - ROSACE_effects(DMSO) for the same variant (the article's gamma); log2 fitness scale; each inhibitor dosed at its own IC50, so inhibitors are modelled separately",
                             "inhibitors": D["drugs"]},
    "changes_from_v1": ["glycine side-chain distance: explicit deterministic replacement (L_min_sidechain_dist_A := L_min_heavy_dist_A, plus an indicator column); the exporter's no-imputation rule is unchanged, this is an analysis-side fill", "ridge alpha chosen by position-weighted validation loss over inner position blocks with the same buffer; ties go to the larger alpha; no usable inner fold -> alpha 10",
                        "empty or insufficient folds defined before scoring (below); the ESM comparison becomes part of the single pass and all six feature sets share identical rows", "GBM: no early stopping, no validation split, fixed seed", "'no improvement' replaced by 'criterion not met'", "family agreement and the engineering outcome fixed explicitly"],
    "rows": {"exclusion_classes": ["no eligible no-drug baseline", "structural block masked by the exporter (reason retained in the exporter coverage)", "mutation letter not one of the 20 residues (three variants labelled missense with mutation X)", "no ESM score (none among eligible rows at freeze time)"],
             "counts_inspected_before_freezing_no_model_fitted": {d: a["counts"] for d, a in data.items()}, "positions_per_drug": {d: len(set(a["pos"].tolist())) for d, a in data.items()}},
    "feature_sets": D["feature_sets"], "blocks": {b: ev.BLOCK_COLS[b] for b in ev.BLOCK_COLS},
    "esm_input": {"file": "Machine_Learning/data/all_features_all_data_pre_binning.csv (fraser-lab/MET_kinase_Inhibitor_DMS, commit c2052b97)", "bytes": (HERE / "data_raw/all_features_all_data_pre_binning.csv").stat().st_size, "sha256": ev.sha256_file(HERE / "data_raw/all_features_all_data_pre_binning.csv"),
                  "use": "only the pos_mut and score columns; score is the ESM-1b log-likelihood ratio the article says was obtained from the esm-variants web server; the same file holds fitness, resistance and the authors' structure-derived columns, none of which is read",
                  "checks": {"keys": len(esm), "possible_keys": 287 * 19, "score_invariant_across_the_nine_inhibitor_rows_per_key": True, "dms_position_range": [min(k[0] for k in esm), max(k[0] for k in esm)], "spearman_with_blosum62_over_keys_feature_only_plausibility": rho},
                  "limits": "provenance beyond the article statement and these consistency checks was not independently verified; the model behind the scores was not re-run"},
    "folds": {"outer": f"{D['outer_blocks']} contiguous blocks of positions (sizes 29 x 7 and 28 x 3), identical for all feature sets, inhibitors and models", "buffer_residues": D["buffer"], "inner": f"{D['inner_blocks']} contiguous blocks of the sorted training positions with the same buffer, ridge alpha only",
              "grouping": "all substitutions and all inhibitors at a position are in the same fold", "empty_or_insufficient_cells": f"a held-out block with no eligible row for an inhibitor is recorded as empty_test and contributes nothing; fewer than {D['min_train_rows']} training rows -> every feature set and model of that cell uses the training mean (flagged, identical for all sets); zero training rows -> the cell is flagged no_training_rows and the evaluation is reported as a gate failure with partial results; no inhibitor or fold is dropped silently. Inner folds with an empty validation set or fewer than {D['min_inner_train_rows']} training rows are skipped.",
              "limits": "blocking is in sequence only; residues close in space but far in sequence remain correlated; shared controls and libraries make inhibitors correlated"},
    "models": {"ridge": f"per-inhibitor ridge; standardization from training rows only (inside every outer and inner fit); alpha in {D['alphas']}", "gbm": f"per-inhibitor HistGradientBoostingRegressor({D['gbm']}); no validation split of any kind", "threads": "1 (environment set before numpy is imported)"},
    "metric": {"primary": "position-weighted mean squared error of gamma: per inhibitor, squared errors averaged per position, then over positions; then equal-weight mean over the six inhibitors; Delta = MSE(a) - MSE(b), negative favours a", "secondary": ["relative MSE reduction (MSE(b)-MSE(a))/MSE(b)", "per-block and per-inhibitor Delta", "per-inhibitor Pearson and Spearman correlation of out-of-fold predictions with gamma (pooled rows)"],
               "sensitivity": "same out-of-fold predictions re-scored on rows whose approximate sd_gamma is at or below the median (no refit)"},
    "comparisons": [{"id": c, "a": a, "b": b} for c, a, b in D["comparisons"]], "primary_question": "C4 (F3 vs F2): does the angular encoding beat the Cα distances computed from the same nine coordinates?  C3 (F3 vs F1): does it add beyond ordinary site context? All six comparisons are reported; none is selected afterwards.",
    "classification": {"criterion_met": "Delta < 0 AND the upper end of the 95% moving-block-bootstrap percentile interval < 0 AND Delta < 0 in at least 8 of 10 outer blocks AND in at least 5 of 6 inhibitors", "criterion_not_met": "Delta >= 0, or Delta < 0 in at most 5 of 10 blocks", "inconclusive": "otherwise",
                       "bootstrap": f"blocks of {D['bootstrap']['block_length']} consecutive positions, all inhibitors and models resampled together, {D['bootstrap']['resamples']} resamples, seed {D['bootstrap']['seed']}; descriptive, conditional on this dataset and procedure, coverage not assumed",
                       "no_equivalence": "an inconclusive interval or a small point difference is never read as equivalence"},
    "family_rule": {"both_met": "criterion met in both families", "one_met": "criterion met in one family only: family-dependent, not supported", "neither": "criterion not met in either family (not met or inconclusive)"},
    "engineering_outcome": {"threshold": "5% relative MSE reduction (provisional engineering continuation threshold; not demonstrated practical or commercial usefulness)", "continue_angular": "C3 criterion met in both families with relative reduction >= 5% in both, and C4 criterion met in both families", "continue_geometry": "C3 as above but C4 not met in both families",
                            "do_not_continue": "otherwise (not equivalence, not proof of no effect)"},
    "uncertainty_context": "approximate sd_gamma = sqrt(sd_inh^2 + sd_DMSO^2) assumes zero covariance between estimates that share controls; the covariance is unknown, so this is context and a row filter for one sensitivity analysis, not a measurement-error model and not a detectable-improvement bound",
    "external_baseline": "The authors report correlations between predicted and observed fitness scores. Our endpoint is gamma (inhibitor minus no-drug score) on no-drug-eligible variants, evaluated by position blocks with a buffer, per inhibitor, from wild-type-template features; theirs is a pooled model trained with predicted mutant-complex structure features and tested on held-out inhibitors, 20% of random positions and two amino acids. The endpoints and splits differ, so no numerical comparison or superiority ranking is made; the ESM arm (F4, F5) is the only external-model comparison, and their model was not recreated.",
    "resources": {"cpu_seconds_cap": D["cpu_seconds_cap"], "rss_mb_cap": D["rss_mb_cap"], "accounting": "RUSAGE_SELF plus RUSAGE_CHILDREN for the whole process; no child process is started; RLIMIT_CPU backstop; progress saved after every (inhibitor, fold) unit; a restart resumes with the frozen design and never recomputes completed units; budget exhaustion preserves partial results and ends the evaluation as partial"},
    "outputs": ["results/results.json", "results/oof_predictions.csv.gz (predictions only, no outcome values)", "results/folds.csv", "logs/run.log", "REPORT.md and two figures generated from these records"],
    "design": D,
}
(HERE / "plan.json").write_text(json.dumps(plan, indent=1, default=list))
md = f"""# MET exploratory evaluation: corrected plan v2 (generated from `plan.json`; frozen with the evaluator before scoring)

{plan['authorization']}

**Target population:** {plan['target_population']}

## Endpoint and rows
- {plan['construct_and_endpoint']['construct']}. Endpoint gamma = {plan['construct_and_endpoint']['gamma']}. Inhibitors: {', '.join(D['drugs'])}.
- Rows used per inhibitor (inspected before freezing; no model fitted): {'; '.join(f"{d} {a['counts']['used']} of {a['counts']['eligible']} eligible" for d, a in data.items())}. Exclusion classes: {'; '.join(plan['rows']['exclusion_classes'])}.
- Feature sets: {json.dumps(D['feature_sets'])}. Glycine side-chain distance: explicit deterministic replacement by the heavy-atom distance plus an indicator.
- ESM input: {plan['esm_input']['use']}. Checks: {json.dumps(plan['esm_input']['checks'])}. Limits: {plan['esm_input']['limits']}.

## Folds, models, metric
- Outer: {plan['folds']['outer']}; buffer {D['buffer']} residues; inner: {plan['folds']['inner']}; grouping: {plan['folds']['grouping']}. Cells: {plan['folds']['empty_or_insufficient_cells']} Limits: {plan['folds']['limits']}.
- Models: ridge: {plan['models']['ridge']}; GBM: {plan['models']['gbm']}.
- Primary metric: {plan['metric']['primary']}. Secondary: {'; '.join(plan['metric']['secondary'])}. Sensitivity: {plan['metric']['sensitivity']}.

## Comparisons and decisions
- {plan['primary_question']} Comparisons: {', '.join(f"{c[0]} = {c[1]} vs {c[2]}" for c in D['comparisons'])}.
- Criterion met: {plan['classification']['criterion_met']}. Not met: {plan['classification']['criterion_not_met']}. Otherwise inconclusive. {plan['classification']['no_equivalence']}. Bootstrap: {plan['classification']['bootstrap']}.
- Families: {plan['family_rule']['both_met']}; {plan['family_rule']['one_met']}; {plan['family_rule']['neither']}.
- Engineering outcome ({plan['engineering_outcome']['threshold']}): continue (angular): {plan['engineering_outcome']['continue_angular']}; continue (geometry): {plan['engineering_outcome']['continue_geometry']}; otherwise: {plan['engineering_outcome']['do_not_continue']}.
- Uncertainty context: {plan['uncertainty_context']}.
- External baseline: {plan['external_baseline']}

## Resources and recovery
{plan['resources']['accounting']}. Caps: {D['cpu_seconds_cap']} CPU-seconds, {D['rss_mb_cap']} MB resident memory.

Changes from v1: {'; '.join(plan['changes_from_v1'])}.
"""
(HERE / "PLAN.md").write_text(md)
print("plan v2 written; spearman esm~blosum", round(rho, 3))
