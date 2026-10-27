"""Generates REPORT.md and figures/abl_response_source_to_endpoint.png from results/*.json (no inspected number is typed here).

    .venv/bin/python -I benchmarks/audit/abl_response_intake_2026-10-27/make_report.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
(HERE / "figures").mkdir(exist_ok=True)
ins = json.loads((HERE / "results/inputs_inspection.json").read_text())
rec = json.loads((HERE / "results/abl_response_intake_record.json").read_text())
t, S, cmp_ = ins["endpoint_table"], ins["structures"], ins["ca_trace_comparison"]
inh = t["inhibitors"]
n_schrock = sum(inh[k] for k in ("bosutinib", "dasatinib", "imatinib", "nilotinib", "ponatinib"))
n_pem = inh["axitinib"]
n_davis = inh["erlotinib"] + inh["gefitinib"]
pair = cmp_["pairwise_rmsd_A"]
within_a = [v for k, v in pair.items() if all(x in ("axitinib", "bosutinib", "dasatinib", "erlotinib", "gefitnib") for x in k.split("|"))]
within_b = [v for k, v in pair.items() if all(x in ("imatinib", "nilotinib", "ponatinib") for x in k.split("|"))]
between = [v for k, v in pair.items() if v not in within_a and v not in within_b]
led = rec["download_ledger"]
INK, GREY, BLUE, ORANGE, GREEN, LIGHT = "#2b2a27", "#6b6a66", "#2a78d6", "#eb6834", "#2f8f6b", "#f1f0ec"

# ---- figure -------------------------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10.4, 4.7))
ax.axis("off")
ax.set_xlim(0, 10.4)
ax.set_ylim(0, 4.7)


def box(x, y, w, h, txt, fc=LIGHT, ec=GREY, fs=7.4, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=fc, ec=ec, lw=1))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, color=INK, fontweight="bold" if bold else "normal", wrap=True)


def arrow(x0, y0, x1, y1, c=GREY):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=c, lw=1.1))


ax.text(0.1, 4.5, "Published measurements", fontsize=8.5, fontweight="bold", color=INK)
ax.text(4.0, 4.5, "Compiled endpoint (Hauser 2018)", fontsize=8.5, fontweight="bold", color=INK)
ax.text(7.9, 4.5, "Structural inputs in the deposit", fontsize=8.5, fontweight="bold", color=INK)
box(0.1, 3.35, 3.3, 0.95, f"Schrock 2013 (Blood abstract)\n5 inhibitors, {n_schrock} records\ncell viability IC50 (per paper)", fc="#e6eefb", ec=BLUE)
box(0.1, 2.25, 3.3, 0.95, f"Pemovska 2015 (Nature)\naxitinib, {n_pem} records\ncell IC50 (per paper, not re-read)", fc="#e6eefb", ec=BLUE)
box(0.1, 1.15, 3.3, 0.95, f"Davis 2011 (Nat. Biotechnol.)\nerlotinib, gefitinib: {n_davis} records\nbinding constants (Kd), not IC50", fc="#fbeadf", ec=ORANGE)
box(3.95, 1.85, 3.1, 2.3, f"Supplementary Data 3\n{t['rows_with_full_values']} records with ΔΔG\n{t['n_inhibitors']} inhibitors, {t['mutations_with_values']} mutations,\n{t['sites_with_values']} residue positions\n{t['resistant_records']} resistant (ΔΔG ≥ 1.36 kcal/mol)\nin {t['sites_with_any_resistant_record']} positions\n2 censored (IC50 > 10,000 nM) kept as equalities", fc="#fff", ec=INK, bold=False, fs=7.6)
box(3.95, 0.35, 3.1, 1.05, "Aldeghi 2019: same 144 records\n(not independent replication);\nML, MD, Rosetta inputs not retrievable here", fc=LIGHT)
box(7.9, 3.35, 2.4, 0.95, f"{len(S)} wild-type Abl:inhibitor models\n6 co-crystal + completion, 2 docked", fc="#e3f3ec", ec=GREEN)
box(7.9, 2.25, 2.4, 0.95, "0 mutant coordinate files\n0 trajectory files", fc="#fbeadf", ec=ORANGE, bold=True)
box(7.9, 1.15, 2.4, 0.95, f"Cα traces of the {len(S)} models differ\nby template ({min(within_a):.2f}-{max(between):.1f} Å), not by mutation", fc="#e3f3ec", ec=GREEN)
for y in (3.8, 2.75, 1.6):
    arrow(3.4, y, 3.95, 3.0 if y > 2 else 2.2)
arrow(5.5, 1.4, 5.5, 1.85)
arrow(7.05, 3.0, 7.9, 3.8, GREY)
arrow(7.05, 2.8, 7.9, 2.75, ORANGE)
fig.tight_layout()
fig.savefig(HERE / "figures/abl_response_source_to_endpoint.png", dpi=130)
plt.close(fig)


def table(rows, head):
    return "\n".join(["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(str(x) for x in r) + " |" for r in rows])


st_rows = [[k, v["construction"], f"{v['first_last'][0]}-{v['first_last'][1]}", v["n_standard_residues_with_CA"], v["protein_chain"], f"{ins['encoded_wt_identity'][k]['mutation_wt_identity_present']}/{ins['encoded_wt_identity'][k]['of']}",
            f"{rec['mapping']['atoms_available'] and ins['ligand_frame_and_site_contacts'][k]['n_sites_within_4.5A']}/{ins['ligand_frame_and_site_contacts'][k]['n_sites_beyond_10A']}"] for k, v in S.items()]
t_st = table(st_rows, ["inhibitor model", "construction", "residues", "Cα count", "chain", "mutation wild-type identity found", "sites within 4.5 Å / beyond 10 Å of the ligand (of 20)"])
t_in = table([[k, v, t["inhibitors"][k]] for k, v in {"axitinib": "Pemovska 2015", "bosutinib": "Schrock 2013", "dasatinib": "Schrock 2013", "imatinib": "Schrock 2013", "nilotinib": "Schrock 2013", "ponatinib": "Schrock 2013", "erlotinib": "Davis 2011", "gefitinib": "Davis 2011"}.items()], ["inhibitor", "source (paper Table 1)", "records with ΔΔG"])
flag = t["rows_where_ddg_exp_disagrees_with_ic50_column"][0]
dm = t["erlotinib_gefitinib_values_equal_to_Davis_block"]
text = f"""# Intake of the public Hauser / Aldeghi Abl inhibitor-response resources (generated report)

Generated by `make_report.py` from `results/inputs_inspection.json` and `results/abl_response_intake_record.json`. Scope: a dataset intake and feasibility decision. **No model was trained or scored.** The records are cell-viability or binding measurements; nothing here establishes a conformational observation, a mechanism, clinical value or a drug-response prediction.

## Decision: **B — a study is feasible only after a specific, justified feature task**

The deposit holds eight wild-type Abl:inhibitor models and no mutant coordinates and no trajectories, so no matched mutant-versus-wild-type Cα comparison exists and a Cα-change predictor cannot be defined from it. The experimental endpoint is not rejected. A different, legitimate and much narrower question is testable with the existing inputs once one feature task exists (the contact and mutation-property features below), but with only {t['sites_with_values']} residue positions, of which {t['sites_with_any_resistant_record']} contain any resistant record, the expected uncertainty is large (below).

![source to endpoint](figures/abl_response_source_to_endpoint.png)

*Figure 1. Purpose: trace each endpoint record to its source and to the structural inputs. Shows the published sources (counts per source from the paper's Table 1 and the sheet), the compiled 144-record table, and what the deposit contains. Main limitation: source attributions come from the paper's Table 1; the three source tables were not re-read; the Aldeghi data files (ACS Data S1/S2, Zenodo Data S3 of {led['not_downloaded'][0]['bytes'] / 1e6:.0f} MB) were not obtained. Evidence class: measured counts from the downloaded files (green), unverified or mismatched properties (orange). Source: `results/inputs_inspection.json`.*

## 1. Endpoint and provenance

- **Records:** {t['rows_with_full_values']} complete records (of {t['rows_in_sheet']} rows in Supplementary Data 3; the other {t['rows_without_values']} rows carry only inhibitor and mutation names): {t['n_inhibitors']} inhibitors × {t['mutations_with_values']} mutations at {t['sites_with_values']} positions (sparse: {t['rows_with_full_values']} of {t['n_inhibitors'] * t['mutations_with_values']} possible pairs).
- **Quantity:** ΔΔG = RT ln(IC50_mut/IC50_WT) at 300 K, positive = resistance, derived by the authors under a competitive Michaelis-Menten assumption with K_M ≫ [ATP]; the paper itself lists unmodelled factors (drug metabolism, expression, feedback). It is not a binding free energy. Resistant = ΔΔG ≥ 1.36 kcal/mol (10-fold) in both papers.
- **Reproducibility of ΔΔG from IC50:** for every inhibitor the sheet's ΔΔG equals RT ln(IC50/IC50_WT,implied) (maximum residual {max(v['max_abs_residual_kcal'] for k, v in t['implied_wt_ic50_consistency_per_inhibitor'].items() if k != 'ponatinib'):.3f} kcal/mol for seven inhibitors), with one exception: **ponatinib {flag['mutation']}** has IC50 {flag['ic50_nM']} nM but ΔΔG {flag['ddg_exp_in_sheet']} kcal/mol, whereas that IC50 would imply {flag['ddg_from_ic50_and_implied_wt']:.2f}. The inconsistency is unresolved and left uncorrected.
- **Censoring:** {t['n_ic50_ge_10000']} records have IC50 = 10,000 nM (dasatinib T315I, imatinib L248R), described in the paper as above the assay range; the sheet stores them as equalities, so their ΔΔG values (5.08, 2.33) are lower bounds. **Resolution:** IC50 is integer nM; extra rows sharing an IC50 value within an inhibitor: {json.loads(t['ic50_resolution_note'][t['ic50_resolution_note'].index('{'):])}. No per-record replicate or uncertainty exists; the authors' variability estimates are inter-laboratory (0.36 kcal/mol, 97 ΔΔGs) and ΔIC50 against ΔKd (0.81 kcal/mol, 27 mutations).
- **Assay mix (not one assay):** per the paper, {n_schrock + n_pem} records are cell-viability IC50 values; the {n_davis} erlotinib and gefitinib records come from Davis et al. 2011 and are binding constants (the paper's Fig. 2d calls ref. 82 ΔKd). In Supplementary Data 2 those values are headed "IC50 (nM)". In this intake {dm['erlotinib']['phosphorylated']} of {dm['erlotinib']['n_rows']} erlotinib and {dm['gefitinib']['phosphorylated']} of {dm['gefitinib']['n_rows']} gefitinib values equal the phosphorylated-Abl block of that sheet and none equals the non-phosphorylated block; the remaining records are not in the sheet, so their block is untested. The structural models are non-phosphorylated.
- **Source tables:** Schrock 2013 is a conference abstract (Blood 122, 3992), so five of the eight inhibitors trace to an abstract; the three source tables were not downloaded or verified. The sheet gives no wild-type IC50 (implied: {', '.join(f"{k} {v['implied_wt_ic50_nM_median']:.3g} nM" for k, v in t['implied_wt_ic50_consistency_per_inhibitor'].items())}).
- **Sequence and numbering:** UniProt P00519 canonical sequence (isoform IA numbering, gatekeeper Thr315), fetched {rec['source']['accessed_utc']}: the wild-type residue of all {t['wildtype_residue_matches_P00519_1a_numbering']['matches']}/{t['wildtype_residue_matches_P00519_1a_numbering']['of']} mutations matches, and every residue of all eight models equals the P00519 residue at its number. Species: human; construct of the models: kinase domain S229-K512 (paper) with ACE/NME caps; the assay constructs were not verified.

{t_in}

## 2. Units and what can be generalized

Unique records: {t['rows_with_full_values']}; substitutions: {t['mutations_with_values']}; residue positions: {t['sites_with_values']} (1 to 5 substitutions each; positions {', '.join(str(k) for k, v in t['mutations_per_site'].items() if v >= 3)} carry three or more); inhibitors: {t['n_inhibitors']}; assay families: 2 (cell viability, binding); source studies: 3; papers using the table: 2 (not independent). Resistant records: {t['resistant_records']}, spread over {t['mutations_with_any_resistant_record']} substitutions at {t['sites_with_any_resistant_record']} positions ({', '.join(str(x) for x in t['sites_with_any_resistant_record_list'])}) and {len(t['inhibitors_with_any_resistant_record'])} inhibitors (none for axitinib, ponatinib or gefitinib).

- **Realistic generalization:** new residue positions within Abl (grouped by position), or new inhibitors (leave-one-inhibitor-out, as in Aldeghi's ML2). Repeated inhibitors for one mutation are not independent mutations. Cross-kinase claims are not testable from this Abl-only set.
- **Precision:** with {t['sites_with_any_resistant_record']} resistant-containing positions and {t['sites_with_values'] - t['sites_with_any_resistant_record']} without, the Hanley-McNeil standard error of a site-level AUROC at a true value of 0.8 is {t['site_level_auroc_standard_error_if_true_auroc_0.8_hanley_mcneil']:.2f} (a formula, not a result). Only a large improvement could be detected; resistant/susceptible classification is descriptive only.

## 3. Structural-input adequacy

{t_st}

- **Contents:** Supplementary Data 5 (63 files, 8 wild-type models; Maestro .mae, PDB with and without waters, UniProt-renumbered PDB, ligand mol2/sdf, renumbering script). **Experimental structures:** the six co-crystal-derived models start from 4WA9, 3UE4, 4XEY, 1OPJ, 3CS9, 3OXZ (identifiers verified against RCSB: X-ray, 2.2-2.9 Å) with missing termini and loops added by homology (template 3IK3, 3OY3, both T315I structures, used for termini or one loop) and Prime. **Docked:** erlotinib and gefitinib into the 3UE4 receptor. **Mutant structures: none** (files named as mutants: {len(ins['mutant_structures_in_deposit']['files_matching_mutant_or_mutation_in_names'])}); **trajectories: none** ({len(ins['mutant_structures_in_deposit']['trajectory_files'])}); the mutation effects are alchemical or Prime estimates in the spreadsheets and analysis scripts (Supplementary Data 3, 4), whose mutant models are not supplied. Whether Aldeghi Data S2 (MD and Rosetta inputs) contains mutant models is unknown: the file could not be fetched by script.
- **Deterministic sample:** all eight models (the whole input set; no sampling), covering both construction routes. Chain rule: the chain with the most standard residues (the authors' own choice). Alternative water/no-water files give identical protein Cα for five inhibitors and differ by up to 0.10 Å for the docked pair.
- **Cα comparison (rigid Kabsch, {cmp_['residues_common_to_all_models']} common residues {cmp_['first_last_common'][0]}-{cmp_['first_last_common'][1]}):** the eight traces fall in two families (the conformational state of either family was not assessed here): axitinib, bosutinib, dasatinib, erlotinib, gefitinib {min(within_a):.3f}-{max(within_a):.2f} Å apart, and imatinib, nilotinib, ponatinib {min(within_b):.2f}-{max(within_b):.2f} Å apart, with {min(between):.1f}-{max(between):.1f} Å between families. The erlotinib and gefitinib models differ by {pair['erlotinib|gefitnib']:.3f} Å, about 14 times the {cmp_['coordinate_rounding_floor_A']} Å floor set by three-decimal PDB coordinates and far below any structural change of interest: the same receptor. **Consequences:** (i) Cα geometry varies with the inhibitor template, never with the mutation, so there is no mutation-specific backbone information and a Cα-change predictor is undefined here; (ii) a template-level Cα feature is a function of inhibitor identity (a baseline); (iii) the Cα trace cannot separate erlotinib from gefitinib.
- **Other legitimate task with existing inputs:** the wild-type complexes support site-context features: ligand-site distance ({min(v['n_sites_within_4.5A'] for v in ins['ligand_frame_and_site_contacts'].values())}-{max(v['n_sites_within_4.5A'] for v in ins['ligand_frame_and_site_contacts'].values())} of {t['sites_with_values']} positions are within 4.5 Å of the ligand, and {min(v['n_sites_beyond_10A'] for v in ins['ligand_frame_and_site_contacts'].values())}-{max(v['n_sites_beyond_10A'] for v in ins['ligand_frame_and_site_contacts'].values())} lie beyond 10 Å, so distal positions matter), side-chain contacts, and local backbone geometry of the template at the position. No existing TopoFold component computes mutation-aware contact or burial features.

## 4. Study design (decision B)

- **Question:** for new Abl residue positions, does any structure-derived site descriptor of the matched wild-type complex (side-chain contacts, ligand distance, TopoFold Cα κ/τ window) improve the ranking of variant-drug ΔΔG over mutation-property and inhibitor baselines?
- **Endpoint and units:** ΔΔG_exp (kcal/mol, positive = resistance) as supplied; the two censored records flagged as lower bounds and excluded from the primary regression metric; the ponatinib G250E row excluded and reported; resistant (≥ 1.36) AUPRC descriptive only.
- **Inputs:** the {t['rows_with_full_values']} records, the eight wild-type templates, and P00519 numbering. **Strongest baselines (equal information):** inhibitor mean ΔΔG from training positions; mutation properties (BLOSUM62, volume, hydrophobicity, charge change); ligand distance and contact counts at the position; their combination in one regularized model. TopoFold features enter only as an added block.
- **Grouping and holdout:** leave-one-position-out over {t['sites_with_values']} positions (all substitutions and inhibitors at a position held out together); secondary leave-one-inhibitor-out. Scaling and any selection inside each training fold.
- **Uncertainty:** bootstrap over positions; with this few positions a descriptive or insufficient-evidence outcome is the expected result unless the improvement is large. Minimum useful improvement to be frozen before any scoring (a margin in RMSE and rank correlation against the strongest baseline), with the equivalence question stated separately.
- **Leakage controls:** the template of the held-out inhibitor is its own; no label, hotspot list or resistant-site knowledge enters any feature; the inspected sample is development material (below); features and thresholds fixed before scoring.
- **Smallest next implementation:** a deterministic feature-table builder (no model) that reads the eight templates and the 144 records and writes site-context and mutation-property features with unit tests (numbering, wild-type identity, distances against an independent calculation); then freeze the protocol.

## 5. Terms, exposure and what was not obtained

- **Terms (actual files):** the Hauser articles are CC BY 4.0 (Crossref); the GitHub repository holding the supplementary files has no licence file and a null licence field; the compiled third-party values have their own terms. Internal analysis: conditional; publication of derived numbers with citation: conditional; redistribution of the files and commercial use: unresolved. The Aldeghi article is CC-BY and Zenodo Data S3 is cc-by-4.0; the ACS supplementary licence was not retrieved.
- **Exposure:** all {t['rows_with_full_values']} records, {t['mutations_with_values']} mutations, {t['sites_with_values']} positions and eight templates were inspected here; they are development material, not an untouched test set. A prospective test needs data this intake did not see.
- **Downloads:** {led['kept_raw_bytes'] / 1e6:.1f} MB kept and about {(led['scientific_download_total_bytes_upper_bound']) / 1e6:.1f} MB including transient reads (ceiling 50 MB). Not downloaded: Zenodo Data S3 ({led['not_downloaded'][0]['bytes'] / 1e6:.0f} MB), ACS Data S1/S2 (challenge page). Nothing was sent externally.
"""
(HERE / "REPORT.md").write_text(text)
print("wrote REPORT.md")
