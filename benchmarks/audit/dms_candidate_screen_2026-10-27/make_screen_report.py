"""Writes results/dms_candidate_screen.json (compact evidence table), SCREEN.md, FEATURE_EXPORT_SPEC_MET_DRAFT.md and
figures/dms_screen_decision.png from results/met_screen.json, results/mavedb_00000091-a-1_record.json and the recorded notices.
Nothing here trains or scores a model.   python3 -I make_screen_report.py  (needs matplotlib: use the repository .venv)
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
met = json.loads((HERE / "results/met_screen.json").read_text())
bra = json.loads((HERE / "results/mavedb_00000091-a-1_record.json").read_text())
pc = met["per_condition"]
DRUGS = [k for k in pc if k != "DMSO"]
dm = pc["DMSO"]
tol = [pc[k]["paired_and_not_DMSO_LOF"] for k in DRUGS]
lof_frac = dm["flag_LOF"] / dm["missense_variants_with_effect"]
gof = [pc[k]["flag_GOF"] for k in DRUGS]
sd = met["posterior_sd_of_filtered_scores"]
led = met["download_ledger"]
st = met["structures"]

hauser_closure = {
    "statements": [
        {"source": "Article page, Rights and permissions (https://www.nature.com/articles/s42003-018-0075-x)", "quote": "Open Access This article is licensed under a Creative Commons Attribution 4.0 International License ... The images or other third party material in this article are included in the article's Creative Commons license, unless indicated otherwise in a credit line to the material. If material is not included in the article's Creative Commons license and your intended use is not permitted by statutory regulation or exceeds the permitted use, you will need to obtain permission directly from the copyright holder.", "scope": "the article and its included third-party material"},
        {"source": "Journal policy page (https://www.nature.com/commsbio/editorial-policies/self-archiving-and-license-to-publish)", "quote": "Open access articles in Communications Biology are published under a CC BY-NC-ND ... or CC BY ... licence.", "scope": "licence options for articles; authors are not required to assign copyright"},
        {"source": "Same article page: Supplementary Information link (media.springernature.com .../42003_2018_75_MOESM1_ESM.pdf)", "quote": "(no licence statement attached to the link)", "scope": "the page lists the official supplement without a separate notice"},
        {"source": "Springer Nature licensing and copyright page and the reporting-standards page", "quote": "(searched: no statement on the licence of supplementary files)", "scope": "none found"},
        {"source": "RCSB policies (https://www.rcsb.org/pages/policies)", "quote": "data files contained in the PDB archive are available under the CC0 1.0 Universal (CC0 1.0) Public Domain Dedication", "scope": "the starting PDB entries 4WA9, 3UE4, 4XEY, 1OPJ, 3CS9, 3OXZ only, not the prepared models"},
        {"source": "GitHub repository kehauser/Predicting-resistance-... (commit d5c6c508)", "quote": "licence field null; no licence file", "scope": "the GitHub copies of the supplementary data and the prepared structural inputs (Supplementary Data 5)"}],
    "recorded_status": {"article_text": "CC BY 4.0 (stated)", "official_supplement_PDF": "no separate statement found; unresolved", "GitHub_copies_and_prepared_models": "no licence stated; unresolved", "starting_PDB_entries": "CC0 (wwPDB policy)", "third_party_source_values": "unresolved"},
    "note": "No legal conclusion is drawn. Clearer terms alone would not make the Hauser table a strong validation dataset: the scientific limitations of D90 stand."}

rows = [
    {"candidate": "A. MET kinase domain inhibitor DMS (Estevam et al., eLife 2025; Ba/F3 TPR-MET)", "access_terms": "Article: Creative Commons Attribution (Crossref: CC BY 4.0). Repository fraser-lab/MET_kinase_Inhibitor_DMS (commit c2052b97): LICENSE = MIT, copyright Fraser Lab 2024; article data statement: 'Original data files, analysis, and source code is available at' that repository; raw reads at SRA PRJNA1136906 (not downloaded). The MIT text refers to 'Software and associated documentation files'; whether it is the stated licence of the data tables is not spelled out.",
     "endpoint": f"Rosace effect (log2 fitness scale; wild-type synonymous = 0) per variant and condition; 12 conditions: DMSO + 11 inhibitors, each dosed at its experimentally determined IC50 (not a common exposure scale); IL-3 withdrawal selection. Posterior sd of the filtered scores: median {sd['median_sd_inhibitor_missense']:.2f}, max {sd['max_sd_inhibitor_missense']:.2f} (replicate values not in the files read).",
     "counts": f"measured: {met['mapping']['dms_positions']} positions (UniProt P08581 {met['mapping']['uniprot_residue_range'][0]}-{met['mapping']['uniprot_residue_range'][1]}, wild-type identity verified at all positions), 11 inhibitors, {dm['missense_variants_with_effect']} DMSO missense variants of {met['possible_missense_per_condition']} possible, {min(pc[k]['missense_variants_with_effect'] for k in DRUGS)}-{max(pc[k]['missense_variants_with_effect'] for k in DRUGS)} per inhibitor; {met['missense_variants_measured_in_all_12_conditions']} variants in all 12 conditions; 1 source group (one laboratory, one assay system)",
     "function_control": f"DMSO condition per variant: {dm['flag_LOF']} of {dm['missense_variants_with_effect']} DMSO missense variants ({100 * lof_frac:.0f}%) are flagged loss-of-function, so the drug effect is interpretable only for the {min(tol)}-{max(tol)} DMSO-tolerated variants per inhibitor; expression and folding are not measured separately",
     "structures": f"8 PDB entries downloaded (X-ray, {min(float(v['resolution_A']) for v in st.values()):.2f}-{max(float(v['resolution_A']) for v in st.values()):.2f} Å): experimental ligand complexes for 7 inhibitors (2WGJ, 6SDE, 4EEV, 5T3Q, 3RHK, 4R1V, 3QTI) plus ATP-bound 3DKC; 4 inhibitors (capmatinib, glumetinib, cabozantinib, glesatinib analog) have no experimental structure (docked per the article; poses not in the files read); resolved scanned-span residues per entry {min(max(c['resolved_CA_in_scanned_span'] for c in v['chains'].values()) for v in st.values())}-{max(max(c['resolved_CA_in_scanned_span'] for c in v['chains'].values()) for v in st.values())} of 287; engineered differences from P08581: 3QTI Y1194F, Y1234F; 3DKC Y1194F, Y1234F, Y1235D",
     "evaluation": "unit = variant-drug; group = position (287); dependent across the ~19 substitutions of a position and across drugs (shared libraries and DMSO control); drug confounded with template; single assay and laboratory, so no transfer claim beyond it",
     "exposure": "published XGBoost models in the repository (14 features incl. ESM-1b log-likelihood ratio, stability, distance, conformation; test = two held-out inhibitors, 20% of positions, two amino acids); ESM-1b correlation reported 0.50 (DMSO) against 0.28 (inhibitor-treated). MET structures were used as evaluation entries in KIN-GEN-2 (a different task).",
     "verdict": "qualifies for an exploratory comparison"},
    {"candidate": f"B. MaveDB {bra['urn']} (BRAF vemurafenib resistance)", "access_terms": "MaveDB record: licence field null; no primary publication or DOI in the record",
     "endpoint": "log enrichment, median over synonymous codons, two replicates; 'only a small subset of variants were reported'", "counts": f"measured: {bra['numVariants']} variants, 1 drug, target sequence {len(bra['targetGenes'][0]['targetSequence']['sequence'])} residues",
     "function_control": "none stated", "structures": "not assessed (stopped at access and provenance)", "evaluation": "one drug: template-drug confounding cannot be separated; too few units for a position-grouped comparison", "exposure": "unknown (no publication linked)", "verdict": "does not qualify: one drug, partial reporting, no licence field, no linked source"},
    {"candidate": "C. Persky et al., Nat. Struct. Mol. Biol. 2020 (inhibitor-resistance DMS of several kinases)", "access_terms": "subscription article (Crossref licence element: Springer TDM only; Europe PMC: not open access); no data repository or supplementary table reached within the bounds",
     "endpoint": "not read", "counts": "not read", "function_control": "not read", "structures": "not read", "evaluation": "not assessed", "exposure": "unknown", "verdict": "not evaluated: essential access missing within the bounded search"}]
out = {"schema": "topofold.dms-candidate-screen/1", "date": "2026-10-27", "hauser_permissions_closure": hauser_closure, "candidates": rows, "decision": {"value": "proceed (exploratory)", "dataset": "A. MET kinase domain inhibitor DMS", "conditions": ["use is based on the stated CC BY article licence and the repository's MIT licence file; the owner confirms this basis before any real-data run", "inhibitor-wise analysis only (doses differ per drug); analysis restricted to DMSO-tolerated variants", "primary inhibitors limited to those with experimental ligand-bound structures", "published models and feature files of the authors are excluded as inputs"]},
       "download_ledger": led, "ceiling_bytes": 20_000_000}
(HERE / "results/dms_candidate_screen.json").write_text(json.dumps(out, indent=1))


def table(rs):
    head = ["", "A. MET DMS", "B. MaveDB BRAF/vemurafenib", "C. Persky 2020"]
    keys = [("access_terms", "access and terms"), ("endpoint", "endpoint, assay, uncertainty"), ("counts", "measured counts"), ("function_control", "function control"), ("structures", "structural inputs, mapping"), ("evaluation", "evaluation, dependence, confounding"), ("exposure", "exposure, published baselines"), ("verdict", "verdict")]
    return "\n".join(["| " + " | ".join(head) + " |", "|---|---|---|---|"] + ["| " + lab + " | " + " | ".join(str(r[k]).replace("|", "/") for r in rs) + " |" for k, lab in keys])


md = f"""# Candidate screen for a variant-drug response test of structural site context (generated, 2026-10-27)

Generated by `make_screen_report.py` from `results/*.json`. Task: one bounded screen of at most three alternatives plus the narrow Hauser permissions check. **No model was trained or scored; nobody was contacted.** Downloads: {led['total_bytes_upper_bound'] / 1e6:.1f} MB of 20 MB (structure files, two DMS tables, article XML, publisher and policy pages). Different tasks are not substituted for each other: all three candidates have a drug-exposure endpoint; functional scans without drug exposure were not considered.

## Decision: **Proceed (exploratory) with candidate A, the MET kinase inhibitor DMS**

![decision](figures/dms_screen_decision.png)

*Figure 1. Purpose: show how the three candidates fared against the gates and where the decision falls. Source: `results/dms_candidate_screen.json`. Limitation: the gates are a screen of access, endpoint, controls, inputs and dependence; none of them measures predictive value, and candidate C was stopped at access.*

## Evidence table

{table(rows)}

*Table 1. Counts are measured from the downloaded files (A, B record) except where stated; statements in the exposure row come from the article text and the repository file names (the authors' model files were not read).*

## What the proposed evaluation can and cannot resolve

- **Unit and group:** variant-drug records of the six ATP-competitive inhibitors with experimental complexes, held out by position ({met['mapping']['dms_positions']} positions). Records at one position and across drugs are dependent, so position is the independent unit and drugs are correlated replicates, not independent tests.
- **Can resolve:** whether adding a site-context block (ligand-site distance, contacts, TopoFold or distance-based Cα windows of the template) to drug and mutation-property baselines gives a consistent direction of change in out-of-fold agreement across held-out positions, with position-level resampling as a descriptive interval. Per-variant noise is large (posterior sd median {sd['median_sd_inhibitor_missense']:.2f} on a scale where the authors' resistance cut is 0.5), so only an effect that is large relative to that noise and consistent across positions is interpretable.
- **Cannot resolve:** small gains; resistant/susceptible classification (the GOF-flagged variants number {min(gof)}-{max(gof)} per inhibitor; their position distribution was not examined; descriptive only); transfer to other assays, laboratories or kinases; whether any gain exceeds what the authors' own structural features already give; mutation-specific backbone effects (no mutant coordinates).
- **Confounds recorded:** drug ⇔ template; docked inhibitors without poses; engineered template residues (3QTI, 3DKC); one source group; expression not separated from function.

## Hauser permissions, closed narrowly

{chr(10).join(f"- **{s['source']}:** {s['quote']} *(scope: {s['scope']})*" for s in hauser_closure['statements'])}

Recorded statuses: {json.dumps(hauser_closure['recorded_status'])}. {hauser_closure['note']}

## Decision conditions and next step

Conditions: {'; '.join(out['decision']['conditions'])}. The feature-export specification is `FEATURE_EXPORT_SPEC_MET_DRAFT.md` (draft, unfrozen).
"""
(HERE / "SCREEN.md").write_text(md)

spec = f"""# Draft feature-export specification for the MET inhibitor DMS (unfrozen; not implemented)

Purpose: the smallest exporter whose output can be used for the exploratory comparison in `SCREEN.md`. No model, no feature selection, no outcome reading.

## Rows and files
- **Unit:** (position p in 1..{met['mapping']['dms_positions']}, substituted residue, inhibitor). Inhibitors used for the structure-dependent blocks: crizotinib (2WGJ), savolitinib (6SDE), merestinib (4EEV), AMG-458 (5T3Q), tepotinib (4R1V), NVP-BVU972 (3QTI). Tivantinib (MET-independent activity per the article) and the four docked inhibitors (no poses available) are not exported with structure blocks; this rule uses article statements only.
- **Outputs:** `features.csv` (inputs only) with a `features_manifest.json` (units, definitions, input SHA-256, template construction, missingness counts). Outcome tables are produced by a separate script and are never read by the exporter.

## Mapping rules
- DMS position p = UniProt P08581 residue p + {1058} (wild-type identity verified at all {met['mapping']['dms_positions']} positions; the exporter asserts it from the file hash).
- Template residues are matched by PDB residue number (equal to UniProt numbering in the entries checked); a template residue whose name differs from P08581 is flagged (3QTI 1194 and 1234, 3DKC 1194, 1234, 1235) and its structural features set missing.
- Chain: the chain with the most resolved residues in 1059-1345, ties to the first chain; altloc blank or A; hydrogens and waters ignored; unresolved residue -> missing value plus a missingness flag; no imputation in the exporter.

## Feature blocks (fixed in advance, units explicit)
- **M, mutation properties:** wild-type and substituted residue; volume change (Å^3), hydropathy change (Kyte-Doolittle), charge change, BLOSUM62 score; proline and glycine involvement flags.
- **L, template site context (per inhibitor):** minimum heavy-atom distance between the residue and the ligand (Å); count of residue atoms within 4.5 Å of the ligand; neighbour count within 10 Å of the Cα.
- **T, TopoFold local geometry:** curvature, torsion (sin and cos) of the 9-residue Cα window centred on p from the template (existing `compute_ca_features`); window with any unresolved residue -> missing.
- **C, equal-information Cα baseline:** the 36 pairwise Cα distances of the same window (Å), same residues and frames as T.
- **D, drug baselines:** inhibitor identity and class; no outcome statistic is exported (training means are computed inside folds by the analysis, not here).

## Tests (written before any real-data run)
1. Numbering: offset {1058} reproduces wild-type identity on the real file hash; a shifted offset fails.
2. Template identity: the flagged residues above are exactly those found; no other mismatch.
3. Distances: ligand-site and window distances equal an independent computation (scipy cdist) on synthetic coordinates, and one real template spot-check.
4. Missing data: a synthetic window with a gap yields missing values and a flag, never an imputed value.
5. Outcome independence: output bytes are identical with the outcome file absent, shuffled or altered; the exporter module imports no outcome reader.
6. Invariance: T is unchanged under rigid motion; a mirror image changes torsion sign only (existing core tests reused).
7. Determinism: two runs give identical hashes; no column named like an outcome or a published prediction.

Stopping condition: stop and report if numbering or template-identity checks fail beyond the listed engineered residues, or if the owner does not confirm the licence basis.
"""
(HERE / "FEATURE_EXPORT_SPEC_MET_DRAFT.md").write_text(spec)

# ---- figure ----
INK, GREY, BLUE, ORANGE, GREEN, LIGHT = "#2b2a27", "#6b6a66", "#2a78d6", "#eb6834", "#2f8f6b", "#f1f0ec"
fig, ax = plt.subplots(figsize=(10.2, 4.2))
ax.axis("off")
ax.set_xlim(0, 10.2)
ax.set_ylim(0, 4.2)


def box(x, y, w, h, txt, fc=LIGHT, ec=GREY, fs=7.4, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=fc, ec=ec, lw=1))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, color=INK, fontweight="bold" if bold else "normal")


def arrow(x0, y0, x1, y1, c=GREY):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=c, lw=1.1))


ax.text(0.1, 3.95, "Candidate", fontsize=8.5, fontweight="bold")
ax.text(3.6, 3.95, "Gate outcome", fontsize=8.5, fontweight="bold")
ax.text(7.5, 3.95, "Decision", fontsize=8.5, fontweight="bold")
box(0.1, 2.75, 3.1, 0.9, f"A. MET inhibitor DMS\n{met['mapping']['dms_positions']} positions, 11 inhibitors, 1 source group", fc="#e3f3ec", ec=GREEN)
box(0.1, 1.55, 3.1, 0.9, f"B. MaveDB BRAF / vemurafenib\n{bra['numVariants']} variants, 1 drug", fc="#fbeadf", ec=ORANGE)
box(0.1, 0.35, 3.1, 0.9, "C. Persky 2020 (subscription)\naccess not reached", fc="#fbeadf", ec=ORANGE)
box(3.6, 2.75, 3.4, 0.9, f"terms stated (CC BY + MIT file); DMSO control;\n{min(tol)}-{max(tol)} tolerated variants per inhibitor;\n7 experimental complexes; mapping verified", fc="#fff", ec=GREEN, fs=7.2)
box(3.6, 1.55, 3.4, 0.9, "no licence field, no linked paper,\npartial reporting, one drug", fc="#fff", ec=ORANGE, fs=7.2)
box(3.6, 0.35, 3.4, 0.9, "no data or terms read within bounds", fc="#fff", ec=ORANGE, fs=7.2)
box(7.5, 2.75, 2.6, 0.9, "PROCEED (exploratory)\nconditions in SCREEN.md", fc="#e3f3ec", ec=GREEN, bold=True)
box(7.5, 1.0, 2.6, 1.45, "not selected / not evaluated;\nHauser table stays a pilot\n(terms still unresolved for the\nsupplement and GitHub copies)", fc=LIGHT)
for y in (3.2, 2.0, 0.8):
    arrow(3.2, y, 3.6, y)
arrow(7.0, 3.2, 7.5, 3.2, GREEN)
arrow(7.0, 2.0, 7.5, 1.9, GREY)
arrow(7.0, 0.8, 7.5, 1.3, GREY)
fig.tight_layout()
fig.savefig(HERE / "figures/dms_screen_decision.png", dpi=130)
plt.close(fig)
print("wrote screen outputs")
