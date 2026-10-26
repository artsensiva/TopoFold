"""Generates REPORT.md and the figure from results/*.json (no number is typed here).

    .venv/bin/python benchmarks/audit/kin_gen_2/make_report.py
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
R = HERE / "results"
(HERE / "figures").mkdir(exist_ok=True)
dm = json.loads((R / "dev_model.json").read_text())
dm1 = json.loads((R / "dev_model_run1.json").read_text())
ev = json.loads((R / "eval_results.json").read_text())
src = json.loads((R / "reference_source_check.json").read_text())
INK, GREY, BLUE, ORANGE, GREEN, LIGHT = "#2b2a27", "#6b6a66", "#2a78d6", "#eb6834", "#2f8f6b", "#e4e3df"
COL = {"in": BLUE, "out": ORANGE, "inter": GREEN}
REPS = ("T", "B1", "B2", "B3")
NAME = {"T": "T: TopoFold κ/τ", "B1": "B1: window Cα distances", "B2": "B2: local RMSD", "B3": "B3: all 36 Cα distances"}
e_out = ev["eligibility"]
plt.rcParams.update({"font.size": 8.5, "text.color": INK, "axes.edgecolor": GREY})
fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 4.0), gridspec_kw={"width_ratios": [1.15, 1]})
for grp in ("in", "inter", "out"):
    d = [x for x in dm["entries"] if x["group"] == grp]
    a.scatter([x["D1"] for x in d], [x["D2"] for x in d], s=14, color=COL[grp], alpha=0.35, linewidths=0, label=f"development {grp} (n={len(d)})")
    e = [x for x in ev["entries"] if x["group"] == grp]
    a.scatter([x["D1"] for x in e], [x["D2"] for x in e], s=26, color=COL[grp], marker="D", edgecolors=INK, linewidths=0.5, label=f"evaluation {grp} (n={len(e)})")
for grp, c in dm["groups"]["centroids"].items():
    a.scatter([c[0]], [c[1]], marker="+", s=90, color=INK)
    a.annotate(grp, c, textcoords="offset points", xytext=(5, 5), fontsize=8)
a.set_xlabel("D1 (Å): αC Glu+4 Cα to DFG-Phe Cζ")
a.set_ylabel("D2 (Å): β3 Lys Cα to DFG-Phe Cζ")
a.legend(fontsize=6.5, frameon=False, loc="upper right")
a.set_title("(a) the independent D1D2 reference groups", fontsize=9, loc="left")
vals = [ev["scores"][r]["weighted"]["ba_w"] for r in REPS]
bars = b.bar(range(4), vals, color=[ORANGE if r == "T" else GREY for r in REPS], width=0.62)
for i, v in enumerate(vals):
    b.text(i, v + 0.006, f"{v:.3f}", ha="center", fontsize=7.5)
b.axhline(1.0, color=LIGHT, lw=1)
b.set_ylim(0.8, 1.05)
b.set_xticks(range(4))
b.set_xticklabels(["T", "B1", "B2", "B3"])
b.set_ylabel("weighted balanced accuracy (descriptive)")
n_cells = f"{e_out['n_out']} out entries from {e_out['kinases_with_out']} kinases (rule: at least 10 from at least 3)"
b.set_title(f"(b) evaluation, descriptive only\nINSUFFICIENT EVIDENCE: {e_out['n_out']} out entries,\n{e_out['kinases_with_out']} kinases (rule: at least 10, at least 3)", fontsize=8, loc="left")
for ax in (a, b):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(HERE / "figures/kin_gen_2_summary.png", dpi=130)
plt.close(fig)


def table(rows, head):
    return "\n".join(["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(str(x) for x in r) + " |" for r in rows])


dev_rows = [[k, v["accession"], f"{dm['per_kinase_group_counts'][k]['in']} / {dm['per_kinase_group_counts'][k]['out']} / {dm['per_kinase_group_counts'][k]['inter']}",
             "; ".join(v["problems"])[:90] or ("αC Glu " + str(v.get("alphaC_glu", {}).get("selected")))] for k, v in dm["kinases"].items()]
t_dev = table(dev_rows, ["development kinase", "UniProt", "entries in / out / inter", "anchor note"])
ev_rows = []
for r in ev["processed_kinases"]:
    c = ev["counts_per_kinase"][r]
    f = ev["family_relation"].get(r)
    ev_rows.append([r, ev["kinases"][r]["accession"], ev["kinases"][r]["search_total"], f"{c['in']} / {c['out']} / {c['inter']}",
                    (f"{f['max_identity']:.2f} ({f['nearest_development_kinase']})" + (" related" if f["family_flag"] else "")) if f else "n/a",
                    "; ".join(ev["kinases"][r]["problems"])[:80] or ("αC Glu " + str(ev["kinases"][r].get("alphaC_glu", {}).get("selected")))])
t_ev = table(ev_rows, ["evaluation kinase (rank order)", "UniProt", "RCSB eligible entries", "entries in / out / inter", "kinase-domain identity to the nearest development kinase", "anchor note"])
s_rows = []
for r in REPS:
    s = ev["scores"][r]
    po = {k: f"{v['tpr']:.2f}/{v['tnr']:.2f}" for k, v in s["per_kinase"].items() if v["n_out"] > 0}
    s_rows.append([NAME[r], s["k"], f"{s['weighted']['ba_w']:.3f}", f"{s['pooled_ba']:.3f}", f"{s['auroc']:.3f}", ", ".join(f"{k} {v}" for k, v in po.items()), f"{dm['loko'][r][str(s['k'])]['ba']:.3f}"])
t_s = table(s_rows, ["representation", "k (development)", "weighted BA", "pooled BA", "AUROC of the neighbour score", "TPR/TNR in kinases with out entries", "development leave-one-kinase-out BA (development, not evidence)"])
g = dm["gates"]
t_g = table([[k, "pass" if v["pass"] else "FAIL"] for k, v in g.items()], ["development gate", "result"])
import collections  # noqa: E402
exc = collections.Counter()
for e in ev["exclusions"]:
    reason = e.get("reason", "")
    exc["no PDB-format file" if "PDB-format" in reason else "kinase_roles refused (role mutation or placement)" if "kinase_roles refused" in reason else "window chain break" if "chain break" in reason else "a reference role unobserved" if "unobserved" in reason else "no chain maps" if "no chain maps" in reason else "non-standard residue in the window" if "non-standard" in reason else reason[:50]] += 1
t_x = table(sorted(exc.items()), ["evaluation entries or kinases excluded, by reason", "count"])
res = ev["resources"]
dev_cost_run1 = dm1["preparation_cost"]
text = f"""# KIN-GEN-2: Cα-only TopoFold features against equal-information baselines for DFG-in versus DFG-out groups (generated report)

Generated by `make_report.py` from `results/*.json`. **Evidence scope: structural-label recovery only.** The labels are geometric groups of static public structures; nothing here establishes function, drug response or biological validation.
Protocol `PROTOCOL.md` (frozen; amendments 1 and 2 are dated and were made before any evaluation structure was selected). It replaces the unfrozen KIN-GEN-1 draft.

## Outcome

**Outcome (frozen rule): {ev['outcome']['label'].upper()}.** The evaluation yielded {e_out['n_scored']} scored entries but only {e_out['n_out']} labelled out, from {e_out['kinases_with_out']} kinases ({', '.join(k + ' ' + str(v['out']) for k, v in ev['counts_per_kinase'].items() if v['out'])}); the rule requires at least 10 from at least 3 kinases.
The benchmark therefore supports none of: comparable recovery, an advantage, or a Cα-only limitation. The scores below are descriptive: all four Cα-only representations separate in from out with weighted balanced accuracy {min(ev['scores'][r]['weighted']['ba_w'] for r in REPS):.3f} to {max(ev['scores'][r]['weighted']['ba_w'] for r in REPS):.3f}, and the differences correspond to one or two entries.

![summary](figures/kin_gen_2_summary.png)

*Figure 1. Purpose: show the reference groups used as labels and the descriptive recovery scores. (a) D1 and D2 of every development (small) and evaluation (diamond) entry, coloured by the D1D2 group (fitted on development; evaluation entries assigned to the nearest development centroid, + marks the centroids). (b) weighted balanced accuracy of the four representations on the evaluation entries. Limitation: the labels come from a full-atom rule that the Cα-only representations never see; the rule defines the labels, so it is not a performance result; the evaluation has too few out entries for any comparison; no interval or p-value is given (too few independent kinases). Source: `results/dev_model.json`, `results/eval_results.json`.*

## Reference route (selected before scoring): independent labels, not KinCoRe

Official pages re-read: the KinCoRe website data ("Last updated {next(iter(src['pages'].values()))['last_updated_text']}") state no licence or terms; the GitHub standalone repositories have no licence file; the article is {src['paper']['licence_statement'].split('under ')[-1].strip().rstrip(' .')}, and that licence is not established for the website data or the supplementary datasets; the paper publishes the distances D1 and D2 and the clustering method but no numeric boundaries.
Unresolved and not assumed: {'; '.join(src['assessment']['unresolved'])}. **Route B** was therefore used: D1 and D2 computed here from public PDB coordinates by an independently written implementation, three groups from average-linkage clustering fitted on development data only. The groups are "D1D2 groups", not KinCoRe labels and not biological measurements. No official label, table or code was read or copied.

## Development (preparation; recorded separately)

{t_dev}

Development entries: {len(dm['entries'])} ({dm['groups']['sizes']['in']} in, {dm['groups']['sizes']['out']} out, {dm['groups']['sizes']['inter']} inter); group centroids (D1, D2 in Å): {'; '.join(f"{k} ({v[0]:.1f}, {v[1]:.1f})" for k, v in dm['groups']['centroids'].items())}. k per representation (leave-one-kinase-out): {dm['k_selected']}; primary baseline fixed on development: **{dm['primary_baseline']}**.
CDK2 yielded no development entry (no unique αC Glu: almost all its entries are αC-out).

{t_g}

Gates: G1 core κ/τ against numpy; G2 rigid-motion invariance; G3 Kabsch against MDAnalysis QCP; G4 role residues and shifted-mapping control; G5 controls; G6 group sizes; G8 D1/D2 against MDAnalysis; G7 (disjointness) is checked in the evaluation: {'pass' if ev['G7']['pass'] else 'FAIL'}. Two pre-evaluation amendments are documented in the protocol (float32 limit for T; anchor rules after MAPK14, AURKA and CDK2 failed; control testability).

## Evaluation selection and counts

{t_ev}

{t_x}

## Scores (descriptive; primary baseline {dm['primary_baseline']})

{t_s}

*Table 1. Purpose: the same k-NN procedure for every representation (k from development, decision at half of the neighbours). The primary metric is the weighted balanced accuracy (each kinase × class cell weighted equally within its class). Limitation: {e_out['n_out']} out entries, so one entry changes a rate by {1 / e_out['n_out']:.2f} within its cell; the in class dominates ({e_out['n_in']} entries). The score is not tuned on evaluation data; no window or feature was selected.*

Family relation: kinase-domain identity of each evaluation kinase to its nearest development kinase is in the evaluation table; {sum(1 for v in ev['family_relation'].values() if v['family_flag'])} of {len(ev['family_relation'])} kinases are flagged related (identity at least 0.40); the out entries come from {', '.join(k for k, v in ev['counts_per_kinase'].items() if v['out'])}, {'none of them flagged' if not any(ev['family_relation'][k]['family_flag'] for k, v in ev['counts_per_kinase'].items() if v['out']) else 'some flagged'}.

## What this does and does not show

- **Not established (insufficient evidence):** whether the Cα-only TopoFold representation recovers the D1D2 group better than, as well as, or worse than the baselines; whether Cα-only information limits recovery. The frozen evaluation did not reach 10 out entries from 3 kinases. Hypothesis (not tested here): entries were taken by resolution without stratifying by class, and out-type structures are a small fraction of the kinases ranked first.
- **Descriptive observation:** on this evaluation set all four representations recovered most or all out entries and almost all in entries (Table 1), so a Cα window of nine residues around the DFG motif carries enough information to recover this label in the entries examined, for in-dominated data. This is a descriptive observation, not a comparison.
- **Development note (not evidence):** on development, T had the highest leave-one-kinase-out balanced accuracy ({dm['loko']['T'][str(dm['k_selected']['T'])]['ba']:.3f} against {dm['loko'][dm['primary_baseline']][str(dm['k_selected'][dm['primary_baseline']])]['ba']:.3f} for {dm['primary_baseline']}); development kinases are few and several are related, so this is not a performance estimate.
- The full-atom D1/D2 rule defines the labels; its agreement with them is 1 by construction.

## Resources

Evaluation: {res['evaluation_download_bytes'] / 1e6:.1f} MB of structure files downloaded ({res['files']} files; limit 100 MB), {res['cpu_seconds']:.0f} CPU-s ({res['wall_seconds']:.0f} s wall; limit 2 CPU-hours).
Preparation (recorded separately): development run 1 downloaded {dev_cost_run1['download_bytes'] / 1e6:.1f} MB ({dev_cost_run1['files']} files), runs 2 and 3: see `results/dev_prep_run2.log` and `results/dev_prep_run3.log`; run 1 fitting and gates took {dev_cost_run1['cpu_seconds']:.0f} CPU-s.
"""
(HERE / "REPORT.md").write_text(text)
print("wrote REPORT.md and figure")
