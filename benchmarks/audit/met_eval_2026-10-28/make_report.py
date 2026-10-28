"""REPORT.md and two figures from results/results.json only (no number is typed here). Run with the repository .venv."""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
R = json.loads((HERE / "results/results.json").read_text())
plan = json.loads((HERE / "plan.json").read_text())
INK, GREY, BLUE, ORANGE = "#2b2a27", "#6b6a66", "#2a78d6", "#eb6834"


def table(rows, head):
    return "\n".join(["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(str(x) for x in r) + " |" for r in rows])


res = R["resources"]
md = [f"# MET exploratory evaluation: result (generated from `results/results.json`)\n",
      f"**Status: {R['status']}.** Internal, exploratory, descriptive. Target population: variants passing the measured no-drug eligibility rule; nothing here shows prediction for variants whose baseline function is unknown. Plan v2 and evaluator were frozen (commit 84215b8) before scoring; one scoring attempt. No commercial or redistribution rights are established; the source assay tables are not published.\n"]
cells = R["cells"]
bad = R["non_ok_cells"]
md.append(f"Cells: {len(cells)} inhibitor × fold; non-ok cells: {len(bad)}" + (": " + json.dumps(bad) if bad else "") + ".\n")
rc = R["row_counts"]
md.append("Rows used per inhibitor (eligible → used): " + "; ".join(f"{d} {c['eligible']} → {c['used']} (structure-masked {c['no_structure']}, mutation letter X {c['unknown_mutation']}, no ESM {c['no_esm']})" for d, c in rc.items()) + ".\n")
if "evaluation" in R:
    ev = R["evaluation"]
    comps = ev["comparisons"]
    sm = ev["per_set_model"]
    md.append("## Predictive performance (position-weighted mean squared error of gamma; lower is better)\n")
    md.append(table([[k.split("|")[0], k.split("|")[1], f"{v['overall_mse']:.4f}"] for k, v in sm.items()], ["feature set", "model", "MSE"]) + "\n")
    md.append("![forest](figures/met_eval_comparisons.png)\n")
    md.append("*Figure 1. Purpose: show every pre-declared paired comparison, not a selected one. Points: relative MSE reduction of the richer set over its reference (positive favours the richer set); bars: descriptive 95% moving-block-bootstrap interval (conditional on this dataset and procedure; coverage not assumed); dashed line: the provisional 5% engineering threshold; filled markers: the consistency criterion was met. Limitation: position blocks in sequence, one assay, six correlated inhibitors, wild-type-template features. Source: `results/results.json`.*\n")
    rows = []
    for cid, c in comps.items():
        for m, f in c["families"].items():
            rows.append([cid, c["test"], m, f"{f['delta']:+.4f}", f"{100 * f['relative_mse_reduction']:+.2f}%", f"[{f['interval_95_descriptive'][0]:+.4f}, {f['interval_95_descriptive'][1]:+.4f}]", f"{f['blocks_improving']}/10", f"{f['drugs_improving']}/6", f["classification"]])
    md.append("## Paired comparisons (identical rows and folds)\n")
    md.append(table(rows, ["id", "test (a vs b)", "model", "Delta MSE", "relative reduction", "95% interval (descriptive)", "blocks improving", "inhibitors improving", "classification"]) + "\n")
    md.append("Family agreement: " + "; ".join(f"{cid}: {c['outcome']}" for cid, c in comps.items()) + ".\n")
    md.append(f"**Engineering outcome under the pre-declared rule (provisional 5% threshold; not demonstrated practical or commercial usefulness): {ev['engineering_outcome']}.**\n")
    md.append("![blocks](figures/met_eval_block_consistency.png)\n")
    md.append("*Figure 2. Purpose: show consistency across the ten outer position blocks. Each cell is Delta MSE for one comparison and model family in one block (blue: the richer set is better, orange: worse; colour scale clipped symmetrically). Limitation: blocks are contiguous stretches of sequence with different numbers of eligible variants; a block-level sign is a coarse indicator. Source: `results/results.json`.*\n")
    s = R["sensitivity_low_sd_rows"]
    md.append("## Sensitivity and context\n")
    md.append(f"Re-scoring the same out-of-fold predictions on rows whose approximate sd_gamma is at or below the median ({s['median_sd_gamma_approx']:.3f}; rows kept {json.dumps(s['rows_kept'])}): " + "; ".join(f"{cid} {m}: {s['comparisons'][cid]['families'][m]['classification']} (relative {100 * s['comparisons'][cid]['families'][m]['relative_mse_reduction']:+.1f}%)" for cid in s["comparisons"] for m in ("ridge", "gbm")) + f". Engineering outcome on this subset: {s['engineering_outcome']}.\n")
    md.append(f"Noise context: mean approximate sd_gamma² = {R['noise_context']['mean_approx_sd_gamma_squared']:.3f}. {R['noise_context']['caveat']}.\n")
    cor = R["secondary_correlations_per_drug_pooled_rows"]
    md.append("Secondary per-inhibitor correlations of out-of-fold predictions with gamma (pooled rows; Spearman; these are not comparable with the authors' correlations, which concern a different endpoint and split):\n")
    md.append(table([[k.replace("|", " ")] + [f"{v['spearman']:.2f}" for v in c.values()] for k, c in cor.items()], ["set model"] + list(next(iter(cor.values())).keys())) + "\n")
md.append("## Resources\n")
md.append(f"CPU {res['cpu_seconds_total_process_and_children']:.0f} s for the whole process and children (cap {res['cpu_cap_seconds']} s), wall {res['wall_seconds']:.0f} s, peak resident memory {res['max_rss_mb']:.0f} MB (cap {res['rss_cap_mb']} MB), threads: {res['threads']}, child processes started: {res['child_processes_started']}.\n")
md.append("## Separation of claims\n")
md.append("- **Predictive performance:** the MSE and correlation tables above, for this dataset, these features and this blocked procedure only.\n- **Uncertainty:** the intervals and block/inhibitor counts are descriptive and conditional; shared controls and libraries make inhibitors dependent; sequence blocking does not remove spatial dependence.\n- **Engineering continuation:** the pre-declared rule with the provisional 5% threshold, a prioritisation device only.\n- **Commercial implications:** none are established by this evaluation. " + plan["external_baseline"] + "\n")
(HERE / "REPORT.md").write_text("\n".join(md))
if "evaluation" in R:
    ids = list(comps)
    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    for i, cid in enumerate(ids):
        for j, m in enumerate(("ridge", "gbm")):
            f = comps[cid]["families"][m]
            y = i + (0.14 if j else -0.14)
            lo, hi = f["interval_95_descriptive"]
            rel = f["relative_mse_reduction"]
            scale = f["mse_b"]
            ax.plot([-hi / scale * 100, -lo / scale * 100], [y, y], color=BLUE if j == 0 else ORANGE, lw=1.6)
            ax.plot([100 * rel], [y], "o", color=BLUE if j == 0 else ORANGE, mfc=(BLUE if j == 0 else ORANGE) if f["classification"] == "criterion met" else "white", ms=6)
    ax.axvline(0, color=GREY, lw=0.8)
    ax.axvline(5, color=GREY, lw=0.8, ls="--")
    ax.set_yticks(range(len(ids)))
    ax.set_yticklabels([f"{c}: {comps[c]['test']}" for c in ids], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("relative MSE reduction of the richer set (%)")
    ax.plot([], [], "o", color=BLUE, label="ridge")
    ax.plot([], [], "o", color=ORANGE, label="GBM")
    ax.legend(fontsize=8, frameon=False)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    fig.tight_layout()
    fig.savefig(HERE / "figures/met_eval_comparisons.png", dpi=120)
    plt.close(fig)
    M = np.array([comps[c]["families"][m]["block_deltas"] for c in ids for m in ("ridge", "gbm")])
    lim = np.nanmax(np.abs(M))
    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    im = ax.imshow(-M, cmap="RdBu", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_yticks(range(len(M)))
    ax.set_yticklabels([f"{c} {m}" for c in ids for m in ("ridge", "gbm")], fontsize=7)
    ax.set_xticks(range(10))
    ax.set_xticklabels([str(k + 1) for k in range(10)])
    ax.set_xlabel("outer position block")
    fig.colorbar(im, ax=ax, label="MSE reduction (blue: richer set better)")
    fig.tight_layout()
    fig.savefig(HERE / "figures/met_eval_block_consistency.png", dpi=120)
    plt.close(fig)
print("report written")
