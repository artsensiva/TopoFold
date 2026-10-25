"""Supplementary visualization of the closed ABL-MYR-1b result (2026-10-25): residue-level coverage and the frozen ±2 rule.

    .venv/bin/python benchmarks/audit/abl_myr_1b/plot_residue_coverage_supplement.py

Presentation only. Reads results/ABL-MYR-1b_residue_summary.json (residue numbers) and results/ABL-MYR-1b_results.json
(topology indices and the recorded outcome sets), checks that they agree, and draws the recorded sets. No score,
threshold or window is recomputed; the ±2 neighbourhood is derived from the recorded R_B by the frozen definition
(run_abl_myr.classify_outcome) only to draw it, and is checked against the recorded R_T_outside_N2_R_B and
R_B_outside_N2_R_T. Writes figures/ABL-MYR-1b_residue_coverage_SUPPLEMENT_2026-10-25.png and a provenance record with the
SHA-256 of the inputs, this generator and the output. No frozen file is modified.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / "results/ABL-MYR-1b_results.json"
SUMMARY = HERE / "results/ABL-MYR-1b_residue_summary.json"
FIG_DIR = HERE / "figures"
STEM = "ABL-MYR-1b_residue_coverage_SUPPLEMENT_2026-10-25"
SPAN = (254, 504)            # represented kinase-domain residues (topology), isoform 1b numbering
ALPHA_I = (518, 531)         # alphaI residues named in the source study; not in the topology
AXIS = (250, 535)
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, ORANGE_LIGHT, PURPLE, HATCH = "#2a78d6", "#d9822b", "#f3d7b8", "#7b5ea7", "#b9b8b3"


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while b := f.read(1 << 20):
            h.update(b)
    return h.hexdigest()


def segments(nums) -> list:
    """Contiguous runs of integers as (first, last)."""
    out, cur = [], []
    for n in sorted(set(nums)):
        if cur and n == cur[-1] + 1:
            cur.append(n)
        else:
            if cur:
                out.append((cur[0], cur[-1]))
            cur = [n]
    if cur:
        out.append((cur[0], cur[-1]))
    return out


def neighbourhood(residues, k: int = 2) -> set:
    """N_k(S): every residue within k of a member of S (the frozen N2 definition)."""
    return {j for i in residues for j in range(i - k, i + k + 1)}


def to_resid(indices, index_to_resid: dict) -> list:
    return sorted(index_to_resid[i] for i in indices)


def load_and_check(res: dict, summ: dict) -> dict:
    """Residue-number sets from the summary, checked against the results file (indices -> resid) and the recorded outside sets."""
    idx = {r["index"]: r["resid"] for r in res["residues"]}
    if sorted(idx.values()) != list(range(SPAN[0], SPAN[1] + 1)):
        raise ValueError("topology residues are not the contiguous represented span")
    o = res["outcome"]
    sets = {"R_T": to_resid(o["R_T"], idx), "R_B": to_resid(o["R_B"], idx), "R_phi_psi": to_resid(o["R_phi_psi"], idx)}
    for k, v in sets.items():
        if v != summ[k]["residues"]:
            raise ValueError(f"{k}: results.json (via topology indices) and the residue summary disagree")
    n2_b, n2_t = neighbourhood(sets["R_B"]), neighbourhood(sets["R_T"])
    outside_t = sorted(set(sets["R_T"]) - n2_b)
    outside_b = sorted(set(sets["R_B"]) - n2_t)
    if outside_t != summ["R_T_outside_N2_R_B"]["residues"] or outside_b != summ["R_B_outside_N2_R_T"]["residues"]:
        raise ValueError("the ±2 rule applied to the recorded sets does not reproduce the recorded outside sets")
    if outside_t != to_resid(o["R_T_outside_N2_R_B"], idx) or outside_b != to_resid(o["R_B_outside_N2_R_T"], idx):
        raise ValueError("recorded outside sets differ between results.json and the summary")
    a = summ["alpha_I"]
    if a["represented"] or a["not_represented"] != list(range(ALPHA_I[0], ALPHA_I[1] + 1)):
        raise ValueError("alphaI record changed")
    return {**sets, "N2_R_B": sorted(n2_b), "T_outside": outside_t, "B_outside": outside_b,
            "label": o["label"], "step": o["decision_rule_step"], "span": tuple(summ["topology_residue_range"])}


def draw(d: dict, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    rows = [("TopoFold κ, τ\n(R_T, 34 residues)", d["R_T"], BLUE, 4),
            ("Cα baselines\n(R_B, 40 residues)", d["R_B"], ORANGE, 3),
            ("Baseline ±2 neighbourhood\nN2(R_B), the frozen rule", [r for r in d["N2_R_B"] if SPAN[0] <= r <= SPAN[1]], ORANGE_LIGHT, 2),
            ("TopoFold residues outside\nN2(R_B)", d["T_outside"], BLUE, 1),
            ("φ/ψ comparator (stronger\ninformation; not in the rule)", d["R_phi_psi"], PURPLE, 0)]
    fig, ax = plt.subplots(figsize=(10.0, 3.9), dpi=200)
    ax.axvspan(SPAN[1] + 0.5, AXIS[1], ymin=0, ymax=1, facecolor="none", edgecolor=HATCH, hatch="////", linewidth=0, zorder=0)
    ax.axvspan(ALPHA_I[0] - 0.5, ALPHA_I[1] + 0.5, ymin=0, ymax=1, facecolor="#e9e8e4", edgecolor="none", zorder=0.5, alpha=0.9)
    for _, res, colour, y in rows:
        for a, b in segments(res):
            ax.broken_barh([(a - 0.5, b - a + 1)], (y - 0.32, 0.64), facecolors=colour, edgecolors="none", zorder=2)
    for y in range(5):
        ax.axhline(y, color=GRID, linewidth=0.6, zorder=1)
    ax.text((SPAN[0] + SPAN[1]) / 2, 1, "none (0 residues)", va="center", ha="center", fontsize=8.5, color=INK, zorder=3,
            style="italic")
    (b317,) = d["B_outside"]
    ax.annotate(f"{b317}: baseline-only (no TopoFold residue within ±2)", xy=(b317, 3.3), xytext=(b317 + 4, 3.5),
                fontsize=7.2, color=ORANGE, ha="left", va="center", arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.8))
    ax.plot([SPAN[0] - 0.5, SPAN[1] + 0.5], [4.78, 4.78], color=INK, linewidth=1.4, clip_on=False)
    ax.text((SPAN[0] + SPAN[1]) / 2, 4.98, f"represented span {SPAN[0]}–{SPAN[1]} (251 residues; N, Cα, C atoms)", ha="center",
            va="bottom", fontsize=8, color=INK)
    white = dict(facecolor="white", edgecolor="none", pad=1.0, alpha=0.85)
    ax.text((SPAN[1] + 0.5 + ALPHA_I[0] - 0.5) / 2, 2, "beyond the topology", rotation=90, ha="center", va="center",
            fontsize=7, color=MUTED, bbox=white, zorder=3)
    ax.text((ALPHA_I[0] + ALPHA_I[1]) / 2, 2, f"αI {ALPHA_I[0]}–{ALPHA_I[1]}: unavailable", rotation=90, ha="center",
            va="center", fontsize=7.5, color=INK, zorder=3)
    ax.set_xlim(*AXIS)
    ax.set_ylim(-0.6, 4.7)
    ax.set_yticks([y for *_, y in rows])
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.set_xticks(range(275, 526, 25))
    ax.set_xlabel("residue number (UniProt P00519-2, isoform 1b numbering)", fontsize=8.5)
    ax.tick_params(axis="x", labelsize=8)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.suptitle("ABL-MYR-1b (closed): residues flagged by TopoFold and by the Cα baselines, and the frozen ±2 rule",
                 fontsize=9.5, x=0.01, ha="left", color=INK)
    handles = [Patch(facecolor=BLUE, label="TopoFold"), Patch(facecolor=ORANGE, label="Cα baseline"),
               Patch(facecolor=ORANGE_LIGHT, label="within ±2 of a baseline residue"),
               Patch(facecolor=PURPLE, label="φ/ψ comparator"),
               Patch(facecolor="none", edgecolor=HATCH, hatch="////", label="beyond the topology")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.6, 0.045), ncol=5, fontsize=7.2, frameon=False,
               handletextpad=0.5, columnspacing=1.2)
    fig.text(0.01, 0.005, "Descriptive, post-hoc, simulation-based; values drawn from results/ABL-MYR-1b_results.json "
             "(outcome sets) and ABL-MYR-1b_residue_summary.json.", fontsize=6.8, color=MUTED, ha="left", va="bottom")
    fig.subplots_adjust(left=0.2, right=0.985, top=0.83, bottom=0.27)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)


def main() -> int:
    res, summ = json.loads(RES.read_text()), json.loads(SUMMARY.read_text())
    d = load_and_check(res, summ)
    png = FIG_DIR / f"{STEM}.png"
    draw(d, png)
    import matplotlib
    prov = {"id": "ABL-MYR-1b-SUPPLEMENT-residue-coverage", "date": "2026-10-25",
            "kind": "measured results (recorded sets); schematic elements: the ±2 band and the span marks only",
            "caption_file": f"figures/{STEM}_caption.md",
            "inputs_sha256": {"results/ABL-MYR-1b_results.json": sha256_file(RES),
                              "results/ABL-MYR-1b_residue_summary.json": sha256_file(SUMMARY)},
            "generator": {"path": "plot_residue_coverage_supplement.py", "sha256": sha256_file(Path(__file__).resolve())},
            "output": {"path": f"figures/{STEM}.png", "bytes": png.stat().st_size, "sha256": sha256_file(png)},
            "software": {"python": sys.version.split()[0], "matplotlib": matplotlib.__version__},
            "reproduce": "PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1b/plot_residue_coverage_supplement.py",
            "plotted_counts": {"R_T": len(d["R_T"]), "R_B": len(d["R_B"]), "R_phi_psi": len(d["R_phi_psi"]),
                               "T_outside_N2_R_B": len(d["T_outside"]), "B_outside_N2_R_T": len(d["B_outside"])},
            "checks": ["results.json indices map to resid 254-504 contiguously",
                       "R_T, R_B, R_phi_psi equal between results.json (via indices) and the residue summary",
                       "the frozen ±2 definition applied to the recorded sets reproduces both recorded outside sets",
                       "alphaI 518-531 recorded as not represented"],
            "outcome_label_recorded": d["label"], "decision_rule_step_recorded": d["step"],
            "not_modified": "no frozen specification, report, measurement or result file is written by this generator",
            "storage_note": "hashes and the ~/Downloads copy of the archives are not an independent backup (same device); "
                            "the gap recorded in CLOSEOUT_2026-10-25.md is unchanged"}
    (FIG_DIR / f"{STEM}_provenance.json").write_text(json.dumps(prov, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps({"png_bytes": prov["output"]["bytes"], "counts": prov["plotted_counts"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
