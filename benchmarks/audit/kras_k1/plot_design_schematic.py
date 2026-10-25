"""Schematic of the proposed KRAS K1 design and its decision branches (a DRAFT design, not a result).

    PYTHONPATH=python .venv/bin/python benchmarks/audit/kras_k1/plot_design_schematic.py

Reads protocol_params_DRAFT.json (design parameters) and results/k1_provenance_checks.json plus
../next_app_feasibility/results/feasibility_metadata.json (recorded counts and sizes). Writes
docs/figures/kras_k1/K1_design_and_decisions_DRAFT_2026-10-25.png and a provenance record with the SHA-256 of the inputs,
this generator and the output. No scientific value is typed into the figure.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PARAMS = HERE / "protocol_params_DRAFT.json"
PROV = HERE / "results/k1_provenance_checks.json"
FEAS = ROOT / "benchmarks/audit/next_app_feasibility/results/feasibility_metadata.json"
OUT_DIR = ROOT / "docs/figures/kras_k1"
STEM = "K1_design_and_decisions_DRAFT_2026-10-25"
INK, MUTED = "#0b0b0b", "#52514e"
STEP, GATE, CLOSE, GO = "#dce9f8", "#f6e3c8", "#e6e5e1", "#d3ead5"
EDGE = {"step": "#2a78d6", "gate": "#d9822b", "close": "#8a8984", "go": "#2e8b57"}


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while b := f.read(1 << 20):
            h.update(b)
    return h.hexdigest()


def facts(params: dict, prov: dict, feas: dict) -> dict:
    """Every number shown in the figure, from the parameter file or the recorded results."""
    reps = feas["kras_pantsar_zenodo"]["replicas_per_system"]
    gdp, gtp = len(reps["GDP_G12D"]), len(reps["GTP_G12D"])
    sub = feas["kras_subsets"]["G12D GDP + GTP + topology files"]
    s = prov["g12d_summary"]
    e = params["endpoint"]
    return {"gdp": gdp, "gtp": gtp, "files": sub["files"], "mb": round(sub["bytes"] / 1e6), "natoms": s["natoms"][0],
            "subsets": math.comb(gtp, gdp), "w0": params["window_ns"][0], "w1": params["window_ns"][1],
            "min_gdp": params["min_retained_runs"]["GDP"], "min_gtp": params["min_retained_runs"]["GTP"],
            "pos_min": e["min_positive"], "pos_max_pct": round(100 * e["max_positive_fraction"]), "segs": e["min_segments"],
            "chains": e["chains"], "delta": params["delta_precision"], "corr": params["correspondence_min_precision"],
            "bal_pct": round(100 * params["balanced_subset_min_fraction"]), "jac": params["redundancy_jaccard_min"],
            "mm": e["min_displacement_A"]}


def draw(f: dict, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    fig, ax = plt.subplots(figsize=(11.0, 5.6), dpi=180)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 56)
    ax.axis("off")

    def box(x, y, w, h, text, kind, size=7.4, bold=None):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=0.8",
                                    facecolor={"step": STEP, "gate": GATE, "close": CLOSE, "go": GO}[kind],
                                    edgecolor=EDGE[kind], linewidth=1.1))
        if bold:
            ax.text(x + w / 2, y + h - 1.4, bold, ha="center", va="top", fontsize=size + 0.4, fontweight="bold", color=INK)
            ax.text(x + w / 2, y + h - 4.1, text, ha="center", va="top", fontsize=size, color=INK, linespacing=1.25)
        else:
            ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size, color=INK, linespacing=1.25)

    def arrow(x0, y0, x1, y1, label=None, lx=0.0, ly=1.2):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.0, shrinkA=0, shrinkB=0))
        if label:
            ax.text((x0 + x1) / 2 + lx, (y0 + y1) / 2 + ly, label, fontsize=6.8, color=MUTED, ha="center", va="center")

    top, h = 37, 15
    box(0.5, top, 21, h, f"5US4 (GDP) vs 5USJ (GMPPNP),\n{f['chains']} chains: per-residue Cα\ndisplacement defines the\nendpoint set E. No MD read.",
        "step", bold="1  Endpoint from crystals")
    box(25.5, top, 18, h, f"{f['pos_min']} ≤ |E| ≤ {f['pos_max_pct']} % of scored\nresidues, and ≥ {f['segs']}\nseparate segments?", "gate", bold="Gate: endpoint")
    box(47.5, top, 22, h, f"{f['files']} files, {f['mb']} MB, MD5-checked.\nG12D: {f['gdp']} GDP, {f['gtp']} GTP runs;\nN, Cα, C only ({f['natoms']} atoms).",
        "step", bold="2  Fetch MD")
    box(73.5, top, 26, h, f"identity table; window {f['w0']}–{f['w1']} ns;\n≥ {f['min_gdp']} GDP and ≥ {f['min_gtp']} GTP runs kept;\nno PBC jumps in the backbone.",
        "gate", bold="3  Gates on the MD files")
    arrow(21.7, top + h / 2, 25.2, top + h / 2)
    arrow(43.8, top + h / 2, 47.2, top + h / 2, "yes", ly=1.6)
    arrow(69.8, top + h / 2, 73.2, top + h / 2)
    box(25.5, 24, 18, 9.5, "Close the application:\nthe endpoint cannot\ndiscriminate methods.\nNo MD downloaded.", "close", size=7)
    arrow(34.5, top, 34.5, 33.9, "no", lx=1.4, ly=0)
    box(83.0, 24, 16.5, 9.5, "Stop, record the cause.\nNo rescue: no other\nchain, window, threshold\nor system.", "close", size=6.8)
    arrow(91.2, top, 91.2, 33.9, "any gate fails", lx=-6.2, ly=0)
    # pass path: from gate 3 down and left to step 4
    ax.plot([77.0, 77.0], [top, 22.6], color=MUTED, lw=1.0)
    ax.plot([77.0, 14.0], [22.6, 22.6], color=MUTED, lw=1.0)
    arrow(14.0, 22.6, 14.0, 21.4)
    ax.text(60.0, 23.5, "all gates pass", fontsize=6.8, color=MUTED, ha="center", va="center")

    box(0.5, 2.0, 27, 19, f"5 methods ranked per residue: κ, τ\n(TopoFold); d(i,i+2), d(i,i+3), RMSF\n(Cα baselines). φ/ψ is shown apart.\nTop-k residues, k = |E|, equal for all;\nprecision at k against E.\nChecks: {f['subsets']} balanced GTP subsets,\nleave-one-run-out, two time halves.",
        "step", size=7, bold="4  Compare at equal budget")
    outs = [("C1", f"No method reaches\nprecision ≥ {f['corr']}", "Close: correspondence\nnot demonstrated", "close", 14.5),
            ("C2", f"Gap < {f['delta']}, or the gain is\nredundant with its matched\nbaseline or in one segment", "Close: no separate\nadvantage for TopoFold", "close", 14.5),
            ("C4", f"Gap ≥ {f['delta']}, but not the\nsame sign across subsets,\nhalves or leave-one-out", "Close: inconclusive\n(unstable)", "close", 14.5),
            ("C3", f"All of: gap ≥ {f['delta']}; same sign in\n≥ {f['bal_pct']} % of subsets, both halves and\nevery leave-one-out run; Jaccard\nwith matched baseline < {f['jac']};\n≥ {f['segs']} endpoint segments",
             "Justifies a separately frozen\nconfirmation on an independent\nsystem. Not a product claim.", "go", 22.0)]
    x = 29.8
    for tag, cond, res, kind, w in outs:
        box(x, 2.0, w, 19, f"{cond}\n\n→ {res}", kind, size=6.5, bold=f"Outcome {tag}")
        x += w + 0.9
    arrow(27.8, 11.5, 29.5, 11.5)
    ax.text(0.5, 54.6, "KRAS K1 (proposed, not frozen): crystal-correspondence test of TopoFold against Cα baselines, "
            "G12D GDP vs GTP", fontsize=9.6, color=INK, va="center", ha="left", fontweight="bold")
    ax.text(0.5, 0.1, "Schematic of a proposed design, not a result. Outcomes are evaluated in the order C1, C2, then C3 or C4. "
            "Parameters: protocol_params_DRAFT.json; counts: recorded provenance checks.", fontsize=6.6, color=MUTED, va="center", ha="left")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def main() -> int:
    params, prov, feas = (json.loads(p.read_text()) for p in (PARAMS, PROV, FEAS))
    f = facts(params, prov, feas)
    png = OUT_DIR / f"{STEM}.png"
    draw(f, png)
    import matplotlib
    rec = {"id": "KRAS-K1-design-schematic", "date": "2026-10-25",
           "kind": "schematic of a proposed design (no measured result is plotted)",
           "caption_file": f"{STEM}_caption.md",
           "inputs_sha256": {str(p.relative_to(ROOT)): sha256_file(p) for p in (PARAMS, PROV, FEAS)},
           "generator": {"path": str(Path(__file__).resolve().relative_to(ROOT)), "sha256": sha256_file(Path(__file__).resolve())},
           "output": {"path": f"docs/figures/kras_k1/{STEM}.png", "bytes": png.stat().st_size, "sha256": sha256_file(png)},
           "software": {"python": sys.version.split()[0], "matplotlib": matplotlib.__version__},
           "reproduce": "PYTHONPATH=python .venv/bin/python benchmarks/audit/kras_k1/plot_design_schematic.py",
           "numbers_shown": f}
    (OUT_DIR / f"{STEM}_provenance.json").write_text(json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps({"png_bytes": rec["output"]["bytes"], **f}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
