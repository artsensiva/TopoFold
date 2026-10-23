"""Presentation-only redraw of the KRAS-STRUCT-2 figure from the frozen outputs (no metric is recomputed).

    .venv/bin/python benchmarks/audit/kras_struct_g12d/plot_figure.py

Reads results/KRAS-STRUCT-2_table.csv and the adapter exports in build/ (for the Calpha trace). The frozen runner's
own figure, results/KRAS-STRUCT-2_figure.png, is kept unchanged; this writes results/KRAS-STRUCT-2_figure_v2.png.
"""

import csv
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

HERE = Path(__file__).resolve().parent
INK, MUTED, GRID, BLUE, BAND, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6", "#c7dcf5", "#fcfcfb"
UNITS = {"d_kappa": "rad", "d_tau": "rad", "d_i_i2": "Å", "d_i_i3": "Å", "rmsd3": "Å", "rmsd4": "Å"}
TITLES = {"d_kappa": "|Δκ|, 3-residue", "d_tau": "|Δτ| (circular), 4-residue", "d_i_i2": "|Δd(i−1,i+1)|",
          "d_i_i3": "|Δd(i−1,i+2)|", "rmsd3": "local Cα RMSD, 3-residue", "rmsd4": "local Cα RMSD, 4-residue"}


def kabsch_fit(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - cP).T @ (Q - cQ))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return (P - cP) @ (Vt.T @ np.diag([1.0, 1.0, d]) @ U.T).T + cQ


def f(v):
    return np.nan if v == "" else float(v)


def main():
    rows = list(csv.DictReader(open(HERE / "results/KRAS-STRUCT-2_table.csv")))
    resid = np.array([int(r["resid"]) for r in rows])
    P = np.load(HERE / "build/5US4_A.npz")["coordinates_angstrom"][0]
    Q = np.load(HERE / "build/5USJ_A.npz")["coordinates_angstrom"][0]
    assert [r["resid"] for r in json.loads((HERE / "build/5USJ_A.json").read_text())["residues"]] == list(resid)
    Pf = kabsch_fit(P, Q)
    c = Q.mean(0)
    Vt = np.linalg.svd(Q - c)[2]
    a2, b2 = (Pf - c) @ Vt[:2].T, (Q - c) @ Vt[:2].T
    r3 = np.array([f(r["rmsd3"]) for r in rows])

    fig = plt.figure(figsize=(14, 8.6), dpi=150, facecolor=SURF)
    gs = fig.add_gridspec(3, 3, width_ratios=[2.1, 1, 1], wspace=0.25, hspace=0.55)
    ax0 = fig.add_subplot(gs[:, 0])
    cmap = LinearSegmentedColormap.from_list("seq", ["#cfe0f6", "#2a78d6", "#0b2f66"])
    segs, vals = [], []
    for i in range(len(resid) - 1):
        if resid[i + 1] - resid[i] == 1:
            segs.append([b2[i], b2[i + 1]])
            vals.append(np.nanmean([r3[i], r3[i + 1]]) if np.isfinite([r3[i], r3[i + 1]]).any() else 0.0)
    ax0.plot(*a2[:59].T, color="#9d9c96", lw=1.1, zorder=1, label="5US4 A (GDP), superposed")
    ax0.plot(*a2[59:].T, color="#9d9c96", lw=1.1, zorder=1)
    lc = LineCollection(segs, cmap=cmap, array=np.array(vals), linewidths=3.0, zorder=2, capstyle="round")
    lc.set_clim(0, max(vals))
    ax0.add_collection(lc)
    ax0.plot(*b2[58:60].T, ls=":", color=MUTED, lw=1.2, label="60–67 unobserved (excluded)")
    for r, off in ((1, (-14, -4)), (30, (6, -2)), (36, (6, 4)), (59, (6, 0)), (68, (6, 0)), (106, (6, 2)), (166, (-24, 4))):
        j = int(np.where(resid == r)[0][0])
        ax0.annotate(str(r), b2[j], fontsize=8, color=INK, xytext=off, textcoords="offset points")
    ax0.set_aspect("equal")
    ax0.set_xticks([])
    ax0.set_yticks([])
    for s in ax0.spines.values():
        s.set_visible(False)
    ax0.legend(fontsize=7.5, frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
    ax0.set_anchor("N")
    cb = fig.colorbar(lc, ax=ax0, orientation="horizontal", fraction=0.035, pad=0.01, aspect=40)
    cb.set_label("5USJ A trace coloured by 3-residue local Cα RMSD vs 5US4 A (Å)", fontsize=8, color=INK)
    cb.ax.tick_params(labelsize=7, colors=MUTED)
    cb.outline.set_visible(False)
    ax0.set_title("A  Primary pair after global Cα superposition\n(projection on the two principal axes of 5USJ A)",
                  fontsize=9, loc="left", color=INK, pad=22)

    for j, k in enumerate(("d_kappa", "d_tau", "d_i_i2", "d_i_i3", "rmsd3", "rmsd4")):
        ax = fig.add_subplot(gs[j // 2, 1 + j % 2])
        ok = np.array([r[k] != "" for r in rows])
        x = resid[ok].astype(float)
        col = lambda name: np.array([f(r[name]) for r in rows])[ok]  # noqa: E731
        for n_seg, seg in enumerate(np.split(np.arange(len(x)), np.where(np.diff(x) > 1)[0] + 1)):
            first = n_seg == 0
            ax.fill_between(x[seg], col(f"{k}_cross_min")[seg], col(f"{k}_cross_max")[seg], color=BAND, lw=0,
                            label="range over 12 GDP×GNP chain pairings" if first else None)
            ax.plot(x[seg], col(f"{k}_copy_max")[seg], color=MUTED, lw=0.9, ls=(0, (3, 2)),
                    label="max over 16 same-entry copy pairs" if first else None)
            ax.plot(x[seg], col(k)[seg], color=BLUE, lw=1.4, label="primary: 5US4 A vs 5USJ A" if first else None)
        ax.axvspan(59.5, 67.5, color=GRID, lw=0)
        ax.set_title(f"{'BCDEFG'[j]}  {TITLES[k]} ({UNITS[k]})", fontsize=8.5, loc="left", color=INK)
        ax.tick_params(labelsize=7, colors=MUTED)
        ax.grid(axis="y", color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(MUTED)
        ax.set_xlim(0, 168)
        ax.set_ylim(bottom=0)
        if j >= 4:
            ax.set_xlabel("window centre residue (κ: i; τ: hinge i→i+1)", fontsize=7.5, color=MUTED)
        if j == 0:
            handles, labels = ax.get_legend_handles_labels()
    fig.suptitle("KRAS-STRUCT-2 · KRAS G12D crystal structures, GDP (5US4) vs GMPPNP (5USJ), no pharmacological inhibitor "
                 "modelled · static and descriptive; grey band 60–67 unobserved", fontsize=10, x=0.01, ha="left", color=INK)
    fig.legend(handles, labels, fontsize=7.5, frameon=False, loc="upper right", bbox_to_anchor=(0.99, 0.965), ncol=3)
    fig.savefig(HERE / "results/KRAS-STRUCT-2_figure_v2.png", bbox_inches="tight", facecolor=SURF)


if __name__ == "__main__":
    main()
