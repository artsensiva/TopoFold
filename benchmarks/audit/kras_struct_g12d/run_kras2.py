"""KRAS-STRUCT-2: descriptive comparison of GDP-bound (5US4) and GMPPNP-bound (5USJ) KRAS G12D crystal structures.

    .venv/bin/python benchmarks/audit/kras_struct_g12d/run_kras2.py freeze
    .venv/bin/python benchmarks/audit/kras_struct_g12d/run_kras2.py run

Reuses the KRAS-STRUCT-1 adapter call, metric definitions and QCP window RMSD unchanged (imported, not copied).
Static structures: direct geometric differences only (see run_spec.json, permitted_interpretation).
"""

import csv
import datetime as dt
import hashlib
import itertools
import json
import resource
import statistics
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/kras_struct_cmp"))
import run_kras as k1  # noqa: E402
from topofold_inputs.compare import spearman  # noqa: E402
from topofold_inputs.features import export_features, features_from_export  # noqa: E402

SPEC = HERE / "run_spec.json"
FREEZE = HERE / "results/KRAS-STRUCT-2_freeze.json"
CODE = ["benchmarks/audit/kras_struct_g12d/run_kras2.py"] + k1.CODE
METRICS = k1.METRICS
FOUR = ("d_tau", "d_i_i3", "rmsd4")
UNITS = {"d_kappa": "rad", "d_tau": "rad", "d_i_i2": "Å", "d_i_i3": "Å", "rmsd3": "Å", "rmsd4": "Å"}
RESIDS = list(range(1, 60)) + list(range(68, 167))
S1_TABLE = ROOT / "benchmarks/audit/kras_struct_cmp/results/KRAS-STRUCT-1_table.csv"


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def freeze():
    if FREEZE.exists():
        print("freeze exists")
        return 1
    rec = {"id": "KRAS-STRUCT-2", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "spec_sha256": sha256(SPEC), "code_sha256": {c: sha256(ROOT / c) for c in CODE},
           "struct1_table_sha256": sha256(S1_TABLE)}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=1))
    return 0


def extract_chain(spec, state, chain, build):
    corr = build / f"corr_{chain}.json"
    corr.write_text(json.dumps([{"chain": chain, "resid": r, "label": f"KRAS:{r}"} for r in RESIDS]))
    run_id = f"{spec['structures'][state]['pdb']}_{chain}"
    rc, err = k1.adapter(spec["structures"][state]["file"], run_id, build, corr, "reject")
    if rc != 0:
        raise SystemExit(f"STOP: {run_id}: adapter rejected ({err})")
    feats, rep = features_from_export(build / f"{run_id}.npz", build / f"{run_id}.json")
    export_features(feats, build / f"{run_id}.npz", build / f"{run_id}.json", build, run_id)
    if [r["resid"] for r in rep["residues"]] != RESIDS or not all(r["status"] == "present" for r in rep["residues"]):
        raise SystemExit(f"STOP: {run_id}: mapping check failed")
    return {"rep": rep, "feats": feats, "coords": np.load(build / f"{run_id}.npz")["coordinates_angstrom"],
            "npz": build / f"{run_id}.npz", "id": run_id}


def top10(m, k):
    return [i for i in sorted((i for i, v in enumerate(m[k]) if v is not None), key=lambda i: -m[k][i])[:10]]


def window_label(i, four):
    return f"{RESIDS[i - 1]}-{RESIDS[i + (2 if four else 1)]}"


def kabsch_fit(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - cP).T @ (Q - cQ))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return (P - cP) @ R.T + cQ


def run():
    spec = json.loads(SPEC.read_text())
    frz = json.loads(FREEZE.read_text())
    bad = [c for c, v in frz["code_sha256"].items() if sha256(ROOT / c) != v]
    if frz["spec_sha256"] != sha256(SPEC) or bad or frz["struct1_table_sha256"] != sha256(S1_TABLE):
        raise SystemExit(f"STOP: spec, code or STRUCT-1 table differs from the freeze: {bad}")
    for st in spec["structures"].values():
        if sha256(ROOT / st["file"]) != st["sha256"]:
            raise SystemExit(f"STOP: {st['pdb']} hash differs")
    t0 = time.time()
    build = HERE / "build"
    build.mkdir(exist_ok=True)
    ch = {("GDP", c): extract_chain(spec, "GDP", c, build) for c in spec["structures"]["GDP"]["chains"]}
    ch.update({("GNP", c): extract_chain(spec, "GNP", c, build) for c in spec["structures"]["GNP"]["chains"]})
    names = [r["resname"] for r in ch[("GDP", "A")]["rep"]["residues"]]
    if any([r["resname"] for r in v["rep"]["residues"]] != names for v in ch.values()):
        raise SystemExit("STOP: residue names differ between chains")
    breaks = {v["id"]: v["rep"]["breaks"] for v in ch.values()}

    def pair(a, b, tag):
        return k1.metrics({"GDP": ch[a], "GNP": ch[b]}, build, tag)

    cross = {f"{a}-{b}": pair(("GDP", a), ("GNP", b), f"x_{a}{b}")
             for a in spec["structures"]["GDP"]["chains"] for b in spec["structures"]["GNP"]["chains"]}
    copies = {}
    for st in ("GDP", "GNP"):
        for a, b in itertools.combinations(spec["structures"][st]["chains"], 2):
            copies[f"{spec['structures'][st]['pdb']}:{a}-{b}"] = pair((st, a), (st, b), f"c_{st}{a}{b}")
    prim = cross["A-A"][0]
    n = len(RESIDS)

    def vals(store, k, i):
        return [m[k][i] for m, _, _ in store.values() if m[k][i] is not None]

    rows = []
    for i in range(n):
        row = {"index": i, "resid": RESIDS[i], "resname": names[i],
               "window3": window_label(i, False) if 1 <= i <= n - 2 else "",
               "window4": window_label(i, True) if 1 <= i <= n - 3 else ""}
        for k in METRICS:
            xs, cs = vals(cross, k, i), vals(copies, k, i)
            row[k] = prim[k][i]
            row[f"{k}_cross_min"] = min(xs) if xs else None
            row[f"{k}_cross_median"] = statistics.median(xs) if xs else None
            row[f"{k}_cross_max"] = max(xs) if xs else None
            row[f"{k}_copy_max"] = max(cs) if cs else None
        rows.append(row)
    excluded = {"three_residue": [], "four_residue": []}
    for i in range(n):
        for four, key in ((False, "three_residue"), (True, "four_residue")):
            lo, hi = i - 1, i + (2 if four else 1)
            if prim["d_tau" if four else "d_kappa"][i] is None and 0 <= lo and hi < n:
                excluded[key].append({"window": window_label(i, four),
                                      "cause": "numbering_discontinuity 59->68" if RESIDS[hi] - RESIDS[lo] != hi - lo
                                      else "invalid core output"})
    recurrence = {}
    for k in METRICS:
        cnt = {}
        for m, _, _ in cross.values():
            for i in top10(m, k):
                cnt[i] = cnt.get(i, 0) + 1
        recurrence[k] = [{"window": window_label(i, k in FOUR), "n_pairings_top10": c, "primary_value": prim[k][i],
                          "in_primary_top10": i in top10(prim, k)}
                         for i, c in sorted(cnt.items(), key=lambda x: (-x[1], x[0]))]
    pairs_c = [("d_kappa", "d_i_i2"), ("d_kappa", "rmsd3"), ("d_i_i2", "rmsd3"),
               ("d_tau", "d_i_i3"), ("d_tau", "rmsd4"), ("d_i_i3", "rmsd4")]

    def rho(m, x, y):
        p = [(m[x][i], m[y][i]) for i in range(n) if m[x][i] is not None and m[y][i] is not None]
        return spearman([a for a, _ in p], [b for _, b in p]), len(p)

    agreement = {}
    for x, y in pairs_c:
        r0, nw = rho(prim, x, y)
        rs = [rho(m, x, y)[0] for m, _, _ in cross.values()]
        agreement[f"{x}~{y}"] = {"primary": r0, "n_windows": nw, "min_over_12": min(rs), "max_over_12": max(rs)}
    # cross-study comparison with the frozen KRAS-STRUCT-1 table
    s1 = {int(r["resid"]): r for r in csv.DictReader(open(S1_TABLE))}
    cross_study = {}
    s2_by_resid = {RESIDS[i]: i for i in range(n)}
    for k in METRICS:
        common = [(float(s1[r][k]), prim[k][i]) for r, i in s2_by_resid.items()
                  if r in s1 and s1[r][k] != "" and prim[k][i] is not None]
        s1_top = sorted(((float(s1[r][k]), r) for r in s2_by_resid if r in s1 and s1[r][k] != ""
                         and prim[k][s2_by_resid[r]] is not None), reverse=True)[:10]
        s2_top = sorted(((prim[k][i], r) for r, i in s2_by_resid.items() if r in s1 and s1[r][k] != ""
                         and prim[k][i] is not None), reverse=True)[:10]
        cross_study[k] = {"n_shared_windows": len(common),
                          "spearman": spearman([a for a, _ in common], [b for _, b in common]),
                          "struct1_top10_centres": sorted(r for _, r in s1_top),
                          "struct2_top10_centres": sorted(r for _, r in s2_top),
                          "overlap": sorted({r for _, r in s1_top} & {r for _, r in s2_top})}
    ru, rc = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    out = {"id": "KRAS-STRUCT-2", "freeze_utc": frz["frozen_utc"], "spec_sha256": frz["spec_sha256"],
           "mapping": {"n_residues": n, "residues": "1-59, 68-166", "residue_names_identical_all_chains": True,
                       "breaks": breaks},
           "global_rmsd_A": {"cross": {p: v[1] for p, v in cross.items()}, "copies": {p: v[1] for p, v in copies.items()},
                             "n_ca": cross["A-A"][2]},
           "primary_top10": {k: [{"window": window_label(i, k in FOUR), "value": prim[k][i]} for i in top10(prim, k)]
                             for k in METRICS},
           "recurrence_over_12_pairings": recurrence, "rank_agreement": agreement, "exclusions": excluded,
           "cross_study_vs_struct1": cross_study,
           "resources": {"cpu_s": ru.ru_utime + ru.ru_stime + rc.ru_utime + rc.ru_stime, "wall_s": time.time() - t0},
           "core": json.loads((build / "5US4_A_features.json").read_text())["provenance"]["core"]}
    (HERE / "results/KRAS-STRUCT-2_results.json").write_text(json.dumps(out, indent=1) + "\n")
    with open(HERE / "results/KRAS-STRUCT-2_table.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else ("" if v is None else v)) for k, v in r.items()})
    figure(rows, ch, prim)
    (HERE / "results/KRAS-STRUCT-2_report.md").write_text(report(out))
    print(json.dumps({k: out[k] for k in ("global_rmsd_A", "rank_agreement", "cross_study_vs_struct1", "resources")},
                     indent=1, default=str))
    return 0


def figure(rows, ch, prim):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    ink, muted, grid, blue = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6"
    fig = plt.figure(figsize=(13, 8.2), dpi=150)
    gs = fig.add_gridspec(3, 3, width_ratios=[1.15, 1, 1], wspace=0.28, hspace=0.45)
    ax0 = fig.add_subplot(gs[:, 0])
    P, Q = ch[("GDP", "A")]["coords"][0], ch[("GNP", "A")]["coords"][0]
    Pf = kabsch_fit(P, Q)
    c = Q.mean(0)
    _, _, Vt = np.linalg.svd(Q - c)
    proj = lambda X: (X - c) @ Vt[:2].T  # noqa: E731
    a2, b2 = proj(Pf), proj(Q)
    seg_ok = [RESIDS[i + 1] - RESIDS[i] == 1 for i in range(len(RESIDS) - 1)]
    r3 = np.array([v if v is not None else np.nan for v in prim["rmsd3"]])
    segs, cols = [], []
    for i in range(len(RESIDS) - 1):
        if seg_ok[i]:
            segs.append([b2[i], b2[i + 1]])
            cols.append(np.nanmean([r3[i], r3[i + 1]]))
    ax0.plot(*a2.T, color="#b9b8b2", lw=1.0, zorder=1)
    lc = LineCollection(segs, cmap="Blues", array=np.nan_to_num(np.array(cols), nan=0.0), linewidths=2.4, zorder=2)
    lc.set_clim(0, np.nanmax(cols))
    ax0.add_collection(lc)
    gi = RESIDS.index(59)
    ax0.plot([b2[gi, 0], b2[gi + 1, 0]], [b2[gi, 1], b2[gi + 1, 1]], ls=":", color=muted, lw=1)
    for r in (1, 30, 59, 68, 106, 166):
        j = RESIDS.index(r)
        ax0.annotate(str(r), b2[j], fontsize=7.5, color=ink, xytext=(3, 3), textcoords="offset points")
    ax0.set_aspect("equal")
    ax0.set_xticks([]), ax0.set_yticks([])
    for s in ax0.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(lc, ax=ax0, orientation="horizontal", fraction=0.04, pad=0.02)
    cb.set_label("three-residue Cα RMSD, 5US4 A vs 5USJ A (Å)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax0.set_title("A  Superposed Cα traces (projection on the two\nprincipal axes of 5USJ A); grey = 5US4 A;\n"
                  "dotted = unobserved 60–67 (excluded)", fontsize=8.5, loc="left", color=ink)
    for j, k in enumerate(("d_kappa", "d_tau", "d_i_i2", "d_i_i3", "rmsd3", "rmsd4")):
        ax = fig.add_subplot(gs[j // 2, 1 + j % 2])
        x = [r["resid"] for r in rows if r[k] is not None]
        get = lambda f: np.array([np.nan if r[f] is None else r[f] for r in rows if r[k] is not None], float)  # noqa: E731
        xs = np.array(x, float)
        brk = np.where(np.diff(xs) > 1)[0] + 1
        for seg in np.split(np.arange(len(xs)), brk):
            ax.fill_between(xs[seg], get(f"{k}_cross_min")[seg], get(f"{k}_cross_max")[seg], color="#c7dcf5",
                            lw=0, label="range, 12 GDP×GNP chain pairings" if seg[0] == 0 else None)
            ax.plot(xs[seg], get(f"{k}_copy_max")[seg], color=muted, lw=1, ls="--",
                    label="max, 16 same-entry copy pairs" if seg[0] == 0 else None)
            ax.plot(xs[seg], get(k)[seg], color=blue, lw=1.4, label="primary 5US4 A vs 5USJ A" if seg[0] == 0 else None)
        ax.axvspan(59.5, 67.5, color=grid, lw=0)
        ax.set_title(f"{'BCDEFG'[j]}  {k} ({UNITS[k]})", fontsize=8.5, loc="left", color=ink)
        ax.tick_params(labelsize=7, colors=muted)
        ax.grid(axis="y", color=grid, lw=0.6)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.set_xlim(0, 168)
        if j >= 4:
            ax.set_xlabel("window centre residue", fontsize=7.5, color=muted)
        if j == 0:
            ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    fig.suptitle("KRAS-STRUCT-2: KRAS G12D GDP (5US4) vs GMPPNP (5USJ) crystal structures — static, descriptive; "
                 "shaded 60–67 unobserved", fontsize=10, x=0.01, ha="left", color=ink)
    fig.savefig(HERE / "results/KRAS-STRUCT-2_figure.png", bbox_inches="tight", facecolor="#fcfcfb")
    plt.close(fig)


def report(o):
    L = ["# KRAS-STRUCT-2: GDP-bound (5US4) vs GMPPNP-bound (5USJ) KRAS G12D, no pharmacological inhibitor modelled", "",
         "Generated from `KRAS-STRUCT-2_results.json` and `KRAS-STRUCT-2_table.csv`. **Static structures: direct "
         "geometric differences only.** Confounded by crystal form (H 3, two copies vs P 1, six copies), "
         "crystallization pH and additives (9.1 with a compound present, no density reported, vs 5.5), resolution "
         "(1.83 vs 1.94 Å), unobserved Switch-II residues 60–67 (excluded) and coordinate error. Crystal copies are "
         "not independent replicates. No ensemble, population, kinetic, functional-state or treatment-response claim.",
         "", f"- Mapping: residues {o['mapping']['residues']} ({o['mapping']['n_residues']}), identical in all 8 chains.",
         "- Global Cα RMSD (QCP, {} Cα): primary A-A {:.3f} Å; 12 GDP×GNP pairings {:.3f}–{:.3f} Å; 16 same-entry "
         "copy pairs {:.3f}–{:.3f} Å.".format(
             o["global_rmsd_A"]["n_ca"], o["global_rmsd_A"]["cross"]["A-A"], min(o["global_rmsd_A"]["cross"].values()),
             max(o["global_rmsd_A"]["cross"].values()), min(o["global_rmsd_A"]["copies"].values()),
             max(o["global_rmsd_A"]["copies"].values())),
         "- Excluded windows: " + "; ".join(f"{k}: " + ", ".join(f"{e['window']} ({e['cause']})" for e in v)
                                            for k, v in o["exclusions"].items()), "",
         "## Largest local differences, primary pair, and recurrence over the 12 chain pairings (descriptive)", "",
         "| Metric | Unit | Ten largest primary windows (value; n/12 pairings with this window in their ten largest) |",
         "|---|---|---|"]
    for k in METRICS:
        rec = {r["window"]: r["n_pairings_top10"] for r in o["recurrence_over_12_pairings"][k]}
        L.append(f"| {k} | {UNITS[k]} | " + ", ".join(f"{t['window']} ({t['value']:.3f}; {rec.get(t['window'], 0)}/12)"
                                                      for t in o["primary_top10"][k]) + " |")
    L += ["", "## Agreement between equal-information metrics (Spearman over windows; descriptive)", "",
          "| Pair | Primary | Range over 12 pairings | Windows |", "|---|---|---|---|"]
    for k, v in o["rank_agreement"].items():
        L.append(f"| {k} | {v['primary']:.3f} | {v['min_over_12']:.3f}–{v['max_over_12']:.3f} | {v['n_windows']} |")
    L += ["", "## Cross-study comparison with KRAS-STRUCT-1 (frozen table; descriptive)", "",
          "Not attributable to inhibitor binding or G12D: these factors are not independently controlled across the "
          "studies, which also differ in crystal forms and conditions.", "",
          "| Metric | Spearman (shared windows) | Overlapping ten-largest window centres | STRUCT-1 ten largest | STRUCT-2 ten largest |",
          "|---|---|---|---|---|"]
    for k, v in o["cross_study_vs_struct1"].items():
        L.append(f"| {k} | {v['spearman']:.3f} ({v['n_shared_windows']}) | {v['overlap']} | "
                 f"{v['struct1_top10_centres']} | {v['struct2_top10_centres']} |")
    L += ["", "![KRAS-STRUCT-2 figure](KRAS-STRUCT-2_figure.png)", "",
          f"Resources: {o['resources']['cpu_s']:.1f} CPU-s, {o['resources']['wall_s']:.1f} s wall.", ""]
    return "\n".join(L)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("freeze", "run"):
        print(__doc__)
        sys.exit(2)
    sys.exit(freeze() if sys.argv[1] == "freeze" else run())
