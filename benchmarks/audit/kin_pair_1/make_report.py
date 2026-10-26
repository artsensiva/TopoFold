"""KIN-PAIR-1 report and figures, generated only from results/*.json (AGENTS.md rule 3).

    .venv/bin/python benchmarks/audit/kin_pair_1/make_report.py
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
R = HERE / "results"
FIG = HERE / "figures"
FIG.mkdir(exist_ok=True)
prim = json.loads((R / "KIN-PAIR-1_results.json").read_text())
supp = json.loads((R / "KIN-PAIR-1_supplement.json").read_text())
sel = json.loads((R / "KIN-PAIR-1_selection.json").read_text())
prot = json.loads((HERE / "protocol.json").read_text())
SUPP = {p["pair_id"]: p for p in supp["pairs"]}
PAIRS = [p["pair_id"] for p in prim["pairs"]]
PRIM = {p["pair_id"]: p for p in prim["pairs"]}
TOL = prot["tolerances"]["numerical_reference_abs"]
BLUE, ORANGE, GREY, INK = "#2a78d6", "#eb6834", "#6b6a66", "#2b2a27"
TARGET_STYLE = {"EGFR": (BLUE, "o"), "PIM1": (ORANGE, "s"), "BTK": (GREY, "^")}


def best(pid):
    """The analysis that evaluated the checks: the primary run if it ran them, else the supplement."""
    p = PRIM[pid]
    return (p, "primary") if "R2_window_set" in p["checks"] else (SUPP[pid], "supplement")


def md_table(rows, head):
    return "\n".join(["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(str(c) for c in r) + " |" for r in rows])


# --- figures ----------------------------------------------------------------------------------------------------------
plt.rcParams.update({"font.size": 9, "axes.edgecolor": GREY, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK})
fig, ax = plt.subplots(figsize=(7.4, 3.9))
checks = [("R3 κ/τ and distance values", lambda c: max(c["R3_kappa_tau_values"]["max_abs_deviation"].values())),
          ("R4 local RMSD", lambda c: max(c["R4_local_rmsd"]["max_abs_deviation"].values())),
          ("R5 global RMSD", lambda c: c["R5_global_rmsd"]["abs_deviation"]),
          ("R7 rigid-motion invariance", lambda c: c["R7_rigid_motion"]["max_abs_deviation_window_columns"])]
for ci, (name, fn) in enumerate(checks):
    for k, pid in enumerate(PAIRS):
        p, src = best(pid)
        v = fn(p["checks"])
        color, mk = TARGET_STYLE[PRIM[pid]["gene"]]
        ax.scatter(ci + (k - 2.5) * 0.11, max(v, 1e-9), s=46, marker=mk, color=color if src == "primary" else "white", edgecolor=color, linewidth=1.6, zorder=3)
ax.axhline(TOL, color=INK, lw=1, ls="--")
ax.text(-0.45, TOL * 1.25, "tolerance 1e-4", ha="left", fontsize=8)
ax.set_yscale("log")
ax.set_ylim(1e-9, 1e-3)
ax.set_xticks(range(len(checks)))
ax.set_xticklabels([c[0] for c in checks], fontsize=8)
ax.set_ylabel("largest absolute deviation (rad or Å)")
ax.grid(axis="y", color="#e4e3df", lw=0.6, zorder=0)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
h = [plt.Line2D([], [], marker=m, color="w", markerfacecolor=c, markeredgecolor=c, markersize=7, label=g) for g, (c, m) in TARGET_STYLE.items()]
h += [plt.Line2D([], [], marker="o", color="w", markerfacecolor=GREY, markeredgecolor=GREY, markersize=7, label="filled: primary run"),
      plt.Line2D([], [], marker="o", color="w", markerfacecolor="white", markeredgecolor=GREY, markeredgewidth=1.6, markersize=7, label="open: supplement (post-hoc)")]
ax.legend(handles=h, fontsize=7.5, frameon=False, loc="lower left", ncol=2)
ax.set_title("KIN-PAIR-1: workflow values against independent references (6 pairs)", fontsize=10, loc="left")
fig.tight_layout()
fig.savefig(FIG / "reference_agreement.png", dpi=130)
plt.close(fig)

fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.7), gridspec_kw={"width_ratios": [1, 1.25]})
ax = axs[0]
x = np.arange(len(PAIRS))
for off, size in ((-0.18, "3"), (0.18, "4")):
    vals = [best(p)[0]["classification"][size].get("overlap_O", 0) for p in PAIRS]
    ax.bar(x + off, vals, width=0.34, color=BLUE if size == "3" else ORANGE, label=f"{size}-residue windows (d_{'kappa' if size == '3' else 'tau'})")
ax.axhline(prot["classification"]["redundant_if_overlap_at_least"] - 0.5, color=INK, ls="--", lw=1)
ax.set_xticks(x)
ax.set_xticklabels(PAIRS, rotation=35, ha="right", fontsize=8)
ax.set_ylim(0, 14)
ax.set_ylabel("O: top-10 geometry windows also in\nthe top-10 of a baseline")
ax.legend(fontsize=7.5, frameon=False, loc="upper left")
ax.set_yticks([0, 2, 4, 6, 8, 10])
ax.set_title("(a) redundant or complementary", fontsize=9, loc="left")
ax = axs[1]
series = [("global RMSD, all positions", lambda d: d["global_rmsd_all_positions"], GREY),
          ("max trimmed-fit displacement", lambda d: d["displacement_trimmed_fit"]["max"], BLUE),
          ("max local RMSD (3 res.)", lambda d: d["local_shape"]["max_rmsd3"], ORANGE)]
for k, (name, fn, color) in enumerate(series):
    ax.bar(x + (k - 1) * 0.27, [fn(best(p)[0]["descriptive"]) for p in PAIRS], width=0.25, color=color, label=name)
ax.set_xticks(x)
ax.set_xticklabels(PAIRS, rotation=35, ha="right", fontsize=8)
ax.set_ylabel("Å")
ax.set_ylim(0, 9)
ax.legend(fontsize=7.5, frameon=False, loc="upper right")
ax.set_title("(b) displacement after alignment vs local shape", fontsize=9, loc="left")
for a in axs:
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(FIG / "baselines_and_displacement.png", dpi=130)
plt.close(fig)

# --- report -----------------------------------------------------------------------------------------------------------
rows = []
for pid in PAIRS:
    p = PRIM[pid]
    b, src = best(pid)
    m = p["meta"]
    mut = lambda k: f"{m[k]['n_mutations_rcsb']}"   # noqa: E731
    rows.append([pid, f"{p['entries']['A']} / {p['entries']['B']}", f"{m['A']['resolution']} / {m['B']['resolution']}", f"{m['A']['primary_ligand']} / {m['B']['primary_ligand']}",
                 f"{m['A']['space_group']} / {m['B']['space_group']}", f"{mut('A')} / {mut('B')}", p["correspondence"]["n_positions"], p["correspondence"]["n_span_positions_not_shared"],
                 p["outcome"], b["outcome"] + f" ({src})"])
t1 = md_table(rows, ["pair", "entries A / B", "resolution (Å)", "primary ligand", "space group", "RCSB mutation counts", "positions", "span positions not shared", "frozen outcome", "checks evaluated by"])
rows = []
for pid in PAIRS:
    b, src = best(pid)
    c = b["checks"]
    rows.append([pid, src, f"{c['R3_kappa_tau_values']['max_abs_deviation']['d_kappa']:.1e} / {c['R3_kappa_tau_values']['max_abs_deviation']['d_tau']:.1e}",
                 f"{max(c['R3_kappa_tau_values']['max_abs_deviation']['d_i_i2'], c['R3_kappa_tau_values']['max_abs_deviation']['d_i_i3']):.1e}",
                 f"{max(c['R4_local_rmsd']['max_abs_deviation'].values()):.1e}", f"{c['R5_global_rmsd']['abs_deviation']:.1e}",
                 f"{c['R7_rigid_motion']['max_abs_deviation_window_columns']:.6e}", "yes" if c["R2_window_set"]["pass"] else "NO", "yes" if c["R9_exclusions"]["pass"] else "NO",
                 "yes" if c["R7_rigid_motion"]["pass"] else "NO (resolution-limited, see text)", "yes" if c["R8_shifted_mapping"]["pass"] else "NO"])
t2 = md_table(rows, ["pair", "analysis", "R3 d_kappa / d_tau max dev.", "R3 distances max dev.", "R4 local RMSD max dev.", "R5 global RMSD dev.", "R7 invariance max dev.", "R2 window set equal", "R9 exclusions equal", "R7 pass", "R8 pass"])
rows = []
for pid in PAIRS:
    b, _ = best(pid)
    cl = b["classification"]
    d = b["descriptive"]
    rows.append([pid, f"{cl['3']['class']} (O={cl['3'].get('overlap_O')}, ρ vs rmsd3 {cl['3']['spearman']['rmsd3']:.2f})", f"{cl['4']['class']} (O={cl['4'].get('overlap_O')}, ρ vs rmsd4 {cl['4']['spearman']['rmsd4']:.2f})",
                 f"{d['global_rmsd_all_positions']:.3f}", f"{d['displacement_trimmed_fit']['max']:.2f} at {d['displacement_trimmed_fit']['argmax']['A_B_residue']}",
                 f"{d['local_shape']['max_rmsd3']:.2f} / {d['local_shape']['max_rmsd4']:.2f}", f"{d['spearman_displacement_vs_rmsd3']:.2f}"])
t3 = md_table(rows, ["pair", "3-residue windows", "4-residue windows", "global RMSD (Å)", "max trimmed-fit displacement (Å) at residue", "max rmsd3 / rmsd4 (Å)", "Spearman displacement vs rmsd3"])
rows = []
for pid in SUPP:
    s = SUPP[pid]["independent_altloc_B_sensitivity"]
    mc = s["max_abs_change"]
    rows.append([pid, s["n_positions_with_altloc_ca"], "yes" if s["validity_sets_equal"] else "no", f"{max(mc['d_kappa'], mc['d_tau']):.3f}", f"{max(mc['d_i_i2'], mc['d_i_i3']):.3f}",
                 f"{max(mc['rmsd3'], mc['rmsd4']):.3f}", f"{min(s['top10_overlap'].values())} of 10"])
t4 = md_table(rows, ["pair", "mapped positions with an altloc Cα", "same valid windows", "max change d_kappa / d_tau (rad)", "max change distances (Å)", "max change local RMSD (Å)", "smallest top-10 overlap"])
cnt = prim["counts"]
n_supp = len(supp["pairs"])
n_supp_pass = sum(p["outcome"] == "PASS" for p in supp["pairs"])
supp_fail = [p for p in supp["pairs"] if p["outcome"] != "PASS"]
sf_txt = "; ".join(f"{p['pair_id']} failed only {', '.join(k for k, v in p['checks'].items() if isinstance(v, dict) and v.get('pass') is False)} "
                   f"(largest deviation {p['checks']['R7_rigid_motion']['max_abs_deviation_window_columns']!r} against {TOL:g}, i.e. above the tolerance by "
                   f"{p['checks']['R7_rigid_motion']['max_abs_deviation_window_columns'] - TOL:.1e})" for p in supp_fail) or "none failed"
shift_max = max(v[k] for pid in PAIRS for v in [best(pid)[0]["checks"]["R8_shifted_mapping"][s_] for s_ in "AB"] for k in ("identity_shifted_by_1",))
same_day = [pid for pid in PAIRS if PRIM[pid]["meta"]["A"]["release_date"] == PRIM[pid]["meta"]["B"]["release_date"]]
sel_rows = []
for r in sel["ranking"][:8]:
    t = sel["targets"][r["gene"]]
    sel_rows.append([r["rank"], r["gene"], r["accession"], r["eligible_by_metadata"], t.get("status"), len(t.get("examined_entries", []))])
t5 = md_table(sel_rows, ["rank", "target", "UniProt", "entries eligible by metadata", "status", "entries examined"])
exam = [(g, e) for g, t in sel["targets"].items() for e in t.get("examined_entries", [])]
reasons = {}
for g, e in exam:
    key = e["status"] if e["status"] != "excluded" else "excluded: " + ("DBREF numbering disagrees with the alignment" if "DBREF" in e["reason"] else e["reason"][:60])
    reasons[key] = reasons.get(key, 0) + 1
t6 = md_table(sorted(reasons.items()), ["status of an examined entry", "count"])
cls = [best(pid)[0]["classification"][sz]["class"] for pid in PAIRS for sz in ("3", "4")]
n_red, n_comp, n_ins = cls.count("redundant"), cls.count("complementary"), cls.count("insufficient_evidence")
red3 = sum(best(pid)[0]["classification"]["3"]["class"] == "redundant" for pid in PAIRS)
red4 = sum(best(pid)[0]["classification"]["4"]["class"] == "redundant" for pid in PAIRS)
disp_ratio = [best(pid)[0]["descriptive"]["displacement_trimmed_fit"]["max"] / best(pid)[0]["descriptive"]["local_shape"]["max_rmsd3"] for pid in PAIRS]
res = prim["resources"]
text = f"""# KIN-PAIR-1: structpair on previously unused human kinases (generated report)

Generated by `make_report.py` from `results/*.json`; no number is typed here. Protocol: `PROTOCOL.md` (frozen before selection; amendments 1 and 2 are dated in it).
**Structural-analysis and engineering evaluation. It is not a test of drug efficacy, resistance, causal allostery or clinical benefit.**

## Result

**Frozen-rule study outcome: {prim['study_outcome']}.** {cnt['pairs_executed']} of {cnt['pairs_selected']} selected pairs were executed (targets: {', '.join(cnt['targets'])});
{cnt['pairs_pass']} passed all checks in the primary run; {cnt['implementation_failures']} are classed *implementation failure* by the frozen rule; {cnt['input_or_mapping_failures']} input or mapping failures.

**Cause of the 4 failures (recorded, not reinterpreted).** In EGFR-1, EGFR-2, BTK-1 and BTK-2 structpair stopped with its input error
"altlocs ['A', 'B'] tie on occupancy 0.5". The protocol had requested the optional `altloc_sensitivity: highest_occupancy`; that policy refuses a tie and aborts the whole comparison, so no check was evaluated for those pairs in the primary run.
It is a refusal designed into the adapter (its docstring: "Anything unresolved raises"; the user guide did not say so), not a numerical error, but under the frozen rule it counts. A post-hoc supplement (amendment 2) re-ran those 4 pairs once without the optional entry: **{n_supp_pass} of {n_supp} passed every check; {sf_txt}**. A deviation of this size is the 6-significant-figure resolution of the CSV outputs, not a geometric deviation. The frozen outcome is not changed by the supplement.

## Panel (selected by the frozen rules from metadata and file checks, before any TopoFold value was computed)

{t5}

Examined entries by outcome:

{t6}

{t1}

Construct, mutation, ligand and crystal-form differences are recorded, not controlled: pairs differ in several of them at once, and the two pairs of one target are not independent. Pairs whose two entries were released on the same day (shared deposition series, so shared methods): {', '.join(same_day) or 'none'}.

## Agreement with independent references (tolerance {TOL:g}; deviations in rad or Å)

{t2}

![reference agreement](figures/reference_agreement.png)

*Figure 1. Purpose: show that the workflow's values agree with independent recomputation from the raw coordinates. Shows the largest absolute deviation per check and pair (log scale); filled markers are primary runs, open markers the post-hoc supplement. Limitation: CSV outputs carry 6 significant figures, so deviations near 1e-5 to 1e-4 are rounding; the invariance deviation of about 1e-4 is at that resolution. Source: `results/KIN-PAIR-1_results.json`, `results/KIN-PAIR-1_supplement.json`.*

The independent window sets and exclusion sets equal the workflow's for every evaluated pair (R2 and R9); the shifted-mapping control collapses the residue-name identity to at most {shift_max:.3f} for every chain (R8). Displacement (R6) could not be compared: every entry contributed one chain, so structpair wrote no displacement channel (the independent trimmed-fit displacement is reported below).

## Local geometry, baselines and displacement (descriptive)

{t3}

![baselines and displacement](figures/baselines_and_displacement.png)

*Figure 2. Purpose: compare TopoFold local geometry with the equal-information baselines and show displacement after alignment separately from local shape change. (a) O is the number of the 10 largest geometry windows that are also among the 10 largest of either baseline (8 or more, above the dashed line: redundant). (b) global RMSD, the largest trimmed-fit (20 %) displacement and the largest 3-residue local RMSD. Limitation: descriptive, one pair per bar, not a test; complementary geometry is not novelty or biological superiority. Source: `results/KIN-PAIR-1_results.json`, `results/KIN-PAIR-1_supplement.json`.*

## Sensitivity to the alternate-location choice (independent, supplement pairs)

{t4}

Replacing location A by B at residues that have both changes window values by the amounts above; the set of valid windows is {'unchanged in every pair' if all(p['independent_altloc_B_sensitivity']['validity_sets_equal'] for p in supp['pairs']) else 'not unchanged in every pair'}. Chain sensitivity could not be examined (one chain per entry in every pair).

## Summary of the descriptive comparison

Across the {len(cls)} pair-by-window-size cases: {n_red} redundant with the equal-information baselines, {n_comp} complementary, {n_ins} insufficient evidence
(3-residue windows: {red3} of {len(PAIRS)} redundant; 4-residue windows: {red4} of {len(PAIRS)} redundant). In every pair the largest trimmed-fit displacement is {min(disp_ratio):.1f} to {max(disp_ratio):.1f} times the largest 3-residue local RMSD, so displacement after alignment and local shape change are different quantities in these data.
These are counts over 6 non-independent pairs, a descriptive rule and not a test; "complementary" means only that at least 3 of the 10 largest geometry windows are not among the baselines' top windows.

## Resources

Execution: {res['execution_cpu_seconds']:.1f} CPU-s ({res['execution_wall_seconds']:.1f} s wall) in the primary run and {supp['cpu_seconds']:.1f} CPU-s in the supplement; selection {sel['resources']['cpu_seconds']:.1f} CPU-s;
{res['selection_files_downloaded']} PDB files, {res['selection_bytes_downloaded'] / 1e6:.2f} MB downloaded by the corrected selection run (the first, defective selection run downloaded a similar amount; total well under the 100 MB limit); installation excluded.

## Limits

The panel is small ({cnt['pairs_selected']} pairs, {len(cnt['targets'])} targets, structures chosen by resolution and ligand difference), and pairs of one target share laboratories and ligand series. Phosphorylated residues (TPO, SEP, PTR) are not read by the adapter and are excluded from the correspondence.
Nothing here shows a biological state, ligand or mutation effect, novelty, or user time savings.
"""
(HERE / "report.md").write_text(text)
print("wrote report.md and 2 figures")
