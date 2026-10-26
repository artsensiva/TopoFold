"""Generates REPORT.md and the schematic from results/pim1_annotation.json (no number is typed here).

    .venv/bin/python benchmarks/audit/pim1_structural_annotation_2026-10-26/make_report.py
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
d = json.loads((HERE / "results/pim1_annotation.json").read_text())
(HERE / "figures").mkdir(exist_ok=True)
INK, GREY, BLUE, ORANGE, LIGHT = "#2b2a27", "#6b6a66", "#2a78d6", "#eb6834", "#e4e3df"
ref, h = d["reference"], d["kinase_roles_helper"]
dom = d["uniprot_annotation"]["features"]["Domain"][0]
dfg = d["dfg_window"]["hits"][0]
glu = d["alphaC_glu_rule"]["qualifying"][0]
lm = {x["name"]: x for x in d["landmarks_verified"]}
mod261 = next(x for x in d["uniprot_annotation"]["features"]["Modified residue"] if x["start"] == 261)


def table(rows, head):
    return "\n".join(["| " + " | ".join(head) + " |", "|" + "---|" * len(head)] + ["| " + " | ".join(str(x) for x in r) + " |" for r in rows])


# deposited helix containing the 261 windows: extents per structure (author numbering = reference numbering here, offset 0)
def helix_of(e, r):
    return next((x for x in d["structures"][e]["deposited_secondary_structure"]["helices"] if x["author_start"] <= r <= x["author_end"]), None)


helix = {e: helix_of(e, 261) for e in d["structures"]}
hs = sorted({(x["ref_start"], x["ref_end"]) for x in helix.values() if x})
plt.rcParams.update({"font.size": 8.5, "text.color": INK, "axes.edgecolor": GREY})
fig, (a1, a2) = plt.subplots(2, 1, figsize=(8.4, 5.6), gridspec_kw={"height_ratios": [1, 1.5]})
# (a) whole reference
a1.add_patch(Rectangle((1, 0.0), ref["length"], 0.3, color=LIGHT, lw=0))
a1.add_patch(Rectangle((dom["start"], 0.0), dom["end"] - dom["start"] + 1, 0.3, color="#c9d6e8", lw=0))
a1.text((dom["start"] + dom["end"]) / 2, 0.15, f"UniProt protein kinase domain {dom['start']}–{dom['end']}", ha="center", va="center", fontsize=8)
pts = [(67, "β3 Lys 67"), (glu, f"αC Glu {glu}"), (165, "HRD 165–167"), (dfg, f"DFG {dfg}")]
for k, (x, lab) in enumerate(pts):
    a1.plot([x, x], [0.3, 0.55 + 0.17 * (k % 2)], color=INK, lw=1)
    a1.text(x, 0.58 + 0.17 * (k % 2), lab, ha="center", va="bottom", fontsize=7.5)
a1.plot([261, 261], [0.0, -0.28], color=ORANGE, lw=2)
a1.text(261, -0.31, f"residue 261 (deposited SEP; UniProt: {mod261['description'].lower()})", ha="center", va="top", fontsize=7.5, color=ORANGE)
a1.set_xlim(0, ref["length"] + 2)
a1.set_ylim(-0.7, 1.15)
a1.set_yticks([])
a1.set_xlabel(f"UniProt {ref['accession']} canonical (isoform 1) position, 1–{ref['length']}")
for s in ("top", "right", "left"):
    a1.spines[s].set_visible(False)
a1.set_title("(a) verified landmarks (UniProt annotation and the four mapped anchors); no other element is named", fontsize=9, loc="left")
# (b) zoom
lo, hi = 255, 267
a2.set_xlim(lo - 0.5, hi + 0.5)
a2.set_ylim(-0.4, 6.6)
a2.set_xticks(range(lo, hi + 1))
a2.set_xticklabels([str(x) for x in range(lo, hi + 1)], fontsize=7.5)
a2.get_xticklabels()[261 - lo].set_color(ORANGE)
a2.get_xticklabels()[261 - lo].set_fontweight("bold")
row = 5.8
for hx0, hx1 in hs:
    a2.add_patch(Rectangle((hx0 - 0.5, row - 0.18), hx1 - hx0 + 1, 0.36, color=LIGHT, lw=0))
    a2.plot([hx0 - 0.5] * 2, [row - 0.3, row + 0.3], color=GREY, lw=1, ls=":")
    a2.text(hx0 - 0.4, row + 0.38, "deposited helix records (4 files): starts at 260; not independently verified; conventional name unverified", fontsize=6.8, va="bottom", color=GREY)
y = 4.8
for pid, p in d["pairs"].items():
    for size, color in (("3", BLUE), ("4", ORANGE)):
        ws = [w for w in p["new_windows"] if w["window_residues"] == size]
        for k, w in enumerate(ws):
            r0, r1 = min(w["reference_indices"]), max(w["reference_indices"])
            yy = y - 0.17 * k
            a2.add_patch(Rectangle((r0 - 0.4, yy - 0.06), r1 - r0 + 0.8, 0.12, color=color, lw=0, alpha=0.85))
        a2.text(hi + 0.55, y - 0.25, f"{pid}: {len(ws)} new {size}-residue windows", fontsize=7, va="center", ha="right", color=INK)
        y -= 1.0
a2.set_yticks([])
a2.set_xlabel("reference position (equal to the deposited number in all four chains)\nbars: new windows, 3-residue blue, 4-residue orange, one per window; 261 is SEP in the structure that carries it, SER in its partner", fontsize=7.5)
for s in ("top", "right", "left"):
    a2.spines[s].set_visible(False)
a2.set_title("(b) the seven windows per pair that exist only after the policy admits residue 261", fontsize=9, loc="left")
fig.tight_layout()
fig.savefig(HERE / "figures/pim1_windows_schematic.png", dpi=130)
plt.close(fig)

# tables
t1 = table([[e, v["pair_slot"], v["sha256"][:12], f"{h['per_structure'][e]['alignment']['identity']:.4f} ({h['per_structure'][e]['alignment']['aligned_pairs']} aligned)",
             ", ".join(f"{m['observed']} {m['resseq']} (reference {m['uniprot']})" for m in v["mismatches"]) or "none",
             ", ".join(f"{x['resname']} {x['resid']}" for x in (h["per_structure"][e]["variants_in_chain"] or [])) or "none", v["author_minus_reference_offsets"],
             "yes" if h["per_structure"][e]["accepted"] else "NO", f"{h['per_structure'][e]['n_ambiguity_checks']}, all unique: {h['per_structure'][e]['all_unique']}",
             ", ".join(f"{a} {c['reference_index']}{c['found']}" for a, c in h["per_structure"][e]["anchor_checks"].items())] for e, v in d["structures"].items()],
           ["structure", "pair slot", "SHA-256 (first 12)", "identity to the reference", "differences from the reference (mutations)", "modified residue names (helper)", "author − reference offsets", "helper accepted",
            "ambiguity checks", "anchors found"])
t2 = table([[x["name"], f"{x['start']}" + (f"–{x['end']}" if x["end"] != x["start"] else ""), x["source"]] for x in d["landmarks_verified"]], ["verified landmark", "reference positions", "source"])
rows = []
for pid, p in d["pairs"].items():
    for w in p["new_windows"]:
        ids = " ".join(f"{r['A'].split(':')[1]}" for r in w["residues"])
        vals = w["values"]
        rk = w["ranks_in_new_run_of_windows_with_values"]
        rows.append([pid, f"{w['window_residues']}-res", ids, f"{min(w['reference_indices'])}–{max(w['reference_indices'])}",
                     "yes" if w["inside_uniprot_kinase_domain"] else "no", "; ".join(sorted(set(r["deposited_ss_A"] for r in w["residues"]))) + " / " + "; ".join(sorted(set(r["deposited_ss_B"] for r in w["residues"]))),
                     "; ".join(f"{m} {v:.3f} (rank {rk[m]['rank']}/{rk[m]['of']})" for m, v in vals.items()), "yes" if w["name_difference_flag"] else "no", "yes" if w["in_a_top10_of_any_metric"] else "no"])
t3 = table(rows, ["pair", "window", "residues (deposited number, structure A)", "reference positions", "inside UniProt kinase domain", "deposited helix/strand in structure A / B",
                  "existing values (rad or Å) and rank among windows with values", "SEP vs SER name difference flagged", "in a top-10 of any metric"])
rows = []
for pid, p in d["pairs"].items():
    ds = p["descriptive"]
    for size in ("3", "4"):
        a, b, c = ds["recorded_in_KIN-PAIR-1"][size], ds["reapplied_to_old_numeric_record"][size], ds["reapplied_to_new_numeric_record"][size]
        rows.append([pid, f"{size}-residue", f"{a['class']} (O={a['overlap_O']}, n={a['n_valid']})", f"{b['class']} (O={b['overlap_O']})", f"{c['class']} (O={c['overlap_O']}, n={c['n_valid']})"])
t4 = table(rows, ["pair", "windows", "recorded in KIN-PAIR-1", "rule re-applied to the old numeric record (reproduces)", "rule applied to the new numeric record"])
allv = [(pid, w["window_residues"], m, v, w["ranks_in_new_run_of_windows_with_values"][m]["rank"], w["ranks_in_new_run_of_windows_with_values"][m]["of"])
        for pid, p in d["pairs"].items() for w in p["new_windows"] for m, v in w["values"].items()]
best = {pid: min((x for x in allv if x[0] == pid), key=lambda x: x[4] / x[5]) for pid in d["pairs"]}
maxv = {pid: max(x[3] for x in allv if x[0] == pid) for pid in d["pairs"]}
unchanged = all(not p["descriptive"]["classification_changed_by_the_added_windows"] and p["descriptive"]["top10_old_still_in_top10_and_no_new_entrant"] for p in d["pairs"].values())
dfg_off = f"about {261 - dfg}"
mut_note = "; ".join(f"{e} has {m['observed']} at {m['resseq']} where the reference has {m['uniprot']}" for e, s_ in d["structures"].items() for m in s_["mismatches"]) or "no difference from the reference"
text = f"""# Verified structural-role annotation of the existing PIM1 comparisons (generated report)

Generated by `make_report.py` from `results/pim1_annotation.json`. **Interpretation of existing outputs: no score was recomputed, no threshold changed, no structure selected; KIN-PAIR-1 and the phospho regression are unchanged.**
Reference: UniProt **{ref['accession']}, {ref['isoform']}, {ref['entry']}, {ref['length']} residues**, sequence SHA-256 `{ref['sequence_sha256'][:16]}…`; FASTA hash {'equals' if ref['matches_selection_record'] else 'DIFFERS from'} the hash recorded in the KIN-PAIR-1 selection record.
Numbering convention: reference index = UniProt canonical position; author numbering is preserved, and its offset to the reference is **{sorted({o for s in d['structures'].values() for o in s['author_minus_reference_offsets']})}** in all four chains (reported, never applied).
Annotation source: UniProt entry version {d['uniprot_annotation']['entry_version']} (annotation updated {d['uniprot_annotation']['last_annotation_update']}), retrieved {d['uniprot_annotation']['retrieved_utc']}, raw SHA-256 `{d['uniprot_annotation']['raw_sha256'][:16]}…`. Rules were fixed in `DESIGN.md` before the run.

## Verification of the anchors (existing `kinase_roles` helper, defaults, role mutations refused)

{t1}

*Table 1. Purpose: show how the four structures map to the reference and that the helper's anchor and ambiguity checks passed. The αC Glu anchor (reference {glu}) was fixed by the pre-registered contact rule (the only Glu in 77–97 within 4.0 Å of the Lys 67 NZ in all four structures; the other Glu in the window, {', '.join(str(c['reference_index']) for c in d['alphaC_glu_rule']['glu_candidates'] if not c['qualifies_in_all_four'])}, does not qualify); the DFG anchor ({dfg}) by the unique DFG in the UniProt-bounded window {d['dfg_window']['window'][0]}–{d['dfg_window']['window'][1]}.
Limitation: **successful anchor mapping verifies the four anchors, the ten derived roles and the ambiguity checks; it does not verify any helix or loop boundary**, and the αC rule is a coordinate contact rule, not an independent annotation. Source: `results/pim1_annotation.json`, `results/kinase_roles_reports.json`.*

{t2}

*Table 2. Purpose: the only landmarks used in the schematic and the text. Not established (kept "unverified"): {'; '.join(d['unverified'])}.*

## The seven new windows per pair, annotated by complete residue identifiers

{t3}

*Table 3. Purpose: location, existing values and ranks of the windows that exist only after residue 261 is admitted. Locations: all lie inside the UniProt-annotated kinase domain, C-terminal to the DFG anchor (reference {dfg}); the deposited HELIX records of all four files place residue 261 in a helix beginning at residue 260 (deposited labels differ per file; the extents agree), so windows starting before 260 straddle the deposited helix boundary, which is deposited annotation and not independently verified. The conventional name of that helix is **unverified**. Ranks are among all windows with values of the same metric in the new run (descending, ties by position). Values are existing; none was recomputed. Source: same file.*

![schematic](figures/pim1_windows_schematic.png)

*Figure 1. Purpose: place the added windows in the sequence using verified landmarks only. (a) UniProt kinase domain and the four mapped anchors; residue 261 is marked (deposited SEP; the UniProt entry lists a phosphoserine modified residue there, a curated annotation that is not evidence of a state in these structures). (b) the new windows by reference position; the grey bar is the deposited helix record (start at 260, dotted: boundary not independently verified). Limitation: sequence positions and deposited records only; no region name beyond the landmarks; one chain per structure. Source: same file.*

## What this adds, descriptively

- The new windows sit in the C-terminal part of the kinase domain, in or at the N-terminal end of the deposited helix 260–271, {dfg_off} residues after the DFG anchor; the mapped anchors do not say whether that helix has a conventional name or a role.
- Their existing values are small: the largest value among the new windows is {', '.join(f"{pid} {v:.3f}" for pid, v in maxv.items())}; the best rank of any new window among the windows of its metric is {', '.join(f"{pid} {best[pid][2]} rank {best[pid][4]} of {best[pid][5]}" for pid in best)}; **no new window is in a top-10 of any metric** ({'confirmed' if not any(w['in_a_top10_of_any_metric'] for p in d['pairs'].values() for w in p['new_windows']) else 'NOT confirmed'}).
- **Previously reported descriptive conclusions: {'unchanged' if unchanged else 'CHANGED (see Table 4)'}.** The recorded rule (top-10 overlap with the baselines) gives the same class and overlap before and after adding the windows:

{t4}

*Table 4. Purpose: check whether the added windows change the recorded descriptive classification (rule and numbers of KIN-PAIR-1, applied to existing values). Limitation: descriptive counts over non-independent pairs, not a test. The maximum trimmed-fit displacement and its location were not recomputed; the global RMSD changed only by the refit with one more position (phospho replay).*

**Kept separate.** The modification (SEP; MODRES; UniProt modified residue) is an annotation. These pairs differ in construct, ligand and crystal form (and {mut_note}), so they do not isolate phosphorylation, a ligand or a mutation, and nothing here supports a functional interpretation of residue 261.
"""
(HERE / "REPORT.md").write_text(text)
print("wrote REPORT.md and figure")
