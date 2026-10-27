"""Endpoint cohort table (correction 2026-10-27): separates assay types per record and cross-checks the GitHub sheet against the
official Supplementary Information PDF and the stated ranges of the Schrock 2013 abstract. No model is trained or scored.

    .venv/bin/python -I benchmarks/audit/abl_response_intake_2026-10-27/build_cohort.py
Needs pdftotext (poppler) for the SI PDF; raises if the raw files are missing (no substitute is synthesised).
"""
import collections
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from xlsx_min import read_xlsx  # noqa: E402

RAW = HERE / "data_raw"
SI_PDF = RAW / "corrections_2026-10-27/hauser2018_SI_MOESM1.pdf"
SD3 = RAW / "hauser/Supplementary Data 3.xlsx"
SD2 = RAW / "hauser/Supplementary Data 2.xlsx"
for p in (SI_PDF, SD3, SD2, RAW / "P00519.fasta"):
    if not p.is_file():
        raise FileNotFoundError(p)
RT = 0.0019872041 * 300.0
TEMPLATE = {"axitinib": "4WA9", "bosutinib": "3UE4", "dasatinib": "4XEY", "imatinib": "1OPJ", "nilotinib": "3CS9", "ponatinib": "3OXZ", "erlotinib": "docked into the 3UE4 receptor", "gefitinib": "docked into the 3UE4 receptor"}
SCHROCK = ("bosutinib", "dasatinib", "imatinib", "nilotinib", "ponatinib")
# Literature-stated values (Schrock et al., Blood 122, 3992 (2013), abstract text retrieved from Crossref metadata): assay and IC50 ranges in nM.
ABSTRACT = {"assay": "engineered Ba/F3 cells, cell viability at 72 hours", "n_single_mutants": 21,
            "ic50_range_nM_excluding_T315I": {"ponatinib": (3, 16), "imatinib": (201, 10000), "nilotinib": (12, 784), "dasatinib": (2, 104), "bosutinib": (40, 1280)},
            "T315I_statement": "IC50s for the other TKIs excluding ponatinib: >4000 for all (T315I)"}

si = subprocess.run(["pdftotext", "-layout", str(SI_PDF), "-"], capture_output=True, text=True, check=True).stdout.splitlines()
# --- SI Tables 2-7: per inhibitor, wild-type IC50 and per-mutation IC50 and experimental ΔΔG
heads = [(i, re.match(r"\s*Supplementary Table (\d+)\. (\w+):", l)) for i, l in enumerate(si)]
heads = [(i, m) for i, m in heads if m and m.group(2).lower() in TEMPLATE]
si_tab = {}
for k, (i, m) in enumerate(heads):
    end = heads[k + 1][0] if k + 1 < len(heads) else len(si)
    # the next 'Supplementary Table' heading of any kind ends the block
    for j in range(i + 1, end):
        if re.match(r"\s*Supplementary Table \d+\.", si[j]):
            end = j
            break
    inh, rows, wt = m.group(2).lower(), {}, None
    for l in si[i:end]:
        w = re.match(r"\s*wild-type\s+(\S+)", l)
        if w:
            wt = float(w.group(1))
        r = re.match(r"\s*([A-Z]\d{3}[A-Z])\s+(\S+)\s+(-?\d+\.\d+)\s", l)
        if r:
            rows[r.group(1)] = (float(r.group(2)), float(r.group(3)))
    si_tab[inh] = {"table": int(m.group(1)), "wt_ic50_nM": wt, "rows": rows}
# --- SI Table 1: ΔΔG per mutation and inhibitor (columns: axitinib bosutinib dasatinib imatinib nilotinib ponatinib gefitinib erlotinib)
order = ["axitinib", "bosutinib", "dasatinib", "imatinib", "nilotinib", "ponatinib", "gefitinib", "erlotinib"]
t1 = {}
start = next(i for i, l in enumerate(si) if re.match(r"\s*Supplementary Table 1\.", l))
stop = next(i for i, l in enumerate(si) if re.match(r"\s*Supplementary Table 2\.", l))
for l in si[start:stop]:
    tok = l.split()
    if tok and re.fullmatch(r"[A-Z]\d{3}[A-Z]", tok[0]) and len(tok) >= 9:
        t1[tok[0]] = {inh: (None if v == "nd" else float(v)) for inh, v in zip(order, tok[1:9])}
# --- GitHub sheet
rows = read_xlsx(SD3)["Sheet1"]
full = [r for r in rows[1:] if len(r) == 11]
sd2 = read_xlsx(SD2)["Davis et al. data"]
davis_wt = {}
for r in sd2:
    if len(r) > 3 and r[3] == "ABL1-nonphosphorylated":
        davis_wt = {"erlotinib": {"non_phosphorylated": float(r[9]), "phosphorylated": float(r[17])}, "gefitinib": {"non_phosphorylated": float(r[10]), "phosphorylated": float(r[18])}}
uni = "".join(l.strip() for l in open(RAW / "P00519.fasta") if not l.startswith(">"))


def implied_wt(inh):
    rr = [r for r in full if r[0] == inh]
    return float(np.median([r[2] * np.exp(-r[3] / RT) for r in rr]))


records = []
for r in full:
    inh, mut, ic50, ddg = r[0], r[1], r[2], r[3]
    pos = int(mut[1:-1])
    sirow = si_tab.get(inh, {}).get("rows", {}).get(mut)
    wt_si = si_tab.get(inh, {}).get("wt_ic50_nM")
    if inh in SCHROCK:
        cohort, study = "primary_candidate_bafs_viability", "Schrock et al. 2013 (Blood abstract; underlying figure not obtained)"
        assay = ABSTRACT["assay"] + " (from the abstract)"
        qty, wt_ref, wt_basis = "IC50 (nM)", wt_si, "stated in Hauser SI Table %d; original WT measurement not independently seen" % si_tab[inh]["table"]
    elif inh == "axitinib":
        cohort, study = "secondary_axitinib_cell_ic50", "Pemovska et al. 2015 (Nature; paywalled, source table not obtained)"
        assay, qty, wt_ref, wt_basis = "unresolved (not verified against the source)", "IC50 (nM), assay type unverified", wt_si, "stated in Hauser SI Table %d; original WT measurement not seen" % si_tab[inh]["table"]
    else:
        cohort, study = "secondary_binding_davis_kd", "Davis et al. 2011 (Nat. Biotechnol.; paywalled, supplementary table not obtained)"
        assay, qty = "competition binding (kinase panel), phosphorylated Abl block per this check; construct not verified", "binding constant; headed 'IC50 (nM)' in Hauser sheets"
        wt_ref = round(implied_wt(inh), 3)
        wt_basis = "inferred from the sheet; equals the phosphorylated-Abl wild-type value %.0f in the Hauser reproduction of the Davis data (not the Davis paper)" % davis_wt[inh]["phosphorylated"]
    recomputed = RT * np.log(ic50 / wt_ref) if wt_ref else None
    flags = []
    if inh in SCHROCK:
        lo, hi = ABSTRACT["ic50_range_nM_excluding_T315I"][inh]
        if mut != "T315I" and not (lo <= ic50 <= hi):
            flags.append("ic50_outside_range_stated_in_source_abstract")
        if mut == "T315I" and inh != "ponatinib":
            flags.append("abstract_reports_only_'>4000'; numeric value comes from the compilation")
    if ic50 >= 10000:
        flags.append("at_10000_nM_assay_ceiling_stored_as_equality")
    if recomputed is not None and abs(recomputed - ddg) > 0.05:
        flags.append("ddg_in_sheet_disagrees_with_ic50_and_wt")
    sit1 = t1.get(mut, {}).get(inh)
    if sit1 is None or abs(sit1 - ddg) > 0.006:
        flags.append("ddg_differs_between_github_sheet_and_SI_table1")
    if sirow and (sirow[0] != ic50 or abs(sirow[1] - ddg) > 0.006):
        flags.append("ic50_or_ddg_differs_between_github_sheet_and_SI_table")
    excl = None
    if "ic50_outside_range_stated_in_source_abstract" in flags and "ddg_in_sheet_disagrees_with_ic50_and_wt" in flags:
        excl = "unresolved discrepancy: IC50 outside the source abstract's stated range and inconsistent with the ΔΔG (also in the official SI)"
    records.append({"id": f"{inh}:{mut}", "inhibitor": inh, "mutation": mut, "position": pos, "wt_aa": mut[0], "mut_aa": mut[-1], "wt_matches_P00519": uni[pos - 1] == mut[0],
                    "template": TEMPLATE[inh], "cohort": cohort, "source_study": study, "assay_construct": assay, "quantity": qty, "units": "nM", "value": ic50,
                    "wt_reference_nM": wt_ref, "wt_reference_basis": wt_basis, "ddg_sheet_kcal_mol": ddg, "ddg_SI_table1_kcal_mol": sit1, "ddg_recomputed_kcal_mol": None if recomputed is None else round(float(recomputed), 3),
                    "flags": flags, "exclusion": excl, "resistant_ge_1.36": ddg >= 1.36})
prim = [r for r in records if r["cohort"] == "primary_candidate_bafs_viability"]
prim_kept = [r for r in prim if r["exclusion"] is None]
flagged = [r for r in prim_kept if any(f.startswith(("abstract_reports", "at_10000")) for f in r["flags"])]
sens = [r for r in prim_kept if r not in flagged]


def summ(rs):
    pos = collections.Counter(r["position"] for r in rs)
    return {"records": len(rs), "substitutions": len({r["mutation"] for r in rs}), "positions": len(pos), "drugs": len({r["inhibitor"] for r in rs}), "templates": len({r["template"] for r in rs}),
            "source_groups": len({r["source_study"] for r in rs}), "records_per_position": dict(sorted(pos.items())), "resistant_records": sum(r["resistant_ge_1.36"] for r in rs),
            "positions_with_resistant_record": sorted({r["position"] for r in rs if r["resistant_ge_1.36"]})}


# ranges of the abstract against the sheet; wild-type values against the SI
ab_chk = {}
for inh in SCHROCK:
    v = [r["value"] for r in prim if r["inhibitor"] == inh and r["mutation"] != "T315I"]
    lo, hi = ABSTRACT["ic50_range_nM_excluding_T315I"][inh]
    ab_chk[inh] = {"sheet_min_max_excluding_T315I": [min(v), max(v)], "abstract_range": [lo, hi], "wt_in_SI_nM": si_tab[inh]["wt_ic50_nM"], "implied_wt_from_sheet_nM": round(implied_wt(inh), 2)}
si_vs_sheet = {"six_inhibitors_ic50_equal": sum(1 for r in records if r["inhibitor"] in si_tab and r["mutation"] in si_tab[r["inhibitor"]]["rows"] and si_tab[r["inhibitor"]]["rows"][r["mutation"]][0] == r["value"]),
               "six_inhibitors_records": sum(1 for r in records if r["inhibitor"] in si_tab),
               "table1_ddg_equal_within_0.006": sum(1 for r in records if t1.get(r["mutation"], {}).get(r["inhibitor"]) is not None and abs(t1[r["mutation"]][r["inhibitor"]] - r["ddg_sheet_kcal_mol"]) <= 0.006), "all_records": len(records)}
out = {"schema": "topofold.abl-endpoint-cohort/1", "date": "2026-10-27", "si_wt_ic50_nM": {k: v["wt_ic50_nM"] for k, v in si_tab.items()}, "inputs": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (SI_PDF, SD3, SD2, RAW / "P00519.fasta")},
       "schrock_abstract_statements": ABSTRACT, "abstract_vs_sheet": ab_chk, "si_vs_github_sheet": si_vs_sheet, "davis_wt_in_hauser_reproduction_nM": davis_wt,
       "cohorts": {c: summ([r for r in records if r["cohort"] == c]) for c in sorted({r["cohort"] for r in records})},
       "primary_candidate_after_exclusion": summ(prim_kept), "primary_flagged_censored_or_bound_records": [r["id"] for r in flagged], "primary_sensitivity_subset_without_flagged": summ(sens),
       "excluded": [{"id": r["id"], "reason": r["exclusion"]} for r in records if r["exclusion"]],
       "records_with_unresolved_flags": {r["id"]: r["flags"] for r in records if r["flags"]}, "records": records}
(HERE / "results/endpoint_cohort_2026-10-27.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: out[k] for k in ("si_vs_github_sheet", "abstract_vs_sheet", "excluded", "primary_flagged_censored_or_bound_records")}, indent=0)[:2500])
for k in ("cohorts", "primary_candidate_after_exclusion", "primary_sensitivity_subset_without_flagged"):
    print(k, json.dumps(out[k])[:1500])
