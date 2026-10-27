"""Measured counts and structural-input checks for the MET kinase inhibitor DMS candidate (Estevam et al., eLife 2025).
Read-only; no model is trained or scored. Run: python3 -I screen_met.py   (writes results/met_screen.json)
"""
import collections
import csv
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
RAW = HERE / "data_raw"
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}
uni = "".join(l.strip() for l in open(RAW / "P08581.fasta") if not l.startswith(">"))
rows = list(csv.DictReader(open(RAW / "met_dms/WT_rosace_effect_all.tsv"), delimiter="\t"))
OFFSET = 1058  # DMS position p is UniProt P08581 residue p + 1058 (verified below)
seq = {}
for r in rows:
    seq[int(r["position"])] = r["wildtype"]
mapping_ok = all(uni[p + OFFSET - 1] == a for p, a in seq.items())
MISS = [r for r in rows if r["type"] == "missense"]


def num(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


by = collections.defaultdict(dict)
for r in MISS:
    v = num(r["ROSACE_effects"])
    by[r["inhibitor"]][(int(r["position"]), r["mutation"])] = {"effect": v, "gof": r["ROSACE.GOF"] == "TRUE", "lof": r["ROSACE.LOF"] == "TRUE", "alt": r["ROSACE_alt"] == "TRUE"}
dmso = by["DMSO"]
per = {}
for inh, d in sorted(by.items()):
    ok = {k: v for k, v in d.items() if v["effect"] is not None}
    paired = [k for k in ok if k in dmso and dmso[k]["effect"] is not None]
    per[inh] = {"missense_variants_with_effect": len(ok), "positions": len({k[0] for k in ok}), "substitution_types": len({k[1] for k in ok}), "paired_with_DMSO": len(paired),
                "flag_GOF": sum(v["gof"] for v in ok.values()), "flag_LOF": sum(v["lof"] for v in ok.values()), "flag_alt": sum(v["alt"] for v in ok.values()),
                "DMSO_LOF_among_paired": sum(dmso[k]["lof"] for k in paired), "paired_and_not_DMSO_LOF": sum(1 for k in paired if not dmso[k]["lof"])}
n_pos = len({int(r["position"]) for r in rows})
all_conditions = set.intersection(*[{k for k, v in d.items() if v["effect"] is not None} for d in by.values()])
possible = n_pos * 19
# structures
PDB = {"2WGJ": "crizotinib", "6SDE": "savolitinib", "4EEV": "merestinib", "5T3Q": "AMG-458", "3RHK": "tivantinib", "4R1V": "tepotinib", "3QTI": "NVP-BVU972", "3DKC": "ATP-bound (no inhibitor)"}
st = {}
for p in PDB:
    f = RAW / "pdb" / f"{p}.pdb"
    txt = f.read_text().splitlines()
    res = next((l.split(":")[1].strip() for l in txt if l.startswith("REMARK   2 RESOLUTION")), None) if False else None
    resolution = next((l[23:30].strip() for l in txt if l.startswith("REMARK   2 RESOLUTION")), None)
    method = next((l[10:].strip() for l in txt if l.startswith("EXPDTA")), None)
    chains, hets = collections.defaultdict(dict), collections.Counter()
    for l in txt:
        if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[16] in " A":
            chains[l[21]][int(l[22:26])] = THREE.get(l[17:20].strip(), "X")
        elif l.startswith("HETATM") and l[17:20].strip() not in ("HOH",):
            hets[l[17:20].strip()] += 1
    ch = {}
    for c, d in chains.items():
        inspan = {n: a for n, a in d.items() if 1059 <= n <= 1345}
        ch[c] = {"resolved_CA_in_scanned_span": len(inspan), "wt_identity_matches_P08581_at_pdb_number": sum(1 for n, a in inspan.items() if uni[n - 1] == a), "first_last_in_span": [min(inspan), max(inspan)] if inspan else None}
    mism = {c: sorted((n, a, uni[n - 1]) for n, a in d.items() if 1059 <= n <= 1345 and uni[n - 1] != a) for c, d in chains.items()}
    st[p] = {"residues_differing_from_P08581_in_span": {c: v for c, v in mism.items() if v}, "ligand_in_title": PDB[p], "method": method, "resolution_A": resolution, "chains": ch, "hetero_groups": dict(hets), "sha256": hashlib.sha256(f.read_bytes()).hexdigest()}
docked = ["capmatinib", "glumetinib", "cabozantinib", "glesatinib analog"]
import statistics
filt = list(csv.DictReader(open(RAW / "met_dms/met_scores_filtered.tsv"), delimiter="\t"))
sds = [float(r["sd"]) for r in filt if r["type"] == "missense" and r["key"] != "DMSO"]
raw_bytes = sum(p.stat().st_size for p in RAW.rglob("*") if p.is_file())
work_bytes = sum(p.stat().st_size for p in (HERE / "work").glob("*") if p.is_file())
out = {"schema": "topofold.dms-candidate-met/1", "files": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [RAW / "met_dms/WT_rosace_effect_all.tsv", RAW / "met_dms/LICENSE", RAW / "P08581.fasta"]},
       "mapping": {"dms_positions": n_pos, "uniprot_residue_range": [1 + OFFSET, n_pos + OFFSET], "wild_type_identity_all_positions": mapping_ok},
       "conditions": sorted(by), "possible_missense_per_condition": possible, "missense_variants_measured_in_all_12_conditions": len(all_conditions), "per_condition": per,
       "synonymous_rows": sum(1 for r in rows if r["type"] == "synonymous"), "rows_total": len(rows), "rows_with_non_numeric_effect": sum(1 for r in rows if num(r["ROSACE_effects"]) is None),
       "inhibitor_dosing": "each inhibitor at its experimentally determined IC50 (article text): effects of different drugs are not on a common exposure scale",
       "posterior_sd_of_filtered_scores": {"file": "met_scores_filtered.tsv", "rows": len(filt), "median_sd_inhibitor_missense": statistics.median(sds), "max_sd_inhibitor_missense": max(sds), "note": "the filtered file keeps fewer variants than the Rosace effect file; replicate-level values are not in the repository files read"},
       "download_ledger": {"raw_files_bytes": raw_bytes, "publisher_and_policy_pages_bytes": work_bytes, "total_bytes_upper_bound": raw_bytes + work_bytes, "ceiling_bytes": 20_000_000, "not_downloaded": ["SRA reads PRJNA1136906", "GitHub Machine_Learning folder (58.6 MB)", "EX_rosace_effect_all.tsv and the unfiltered score tables"]},
       "structures": st, "docked_without_experimental_structure_per_article": docked, "independent_source_groups": 1}
(HERE / "results/met_screen.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: out[k] for k in ("mapping", "possible_missense_per_condition", "missense_variants_measured_in_all_12_conditions", "rows_with_non_numeric_effect")}))
for k, v in per.items():
    print(k, v)
for p, v in st.items():
    print(p, v["method"], v["resolution_A"], {c: (d["resolved_CA_in_scanned_span"], d["wt_identity_matches_P08581_at_pdb_number"]) for c, d in v["chains"].items()}, list(v["hetero_groups"])[:6])
