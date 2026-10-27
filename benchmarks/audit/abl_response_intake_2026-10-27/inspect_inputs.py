"""Read-only inspection of the Hauser et al. 2018 structural inputs (Supplementary Data 5) and endpoint table (Supplementary Data 3).

Writes results/inputs_inspection.json. No model is trained or scored. Run: python3 -I inspect_inputs.py
Dependencies: numpy and the repository's kinpair.rmsd via kin_gen_2.kingen2 (Kabsch).
"""
import collections
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "kin_gen_2"))
sys.path.insert(0, str(HERE.parent / "kin_pair_1"))
from xlsx_min import read_xlsx  # noqa: E402
import kingen2  # noqa: E402

S5 = HERE / "data_raw/hauser_s5_extracted/Supplementary_Data_5"
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "ASH": "D", "GLH": "E", "LYN": "K"}
KINASE_DOMAIN = (229, 512)  # per the paper: ABL1 isoform IA, S229-K512
METHOD = {"axitinib": "co-crystal 4WA9 + homology/Prime completion", "bosutinib": "co-crystal 3UE4 + homology/Prime completion", "dasatinib": "co-crystal 4XEY + homology/Prime completion",
          "imatinib": "co-crystal 1OPJ + Prime completion", "nilotinib": "co-crystal 3CS9 + homology/Prime completion", "ponatinib": "co-crystal 3OXZ + homology/Prime completion",
          "erlotinib": "docked into 3UE4 receptor (Glide-SP)", "gefitnib": "docked into 3UE4 receptor (Glide-SP)"}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_pdb(p):
    res, order = {}, []
    het = collections.Counter()
    for ln in Path(p).read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            name, rn, ch = ln[12:16].strip(), ln[17:20].strip(), ln[21]
            num, ic = int(ln[22:26]), ln[26]
            key = (ch, num, ic)
            xyz = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            el = (ln[76:78].strip() or name[0])
            if key not in res:
                res[key] = {"rn": rn, "atoms": {}, "het": ln.startswith("HETATM"), "alt": set()}
                order.append(key)
            if ln[16] not in " A":
                continue
            res[key]["atoms"][name] = (xyz, el)
            if ln[16] != " ":
                res[key]["alt"].add(ln[16])
    return res, order


def seq_of(res, order, chain):
    return [(k[1], THREE.get(res[k]["rn"], "X"), res[k]["rn"]) for k in order if k[0] == chain and res[k]["rn"] in THREE and "CA" in res[k]["atoms"]]


uni = "".join(l.strip() for l in open(HERE / "data_raw/P00519.fasta") if not l.startswith(">"))
# endpoint table
rows = read_xlsx(HERE / "data_raw/hauser/Supplementary Data 3.xlsx")["Sheet1"]
full = [r for r in rows[1:] if len(r) == 11]
tki_n = collections.Counter(r[0] for r in full)
muts = sorted({r[1] for r in full}, key=lambda m: int(m[1:-1]))
mut_ok = {m: (uni[int(m[1:-1]) - 1] == m[0]) for m in muts}
sites = sorted({int(m[1:-1]) for m in muts})
sites_all_sheet = sorted({int(r[1][1:-1]) for r in rows[1:]})
per_site = collections.Counter(int(m[1:-1]) for m in muts)
table = {
    "rows_in_sheet": len(rows) - 1, "rows_with_full_values": len(full), "rows_without_values": len(rows) - 1 - len(full),
    "inhibitors": dict(tki_n), "n_inhibitors": len(tki_n), "mutations_with_values": len(muts), "sites_with_values": len(sites), "sites_in_sheet_all_rows": len(sites_all_sheet),
    "mutations_per_site": dict(sorted(per_site.items())), "wildtype_residue_matches_P00519_1a_numbering": {"matches": sum(mut_ok.values()), "of": len(mut_ok), "mismatch": [m for m, v in mut_ok.items() if not v]},
    "ic50_nM_min": min(r[2] for r in full), "ic50_nM_max": max(r[2] for r in full), "n_ic50_at_10000": sum(1 for r in full if r[2] == 10000.0),
    "n_ic50_ge_10000": sum(1 for r in full if r[2] >= 10000.0), "ddg_exp_min": min(r[3] for r in full), "ddg_exp_max": max(r[3] for r in full),
    "n_resistant_ge_1.36": sum(1 for r in full if r[3] >= 1.36), "unique_ic50_values_per_inhibitor": {t: len({r[2] for r in full if r[0] == t}) for t in tki_n},
}
# consistency of DDG_exp with IC50: DDG = RT ln(IC50_mut / IC50_WT) (positive = resistance). The WT IC50 is not in the sheet; per inhibitor the implied WT IC50 is the median of IC50 * exp(-DDG/RT)
RT = 0.0019872041 * 300.0
imp, flagged = {}, []
for t in tki_n:
    rr = [r for r in full if r[0] == t]
    wt = float(np.median([r[2] * np.exp(-r[3] / RT) for r in rr]))
    resid = {r[1]: RT * np.log(r[2] / wt) - r[3] for r in rr}
    imp[t] = {"implied_wt_ic50_nM_median": wt, "max_abs_residual_kcal": float(max(abs(v) for v in resid.values())), "n_rows_residual_gt_0.05": int(sum(abs(v) > 0.05 for v in resid.values()))}
    flagged += [{"inhibitor": t, "mutation": m, "ic50_nM": next(r[2] for r in rr if r[1] == m), "ddg_exp_in_sheet": next(r[3] for r in rr if r[1] == m), "ddg_from_ic50_and_implied_wt": float(RT * np.log(next(r[2] for r in rr if r[1] == m) / wt)), "residual_kcal": float(v)} for m, v in resid.items() if abs(v) > 0.05]
table["implied_wt_ic50_consistency_per_inhibitor"] = imp
table["rows_where_ddg_exp_disagrees_with_ic50_column"] = flagged
table["ic50_resolution_note"] = "IC50 values are given to integer nM in most rows; rows sharing one IC50 per inhibitor: " + json.dumps({t: sum(c - 1 for c in collections.Counter(r[2] for r in full if r[0] == t).values() if c > 1) for t in tki_n})


# provenance of the erlotinib/gefitinib values: compare with the Davis et al. block of Supplementary Data 2 (non-phosphorylated vs phosphorylated columns)
dv = read_xlsx(HERE / "data_raw/hauser/Supplementary Data 2.xlsx")["Davis et al. data"]
col = {"erlotinib": 9, "gefitinib": 10}
davis = {"phosphorylated": {}, "non_phosphorylated": {}}
for r in dv[1:]:
    if len(r) > 3 and r[3] and "ABL1(" in str(r[3]):
        m = re.match(r"ABL1\(([A-Z]\d+[A-Z])\)", r[3])
        if m:
            for inh, c in col.items():
                for blk, off in (("non_phosphorylated", 0), ("phosphorylated", 8)):
                    v = r[c + off] if len(r) > c + off else None
                    if v not in (None, ""):
                        davis[blk][(inh, m.group(1))] = float(v)
dmatch = {}
for inh in col:
    rr = [(r[1], r[2]) for r in full if r[0] == inh]
    dmatch[inh] = {blk: sum(1 for m, v in rr if davis[blk].get((inh, m)) == v) for blk in davis}
    dmatch[inh]["n_rows"] = len(rr)
table["erlotinib_gefitinib_values_equal_to_Davis_block"] = dmatch
# site-level structure of the labels (resistant = DDG_exp >= 1.36 kcal/mol, the threshold used by both papers)
bys = collections.defaultdict(list)
for r in full:
    bys[int(r[1][1:-1])].append(r[3] >= 1.36)
n_pos_sites = sum(1 for v in bys.values() if any(v))
table["resistant_records"] = sum(1 for r in full if r[3] >= 1.36)
table["sites_with_any_resistant_record"] = n_pos_sites
table["sites_with_any_resistant_record_list"] = sorted(s_ for s_, v in bys.items() if any(v))
table["mutations_with_any_resistant_record"] = len({r[1] for r in full if r[3] >= 1.36})
table["inhibitors_with_any_resistant_record"] = sorted({r[0] for r in full if r[3] >= 1.36})
A, n1, n0 = 0.8, n_pos_sites, len(bys) - n_pos_sites
Q1, Q2 = A / (2 - A), 2 * A * A / (1 + A)
table["site_level_auroc_standard_error_if_true_auroc_0.8_hanley_mcneil"] = float(np.sqrt((A * (1 - A) + (n1 - 1) * (Q1 - A * A) + (n0 - 1) * (Q2 - A * A)) / (n1 * n0)))

# structures
files = sorted((S5 / "renumbered").glob("*_nowat-uniprot.pdb")) + [S5 / "renumbered/cAbl-dasatinib-uniprot.pdb"]
files = sorted(set(files))
struct, cas = {}, {}
for f in files:
    inh = f.name.split("-")[1].split("_")[0]
    res, order = read_pdb(f)
    chains = collections.Counter(k[0] for k in order)
    prot_chain = max(chains, key=lambda c: sum(1 for k in order if k[0] == c and res[k]["rn"] in THREE))
    sq = seq_of(res, order, prot_chain)
    nums = [n for n, _, _ in sq]
    mism = [(n, a, uni[n - 1]) for n, a, _ in sq if not (0 < n <= len(uni)) or uni[n - 1] != a]
    gaps = [(a, b) for a, b in zip(nums, nums[1:]) if b != a + 1]
    sitecover = {}
    for s in sites:
        k = (prot_chain, s, " ")
        sitecover[s] = bool(k in res and "CA" in res[k]["atoms"])
    nonstd = sorted({res[k]["rn"] for k in order if k[0] == prot_chain and res[k]["rn"] not in THREE})
    hetero = sorted({res[k]["rn"] for k in order if k[0] != prot_chain or res[k]["het"]})
    heavy = sum(1 for k in order if k[0] == prot_chain for a, (x, el) in res[k]["atoms"].items() if el != "H")
    # ligand-contact counts (heavy atoms within 4.5 A of ligand) per mutation site, from the wild-type complex only
    lig = [xyz for k in order if k[0] != prot_chain for a, (xyz, el) in res[k]["atoms"].items() if el != "H"]
    mind = {}
    if lig:
        L = np.array(lig)
        for s in sites:
            k = (prot_chain, s, " ")
            if k in res:
                P = np.array([xyz for a, (xyz, el) in res[k]["atoms"].items() if el != "H"])
                mind[s] = float(np.sqrt(((P[:, None] - L[None]) ** 2).sum(-1)).min())
    struct[inh] = {
        "file": str(f.relative_to(HERE)), "sha256": sha(f), "protein_chain": prot_chain, "n_chains": dict(chains), "construction": METHOD[inh],
        "n_standard_residues_with_CA": len(sq), "first_last": [nums[0], nums[-1]], "gaps_in_numbering": gaps, "mismatches_to_P00519_1a": mism,
        "nonstandard_protein_residue_names": nonstd, "hetero_groups_and_ligand_residues": hetero, "protein_heavy_atoms": heavy, "has_hydrogens": any(el == "H" for k in order for a, (x, el) in res[k]["atoms"].items()),
        "mutation_sites_with_CA": sum(sitecover.values()), "mutation_sites_total": len(sites), "min_ligand_distance_A_per_site": mind,
    }
    cas[inh] = {n: np.array(res[(prot_chain, n, " ")]["atoms"]["CA"][0]) for n, _, _ in sq}
# encoded wild-type identity of every mutation in every structure
enc = {}
for inh, s in struct.items():
    res, order = read_pdb(S5 / "renumbered" / Path(s["file"]).name)
    ok = 0
    tot = 0
    for m in muts:
        k = (s["protein_chain"], int(m[1:-1]), " ")
        tot += 1
        ok += int(k in res and THREE.get(res[k]["rn"]) == m[0])
    enc[inh] = {"mutation_wt_identity_present": ok, "of": tot}
# Ca-trace comparison across the eight wild-type models on the residues common to all (rigid Kabsch)
common = sorted(set.intersection(*[set(c) for c in cas.values()]))
names = sorted(cas)
M = {n: np.array([cas[n][i] for i in common]) for n in names}
pair = {}
for i, a in enumerate(names):
    for b in names[i + 1:]:
        pair[f"{a}|{b}"] = float(kingen2.rmsd(M[a], M[b]))
# precision floor: PDB-format coordinates carry 3 decimals -> rounding error per axis <= 0.0005 A; RMSD floor taken as 0.001 A
vals = np.array(list(pair.values()))
# wat vs nowat and receptors vs renumbered: identical protein Ca?
det = {}
for inh in ("axitinib", "bosutinib", "imatinib", "nilotinib", "ponatinib", "erlotinib", "gefitnib"):
    a = S5 / "renumbered" / f"cAbl-{inh}_nowat-uniprot.pdb"
    b = S5 / "renumbered" / f"cAbl-{inh}_wat-uniprot.pdb"
    ra, oa = read_pdb(a)
    rb, ob = read_pdb(b)
    ch = struct[inh]["protein_chain"]
    ka = [k for k in oa if k[0] == ch and "CA" in ra[k]["atoms"]]
    kb = [k for k in ob if k[0] == ch and "CA" in rb[k]["atoms"]]
    det[inh] = {"same_residues": ka == kb, "max_abs_coordinate_difference_A": float(max(np.abs(np.array(ra[k]["atoms"]["CA"][0]) - np.array(rb[k]["atoms"]["CA"][0])).max() for k in ka)) if ka == kb else None}
def mol2_heavy(p):
    atoms, on = [], False
    for ln in Path(p).read_text().splitlines():
        if ln.startswith("@<TRIPOS>"):
            on = ln.startswith("@<TRIPOS>ATOM")
            continue
        if on and ln.strip():
            f = ln.split()
            if not f[5].startswith("H"):
                atoms.append((float(f[2]), float(f[3]), float(f[4])))
    return np.array(atoms)


LIGFILE = {"axitinib": "axitinib", "bosutinib": "bosutinib", "dasatinib": "dasatinib", "imatinib": "imatinib", "nilotinib": "nilotinib", "ponatinib": "ponatinib", "erlotinib": "erlotinib_nowat", "gefitnib": "gefitinib_nowat"}
ligframe = {}
for inh, s_ in struct.items():
    L = mol2_heavy(S5 / "ligands" / (LIGFILE[inh] + ".mol2"))
    res, order = read_pdb(S5 / "renumbered" / Path(s_["file"]).name)
    P = np.array([xyz for k in order if k[0] == s_["protein_chain"] for a, (xyz, el) in res[k]["atoms"].items() if el != "H"])
    d = np.sqrt(((L[:, None] - P[None]) ** 2).sum(-1))
    sited = {}
    for s0 in sites:
        k = (s_["protein_chain"], s0, " ")
        Q = np.array([xyz for a, (xyz, el) in res[k]["atoms"].items() if el != "H"])
        sited[s0] = float(np.sqrt(((L[:, None] - Q[None]) ** 2).sum(-1)).min())
    ligframe[inh] = {"ligand_heavy_atoms": int(len(L)), "min_ligand_protein_heavy_atom_distance_A": float(d.min()), "sites_within_4.5A_of_ligand": sorted(k for k, v in sited.items() if v <= 4.5), "n_sites_within_4.5A": int(sum(v <= 4.5 for v in sited.values())), "n_sites_beyond_10A": int(sum(v > 10 for v in sited.values()))}
    s_["min_ligand_distance_A_per_site"] = sited
out = {
    "ligand_frame_and_site_contacts": ligframe,
    "schema": "topofold.abl-response-inputs-inspection/1", "endpoint_table": table, "structures": struct, "encoded_wt_identity": enc,
    "ca_trace_comparison": {"residues_common_to_all_models": len(common), "first_last_common": [common[0], common[-1]], "pairwise_rmsd_A": pair, "min_A": float(vals.min()), "max_A": float(vals.max()), "median_A": float(np.median(vals)),
                            "coordinate_rounding_floor_A": 0.001, "note": "Only wild-type models exist in the deposit; no mutant coordinates, so no mutant-versus-wild-type Ca comparison is possible from it."},
    "wat_vs_nowat_identical_protein_CA": det,
    "mutant_structures_in_deposit": {"files_matching_mutant_or_mutation_in_names": [p.name for p in S5.rglob("*") if re.search(r"mut|[A-Z]\d{3}[A-Z]", p.name)], "trajectory_files": [p.name for p in S5.rglob("*") if p.suffix in (".xtc", ".trr", ".dcd", ".nc", ".h5")]},
}
(HERE / "results").mkdir(exist_ok=True)
(HERE / "results/inputs_inspection.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: out[k] for k in ("ligand_frame_and_site_contacts",)}, indent=0)[:2500])
print(json.dumps({k: table[k] for k in ("implied_wt_ic50_consistency_per_inhibitor","rows_where_ddg_exp_disagrees_with_ic50_column","ic50_resolution_note")}, indent=0))
for k, v in struct.items():
    print(k, v["protein_chain"], v["n_chains"], v["first_last"], v["gaps_in_numbering"], v["mismatches_to_P00519_1a"], v["nonstandard_protein_residue_names"], v["hetero_groups_and_ligand_residues"][:4], v["mutation_sites_with_CA"], v["has_hydrogens"])
