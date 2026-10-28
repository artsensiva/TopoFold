"""MET kinase domain: wild-type template site-context feature exporter (contract: CONTRACT.json).

Reads structural inputs (six experimental PDB complexes) and the pinned UniProt P08581 sequence only. It never reads, imports
or accepts an assay or outcome table. Features describe the wild-type template; they do not measure mutation-induced change.

    python3 met_export.py freeze      # writes freeze_exporter.json (hashes of code, contract, pinned data, tests)
    python3 met_export.py run         # verifies the freeze, then writes results/ (manifest, coverage, feature tables)

Reuse: the nine-Cα curvature and torsion come from the TopoFold core through `kingen2.core_kappa_tau` (checked against the
frozen numpy references); the PDB reader here is new because the existing stdlib reader drops ligands, insertion codes and
non-backbone residues, and the MDAnalysis adapters are trajectory tools that are not installed in the main environment.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "kin_gen_2"))
import kingen2  # noqa: E402  (core_kappa_tau: TopoFold core kappa/tau with numpy references)

CONTRACT_PATH = HERE / "CONTRACT.json"
AA20 = "ACDEFGHIKLMNPQRSTVWY"
ONE_TO_THREE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}
BACKBONE = ("N", "CA", "C", "O")
SIDE = {"A": ("CB",), "R": ("CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"), "N": ("CB", "CG", "OD1", "ND2"), "D": ("CB", "CG", "OD1", "OD2"), "C": ("CB", "SG"), "Q": ("CB", "CG", "CD", "OE1", "NE2"),
        "E": ("CB", "CG", "CD", "OE1", "OE2"), "G": (), "H": ("CB", "CG", "ND1", "CD2", "CE1", "NE2"), "I": ("CB", "CG1", "CG2", "CD1"), "L": ("CB", "CG", "CD1", "CD2"), "K": ("CB", "CG", "CD", "CE", "NZ"),
        "M": ("CB", "CG", "SD", "CE"), "F": ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"), "P": ("CB", "CG", "CD"), "S": ("CB", "OG"), "T": ("CB", "OG1", "CG2"),
        "W": ("CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"), "Y": ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"), "V": ("CB", "CG1", "CG2")}
VOLUME = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8, "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0}
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}
CA_CA_BREAK = 4.2
CONTACT_A = 4.5


def sha256_file(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_contract() -> dict:
    return json.loads(CONTRACT_PATH.read_text())


# --- reference sequence and mutation properties --------------------------------------------------------------------------------

def read_fasta(path) -> str:
    return "".join(l.strip() for l in Path(path).read_text().splitlines() if not l.startswith(">"))


def read_blosum62(path) -> dict:
    rows = [l.split() for l in Path(path).read_text().splitlines() if l.strip() and not l.startswith("#")]
    cols = rows[0]
    return {(r[0], c): int(v) for r in rows[1:] for c, v in zip(cols, r[1:])}


def mutation_table(blosum: dict) -> list:
    out = []
    for wt in AA20:
        for mut in AA20:
            if mut == wt:
                continue
            out.append({"wt": wt, "mut": mut, "vol_wt": VOLUME[wt], "vol_mut": VOLUME[mut], "dvol": round(VOLUME[mut] - VOLUME[wt], 4), "kd_wt": KD[wt], "kd_mut": KD[mut], "dkd": round(KD[mut] - KD[wt], 4),
                        "charge_wt": CHARGE.get(wt, 0), "charge_mut": CHARGE.get(mut, 0), "dcharge": CHARGE.get(mut, 0) - CHARGE.get(wt, 0), "blosum62": blosum[(wt, mut)],
                        "to_pro": int(mut == "P"), "from_pro": int(wt == "P"), "to_gly": int(mut == "G"), "from_gly": int(wt == "G")})
    return out


def row_universe(ref: str, contract: dict) -> list:
    lo, hi = contract["reference"]["span_uniprot"]
    off = contract["reference"]["dms_position_offset"]
    rows = []
    for p in range(lo, hi + 1):
        wt = ref[p - 1]
        for mut in AA20:
            if mut == wt:
                continue
            for d in contract["drugs"]:
                rows.append({"row_id": f"{d['code']}:{wt}{p}{mut}", "drug": d["code"], "uniprot_position": p, "dms_position": p - off, "wt": wt, "mut": mut})
    return rows


# --- structure reading ----------------------------------------------------------------------------------------------------------

def parse_pdb(text: str) -> tuple[dict, list]:
    """First model only, heavy atoms only. Residue key = (chain, resseq, icode). atoms: name -> [(altloc, xyz, occupancy)]."""
    res, order = {}, []
    for ln in text.splitlines():
        if ln.startswith("ENDMDL"):
            break
        if ln[:6] not in ("ATOM  ", "HETATM"):
            continue
        element = (ln[76:78].strip() if len(ln) >= 78 else "") or ln[12:16].strip().lstrip("0123456789")[:1]
        if element.upper() in ("H", "D"):
            continue
        key = (ln[21], int(ln[22:26]), ln[26])
        if key not in res:
            res[key] = {"resname": ln[17:20].strip(), "hetero": ln.startswith("HETATM"), "atoms": {}}
            order.append(key)
        xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        occ = float(ln[54:60]) if ln[54:60].strip() else 1.0
        res[key]["atoms"].setdefault(ln[12:16].strip(), []).append((ln[16], xyz, occ))
    return res, order


def atom(res: dict, name: str):
    """One conformer rule: blank location, else 'A', else missing."""
    cand = res["atoms"].get(name)
    if not cand:
        return None
    for want in (" ", "A"):
        for alt, xyz, _ in cand:
            if alt == want:
                return xyz
    return None


def has_altloc(res: dict) -> bool:
    return any(alt != " " for v in res["atoms"].values() for alt, _, _ in v)


def choose_chain(res: dict, order: list, span: tuple) -> str:
    """Contract chain_rule: most residues with a Cα, blank insertion code, number inside the span; ties to the first chain letter."""
    cnt = {}
    for k in order:
        if span[0] <= k[1] <= span[1] and k[2] == " " and atom(res[k], "CA") is not None:
            cnt[k[0]] = cnt.get(k[0], 0) + 1
    if not cnt:
        raise ValueError("no chain with Cα in the span")
    best = max(cnt.values())
    return sorted(c for c, n in cnt.items() if n == best)[0]


def find_ligand(res: dict, order: list, chain: str, component: str) -> tuple:
    keys = [k for k in order if res[k]["resname"] == component and k[0] == chain and res[k]["hetero"]]
    if len(keys) != 1:
        raise ValueError(f"expected exactly one {component} in chain {chain}, found {len(keys)}")
    k = keys[0]
    xyz = [a for name in res[k]["atoms"] for a in [atom(res[k], name)] if a is not None]
    if not xyz:
        raise ValueError("ligand without coordinates")
    return k, np.array(xyz)


def residue_table(res: dict, order: list, chain: str, ref: str, lo: int, hi: int) -> dict:
    """Status of every residue number lo..hi in the chosen chain (see CONTRACT residue_rules)."""
    ic_numbers = {k[1] for k in order if k[0] == chain and k[2] != " "}
    out = {}
    for q in range(lo, hi + 1):
        key = (chain, q, " ")
        info = {"status": "ok", "resname": None, "side_chain_complete": None, "altloc": False, "ca": None}
        if q in ic_numbers:
            info["status"] = "insertion_code"
        elif key not in res:
            info["status"] = "missing"
        else:
            r = res[key]
            info["resname"], info["altloc"] = r["resname"], has_altloc(r)
            if r["resname"] != ONE_TO_THREE[ref[q - 1]]:
                info["status"] = "engineered"
            else:
                ca = atom(r, "CA")
                if ca is None:
                    info["status"] = "missing_CA"
                else:
                    info["ca"] = ca
                    info["side_chain_complete"] = all(atom(r, n) is not None for n in BACKBONE + SIDE[ref[q - 1]])
        out[q] = info
    for q in range(lo, hi):
        a, b = out[q], out[q + 1]
        a["break_after"] = bool(a["ca"] is not None and b["ca"] is not None and float(np.linalg.norm(a["ca"] - b["ca"])) > CA_CA_BREAK)
    out[hi]["break_after"] = False
    return out


# --- features -------------------------------------------------------------------------------------------------------------------

KAPPA_COLS = [f"kappa_v{i}" for i in range(1, 8)]
TAU_COLS = [f"tau_sin_h{i}" for i in range(1, 7)] + [f"tau_cos_h{i}" for i in range(1, 7)]
IU = np.triu_indices(9, 1)
DIST_COLS = [f"d_{i}_{j}" for i, j in zip(*IU)]
L_COLS = ["L_min_heavy_dist_A", "L_min_sidechain_dist_A", "L_n_contact_atoms_4p5", "L_ca_min_dist_A"]
NAN = float("nan")


REASON = {"missing": "missing_residue", "engineered": "engineered_residue", "insertion_code": "insertion_code", "missing_CA": "missing_CA"}


def window_block(table: dict, p: int, half: int = 4):
    """(kappa, tau, dists, reason, first, last). Masked if any residue of the window is not ok or a break lies inside it."""
    qs = list(range(p - half, p + half + 1))
    reasons = []
    for q in qs:
        s = table[q]["status"] if q in table else "missing"
        if s != "ok":
            reasons.append(f"{REASON[s]}:{q}")
    for q in qs[:-1]:
        if q in table and table[q]["break_after"]:
            reasons.append(f"chain_break:{q}-{q + 1}")
    if reasons:
        return None, None, None, ";".join(reasons), qs[0], qs[-1]
    ca = np.array([table[q]["ca"] for q in qs])
    kap, tau, rk, rt = kingen2.core_kappa_tau(ca)
    if not (np.isfinite(kap).all() and np.isfinite(tau).all()):
        return None, None, None, "degenerate_geometry", qs[0], qs[-1]
    ref_k, ref_t = np.asarray(rk, float), np.asarray(rt, float)
    if np.abs(kap - ref_k).max() > 1e-4 or np.abs(np.sin(tau) - np.sin(ref_t)).max() > 1e-4:
        raise ValueError("core kappa/tau differ from the numpy reference beyond 1e-4")
    d = np.linalg.norm(ca[:, None] - ca[None], axis=-1)[IU]
    return kap, tau, d, "ok", qs[0], qs[-1]


def ligand_block(res: dict, chain: str, p: int, wt: str, status: dict, lig: np.ndarray):
    """(values[4], reason). Masked if residue p is not ok or has an incomplete side chain."""
    if status["status"] != "ok":
        return [NAN] * 4, f"{REASON[status['status']]}:{p}"
    if status["side_chain_complete"] is False:
        return [NAN] * 4, f"incomplete_side_chain:{p}"
    r = res[(chain, p, " ")]
    names = [n for n in BACKBONE + SIDE[wt]]
    pts = {n: atom(r, n) for n in names}
    allp = np.array([v for v in pts.values()])
    dall = np.linalg.norm(allp[:, None] - lig[None], axis=-1)
    sc = [pts[n] for n in SIDE[wt]]
    if sc:
        dsc = float(np.linalg.norm(np.array(sc)[:, None] - lig[None], axis=-1).min())
        reason = "ok"
    else:
        dsc, reason = NAN, "no_side_chain"
    ca_d = float(np.linalg.norm(pts["CA"] - lig, axis=-1).min())
    return [float(dall.min()), dsc, float((dall.min(axis=1) <= CONTACT_A).sum()), ca_d], reason


def build_entry(drug: dict, pdb_text: str, ref: str, contract: dict) -> dict:
    half = contract["window"]["half_width"]
    lo_s, hi_s = contract["reference"]["span_uniprot"]
    lo, hi = lo_s - half, hi_s + half
    res, order = parse_pdb(pdb_text)
    chain = choose_chain(res, order, (lo, hi))
    lkey, lig = find_ligand(res, order, chain, drug["inhibitor_component"])
    table = residue_table(res, order, chain, ref, lo, hi)
    rows = []
    for p in range(lo_s, hi_s + 1):
        wt = ref[p - 1]
        kap, tau, dist, t_reason, first, last = window_block(table, p, half)
        lvals, l_reason = ligand_block(res, chain, p, wt, table[p], lig)
        row = {"drug": drug["code"], "uniprot_position": p, "dms_position": p - contract["reference"]["dms_position_offset"], "wt": wt, "residue_status": table[p]["status"],
               "pdb": drug["pdb"], "chain": chain, "ligand_component": drug["inhibitor_component"], "ligand_resseq": lkey[1], "article_class": drug["article_class"],
               "window_first": first, "window_last": last, "altloc_in_window": sum(1 for q in range(first, last + 1) if q in table and table[q]["altloc"]), "L_reason": l_reason, "T_reason": t_reason, "C_reason": t_reason}
        row.update(dict(zip(L_COLS, lvals)))
        row.update(dict(zip(KAPPA_COLS, kap if kap is not None else [NAN] * 7)))
        row.update(dict(zip(TAU_COLS, list(np.sin(tau)) + list(np.cos(tau)) if tau is not None else [NAN] * 12)))
        row.update(dict(zip(DIST_COLS, dist if dist is not None else [NAN] * 36)))
        rows.append(row)
    flagged = {s: [q for q in range(lo, hi + 1) if table[q]["status"] == s] for s in ("missing", "insertion_code", "engineered", "missing_CA")}
    return {"rows": rows, "chain": chain, "ligand": {"component": drug["inhibitor_component"], "chain": lkey[0], "resseq": lkey[1], "icode": lkey[2], "heavy_atoms": int(len(lig))},
            "flagged_residues": flagged, "engineered_detail": {q: table[q]["resname"] for q in flagged["engineered"]},
            "incomplete_side_chain_residues": [q for q in range(lo_s, hi_s + 1) if table[q]["side_chain_complete"] is False], "altloc_residues": [q for q in range(lo, hi + 1) if table[q]["altloc"]],
            "chain_breaks_after": [q for q in range(lo, hi) if table[q]["break_after"]], "residues_considered": [lo, hi]}


def join_flat(universe: list, mutation: list, structure: list):
    """Row-level view used by analyses: universe x mutation properties x structural descriptors (no outcome columns)."""
    m = {(r["wt"], r["mut"]): r for r in mutation}
    s = {(r["drug"], r["uniprot_position"]): r for r in structure}
    for u in universe:
        out = dict(u)
        out.update({k: v for k, v in m[(u["wt"], u["mut"])].items() if k not in ("wt", "mut")})
        out.update({k: v for k, v in s[(u["drug"], u["uniprot_position"])].items() if k not in ("drug", "uniprot_position", "dms_position", "wt")})
        yield out


# --- export, manifest, freeze ----------------------------------------------------------------------------------------------------

FROZEN = ["met_export.py", "CONTRACT.json", "data/BLOSUM62.txt", "data/P08581.fasta", "../../../tests/test_met_export.py"]


def write_csv(path: Path, rows: list, cols: list):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(cols)
        for r in rows:
            w.writerow([("" if isinstance(r[c], float) and np.isnan(r[c]) else (f"{r[c]:.8g}" if isinstance(r[c], float) else r[c])) for c in cols])


def freeze():
    rec = {"schema": "topofold.met-export-freeze/1", "files": {f: sha256_file(HERE / f) for f in FROZEN}}
    (HERE / "freeze_exporter.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))


def check_freeze():
    fz = json.loads((HERE / "freeze_exporter.json").read_text())
    bad = [f for f, h in fz["files"].items() if sha256_file(HERE / f) != h]
    if bad:
        raise SystemExit(f"frozen file differs: {bad}")
    return fz


def run(outdir: Path | None = None):
    fz = check_freeze()
    contract = load_contract()
    ref = read_fasta(HERE / contract["reference"]["file"])
    if hashlib.sha256((HERE / contract["reference"]["file"]).read_bytes()).hexdigest() != contract["reference"]["sha256"]:
        raise SystemExit("reference sequence hash differs from the contract")
    outdir = outdir or HERE / "results"
    outdir.mkdir(exist_ok=True)
    structure, entries = [], {}
    for d in contract["drugs"]:
        f = HERE / d["input_file"]
        if not f.is_file():
            raise FileNotFoundError(f)
        if sha256_file(f) != d["input_sha256"]:
            raise SystemExit(f"input hash differs for {d['pdb']}")
        e = build_entry(d, f.read_text(), ref, contract)
        entries[d["code"]] = {k: v for k, v in e.items() if k != "rows"} | {"pdb": d["pdb"], "input_sha256": d["input_sha256"]}
        structure += e["rows"]
    mutation = mutation_table(read_blosum62(HERE / "data/BLOSUM62.txt"))
    universe = row_universe(ref, contract)
    s_cols = ["drug", "uniprot_position", "dms_position", "wt", "residue_status", "pdb", "chain", "ligand_component", "ligand_resseq", "article_class", "window_first", "window_last", "altloc_in_window"] + L_COLS + ["L_reason"] + KAPPA_COLS + TAU_COLS + ["T_reason"] + DIST_COLS + ["C_reason"]
    m_cols = ["wt", "mut", "vol_wt", "vol_mut", "dvol", "kd_wt", "kd_mut", "dkd", "charge_wt", "charge_mut", "dcharge", "blosum62", "to_pro", "from_pro", "to_gly", "from_gly"]
    write_csv(outdir / "features_structure.csv", structure, s_cols)
    write_csv(outdir / "features_mutation.csv", mutation, m_cols)
    write_csv(outdir / "row_universe.csv", universe, ["row_id", "drug", "uniprot_position", "dms_position", "wt", "mut"])
    cov = {}
    for d in contract["drugs"]:
        rs = [r for r in structure if r["drug"] == d["code"]]
        cov[d["code"]] = {}
        for blk, col in (("L", "L_reason"), ("T", "T_reason"), ("C", "C_reason")):
            ok = sum(r[col] == "ok" or (blk == "L" and r[col] == "no_side_chain") for r in rs)
            cats = {}
            for r in rs:
                if r[col] not in ("ok",):
                    for part in r[col].split(";"):
                        k = part.split(":")[0]
                        cats[k] = cats.get(k, 0) + 1
            cov[d["code"]][blk] = {"positions": len(rs), "ok": ok, "masked": len(rs) - ok, "reasons_incl_multiple_per_position": cats}
    with open(outdir / "coverage.csv", "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["drug", "pdb", "block", "positions", "ok", "masked", "reasons"])
        for d in contract["drugs"]:
            for blk in ("L", "T", "C"):
                c = cov[d["code"]][blk]
                w.writerow([d["code"], d["pdb"], blk, c["positions"], c["ok"], c["masked"], json.dumps(c["reasons_incl_multiple_per_position"], sort_keys=True)])
    man = {"schema": "topofold.met-export-manifest/1", "contract_version": contract["contract_version"], "scope": contract["scope"], "freeze": fz["files"], "units": contract["units"],
           "reference": contract["reference"], "row_universe": {"rows": len(universe), "positions": contract["reference"]["span_uniprot"][1] - contract["reference"]["span_uniprot"][0] + 1, "substitutions_per_position": 19, "drugs": len(contract["drugs"])},
           "tables": {n: {"rows": sum(1 for _ in open(outdir / n)) - 1, "sha256": sha256_file(outdir / n)} for n in ("features_structure.csv", "features_mutation.csv", "row_universe.csv", "coverage.csv")},
           "entries": entries, "coverage": cov, "block_definitions": contract["blocks"], "reason_vocabulary": contract["reason_vocabulary"], "source_terms": contract["source_terms"],
           "wider_engineered_template_effects_not_removed_by_masking": True}
    (outdir / "manifest.json").write_text(json.dumps(man, indent=1))
    print("exported", {k: v["rows"] for k, v in man["tables"].items()})
    return man


if __name__ == "__main__":
    {"freeze": freeze, "run": run}[sys.argv[1]]()
