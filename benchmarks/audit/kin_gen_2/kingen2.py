"""KIN-GEN-2 library: independent D1D2 reference groups, the 9-residue Cα window representations, one k-NN scoring procedure, weighted balanced accuracy.

Protocol: PROTOCOL.md / protocol.json (frozen before use). Reuses the existing `kinase_roles` helper (mapping), the TopoFold core through `compute_ca_features`
and `kinpair` (frozen KIN-PAIR-1 numerical references). Evidence scope: structural-label recovery only.
"""

from __future__ import annotations

import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/kin_pair_1"))
import kinpair as kp  # noqa: E402  (kappa, tau references; kabsch; frozen)
from topofold_inputs.features import compute_ca_features  # noqa: E402
from topofold_inputs.kinase_roles import MotifWindow, Numbering, ObservedResidue, Reference, RoleMappingError, _align, map_kinase_roles  # noqa: E402

STD20 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M",
       "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "MSE": "M"}
WINDOW_BEFORE, WINDOW_AFTER = 3, 5            # DFG-Asp - 3 .. DFG-Asp + 5: 9 residues
N_WIN = WINDOW_BEFORE + WINDOW_AFTER + 1
CA_CA_MAX, PEPTIDE_BOND_MAX = 4.2, 2.0
REPS = ("T", "B1", "B2", "B3")


def sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


# --- parsing (fixed columns, first model, altloc blank then A) ----------------------------------------------------------------

def parse_chains(text: str) -> dict:
    """{chain: [residue dict]} in file order for residues with N, CA and C atoms (amino-acid-like). Atoms by blank_then_A."""
    chains, index, ended = {}, {}, False
    for ln in text.splitlines():
        if ln.startswith("ENDMDL"):
            ended = True
        if ended or ln[:6] not in ("ATOM  ", "HETATM"):
            continue
        key = (ln[21], int(ln[22:26]), ln[26].strip())
        r = index.get(key)
        if r is None:
            r = index[key] = {"chain": key[0], "resid": key[1], "icode": key[2], "resname": ln[17:20].strip(), "atoms_all": {}}
            chains.setdefault(key[0], []).append(r)
        if ln[17:20].strip() != r["resname"]:
            continue
        r["atoms_all"].setdefault(ln[12:16].strip(), []).append((ln[16].strip(), np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    out = {}
    for ch, rs in chains.items():
        for r in rs:
            r["atoms"] = {n: _pick(v) for n, v in r.pop("atoms_all").items() if _pick(v) is not None}
        keep = [r for r in rs if all(a in r["atoms"] for a in ("N", "CA", "C"))]
        for i, r in enumerate(keep):
            prev = keep[i - 1] if i else None
            r["bonded_prev"] = bool(prev and float(np.linalg.norm(prev["atoms"]["C"][1] - r["atoms"]["N"][1])) <= PEPTIDE_BOND_MAX)
        out[ch] = keep
    return out


def _pick(recs):
    for want in ("", "A"):
        for alt, p in recs:
            if alt == want:
                return (alt, p)
    return None


def xyz(res, name):
    a = res["atoms"].get(name)
    return None if a is None else a[1]


def observed(res_list):
    return [ObservedResidue(r["resid"], r["icode"], r["resname"], r["bonded_prev"] if i else None) for i, r in enumerate(res_list)]


def letters(res_list):
    return "".join(ONE.get(r["resname"], "X") for r in res_list)


def pre_map(ref: str, res_list):
    """Alignment mapping reference index (1-based) -> residue dict, by the existing local-alignment function; plus identity and aligned count."""
    obs = letters(res_list)
    at_break = np.array([False] + [not r["bonded_prev"] for r in res_list[1:]], dtype=bool)
    score, pairs = _align(ref, obs, at_break)
    ident = sum(ref[i] == obs[j] for i, j in pairs)
    return {i + 1: res_list[j] for i, j in pairs}, (ident / len(pairs) if pairs else 0.0), len(pairs)


# --- UniProt-derived anchors ---------------------------------------------------------------------------------------------------

def uniprot_info(acc: str, cache: Path) -> dict:
    """Anchors by the rules of amendment_1.json (UniProt Protein kinase domain, ATP-binding Lys, motif windows measured from the Lys)."""
    f = cache / f"{acc}.json"
    if not f.is_file():
        cache.mkdir(parents=True, exist_ok=True)
        f.write_bytes(urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=60).read())
    raw = f.read_bytes()
    u = json.loads(raw)
    ref = u["sequence"]["value"]
    feats = u["features"]
    info = {"accession": acc, "entry": u["uniProtkbId"], "length": len(ref), "sequence": ref, "sequence_sha256": sha256(ref.encode()), "raw_sha256": sha256(raw),
            "entry_version": u["entryAudit"]["entryVersion"], "problems": []}
    atp_lys = [x["location"]["start"]["value"] for x in feats if x["type"] == "Binding site" and x.get("ligand", {}).get("name") == "ATP"
               and x["location"]["start"]["value"] == x["location"]["end"]["value"] and ref[x["location"]["start"]["value"] - 1] == "K"]
    pairs = []
    for d in feats:
        if d["type"] == "Domain" and "protein kinase" in d.get("description", "").lower():
            d0, d1 = d["location"]["start"]["value"], d["location"]["end"]["value"]
            ks = [p for p in atp_lys if d0 + 20 <= p <= d0 + 40]
            if len(ks) == 1:
                pairs.append((d0, d1, ks[0], d.get("description", "")))
    if len(pairs) != 1:
        info["problems"].append(f"{len(pairs)} Protein kinase domains with exactly one ATP-binding Lys at domain start + 20..+40 (need exactly one)")
        return info
    d0, d1, k, name = pairs[0]
    hrd = [i for i in range(k + 70, k + 126) if ref[i - 1:i + 2] in ("HRD", "YRD")]
    dfg = [i for i in range(k + 95, k + 141) if ref[i - 1:i + 2] == "DFG"]
    info.update({"domain": [d0, d1], "domain_name": name, "beta3_lys": k, "hrd_window": [k + 70, k + 125], "dfg_window": [k + 95, k + 140], "hrd_hits": hrd, "dfg_hits": dfg, "domain_sequence": ref[d0 - 1:d1]})
    if len(hrd) != 1:
        info["problems"].append(f"{len(hrd)} HRD/YRD motifs in the window {k + 70}-{k + 125} (need exactly one)")
    if len(dfg) != 1:
        info["problems"].append(f"{len(dfg)} DFG motifs in the window {k + 95}-{k + 140} (need exactly one)")
    if len(hrd) == 1:
        info["hrd_h"] = hrd[0]
    return info


def alphac_glu(info: dict, mapped: list) -> dict:
    """mapped: list of {reference index -> residue dict} for the kinase's entries (one chain each). Contact rule of protocol.json."""
    k = info["beta3_lys"]
    counts = {}
    for m in mapped:
        lys = m.get(k)
        nz = xyz(lys, "NZ") if lys else None
        if nz is None:
            continue
        for idx in range(k + 12, k + 26):
            if info["sequence"][idx - 1] != "E" or idx not in m:
                continue
            d = [float(np.linalg.norm(xyz(m[idx], o) - nz)) for o in ("OE1", "OE2") if xyz(m[idx], o) is not None]
            if d and min(d) <= 4.5:
                counts[idx] = counts.get(idx, 0) + 1
    order = sorted(counts.items(), key=lambda t: (-t[1], t[0]))
    ok = bool(order) and order[0][1] >= 3 and (len(order) == 1 or order[0][1] - order[1][1] >= 2)
    return {"counts": counts, "selected": order[0][0] if ok else None, "reason": "" if ok else f"no unique contact Glu (counts {counts})"}


# --- one entry: roles, window, reference distances --------------------------------------------------------------------------

def analyse_chain(info: dict, glu: int, res_list: list) -> dict:
    """Existing helper mapping with the kinase anchors, then the window and D1, D2. Returns {'ok': bool, 'reason': str, ...}."""
    ref = Reference(info["accession"], info["entry"], info["sequence"], f"UniProt {info['accession']} v{info['entry_version']}")
    anchors = {"beta3_Lys": info["beta3_lys"], "alphaC_Glu": glu, "HRD_H": MotifWindow(*info["hrd_window"]), "DFG_Asp": MotifWindow(*info["dfg_window"])}
    try:
        m = map_kinase_roles(observed(res_list), res_list[0]["chain"], ref, Numbering.reference_index(), anchors)
    except RoleMappingError as e:
        return {"ok": False, "reason": "kinase_roles refused: " + "; ".join(e.problems)[:240]}
    except Exception as e:     # noqa: BLE001
        return {"ok": False, "reason": f"mapping error: {str(e)[:160]}"}
    key = lambda spec: (spec.resid, spec.icode)   # noqa: E731
    pos = {(r["resid"], r["icode"]): i for i, r in enumerate(res_list)}
    roles = m.roles
    if any(roles[r] is None for r in ("beta3_Lys", "alphaC_Glu_plus4", "DFG_Asp", "DFG_Phe")):
        return {"ok": False, "reason": "a reference role is unobserved in this chain"}
    jd = pos[key(roles["DFG_Asp"])]
    lo, hi = jd - WINDOW_BEFORE, jd + WINDOW_AFTER
    if lo < 0 or hi >= len(res_list):
        return {"ok": False, "reason": "window runs past the chain ends"}
    win = res_list[lo:hi + 1]
    if not all(r["bonded_prev"] for r in win[1:]):
        return {"ok": False, "reason": "window has a chain break (no peptide bond between consecutive residues)"}
    ca = np.array([xyz(r, "CA") for r in win])
    if np.linalg.norm(np.diff(ca, axis=0), axis=1).max() > CA_CA_MAX:
        return {"ok": False, "reason": "window Cα-Cα distance above 4.2 Å"}
    nums = [(r["resid"], r["icode"]) for r in win]
    if any(not _consec(a, b) for a, b in zip(nums, nums[1:])):
        return {"ok": False, "reason": "window numbering is not consecutive"}
    odd = [r["resname"] for r in win if r["resname"] not in STD20 and r["resname"] != "MSE"]
    if odd:
        return {"ok": False, "reason": f"non-standard residue name(s) in the window: {sorted(set(odd))}"}
    by = lambda role: res_list[pos[key(roles[role])]]   # noqa: E731
    k_ca, e4_ca, phe_cz = xyz(by("beta3_Lys"), "CA"), xyz(by("alphaC_Glu_plus4"), "CA"), xyz(by("DFG_Phe"), "CZ")
    if k_ca is None or e4_ca is None or phe_cz is None:
        return {"ok": False, "reason": "Cα of the β3 Lys or αC Glu+4, or Cζ of the DFG Phe, is not observed"}
    return {"ok": True, "reason": "", "window_ca": ca.tolist(), "window_ids": [f"{r['chain']}:{r['resid']}{r['icode']}" for r in win],
            "window_resnames": [r["resname"] for r in win], "D1": float(np.linalg.norm(e4_ca - phe_cz)), "D2": float(np.linalg.norm(k_ca - phe_cz)),
            "anchor_residues": {k: f"{res_list[pos[key(v)]]['chain']}:{v.resid}{v.icode}" for k, v in roles.items() if v is not None}, "identity": m.report["alignment"]["identity"]}


def _consec(a, b):
    if b[0] == a[0]:
        return b[1] == ("A" if a[1] == "" else chr(ord(a[1]) + 1))
    return b[0] == a[0] + 1 and b[1] == ""


# --- representations --------------------------------------------------------------------------------------------------------

def core_kappa_tau(ca: np.ndarray):
    """kappa at the 7 interior vertices and tau at the 6 hinges, from the TopoFold core (compute_ca_features), plus the numpy reference."""
    X = np.asarray(ca, float)[None]
    f = compute_ca_features(X, np.ones(N_WIN, bool), [(0, N_WIN)])
    kap, tau = f.kappa[0, 1:N_WIN - 1], f.tau[0, 1:N_WIN - 2]
    ref_k = np.array([kp.kappa(*X[0, i - 1:i + 2]) for i in range(1, N_WIN - 1)], float)
    ref_t = np.array([kp.tau(*X[0, i - 1:i + 3]) for i in range(1, N_WIN - 2)], float)
    return kap, tau, ref_k, ref_t


def raw_features(ca: np.ndarray) -> dict:
    ca = np.asarray(ca, float)
    kap, tau, _, _ = core_kappa_tau(ca)
    d2 = [np.linalg.norm(ca[i + 1] - ca[i - 1]) for i in range(1, N_WIN - 1)]
    d3 = [np.linalg.norm(ca[i + 2] - ca[i - 1]) for i in range(1, N_WIN - 2)]
    iu = np.triu_indices(N_WIN, 1)
    allp = np.linalg.norm(ca[:, None] - ca[None], axis=-1)[iu]
    return {"T": np.concatenate([kap, np.sin(tau), np.cos(tau)]), "B1": np.array(d2 + d3), "B3": allp}


def rmsd(a, b):
    return kp.rmsd(np.asarray(a, float), np.asarray(b, float))


# --- reference groups and classifier ---------------------------------------------------------------------------------------------

def fit_groups(points: np.ndarray) -> dict:
    """Average-linkage hierarchical clustering of (D1, D2) into three groups; named by the rule of protocol.json."""
    from scipy.cluster.hierarchy import fcluster, linkage
    lab = fcluster(linkage(points, method="average"), 3, criterion="maxclust")
    ids = sorted(set(lab))
    mean = {c: points[lab == c].mean(0) for c in ids}
    n = {c: int((lab == c).sum()) for c in ids}
    out_c = max(ids, key=lambda c: mean[c][0])
    rest = [c for c in ids if c != out_c]
    in_c = max(rest, key=lambda c: n[c])
    inter_c = [c for c in rest if c != in_c][0] if len(rest) == 2 else None
    names = {out_c: "out", in_c: "in", **({inter_c: "inter"} if inter_c is not None else {})}
    cent = {names[c]: mean[c].tolist() for c in ids}
    return {"centroids": cent, "sizes": {names[c]: n[c] for c in ids}, "ranges": {names[c]: {"D1": [float(points[lab == c][:, 0].min()), float(points[lab == c][:, 0].max())],
                                                                                         "D2": [float(points[lab == c][:, 1].min()), float(points[lab == c][:, 1].max())]} for c in ids},
            "labels": [names[c] for c in lab]}


def assign_group(d1d2, groups):
    c = groups["centroids"]
    return min(c, key=lambda g: float(np.hypot(d1d2[0] - c[g][0], d1d2[1] - c[g][1])))


class Model:
    """Development-only classifier: scaling, neighbours, k. Holds development ids only (leakage control)."""

    def __init__(self, dev: list, y: list, k: dict | None = None):
        self.ids = [d["id"] for d in dev]
        self.y = np.array([1 if v == "out" else 0 for v in y])
        self.dev_ca = [np.asarray(d["window_ca"], float) for d in dev]
        raws = [raw_features(c) for c in self.dev_ca]
        self.mean, self.std, self.F = {}, {}, {}
        for r in ("T", "B1", "B3"):
            M = np.array([x[r] for x in raws])
            self.mean[r], self.std[r] = M.mean(0), M.std(0) + 1e-9
            self.F[r] = (M - self.mean[r]) / self.std[r]
        self.k = k or {}

    def dist_to_dev(self, rep, ca, exclude=None):
        ca = np.asarray(ca, float)
        if rep == "B2":
            d = np.array([rmsd(ca, c) for c in self.dev_ca])
        else:
            v = (raw_features(ca)[rep] - self.mean[rep]) / self.std[rep]
            d = np.linalg.norm(self.F[rep] - v, axis=1)
        if exclude is not None:
            d = np.where(exclude, np.inf, d)
        return d

    def score(self, rep, ca, k, exclude=None):
        d = self.dist_to_dev(rep, ca, exclude)
        nn = np.argsort(d, kind="stable")[:k]
        return float(self.y[nn].mean())


def balanced_accuracy(y, pred):
    y, pred = np.asarray(y), np.asarray(pred)
    tpr = float((pred[y == 1] == 1).mean()) if (y == 1).any() else float("nan")
    tnr = float((pred[y == 0] == 0).mean()) if (y == 0).any() else float("nan")
    return (tpr + tnr) / 2, tpr, tnr


def loko(model: Model, kinases: list, reps=REPS, k_grid=(1, 3, 5, 7)) -> dict:
    """Leave-one-kinase-out on development entries: the held-out kinase's entries never act as neighbours."""
    kin = np.array(kinases)
    res = {}
    for rep in reps:
        res[rep] = {}
        for k in k_grid:
            pred = []
            for i, ca in enumerate(model.dev_ca):
                pred.append(1 if model.score(rep, ca, k, exclude=(kin == kin[i])) >= 0.5 else 0)
            ba, tpr, tnr = balanced_accuracy(model.y, pred)
            res[rep][k] = {"ba": ba, "tpr": tpr, "tnr": tnr}
    return res


def weighted_ba(rows: list) -> dict:
    """rows: dicts with kinase, y (1 out / 0 in), pred. Each (kinase, class) cell has equal weight within its class."""
    out = {}
    for cls, name in ((1, "tpr"), (0, "tnr")):
        cells = {}
        for r in rows:
            if r["y"] == cls:
                cells.setdefault(r["kinase"], []).append(1.0 if r["pred"] == cls else 0.0)
        out[name] = float(np.mean([np.mean(v) for v in cells.values()])) if cells else float("nan")
        out[f"n_kinases_{name}"] = len(cells)
    out["ba_w"] = (out["tpr"] + out["tnr"]) / 2
    return out


def auroc(y, s):
    y, s = np.asarray(y), np.asarray(s, float)
    pos, neg = s[y == 1], s[y == 0]
    if not len(pos) or not len(neg):
        return float("nan")
    return float(np.mean([(p > q) + 0.5 * (p == q) for p in pos for q in neg]))


def check_freeze(freeze_file: Path) -> list:
    """Hash check of a freeze record. PROTOCOL.md may only have gained appended, dated amendments (the frozen text, minus its trailing '(none)', must be a prefix)."""
    import subprocess
    fr = json.loads(freeze_file.read_text())
    bad = []
    for f, h in fr["sha256"].items():
        cur = (HERE / f).read_bytes()
        if sha256(cur) == h:
            continue
        if f == "PROTOCOL.md":
            commit = subprocess.run(["git", "log", "--format=%h", "--diff-filter=A", "--", str(freeze_file.relative_to(ROOT))], capture_output=True, text=True, cwd=ROOT).stdout.split()[0]
            frozen = subprocess.run(["git", "show", f"{commit}:{(HERE / f).relative_to(ROOT)}"], capture_output=True, cwd=ROOT).stdout
            base = frozen.rstrip()
            base = base[: -len(b"(none)")].rstrip() if base.endswith(b"(none)") else base
            if sha256(frozen) == h and cur.startswith(base):
                continue
        bad.append(f)
    return bad
