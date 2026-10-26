"""KIN-PAIR-1 helpers: independent PDB parsing, sequence mapping checks and independent numerical references.

Everything here is written independently of the structpair workflow under test, which reads coordinates through
MDAnalysis and the TopoFold core. It uses a fixed-column PDB parser, numpy only, and the existing local-alignment function
of `topofold_inputs.kinase_roles` (reuse; it is the mapping helper developed for KIN-MAP-1). No TopoFold score is computed here.
Protocol: PROTOCOL.md (frozen before selection or execution).
"""

from __future__ import annotations

import hashlib
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs.kinase_roles import _align  # noqa: E402  (existing local alignment with chain-break aware gap costs)

KINPAIR_VERSION = "1.0.0"
THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
             "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V", "MSE": "M", "HYP": "P", "CME": "C", "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H",
             "HSE": "H", "HSP": "H"}
# Residue names that the Cα adapter's `protein` selection keeps (MDAnalysis 2.10.0 ProteinSelection.prot_res, 104 names).
# A residue with another name (for example TPO, SEP, PTR) is not read by the adapter: for this study it is "not read".
PROT_RES = frozenset("""ACE ALA ALAD ARG ARGN ASF ASH ASN ASN1 ASP ASPH CALA CARG CASF CASN CASP CCYS CCYX CGLN CGLU CGLY CHID CHIE
CHIP CILE CLEU CLYS CME CMET CPHE CPRO CSER CTHR CTRP CTYR CVAL CYM CYS CYS1 CYS2 CYSH CYX DAB GLH GLN GLU GLUH GLY HID HIE HIP
HIS HIS1 HIS2 HISA HISB HISD HISE HISH HSD HSE HSP HYP ILE LEU LYN LYS LYSH MET MSE NALA NARG NASN NASP NCYS NCYX NGLN NGLU
NGLY NHID NHIE NHIP NILE NLEU NLYS NMET NPHE NPRO NSER NTHR NTRP NTYR NVAL ORN PGLU PHE PRO QLN SER THR TRP TYR VAL""".split())
CA_CA_MAX = 4.2                 # Å, the adapter's documented break threshold
PEPTIDE_BOND_MAX = 2.0          # Å, C(i-1)-N(i)
COLLINEAR_EPS = 1e-9            # core GEOMETRY_EPSILON


# --- parsing ---------------------------------------------------------------------------------------------------------

@dataclass
class Residue:
    chain: str
    resseq: int
    icode: str
    resname: str
    letter: str
    ca: np.ndarray | None          # Cα by the altloc policy, or None
    n: np.ndarray | None
    c: np.ndarray | None
    read: bool                     # True if the adapter would read this residue (name in PROT_RES) and a Cα exists
    altloc_ca: bool                # the Cα has alternate locations in the file

    @property
    def key(self):
        return (self.resseq, self.icode)


@dataclass
class Parsed:
    text_sha256: str
    n_models: int
    residues: dict = field(default_factory=dict)     # chain -> [Residue] in file order
    dbref: list = field(default_factory=list)        # dicts
    seqadv: list = field(default_factory=list)       # dicts
    modres: dict = field(default_factory=dict)       # resname -> parent resname
    space_group: str = ""
    n_hetero_polymer_residues: int = 0


def _f(line, a, b):
    return line[a:b].strip()


def parse_pdb(text: str, altloc_policy: str = "blank_then_A") -> Parsed:
    """Fixed-column parse of the first model. altloc_policy: 'blank_then_A' or 'highest_occupancy' (ties: first listed)."""
    if altloc_policy not in ("blank_then_A", "highest_occupancy"):
        raise ValueError("unknown altloc policy")
    P = Parsed(hashlib.sha256(text.encode()).hexdigest(), 0)
    db1 = {}
    atoms: dict = {}          # (chain, resseq, icode) -> {"resname", "order", "recs": {name: [(alt, occ, xyz)]}}
    order: dict = {}
    in_model = True
    for line in text.splitlines():
        rec = line[:6]
        if rec == "MODEL ":
            P.n_models += 1
            if P.n_models > 1:
                in_model = False
            continue
        if rec == "ENDMDL":
            in_model = False                                  # only the first model is read
            continue
        if rec == "MODRES":
            P.modres[_f(line, 12, 15)] = _f(line, 24, 27)
        elif rec == "CRYST1":
            P.space_group = line[55:66].strip()
        elif rec == "DBREF ":
            P.dbref.append({"chain": line[12], "seq_begin": int(line[14:18]), "ins_begin": line[18].strip(),
                            "seq_end": int(line[20:24]), "ins_end": line[24].strip(), "database": _f(line, 26, 32),
                            "accession": _f(line, 33, 41), "id_code": _f(line, 42, 54), "db_begin": int(line[55:60]),
                            "db_end": int(line[62:67])})
        elif rec == "DBREF1":
            db1[line[12]] = {"chain": line[12], "seq_begin": int(line[14:18]), "ins_begin": line[18].strip(),
                             "seq_end": int(line[20:24]), "ins_end": line[24].strip(), "database": _f(line, 26, 32),
                             "id_code": _f(line, 47, 67)}
        elif rec == "DBREF2":
            d = db1.pop(line[12], None)
            if d:
                d.update({"accession": _f(line, 18, 40), "db_begin": int(line[45:55]), "db_end": int(line[57:67])})
                P.dbref.append(d)
        elif rec == "SEQADV":
            P.seqadv.append({"resname": _f(line, 12, 15), "chain": line[16], "seqnum": line[18:22].strip(),
                             "icode": line[22].strip(), "database": _f(line, 24, 28), "accession": _f(line, 29, 38),
                             "db_resname": _f(line, 39, 42), "conflict": line[49:70].strip()})
        elif rec in ("ATOM  ", "HETATM") and in_model:
            resname = line[17:20].strip()
            chain, icode = line[21], line[26].strip()
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            key = (chain, resseq, icode)
            if key not in atoms:
                atoms[key] = {"resname": resname, "recs": {}, "rec": rec}
                order.setdefault(chain, []).append(key)
            if atoms[key]["resname"] != resname:      # point-mutation microheterogeneity: keep the first name listed
                continue
            name = line[12:16].strip()
            if name in ("CA", "N", "C"):
                occ = float(line[54:60]) if line[54:60].strip() else 1.0
                atoms[key]["recs"].setdefault(name, []).append(
                    (line[16].strip(), occ, np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    for chain, keys in order.items():
        out = []
        for key in keys:
            a = atoms[key]
            resname = a["resname"]
            parent = P.modres.get(resname, resname)
            known = resname in PROT_RES or resname in THREE2ONE or resname in P.modres
            if not known:
                continue                                  # waters, ligands, ions
            letter = THREE2ONE.get(resname) or THREE2ONE.get(parent) or "X"
            sel = {n: _choose(a["recs"].get(n, []), altloc_policy) for n in ("CA", "N", "C")}
            ca_recs = a["recs"].get("CA", [])
            alt = len({r[0] for r in ca_recs}) > 1 or any(r[0] for r in ca_recs)
            out.append(Residue(chain, key[1], key[2], resname, letter, sel["CA"], sel["N"], sel["C"],
                               read=(resname in PROT_RES and sel["CA"] is not None), altloc_ca=alt))
        P.residues[chain] = out
    return P


def _choose(recs, policy):
    if not recs:
        return None
    if policy == "blank_then_A":
        for want in ("", "A"):
            for alt, _, xyz in recs:
                if alt == want:
                    return xyz
        return None
    best = max(recs, key=lambda r: r[1])                 # first listed among equal occupancy (max returns the first)
    return best[2]


# --- mapping to the UniProt reference ---------------------------------------------------------------------------------

def map_chain(P: Parsed, chain: str, ref_seq: str, accession: str) -> dict:
    """Align the chain's residues to the UniProt sequence and check the mapping against the deposited identifiers.

    Returns a dict with `index` ({(resseq, icode): 1-based UniProt index}), identity, coverage, mismatches (with their
    SEQADV documentation), DBREF agreement, and `problems` (empty if the mapping passes)."""
    res = P.residues.get(chain, [])
    problems = []
    if not res:
        return {"index": {}, "problems": [f"chain {chain} has no polymer residues"], "n_residues": 0}
    obs = "".join(r.letter for r in res)
    bonded = [False] * len(res)
    for j in range(1, len(res)):
        a, b = res[j - 1].c, res[j].n
        bonded[j] = a is not None and b is not None and float(np.linalg.norm(a - b)) <= PEPTIDE_BOND_MAX
    at_break = np.array([not b for b in bonded], dtype=bool)
    score, pairs = _align(ref_seq, obs, at_break)
    index = {res[j].key: i + 1 for i, j in pairs}
    ident = sum(ref_seq[i] == obs[j] for i, j in pairs)
    n_al = len(pairs)
    sa = {(d["chain"], d["seqnum"], d["icode"]): d for d in P.seqadv}
    mism = []
    for i, j in pairs:
        if ref_seq[i] != obs[j]:
            r = res[j]
            doc = sa.get((chain, str(r.resseq), r.icode))
            mism.append({"resseq": r.resseq, "icode": r.icode, "observed": r.resname, "uniprot": ref_seq[i],
                         "uniprot_index": i + 1, "seqadv": None if doc is None else doc["conflict"] or "listed"})
    segs = [d for d in P.dbref if d["chain"] == chain and d["accession"] == accession]
    if not segs:
        problems.append(f"no DBREF record of chain {chain} names accession {accession}")
    n_in, n_disagree = 0, []
    for j, r in enumerate(res):
        if r.key not in index or r.icode:
            continue
        for d in segs:
            if d["seq_begin"] <= r.resseq <= d["seq_end"]:
                n_in += 1
                exp = d["db_begin"] + (r.resseq - d["seq_begin"])
                if exp != index[r.key]:
                    n_disagree.append((r.resseq, exp, index[r.key]))
                break
    undocumented = [m for m in mism if m["seqadv"] is None]
    identity = ident / n_al if n_al else 0.0
    coverage = n_al / len(res)
    if identity < 0.98:
        problems.append(f"identity to the UniProt sequence over aligned residues is {identity:.3f} (< 0.98)")
    if coverage < 0.80:
        problems.append(f"only {coverage:.2f} of the chain's residues align to the UniProt sequence (< 0.80)")
    if undocumented:
        problems.append("residue-name differences without a SEQADV record: " + ", ".join(
            f"{m['observed']}{m['resseq']}{m['icode']}->{m['uniprot']}" for m in undocumented))
    if n_disagree:
        problems.append(f"{len(n_disagree)} aligned residues disagree with the DBREF numbering, first {n_disagree[0]}")
    if segs and n_in == 0:
        problems.append("no aligned residue lies inside the DBREF segment of this chain")
    return {"index": index, "n_residues": len(res), "n_aligned": n_al, "identity": identity, "coverage": coverage,
            "mismatches": mism, "dbref_segments": [{k: d[k] for k in ("seq_begin", "seq_end", "db_begin", "db_end")} for d in segs],
            "n_checked_against_dbref": n_in, "n_dbref_disagreements": len(n_disagree), "problems": problems,
            "n_unaligned": len(res) - n_al, "alignment_score": int(score)}


def shifted_identity(P: Parsed, chain: str, ref_seq: str, mapping: dict, shift: int = 1) -> float:
    """Residue-name identity when every mapped residue is moved by `shift` UniProt positions (the corruption control)."""
    res = {r.key: r for r in P.residues[chain]}
    ok = tot = 0
    for key, idx in mapping["index"].items():
        k = idx - 1 + shift
        if 0 <= k < len(ref_seq):
            tot += 1
            ok += res[key].letter == ref_seq[k]
    return ok / tot if tot else float("nan")


def correspondence(PA, ma, chainA, PB, mb, chainB) -> dict:
    """Pairs (A residue, B residue) of residues the adapter reads in both chains and that map to the same UniProt index."""
    ra = {r.key: r for r in PA.residues[chainA]}
    rb = {r.key: r for r in PB.residues[chainB]}
    ia = {idx: k for k, idx in ma["index"].items() if ra[k].read}
    ib = {idx: k for k, idx in mb["index"].items() if rb[k].read}
    shared = sorted(set(ia) & set(ib))
    pairs = [(ia[i], ib[i], i) for i in shared]
    span = range(max(min(ma["index"].values()), min(mb["index"].values())), min(max(ma["index"].values()), max(mb["index"].values())) + 1)
    unread = [i for i in span if i not in set(shared)]
    return {"pairs": pairs, "n_shared": len(shared), "n_span_not_shared": len(unread),
            "span": [span.start, span.stop - 1] if len(span) else None}


# --- independent numerical references ---------------------------------------------------------------------------------

def kabsch(P, Q):
    """Proper rotation R and translation so that R @ p + t best fits q (least squares); returns R, t."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def rmsd(P, Q):
    R, t = kabsch(P, Q)
    return float(np.sqrt(np.mean(np.sum((P @ R.T + t - Q) ** 2, axis=1))))


def kappa(p0, p1, p2):
    """Turning angle between consecutive bond vectors; None if a vector or the cross product is degenerate."""
    t0, t1 = p1 - p0, p2 - p1
    n0, n1 = np.linalg.norm(t0), np.linalg.norm(t1)
    if n0 < COLLINEAR_EPS or n1 < COLLINEAR_EPS:
        return None
    u0, u1 = t0 / n0, t1 / n1
    if np.linalg.norm(np.cross(u0, u1)) < COLLINEAR_EPS:
        return None
    return float(math.acos(max(-1.0, min(1.0, float(u0 @ u1)))))


def tau(p0, p1, p2, p3):
    """Signed angle between consecutive osculating planes about the middle bond (the core's definition)."""
    t = [p1 - p0, p2 - p1, p3 - p2]
    n = [np.linalg.norm(v) for v in t]
    if min(n) < COLLINEAR_EPS:
        return None
    u = [v / k for v, k in zip(t, n)]
    c0, c1 = np.cross(u[0], u[1]), np.cross(u[1], u[2])
    if np.linalg.norm(c0) < COLLINEAR_EPS or np.linalg.norm(c1) < COLLINEAR_EPS:
        return None
    b0, b1 = c0 / np.linalg.norm(c0), c1 / np.linalg.norm(c1)
    return float(math.atan2(float(np.cross(b0, b1) @ u[1]), float(b0 @ b1)))


def _next_icode(ic):
    return "A" if ic == "" else (chr(ord(ic) + 1) if len(ic) == 1 and "A" <= ic < "Z" else None)


def consecutive(a: Residue, b: Residue) -> bool:
    """The adapter's documented numbering adjacency: 52->52A->52B, 52->53, 52B->53; every other jump is a break."""
    if b.resseq == a.resseq:
        return b.icode == _next_icode(a.icode)
    return b.resseq == a.resseq + 1 and b.icode == ""


def independent_windows(resA: list, resB: list) -> dict:
    """Independent window metrics for aligned residue lists (same length; position j <-> resA[j], resB[j]).

    Per structure, a link between positions j and j+1 is continuous if the numbering is consecutive and the Cα-Cα distance
    is at most 4.2 Å; a run is a maximal chain of continuous links. The TopoFold features are computed only inside runs of
    at least 4 positions (a shorter run gives no κ or τ). A three-residue window j (positions j-1..j+1) and a
    four-residue window j (j-1..j+2) are valid if all their positions lie in one run of at least 4 positions in each
    structure. (Splitting of runs at collinear or coincident Cα triples is not modelled; such windows would appear as
    disagreements and be reported.)"""
    n = len(resA)
    XA = np.array([r.ca for r in resA])
    XB = np.array([r.ca for r in resB])

    def runs(res, X):
        link = [consecutive(res[j], res[j + 1]) and np.linalg.norm(X[j + 1] - X[j]) <= CA_CA_MAX for j in range(n - 1)]
        rid, length, cur = [0] * n, {}, 0
        for j in range(n):
            if j and not link[j - 1]:
                cur += 1
            rid[j] = cur
            length[cur] = length.get(cur, 0) + 1
        return link, rid, length

    linkA, ridA, lenA = runs(resA, XA)
    linkB, ridB, lenB = runs(resB, XB)

    def ok(lo, hi):
        return (len({ridA[k] for k in range(lo, hi + 1)}) == 1 and len({ridB[k] for k in range(lo, hi + 1)}) == 1
                and lenA[ridA[lo]] >= 4 and lenB[ridB[lo]] >= 4)

    link = [linkA[j] and linkB[j] for j in range(n - 1)]
    out = {"n": n, "valid3": [False] * n, "valid4": [False] * n, "d_kappa": [None] * n, "d_i_i2": [None] * n,
           "rmsd3": [None] * n, "d_tau": [None] * n, "d_i_i3": [None] * n, "rmsd4": [None] * n, "links_broken": [j for j in range(n - 1) if not link[j]]}
    for j in range(1, n - 1):
        if ok(j - 1, j + 1):
            ka, kb = kappa(*XA[j - 1:j + 2]), kappa(*XB[j - 1:j + 2])
            if ka is None or kb is None:
                continue
            out["valid3"][j] = True
            out["d_kappa"][j] = abs(kb - ka)
            out["d_i_i2"][j] = abs(float(np.linalg.norm(XB[j + 1] - XB[j - 1]) - np.linalg.norm(XA[j + 1] - XA[j - 1])))
            out["rmsd3"][j] = rmsd(XA[j - 1:j + 2], XB[j - 1:j + 2])
    for j in range(1, n - 2):
        if ok(j - 1, j + 2):
            ta, tb = tau(*XA[j - 1:j + 3]), tau(*XB[j - 1:j + 3])
            if ta is None or tb is None:
                continue
            d = tb - ta
            out["valid4"][j] = True
            out["d_tau"][j] = abs(math.atan2(math.sin(d), math.cos(d)))
            out["d_i_i3"][j] = abs(float(np.linalg.norm(XB[j + 2] - XB[j - 1]) - np.linalg.norm(XA[j + 2] - XA[j - 1])))
            out["rmsd4"][j] = rmsd(XA[j - 1:j + 3], XB[j - 1:j + 3])
    out["global_rmsd"] = rmsd(XA, XB)
    return out


def trimmed_displacement(XA, XB, trim=0.2):
    """Cα displacement after a trimmed fit: fit all, drop the floor(trim*n) most displaced (ties: lower position first),
    refit once; the displacement of every position under the refit. Also returns the all-position fit displacement."""
    R, t = kabsch(XA, XB)
    d0 = np.linalg.norm(XA @ R.T + t - XB, axis=1)
    k = int(math.floor(trim * len(XA)))
    order = sorted(range(len(XA)), key=lambda i: (-d0[i], i))
    keep = sorted(order[k:])
    R2, t2 = kabsch(XA[keep], XB[keep])
    return np.linalg.norm(XA @ R2.T + t2 - XB, axis=1), d0


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    rx = _rank(x)
    ry = _rank(y)
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def _rank(v):
    order = np.argsort(v, kind="stable")
    r = np.empty(len(v))
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        r[order[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return r


def top_k(values: list, k: int = 10) -> list:
    """Positions of the k largest values (ties: lower position first); None values never rank."""
    idx = [i for i, v in enumerate(values) if v is not None]
    return sorted(sorted(idx, key=lambda i: (-values[i], i))[:k])


def classify_pair(win: dict, k: int = 10, min_windows: int = 50) -> dict:
    """Redundant / complementary / insufficient, by the protocol rule (per window size).

    O = how many of the k top windows of the geometry metric (d_kappa or d_tau) are also among the k top windows of either
    equal-information baseline (Cα-distance difference or local RMSD). O >= 8: redundant with the baselines. O <= 7: at least
    3 of the geometry metric's top windows are not top windows of any baseline: complementary descriptive information.
    Fewer than `min_windows` valid windows: insufficient evidence."""
    out = {}
    for size, geo, bases in ((3, "d_kappa", ("d_i_i2", "rmsd3")), (4, "d_tau", ("d_i_i3", "rmsd4"))):
        valid = [i for i, v in enumerate(win[geo]) if v is not None]
        if len(valid) < min_windows:
            out[str(size)] = {"class": "insufficient_evidence", "n_valid": len(valid)}
            continue
        g = top_k(win[geo], k)
        union = set()
        for b in bases:
            union |= set(top_k(win[b], k))
        o = len(set(g) & union)
        rho = {b: spearman([win[geo][i] for i in valid], [win[b][i] for i in valid]) for b in bases}
        out[str(size)] = {"class": "redundant" if o >= 8 else "complementary", "overlap_O": o, "k": k, "n_valid": len(valid),
                          "spearman": rho, "top_geometry": g, "top_baselines": {b: top_k(win[b], k) for b in bases}}
    return out
