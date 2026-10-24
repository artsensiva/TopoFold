"""ABL-MYR-1b comparability check (post-hoc; see check_spec.json, hash-frozen before this code reads any study file).

    <python with MDAnalysis> benchmarks/audit/abl_myr_1b/run_comp_check.py freeze
    <python with MDAnalysis> benchmarks/audit/abl_myr_1b/run_comp_check.py run
    <python with MDAnalysis> benchmarks/audit/abl_myr_1b/run_comp_check.py report

Reads only run inputs and parameter files (md.tpr, nvt/npt/production mdp, index.ndx, topol.top, system_em.gro). Opens no
trajectory. GROMACS is not installed, so md.tpr coordinates read by MDAnalysis are validated independently: atom order
against the text-parsed system_em.gro, protein bond lengths, and a Cα superposition on system_em.gro (NumPy Kabsch).
The archive helpers and text parsers are those of the frozen ABL-MYR-1 runner, imported unchanged.
"""

from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import re
import sys
import time
import zipfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ORIG = ROOT / "benchmarks/audit/abl_myr_1"
AMEND = ROOT / "benchmarks/audit/abl_myr_1a"
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ORIG))
import run_abl_myr as R0  # noqa: E402  (frozen runner, unchanged)

SPEC = HERE / "check_spec.json"
FREEZE = HERE / "results/ABL-MYR-1b_check_freeze.json"
RESULTS = "ABL-MYR-1b_check_results.json"
CA_NPZ = "ABL-MYR-1b_production_start_ca.npz"
SELF = Path(__file__).resolve()
TEST = "tests/test_abl_myr_1b_check.py"
PRESERVED = [f"benchmarks/audit/abl_myr_1/{p}" for p in (
    "run_spec.json", "run_abl_myr.py", "g1_member_audit.py", "results/ABL-MYR-1_freeze.json",
    "results/ABL-MYR-1_code_freeze.json", "results/ABL-MYR-1_blocked.json", "results/ABL-MYR-1_g1_member_audit.json")] + [
    f"benchmarks/audit/abl_myr_1a/{p}" for p in (
        "audit_spec.json", "run_init_audit.py", "results/ABL-MYR-1a_audit_freeze.json",
        "results/ABL-MYR-1a_audit_results.json", "results/ABL-MYR-1a_audit_report.md")]
SETS = {"Set1": list(range(11, 21)), "Set2": list(range(21, 31))}
CONDITIONS = {"Myr": "kinase_kinked_Myr", "noMyr": "kinase_kinked_noMyr"}
MDP_FILES = ("nvt.mdp", "npt.mdp", "production_md.mdp")
BOND_MIN, BOND_MAX, CA_RMSD_MAX, IDENT_RMSD = 0.9, 2.3, 1.0, 0.01
SOL_REL_MAX, CHARGE_DIFF_MAX, REF_T = 0.01, 0.5, "300"
SOLUTE_EXTRA = {"ATP", "MG"}
SOLVENT = {"SOL", "NA", "CL"}


class CheckStop(Exception):
    def __init__(self, code, cause):
        self.code, self.cause = code, cause
        super().__init__(f"{code}: {cause}")


# --- pure helpers ------------------------------------------------------------------------------------------------------

def norm_value(v):
    return " ".join(str(v).split()) if v is not None else None


def strip_myr(value: str) -> str:
    out = []
    for tok in value.split():
        parts = [p for p in tok.split("_") if p.upper() != "MYRP"]
        out.append("_".join(parts))
    return " ".join(out)


def compare_mdp(myr_text: str, nomyr_text: str) -> dict:
    va = {k: norm_value(v) for k, v in R0.mdp_values(myr_text).items()}
    vb = {k: norm_value(v) for k, v in R0.mdp_values(nomyr_text).items()}
    diff = {}
    for k in sorted(set(va) | set(vb)):
        a, b = va.get(k), vb.get(k)
        if a == b:
            continue
        cls = "inherent_group_name" if (a is not None and b is not None and "MYRP" in a.upper()
                                        and strip_myr(a) == b) else "potential_confound"
        diff[k] = {"Myr": a, "noMyr": b, "class": cls}
    return {"differences": diff, "potential_confounds": sorted(k for k, d in diff.items() if d["class"] == "potential_confound"),
            "values_Myr": va, "values_noMyr": vb}


def ref_t_ok(v) -> bool:
    """Every ref_t entry equals 300 K numerically."""
    try:
        return v is not None and len(v.split()) > 0 and all(float(t) == float(REF_T) for t in v.split())
    except ValueError:
        return False


def parse_ndx(text: str) -> dict:
    groups, cur = {}, None
    for ln in text.splitlines():
        s = ln.strip()
        m = re.match(r"^\[\s*(.+?)\s*\]$", s)
        if m:
            cur = m.group(1)
            if cur in groups:
                raise ValueError(f"duplicate index group {cur!r}")
            groups[cur] = []
        elif s and cur is not None:
            groups[cur] += [int(x) - 1 for x in s.split()]
        elif s:
            raise ValueError("index entries before the first group")
    return {k: np.array(v, dtype=np.int64) for k, v in groups.items()}


def tc_group_check(groups: dict, tc_names: list, atom_resnames, atom_resindex, protein_resnames: set, condition: str) -> dict:
    n = len(atom_resnames)
    out = {"tc_grps": tc_names, "groups": {}, "ok": True, "problems": []}
    seen = np.zeros(n, dtype=np.int64)
    atom_resindex = np.asarray(atom_resindex)
    ri_u, first = np.unique(atom_resindex, return_index=True)
    name_by_ri = {int(r): str(atom_resnames[i]) for r, i in zip(ri_u, first)}
    for g in tc_names:
        if g not in groups:
            out["ok"] = False
            out["problems"].append(f"tc group {g!r} not in index.ndx")
            continue
        idx = groups[g]
        if idx.size and (idx.min() < 0 or idx.max() >= n):
            out["ok"] = False
            out["problems"].append(f"tc group {g!r} has indices outside the system")
            continue
        seen[idx] += 1
        res = collections.Counter(name_by_ri[int(r)] for r in np.unique(atom_resindex[idx]))
        kinds = set(res)
        is_solvent = bool(kinds) and kinds <= SOLVENT
        allowed = SOLVENT if is_solvent else (protein_resnames | SOLUTE_EXTRA | ({"MYRP"} if condition == "Myr" else set()))
        bad = sorted(kinds - allowed)
        out["groups"][g] = {"atoms": int(idx.size), "residue_counts": dict(sorted(res.items())), "role": "solvent" if is_solvent
                            else "solute", "outside_rule": bad}
        if bad:
            out["ok"] = False
            out["problems"].append(f"group {g!r} contains {bad}")
    if out["ok"]:
        if (seen != 1).any():
            out["ok"] = False
            out["problems"].append(f"tc groups do not partition the system: {int((seen == 0).sum())} atoms in no group, "
                                   f"{int((seen > 1).sum())} in several")
        roles = [v["role"] for v in out["groups"].values()]
        if sorted(roles) != ["solute", "solvent"]:
            out["ok"] = False
            out["problems"].append(f"expected one solute and one solvent group, got {roles}")
        myr_in = any("MYRP" in v["residue_counts"] for v in out["groups"].values())
        if myr_in != (condition == "Myr"):
            out["ok"] = False
            out["problems"].append(f"MYRP {'absent' if condition == 'Myr' else 'present'} in the {condition} tc groups")
    return out


def composition_check(counts_myr: dict, counts_nomyr: dict, q_myr: float, q_nomyr: float) -> dict:
    diffs, inherent, confounds = {}, [], []
    for k in sorted(set(counts_myr) | set(counts_nomyr)):
        a, b = counts_myr.get(k, 0), counts_nomyr.get(k, 0)
        if a == b:
            continue
        diffs[k] = [a, b]
        if k == "MYRP" and b == 0 and a > 0:
            inherent.append("MYRP present only with Myr")
        elif k == "SOL" and b and abs(a - b) / b <= SOL_REL_MAX:
            inherent.append(f"SOL differs by {a - b} ({abs(a - b) / b:.4%} of noMyr)")
        else:
            confounds.append(f"{k}: {a} vs {b}")
    dq = abs(q_myr - q_nomyr)
    if dq > CHARGE_DIFF_MAX:
        confounds.append(f"net charge differs by {dq:.3f} e")
    return {"differences": diffs, "net_charge_e": {"Myr": round(q_myr, 6), "noMyr": round(q_nomyr, 6)},
            "net_charge_difference_e": round(dq, 6), "inherent": inherent, "potential_confounds": confounds,
            "ok": not confounds}


def kabsch_rmsd(P, Q) -> float:
    """RMSD after optimal proper superposition of Q onto P (NumPy SVD; independent of the MDAnalysis reader)."""
    P, Q = np.asarray(P, float), np.asarray(Q, float)
    p, q = P - P.mean(0), Q - Q.mean(0)
    u, _, vt = np.linalg.svd(q.T @ p)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1.0, 1.0, d]) @ vt
    return float(np.sqrt(np.mean(np.sum((q @ r - p) ** 2, axis=1))))


def radius_of_gyration(X) -> float:
    X = np.asarray(X, float)
    return float(np.sqrt(np.mean(np.sum((X - X.mean(0)) ** 2, axis=1))))


def whole_chain(xyz, dims):
    """Make a chain trace continuous by sequential minimum image (triclinic boxes supported)."""
    from MDAnalysis.lib.distances import minimize_vectors
    X = np.array(xyz, dtype=np.float64)
    for i in range(1, len(X)):
        d = minimize_vectors((X[i] - X[i - 1])[None].astype(np.float32), np.asarray(dims, np.float32))[0]
        X[i] = X[i - 1] + d
    return X


def bond_stats(pos, bonds, dims) -> dict:
    from MDAnalysis.lib.distances import calc_bonds
    if len(bonds) == 0:
        return {"n_bonds": 0, "min_A": None, "max_A": None, "outside": 0, "ok": False}
    d = calc_bonds(pos[bonds[:, 0]].astype(np.float64), pos[bonds[:, 1]].astype(np.float64), box=np.asarray(dims, np.float64))
    out = int(((d < BOND_MIN) | (d > BOND_MAX)).sum())
    return {"n_bonds": int(len(d)), "min_A": round(float(d.min()), 4), "max_A": round(float(d.max()), 4),
            "outside": out, "ok": out == 0}


def gro_box_dims(path) -> np.ndarray:
    from MDAnalysis.lib.mdamath import triclinic_box
    v = [float(x) * 10.0 for x in Path(path).read_text().splitlines()[-1].split()]
    if len(v) == 3:
        return np.array([v[0], v[1], v[2], 90.0, 90.0, 90.0])
    if len(v) != 9:
        raise ValueError(f"{path}: box line has {len(v)} values")
    a, b, c = np.array([v[0], v[3], v[4]]), np.array([v[5], v[1], v[6]]), np.array([v[7], v[8], v[2]])
    return np.asarray(triclinic_box(a, b, c), float)


def groups_distinct(n: int, identical_pairs) -> int:
    parent = list(range(n))

    def f(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i, j in identical_pairs:
        parent[f(i)] = f(j)
    return len({f(i) for i in range(n)})


def summarize(values) -> dict:
    v = np.asarray(values, float)
    return {"n": int(v.size), "min": float(v.min()), "median": float(np.median(v)), "max": float(v.max())} if v.size else {"n": 0}


def decide(k: dict) -> dict:
    names = ["K1_usable_coordinates", "K2_thermostat_comparable", "K3_tc_groups", "K4_composition", "K5_replicates_not_copies"]
    failed = [n for n in names if not k[n]]
    return {"criteria": {n: bool(k[n]) for n in names}, "failed": failed,
            "outcome": "comparable_for_1b" if not failed else "not_comparable"}


# --- MDAnalysis loader ---------------------------------------------------------------------------------------------------

def load_tpr(path) -> dict:
    import warnings

    import MDAnalysis as mda
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = mda.Universe(str(path))
        pos = np.asarray(u.atoms.positions, dtype=np.float64)
        dims = None if u.dimensions is None else np.asarray(u.dimensions, dtype=np.float64)
        try:
            vel = np.asarray(u.atoms.velocities, dtype=np.float64)
        except Exception:
            vel = None
        prot = u.select_atoms("protein")
        pidx = set(prot.indices.tolist())
        b = np.asarray(u.atoms.bonds.indices, dtype=np.int64) if hasattr(u.atoms, "bonds") else np.zeros((0, 2), np.int64)
        pb = b[[(i in pidx and j in pidx) for i, j in b]] if len(b) else b
        names = np.array([str(x) for x in u.atoms.names])
        resn = np.array([str(x) for x in u.atoms.resnames])
        ca_idx = np.flatnonzero(np.isin(resn, sorted(R0.AMINO)) & (names == "CA"))   # the same residue set as the gro parser
        heavy = u.select_atoms("protein and not name H*")
        heavy_whole = heavy.unwrap(compound="fragments", reference="com", inplace=False) if dims is not None else None
        counts = collections.Counter(str(r) for r in u.residues.resnames)
        prot_res = set(str(x) for x in prot.residues.resnames)
        return {"n_atoms": int(u.atoms.n_atoms), "names": names.tolist(), "resnames": resn,
                "resindex": np.asarray(u.atoms.resindices), "protein_resnames": prot_res,
                "protein_residues_outside_gro_set": sorted(prot_res - set(R0.AMINO)), "positions": pos, "dims": dims,
                "velocities": vel, "protein_bonds": pb, "ca_positions": pos[ca_idx], "ca_resnames": resn[ca_idx].tolist(),
                "heavy_whole": None if heavy_whole is None else np.asarray(heavy_whole, float),
                "protein_vel": None if vel is None else vel[prot.indices], "residue_counts": dict(sorted(counts.items())),
                "net_charge": float(np.sum(u.atoms.charges))}


# --- check -------------------------------------------------------------------------------------------------------------

def member_plan() -> dict:
    plan = {}
    for c, model in CONDITIONS.items():
        m = {}
        for s, reps in SETS.items():
            for f in MDP_FILES + ("index.ndx", "topol.top", "system_em.gro"):
                m[f"{f}:{s}"] = f"{s}/{f}"
            for r in reps:
                m[f"tpr:{s}:{r}"] = f"{s}/{r}/md.tpr"
        plan[f"{model}.zip"] = m
    return plan


def validate_run(rec: dict, em_names: list, em_ca_A, em_ca_resnames: list, em_dims) -> dict:
    out = {"positions_finite": bool(np.isfinite(rec["positions"]).all()), "box_present": rec["dims"] is not None
           and bool(np.all(np.asarray(rec["dims"])[:3] > 0))}
    out["atom_order_matches_em"] = rec["names"] == em_names
    out["bonds"] = bond_stats(rec["positions"], rec["protein_bonds"], rec["dims"]) if out["box_present"] else {"ok": False}
    ca_ok = rec["ca_resnames"] == em_ca_resnames and len(em_ca_A) == len(rec["ca_positions"]) and out["box_present"]
    if ca_ok:
        mine = whole_chain(rec["ca_positions"], rec["dims"])
        ref = whole_chain(em_ca_A, em_dims)
        out["ca_rmsd_to_em_A"] = round(kabsch_rmsd(ref, mine), 6)
        out["ca_whole"] = mine
    else:
        out["ca_rmsd_to_em_A"] = None
    out["ok"] = bool(out["positions_finite"] and out["box_present"] and out["atom_order_matches_em"] and out["bonds"]["ok"]
                     and out["ca_rmsd_to_em_A"] is not None and out["ca_rmsd_to_em_A"] <= CA_RMSD_MAX)
    return out


def run_check(orig_spec: dict, data_dir: Path, extract_dir: Path, spec: dict, loader=load_tpr, log=print) -> tuple:
    t0, cpu0 = time.time(), R0.cpu_seconds()
    budget = spec["budget"]
    res = {"id": spec["id"], "post_hoc": True}
    ca_store = {}
    try:
        for a in CONDITIONS.values():
            want, p = orig_spec["inputs"][f"{a}.zip"], Path(data_dir) / f"{a}.zip"
            if not p.is_file() or p.stat().st_size != want["bytes"] or R0.sha256_file(p) != want["sha256"]:
                raise CheckStop("S1", f"{a}.zip is missing or differs from the frozen size or SHA-256")
        paths, files, total = {}, {}, 0
        for arch, keys in member_plan().items():
            with zipfile.ZipFile(Path(data_dir) / arch) as zf:
                names = zf.namelist()
                chosen = {}
                for k, suf in keys.items():
                    try:
                        chosen[k] = R0.resolve_member(names, suf)
                    except R0.GateError as e:
                        raise CheckStop("S2", e.message)
                total += sum(zf.getinfo(x).file_size for x in chosen.values())
                if total > budget["extracted_bytes_max"]:
                    raise CheckStop("S4", f"planned extraction {total} B exceeds the budget")
                for k, member in chosen.items():
                    e = R0.safe_extract(zf, member, Path(extract_dir) / Path(arch).stem)
                    paths[(arch, k)] = e["path"]
                    files[f"{arch}:{k}"] = {"member": member, "bytes": e["bytes"], "sha256": e["sha256"]}
        res["extracted"] = {"files": len(files), "bytes": sum(f["bytes"] for f in files.values()), "sha256": files}
        P = lambda c, k: paths[(f"{CONDITIONS[c]}.zip", k)]  # noqa: E731
        try:
            mdp = {s: {f: compare_mdp(Path(P("Myr", f"{f}:{s}")).read_text(), Path(P("noMyr", f"{f}:{s}")).read_text())
                       for f in MDP_FILES} for s in SETS}
            ndx = {s: {c: parse_ndx(Path(P(c, f"index.ndx:{s}")).read_text()) for c in CONDITIONS} for s in SETS}
            em = {}
            for s in SETS:
                for c in CONDITIONS:
                    atoms = R0.parse_gro(P(c, f"system_em.gro:{s}"))
                    prot = R0.protein_atoms(atoms)
                    em[(s, c)] = {"names": [a[1] for a in atoms], "ca": np.array([a[2] for a in prot if a[1] == "CA"]) * 10.0,
                                  "ca_resnames": [a[0] for a in prot if a[1] == "CA"], "dims": gro_box_dims(P(c, f"system_em.gro:{s}"))}
        except (ValueError, IndexError) as e:
            raise CheckStop("S3", f"cannot parse a parameter or structure file: {type(e).__name__}: {e}")
        # R4 thermostat and coupling
        r4, k2 = {}, True
        for s in SETS:
            r4[s] = {}
            for f in MDP_FILES:
                m = mdp[s][f]
                va, vb = m["values_Myr"], m["values_noMyr"]
                ngr = [len((va.get("tc_grps") or "").split()), len((vb.get("tc_grps") or "").split())]
                rt = [ref_t_ok(va.get("ref_t")), ref_t_ok(vb.get("ref_t"))]
                stated = {c: {k: v.get(k) for k in ("integrator", "tcoupl", "tc_grps", "tau_t", "ref_t", "pcoupl", "tau_p",
                                                    "ref_p", "compressibility", "define", "gen_vel", "gen_seed", "gen_temp",
                                                    "continuation", "dt", "nsteps", "constraints", "refcoord_scaling")}
                          for c, v in (("Myr", va), ("noMyr", vb))}
                r4[s][f] = {"differences": m["differences"], "potential_confounds": m["potential_confounds"],
                            "ref_t_all_300": rt, "n_tc_groups": ngr, "stated": stated}
                k2 &= not m["potential_confounds"] and all(rt) and ngr[0] == ngr[1] and ngr[0] > 0
        # R1, R2: run inputs
        r1, recs = {}, {}
        for s, reps in SETS.items():
            for c in CONDITIONS:
                ref = em[(s, c)]
                for r in reps:
                    key = f"{s}:{c}:{r}"
                    try:
                        rec = loader(P(c, f"tpr:{s}:{r}"))
                    except Exception as e:
                        r1[key] = {"readable": False, "error": f"{type(e).__name__}: {str(e)[:200]}", "ok": False}
                        continue
                    v = validate_run(rec, ref["names"], ref["ca"], ref["ca_resnames"], ref["dims"])
                    ca_store[key] = v.pop("ca_whole", None)
                    v.update(readable=True, n_atoms=rec["n_atoms"], velocities_present=rec["velocities"] is not None,
                             protein_residues_outside_gro_set=rec.get("protein_residues_outside_gro_set", []),
                             box_A=None if rec["dims"] is None else [round(float(x), 4) for x in rec["dims"]],
                             ca_rg_A=None if ca_store[key] is None else round(radius_of_gyration(ca_store[key]), 6))
                    if rec["dims"] is not None:
                        from MDAnalysis.lib.mdamath import box_volume
                        v["box_volume_A3"] = round(float(box_volume(rec["dims"])), 3)
                    r1[key] = v
                    recs[key] = rec
                if R0.cpu_seconds() - cpu0 > budget["cpu_seconds_max"]:
                    raise CheckStop("S4", "CPU budget exceeded")
        k1 = len(r1) == 40 and all(v.get("ok") for v in r1.values())
        res["R1_usable_coordinates"] = {"per_run": r1, "runs_ok": sum(bool(v.get("ok")) for v in r1.values()), "runs": len(r1)}
        # R2 pairwise geometry
        def pair_stats(keys_a, keys_b, same):
            ca, hv, ident = [], [], []
            for i, a in enumerate(keys_a):
                for j, b in enumerate(keys_b):
                    if same and j <= i:
                        continue
                    if ca_store.get(a) is None or ca_store.get(b) is None:
                        continue
                    x = kabsch_rmsd(ca_store[a], ca_store[b])
                    ca.append(x)
                    ha, hb = recs[a]["heavy_whole"], recs[b]["heavy_whole"]
                    if ha is not None and hb is not None and len(ha) == len(hb):
                        hv.append(kabsch_rmsd(ha, hb))
                    if x < IDENT_RMSD:
                        ident.append((i, j))
            return {"ca_rmsd_A": summarize(ca), "heavy_rmsd_A": summarize(hv), "pairs_below_0.01A": len(ident)}, ident
        r2, k5 = {}, True
        for s, reps in SETS.items():
            r2[s] = {}
            keys = {c: [f"{s}:{c}:{r}" for r in reps if f"{s}:{c}:{r}" in recs] for c in CONDITIONS}
            for c in CONDITIONS:
                st, ident = pair_stats(keys[c], keys[c], True)
                vel = [recs[k]["protein_vel"] for k in keys[c]]
                vh = [None if v is None else hashlib.sha256(np.round(v, 4).tobytes()).hexdigest() for v in vel]
                st["distinct_starts"] = groups_distinct(len(keys[c]), ident)
                st["velocities_present"] = all(h is not None for h in vh) and bool(vh)
                st["distinct_velocity_sets"] = len(set(vh)) if st["velocities_present"] else None
                r2[s][f"within_{c}"] = st
                k5 &= bool(keys[c]) and (st["distinct_starts"] >= 2 or (st["velocities_present"] and st["distinct_velocity_sets"] >= 2))
            r2[s]["between_Myr_noMyr"] = pair_stats(keys["Myr"], keys["noMyr"], False)[0]
        r2["cross_set"] = {c: pair_stats([f"Set1:{c}:{r}" for r in SETS["Set1"]], [f"Set2:{c}:{r}" for r in SETS["Set2"]],
                                         False)[0] for c in CONDITIONS}
        res["R2_production_start_geometry"] = r2
        # R3 tc groups and R5 composition (from the first readable run input of each set and condition)
        r3, r5, k3, k4 = {}, {}, True, True
        for s, reps in SETS.items():
            r3[s], comp = {}, {}
            for c in CONDITIONS:
                ks = [f"{s}:{c}:{r}" for r in reps if f"{s}:{c}:{r}" in recs]
                if not ks:
                    k3 = k4 = False
                    continue
                rec = recs[ks[0]]
                tc = (mdp[s]["production_md.mdp"]["values_" + c].get("tc_grps") or "").split()
                r3[s][c] = tc_group_check(ndx[s][c], tc, rec["resnames"], rec["resindex"], rec["protein_resnames"], c)
                k3 &= r3[s][c]["ok"]
                distinct = {json.dumps(recs[k]["residue_counts"], sort_keys=True) for k in ks}
                charges = sorted({round(recs[k]["net_charge"], 6) for k in ks})
                comp[c] = {"residue_counts": rec["residue_counts"], "net_charge": rec["net_charge"],
                           "distinct_compositions_across_runs": len(distinct), "net_charges_across_runs": charges,
                           "n_atoms": rec["n_atoms"]}
                k4 &= len(distinct) == 1
            if len(comp) == 2:
                cc = composition_check(comp["Myr"]["residue_counts"], comp["noMyr"]["residue_counts"],
                                       comp["Myr"]["net_charge"], comp["noMyr"]["net_charge"])
                r5[s] = {"per_condition": comp, "comparison": cc}
                k4 &= cc["ok"]
        res["R3_tc_groups"], res["R4_thermostat_and_coupling"], res["R5_composition"] = r3, r4, r5
        # R6 classification
        cls = {"inherent_to_Myr_model": [], "potential_confound": [], "not_separable_from_Myr_presence": []}
        for s in SETS:
            for f in MDP_FILES:
                for k, d in r4[s][f]["differences"].items():
                    (cls["inherent_to_Myr_model"] if d["class"] == "inherent_group_name" else cls["potential_confound"]).append(
                        f"{s} {f} {k}: {d['Myr']} | {d['noMyr']}")
            if s in r5:
                cls["inherent_to_Myr_model"] += [f"{s} {x}" for x in r5[s]["comparison"]["inherent"]]
                cls["potential_confound"] += [f"{s} {x}" for x in r5[s]["comparison"]["potential_confounds"]]
            for c in CONDITIONS:
                if c in r3.get(s, {}):
                    cls["potential_confound"] += [f"{s} {c} tc groups: {p}" for p in r3[s][c]["problems"]]
            cls["not_separable_from_Myr_presence"].append(f"{s}: production-start geometry differences between Myr and noMyr "
                                                          f"(see R2 between_Myr_noMyr)")
        res["R6_classification"] = cls
        res["decision"] = decide({"K1_usable_coordinates": k1, "K2_thermostat_comparable": k2, "K3_tc_groups": k3,
                                  "K4_composition": k4, "K5_replicates_not_copies": k5})
        res["outcome"] = res["decision"]["outcome"]
        if R0.cpu_seconds() - cpu0 > budget["cpu_seconds_max"]:
            raise CheckStop("S4", "CPU budget exceeded")
    except CheckStop as e:
        res["outcome"] = "check_incomplete"
        res["stop"] = {"code": e.code, "cause": e.cause}
    res["resources"] = {"cpu_s": R0.cpu_seconds() - cpu0, "wall_s": time.time() - t0, "download_bytes": 0}
    return res, ca_store


def _group_text(g: str, x: dict) -> str:
    rc = x["residue_counts"]
    if x["role"] == "solute":
        lig = {k: n for k, n in rc.items() if k in ("ATP", "MG", "MYRP")}
        return f"{g}: {x['atoms']} atoms, {sum(rc.values()) - sum(lig.values())} protein residues + {lig}"
    return f"{g}: {x['atoms']} atoms, {rc}"


def render_report(res: dict) -> str:
    L = ["# ABL-MYR-1b comparability check (generated from ABL-MYR-1b_check_results.json)", "",
         "Post-hoc check of run inputs and parameter files; no trajectory opened. ABL-MYR-1 (G1 stop) and ABL-MYR-1a (V4) remain "
         "valid. Not evidence about Ayaz WT/M472I, HDX, drug response or biological function; no causal Myr claim.", "",
         f"**Outcome: {res['outcome']}**" + (f" (stop {res['stop']['code']}: {res['stop']['cause']})" if "stop" in res else ""), ""]
    if "decision" not in res:
        return "\n".join(L) + "\n"
    L += ["| Criterion | Pass |", "|---|---|"] + [f"| {k} | {v} |" for k, v in res["decision"]["criteria"].items()] + [""]
    r1 = res["R1_usable_coordinates"]
    rms = [v["ca_rmsd_to_em_A"] for v in r1["per_run"].values() if v.get("ca_rmsd_to_em_A") is not None]
    bmax = [v["bonds"]["max_A"] for v in r1["per_run"].values() if v.get("bonds", {}).get("max_A") is not None]
    bmin = [v["bonds"]["min_A"] for v in r1["per_run"].values() if v.get("bonds", {}).get("min_A") is not None]
    L += [f"R1: {r1['runs_ok']} of {r1['runs']} run inputs usable and validated. Atom order equal to system_em.gro in "
          f"{sum(bool(v.get('atom_order_matches_em')) for v in r1['per_run'].values())} runs; protein bond lengths "
          f"{min(bmin) if bmin else None}–{max(bmax) if bmax else None} Å; Cα RMSD to system_em.gro "
          f"{min(rms) if rms else None}–{max(rms) if rms else None} Å (limit 1.0 Å); velocities present in "
          f"{sum(bool(v.get('velocities_present')) for v in r1['per_run'].values())} runs.", ""]
    L += ["R2 production-start Cα RMSD after superposition (Å; min / median / max, pairs):", ""]
    for s in ("Set1", "Set2"):
        for k, v in res["R2_production_start_geometry"][s].items():
            c = v["ca_rmsd_A"]
            extra = (f"; distinct starts {v['distinct_starts']}, distinct velocity sets {v['distinct_velocity_sets']}"
                     if "distinct_starts" in v else "")
            L.append(f"- {s} {k}: {c.get('min', 0):.3f} / {c.get('median', 0):.3f} / {c.get('max', 0):.3f} ({c['n']} pairs); "
                     f"heavy-atom median {v['heavy_rmsd_A'].get('median', float('nan')):.3f}{extra}")
    for c, v in res["R2_production_start_geometry"]["cross_set"].items():
        L.append(f"- Set1 vs Set2 {c}: Cα median {v['ca_rmsd_A'].get('median', float('nan')):.3f} ({v['ca_rmsd_A']['n']} pairs)")
    L += ["", "R3 tc groups:", ""]
    for s, d in res["R3_tc_groups"].items():
        for c, v in d.items():
            L.append(f"- {s} {c}: {v['tc_grps']} ok={v['ok']}; " + "; ".join(_group_text(g, x) for g, x in v["groups"].items())
                     + (f"; problems {v['problems']}" if v["problems"] else ""))
    L += ["", "R4 thermostat and coupling (differences between conditions):", ""]
    for s, d in res["R4_thermostat_and_coupling"].items():
        for f, v in d.items():
            st = v["stated"]["Myr"]
            L.append(f"- {s} {f}: tcoupl {st['tcoupl']}, ref_t {st['ref_t']}, tau_t {st['tau_t']}, pcoupl {st['pcoupl']}, "
                     f"gen_vel {st['gen_vel']}, gen_seed {st['gen_seed']}, define {st['define']}; differences "
                     + (", ".join(f"{k} ({x['class']})" for k, x in v["differences"].items()) or "none"))
    L += ["", "R5 composition:", ""]
    for s, d in res["R5_composition"].items():
        cc = d["comparison"]
        L.append(f"- {s}: differences {cc['differences']}; net charge {cc['net_charge_e']}; inherent {cc['inherent']}; "
                 f"potential confounds {cc['potential_confounds']}")
    L += ["", "R6 classification:", ""] + [f"- {k}: {v}" for k, v in res["R6_classification"].items()] + [""]
    return "\n".join(L) + "\n"


# --- freeze / run / report ---------------------------------------------------------------------------------------------

def _h(p):
    return R0.sha256_file(ROOT / p)


def cmd_freeze(args) -> int:
    if FREEZE.exists():
        print("freeze exists")
        return 1
    if (ORIG / "data_dryad/extracted_1b").exists():
        raise SystemExit("STOP: extracted_1b exists; the check freeze must precede extraction")
    rec = {"id": json.loads(SPEC.read_text())["id"], "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "spec_sha256": R0.sha256_file(SPEC), "code": {"path": str(SELF.relative_to(ROOT)), "sha256": R0.sha256_file(SELF)},
           "test": {TEST: _h(TEST)}, "preserved_sha256": {p: _h(p) for p in PRESERVED},
           "frozen_runner_matches_original_code_freeze": R0.sha256_file(ORIG / "run_abl_myr.py") == json.loads(
               (ORIG / "results/ABL-MYR-1_code_freeze.json").read_text())["runner"]["sha256"],
           "extracted_before_freeze": [], "gromacs_available": False}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps({k: rec[k] for k in ("frozen_utc", "spec_sha256", "code")}, indent=1))
    return 0


def verify_freeze() -> dict:
    if not FREEZE.exists():
        raise SystemExit("STOP: no check freeze record")
    f = json.loads(FREEZE.read_text())
    bad = [n for n, ok in (("check_spec.json", R0.sha256_file(SPEC) == f["spec_sha256"]),
                           ("run_comp_check.py", R0.sha256_file(SELF) == f["code"]["sha256"])) if not ok]
    bad += [p for p, h in f["preserved_sha256"].items() if _h(p) != h]
    if bad:
        raise SystemExit(f"STOP: frozen or preserved files differ: {bad}")
    return f


def cmd_run(args) -> int:
    f = verify_freeze()
    spec = json.loads(SPEC.read_text())
    res, ca = run_check(json.loads((ORIG / "run_spec.json").read_text()), ORIG / "data_dryad", ORIG / "data_dryad/extracted_1b", spec)
    res["check_spec_sha256"], res["check_code_sha256"] = f["spec_sha256"], f["code"]["sha256"]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / CA_NPZ, **{k.replace(":", "_"): v for k, v in ca.items() if v is not None})
    res["production_start_ca_npz"] = {"file": CA_NPZ, "sha256": R0.sha256_file(out / CA_NPZ), "units": "Å", "atoms": "protein Cα, whole"}
    (out / RESULTS).write_text(json.dumps(R0.scrub(res), indent=1, default=R0._json_default) + "\n")
    print(res["outcome"], json.dumps(res["resources"]))
    return 0


def cmd_report(args) -> int:
    res = json.loads((Path(args.out) / RESULTS).read_text())
    (Path(args.out) / "ABL-MYR-1b_check_report.md").write_text(render_report(res))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("freeze")
    for name in ("run", "report"):
        p = sub.add_parser(name)
        p.add_argument("--out", default=str(HERE / "results"))
    args = ap.parse_args(argv)
    return {"freeze": cmd_freeze, "run": cmd_run, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
