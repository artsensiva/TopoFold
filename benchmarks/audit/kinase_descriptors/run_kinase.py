"""KIN-DESC-1: which representations separate stated DFG / A-loop states on four public Abl chains.

    .venv/bin/python benchmarks/audit/kinase_descriptors/run_kinase.py freeze
    .venv/bin/python benchmarks/audit/kinase_descriptors/run_kinase.py run

Representation capability only (see run_spec.json, permitted_interpretation). Atoms are read with a fixed-column PDB
reader (first model, blank or 'A' altloc); kappa/tau come from the TopoFold core via topofold_inputs.features.
"""

import csv
import datetime as dt
import hashlib
import itertools
import json
import math
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs.features import compute_ca_features  # noqa: E402

SPEC = HERE / "run_spec.json"
FREEZE = HERE / "results/KIN-DESC-1_freeze.json"
CODE = ["benchmarks/audit/kinase_descriptors/run_kinase.py", "python/topofold_inputs/features.py"]
RANGE = range(240, 501)


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def freeze():
    if FREEZE.exists():
        print("freeze exists")
        return 1
    FREEZE.parent.mkdir(exist_ok=True)
    rec = {"id": "KIN-DESC-1", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "spec_sha256": sha256(SPEC), "code_sha256": {c: sha256(ROOT / c) for c in CODE}}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=1))
    return 0


def read_atoms(path, chain, offset):
    """{(resid_1a, atom_name): xyz} for one chain, first model, altloc blank or A."""
    out = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM") or line[21] != chain or line[16] not in " A":
            continue
        r = int(line[22:26]) + offset
        out[(r, line[12:16].strip())] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return out


def dist(a, b):
    return float(np.linalg.norm(a - b))


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    return float(math.atan2(np.dot(np.cross(b1n, v), w), np.dot(v, w)))


def kabsch_rmsd(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - cP).T @ (Q - cQ))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return float(np.sqrt(np.mean(np.sum(((P - cP) @ R.T + cQ - Q) ** 2, axis=1))))


def established(at):
    g = lambda r, a: at[(r, a)]  # noqa: E731
    return {"D1": dist(g(290, "CA"), g(382, "CZ")), "D2": dist(g(271, "CA"), g(382, "CZ")),
            "alphaC_CbCb": dist(g(286, "CB"), g(271, "CB")), "Aloop_Hbond_NO": dist(g(386, "N"), g(360, "O")),
            "Phe382_chi1": dihedral(g(382, "N"), g(382, "CA"), g(382, "CB"), g(382, "CG")),
            **{f"phi{r}": dihedral(g(r - 1, "C"), g(r, "N"), g(r, "CA"), g(r, "C")) for r in (380, 381, 382)},
            **{f"psi{r}": dihedral(g(r, "N"), g(r, "CA"), g(r, "C"), g(r + 1, "N")) for r in (380, 381, 382)}}


ANGULAR = {"Phe382_chi1", "phi380", "phi381", "phi382", "psi380", "psi381", "psi382", "kappa", "tau"}


def diff(a, b, angular):
    d = abs(a - b)
    return abs(math.atan2(math.sin(a - b), math.cos(a - b))) if angular else d


def separation(values, groups, angular, floor):
    """values: {chain: x}; groups: [list, list]. Returns between, within, separates."""
    between = min(diff(values[i], values[j], angular) for i in groups[0] for j in groups[1])
    within = max([diff(values[i], values[k], angular) for g in groups for i, k in itertools.combinations(g, 2)] or [0.0])
    return between, within, bool(between > within and between >= floor)


def pair_separation(pairval, groups, floor):
    between = min(pairval[frozenset((i, j))] for i in groups[0] for j in groups[1])
    within = max([pairval[frozenset((i, k))] for g in groups for i, k in itertools.combinations(g, 2)] or [0.0])
    return between, within, bool(between > within and between >= floor)


def run():
    spec = json.loads(SPEC.read_text())
    frz = json.loads(FREEZE.read_text())
    bad = [c for c, v in frz["code_sha256"].items() if sha256(ROOT / c) != v]
    if frz["spec_sha256"] != sha256(SPEC) or bad:
        raise SystemExit(f"STOP: spec or code differs from the freeze: {bad}")
    t0 = time.time()
    chains, ca, est = {}, {}, {}
    for pid, s in spec["structures"].items():
        if sha256(ROOT / s["file"]) != s["sha256"]:
            raise SystemExit(f"STOP: {pid} hash differs")
        for c in s["chains"]:
            key = f"{pid}:{c}"
            at = read_atoms(ROOT / s["file"], c, s["numbering_offset_to_cAbl_1a"])
            missing = [r for r in RANGE if (r, "CA") not in at]
            if missing:
                raise SystemExit(f"STOP: {key} lacks CA at {missing[:5]}")
            chains[key] = at
            ca[key] = np.array([at[(r, "CA")] for r in RANGE])
            est[key] = established(at)
    keys = list(chains)
    resid = list(RANGE)
    n = len(resid)
    feats = {}
    for k in keys:
        f = compute_ca_features(ca[k][None], np.ones(n, bool), [[0, n]])
        feats[k] = {"kappa": f.kappa[0], "tau": f.tau[0], "kv": f.kappa_valid[0], "tv": f.tau_valid[0]}
    # per-window scalar metrics
    win = {}
    for i in range(1, n - 1):
        win[("kappa", resid[i])] = {k: float(feats[k]["kappa"][i]) if feats[k]["kv"][i] else None for k in keys}
        win[("d_i_i2", resid[i])] = {k: dist(ca[k][i - 1], ca[k][i + 1]) for k in keys}
    for i in range(1, n - 2):
        win[("tau", resid[i])] = {k: float(feats[k]["tau"][i]) if feats[k]["tv"][i] else None for k in keys}
        win[("d_i_i3", resid[i])] = {k: dist(ca[k][i - 1], ca[k][i + 2]) for k in keys}
    pair = {}
    for i in range(1, n - 1):
        pair[("rmsd3", resid[i])] = {frozenset((a, b)): kabsch_rmsd(ca[a][i - 1:i + 2], ca[b][i - 1:i + 2])
                                     for a, b in itertools.combinations(keys, 2)}
    for i in range(1, n - 2):
        pair[("rmsd4", resid[i])] = {frozenset((a, b)): kabsch_rmsd(ca[a][i - 1:i + 3], ca[b][i - 1:i + 3])
                                     for a, b in itertools.combinations(keys, 2)}
    span = {"kappa": (1, 1), "d_i_i2": (1, 1), "rmsd3": (1, 1), "tau": (1, 2), "d_i_i3": (1, 2), "rmsd4": (1, 2)}
    proxies = {k: {"proxy_D1": dist(chains[k][(290, "CA")], chains[k][(382, "CA")]),
                   "proxy_D2": dist(chains[k][(271, "CA")], chains[k][(382, "CA")]),
                   "proxy_alphaC": dist(chains[k][(286, "CA")], chains[k][(271, "CA")]),
                   "proxy_Aloop": dist(chains[k][(386, "CA")], chains[k][(360, "CA")])} for k in keys}
    regions = {"DFG": (378, 386), "A_loop": (381, 409), "alphaC": (279, 293)}
    floors = spec["decision_rule"]["floors"]
    contrasts = {"DFG": [spec["labels"]["DFG"]["in"], spec["labels"]["DFG"]["out"]],
                 "A_loop": [spec["labels"]["A_loop"]["closed"], spec["labels"]["A_loop"]["open"]]}
    acin = {k: est[k]["alphaC_CbCb"] <= 10.0 for k in keys}
    if len(set(acin.values())) == 2:
        contrasts["alphaC"] = [[k for k in keys if acin[k]], [k for k in keys if not acin[k]]]
    results = {}
    for cname, groups in contrasts.items():
        lo, hi = regions[cname]
        rows = {}
        for d in ("D1", "D2", "alphaC_CbCb", "Aloop_Hbond_NO", "Phe382_chi1", "phi380", "phi381", "phi382",
                  "psi380", "psi381", "psi382"):
            ang = d in ANGULAR
            b, w, s = separation({k: est[k][d] for k in keys}, groups, ang, floors["angle_rad" if ang else "distance_and_rmsd_A"])
            rows[f"established:{d}"] = {"between": b, "within": w, "separates": s}
        for d in ("proxy_D1", "proxy_D2", "proxy_alphaC", "proxy_Aloop"):
            b, w, s = separation({k: proxies[k][d] for k in keys}, groups, False, floors["distance_and_rmsd_A"])
            rows[f"calpha_proxy:{d}"] = {"between": b, "within": w, "separates": s}
        for metric in ("kappa", "tau", "d_i_i2", "d_i_i3", "rmsd3", "rmsd4"):
            a, z = span[metric]
            centres = [r for r in resid if r - a >= lo and r + z <= hi]
            tested, sep, best = 0, [], None
            for r in centres:
                if metric.startswith("rmsd"):
                    if (metric, r) not in pair:
                        continue
                    b, w, s = pair_separation(pair[(metric, r)], groups, floors["distance_and_rmsd_A"])
                else:
                    vals = win.get((metric, r))
                    if vals is None or any(v is None for v in vals.values()):
                        continue
                    ang = metric in ANGULAR
                    b, w, s = separation(vals, groups, ang, floors["angle_rad" if ang else "distance_and_rmsd_A"])
                tested += 1
                if s:
                    sep.append(r)
                if best is None or b - w > best["margin"]:
                    best = {"centre": r, "between": b, "within": w, "margin": b - w}
            rows[f"window:{metric}"] = {"windows_tested": tested, "windows_separating": len(sep),
                                        "separating_centres": sep, "largest_margin_window": best}
        results[cname] = {"groups": groups, "region": [lo, hi], "rows": rows}
    ru = resource.getrusage(resource.RUSAGE_SELF)
    out = {"id": "KIN-DESC-1", "freeze_utc": frz["frozen_utc"], "spec_sha256": frz["spec_sha256"], "chains": keys,
           "established_values": est, "calpha_proxy_values": proxies, "alphaC_in_by_CbCb_le_10A": acin,
           "alphaC_contrast": "evaluated" if "alphaC" in contrasts else "not assessable: all four chains fall in one alphaC class",
           "contrasts": results, "resources": {"cpu_s": ru.ru_utime + ru.ru_stime, "wall_s": time.time() - t0}}
    (HERE / "results/KIN-DESC-1_results.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    with open(HERE / "results/KIN-DESC-1_summary.csv", "w", newline="") as fh:
        wr = csv.writer(fh, lineterminator="\n")
        wr.writerow(["contrast", "representation", "outcome", "detail"])
        for cname, c in results.items():
            for rep, v in c["rows"].items():
                if rep.startswith("window:"):
                    outc = "separates in >= 1 tested element" if v["windows_separating"] else "separates in none"
                    det = f"{v['windows_separating']}/{v['windows_tested']} windows; centres {v['separating_centres']}"
                else:
                    outc = "separates" if v["separates"] else "does not separate"
                    det = f"between {v['between']:.3f}, within {v['within']:.3f}"
                wr.writerow([cname, rep, outc, det])
    print(json.dumps({"alphaC": out["alphaC_contrast"], "established": est, "resources": out["resources"]}, indent=1, default=float))
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("freeze", "run"):
        print(__doc__)
        sys.exit(2)
    sys.exit(freeze() if sys.argv[1] == "freeze" else run())
