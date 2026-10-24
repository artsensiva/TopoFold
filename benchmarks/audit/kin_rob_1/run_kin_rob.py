"""KIN-ROB-1: robustness of named-atom extraction and kinase descriptors on a frozen panel of public kinase structures.

    <python with MDAnalysis> benchmarks/audit/kin_rob_1/run_kin_rob.py freeze
    <python with MDAnalysis> benchmarks/audit/kin_rob_1/run_kin_rob.py run

Engineering check only (see run_spec.json): no state labels, thresholds or classification. The run downloads the panel
(PDB, CC0) and UniProt sequences into git-ignored folders and records their SHA-256.
"""

import csv
import datetime as dt
import hashlib
import json
import math
import re
import resource
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs.atom_adapter import BACKBONE, extract_atoms, role_request  # noqa: E402
from topofold_inputs.ca_adapter import ResidueSpec, extract_ca, load_universe  # noqa: E402
from topofold_inputs.kinase_geometry import DESCRIPTOR_ATOMS, kinase_descriptors  # noqa: E402

SPEC = HERE / "run_spec.json"
FREEZE = HERE / "results/KIN-ROB-1_freeze.json"
CODE = ["benchmarks/audit/kin_rob_1/run_kin_rob.py", "python/topofold_inputs/atom_adapter.py",
        "python/topofold_inputs/kinase_geometry.py", "python/topofold_inputs/ca_adapter.py"]
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
         "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V", "MSE": "M"}
IDEAL = {("N", "CA"): 1.458, ("CA", "C"): 1.525, ("C", "O"): 1.231, ("CA", "CB"): 1.530}


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def freeze():
    if FREEZE.exists():
        print("freeze exists")
        return 1
    FREEZE.parent.mkdir(exist_ok=True)
    rec = {"id": "KIN-ROB-1", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "spec_sha256": sha256(SPEC),
           "code_sha256": {c: sha256(ROOT / c) for c in CODE},
           "panel_files_present_at_freeze": sorted(p.name for p in (HERE / "data_pdb").glob("*")) if (HERE / "data_pdb").exists() else []}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=1))
    return 0


def fetch(url, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        dest.write_bytes(urllib.request.urlopen(url, timeout=60).read())
    return {"path": str(dest.relative_to(ROOT)), "url": url, "bytes": dest.stat().st_size, "sha256": sha256(dest)}


# --- independent fixed-column reading (not MDAnalysis, not the adapters) ------------------------------------------------

def read_fixed(path, chain):
    """{(resid, icode, atom): xyz} first model, altloc blank-or-A; residue names {(resid, icode): name}; insertion codes."""
    atoms, names, icodes = {}, {}, set()
    for line in Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if line.startswith(("ATOM", "HETATM")) and line[21] == chain and line[16] in " A":
            key = (int(line[22:26]), line[26].strip())
            if line.startswith("HETATM") and line[17:20] not in THREE:
                continue
            atoms[key + (line[12:16].strip(),)] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            names[key] = line[17:20]
            if key[1]:
                icodes.add(key)
    return atoms, names, icodes


def remark_missing(path, chain):
    """REMARK 465 residues and REMARK 470 atoms for one chain."""
    miss_res, miss_atoms = set(), {}
    rx465 = re.compile(r"^REMARK 465\s+(?:\d+\s+)?([A-Z0-9]{1,3})\s+(\S)\s*(-?\d+)([A-Z]?)\s*$")
    rx470 = re.compile(r"^REMARK 470\s+(?:\d+\s+)?([A-Z0-9]{1,3})\s+(\S)\s*(-?\d+)([A-Z]?)\s+(.+)$")
    for line in Path(path).read_text().splitlines():
        m = rx465.match(line)
        if m and m.group(2) == chain:
            miss_res.add((int(m.group(3)), m.group(4)))
        m = rx470.match(line)
        if m and m.group(2) == chain and m.group(1) != "RES":
            miss_atoms[(int(m.group(3)), m.group(4))] = set(m.group(5).split())
    return miss_res, miss_atoms


def dist(a, b):
    return float(np.linalg.norm(b - a))


def dihed(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def independent_descriptors(atoms, roles):
    out = {}
    for name, needs in DESCRIPTOR_ATOMS.items():
        pts = [atoms.get((roles[r][0], roles[r][1], a)) for r, a in needs]
        if any(p is None for p in pts):
            out[name] = None
        else:
            out[name] = dist(*pts) if len(pts) == 2 else dihed(*pts)
    return out


def run():
    spec = json.loads(SPEC.read_text())
    frz = json.loads(FREEZE.read_text())
    bad = [c for c, v in frz["code_sha256"].items() if sha256(ROOT / c) != v]
    if frz["spec_sha256"] != sha256(SPEC) or bad:
        raise SystemExit(f"STOP: spec or code differs from the freeze: {bad}")
    t0 = time.time()
    inputs, rows, chains_out = {}, [], []
    for entry in spec["panel"]:
        pdb, acc = entry["pdb"], entry["uniprot"]
        inputs[pdb] = fetch(f"https://files.rcsb.org/download/{pdb}.pdb", HERE / f"data_pdb/{pdb}.pdb")
        inputs[acc] = fetch(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", HERE / f"data_uniprot/{acc}.fasta")
        seq = "".join(l.strip() for l in (HERE / f"data_uniprot/{acc}.fasta").read_text().splitlines() if not l.startswith(">"))
        pos = spec["role_positions_uniprot"][acc]
        # C1: role identity in UniProt
        aa = {r: (seq[p - 1] if 0 < p <= len(seq) else "?") for r, p in pos.items()}
        hrd = seq[pos["XHRD"]: pos["XHRD"] + 3]
        c1 = (aa["beta3_Lys"] == "K" and aa["alphaC_Glu"] == "E" and aa["DFG_Asp"] + aa["DFG_Phe"] + aa["DFG_Gly"] == "DFG"
              and hrd[0] in "HY" and hrd[1:] == "RD")
        pdb_path = ROOT / inputs[pdb]["path"]
        u = load_universe(pdb_path)
        for chain in entry["chains"]:
            rec = {"pdb": pdb, "chain": chain, "kinase": entry["kinase"], "uniprot": acc, "primary": chain == entry["chains"][0],
                   "C1_uniprot_role_identity": {"pass": c1, "residues": aa, "hrd_motif": hrd}}
            try:
                fixed, names, icodes = read_fixed(pdb_path, chain)
                observed = sorted({(k[0], k[1]) for k in fixed if k[2] == "CA"})
                # C2: offset mapping
                best = None
                for off in range(-700, 701):
                    plain = [(r, names[(r, ic)]) for r, ic in observed if not ic]
                    n = sum(1 for r, nm in plain if 0 < r + off <= len(seq) and THREE.get(nm, "X") == seq[r + off - 1])
                    if best is None or n > best[1]:
                        best = (off, n, len(plain))
                off, nmatch, nobs = best
                identity = nmatch / nobs if nobs else 0.0
                mism = [(r, names[(r, "")], seq[r + off - 1] if 0 < r + off <= len(seq) else "-") for r, ic in observed
                        if not ic and not (0 < r + off <= len(seq) and THREE.get(names[(r, "")], "X") == seq[r + off - 1])]
                role_pdb = {r: (p - off, "") for r, p in pos.items()}
                icode_in_roles = sorted(k for k in icodes if any(k[0] == v[0] for v in role_pdb.values()))
                rec["C2_numbering_mapping"] = {"pass": identity >= 0.95 and not icode_in_roles, "offset_uniprot_minus_pdb": off,
                                               "identity": identity, "observed_residues": nobs, "insertion_codes_in_chain": len(icodes),
                                               "insertion_codes_in_role_residues": icode_in_roles, "mismatches": mism[:30],
                                               "n_mismatches": len(mism)}
                # C3: residue names at role positions vs UniProt
                c3 = {}
                for r, (rid, ic) in role_pdb.items():
                    nm = names.get((rid, ic))
                    c3[r] = {"pdb_resid": rid, "pdb_resname": nm, "uniprot_aa": aa[r],
                             "agrees": nm is not None and THREE.get(nm, "X") == aa[r]}
                rec["C3_role_residue_identity_in_structure"] = {"all_agree": all(v["agrees"] for v in c3.values()), "roles": c3}
                roles = {r: ResidueSpec(chain, rid, ic) for r, (rid, ic) in role_pdb.items()}
                req = role_request(roles)
                # C4 + C6: extraction under both altloc policies
                runs = {pol: extract_atoms(u, f"{pdb}{chain}_{pol}", req, on_missing="mark", altloc=pol)
                        for pol in ("blank_then_A", "highest_occupancy")}
                ta = runs["blank_then_A"]
                st = {(a["resid"], a["icode"], a["name"]): a["status"] for a in ta.atoms}
                cov = sum(1 for v in st.values() if v == "present") / len(st)
                desc = {pol: kinase_descriptors(t, roles) for pol, t in runs.items()}
                rec["C4_required_atom_coverage"] = {"atoms_requested": len(st), "present_fraction": cov,
                                                    "missing": sorted(f"{k[0]}{k[1]}:{k[2]}={v}" for k, v in st.items() if v != "present"),
                                                    "descriptors_nan": {k: v for k, v in desc["blank_then_A"]["missing"].items()}}
                # C5: compare missing atoms/residues with REMARK 465/470
                m465, m470 = remark_missing(pdb_path, chain)
                role_res = {(rid, ic) for rid, ic in role_pdb.values()}
                ours_res = {(a["resid"], a["icode"]) for a in ta.atoms if a["status"] == "residue_not_found"}
                ours_atoms = {}
                for a in ta.atoms:
                    if a["status"] == "atom_not_found":
                        ours_atoms.setdefault((a["resid"], a["icode"]), set()).add(a["name"])
                want_res = role_res & m465
                want_atoms = {k: {n for n in v if any(n in names_ for (s, names_) in req if s.key == (chain, k[0], k[1]))}
                              for k, v in m470.items() if k in role_res}
                want_atoms = {k: v for k, v in want_atoms.items() if v}
                c5 = ours_res == want_res and ours_atoms == want_atoms
                rec["C5_missing_atom_consistency"] = {"pass": c5, "adapter_missing_residues": sorted(ours_res),
                                                      "remark465_role_residues": sorted(want_res),
                                                      "adapter_missing_atoms": {f"{k[0]}{k[1]}": sorted(v) for k, v in ours_atoms.items()},
                                                      "remark470_role_atoms": {f"{k[0]}{k[1]}": sorted(v) for k, v in want_atoms.items()}}
                # C6: altlocs, residue consistency, descriptor change
                mixed = []
                alt_res = set()
                for pol, t in runs.items():
                    per = {}
                    for a in t.atoms:
                        if a["altloc"]:
                            per.setdefault((a["resid"], a["icode"]), set()).add(a["altloc"])
                            alt_res.add((a["resid"], a["icode"]))
                    mixed += [(pol, k, sorted(v)) for k, v in per.items() if len(v) > 1]
                dmax = {}
                for k in DESCRIPTOR_ATOMS:
                    a, b = desc["blank_then_A"]["values"][k][0], desc["highest_occupancy"]["values"][k][0]
                    if np.isfinite(a) and np.isfinite(b):
                        dmax[k] = abs(math.atan2(math.sin(a - b), math.cos(a - b))) if len(DESCRIPTOR_ATOMS[k]) == 4 else abs(a - b)
                rec["C6_alternate_locations"] = {"pass": not mixed, "role_residues_with_altlocs": sorted(f"{r}{i}" for r, i in alt_res),
                                                 "mixed_conformers": mixed, "max_descriptor_change_between_policies": dmax}
                # C7: independent parser
                ind = independent_descriptors(fixed, role_pdb)
                diffs = {}
                for k, v in ind.items():
                    ours = float(desc["blank_then_A"]["values"][k][0])
                    if v is None or not np.isfinite(ours):
                        diffs[k] = None if (v is None) == (not np.isfinite(ours)) else "availability differs"
                    else:
                        diffs[k] = abs(math.atan2(math.sin(ours - v), math.cos(ours - v))) if len(DESCRIPTOR_ATOMS[k]) == 4 else abs(ours - v)
                num = [d for d in diffs.values() if isinstance(d, float)]
                rec["C7_independent_parser"] = {"pass": all(not isinstance(d, str) for d in diffs.values()) and (max(num) if num else 0.0) <= 1e-3,
                                                "max_abs_difference": max(num) if num else None, "n_compared": len(num),
                                                "availability_mismatches": [k for k, d in diffs.items() if isinstance(d, str)]}
                # C8: Calpha identity with extract_ca
                uniq = sorted({(rid, ic) for rid, ic in role_pdb.values()})
                c8 = True
                for pol in ("blank_then_A", "highest_occupancy"):
                    specs = [ResidueSpec(chain, rid, ic) for rid, ic in uniq]
                    ca = extract_ca(u, "ca", correspondence=specs, on_missing="mark", altloc=pol)
                    at = extract_atoms(u, "at", [(s, ["CA"]) for s in specs], on_missing="mark", altloc=pol)
                    c8 &= bool(np.array_equal(ca.coordinates, at.coordinates, equal_nan=True))
                rec["C8_calpha_identity"] = {"pass": c8, "residues": len(uniq)}
                # C9: determinism
                again = extract_atoms(u, "again", req, on_missing="mark", altloc="blank_then_A")
                rec["C9_determinism"] = {"pass": bool(np.array_equal(again.coordinates, ta.coordinates, equal_nan=True))}
                # C10: geometry plausibility of role residues
                bb = extract_atoms(u, "bb", [(ResidueSpec(chain, rid, ic), list(BACKBONE) + ["CB"]) for rid, ic in uniq],
                                   on_missing="mark", altloc="blank_then_A")
                warn = []
                for rid, ic in uniq:
                    for (x, y), ideal in IDEAL.items():
                        p, q = bb.xyz(chain, rid, x, ic)[0], bb.xyz(chain, rid, y, ic)[0]
                        if np.isfinite(p).all() and np.isfinite(q).all() and abs(dist(p, q) - ideal) > 0.10:
                            warn.append(f"{rid}{ic} {x}-{y} {dist(p, q):.3f}")
                phe = role_pdb["DFG_Phe"]
                cg, cz = fixed.get((phe[0], phe[1], "CG")), fixed.get((phe[0], phe[1], "CZ"))
                if cg is not None and cz is not None and abs(dist(cg, cz) - 2.78) > 0.10:
                    warn.append(f"{phe[0]} CG-CZ {dist(cg, cz):.3f}")
                rec["C10_geometry_plausibility"] = {"warnings": warn}
                # outcome
                fail = (not c1 or not rec["C2_numbering_mapping"]["pass"] or not c5 or not rec["C7_independent_parser"]["pass"]
                        or not c8 or not rec["C9_determinism"]["pass"] or mixed)
                excl = bool(desc["blank_then_A"]["missing"]) or not rec["C3_role_residue_identity_in_structure"]["all_agree"]
                rec["outcome"] = "fail" if fail else ("pass_with_exclusions" if excl else "pass")
                rec["descriptors_blank_then_A"] = {k: (None if not np.isfinite(v[0]) else float(v[0]))
                                                   for k, v in desc["blank_then_A"]["values"].items()}
                chains_out.append(rec)
                rows.append({"pdb": pdb, "chain": chain, "kinase": entry["kinase"], "outcome": rec["outcome"],
                             "C1": c1, "C2_identity": round(identity, 4), "C2_offset": off, "C3_roles_agree": rec["C3_role_residue_identity_in_structure"]["all_agree"],
                             "C4_coverage": round(cov, 4), "C5": c5, "C6_altloc_roles": len(alt_res), "C6_mixed": len(mixed),
                             "C7_max_diff": rec["C7_independent_parser"]["max_abs_difference"], "C8": c8, "C9": rec["C9_determinism"]["pass"],
                             "C10_warnings": len(warn), "descriptors_nan": ";".join(sorted(desc["blank_then_A"]["missing"]))})
            except Exception as e:   # an adapter refusal (e.g. altloc tie) is a recorded failure, not a crash
                rec["outcome"] = "fail"
                rec["error"] = f"{type(e).__name__}: {e}"
                chains_out.append(rec)
                rows.append({"pdb": pdb, "chain": chain, "kinase": entry["kinase"], "outcome": "fail", "C1": c1,
                             "error": rec["error"][:200]})
    ru, rc = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    out = {"id": "KIN-ROB-1", "freeze_utc": frz["frozen_utc"], "spec_sha256": frz["spec_sha256"], "inputs": inputs,
           "chains": chains_out, "summary": {o: sum(1 for c in chains_out if c["outcome"] == o) for o in ("pass", "pass_with_exclusions", "fail")},
           "resources": {"cpu_s": ru.ru_utime + ru.ru_stime + rc.ru_utime + rc.ru_stime, "wall_s": time.time() - t0,
                         "download_bytes": sum(v["bytes"] for v in inputs.values())}}
    (HERE / "results/KIN-ROB-1_results.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
    with open(HERE / "results/KIN-ROB-1_summary.csv", "w", newline="") as fh:
        fields = list(dict.fromkeys(k for r in rows for k in r))
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(json.dumps(out["summary"]), json.dumps(out["resources"]))
    for r in rows:
        print(r)
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("freeze", "run"):
        print(__doc__)
        sys.exit(2)
    sys.exit(freeze() if sys.argv[1] == "freeze" else run())
