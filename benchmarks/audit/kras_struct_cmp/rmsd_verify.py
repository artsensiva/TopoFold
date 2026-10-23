"""KRAS-STRUCT-1 closeout: independent check of the mapped global Calpha RMSD (does not modify frozen artifacts).

    .venv/bin/python benchmarks/audit/kras_struct_cmp/rmsd_verify.py

Independent of the run's pipeline: a minimal fixed-column PDB reader (not the adapter, not MDAnalysis) and a numpy
SVD/Kabsch superposition with the determinant correction (proper rotation only). Declared before running:
  1. primary: chain A residues 1-168, Calpha, altloc blank-then-A (the frozen KRAS-STRUCT-1 mapping);
  2. altloc sensitivity: highest occupancy;
  3. per-residue deviations after the primary superposition (distribution and largest values);
  4. context only, not a reproduction: an emulation of iterative outlier rejection in the style of PyMOL `align`
     (documented defaults cutoff 2.0 "in RMS", 5 cycles; RMSD reported after rejection). Interpreted here as rejecting
     pairs whose deviation exceeds 2.0 x the current RMSD. PyMOL is not installed and this emulation is not verified
     against it; PyMOL also performs its own sequence alignment and atom selection.
No other subsets are evaluated.
"""

import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / "results/KRAS-STRUCT-1_rmsd_verification.json"


def read_ca(path, chain, resids, altloc_policy):
    """{resid: (resname, xyz)} for Calpha atoms of the first model, fixed PDB columns."""
    cand = {}
    for line in Path(path).read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")) or line[12:16].strip() != "CA" or line[21] != chain:
            continue
        if line[26] != " ":
            raise ValueError(f"insertion code in {path}: {line}")
        r = int(line[22:26])
        if r in resids:
            cand.setdefault(r, []).append((line[16], float(line[54:60]), line[17:20],
                                           np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    out = {}
    for r in resids:
        alts = cand[r]
        if altloc_policy == "blank_then_A":
            pick = next((a for a in alts if a[0] == " "), None) or next(a for a in alts if a[0] == "A")
        else:
            occ = sorted(alts, key=lambda a: -a[1])
            if len(occ) > 1 and occ[0][1] == occ[1][1]:
                raise ValueError(f"occupancy tie at {r}")
            pick = occ[0]
        out[r] = (pick[2], pick[3])
    return out


def kabsch(P, Q):
    """Proper rotation R and translation minimizing |(P - cP) R - (Q - cQ)|; returns rmsd, per-atom deviations."""
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    fit = (P - cP) @ R.T + cQ
    dev = np.linalg.norm(fit - Q, axis=1)
    return float(np.sqrt(np.mean(dev ** 2))), dev, R


def self_test():
    rng = np.random.default_rng(0)
    P = rng.normal(size=(30, 3))
    a = 0.7
    Rz = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
    rot_ok = kabsch(P, P @ Rz.T + [1, 2, 3])[0] < 1e-12
    mirror = P * [1, 1, -1]
    r_m, _, R = kabsch(P, mirror)
    return {"recovers_rotation": bool(rot_ok), "mirror_rmsd_positive": bool(r_m > 0.1),
            "det_R_on_mirror_case": float(np.linalg.det(R))}


def iterative(P, Q, cutoff=2.0, cycles=5):
    keep = np.ones(len(P), bool)
    for _ in range(cycles):
        r, dev, _ = kabsch(P[keep], Q[keep])
        new = keep.copy()
        new[keep] = dev <= cutoff * r
        if new.sum() == keep.sum() or new.sum() < 3:
            break
        keep = new
    r, _, _ = kabsch(P[keep], Q[keep])
    return r, keep


def main():
    spec = json.loads((HERE / "run_spec.json").read_text())
    resids = list(range(1, 169))
    files = {s: ROOT / spec["structures"][s]["file"] for s in ("GDP", "GNP")}
    for s, f in files.items():
        assert hashlib.sha256(f.read_bytes()).hexdigest() == spec["structures"][s]["sha256"], s
    res = {"id": "KRAS-STRUCT-1 RMSD verification", "self_test": self_test(), "policies": {}}
    for pol in ("blank_then_A", "highest_occupancy"):
        a = read_ca(files["GDP"], "A", resids, pol)
        b = read_ca(files["GNP"], "A", resids, pol)
        assert all(a[r][0] == b[r][0] for r in resids), "residue names differ"
        P = np.array([a[r][1] for r in resids])
        Q = np.array([b[r][1] for r in resids])
        r, dev, _ = kabsch(P, Q)
        ent = {"n_ca": len(resids), "rmsd_A": r,
               "deviation_quantiles_A": {q: float(np.quantile(dev, float(q))) for q in ("0.5", "0.75", "0.9", "0.95", "1.0")},
               "n_deviation_gt_1A": int((dev > 1.0).sum()),
               "largest_deviations": [{"resid": resids[i], "resname": a[resids[i]][0], "dev_A": float(dev[i])}
                                      for i in np.argsort(-dev)[:12]]}
        if pol == "blank_then_A":
            ri, keep = iterative(P, Q)
            ent["context_iterative_rejection_emulation"] = {
                "rmsd_A": ri, "n_kept": int(keep.sum()), "rejected_resids": [resids[i] for i in np.where(~keep)[0]],
                "status": "illustrative emulation of a documented PyMOL default; not verified against PyMOL"}
            ent["pipeline_value_mdanalysis_qcp_A"] = json.loads(
                (HERE / "results/KRAS-STRUCT-1_results.json").read_text())["global_rmsd_A"]
            ent["abs_difference_vs_pipeline_A"] = abs(r - ent["pipeline_value_mdanalysis_qcp_A"])
        res["policies"][pol] = ent
    res["paper"] = {
        "statement": "They superimpose with a low RMSD (0.3 Å) and show minor differences in Switch-I and -II "
                     "conformations (Supplementary Fig. 5e).",
        "source": "Bonsor et al., Nat Commun 17:1614 (2026), Results; PMC12905138",
        "definition_given": "none: no atom set, residue range, software or outlier rejection is stated for this value "
                            "(the article's only atom specification, 'over Calpha atoms', belongs to a different "
                            "comparison of cryo-EM complexes); Methods name PyMOL for figures only",
        "referent": "the preceding sentence mentions existing KRAS-G12D-MRTX1133 structures, so 'They' most plausibly "
                    "means the two new structures but is grammatically ambiguous",
        "not_inspected": "Supplementary Information PDF (37.2 MB, exceeds the 25 MB download budget); the peer-review "
                         "file (0.3 MB) was searched and contains no RMSD statement"}
    OUT.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
