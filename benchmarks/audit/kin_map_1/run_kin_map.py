"""KIN-MAP-1: acceptance check of the kinase residue-role mapping helper on the frozen KIN-ROB-1 panel.

    <python with MDAnalysis> benchmarks/audit/kin_map_1/run_kin_map.py freeze
    <python with MDAnalysis> benchmarks/audit/kin_map_1/run_kin_map.py run

Engineering check only (see run_spec.json). Reads the KIN-ROB-1 downloads (git-ignored; nothing is downloaded) and the
KIN-ROB-1 / KIN-ROB-1a records without modifying them.
"""

import csv
import datetime as dt
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs.atom_adapter import extract_atoms  # noqa: E402
from topofold_inputs.ca_adapter import load_universe  # noqa: E402
from topofold_inputs.kinase_geometry import kinase_descriptors  # noqa: E402
from topofold_inputs.kinase_roles import (Numbering, Reference, RoleMappingError, map_kinase_roles,  # noqa: E402
                                          observed_chain)

SPEC = HERE / "run_spec.json"
FREEZE = HERE / "results/KIN-MAP-1_freeze.json"
ROB, ROB_A = ROOT / "benchmarks/audit/kin_rob_1", ROOT / "benchmarks/audit/kin_rob_1a"
CODE = ["benchmarks/audit/kin_map_1/run_kin_map.py", "python/topofold_inputs/kinase_roles.py",
        "tests/fixtures/kin_rob_1_panel_observed_chains.json"]
UNCHANGED_CODE = ["python/topofold_inputs/ca_adapter.py", "python/topofold_inputs/atom_adapter.py",
                  "python/topofold_inputs/kinase_geometry.py"]
FROZEN = [f"benchmarks/audit/{d}/{f}" for d in ("kin_rob_1", "kin_rob_1a") for f in
          ("run_spec.json", "run_kin_rob.py", "README.md", "results/KIN-ROB-1_freeze.json",
           "results/KIN-ROB-1_results.json", "results/KIN-ROB-1_summary.csv")]


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def freeze():
    if FREEZE.exists():
        print("freeze exists")
        return 1
    FREEZE.parent.mkdir(exist_ok=True)
    rec = {"id": "KIN-MAP-1", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "spec_sha256": sha256(SPEC),
           "code_sha256": {c: sha256(ROOT / c) for c in CODE},
           "frozen_artifacts_sha256": {c: sha256(ROOT / c) for c in FROZEN}}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=1))
    return 0


def anchors_of(pos):
    return {"beta3_Lys": pos["beta3_Lys"], "alphaC_Glu": pos["alphaC_Glu"], "HRD_H": pos["XHRD"] + 1,
            "DFG_Asp": pos["DFG_Asp"]}


def recorded_roles(rec):
    c3 = rec["C3_role_residue_identity_in_structure"]["roles"]
    gone = {int(s.split(":")[0]) for s in rec["C4_required_atom_coverage"]["missing"] if "residue_not_found" in s}
    return {r: (None if v["pdb_resid"] in gone else [v["pdb_resid"], ""]) for r, v in c3.items()}


def run():
    spec, frz = json.loads(SPEC.read_text()), json.loads(FREEZE.read_text())
    bad = [c for c, v in frz["code_sha256"].items() if sha256(ROOT / c) != v]
    if frz["spec_sha256"] != sha256(SPEC) or bad:
        raise SystemExit(f"STOP: spec or code differs from the freeze: {bad}")
    t0 = time.time()
    rob_spec = json.loads((ROB / "run_spec.json").read_text())
    rob_res = json.loads((ROB / "results/KIN-ROB-1_results.json").read_text())
    rob_a_pos = json.loads((ROB_A / "run_spec.json").read_text())["role_positions_uniprot"]["P05132"]
    rob_a_rec = json.loads((ROB_A / "results/KIN-ROB-1_results.json").read_text())["chains"][0]
    rob_frz = json.loads((ROB / "results/KIN-ROB-1_freeze.json").read_text())
    bad_inputs = [k for k, v in rob_res["inputs"].items() if sha256(ROOT / v["path"]) != v["sha256"]]
    if bad_inputs:
        raise SystemExit(f"STOP: KIN-ROB-1 inputs differ from their records: {bad_inputs}")
    fixture = {(c["pdb"], c["chain"]): c["residues"] for c in
               json.loads((ROOT / "tests/fixtures/kin_rob_1_panel_observed_chains.json").read_text())["chains"]}
    records = {(c["pdb"], c["chain"]): c for c in rob_res["chains"]}
    pos_all = rob_spec["role_positions_uniprot"]

    m7 = {"unchanged_code": {c: sha256(ROOT / c) == rob_frz["code_sha256"][c] for c in UNCHANGED_CODE},
          "frozen_artifacts": {c: sha256(ROOT / c) == v for c, v in frz["frozen_artifacts_sha256"].items()}}
    m7["pass"] = all(m7["unchanged_code"].values()) and all(m7["frozen_artifacts"].values())
    cases, rows, m5 = [], [], {}
    for entry in rob_spec["panel"]:
        pdb, acc = entry["pdb"], entry["uniprot"]
        u = load_universe(ROB / f"data_pdb/{pdb}.pdb")
        fasta = (ROB / f"data_uniprot/{acc}.fasta").read_text()
        ref = Reference.from_fasta(fasta, acc, f"{acc}-1", f"{rob_res['inputs'][acc]['url']} sha256 "
                                                          f"{rob_res['inputs'][acc]['sha256']}")
        for chain in entry["chains"]:
            obs = observed_chain(u, chain)
            m5[f"{pdb} {chain}"] = [[o.resid, o.icode, o.resname, o.bonded_to_previous] for o in obs] == fixture[(pdb, chain)]
            pka = pdb == "1ATP"
            want_pos = rob_a_pos if pka else pos_all[acc]
            want_res = recorded_roles(rob_a_rec if pka else records[(pdb, chain)])
            want_desc = (rob_a_rec if pka else records[(pdb, chain)])["descriptors_blank_then_A"]
            author = {a: want_res[r][0] for a, r in (("beta3_Lys", "beta3_Lys"), ("alphaC_Glu", "alphaC_Glu"),
                                                      ("DFG_Asp", "DFG_Asp"))}
            author["HRD_H"] = want_res["XHRD"][0] + 1
            todo = [("frozen_spec", Numbering.reference_index("as labelled in the KIN-ROB-1 spec (UniProt)"),
                     anchors_of(pos_all[acc]), not pka),
                    ("structure_author", Numbering("structure_author"), author, True)]
            if pka:
                todo += [("pka_kin_rob_1a", Numbering.reference_index("UniProt P05132 (KIN-ROB-1a)"),
                          anchors_of(rob_a_pos), True),
                         ("pka_literature_declared", Numbering("reference", 1, "PKA literature numbering (no "
                                                               "initiator Met)"), anchors_of(pos_all[acc]), True),
                         ("pka_structure_author_original", Numbering("structure_author"),
                          {a: v - 1 for a, v in author.items()}, False)]
            for case, numbering, anchors, expect_accept in todo:
                out = {"case": case, "pdb": pdb, "chain": chain, "uniprot": acc, "anchors": anchors,
                       "numbering": {"scheme": numbering.scheme, "reference_minus_label": numbering.reference_minus_label,
                                     "description": numbering.description}, "expected_accepted": expect_accept}
                try:
                    m = map_kinase_roles(obs, chain, ref, numbering, anchors)
                    rep, accepted = m.report, True
                except RoleMappingError as e:
                    m, rep, accepted = None, e.report, False
                out["accepted"] = accepted
                out["M1_decision"] = accepted == expect_accept
                if accepted:
                    got_pos = {r: v["reference_index"] for r, v in rep["roles"].items()}
                    got_res = {r: (None if s is None else [s.resid, s.icode]) for r, s in m.roles.items()}
                    out["M2_reference_positions"] = got_pos == want_pos
                    out["M3_structure_residues"] = got_res == want_res
                    traj = extract_atoms(u, f"{pdb}{chain}", m.atom_request(), on_missing="mark", altloc="blank_then_A")
                    d = kinase_descriptors(traj, m.residue_specs())["values"]
                    mism = [k for k, v in want_desc.items()
                            if not (np.isnan(d[k][0]) if v is None else float(d[k][0]) == v)]
                    out["M4_descriptors"] = {"pass": not mism, "mismatches": mism,
                                             "nan": sorted(k for k, v in d.items() if np.isnan(v[0]))}
                else:
                    motif = [p for p in rep["problems"] if "expected" in p]
                    out["M6_refusal_content"] = {
                        "motif_problems": len(motif),
                        "diagnostic_offsets": rep.get("diagnostic_offsets_not_applied", {}).get("d"),
                        "pass": len(motif) == 4 and (numbering.scheme != "reference"
                                                     or rep["diagnostic_offsets_not_applied"]["d"] == [1])}
                out["problems"] = rep["problems"]
                out["alignment"] = {k: rep["alignment"][k] for k in
                                    ("aligned_pairs", "identity", "reference_span", "mutations", "insertions",
                                     "unaligned_terminal", "unobserved_reference_ranges",
                                     "absent_between_bonded_residues", "chain_breaks_before", "bonding_unknown",
                                     "insertion_code_residues", "variant_residue_names",
                                     "author_numbering_diagnostic")}
                out["roles"] = rep.get("roles")
                checks = [out["M1_decision"], out.get("M2_reference_positions", True),
                          out.get("M3_structure_residues", True), out.get("M4_descriptors", {}).get("pass", True),
                          out.get("M6_refusal_content", {}).get("pass", True)]
                out["pass"] = all(checks)
                cases.append(out)
                aln = rep["alignment"]
                rows.append({"case": case, "pdb": pdb, "chain": chain, "expected": "accept" if expect_accept else
                             "refuse", "decision": "accepted" if accepted else "refused", "pass": out["pass"],
                             "M2": out.get("M2_reference_positions", ""), "M3": out.get("M3_structure_residues", ""),
                             "M4": out.get("M4_descriptors", {}).get("pass", ""),
                             "M6": out.get("M6_refusal_content", {}).get("pass", ""),
                             "identity": round(aln["identity"], 4), "mutations": len(aln["mutations"]),
                             "unobserved_ranges": len(aln["unobserved_reference_ranges"]),
                             "unobserved_roles": ";".join(r for r, v in (rep.get("roles") or {}).items()
                                                          if v["status"] == "unobserved"),
                             "author_offsets": json.dumps(aln["author_numbering_diagnostic"]
                                                          ["reference_index_minus_author_resid"])})
    ru = resource.getrusage(resource.RUSAGE_SELF)
    out = {"id": "KIN-MAP-1", "freeze_utc": frz["frozen_utc"], "spec_sha256": frz["spec_sha256"],
           "M5_observed_chain": {"pass": all(m5.values()), "per_chain": m5}, "M7_frozen_artifacts_and_code": m7,
           "summary": {"cases": len(cases), "cases_pass": sum(c["pass"] for c in cases),
                       "accepted": sum(c["accepted"] for c in cases),
                       "refused": sum(not c["accepted"] for c in cases),
                       "M5_chains_pass": sum(m5.values()), "M5_chains": len(m5), "M7_pass": m7["pass"]},
           "cases": cases, "resources": {"cpu_s": ru.ru_utime + ru.ru_stime, "wall_s": time.time() - t0,
                                         "download_bytes": 0}}
    (HERE / "results/KIN-MAP-1_results.json").write_text(json.dumps(out, indent=1) + "\n")
    with open(HERE / "results/KIN-MAP-1_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
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
