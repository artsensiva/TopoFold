"""Verified structural-role annotation of the existing PIM1 comparisons (rules: DESIGN.md, fixed before this ran).

    .venv/bin/python benchmarks/audit/pim1_structural_annotation_2026-10-26/annotate_pim1.py

Reads local structures (hash-verified), the recorded reference sequence, the existing comparison outputs and one small UniProt annotation record;
runs the existing kinase_roles helper through the adapter interpreter. Recomputes no comparison score and changes no threshold.
Output: results/pim1_annotation.json.
"""

import csv
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
KP = ROOT / "benchmarks/audit/kin_pair_1"
PH = ROOT / "benchmarks/audit/structpair_phospho_support_2026-10-26/results"
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(KP))
sys.path.insert(0, str(ROOT / "benchmarks/audit/modified_residue_inventory_2026-10-26"))
import kinpair as kp  # noqa: E402
import inventory as inv  # noqa: E402  (scan, pick only)

ADAPTER = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
ENTRIES = {"8R0Q": "PIM1-1:A", "6NO9": "PIM1-1:B", "6VRV": "PIM1-2:A", "6L17": "PIM1-2:B"}
ACC, ISO = "P11309", "P11309-1"
WIN_LO, WIN_HI, CONTACT = 77, 97, 4.0           # DESIGN.md: alpha-C Glu rule
DFG_WINDOW = (168, 290)
HRD_H, BETA3 = 165, 67
M3, M4 = ("d_kappa", "d_i_i2", "rmsd3"), ("d_tau", "d_i_i3", "rmsd4")


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def helices_strands(text, chain, index):
    """Deposited HELIX and SHEET records of the chain, with the author range and its reference range where mapped."""
    out = {"helices": [], "strands": []}
    for ln in text.splitlines():
        if ln.startswith("HELIX ") and ln[19] == chain:
            a, b = int(ln[21:25]), int(ln[33:37])
            out["helices"].append({"id": ln[11:14].strip(), "author_start": a, "author_end": b, "ref_start": index.get((a, ln[25].strip())),
                                   "ref_end": index.get((b, ln[37].strip())), "length": int(ln[71:76])})
        if ln.startswith("SHEET ") and ln[21] == chain:
            a, b = int(ln[22:26]), int(ln[33:37])
            out["strands"].append({"id": f"{ln[11:14].strip()}/{ln[7:10].strip()}", "author_start": a, "author_end": b, "ref_start": index.get((a, ln[26].strip())),
                                   "ref_end": index.get((b, ln[37].strip()))})
    return out


def ss_of(ss, resseq):
    for h in ss["helices"]:
        if h["author_start"] <= resseq <= h["author_end"]:
            return f"helix {h['id']}"
    for s in ss["strands"]:
        if s["author_start"] <= resseq <= s["author_end"]:
            return f"strand {s['id']}"
    return "no deposited helix or strand"


def positions_windows(rec):
    pos = [f"{p['A_residue']}|{p['B_residue']}" for p in rec["positions"]]
    n = len(pos)
    w = {}
    for j in range(1, n - 1):
        if rec["metrics"]["d_kappa"]["values"][j] is not None:
            w[("3", tuple(pos[k] for k in range(j - 1, j + 2)))] = j
    for j in range(1, n - 2):
        if rec["metrics"]["d_tau"]["values"][j] is not None:
            w[("4", tuple(pos[k] for k in range(j - 1, j + 3)))] = j
    return w


def rank_of(vals, j):
    order = sorted((i for i, v in enumerate(vals) if v is not None), key=lambda i: (-vals[i], i))
    return order.index(j) + 1, len(order)


def main():
    sel = json.loads((KP / "results/KIN-PAIR-1_selection.json").read_text())
    panel = {p["pair_id"]: p for p in sel["panel"]}
    res = {"id": "PIM1-STRUCTURAL-ANNOTATION-2026-10-26", "design": "DESIGN.md (anchors and region-label rules fixed before this ran)",
           "scope": "interpretation of existing comparison outputs; no score recomputed, no threshold changed, no structure selected; these pairs do not isolate phosphorylation, ligand or mutation effects"}
    # reference and primary-source annotation
    raw = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{ACC}.json", timeout=60).read()
    u = json.loads(raw)
    fasta = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{ACC}.fasta", timeout=60).read().decode()
    ref = u["sequence"]["value"]
    assert "".join(fasta.splitlines()[1:]) == ref
    rec_hash = sel["targets"]["PIM1"]["uniprot_fasta_sha256"]
    feats = {}
    for f in u["features"]:
        if f["type"] in ("Domain", "Active site", "Binding site", "Modified residue"):
            feats.setdefault(f["type"], []).append({"start": f["location"]["start"]["value"], "end": f["location"]["end"]["value"], "description": f.get("description", ""),
                                                     "evidence": [e.get("source", "") + ":" + e.get("id", "") for e in f.get("evidences", [])]})
    res["reference"] = {"accession": ACC, "isoform": ISO + " (canonical)", "entry": u["uniProtkbId"], "organism": u["organism"]["scientificName"], "length": len(ref),
                        "sequence_sha256": sha(ref.encode()), "fasta_sha256_now": sha(fasta.encode()), "fasta_sha256_in_selection_record": rec_hash,
                        "matches_selection_record": sha(fasta.encode()) == rec_hash,
                        "numbering_convention": "reference index = UniProt canonical position (Numbering.reference_index); author numbering is preserved and its offset reported, never applied"}
    res["uniprot_annotation"] = {"url": f"https://rest.uniprot.org/uniprotkb/{ACC}.json", "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                                 "raw_sha256": sha(raw), "entry_version": u["entryAudit"]["entryVersion"], "last_annotation_update": u["entryAudit"]["lastAnnotationUpdateDate"],
                                 "features": feats, "note": "modified-residue and binding-site annotations are literature or proteomics curated annotations, not evidence of a functional state in these structures"}
    # structures
    P, text, scan_, maps, ss = {}, {}, {}, {}, {}
    res["structures"] = {}
    for e in ENTRIES:
        pair = next(p for p in panel.values() if e in (p["A"], p["B"]))
        k = "A" if pair["A"] == e else "B"
        f = KP / "data_pdb" / f"{e}.pdb"
        raw_f = f.read_bytes()
        assert sha(raw_f) == pair[f"sha256_{k}"], e
        text[e] = raw_f.decode("latin-1")
        P[e] = kp.parse_pdb(text[e], "blank_then_A")
        scan_[e] = inv.scan(text[e])[0]["A"]
        maps[e] = kp.map_chain(P[e], "A", ref, ACC)
        ss[e] = helices_strands(text[e], "A", maps[e]["index"])
        offs = sorted({k_[0] - v for k_, v in maps[e]["index"].items() if k_[1] == ""})
        res["structures"][e] = {"sha256": sha(raw_f), "pair_slot": ENTRIES[e], "chain": "A", "identity": maps[e]["identity"], "n_aligned": maps[e]["n_aligned"],
                                "author_minus_reference_offsets": offs, "mismatches": maps[e]["mismatches"], "dbref": maps[e]["dbref_segments"], "mapping_problems": maps[e]["problems"],
                                "deposited_secondary_structure": ss[e]}
    # alpha-C Glu contact rule (DESIGN.md)
    cands = []
    for idx in range(WIN_LO, WIN_HI + 1):
        if ref[idx - 1] != "E":
            continue
        dists, ok = {}, True
        for e in ENTRIES:
            inv_map = {v: k_ for k_, v in maps[e]["index"].items()}
            glu = next((r for r in scan_[e] if (r["resseq"], r["icode"]) == inv_map.get(idx)), None)
            lys = next((r for r in scan_[e] if (r["resseq"], r["icode"]) == inv_map.get(BETA3)), None)
            nz = inv.pick(lys["atoms"].get("NZ")) if lys else None
            d = [float(np.linalg.norm(inv.pick(glu["atoms"].get(o)) - nz)) for o in ("OE1", "OE2") if glu and nz is not None and inv.pick(glu["atoms"].get(o)) is not None]
            dists[e] = min(d) if d else None
            ok = ok and dists[e] is not None and dists[e] <= CONTACT
        cands.append({"reference_index": idx, "distances_angstrom": dists, "qualifies_in_all_four": ok})
    qual = [c["reference_index"] for c in cands if c["qualifies_in_all_four"]]
    res["alphaC_glu_rule"] = {"window": [WIN_LO, WIN_HI], "contact_angstrom": CONTACT, "glu_candidates": cands, "qualifying": qual, "established": len(qual) == 1}
    dfg_hits = [i for i in range(DFG_WINDOW[0], DFG_WINDOW[1] - 1) if ref[i - 1:i + 2] == "DFG"]
    res["dfg_window"] = {"window": list(DFG_WINDOW), "hits": dfg_hits}
    helper = {"ran": False}
    if len(qual) == 1 and len(dfg_hits) == 1:
        cfg = {"accession": ACC, "isoform": ISO, "sequence": ref, "source": f"UniProt {ACC} entry version {u['entryAudit']['entryVersion']}",
               "files": {e: str(KP / "data_pdb" / f"{e}.pdb") for e in ENTRIES},
               "anchors": {"beta3_Lys": BETA3, "alphaC_Glu": qual[0], "HRD_H": HRD_H, "DFG_Asp": list(DFG_WINDOW)}}
        (HERE / "work").mkdir(exist_ok=True)
        (HERE / "work/roles_config.json").write_text(json.dumps(cfg))
        r = subprocess.run([ADAPTER, str(HERE / "run_roles.py"), str(HERE / "work/roles_config.json"), str(HERE / "results/kinase_roles_reports.json")], capture_output=True, text=True,
                           env=dict(os.environ, PYTHONPATH=str(ROOT / "python")))
        helper = {"ran": True, "returncode": r.returncode, "stderr": r.stderr[-400:], "anchors_declared": cfg["anchors"], "allow_role_mutations": False}
        if r.returncode == 0:
            rep = json.loads((HERE / "results/kinase_roles_reports.json").read_text())
            helper["per_structure"] = {e: {"accepted": v["accepted"], "problems": v.get("problems", []),
                                           "anchor_checks": v["report"].get("anchor_checks"), "anchor_order": v["report"].get("anchor_order"),
                                           "roles": v["report"].get("roles"), "n_ambiguity_checks": len(v["report"].get("ambiguity_checks", {})),
                                           "all_unique": all(a["unique"] for a in v["report"].get("ambiguity_checks", {}).values()),
                                           "alignment": {k_: v["report"]["alignment"].get(k_) for k_ in ("identity", "aligned_pairs", "n_variants")} if "alignment" in v["report"] else None,
                                           "variants_in_chain": v["report"].get("alignment", {}).get("variants")} for e, v in rep.items()}
            helper["all_accepted"] = all(v["accepted"] for v in rep.values())
    else:
        helper["reason_not_run"] = f"alpha-C Glu rule qualifying={qual}, DFG hits={dfg_hits}"
    res["kinase_roles_helper"] = helper
    anchors_ok = helper.get("all_accepted", False)
    # landmarks (verified only)
    dom = feats["Domain"][0]
    landmarks = [{"name": "UniProt protein kinase domain", "start": dom["start"], "end": dom["end"], "source": "UniProt P11309 feature (PROSITE PRU00159)", "kind": "range"},
                 {"name": "active site Asp (UniProt, proton acceptor); HRD Asp", "start": 167, "end": 167, "source": "UniProt active site; helper checks HRD at 165-167", "kind": "point"},
                 {"name": "beta3 Lys (UniProt ATP binding site)", "start": BETA3, "end": BETA3, "source": "UniProt binding site; helper anchor", "kind": "point"}]
    if anchors_ok:
        roles = helper["per_structure"][next(iter(ENTRIES))]["roles"]
        landmarks += [{"name": "alpha-C Glu anchor", "start": qual[0], "end": qual[0], "source": "contact rule (4 structures) + helper motif check", "kind": "point"},
                      {"name": "DFG Asp anchor (motif DFG)", "start": dfg_hits[0], "end": dfg_hits[0], "source": "unique DFG in UniProt-bounded window; helper", "kind": "point"}]
    res["landmarks_verified"] = landmarks
    res["unverified"] = ["activation loop extent", "N-lobe/C-lobe boundary", "conventional names of the helices after the DFG anchor", "boundaries of every secondary-structure element beyond the deposited records",
                         "any functional role of residue 261 (UniProt lists it as a phosphoserine modified residue; the annotation is not evidence of a state here)"]
    # annotate the existing windows
    res["pairs"] = {}
    old_class = {p["pair_id"]: p["classification"] for p in json.loads((KP / "results/KIN-PAIR-1_results.json").read_text())["pairs"] if "classification" in p}
    for pid in ("PIM1-1", "PIM1-2"):
        new = json.loads((PH / f"runs/{pid}/new/numeric_record.json").read_text())
        old = json.loads((PH / f"runs/{pid}/old_with_numeric_record/numeric_record.json").read_text())
        wn, wo = positions_windows(new), positions_windows(old)
        added = sorted((k for k in wn if k not in wo), key=lambda k: (k[0], wn[k]))
        p = panel[pid]
        ea, eb = p["A"], p["B"]
        wcsv = list(csv.DictReader(open(PH / f"runs/{pid}/new/windows.csv")))
        top_prev = json.loads((PH / "phospho_replay.json").read_text())["affected"][pid]["aggregates"]["top10_changes"]
        rows = []
        for size, ident in added:
            j = wn[(size, ident)]
            metrics = M3 if size == "3" else M4
            residues = []
            refs = []
            for pidstr in ident:
                a_, b_ = pidstr.split("|")
                ka = (int(a_.split(":")[1].rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")), "")
                kb = (int(b_.split(":")[1].rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")), "")
                ia, ib = maps[ea]["index"].get(ka), maps[eb]["index"].get(kb)
                pos = next(q for q in new["positions"] if q["A_residue"] == a_ and q["B_residue"] == b_)
                residues.append({"A": a_, "B": b_, "reference_index_A": ia, "reference_index_B": ib, "A_resname": pos["A_resname"], "B_resname": pos["B_resname"],
                                 "A_modification": pos.get("A_modification"), "B_modification": pos.get("B_modification"),
                                 "deposited_ss_A": ss_of(ss[ea], ka[0]), "deposited_ss_B": ss_of(ss[eb], kb[0])})
                refs.append(ia)
            row = {"window_residues": size, "position_index": j, "residues": residues, "reference_indices": refs, "same_reference_index_in_both": all(r["reference_index_A"] == r["reference_index_B"] for r in residues),
                   "inside_uniprot_kinase_domain": all(dom["start"] <= x <= dom["end"] for x in refs),
                   "relative_to_dfg_anchor": (f"C-terminal to the DFG anchor (reference {dfg_hits[0]}); offsets {[x - dfg_hits[0] for x in refs]}" if anchors_ok else "unverified: anchors not mapped"),
                   "values": {m: new["metrics"][m]["values"][j] for m in metrics}, "units": {m: new["metrics"][m]["unit"] for m in metrics},
                   "ranks_in_new_run_of_windows_with_values": {m: dict(zip(("rank", "of"), rank_of(new["metrics"][m]["values"], j))) for m in metrics},
                   "name_difference_flag": wcsv[j][f"mismatch_in_window{size}"] == "True", "contains_modified_residue": wcsv[j][f"modified_in_window{size}"] == "True"}
            row["in_a_top10_of_any_metric"] = any(row["ranks_in_new_run_of_windows_with_values"][m]["rank"] <= 10 for m in metrics)
            rows.append(row)
        # descriptive classification: the recorded rule applied to the existing numeric values
        def winlists(rec):
            return {m: rec["metrics"][m]["values"] for m in M3 + M4}
        cn, co = kp.classify_pair(winlists(new)), kp.classify_pair(winlists(old))
        slim = lambda c: {s: {"class": v["class"], "overlap_O": v.get("overlap_O"), "n_valid": v["n_valid"]} for s, v in c.items()}   # noqa: E731
        core = lambda c: {s: (v["class"], v.get("overlap_O")) for s, v in c.items()}   # noqa: E731  (the conclusion: class and overlap O, not the number of windows)
        res["pairs"][pid] = {"entries": [ea, eb], "n_new_windows": {"3": sum(r["window_residues"] == "3" for r in rows), "4": sum(r["window_residues"] == "4" for r in rows)}, "new_windows": rows,
                             "descriptive": {"recorded_in_KIN-PAIR-1": slim(old_class[pid]), "reapplied_to_old_numeric_record": slim(co), "reapplied_to_new_numeric_record": slim(cn),
                                             "old_reapplication_reproduces_recorded": core(co) == core(old_class[pid]), "classification_changed_by_the_added_windows": core(cn) != core(co),
                                             "top10_old_still_in_top10_and_no_new_entrant": all(t["top10_old_in_new_top10"] == 10 and not t["top10_new_not_in_old_top10"] for t in top_prev.values()),
                                             "not_recomputed": "maximum trimmed-fit displacement and its location; global RMSD changed only by the refit with one more position (recorded in the phospho replay)"}}
    (HERE / "results/pim1_annotation.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps({"alphaC": res["alphaC_glu_rule"]["qualifying"], "dfg": dfg_hits, "helper_accepted": helper.get("all_accepted"), "helper_ran": helper["ran"],
                      "n_new": {k: v["n_new_windows"] for k, v in res["pairs"].items()}}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
