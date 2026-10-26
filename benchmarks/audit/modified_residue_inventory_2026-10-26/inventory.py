"""Modified-residue coverage inventory of the KIN-PAIR-1 inputs (read only; no production code, no scoring).

    .venv/bin/python benchmarks/audit/modified_residue_inventory_2026-10-26/inventory.py

Question: do the installed adapter and the study's correspondence rules omit observed non-standard (including phosphorylated) amino
acids, and which 3- and 4-residue windows does that prevent us from evaluating?

Inputs: the 12 original PDB files of the six KIN-PAIR-1 pairs (`kin_pair_1/data_pdb/`, SHA-256 verified against the selection record) and the
study's own outputs (`results/pairs/` or `results/pairs_supplement/`). The three UniProt FASTA files are re-read and checked against the
hashes in the selection record (they are the study's reference sequences, not new data). Residues are found from the coordinate records,
with MODRES used as annotation only. Inclusion by the installed adapter is taken from running the real adapter command (no selection,
no correspondence) in a diagnostic. No window metric is computed here; windows are counted from the study's recorded exclusion rows and
from continuity of the deposited coordinates. Output: results/modified_residue_inventory.json.
"""

import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
KP = ROOT / "benchmarks/audit/kin_pair_1"
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(KP))
import kinpair as kp  # noqa: E402  (frozen helpers, read only: parsing, UniProt mapping, correspondence)

ADAPTER = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
STD20 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
PEPTIDE_BOND_MAX = 2.0
CA_CA_MAX = 4.2
OUT = HERE / "results/modified_residue_inventory.json"
WORK = HERE / "work"


def sha(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


# --- coordinate-record scan (independent of the adapter and of MODRES) -------------------------------------------------

def scan(text: str):
    """First model only. Per chain, the residues in file order with their record type, names, and N/CA/C/O coordinates
    chosen by the study's altloc rule (blank, else A)."""
    model_records = sum(1 for ln in text.splitlines() if ln.startswith("MODEL "))
    chains, index, ended = {}, {}, False
    for ln in text.splitlines():
        if ln.startswith("ENDMDL"):
            ended = True
        if ended or ln[:6] not in ("ATOM  ", "HETATM"):
            continue
        key = (ln[21], int(ln[22:26]), ln[26].strip())
        r = index.get(key)
        if r is None:
            r = index[key] = {"chain": key[0], "resseq": key[1], "icode": key[2], "resname": ln[17:20].strip(), "record": ln[:6].strip(),
                              "atoms": {}, "altloc_labels": set()}
            chains.setdefault(key[0], []).append(r)
        if ln[17:20].strip() != r["resname"]:
            r.setdefault("other_names", set()).add(ln[17:20].strip())
            continue
        name = ln[12:16].strip()
        alt = ln[16].strip()
        if alt:
            r["altloc_labels"].add(alt)
        r["atoms"].setdefault(name, []).append((alt, np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])))
    return chains, model_records


def pick(recs):
    for want in ("", "A"):
        for alt, xyz in recs or []:
            if alt == want:
                return xyz
    return None


def bonded(prev, nxt):
    c, n = pick(prev["atoms"].get("C")), pick(nxt["atoms"].get("N"))
    if c is None or n is None:
        return None
    return float(np.linalg.norm(c - n))


def parse_header(text: str):
    modres, remark465 = {}, set()
    in465 = False
    for ln in text.splitlines():
        if ln.startswith("MODRES"):
            modres[(ln[16], int(ln[18:22]), ln[22].strip())] = {"resname": ln[12:15].strip(), "parent": ln[24:27].strip(), "comment": ln[29:70].strip()}
        if ln.startswith("REMARK 465"):
            m = re.match(r"REMARK 465\s+([A-Z0-9]{3})\s+([A-Za-z0-9])\s+(-?\d+)\s*([A-Za-z]?)\s*$", ln)
            if m:
                remark465.add((m.group(2), int(m.group(3)), m.group(4)))
    return modres, remark465


# --- the installed adapter, as a diagnostic ---------------------------------------------------------------------------

def adapter_residues(path: Path, tag: str):
    out = WORK / tag
    shutil.rmtree(out, ignore_errors=True)
    out.mkdir(parents=True)
    r = subprocess.run([ADAPTER, "-m", "topofold_inputs", str(path), "--run-id", tag, "--out", str(out), "--altloc", "blank_then_A"],
                       capture_output=True, text=True, env=dict(os.environ, PYTHONPATH=str(ROOT / "python")))
    if r.returncode != 0:
        raise RuntimeError(f"adapter diagnostic failed for {path.name}: {r.stderr[-400:]}")
    rep = json.loads((out / f"{tag}.json").read_text())
    return {(x["chain"], x["resid"], x.get("icode") or ""): {"resname": x["resname"], "status": x["status"]} for x in rep["residues"]}, \
        rep["provenance"]["software"], rep["provenance"]["adapter"]


def has_backbone(r):
    return all(pick(r["atoms"].get(a)) is not None for a in ("N", "CA", "C"))


def chemcomp(code, cache={}):
    if code not in cache:
        try:
            d = json.load(urllib.request.urlopen(f"https://data.rcsb.org/rest/v1/core/chemcomp/{code}", timeout=30))
            cache[code] = {"name": d["chem_comp"].get("name"), "type": d["chem_comp"].get("type"), "parent": d["chem_comp"].get("mon_nstd_parent_comp_id")}
        except Exception as e:     # noqa: BLE001
            cache[code] = {"error": str(e)[:80]}
    return cache[code]


# --- main ------------------------------------------------------------------------------------------------------------

def main():
    sel = json.loads((KP / "results/KIN-PAIR-1_selection.json").read_text())
    panel = sel["panel"]
    files = {}
    for p in panel:
        for k in "AB":
            files.setdefault(p[k], {"sha256": p[f"sha256_{k}"], "chain": p[f"chain_{k}"], "target_chains": p[f"meta_{k}"]["target_chains"], "pairs": []})["pairs"].append(f"{p['pair_id']}:{k}")
    by_hash = {}
    for e, f in files.items():
        raw = (KP / "data_pdb" / f"{e}.pdb").read_bytes()
        assert sha(raw) == f["sha256"], f"{e}: hash differs from the selection record"
        by_hash.setdefault(f["sha256"], []).append(e)
    assert all(len(v) == 1 for v in by_hash.values())
    refs = {}
    for g in {p["gene"] for p in panel}:
        acc = next(p["accession"] for p in panel if p["gene"] == g)
        txt = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=60).read().decode()
        assert sha(txt.encode()) == sel["targets"][g]["uniprot_fasta_sha256"], f"{g}: UniProt FASTA differs from the selection record"
        refs[acc] = "".join(txt.splitlines()[1:])
    parsed, scanned, header, adapter, soft, ADAPTER_META = {}, {}, {}, {}, None, None

    def npoly(e, chain):
        """Residue records with N, CA and C atoms (amino-acid-like), excluding waters, ions and ligands."""
        return sum(1 for r in scanned[e][chain] if has_backbone(r))
    for e, f in files.items():
        text = (KP / "data_pdb" / f"{e}.pdb").read_text(encoding="latin-1")
        parsed[e] = kp.parse_pdb(text, "blank_then_A")
        scanned[e], f["model_records"] = scan(text)
        header[e] = parse_header(text)
        adapter[e], soft, ADAPTER_META = adapter_residues(KP / "data_pdb" / f"{e}.pdb", f"diag_{e}")
    # study pairs: mapping, correspondence, outputs of record
    pair_rows, residues_out = [], []
    corr_keys = {}                      # (entry, chain) -> {key: (pair_id, uniprot index)}
    pair_data = {}
    for p in panel:
        pid, acc = p["pair_id"], p["accession"]
        ref = refs[acc]
        mp = {k: kp.map_chain(parsed[p[k]], p[f"chain_{k}"], ref, acc) for k in "AB"}
        corr = kp.correspondence(parsed[p["A"]], mp["A"], p["chain_A"], parsed[p["B"]], mp["B"], p["chain_B"])
        out_dir = KP / "results/pairs" / pid
        if not (out_dir / "residues.csv").is_file():
            out_dir = KP / "results/pairs_supplement" / pid
        rows = list(csv.DictReader(open(out_dir / "residues.csv")))
        mine = [(f"{p['chain_A']}:{a[0]}{a[1]}", f"{p['chain_B']}:{b[0]}{b[1]}") for a, b, _ in corr["pairs"]]
        assert mine == [(r["A_residue"], r["B_residue"]) for r in rows], f"{pid}: recomputed correspondence differs from the study's residues.csv"
        pair_data[pid] = {"p": p, "mp": mp, "corr": corr, "dir": out_dir, "rows": rows}
        for k, side in (("A", 0), ("B", 1)):
            corr_keys.setdefault((p[k], p[f"chain_{k}"]), {})[pid] = {pr[side]: pr[2] for pr in corr["pairs"]}
    # relevant residues: non-standard names with backbone atoms, or MODRES-annotated, in the selected chains
    names_seen = set()
    for e, f in sorted(files.items()):
        chain = f["chain"]
        modres, rem465 = header[e]
        res = scanned[e][chain]
        for i, r in enumerate(res):
            has_bb = all(pick(r["atoms"].get(a)) is not None for a in ("N", "CA", "C"))
            key = (chain, r["resseq"], r["icode"])
            ann = modres.get(key)
            if r["resname"] in STD20 or not (has_bb or ann):
                continue
            names_seen.add(r["resname"])
            prev = res[i - 1] if i else None
            nxt = res[i + 1] if i + 1 < len(res) else None
            d_prev = bonded(prev, r) if prev else None
            d_next = bonded(r, nxt) if nxt else None
            member = (d_prev is not None and d_prev <= PEPTIDE_BOND_MAX) or (d_next is not None and d_next <= PEPTIDE_BOND_MAX) or bool(ann)
            ad = adapter[e].get(key)
            studies = []
            for pid, pdx in pair_data.items():
                for k in "AB":
                    if pdx["p"][k] == e and pdx["p"][f"chain_{k}"] == chain:
                        incl = key in corr_keys[(e, chain)][pid]
                        studies.append({"pair_id": pid, "side": k, "in_correspondence": incl,
                                        "uniprot_index": corr_keys[(e, chain)][pid].get(key)})
            mapped = None
            for pid, pdx in pair_data.items():
                for k in "AB":
                    if pdx["p"][k] == e and pdx["p"][f"chain_{k}"] == chain:
                        mapped = pdx["mp"][k]["index"].get((r["resseq"], r["icode"]))
            reasons = []
            if not has_bb or pick(r["atoms"].get("CA")) is None:
                category, reasons = 1, ["no Cα coordinate in the records (partial or absent coordinates)"]
            elif not member:
                category, reasons = 4, ["has backbone atoms but is not peptide-bonded to a neighbour and has no MODRES record: polymer membership ambiguous"]
            elif ad is None:
                category, reasons = 2, ["observed in the coordinate records; the installed adapter's `protein` selection does not contain this residue name, so it is absent from the adapter output"]
            elif mapped is None:
                category, reasons = 4, ["read by the adapter but not aligned to the reference sequence: mapping ambiguous"]
            elif not any(s["in_correspondence"] for s in studies):
                category, reasons = 3, ["read by the adapter, aligned to the reference, but not in the study correspondence (position not shared with the partner)"]
            else:
                category, reasons = 0, ["included"]
            cc = chemcomp(r["resname"])
            residues_out.append({
                "input_id": e, "sha256": f["sha256"], "model": 1, "model_records_in_file": f["model_records"], "chain": chain, "resseq": r["resseq"], "icode": r["icode"],
                "resname": r["resname"], "record_type": r["record"], "modres": ann, "modres_annotated": ann is not None,
                "chem_comp": cc, "parent_documented": (ann or {}).get("parent") or cc.get("parent"),
                "atoms_observed": {a: pick(r["atoms"].get(a)) is not None for a in ("N", "CA", "C", "O")},
                "altloc_labels": sorted(r["altloc_labels"]), "neighbour_backbone": {
                    "previous_residue_C_observed": prev is not None and pick(prev["atoms"].get("C")) is not None,
                    "next_residue_N_observed": nxt is not None and pick(nxt["atoms"].get("N")) is not None,
                    "C_prev_to_N_angstrom": d_prev, "C_to_N_next_angstrom": d_next},
                "polymer_member": member, "adapter_included": ad is not None, "adapter_status": None if ad is None else ad["status"],
                "adapter_resname_seen": None if ad is None else ad["resname"], "uniprot_index": mapped, "study_inclusion": studies,
                "category": category, "category_label": {0: "included", 1: "genuinely unobserved or partial coordinates", 2: "observed residue excluded by the parser",
                                                         3: "excluded by the study specification", 4: "ambiguous mapping or polymer membership"}[category],
                "reason": "; ".join(reasons), "in_remark_465": key in rem465,
                "phosphorylated_by_documented_component": bool(re.search("phospho", (cc.get("name") or ""), re.I)),
                "region": "unverified: no verified kinase role mapping exists for this target; deposited number and reference index are reported instead"})
    # MODRES-annotated residues with no coordinate record at all (in the selected chains)
    unobserved_annotated = []
    for e, f in sorted(files.items()):
        modres, _ = header[e]
        have = {(r["chain"], r["resseq"], r["icode"]) for r in scanned[e].get(f["chain"], [])}
        for key, a in modres.items():
            if key[0] == f["chain"] and key not in have:
                unobserved_annotated.append({"input_id": e, "sha256": f["sha256"], "chain": key[0], "resseq": key[1], "icode": key[2], **a,
                                             "category": 1, "reason": "annotated in MODRES but no ATOM/HETATM coordinate record"})
    # windows: gaps between consecutive correspondence positions, classified, and the windows they exclude
    for pid, pdx in pair_data.items():
        p, mp, corr = pdx["p"], pdx["mp"], pdx["corr"]
        ra = {r.key: r for r in parsed[p["A"]].residues[p["chain_A"]]}
        rb = {r.key: r for r in parsed[p["B"]].residues[p["chain_B"]]}
        pairs = corr["pairs"]
        n = len(pairs)
        inv = {k: {idx: key for key, idx in mp[k]["index"].items()} for k in "AB"}
        win_rows = list(csv.DictReader(open(pdx["dir"] / "windows.csv")))
        excl = list(csv.DictReader(open(pdx["dir"] / "exclusions.csv")))
        ex3 = {int(r["position"]) for r in excl if r["window_size"] == "3"}
        ex4 = {int(r["position"]) for r in excl if r["window_size"] == "4"}
        valid3 = sum(1 for r in win_rows if r["d_kappa_rad"] != "")
        valid4 = sum(1 for r in win_rows if r["d_tau_rad"] != "")
        gaps = []
        for j in range(n - 1):
            (a0, b0, u0), (a1, b1, u1) = pairs[j], pairs[j + 1]
            broken = not (kp.consecutive(ra[a0], ra[a1]) and kp.consecutive(rb[b0], rb[b1])) or \
                float(np.linalg.norm(ra[a0].ca - ra[a1].ca)) > CA_CA_MAX or float(np.linalg.norm(rb[b0].ca - rb[b1].ca)) > CA_CA_MAX
            if not broken:
                continue
            skipped = []
            for u in range(u0 + 1, u1):
                row = {"uniprot_index": u}
                for k, R, src in (("A", ra, "A"), ("B", rb, "B")):
                    key = inv[k].get(u)
                    if key is None:
                        row[k] = {"state": "no residue aligned to this reference index in the deposited chain", "remark_465": None}
                    else:
                        r_ = R[key]
                        other = (pairs[j][1 if k == "A" else 0])
                        if r_.ca is None:
                            row[k] = {"deposited": f"{p['chain_' + k]}:{key[0]}{key[1]}", "resname": r_.resname, "state": "residue record without Cα"}
                        elif not r_.read:
                            row[k] = {"deposited": f"{p['chain_' + k]}:{key[0]}{key[1]}", "resname": r_.resname, "state": "observed, not read by the adapter (non-standard residue name)",
                                      "modified": True}
                        else:
                            row[k] = {"deposited": f"{p['chain_' + k]}:{key[0]}{key[1]}", "resname": r_.resname, "state": "observed and read by the adapter in this chain"}
                skipped.append(row)
            # recoverable if every skipped index would be observed with Cα in both chains and the links continuous
            all_present = all(s["A"].get("state", "").startswith(("observed and read", "observed, not read")) and s["B"].get("state", "").startswith(("observed and read", "observed, not read")) for s in skipped) and bool(skipped)
            has_modified = any(s[k].get("modified") for s in skipped for k in "AB")
            only_modified = bool(skipped) and all(any(s[k].get("modified") for k in "AB") and all(s[k].get("state", "").startswith(("observed and read", "observed, not read")) for k in "AB") for s in skipped)
            gaps.append({"between_positions": [j, j + 1], "between_uniprot_indices": [u0, u1], "between_deposited_A": [f"{p['chain_A']}:{a0[0]}{a0[1]}", f"{p['chain_A']}:{a1[0]}{a1[1]}"],
                         "between_deposited_B": [f"{p['chain_B']}:{b0[0]}{b0[1]}", f"{p['chain_B']}:{b1[0]}{b1[1]}"], "skipped": skipped,
                         "contains_observed_unread_modified_residue": has_modified, "all_skipped_residues_observed_in_both_chains": all_present,
                         "gap_consists_only_of_observed_modified_residues": only_modified,
                         "windows3_centres": [c for c in (j, j + 1) if 1 <= c <= n - 2], "windows4_indices": [c for c in (j - 1, j, j + 1) if 1 <= c <= n - 3]})
        attr3 = sorted({c for g in gaps if g["contains_observed_unread_modified_residue"] for c in g["windows3_centres"]})
        attr4 = sorted({c for g in gaps if g["contains_observed_unread_modified_residue"] for c in g["windows4_indices"]})
        only3 = sorted({c for g in gaps if g["gap_consists_only_of_observed_modified_residues"] for c in g["windows3_centres"]})
        only4 = sorted({c for g in gaps if g["gap_consists_only_of_observed_modified_residues"] for c in g["windows4_indices"]})
        union3 = sorted({c for g in gaps for c in g["windows3_centres"]})
        union4 = sorted({c for g in gaps for c in g["windows4_indices"]})
        pair_rows.append({"pair_id": pid, "entries": [p["A"], p["B"]], "n_positions": n, "windows_with_values": {"3": valid3, "4": valid4},
                          "windows_excluded_by_study_output": {"3": len(ex3), "4": len(ex4)},
                          "windows_in_gap_links_recomputed": {"3": len(union3), "4": len(union4)},
                          "gap_links_match_study_exclusions": set(union3) == ex3 and set(union4) == ex4,
                          "windows_touching_a_gap_with_an_observed_unread_modified_residue": {"3": len(attr3), "4": len(attr4), "positions3": attr3, "positions4": attr4},
                          "windows_in_gaps_made_only_of_observed_modified_residues": {"3": len(only3), "4": len(only4), "positions3": only3, "positions4": only4},
                          "n_gaps": len(gaps), "gaps": gaps, "study_outputs": str(pdx["dir"].relative_to(ROOT))})
    # other polymer chains of the same files (counts only; the study used the selected chains)
    other = {}
    for e, f in sorted(files.items()):
        rows = []
        for ch, res in sorted(scanned[e].items()):
            if ch == f["chain"]:
                continue
            nonstd = [r for r in res if r["resname"] not in STD20 and all(pick(r["atoms"].get(a)) is not None for a in ("N", "CA", "C"))]
            n_std = sum(r["resname"] in STD20 for r in res)
            if n_std >= 20 or nonstd:
                rows.append({"chain": ch, "n_standard_residue_records": n_std,
                             "non_standard_residues_with_backbone": [{"resname": r["resname"], "resseq": r["resseq"], "icode": r["icode"], "record": r["record"],
                                                                      "modres_annotated": (ch, r["resseq"], r["icode"]) in header[e][0]} for r in nonstd]})
        other[e] = rows
    # collection-level counts
    rel = residues_out
    by_cat = {}
    for r in rel:
        by_cat[r["category"]] = by_cat.get(r["category"], 0) + 1
    obs_unread = [r for r in rel if r["category"] == 2]
    lost3 = sum(p["windows_touching_a_gap_with_an_observed_unread_modified_residue"]["3"] for p in pair_rows)
    lost4 = sum(p["windows_touching_a_gap_with_an_observed_unread_modified_residue"]["4"] for p in pair_rows)
    n_res_selected = sum(npoly(e, f["chain"]) for e, f in files.items())
    if obs_unread and (lost3 + lost4) > 0:
        rec = "A narrowly scoped modified-residue support task is warranted."
        rec_code = "narrow_support_task"
    elif not rel or all(r["category"] in (0, 1) for r in rel):
        rec = "No coverage problem is demonstrated in this collection."
        rec_code = "no_demonstrated_problem"
    else:
        rec = "Insufficient evidence."
        rec_code = "insufficient_evidence"
    out = {"id": "MODIFIED-RESIDUE-INVENTORY-2026-10-26", "scope": "the 12 original KIN-PAIR-1 input files (six pairs), selected chains only; frequencies apply to this local, selected collection and not to kinases generally",
           "no_new_scores": "no window metric was computed; windows are counted from the study's recorded exclusion rows and from deposited-coordinate continuity",
           "adapter": {"python": ADAPTER, "software": soft, "adapter": ADAPTER_META, "diagnostic": "installed adapter command run on each file with its `protein` selection, altloc blank_then_A, no correspondence"},
           "collection": {"n_pair_slots": 2 * len(panel), "n_files": len(files), "n_unique_sha256": len(by_hash),
                          "files": [{"input_id": e, "sha256": f["sha256"], "verified_against_selection_record": True, "selected_chain": f["chain"], "pairs": f["pairs"],
                                     "model_records_in_file": f["model_records"], "polymer_residues_in_selected_chain": npoly(e, f["chain"]),
                                     "residues_absent_from_adapter_output_in_selected_chain": sum(1 for r in scanned[e][f["chain"]] if (f["chain"], r["resseq"], r["icode"]) not in adapter[e] and r["resname"] in STD20 | {"SEP", "TPO", "PTR"} | set(names_seen)),
                                     "modres_records": len(header[e][0])} for e, f in sorted(files.items())],
                          "n_polymer_residue_records_in_selected_chains": n_res_selected},
           "relevant_residues": rel, "annotated_without_coordinates": unobserved_annotated, "other_chains_in_the_same_files": other,
           "counts": {"relevant_residues_by_category": {str(k): v for k, v in sorted(by_cat.items())}, "n_relevant_residues": len(rel),
                      "n_observed_unread_by_adapter": len(obs_unread), "n_with_modres": sum(r["modres_annotated"] for r in rel),
                      "n_without_modres_annotation": sum(not r["modres_annotated"] for r in rel), "names": sorted(names_seen),
                      "phosphorylated_by_documented_component": sorted({(r["input_id"], r["resname"], r["resseq"]) for r in rel if r["phosphorylated_by_documented_component"]}),
                      "windows_touching_gaps_with_observed_unread_modified_residue": {"3": lost3, "4": lost4},
                      "windows_in_gaps_made_only_of_observed_modified_residues": {"3": sum(p["windows_in_gaps_made_only_of_observed_modified_residues"]["3"] for p in pair_rows),
                                                                                    "4": sum(p["windows_in_gaps_made_only_of_observed_modified_residues"]["4"] for p in pair_rows)},
                      "windows_with_values_total": {"3": sum(p["windows_with_values"]["3"] for p in pair_rows), "4": sum(p["windows_with_values"]["4"] for p in pair_rows)},
                      "windows_excluded_by_study_output_total": {"3": sum(p["windows_excluded_by_study_output"]["3"] for p in pair_rows), "4": sum(p["windows_excluded_by_study_output"]["4"] for p in pair_rows)}},
           "pairs": pair_rows,
           "decision_rule": "narrow support task if at least one observed modified residue with backbone is excluded by the adapter and at least one window touches a gap that contains it; no problem if no modified residue is lost; otherwise insufficient evidence",
           "recommendation": {"code": rec_code, "text": rec}}
    OUT.write_text(json.dumps(out, indent=1, default=lambda o: sorted(o) if isinstance(o, set) else str(o)) + "\n")
    print(json.dumps({"counts": out["counts"], "recommendation": out["recommendation"]}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
