"""Descriptive demonstration of structpair 0.3.0 on the locally held Abl structures: 1OPJ chains A, B versus 1OPK chain A.

    .venv/bin/python benchmarks/audit/structpair_abl_demo_2026-10-25/run_demo.py [--adapter-python PATH]

Not an experiment and not evidence about imatinib, a nucleotide, a mutation or resistance: the two entries differ in species
numbering source, construct, inhibitor, an engineered mutation and crystal form (see the record). The correspondence rule is
fixed here before any comparison is computed:
  * positions are the PDB author numbers 248-534 shared by both entries (the 1OPJ construct range);
  * a position is mapped if it has an observed C-alpha in all three chains (1OPJ A, 1OPJ B, 1OPK A), from the coordinates and
    cross-checked against REMARK 465;
  * a residue-name difference between chains is accepted only if SEQADV of an entry documents it as an engineered mutation (then
    it stays mapped and is flagged); any other difference stops the demonstration;
  * no alternate location or insertion code may occur at a mapped C-alpha (the policy is `reject`).
Inputs are the local copies used by KIN-DESC-1 (SHA-256 checked against that record) and are cross-checked against RCSB entry
metadata. No download of coordinates, raw data or trajectories. Writes results/abl_demo_record.json and the compact outputs.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import resource
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = ROOT / "benchmarks/audit/kinase_descriptors/data_pdb"
KIN_SPEC = ROOT / "benchmarks/audit/kinase_descriptors/run_spec.json"
OUT = HERE / "results"
RANGE = (248, 534)
ENTRIES = {"A": ("1OPJ", ["A", "B"], "A"), "B": ("1OPK", ["A"], "A")}
KEEP = ("chain_variation.csv", "chains.csv", "chain_variation.png", "windows.csv", "residues.csv", "exclusions.csv", "report.md",
        "manifest.json", "summary.json")


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while b := f.read(1 << 20):
            h.update(b)
    return h.hexdigest()


def header_record(path: Path) -> dict:
    """Facts read from the deposited header records of one PDB file (no interpretation)."""
    lines = path.read_text().splitlines()
    rec = lambda tag: [ln for ln in lines if ln.startswith(tag)]  # noqa: E731
    compnd = " ".join(ln[10:].strip() for ln in rec("COMPND"))
    source = " ".join(ln[10:].strip() for ln in rec("SOURCE"))
    pick = lambda text, key: (re.search(key + r":\s*([^;]*);?", text) or [None, None])[1]  # noqa: E731
    seqadv = [{"resname": ln[12:15].strip(), "chain": ln[16], "resid": int(ln[18:22]), "uniprot_resname": ln[39:42].strip(),
               "uniprot_resid": ln[43:48].strip(), "comment": ln[49:].strip()} for ln in rec("SEQADV")]
    dbref = [{"chain": ln[12], "pdb_range": [int(ln[14:18]), int(ln[20:24])], "db": ln[26:32].strip(), "accession": ln[33:41].strip(),
              "db_id": ln[42:54].strip(), "db_range": [int(ln[55:60]), int(ln[62:67])]} for ln in rec("DBREF")]
    het = [{"name": ln[7:10].strip(), "chain": ln[12], "number": int(ln[13:17])} for ln in rec("HET   ")]
    hetnam = {}
    for ln in rec("HETNAM"):
        hetnam[ln[11:14].strip()] = (hetnam.get(ln[11:14].strip(), "") + " " + ln[15:].strip()).strip()
    hetsyn = {}
    for ln in rec("HETSYN"):
        hetsyn[ln[11:14].strip()] = (hetsyn.get(ln[11:14].strip(), "") + " " + ln[15:].strip()).strip()
    cryst = rec("CRYST1")[0]
    missing = {}
    for ln in lines:
        if ln.startswith("REMARK 465") and re.match(r"REMARK 465\s+[A-Z]{3}\s+[A-Z]\s+-?\d+", ln):
            f = ln.split()
            missing.setdefault(f[3], []).append(int(f[4]))
    return {"title": " ".join(ln[10:].strip() for ln in rec("TITLE")), "method": " ".join(ln[10:].strip() for ln in rec("EXPDTA")),
            "resolution_A": float(re.search(r"RESOLUTION\.\s+([\d.]+)", " ".join(rec("REMARK   2 RESOLUTION")))[1]),
            "molecule": pick(compnd, "MOLECULE"), "fragment": pick(compnd, "FRAGMENT"), "chains_declared": pick(compnd, "CHAIN"),
            "engineered": pick(compnd, "ENGINEERED"), "mutation_flag": pick(compnd, "MUTATION"),
            "organism": pick(source, "ORGANISM_SCIENTIFIC"), "taxid": pick(source, "ORGANISM_TAXID"), "gene": pick(source, "GENE"),
            "expression_system": pick(source, "EXPRESSION_SYSTEM"), "dbref": dbref, "seqadv": seqadv, "het": het, "hetnam": hetnam, "hetsyn": hetsyn,
            "space_group": cryst[55:66].strip(), "cell": cryst[6:54].split(), "remark_465_missing": missing,
            "doi": next((ln.split("DOI")[1].strip() for ln in rec("JRNL        DOI")), None)}


def calpha(path: Path) -> dict:
    """{chain: {resid: (resname, altloc, icode)}} of ATOM C-alpha records."""
    out = {}
    for ln in path.read_text().splitlines():
        if ln.startswith("ATOM  ") and ln[12:16] == " CA ":
            out.setdefault(ln[21], {})[int(ln[22:26])] = (ln[17:20].strip(), ln[16], ln[26])
    return out


def ranges(nums) -> str:
    out, cur = [], []
    for n in sorted(nums):
        if cur and n == cur[-1] + 1:
            cur.append(n)
        else:
            if cur:
                out.append(cur)
            cur = [n]
    if cur:
        out.append(cur)
    return ",".join(f"{c[0]}-{c[-1]}" if len(c) > 1 else str(c[0]) for c in out)


def rcsb(pid: str) -> dict:
    def get(url):
        for attempt in range(4):                                   # the RCSB API occasionally drops a connection
            try:
                with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "TopoFold-demo/1"}), timeout=60) as r:
                    return json.loads(r.read())
            except OSError:
                if attempt == 3:
                    raise
                time.sleep(2 * (attempt + 1))
    e = get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
    ent = get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/1")
    comps = sorted(get(f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pid}/{i}")["pdbx_entity_nonpoly"]["comp_id"]
                   for i in e["rcsb_entry_container_identifiers"].get("non_polymer_entity_ids", []))
    ids = ent["rcsb_polymer_entity_container_identifiers"]
    return {"title": e["struct"]["title"], "resolution_A": (e["rcsb_entry_info"].get("resolution_combined") or [None])[0],
            "space_group": e["symmetry"]["space_group_name_H_M"], "nonpolymer_comp_ids": comps,
            "nonpolymer_instances": e["rcsb_entry_info"].get("deposited_nonpolymer_entity_instance_count"),
            "organism": [o.get("ncbi_scientific_name") for o in ent.get("rcsb_entity_source_organism", [])],
            "uniprot": [r.get("database_accession") for r in ids.get("reference_sequence_identifiers", []) if r.get("database_name") == "UniProt"],
            "mutation_count": ent["entity_poly"].get("rcsb_mutation_count"), "chains": ids.get("auth_asym_ids")}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--adapter-python", default=os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python"))
    a = ap.parse_args(argv)
    kin = json.loads(KIN_SPEC.read_text())["structures"]
    inputs = {}
    for key, (pid, chains, prim) in ENTRIES.items():
        path = DATA / f"{pid}.pdb"
        sha = sha256_file(path)
        if sha != kin[pid]["sha256"]:
            raise SystemExit(f"STOP: {pid} differs from the KIN-DESC-1 record")
        hd = header_record(path)
        inputs[key] = {"pdb": pid, "path": str(path.relative_to(ROOT)), "sha256": sha, "sha256_equals_kin_desc_1_record": True,
                       "chains_selected": chains, "primary_chain": prim, "header": hd, "rcsb": rcsb(pid)}
    # ---- input checks and the correspondence (rule fixed in the docstring) ----
    ca = {pid: calpha(DATA / f"{pid}.pdb") for pid in ("1OPJ", "1OPK")}
    chains = [("1OPJ", "A"), ("1OPJ", "B"), ("1OPK", "A")]
    lo, hi = RANGE
    obs = {c: {r for r in ca[c[0]][c[1]] if lo <= r <= hi} for c in chains}
    mapped = sorted(set.intersection(*obs.values()))
    checks = {}
    miss465 = {c: set(inputs["A" if c[0] == "1OPJ" else "B"]["header"]["remark_465_missing"].get(c[1], [])) for c in chains}
    checks["remark_465_consistent_with_coordinates"] = all(not (miss465[c] & set(ca[c[0]][c[1]])) for c in chains)
    checks["no_altloc_or_icode_at_mapped_calpha"] = all(ca[p][ch][r][1] == " " and ca[p][ch][r][2] == " " for p, ch in chains for r in mapped)
    names = {r: {c: ca[c[0]][c[1]][r][0] for c in chains} for r in mapped}
    diff = [r for r in mapped if len(set(names[r].values())) > 1]
    documented = {(s["resid"]) for k in inputs.values() for s in k["header"]["seqadv"] if "ENGINEERED" in s["comment"].upper()}
    undocumented = [r for r in diff if r not in documented]
    checks["residue_name_differences"] = {str(r): {f"{p}:{ch}": n for (p, ch), n in names[r].items()} for r in diff}
    checks["differences_all_documented_as_engineered_in_SEQADV"] = not undocumented
    checks["same_species_and_gene"] = (inputs["A"]["header"]["organism"], inputs["A"]["header"]["gene"]) == (inputs["B"]["header"]["organism"], inputs["B"]["header"]["gene"])
    checks["same_uniprot_accession_in_DBREF"] = {inputs[k]["pdb"]: sorted({d["accession"] for d in inputs[k]["header"]["dbref"]}) for k in "AB"}
    checks["same_numbering_offset_to_database"] = all(d["pdb_range"][0] - d["db_range"][0] == 19 for k in "AB" for d in inputs[k]["header"]["dbref"])
    checks["rcsb_agrees_with_local_records"] = {inputs[k]["pdb"]: {
        "space_group": inputs[k]["rcsb"]["space_group"].replace(" ", "") == inputs[k]["header"]["space_group"].replace(" ", ""),
        "resolution": abs((inputs[k]["rcsb"]["resolution_A"] or 0) - inputs[k]["header"]["resolution_A"]) < 1e-6,
        "ligands": sorted({h["name"] for h in inputs[k]["header"]["het"]}) == inputs[k]["rcsb"]["nonpolymer_comp_ids"],
        "organism": [o.upper() for o in inputs[k]["rcsb"]["organism"]] == [inputs[k]["header"]["organism"]],
        "ligand_instances": inputs[k]["rcsb"]["nonpolymer_instances"] == len(inputs[k]["header"]["het"]),
        "uniprot": inputs[k]["rcsb"]["uniprot"] == sorted({d["accession"] for d in inputs[k]["header"]["dbref"]}),
        "engineered_mutation_count": inputs[k]["rcsb"]["mutation_count"] == sum("ENGINEERED" in x["comment"].upper() for x in inputs[k]["header"]["seqadv"]),
        "chains": sorted(inputs[k]["rcsb"]["chains"]) == sorted(inputs[k]["header"]["chains_declared"].replace(" ", "").split(","))} for k in "AB"}
    ok = (checks["remark_465_consistent_with_coordinates"] and checks["no_altloc_or_icode_at_mapped_calpha"] and checks["same_species_and_gene"]
          and checks["differences_all_documented_as_engineered_in_SEQADV"] and checks["same_numbering_offset_to_database"]
          and all(all(v.values()) for v in checks["rcsb_agrees_with_local_records"].values()))
    if not ok or len(mapped) < 100:
        OUT.mkdir(exist_ok=True)
        (OUT / "abl_demo_STOPPED.json").write_text(json.dumps({"inputs": inputs, "checks": checks, "mapped": len(mapped)}, indent=1, default=str) + "\n")
        raise SystemExit("STOP: an input check failed; the demonstration is not run (see results/abl_demo_STOPPED.json)")
    # ---- structpair spec and run ----
    spec = {"name": "Descriptive demonstration: 1OPJ chains A,B vs 1OPK chain A (mouse Abl1, different inhibitor, construct, engineered D363N, crystal form)",
            "structures": {"A": {"label": "1OPJ (imatinib; kinase domain)", "file": f"../../kinase_descriptors/data_pdb/1OPJ.pdb", "chains": ["A", "B"],
                                 "primary_chain": "A", "altloc": "reject"},
                           "B": {"label": "1OPK (PD166326; SH3-SH2-kinase)", "file": "../../kinase_descriptors/data_pdb/1OPK.pdb", "chains": ["A"],
                                 "primary_chain": "A", "altloc": "reject"}},
            "correspondence": {"same_numbering": ranges(mapped)}, "on_missing": "reject", "resname_mismatch": "flag", "top_n": 8}
    spec_path = HERE / "specs/abl_1OPJ_AB_vs_1OPK_A.json"
    spec_path.parent.mkdir(exist_ok=True)
    spec_path.write_text(json.dumps(spec, indent=1) + "\n")
    with tempfile.TemporaryDirectory(prefix="abl_demo_") as t:
        o = Path(t) / "out"
        t0, r0 = time.time(), resource.getrusage(resource.RUSAGE_CHILDREN)
        p = subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(spec_path), "--out", str(o)], capture_output=True, text=True,
                           env=dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=a.adapter_python))
        r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
        if p.returncode != 0:
            raise SystemExit(f"STOP: structpair failed: {p.stderr.strip()}")
        size = lambda d: sum(f.stat().st_size for f in d.rglob("*") if f.is_file())  # noqa: E731
        res = {"wall_s": round(time.time() - t0, 2), "cpu_s": round((r1.ru_utime + r1.ru_stime) - (r0.ru_utime + r0.ru_stime), 2),
               "output_bytes_total": size(o), "output_bytes_without_work": size(o) - size(o / "work")}
        dest = OUT / "run"
        if dest.exists():
            shutil.rmtree(dest)
        dest.mkdir(parents=True)
        for f in KEEP:
            shutil.copy(o / f, dest / f)
    return finish(spec, inputs, checks, mapped, res)


def finish(spec, inputs, checks, mapped, res) -> int:
    import csv
    dest = OUT / "run"
    rows = list(csv.DictReader(open(dest / "chain_variation.csv")))
    disp = [r for r in rows if r["channel"] == "disp"]
    f = lambda k: [float(r[k]) for r in disp if r[k] != ""]  # noqa: E731
    summ = json.loads((dest / "summary.json").read_text())["within_entry_variation"]
    man = json.loads((dest / "manifest.json").read_text())
    stats = {"positions_mapped": len(mapped), "mapped_ranges": spec["correspondence"]["same_numbering"],
             "within_entry_status": summ["status"], "n_pairs": summ["n_pairs"], "pair_counts_per_position": sorted({(r["between_n"], r["withinA_n"], r["withinB_n"]) for r in disp}),
             "disp_between_median_over_positions": {"median": statistics.median(f("between_median")), "max": max(f("between_median"))},
             "disp_withinA_max_over_positions": max(f("withinA_max")),
             "positions_with_between_median_above_5A": ranges(int(r["A_residue"].split(":")[1]) for r in disp if float(r["between_median"]) > 5.0),
             "positions_with_withinA_max_above_2A": ranges(int(r["A_residue"].split(":")[1]) for r in disp if float(r["withinA_max"]) > 2.0),
             "entry_B_within_cells_empty_not_zero": all(r["withinB_min"] == "" and r["withinB_max"] == "" and r["withinB_n"] == "0" for r in rows),
             "entry_B_positions_with_within_value": sum(r["withinB_max"] != "" for r in disp),
             "advisory_configured_channels": list(summ["parameters"]["advisory"]["min_difference"]),
             "advisory_flags_true": sum(r["advisory_flag"] == "True" for r in rows),
             "advisory_flags_assessed": sum(r["advisory_flag"] != "" for r in rows)}
    note = ("KIN-DESC-1's run_spec lists construct differences against UniProt P00519 (human ABL1); the deposited records of 1OPJ and 1OPK name "
            "P00520 (Mus musculus). KIN-DESC-1 is frozen and was not reopened; mouse-human differences in the kinase domain were not assessed here.")
    rec = {"id": "STRUCTPAIR-ABL-DEMO-2026-10-25", "kin_desc_1_note": note, "date": "2026-10-25", "structpair_version": man["software"]["structpair_version"], "schema": man["schema"],
           "kind": "descriptive demonstration; not an experiment; no difference is attributable to imatinib, nucleotide state, mutation or resistance",
           "spec": spec, "inputs": inputs, "input_checks": checks, "result_summary": stats, "resources": res,
           "outputs_sha256": {p.name: sha256_file(p) for p in sorted(dest.iterdir())}}
    (OUT / "abl_demo_record.json").write_text(json.dumps(rec, indent=1, default=str) + "\n")
    (dest / "chain_variation_caption.md").write_text(caption(rec))
    rec["outputs_sha256"] = {p.name: sha256_file(p) for p in sorted(dest.iterdir())}
    (OUT / "abl_demo_record.json").write_text(json.dumps(rec, indent=1, default=str) + "\n")
    print(json.dumps({"mapped": len(mapped), "stats": stats, "resources": res}, indent=1, default=str))
    return 0


def caption(rec: dict) -> str:
    A, B = rec["inputs"]["A"], rec["inputs"]["B"]
    ha, hb = A["header"], B["header"]

    def lig(h):
        seen = []
        for x in h["het"]:
            if x["name"] not in seen:
                seen.append(x["name"])
        return ", ".join(f"{c} (HETSYN {h['hetsyn'][c]})" if c in h["hetsyn"] else f"{c} ({h['hetnam'].get(c, '')[:40]})" for c in seen)

    def construct(h):
        return h["fragment"].replace("KINASE", "kinase").replace("DOMAIN", "domain")

    def missing(h):
        parts = []
        for ch, res in sorted(h["remark_465_missing"].items()):
            inside = [r for r in res if 248 <= r <= 534]
            if inside:
                parts.append(f"chain {ch}: {ranges(inside)}")
        return "; ".join(parts) or "none"
    s = rec["result_summary"]
    return f"""# Figure: 1OPJ chains A, B versus 1OPK chain A (descriptive; structpair {rec['structpair_version']})

![chain variation](chain_variation.png)

**Purpose.** To show the partial status of within-entry crystal-copy variation on a real pair: two chains in one entry, one in the other.

**What it shows.** Per mapped residue ({s['positions_mapped']} positions: {s['mapped_ranges']}), the Cα displacement between 1OPJ and 1OPK after a trimmed fit:
the min–max band and median over the 2 cross pairings, and the largest difference between the two 1OPJ chains (dashed). **1OPK has a single
chain, so its within-entry variation is unavailable and appears as missing, not zero** (empty cells; hatched legend entry). No advisory rule is
configured, so no flag is assessed. For orientation only: the median over positions of the between-entry displacement is {s['disp_between_median_over_positions']['median']:.2f} Å and its largest
value is {s['disp_between_median_over_positions']['max']:.1f} Å; positions above 5 Å are {s['positions_with_between_median_above_5A']} (PDB numbering); the largest
difference between the two 1OPJ chains is {s['disp_withinA_max_over_positions']:.1f} Å, above 2 Å at {s['positions_with_withinA_max_above_2A']}.

**Main limitation: nothing here can be attributed to imatinib, nucleotide state, mutation or resistance.** From the deposited records (cross-checked
with RCSB): both entries are *{ha['organism'].capitalize()}* Abl1 (UniProt P00520, not human), expressed in insect cells, with myristate bound in each.
They differ in the ligands (1OPJ: {lig(ha)}; 1OPK: {lig(hb)}), in the construct ({construct(ha)}, {ha['dbref'][0]['db_range'][0]}–{ha['dbref'][0]['db_range'][1]} in
1OPJ; {construct(hb)}, {hb['dbref'][0]['db_range'][0]}–{hb['dbref'][0]['db_range'][1]} in 1OPK, database numbering), in an engineered mutation (D363N, PDB residue 382, in 1OPK; flagged in
the output), in crystal form (1OPJ {ha['space_group']}, 2 chains; 1OPK {hb['space_group']}, 1 chain) and in resolution ({ha['resolution_A']} and {hb['resolution_A']} Å). Neither entry
contains a nucleotide. Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of them; the dashed
line is a descriptive range, not a noise floor, an uncertainty or an HDX replicate estimate.

**Missing residues (REMARK 465, inside the shared numbering 248–534).** 1OPJ {missing(ha)}; 1OPK {missing(hb)}. The correspondence keeps only residues
observed in all three chains, so it ends at {s['mapped_ranges'].split('-')[-1]}.

**Source data.** `chain_variation.csv`, `chains.csv`, `manifest.json` and `../abl_demo_record.json` (inputs and their SHA-256, header facts, checks).

**Reproduce.** `.venv/bin/python benchmarks/audit/structpair_abl_demo_2026-10-25/run_demo.py` ({rec['resources']['cpu_s']:.0f} CPU-s, {rec['resources']['wall_s']:.0f} s wall).
"""


if __name__ == "__main__":
    sys.exit(main())
