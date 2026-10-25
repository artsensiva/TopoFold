"""Next-application feasibility (2026-10-25): metadata collection only. Nothing is scored and no trajectory is downloaded.

    python3 -I benchmarks/audit/next_app_feasibility/collect_metadata.py [--fetch-dir <new empty directory>]

Collects, from public APIs:
  * the Zenodo record of Pantsar et al. 2018 (KRAS G12X MD; file names, sizes, MD5 as published);
  * two XTC frame headers (184 bytes, HTTP range requests) of 4 trajectories, for atom count and frame interval;
  * one small file of that record, the 38 kB backbone topology `backbone_4-164.pdb`, only to verify the residue mapping
    by sequence against UniProt P01116 isoforms 2A and 2B (into --fetch-dir; read as untrusted text, never executed);
  * RCSB entry metadata for the candidate structural endpoints (Mpro multi-temperature 7MHF-7MHQ, Mpro NMR 9PYS/9PYT,
    KRAS 4OBE start, KRAS G12D 5US4/5USJ) and for every PDB ID in the Moonshot RUN metadata already committed;
  * the UniProt P0DTD1 nsp5 chain range.
Writes results/feasibility_metadata.json.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import re
import struct
import sys
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / "results/feasibility_metadata.json"
MOONSHOT_CSV = ROOT / "benchmarks/audit/mpro_dataset_search/fah-metadata.csv"
ZENODO_RECORD = "1346073"          # version DOI 10.5281/zenodo.1346073 (concept 10.5281/zenodo.1346072)
TOPOLOGY = "backbone_4-164.pdb"
ENDPOINT_ENTRIES = {
    "7MHF": "Mpro 100 K multiconformer", "7MHG": "Mpro 240 K multiconformer", "7MHH": "Mpro 277 K multiconformer",
    "7MHI": "Mpro 298 K multiconformer", "7MHJ": "Mpro 298 K high humidity multiconformer",
    "7MHK": "Mpro 310 K multiconformer", "7MHL": "Mpro 100 K ensemble", "7MHM": "Mpro 240 K ensemble",
    "7MHN": "Mpro 277 K ensemble", "7MHO": "Mpro 298 K ensemble", "7MHP": "Mpro 298 K high humidity ensemble",
    "7MHQ": "Mpro 310 K ensemble", "9PYS": "monomeric Mpro NMR RDC, helical domain",
    "9PYT": "monomeric Mpro NMR RDC, catalytic domain", "4OBE": "KRAS GDP, start of the Pantsar et al. MD",
    "5US4": "KRAS G12D GDP (KRAS-STRUCT-2)", "5USJ": "KRAS G12D GMPPNP (KRAS-STRUCT-2)",
}
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
         "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V", "HID": "H", "HIE": "H", "HIP": "H", "HSD": "H", "HSE": "H", "HSP": "H"}


def get(url: str, accept: str = "application/json", tries: int = 3) -> bytes:
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": accept, "User-Agent": "TopoFold-feasibility/1"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except OSError:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))
    raise RuntimeError("unreachable")


def get_json(url: str):
    return json.loads(get(url))


def fasta_seq(text: str) -> str:
    return "".join(ln.strip() for ln in text.splitlines() if ln and not ln.startswith(">"))


def pdb_sequence(text: str) -> list:
    """(resid, one-letter) per residue from ATOM records, in file order; unknown names are 'X'."""
    out, seen = [], set()
    for ln in text.splitlines():
        if ln.startswith(("ATOM", "HETATM")) and len(ln) >= 26:
            key = (ln[21], ln[22:27])
            if key not in seen:
                seen.add(key)
                out.append((int(ln[22:26]), THREE.get(ln[17:20].strip(), "X")))
    return out


def map_to_reference(residues: list, ref: str) -> dict:
    """Exact match of author-numbered residues to a reference sequence at offset 0 (resid = reference index)."""
    mism = [(r, aa, ref[r - 1] if 0 < r <= len(ref) else None) for r, aa in residues if not (0 < r <= len(ref) and ref[r - 1] == aa)]
    return {"n": len(residues), "first": residues[0][0] if residues else None, "last": residues[-1][0] if residues else None,
            "mismatches": [{"resid": r, "topology": a, "reference": b} for r, a, b in mism]}


def rcsb_entry(pid: str) -> dict:
    d = get_json(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
    info = d.get("rcsb_entry_info", {})
    return {"title": d.get("struct", {}).get("title"), "method": (d.get("exptl") or [{}])[0].get("method"),
            "resolution_A": (info.get("resolution_combined") or [None])[0], "models": info.get("deposited_model_count"),
            "nonpolymer_entities": info.get("nonpolymer_entity_count"),
            "released": d.get("rcsb_accession_info", {}).get("initial_release_date", "")[:10]}


def get_range(url: str, start: int, end: int) -> bytes:
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}", "User-Agent": "TopoFold-feasibility/1"})
    with urllib.request.urlopen(req, timeout=60) as r:
        if r.status != 206:
            raise OSError(f"range request not honoured ({r.status})")
        return r.read()


def xtc_frame_header(b: bytes) -> dict:
    """First 92 bytes of an XTC frame (XDR, big-endian): magic, natoms, step, time, box, natoms, precision, ..., nbytes."""
    magic, natoms, step = struct.unpack(">iii", b[:12])
    (time_ps,) = struct.unpack(">f", b[12:16])
    if magic != 1995 or natoms <= 9:
        raise ValueError(f"not a compressed XTC frame header (magic {magic}, natoms {natoms})")
    (nbytes,) = struct.unpack(">i", b[88:92])
    return {"natoms": natoms, "step": step, "time_ps": time_ps, "frame_bytes": 92 + (nbytes + 3) // 4 * 4}


def xtc_probe(record: str, key: str, size: int) -> dict:
    """Two tiny range reads: the first two frame headers give natoms and the frame interval; the frame count is only
    estimated from the file size (frames are compressed to varying sizes)."""
    url = f"https://zenodo.org/records/{record}/files/{key}?download=1"
    h1 = xtc_frame_header(get_range(url, 0, 91))
    h2 = xtc_frame_header(get_range(url, h1["frame_bytes"], h1["frame_bytes"] + 91))
    dt_ps = h2["time_ps"] - h1["time_ps"]
    return {"natoms": h1["natoms"], "first_time_ps": h1["time_ps"], "frame_interval_ps": dt_ps,
            "bytes_read": 184, "estimated_frames_from_size": round(size / h1["frame_bytes"]),
            "estimate_note": "file size / first-frame size; not a frame count"}


def zenodo(record: str) -> dict:
    d = get_json(f"https://zenodo.org/api/records/{record}")
    m = d["metadata"]
    files = sorted(({"key": f["key"], "bytes": f["size"], "checksum": f["checksum"]} for f in d.get("files", [])),
                   key=lambda f: f["key"])
    systems = {}
    for f in files:
        mm = re.match(r"^(GDP|GTP)_(WT|G12[A-Z])_(\d+)\.xtc$", f["key"])
        if mm:
            systems.setdefault(f"{mm.group(1)}_{mm.group(2)}", []).append(int(mm.group(3)))
    return {"doi": d.get("doi"), "concept_doi": d.get("conceptdoi"), "title": m.get("title"),
            "publication_date": m.get("publication_date"), "license": (m.get("license") or {}).get("id"),
            "access_right": m.get("access_right"), "files": files, "n_files": len(files),
            "total_bytes": sum(f["bytes"] for f in files),
            "replicas_per_system": {k: sorted(v) for k, v in sorted(systems.items())}}


def moonshot_subsets(rows: list) -> dict:
    """Download size of candidate RUN subsets in the apo (13430) and holo (13432) monomer His41(0) projects, from S3
    listings (no file is read). 'Full-length' clones are those at least 90 % of the project's median file size: a size
    heuristic for unequal clone lengths, not a frame count."""
    sys.path.insert(0, str(ROOT / "benchmarks/audit/mpro_dataset_search"))
    import list_moonshot_metadata as L
    lst = {p: dict(L.list_keys(f"{L.PREFIX}solute-hdf5-trajectories/{p}/")) for p in ("13430", "13432")}
    med = {p: sorted(v.values())[len(v) // 2] for p, v in lst.items()}
    out = {}
    for site in ("XChem Screen - active site", "XChem Screen - xtal contact"):
        runs = [r for r in rows if r["site_name"] == site and r["pdb_entry"]]
        per = {}
        for r in runs:
            n = r["RUN"].removeprefix("RUN")
            files = {p: {k: s for k, s in lst[p].items() if k.rsplit("/", 1)[-1].startswith(f"run{n}-clone")} for p in lst}
            per[r["RUN"]] = {"pdb_entry": r["pdb_entry"], **{p: {"clones": len(f), "full_length": sum(s >= 0.9 * med[p] for s in f.values()),
                                                                 "bytes": sum(f.values()),
                                                                 "bytes_3_largest": sum(sorted(f.values())[-3:])} for p, f in files.items()}}
        out[site] = {"runs": len(per), "bytes_total": sum(v[p]["bytes"] for v in per.values() for p in lst),
                     "runs_with_3_full_length_clones_each": sum(all(v[p]["full_length"] >= 3 for p in lst) for v in per.values()),
                     "bytes_3_full_length_clones_per_arm": sum(v[p]["bytes_3_largest"] for v in per.values() for p in lst
                                                               if all(v[q]["full_length"] >= 3 for q in lst)),
                     "per_run": per}
    return {"median_file_bytes": med, "subsets": out}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fetch-dir", default=None, help="new empty directory for the one small topology file")
    args = ap.parse_args(argv)
    rec: dict = {"id": "NEXT-APP-FEASIBILITY-1", "collected_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "scope": "metadata only; one 38 kB topology file for the residue mapping and 184 bytes (two frame headers) of "
                          "each of 4 trajectories; no trajectory downloaded, nothing scored"}
    z = zenodo(ZENODO_RECORD)
    rec["kras_pantsar_zenodo"] = z
    p2a = fasta_seq(get("https://rest.uniprot.org/uniprotkb/P01116-1.fasta", "text/plain").decode())
    p2b = fasta_seq(get("https://rest.uniprot.org/uniprotkb/P01116-2.fasta", "text/plain").decode())
    rec["uniprot_kras"] = {"P01116-1 (isoform 2A)": {"length": len(p2a), "sha256": hashlib.sha256(p2a.encode()).hexdigest()},
                           "P01116-2 (isoform 2B)": {"length": len(p2b), "sha256": hashlib.sha256(p2b.encode()).hexdigest()},
                           "first_difference_1based": next((i + 1 for i, (a, b) in enumerate(zip(p2a, p2b)) if a != b), None)}
    if args.fetch_dir:
        d = Path(args.fetch_dir)
        d.mkdir(parents=True, exist_ok=False)
        entry = next(f for f in z["files"] if f["key"] == TOPOLOGY)
        raw = get(f"https://zenodo.org/records/{ZENODO_RECORD}/files/{TOPOLOGY}?download=1", "*/*")
        (d / TOPOLOGY).write_bytes(raw)
        md5 = hashlib.md5(raw).hexdigest()
        res = pdb_sequence(raw.decode("ascii", "replace"))
        atoms = [ln for ln in raw.decode("ascii", "replace").splitlines() if ln.startswith(("ATOM", "HETATM"))]
        rec["kras_topology_check"] = {
            "atoms": len(atoms), "atom_names": sorted({ln[12:16].strip() for ln in atoms}),
            "file": TOPOLOGY, "bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest(), "md5": md5,
            "md5_matches_zenodo": f"md5:{md5}" == entry["checksum"],
            "residues": len(res), "sequence": "".join(a for _, a in res),
            "map_P01116-1": map_to_reference(res, p2a), "map_P01116-2": map_to_reference(res, p2b),
            "local_copy": "outside the repository (--fetch-dir); not committed"}
    sizes = {f["key"]: f["bytes"] for f in z["files"]}
    rec["kras_xtc_probe"] = {k: xtc_probe(ZENODO_RECORD, k, sizes[k])
                             for k in ("GDP_WT_1.xtc", "GTP_WT_1.xtc", "GDP_G12D_1.xtc", "GTP_G12D_1.xtc")}
    nsp5 = [f for f in get_json("https://rest.uniprot.org/uniprotkb/P0DTD1.json?fields=ft_chain")["features"]
            if "nsp5" in f.get("description", "")]
    rec["uniprot_mpro"] = [{"description": f["description"], "start": f["location"]["start"]["value"],
                            "end": f["location"]["end"]["value"]} for f in nsp5]
    rec["endpoint_entries"] = {pid: {"role": role, **rcsb_entry(pid)} for pid, role in ENDPOINT_ENTRIES.items()}
    rows = list(csv.DictReader(open(MOONSHOT_CSV)))
    ids = sorted({r["pdb_entry"] for r in rows if r["pdb_entry"]})
    ms = {pid: rcsb_entry(pid) for pid in ids}
    sites = {}
    for r in rows:
        s = sites.setdefault(r["site_name"], {"runs": 0, "with_pdb_entry": 0})
        s["runs"] += 1
        s["with_pdb_entry"] += bool(r["pdb_entry"])
    rec["moonshot_pdb_entries"] = {"runs": len(rows), "runs_with_pdb_entry": sum(bool(r["pdb_entry"]) for r in rows),
                                   "unique_pdb_entries": len(ids),
                                   "xray_entries": sum(v["method"] == "X-RAY DIFFRACTION" for v in ms.values()),
                                   "entries_with_nonpolymer": sum((v["nonpolymer_entities"] or 0) > 0 for v in ms.values()),
                                   "by_site": sites, "entries": ms}
    rec["moonshot_subsets"] = moonshot_subsets(rows)
    rec["kras_subsets"] = {name: {"files": len(keys), "bytes": sum(sizes[k] for k in keys)} for name, keys in {
        "G12D GDP + GTP + topology files": [k for k in sizes if k.startswith(("GDP_G12D_", "GTP_G12D_")) or k.startswith("backbone_")],
        "WT GDP + GTP + topology files": [k for k in sizes if k.startswith(("GDP_WT_", "GTP_WT_")) or k.startswith("backbone_")],
    }.items()}
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    t = rec.get("kras_topology_check", {})
    print(json.dumps({"zenodo_files": z["n_files"], "zenodo_bytes": z["total_bytes"], "license": z["license"],
                      "replicas": {k: len(v) for k, v in z["replicas_per_system"].items()},
                      "topology_md5_ok": t.get("md5_matches_zenodo"), "residues": t.get("residues"),
                      "mismatch_2A": len(t.get("map_P01116-1", {}).get("mismatches", [])) if t else None,
                      "xtc": {k: (v["natoms"], v["frame_interval_ps"], v["estimated_frames_from_size"]) for k, v in rec["kras_xtc_probe"].items()},
                      "mismatch_2B": len(t.get("map_P01116-2", {}).get("mismatches", [])) if t else None,
                      "moonshot": {k: rec["moonshot_pdb_entries"][k] for k in
                                   ("runs_with_pdb_entry", "unique_pdb_entries", "xray_entries", "entries_with_nonpolymer")}},
                     indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
