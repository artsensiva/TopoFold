"""KRAS K1 provenance checks (2026-10-25): metadata and header fields only; no coordinate is decoded or scored.

    python3 -I benchmarks/audit/kras_k1/check_k1_provenance.py --fetch-dir <new empty directory outside the repository>

1. Fetches S1 Table of Pantsar et al. 2018 (a 24 kB PDF, doi:10.1371/journal.pcbi.1006458.s016), checks its SHA-256, parses
   it with `pdftotext -layout` (poppler) and compares the replica counts per system with the file names of Zenodo record
   10.5281/zenodo.1346073 (read from ../next_app_feasibility/results/feasibility_metadata.json).
2. For the 15 KRAS G12D trajectories, reads by HTTP range request only header fields: the first frame header (atoms, time,
   box, compressed size), the second frame header (frame interval) and a 16 kB tail window in which the last frame header is
   located; the last frame is accepted only if its own size ends exactly at the end of the file. The compressed coordinate
   bytes are never decoded. Frame counts are implied from first and last time and the interval (uniform spacing assumed).
3. Reads the residue labels (never the values) of the authors' PyInteraph files (PyInteraph_data.zip, 123 kB, MD5 checked
   against the Zenodo listing) to record which residue name appears at position 12 and whether any nucleotide or Mg label
   appears, per system: system-level evidence of the mutation identity, not of the nucleotide and not per replica.
Writes results/k1_provenance_checks.json.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import struct
import subprocess
import sys
import urllib.request
import zipfile
from io import BytesIO
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / "results/k1_provenance_checks.json"
FEAS = ROOT / "benchmarks/audit/next_app_feasibility/results/feasibility_metadata.json"
S1_URL = "https://journals.plos.org/ploscompbiol/article/file?type=supplementary&id=10.1371/journal.pcbi.1006458.s016"
S1_SHA256 = "2a047d8c441a8c77189034e16c1b5c29efb1d0136a5ca3748bd4e5db8bdb0376"   # as read 2026-10-25 (also via Europe PMC)
ZENODO = "1346073"
TAIL = 16384
ROW = re.compile(r"^\s*(G12[A-Z]|Wild-type)\s+(GDP|GTP)\s+(\d+)\s+([\d ]+?)\s+(\d+)\s+(\d+)\s*$")


def sha256(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def parse_s1(text: str) -> dict:
    """{system: {'replicas', 'atoms', 'K+', 'Cl-'}} from the pdftotext layout of S1 Table; system = 'GDP_G12D', 'GTP_WT', ..."""
    out = {}
    for ln in text.splitlines():
        m = ROW.match(ln)
        if m:
            mut = "WT" if m.group(1) == "Wild-type" else m.group(1)
            out[f"{m.group(2)}_{mut}"] = {"replicas": int(m.group(3)), "atoms": int(m.group(4).replace(" ", "")),
                                          "K+": int(m.group(5)), "Cl-": int(m.group(6))}
    return out


def frame_fields(b: bytes) -> dict:
    """XTC compressed-frame header (XDR, big-endian): magic, natoms, step, time, 3x3 box, natoms, ..., compressed bytes."""
    magic, natoms, step = struct.unpack(">iii", b[:12])
    (t,) = struct.unpack(">f", b[12:16])
    box = struct.unpack(">9f", b[16:52])
    (nbytes,) = struct.unpack(">i", b[88:92])
    if magic != 1995 or natoms <= 9:
        raise ValueError("not a compressed XTC frame header")
    return {"natoms": natoms, "step": step, "time_ps": t, "box_nm": [box[0], box[4], box[8]],
            "off_diagonal_box_nm": [box[i] for i in (1, 2, 3, 5, 6, 7)], "frame_bytes": 92 + (nbytes + 3) // 4 * 4}


def last_frame(tail: bytes, base: int, size: int, natoms: int) -> dict | None:
    """Last frame header in a tail window whose first byte is file offset `base`; valid only if it ends at EOF."""
    pat = struct.pack(">ii", 1995, natoms)
    best = None
    p = tail.find(pat)
    while p != -1:
        if (base + p) % 4 == 0 and p + 92 <= len(tail):
            try:
                f = frame_fields(tail[p:p + 92])
            except ValueError:
                f = None
            if f and base + p + f["frame_bytes"] == size:
                best = {**f, "offset": base + p}
        p = tail.find(pat, p + 4)
    return best


def get_range(url: str, start: int, end: int) -> bytes:
    req = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}", "User-Agent": "TopoFold-k1-check/1"})
    with urllib.request.urlopen(req, timeout=60) as r:
        if r.status != 206:
            raise OSError(f"range request not honoured ({r.status})")
        return r.read()


def probe(key: str, size: int) -> dict:
    url = f"https://zenodo.org/records/{ZENODO}/files/{key}?download=1"
    f1 = frame_fields(get_range(url, 0, 91))
    f2 = frame_fields(get_range(url, f1["frame_bytes"], f1["frame_bytes"] + 91))
    base = size - TAIL
    fl = last_frame(get_range(url, base, size - 1), base, size, f1["natoms"])
    dt_ps = f2["time_ps"] - f1["time_ps"]
    rec = {"bytes": size, "natoms": f1["natoms"], "first_time_ps": f1["time_ps"], "interval_ps": dt_ps,
           "box_first_frame_nm": f1["box_nm"], "off_diagonal_box_nm": f1["off_diagonal_box_nm"],
           "bytes_read": 92 + 92 + TAIL}
    if fl is None:
        rec["last_frame"] = None
    else:
        rec["last_frame"] = {"time_ps": fl["time_ps"], "step": fl["step"], "natoms": fl["natoms"]}
        rec["implied_frames"] = round((fl["time_ps"] - f1["time_ps"]) / dt_ps) + 1 if dt_ps else None
    return rec


def pyinteraph_identity(zip_bytes: bytes) -> dict:
    """{system: {'residue_12_names': {name: n_labels}, 'nucleotide_or_mg_labels': {...}}} from the label column only."""
    out = {}
    with zipfile.ZipFile(BytesIO(zip_bytes)) as zf:
        for name in sorted(zf.namelist()):
            m = re.match(r"^PyInteraph_data/((?:GDP|GTP)_[A-Za-z0-9]+)/(salt-bridges|hydrophobic-clusters)\.dat$", name)
            if not m:
                continue
            rec = out.setdefault(m.group(1), {"residue_12_names": {}, "nucleotide_or_mg_labels": {}})
            for ln in zf.read(name).decode("utf8", "replace").splitlines():
                label = ln.split("\t")[0]            # the value column is discarded here and never parsed
                for r in re.finditer(r"[A-Z]-(\d+)([A-Z0-9]{3})", label):
                    if r.group(1) == "12":
                        rec["residue_12_names"][r.group(2)] = rec["residue_12_names"].get(r.group(2), 0) + 1
                    if r.group(2) in ("GDP", "GTP", "GNP", "MG2"):
                        rec["nucleotide_or_mg_labels"][r.group(2)] = rec["nucleotide_or_mg_labels"].get(r.group(2), 0) + 1
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fetch-dir", required=True, help="new empty directory outside the repository for the S1 Table PDF")
    args = ap.parse_args(argv)
    d = Path(args.fetch_dir)
    d.mkdir(parents=True, exist_ok=False)
    req = urllib.request.Request(S1_URL, headers={"User-Agent": "TopoFold-k1-check/1"})
    with urllib.request.urlopen(req, timeout=60) as r:
        pdf = r.read()
    if sha256(pdf) != S1_SHA256:
        raise SystemExit(f"STOP: S1 Table differs from the version read on 2026-10-25 ({sha256(pdf)})")
    (d / "s1_table.pdf").write_bytes(pdf)
    text = subprocess.run(["pdftotext", "-layout", str(d / "s1_table.pdf"), "-"], capture_output=True, text=True, check=True).stdout
    s1 = parse_s1(text)
    z = json.loads(FEAS.read_text())["kras_pantsar_zenodo"]
    files = {f["key"]: f["bytes"] for f in z["files"]}
    zen = {k: len(v) for k, v in z["replicas_per_system"].items()}
    compare = {k: {"s1_table": v["replicas"], "zenodo_files": zen.get(k), "equal": zen.get(k) == v["replicas"]} for k, v in s1.items()}
    g12d = sorted(k for k in files if re.match(r"^(GDP|GTP)_G12D_\d+\.xtc$", k))
    traj = {k: probe(k, files[k]) for k in g12d}
    zreq = urllib.request.Request(f"https://zenodo.org/records/{ZENODO}/files/PyInteraph_data.zip?download=1",
                                  headers={"User-Agent": "TopoFold-k1-check/1"})
    with urllib.request.urlopen(zreq, timeout=60) as r:
        pyi = r.read()
    expected_md5 = next(f["checksum"] for f in z["files"] if f["key"] == "PyInteraph_data.zip")
    if f"md5:{hashlib.md5(pyi).hexdigest()}" != expected_md5:
        raise SystemExit("STOP: PyInteraph_data.zip differs from the Zenodo listing")
    identity = pyinteraph_identity(pyi)
    rec = {"id": "KRAS-K1-PROVENANCE-CHECKS", "checked_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "scope": "metadata and XTC header fields only; no coordinate decoded or scored",
           "s1_table": {"source": "doi:10.1371/journal.pcbi.1006458.s016 (CC BY 4.0 article; supplementary PDF)", "url": S1_URL,
                        "sha256": S1_SHA256, "bytes": len(pdf), "systems": s1,
                        "total_replicas": sum(v["replicas"] for v in s1.values()),
                        "replica_start_note": "footnote: an equilibrated system was simulated at 350 K with restraints (K-Ras, "
                                              "crystal waters, Mg, ligand); every 50 ps one frame was picked as the initial "
                                              "structure of each replica"},
           "replica_counts": {"compared": len(compare), "equal": sum(c["equal"] for c in compare.values()), "detail": compare,
                              "zenodo_trajectory_files": sum(zen.values())},
           "pyinteraph_label_identity": {"source": "Zenodo PyInteraph_data.zip (authors' analysis of the same replicas; labels only)",
                                         "md5": expected_md5, "systems": identity},
           "g12d_trajectories": traj,
           "g12d_summary": {
               "files": len(traj),
               "natoms": sorted({t["natoms"] for t in traj.values()}),
               "first_time_ps": sorted({t["first_time_ps"] for t in traj.values()}),
               "interval_ps": sorted({t["interval_ps"] for t in traj.values()}),
               "last_time_ps": {k: (t["last_frame"] or {}).get("time_ps") for k, t in traj.items()},
               "implied_frames": {k: t.get("implied_frames") for k, t in traj.items()},
               "box_first_frame_nm": {k: t["box_first_frame_nm"] for k, t in traj.items()},
               "bytes_read_total": sum(t["bytes_read"] for t in traj.values())}}
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    s = rec["g12d_summary"]
    print(json.dumps({"s1_systems": len(s1), "replica_counts_equal": f"{rec['replica_counts']['equal']}/{rec['replica_counts']['compared']}",
                      "total_replicas_s1": rec["s1_table"]["total_replicas"], "zenodo_files": rec["replica_counts"]["zenodo_trajectory_files"],
                      "g12d": {k: s[k] for k in ("files", "natoms", "first_time_ps", "interval_ps", "bytes_read_total")},
                      "last_frames_found": sum(v is not None for v in s["last_time_ps"].values()),
                      "last_time_ps_set": sorted({v for v in s["last_time_ps"].values() if v is not None}),
                      "implied_frames_set": sorted({v for v in s["implied_frames"].values() if v is not None}),
                      "residue_12_names_g12d": {k: v["residue_12_names"] for k, v in identity.items() if "G12D" in k}}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
