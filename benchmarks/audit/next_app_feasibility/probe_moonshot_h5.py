"""Read only the HDF5 metadata, the time axis and the topology of remote Moonshot Mpro trajectory files (no coordinates).

    <python with h5py> -I benchmarks/audit/next_app_feasibility/probe_moonshot_h5.py

The protein sequence is mapped to UniProt P0DTD1 nsp5 by exact match. h5py opens a remote file through a file-like object that issues HTTP range requests in 64 kB blocks, so only the blocks
that hold the metadata, `time` and `topology` are transferred. A hard cap of 8 MB per file aborts the probe. The
`coordinates` dataset is never read: only its shape, dtype and attributes. Writes results/moonshot_h5_probe.json.
"""

from __future__ import annotations

import datetime as dt
import json
import sys
import urllib.request
from collections import Counter
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "results/moonshot_h5_probe.json"
BUCKET = "https://fah-public-data-covid19-moonshot-dynamics.s3.us-east-2.amazonaws.com"
PREFIX = "SARS-CoV-2-Mpro-Diamond-structures/solute-hdf5-trajectories"
# RUN9 = Mpro-x0354_0A, PDB 5RE9, site 'XChem Screen - active site' (fah-metadata.csv); apo and holo monomer, His41(0)
FILES = {"apo_monomer_13430": f"{PREFIX}/13430/run9-clone0.h5", "holo_monomer_13432": f"{PREFIX}/13432/run9-clone0.h5"}
BLOCK, CAP = 1 << 16, 8 << 20
NSP5 = (3264, 3569)


class RangeFile:
    """Read-only, seekable file object over HTTP range requests, with a block cache and a byte cap."""

    def __init__(self, url: str, cap: int = CAP):
        self.url, self.cap, self.pos, self.fetched, self.cache = url, cap, 0, 0, {}
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            self.size = int(r.headers["Content-Length"])

    def _block(self, i: int) -> bytes:
        if i not in self.cache:
            start, end = i * BLOCK, min((i + 1) * BLOCK, self.size) - 1
            if self.fetched + (end - start + 1) > self.cap:
                raise RuntimeError(f"byte cap {self.cap} reached")
            req = urllib.request.Request(self.url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(req, timeout=60) as r:
                if r.status != 206:
                    raise OSError(f"range request not honoured ({r.status})")
                self.cache[i] = r.read()
            self.fetched += len(self.cache[i])
        return self.cache[i]

    def read(self, n: int = -1) -> bytes:
        if n < 0:
            n = self.size - self.pos
        n = max(0, min(n, self.size - self.pos))
        out, p = bytearray(), self.pos
        while len(out) < n:
            b = self._block(p // BLOCK)
            off = p % BLOCK
            take = b[off:off + n - len(out)]
            out += take
            p += len(take)
        self.pos = p
        return bytes(out)

    def readinto(self, buf) -> int:
        data = self.read(len(buf))
        buf[:len(data)] = data
        return len(data)

    def seek(self, off: int, whence: int = 0) -> int:
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def tell(self) -> int:
        return self.pos

    def seekable(self) -> bool:
        return True

    def readable(self) -> bool:
        return True


def attrs(obj) -> dict:
    out = {}
    for k, v in obj.attrs.items():
        v = v.decode() if isinstance(v, bytes) else v.tolist() if isinstance(v, np.ndarray) else v
        out[k] = v if isinstance(v, (str, int, float, list, bool)) or v is None else str(v)
    return out


def probe(url: str) -> dict:
    f = RangeFile(url)
    with h5py.File(f, "r") as h:
        datasets = {}
        h.visititems(lambda name, o: datasets.__setitem__(name, {"shape": list(o.shape), "dtype": str(o.dtype),
                                                                 "attrs": attrs(o)}) if isinstance(o, h5py.Dataset) else None)
        root = attrs(h)
        t = np.asarray(h["time"][...], dtype=np.float64) if "time" in h else None
        topo = json.loads(h["topology"][0].decode() if h["topology"].shape else h["topology"][()].decode()) if "topology" in h else None
    rec = {"bytes_total": f.size, "bytes_fetched": f.fetched, "root_attrs": root, "datasets": datasets}
    if t is not None:
        d = np.diff(t)
        rec["time"] = {"n_frames": int(t.size), "first": float(t[0]), "last": float(t[-1]),
                       "interval_min": float(d.min()) if d.size else None, "interval_max": float(d.max()) if d.size else None,
                       "units": datasets["time"]["attrs"].get("units")}
    if topo is not None:
        residues = [(c["index"], r["resSeq"], r["name"], len(r["atoms"])) for c in topo["chains"] for r in c["residues"]]
        protein = [r for r in residues if r[2] not in ("HOH", "NA", "CL", "LIG", "UNK", "MOL") and r[3] > 1]
        rec["topology"] = {"chains": len(topo["chains"]), "residues": len(residues), "atoms": sum(r[3] for r in residues),
                           "residue_names_not_standard": sorted({r[2] for r in residues} - STANDARD),
                           "per_chain_resSeq_range": {str(ci): [min(r[1] for r in residues if r[0] == ci),
                                                                max(r[1] for r in residues if r[0] == ci)]
                                                      for ci in sorted({r[0] for r in residues})},
                           "sequence_one_letter": "".join(ONE.get(r[2], "X") for r in protein),
                           "elements": dict(Counter(a.get("element") for c in topo["chains"] for r in c["residues"]
                                                    for a in r["atoms"]))}
    return rec


STANDARD = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
            "THR", "TRP", "TYR", "VAL"}
ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
       "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "CYM": "C", "ASH": "D", "GLH": "E", "LYN": "K"}


def main() -> int:
    out = {"id": "NEXT-APP-FEASIBILITY-1-moonshot-h5", "probed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "scope": "HDF5 metadata, time axis and topology only; coordinates not read", "files": {}}
    with urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/P0DTD1.fasta", timeout=60) as r:
        pp = "".join(ln.strip() for ln in r.read().decode().splitlines() if not ln.startswith(">"))
    nsp5 = pp[NSP5[0] - 1:NSP5[1]]      # UniProt P0DTD1 chain '3C-like proteinase nsp5' 3264-3569
    for name, key in FILES.items():
        rec = {"key": key, **probe(f"{BUCKET}/{key}")}
        s = rec.get("topology", {}).get("sequence_one_letter", "")
        rec["mapping_nsp5"] = {"reference": f"UniProt P0DTD1 {NSP5[0]}-{NSP5[1]} (nsp5 1-{len(nsp5)})", "offset": 0,
                               "compared": len(s), "mismatches": [[i + 1, a, nsp5[i] if i < len(nsp5) else None]
                                                                 for i, a in enumerate(s) if i >= len(nsp5) or nsp5[i] != a],
                               "reference_residues_absent": list(range(len(s) + 1, len(nsp5) + 1))}
        out["files"][name] = rec
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    for name, r in out["files"].items():
        print(name, r["bytes_fetched"], "of", r["bytes_total"], json.dumps(r.get("time")),
              json.dumps({k: v for k, v in r.get("topology", {}).items() if k != "sequence_one_letter"}))
        print(" datasets:", {k: (v["shape"], v["dtype"], v["attrs"]) for k, v in r["datasets"].items()})
    return 0


if __name__ == "__main__":
    sys.exit(main())
