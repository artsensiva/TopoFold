"""Thin bridge from one mdCATH HDF5 replica to the Cα adapter (needs h5py and MDAnalysis).

No new file parser: the topology is the file's own ``pdbProteinAtoms`` text read by MDAnalysis's PDB parser, and
coordinates are the stored ``coords`` array. The bridge only checks that the two agree and records provenance.

* **Mapping.** The PDB topology must have exactly as many atoms as ``coords`` and, atom by atom, the same residue
  number, residue name and chain as the HDF5 ``resid``, ``resname`` and ``chain`` arrays; otherwise it raises.
* **Units.** ``coords`` must carry ``unit == "Angstrom"``; the adapter is told so explicitly
  (``assume_angstrom_if_undeclared``) and the basis is recorded. The stored ``box`` (3 × 3, nm, one per replica,
  first frame) is converted to Å (× 10) and applied to every frame; this is recorded as an approximation.
* **Frames and time.** Source frame indices are preserved. The files store no timestamps; the inferred time
  (index × ``frame_interval_ns``, 1 ns per the dataset paper) is reported separately and labelled as inferred.
* **Frame identity** cannot be verified (coordinates only); the adapter reports ``not_verifiable``.
"""

from __future__ import annotations

import hashlib
import io
import json
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np

from .ca_adapter import AdapterError, extract_ca

BRIDGE_VERSION = "0.1.0"


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _text(x) -> str:
    return x.decode() if isinstance(x, bytes) else str(x)


def mdcath_universe(h5_path, domain: str, temperature: str, replica: str, *, start=None, stop=None, step=None):
    """Return (MDAnalysis Universe with in-memory coordinates for the selected frames, bridge record)."""
    import h5py
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    h5_path = Path(h5_path)
    with h5py.File(h5_path, "r") as f:
        layout = _text(f.attrs.get("layout", ""))
        g = f[domain]
        rep = g[str(temperature)][str(replica)]
        ds = rep["coords"]
        unit = _text(ds.attrs.get("unit", ""))
        if unit != "Angstrom":
            raise AdapterError(f"coords unit {unit!r}; the bridge accepts only 'Angstrom'")
        n_frames_src = int(ds.shape[0])
        if "numFrames" in rep.attrs and int(rep.attrs["numFrames"]) != n_frames_src:
            raise AdapterError("numFrames attribute differs from the coords length")
        frames = np.arange(n_frames_src)[slice(start, stop, step)]
        if frames.size == 0:
            raise AdapterError("frame slice selected no frames")
        coords = np.asarray(ds[frames[0]:frames[-1] + 1][frames - frames[0]], dtype=np.float32)
        box_nm = np.asarray(rep["box"][()], dtype=np.float64) if "box" in rep else None
        pdb_text = _text(g["pdbProteinAtoms"][()])
        h5_resid = np.asarray(g["resid"][()])
        h5_resname = [_text(x) for x in g["resname"][()]]
        h5_chain = [_text(x) for x in g["chain"][()]]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        u = mda.Universe(io.StringIO(pdb_text), format="PDB")
    if u.atoms.n_atoms != coords.shape[1]:
        raise AdapterError(f"topology has {u.atoms.n_atoms} atoms but coords has {coords.shape[1]}")
    mism = [i for i in range(u.atoms.n_atoms)
            if int(u.atoms.resids[i]) != int(h5_resid[i]) or str(u.atoms.resnames[i]) != h5_resname[i]
            or str(u.atoms.chainIDs[i]) != h5_chain[i]]
    if mism:
        i = mism[0]
        raise AdapterError(f"topology/HDF5 mapping differs at atom {i}: "
                           f"({u.atoms.resids[i]}, {u.atoms.resnames[i]}, {u.atoms.chainIDs[i]}) vs "
                           f"({h5_resid[i]}, {h5_resname[i]}, {h5_chain[i]}); {len(mism)} atoms differ")
    dims = None
    box_note = "no box stored"
    if box_nm is not None:
        if box_nm.shape != (3, 3) or np.any(np.abs(box_nm - np.diag(np.diag(box_nm))) > 0):
            raise AdapterError(f"box is not a diagonal 3x3 matrix: {box_nm.tolist()}")
        lengths = np.diag(box_nm) * 10.0
        dims = np.tile(np.concatenate([lengths, [90.0, 90.0, 90.0]]), (len(frames), 1)).astype(np.float32)
        box_note = "stored box (nm, one per replica) x 10 -> Angstrom, applied to every frame (approximation)"
    u.load_new(coords, format=MemoryReader, dimensions=dims)
    record = {
        "bridge": {"name": "topofold_inputs.mdcath_bridge", "version": BRIDGE_VERSION,
                   "source_sha256": _sha256(__file__)},
        "source_file": {"role": "hdf5", "path": str(h5_path.resolve()), "bytes": h5_path.stat().st_size,
                        "sha256": _sha256(h5_path)},
        "dataset": {"layout": layout, "domain": domain, "temperature": str(temperature), "replica": str(replica),
                    "coords": f"/{domain}/{temperature}/{replica}/coords", "n_frames_source": n_frames_src},
        "mapping_check": {"atoms": int(u.atoms.n_atoms), "fields": ["resid", "resname", "chain"], "mismatches": 0},
        "units": {"coords_attr": unit, "box": box_note,
                  "box_angstrom": None if box_nm is None else (np.diag(box_nm) * 10.0).tolist()},
        "frames": {"source_index": frames.tolist(), "time_source": "none in file",
                   "inferred_time_ns": frames.astype(float).tolist(), "inferred_time_basis":
                   "INFERRED: index x 1 ns (mdCATH paper); the file stores no timestamps"},
        "topology_parse_warnings": list(dict.fromkeys(str(w.message) for w in caught)),
    }
    return u, record


def extract_ca_mdcath(h5_path, domain, temperature, replica, run_id, *, start=None, stop=None, step=None, **kwargs):
    """Bridge one replica and run the Cα adapter on it; source frame indices replace the adapter's local ones."""
    u, rec = mdcath_universe(h5_path, domain, temperature, replica, start=start, stop=stop, step=step)
    res = extract_ca(u, run_id, assume_angstrom_if_undeclared=True, **kwargs)
    res.frame_index = np.asarray(rec["frames"]["source_index"], dtype=np.int64)
    res.report["provenance"]["source_files"] = [rec["source_file"]]
    res.report["provenance"]["bridge"] = {k: v for k, v in rec.items() if k != "frames"}
    res.report["frames_inferred_time"] = {"inferred_time_ns": rec["frames"]["inferred_time_ns"],
                                          "basis": rec["frames"]["inferred_time_basis"]}
    res.report["units"]["basis"] = (f"HDF5 coords unit attribute {rec['units']['coords_attr']!r}; "
                                    "declared to the adapter as Angstrom")
    return res


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(prog="python -m topofold_inputs.mdcath_bridge")
    ap.add_argument("h5")
    ap.add_argument("domain")
    ap.add_argument("temperature")
    ap.add_argument("replica")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--start", type=int)
    ap.add_argument("--stop", type=int)
    ap.add_argument("--step", type=int)
    ap.add_argument("--pbc", default="none")
    ap.add_argument("--altloc", default="reject")
    ap.add_argument("--max-ca-ca", type=float, default=4.2)
    a = ap.parse_args(argv)
    try:
        res = extract_ca_mdcath(a.h5, a.domain, a.temperature, a.replica, a.run_id, start=a.start, stop=a.stop,
                                step=a.step, pbc=a.pbc, altloc=a.altloc, max_ca_ca_angstrom=a.max_ca_ca)
    except AdapterError as exc:
        print(f"rejected: {exc}", file=sys.stderr)
        return 2
    npz, js = res.export(a.out)
    s = res.report["summary"]
    print(json.dumps({"run_id": a.run_id, "frames": s["frames"], "residues": s["residues"], "breaks": s["breaks"],
                      "npz": str(npz), "json": str(js)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
