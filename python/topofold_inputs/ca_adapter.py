"""Cα extraction from MDAnalysis Universes with explicit residue identity, gaps and provenance.

Engineering adapter, not a validated analysis. It never parses file formats itself: reading, unit conversion and
format detection are MDAnalysis's. It adds explicit policies on top:

* **Units.** Output coordinates are in Å (MDAnalysis converts every supported reader to Å). The source length unit is
  recorded. Output is float64, but values are only as precise as the source (MDAnalysis positions are float32, and the
  file format may quantize further); conversion to float64 does not restore precision.
* **Identity.** Each output residue is one topology residue, identified by (chain, resid, icode) plus its MDAnalysis
  residue and atom indices. Identity is fixed once and reused for every frame. Correspondence to any external numbering
  is only what the caller supplies; nothing is inferred from residue numbers across proteins.
* **Gaps.** Chain/segment boundaries, residues without Cα, numbering discontinuities and over-long Cα–Cα distances
  are reported as breaks. ``segments`` lists the maximal runs without any break; consumers must not join across them.
* **Ambiguity.** Duplicate (chain, resid, icode) keys, several Cα atoms in a residue that the altloc policy cannot
  resolve, and non-finite coordinates raise :class:`AdapterError` instead of producing coordinates.
* **Periodic boundaries.** ``pbc="none"`` (default) leaves coordinates as stored; ``pbc="make_whole_sequential"``
  applies the WP-B sequential minimum-image rule within each segment, for orthorhombic boxes only.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import os
import platform
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

ADAPTER_VERSION = "0.1.0"
LENGTH_UNIT = "angstrom"
TIME_UNIT = "ps"

ALTLOC_POLICIES = ("reject", "blank_then_A", "highest_occupancy")
PBC_POLICIES = ("none", "make_whole_sequential")

# Nominal resolution by format convention. Not read from the file; recorded so consumers know the claim's basis.
_NOMINAL_RESOLUTION = {
    "PDB": "0.001 Å (fixed-column %8.3f coordinates)",
    "ENT": "0.001 Å (fixed-column %8.3f coordinates)",
    "GRO": "0.01 Å if written with the usual %8.3f nm fields; other precisions exist",
    "XTC": "lossy compression; precision set by the writer (commonly 0.001 nm = 0.01 Å); not read by this adapter",
    "DCD": "float32",
    "TRR": "float32 or float64 as written",
}


class AdapterError(ValueError):
    """Raised for inputs the adapter refuses to convert (ambiguity, non-finite values, unsupported policy)."""


@dataclass(frozen=True)
class ResidueSpec:
    """One requested residue: (chain, resid, icode) in the source structure and the caller's label for it."""

    chain: str
    resid: int
    icode: str = ""
    label: str | None = None

    @classmethod
    def from_mapping(cls, m: Mapping[str, Any]) -> "ResidueSpec":
        return cls(chain=str(m["chain"]), resid=int(m["resid"]), icode=str(m.get("icode", "") or ""),
                   label=None if m.get("label") is None else str(m["label"]))

    @property
    def key(self) -> tuple[str, int, str]:
        return (self.chain, self.resid, self.icode)


@dataclass
class CATrajectory:
    """Cα coordinates of one run with per-residue identity, breaks and provenance."""

    run_id: str
    coordinates: np.ndarray  # (F, N, 3) float64, Å; NaN where ``present`` is False
    present: np.ndarray  # (N,) bool
    frame_index: np.ndarray  # (F,) int64, index in the source trajectory
    time_ps: np.ndarray  # (F,) float64; NaN when the format provides no time
    residues: list[dict[str, Any]]
    breaks: list[dict[str, Any]]
    segments: list[list[int]]  # half-open [start, stop) output-index ranges
    report: dict[str, Any] = field(default_factory=dict)

    def to_report(self) -> dict[str, Any]:
        out = dict(self.report)
        out["run_id"] = self.run_id
        out["shape"] = {"frames": int(self.coordinates.shape[0]), "residues": int(self.coordinates.shape[1])}
        out["residues"] = self.residues
        out["breaks"] = self.breaks
        out["segments_half_open"] = self.segments
        out["frames"] = {"frame_index": self.frame_index.tolist(),
                         "time_ps": [None if math.isnan(t) else float(t) for t in self.time_ps]}
        return out

    def export(self, out_dir: str | os.PathLike, stem: str | None = None) -> tuple[Path, Path]:
        """Write ``<stem>.npz`` (arrays) and ``<stem>.json`` (mapping and provenance, incl. the npz SHA-256)."""
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        stem = stem or self.run_id
        npz_path, json_path = out / f"{stem}.npz", out / f"{stem}.json"
        np.savez(
            npz_path,
            coordinates_angstrom=self.coordinates,
            present=self.present,
            frame_index=self.frame_index,
            time_ps=self.time_ps,
            chain=np.array([r["chain"] for r in self.residues], dtype=str),
            resid=np.array([r["resid"] if r["resid"] is not None else -1 for r in self.residues], dtype=np.int64),
            icode=np.array([r["icode"] or "" for r in self.residues], dtype=str),
            label=np.array([r["label"] or "" for r in self.residues], dtype=str),
        )
        rep = self.to_report()
        rep["export"] = {"npz": npz_path.name, "npz_sha256": _sha256(npz_path),
                         "npz_arrays": {"coordinates_angstrom": "(frames, residues, 3) float64, Å, NaN if absent",
                                        "present": "(residues,) bool", "frame_index": "(frames,) int64",
                                        "time_ps": "(frames,) float64, NaN if not provided by the format",
                                        "chain/resid/icode/label": "(residues,) per output residue; resid -1 if absent"}}
        json_path.write_text(json.dumps(rep, indent=2, sort_keys=False, allow_nan=False) + "\n")
        return npz_path, json_path


def load_universe(topology: str | os.PathLike, *trajectories: str | os.PathLike, **kwargs: Any):
    """Thin local-file entry point: ``MDAnalysis.Universe(topology, *trajectories, **kwargs)``.

    Format support is exactly MDAnalysis's; this function adds nothing beyond recording what was read.
    """
    import MDAnalysis as mda

    for p in (topology, *trajectories):
        if not Path(p).is_file():
            raise FileNotFoundError(f"input file not found: {p}")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        u = mda.Universe(str(topology), *[str(t) for t in trajectories], **kwargs)
    u._topofold_load_warnings = list(dict.fromkeys(str(w.message) for w in caught))  # recorded in the report
    return u


def extract_ca(
    universe,
    run_id: str,
    *,
    selection: str = "protein",
    correspondence: Sequence[ResidueSpec | Mapping[str, Any]] | None = None,
    on_missing: str = "mark",
    chain_attribute: str = "auto",
    altloc: str = "reject",
    pbc: str = "none",
    max_ca_ca_angstrom: float = 4.2,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    assume_angstrom_if_undeclared: bool = False,
) -> CATrajectory:
    """Extract Cα coordinates (Å) for the residues of ``selection``, or for ``correspondence`` in its order.

    Parameters
    ----------
    universe : MDAnalysis.Universe
        Already loaded; see :func:`load_universe` for the file entry point.
    run_id : str
        Caller-chosen identifier of this run, stored with every output.
    selection : str
        MDAnalysis selection defining candidate residues (default ``"protein"``; avoids calcium ions named CA).
    correspondence : sequence of ResidueSpec or mappings with chain, resid, icode, label
        Explicit residue list and labels. Without it, all selected residues are returned in topology order with no
        label. Each entry must match at most one residue; a residue may be requested at most once.
    on_missing : {"mark", "reject"}
        For correspondence entries absent from the structure: return a NaN row with ``present=False``, or raise.
    chain_attribute : {"auto", "chainID", "segid"}
        "auto" uses chainIDs when every selected atom has a non-empty one, else segids.
    altloc : {"reject", "blank_then_A", "highest_occupancy"}
        Policy when a residue has several Cα atoms carrying altloc labels. "blank_then_A" takes the unlabelled atom,
        else "A". "highest_occupancy" takes the single highest-occupancy atom. Anything unresolved raises.
        Several Cα atoms without altloc labels always raise.
    pbc : {"none", "make_whole_sequential"}
        See module docstring.
    max_ca_ca_angstrom : float
        Consecutive Cα–Cα distances above this value in any frame are reported as geometric breaks.
    start, stop, step : int, optional
        Frame slice of the source trajectory.
    assume_angstrom_if_undeclared : bool
        Readers that declare no length unit (e.g. MDAnalysis MemoryReader) are rejected unless this is True, in which
        case the coordinates are taken to be in Å and the assumption is recorded.

    Fixed-topology contract
    -----------------------
    Atom identity comes from the topology and is assumed to be the same, in the same order, in every frame. Output
    order (topology order or ``correspondence``) is a separate, explicit choice. Where the format repeats atom identity
    per frame (multi-model PDB read by MDAnalysis's PDBReader), every model's identity records are checked against the
    first model and any difference raises. For coordinate-only formats (XTC, DCD, TRR, ...) and in-memory readers the
    order cannot be checked; this is reported as ``frame_identity.status = "not_verifiable"``. Arbitrary reordering of
    atoms between frames is not supported.
    """
    if altloc not in ALTLOC_POLICIES:
        raise AdapterError(f"unknown altloc policy {altloc!r}; expected one of {ALTLOC_POLICIES}")
    if pbc not in PBC_POLICIES:
        raise AdapterError(f"unknown pbc policy {pbc!r}; expected one of {PBC_POLICIES}")
    if on_missing not in ("mark", "reject"):
        raise AdapterError(f"unknown on_missing policy {on_missing!r}")
    if not run_id:
        raise AdapterError("run_id must be a non-empty string")
    if not (max_ca_ca_angstrom > 0 and math.isfinite(max_ca_ca_angstrom)):
        raise AdapterError("max_ca_ca_angstrom must be a positive finite number")

    import MDAnalysis as mda

    sel = universe.select_atoms(selection)
    if sel.n_atoms == 0:
        raise AdapterError(f"selection {selection!r} matched no atoms")
    chain_src = _chain_source(universe, sel, chain_attribute)

    # --- residue identity (fixed once, reused for every frame) ---------------------------------------------------
    topo_res = []  # one entry per topology residue in the selection, topology order
    seen_keys: dict[tuple[str, int, str], int] = {}
    altloc_decisions = []
    for res in sel.residues:
        atoms = res.atoms & sel
        chain = _chain_of(atoms, chain_src)
        icode = str(getattr(res, "icode", "") or "")
        key = (chain, int(res.resid), icode)
        if key in seen_keys:
            raise AdapterError(
                f"duplicate residue key chain={key[0]!r} resid={key[1]} icode={key[2]!r} (topology residues "
                f"{seen_keys[key]} and {int(res.ix)}); chain identity or insertion codes may have been lost by the "
                f"format. Use a topology that preserves them, or a different chain_attribute.")
        seen_keys[key] = int(res.ix)
        ca_atom, decision = _choose_ca(atoms.select_atoms("name CA"), key, altloc)
        if decision is not None:
            altloc_decisions.append(decision)
        topo_res.append({"key": key, "residue_ix": int(res.ix), "resname": str(res.resname),
                         "segid": str(res.segment.segid), "atom": ca_atom})

    # --- output order: topology order, or the caller's correspondence --------------------------------------------
    by_key = {r["key"]: r for r in topo_res}
    rows: list[dict[str, Any]] = []
    if correspondence is None:
        for r in topo_res:
            rows.append({"spec": None, "topo": r})
    else:
        specs = [s if isinstance(s, ResidueSpec) else ResidueSpec.from_mapping(s) for s in correspondence]
        used: set[tuple[str, int, str]] = set()
        for s in specs:
            if s.key in used:
                raise AdapterError(f"correspondence requests residue {s.key} more than once")
            used.add(s.key)
            r = by_key.get(s.key)
            if r is None and on_missing == "reject":
                raise AdapterError(f"correspondence residue {s.key} not found in selection {selection!r}")
            rows.append({"spec": s, "topo": r})

    residues, atom_ix = [], []
    for i, row in enumerate(rows):
        s, r = row["spec"], row["topo"]
        ca = None if r is None else r["atom"]
        status = "present" if ca is not None else ("residue_not_found" if r is None else "no_CA_atom")
        key = s.key if r is None else r["key"]
        residues.append({
            "index": i, "chain": key[0], "chain_source": chain_src, "resid": key[1], "icode": key[2],
            "label": None if s is None else s.label,
            "resname": None if r is None else r["resname"], "segid": None if r is None else r["segid"],
            "residue_ix": None if r is None else r["residue_ix"],
            "atom_ix": None if ca is None else int(ca.ix),
            "altloc": None if ca is None else (str(getattr(ca, "altLoc", "") or "")),
            "status": status,
        })
        atom_ix.append(-1 if ca is None else int(ca.ix))
    present = np.array([a >= 0 for a in atom_ix], dtype=bool)
    idx_present = np.array([a for a in atom_ix if a >= 0], dtype=np.int64)
    n_res = len(residues)

    # --- static breaks (topology) --------------------------------------------------------------------------------
    breaks: list[dict[str, Any]] = []
    for i in range(n_res - 1):
        a, b = residues[i], residues[i + 1]
        reasons = []
        if not (a["status"] == "present" and b["status"] == "present"):
            reasons.append("missing_CA")
        if a["chain"] != b["chain"] or (a["segid"] is not None and b["segid"] is not None and a["segid"] != b["segid"]):
            reasons.append("chain_boundary")
        elif not _consecutive_numbering(a, b):
            reasons.append("numbering_discontinuity")
        if reasons:
            breaks.append({"after_index": i, "reasons": reasons})

    # --- units and frame identity --------------------------------------------------------------------------------
    traj = universe.trajectory
    units = dict(getattr(traj, "units", {}) or {})
    unit_note = _check_length_unit(traj, units, assume_angstrom_if_undeclared)
    frame_identity = _check_frame_identity(traj)

    # --- frames ----------------------------------------------------------------------------------------------------
    time_available = units.get("time") is not None
    frames = list(range(len(traj)))[slice(start, stop, step)]
    if not frames:
        raise AdapterError("frame slice selected no frames")
    coords = np.full((len(frames), n_res, 3), np.nan, dtype=np.float64)
    times = np.full(len(frames), np.nan, dtype=np.float64)
    reader_dtype = None
    pbc_stats = {"policy": pbc, "frames_changed": 0, "max_shift_angstrom": 0.0}
    box_seen = False
    static_break_after = {b["after_index"] for b in breaks}
    seg_static = _segments(n_res, static_break_after, present)
    read_warnings: list[str] = []
    for f_out, f_src in enumerate(frames):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ts = traj[f_src]
        read_warnings += [str(w.message) for w in caught if str(w.message) not in read_warnings]
        pos = universe.atoms.positions[idx_present] if idx_present.size else np.empty((0, 3), np.float32)
        reader_dtype = reader_dtype or str(pos.dtype)
        xyz = np.asarray(pos, dtype=np.float64)
        bad = ~np.isfinite(xyz).all(axis=1)
        if bad.any():
            first = int(np.flatnonzero(present)[np.flatnonzero(bad)[0]])
            raise AdapterError(f"non-finite Cα coordinates in source frame {f_src} (output residue {first}, "
                               f"{residues[first]['chain']}:{residues[first]['resid']}{residues[first]['icode']})")
        coords[f_out, present] = xyz
        if time_available:
            times[f_out] = float(ts.time)
        box = ts.dimensions
        if box is not None and not np.any(np.asarray(box[:3], dtype=np.float64)):
            box = None  # MDAnalysis represents a missing box as zeros for some formats
        box_seen = box_seen or box is not None
        if pbc == "make_whole_sequential":
            coords[f_out] = _make_whole_segments(coords[f_out], seg_static, box, f_src, pbc_stats)

    # --- geometric breaks (any frame) --------------------------------------------------------------------------------
    geometric = []
    for i in range(n_res - 1):
        if i in static_break_after:
            continue
        d = np.linalg.norm(coords[:, i + 1] - coords[:, i], axis=1)
        over = np.flatnonzero(d > max_ca_ca_angstrom)
        if over.size:
            geometric.append(i)
            breaks.append({"after_index": i, "reasons": ["ca_ca_distance"],
                           "n_frames": int(over.size), "first_frames": frames_slice(frames, over[:10]),
                           "max_distance_angstrom": float(d.max())})
    breaks.sort(key=lambda b: b["after_index"])
    segments = _segments(n_res, static_break_after | set(geometric), present)

    fmt = _format_name(traj)
    report = {
        "schema": "topofold_inputs.ca_adapter/1",
        "engineering_note": "Engineering output; no biological or scientific validation is implied.",
        "units": {"coordinates": LENGTH_UNIT, "time": TIME_UNIT if time_available else None,
                  "source_length_unit": units.get("length"), "source_time_unit": units.get("time"),
                  "reader_convert_units": getattr(traj, "convert_units", None), "basis": unit_note},
        "frame_identity": frame_identity,
        "output_order": "correspondence" if correspondence is not None else "topology",
        "precision": {"reader_position_dtype": reader_dtype, "output_dtype": "float64",
                      "nominal_source_resolution": _NOMINAL_RESOLUTION.get(fmt, "unknown"),
                      "note": "float64 output does not restore precision lost in the source or the reader"},
        "time_source": "trajectory" if time_available else "not provided by the format (time_ps is NaN)",
        "config": {"selection": selection, "chain_attribute": chain_src, "chain_attribute_requested": chain_attribute,
                   "altloc": altloc, "pbc": pbc, "on_missing": on_missing,
                   "max_ca_ca_angstrom": max_ca_ca_angstrom, "frame_slice": [start, stop, step],
                   "correspondence_supplied": correspondence is not None,
                   "assume_angstrom_if_undeclared": assume_angstrom_if_undeclared},
        "altloc_decisions": altloc_decisions,
        "pbc": {**pbc_stats, "box_present_in_any_frame": box_seen},
        "summary": {"residues": n_res, "present": int(present.sum()), "missing": int((~present).sum()),
                    "breaks": len(breaks), "segments": len(segments), "frames": len(frames)},
        "provenance": {**_provenance(universe, fmt, mda.__version__), "read_warnings": read_warnings},
    }
    return CATrajectory(run_id=run_id, coordinates=coords, present=present,
                        frame_index=np.asarray(frames, dtype=np.int64), time_ps=times,
                        residues=residues, breaks=breaks, segments=segments, report=report)


# --- helpers -------------------------------------------------------------------------------------------------------

def frames_slice(frames: list[int], idx: Iterable[int]) -> list[int]:
    return [int(frames[int(i)]) for i in idx]


def _chain_source(universe, sel, requested: str) -> str:
    has_chain = hasattr(universe.atoms, "chainIDs")
    if requested == "chainID":
        if not has_chain:
            raise AdapterError("chain_attribute='chainID' but the topology has no chainIDs")
        return "chainID"
    if requested == "segid":
        return "segid"
    if requested != "auto":
        raise AdapterError(f"unknown chain_attribute {requested!r}")
    if has_chain and all(str(c).strip() for c in sel.chainIDs):
        return "chainID"
    return "segid"


def _chain_of(atoms, chain_src: str) -> str:
    vals = {str(v) for v in (atoms.chainIDs if chain_src == "chainID" else atoms.segids)}
    if len(vals) != 1:
        raise AdapterError(f"one residue spans several {chain_src} values: {sorted(vals)}")
    return vals.pop()


def _choose_ca(cas, key, policy: str):
    if cas.n_atoms == 0:
        return None, None
    if cas.n_atoms == 1:
        return cas[0], None
    alts = [str(getattr(a, "altLoc", "") or "") for a in cas]
    if not hasattr(cas, "altLocs") or sum(1 for a in alts if a == "") > 1 or len(set(alts)) != len(alts):
        raise AdapterError(f"residue {key}: {cas.n_atoms} Cα atoms without distinguishing altloc labels "
                           f"({alts}); residues may have been merged by the format")
    decision = {"chain": key[0], "resid": key[1], "icode": key[2], "available": alts, "policy": policy}
    if policy == "reject":
        raise AdapterError(f"residue {key}: alternate Cα locations {alts}; altloc policy is 'reject'")
    if policy == "blank_then_A":
        for want in ("", "A"):
            if want in alts:
                chosen = cas[alts.index(want)]
                decision["used"] = want
                return chosen, decision
        raise AdapterError(f"residue {key}: altlocs {alts} contain neither blank nor 'A'")
    occ = np.array([float(a.occupancy) for a in cas])
    best = np.flatnonzero(occ == occ.max())
    if best.size != 1:
        raise AdapterError(f"residue {key}: altlocs {alts} tie on occupancy {occ.max()}")
    decision["used"] = alts[int(best[0])]
    decision["occupancies"] = occ.tolist()
    return cas[int(best[0])], decision


def _next_icode(icode: str) -> str | None:
    """The insertion code that may follow ``icode`` within one residue number: "" -> "A" -> "B" ... -> "Z"."""
    if icode == "":
        return "A"
    if len(icode) == 1 and "A" <= icode < "Z":
        return chr(ord(icode) + 1)
    return None


def _consecutive_numbering(a: Mapping[str, Any], b: Mapping[str, Any]) -> bool:
    """Conservative sequence adjacency from numbering only (distance is never used to establish it).

    Adjacent only if b is the next insertion code of the same residue number (52 -> 52A -> 52B), or the next residue
    number without an insertion code (52 -> 53, 52B -> 53). Skipped, reversed, lowercase or other codes, a next number
    that starts with an insertion code (52 -> 53A) and every other jump are reported as discontinuities.
    """
    if b["resid"] == a["resid"]:
        return b["icode"] == _next_icode(a["icode"])
    return b["resid"] == a["resid"] + 1 and b["icode"] == ""


def _segments(n: int, break_after: set[int], present: np.ndarray) -> list[list[int]]:
    segs, start = [], None
    for i in range(n):
        if not present[i]:
            if start is not None:
                segs.append([start, i])
            start = None
            continue
        if start is None:
            start = i
        if i in break_after or i == n - 1:
            segs.append([start, i + 1])
            start = None
    return segs


def _make_whole_segments(X: np.ndarray, segments, box, frame: int, stats: dict) -> np.ndarray:
    if box is None:
        raise AdapterError(f"pbc='make_whole_sequential' but frame {frame} has no box")
    box = np.asarray(box, dtype=np.float64)
    if not np.allclose(box[3:6], 90.0, atol=1e-6) or not np.all(box[:3] > 0):
        raise AdapterError(f"frame {frame}: box {box.tolist()} is not a valid orthorhombic box; unsupported for make_whole")
    L = box[:3]
    Y = X.copy()
    for s, e in segments:  # never across a break
        for i in range(s + 1, e):
            Y[i] -= L * np.round((Y[i] - Y[i - 1]) / L)
    shift = float(np.nanmax(np.abs(Y - X))) if np.isfinite(X).any() else 0.0
    if shift > 0:
        stats["frames_changed"] += 1
        stats["max_shift_angstrom"] = max(stats["max_shift_angstrom"], shift)
    return Y


_ANGSTROM_NAMES = {"angstrom", "a", "å"}


def _check_length_unit(traj, units: Mapping[str, Any], assume_angstrom: bool) -> str:
    """Return the basis for treating positions as Å, or raise."""
    length = units.get("length")
    convert = getattr(traj, "convert_units", None)
    if length is None:
        if assume_angstrom:
            return "reader declares no length unit; Å assumed by caller (assume_angstrom_if_undeclared=True)"
        raise AdapterError(f"{type(traj).__name__} declares no length unit; pass assume_angstrom_if_undeclared=True "
                           f"only if the coordinates are known to be in Å")
    if str(length).lower() in _ANGSTROM_NAMES:
        return f"source unit {length}; no conversion needed"
    if convert is True:
        return f"source unit {length}; converted to Å by MDAnalysis (convert_units=True)"
    raise AdapterError(f"positions are in {length} (reader convert_units={convert}); load with unit conversion "
                       f"enabled so that coordinates are in Å")


def _check_frame_identity(traj) -> dict[str, Any]:
    """Verify per-frame atom identity where the format records it (multi-model PDB); otherwise disclose."""
    n = len(traj)
    name = type(traj).__name__
    if n == 1:
        return {"status": "single_frame", "reader": name}
    path = getattr(traj, "filename", None)
    if name == "PDBReader" and path and Path(str(path)).is_file():
        models, current = [], None
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                rec = line[:6]
                if rec == "MODEL ":
                    current = []
                    models.append(current)
                elif rec in ("ATOM  ", "HETATM"):
                    if current is None:
                        current = []
                        models.append(current)
                    # identity columns only (name, altloc, resname, chain, resSeq, iCode, segID); not coordinates
                    current.append((line[12:16], line[16], line[17:20], line[21], line[22:26], line[26],
                                    line[72:76].rstrip("\n")))
                elif rec == "ENDMDL":
                    current = None
        if len(models) != n:
            raise AdapterError(f"PDB has {len(models)} models but the reader reports {n} frames; cannot verify identity")
        ref = models[0]
        for k, m in enumerate(models[1:], start=2):
            if m != ref:
                diff = next((i for i, (x, y) in enumerate(zip(ref, m)) if x != y), min(len(ref), len(m)))
                raise AdapterError(f"atom identity differs between model 1 and model {k} at record {diff + 1}; "
                                   f"frames do not share one fixed topology")
        return {"status": "verified", "reader": name, "models_checked": n,
                "fields": ["name", "altloc", "resname", "chain", "resSeq", "iCode", "segID"]}
    return {"status": "not_verifiable", "reader": name,
            "note": "format stores no per-frame atom identity; fixed atom order is assumed (fixed-topology contract)"}


def _format_name(traj) -> str:
    fmt = getattr(traj, "format", None)
    if isinstance(fmt, (list, tuple)):
        fmt = fmt[0] if fmt else None
    return str(fmt).upper() if fmt else type(traj).__name__


def _source_files(universe) -> list[dict[str, Any]]:
    files = []
    top = getattr(universe, "filename", None)
    if top:
        files.append(("topology", top))
    traj = universe.trajectory
    names = getattr(traj, "filenames", None) or [getattr(traj, "filename", None)]
    for n in names:
        if n and all(n != f for _, f in files):
            files.append(("trajectory", n))
    out = []
    for role, name in files:
        p = Path(str(name))
        if p.is_file():
            out.append({"role": role, "path": str(p.resolve()), "bytes": p.stat().st_size, "sha256": _sha256(p)})
        else:
            out.append({"role": role, "path": str(name), "bytes": None, "sha256": None,
                        "note": "not a local file; content not hashed"})
    return out


def _provenance(universe, fmt: str, mda_version: str) -> dict[str, Any]:
    here = Path(__file__)
    return {
        "adapter": {"name": "topofold_inputs.ca_adapter", "version": ADAPTER_VERSION,
                    "source_sha256": _sha256(here)},
        "reader": {"class": type(universe.trajectory).__name__, "format": fmt, "n_frames_source": len(universe.trajectory)},
        "source_files": _source_files(universe),
        "load_warnings": list(getattr(universe, "_topofold_load_warnings", [])),
        "software": {"python": platform.python_version(), "numpy": np.__version__, "MDAnalysis": mda_version},
        "created_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
    }


def _sha256(path: str | os.PathLike) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
