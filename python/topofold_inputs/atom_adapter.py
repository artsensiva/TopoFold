"""Bounded extraction of named backbone and side-chain atoms, with explicit missing-atom and altloc policies.

Engineering capability ("B"), not a validated analysis. It extends, without modifying, the Cα adapter
(:mod:`topofold_inputs.ca_adapter`): residue identity, chain attribution, unit checks, frame-identity checks and
provenance are reused from there, so Cα outputs of that module are unaffected.

What it returns is a flat list of requested atoms, each identified by (chain, resid, icode, atom name), with
coordinates in Å for every frame. Nothing is inferred:

* **Requested atoms only.** The caller names residues (explicit correspondence, required) and atom names per residue.
  :data:`KINASE_DESCRIPTOR_ROLES` lists the atoms the DFG / αC / activation-loop descriptors need; the caller maps each
  role to a residue of its structure. No kinase numbering is built in.
* **Missing atoms.** ``on_missing="mark"`` returns NaN coordinates with a per-atom status
  (``residue_not_found`` or ``atom_not_found``); ``"reject"`` raises :class:`AdapterError`. A missing atom is never filled,
  rebuilt or substituted.
* **Alternate locations, per residue.** ``reject`` raises if any requested atom of a residue has an altloc.
  ``blank_then_A`` takes each atom's blank location, else ``A``. ``highest_occupancy`` chooses one altloc identifier for
  the whole residue (highest mean occupancy over its requested atoms that have altlocs; a tie raises) and each atom takes
  its blank location, else that identifier. One residue never mixes conformers. Every decision is recorded.
* **Units and precision** as in the Cα adapter: Å, float64 output that cannot restore source precision.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .ca_adapter import (ALTLOC_POLICIES, LENGTH_UNIT, TIME_UNIT, AdapterError, ResidueSpec, _check_frame_identity,
                         _check_length_unit, _chain_of, _chain_source, _format_name, _provenance, _sha256)

ATOM_ADAPTER_VERSION = "0.1.0"
BACKBONE = ("N", "CA", "C", "O")

# Atoms each kinase-descriptor role needs (Modi & Dunbrack, PNAS 116:6818, 2019 definitions; see kinase_geometry.py).
# The caller maps every role to one residue of its structure; residue numbering is structure-specific.
KINASE_DESCRIPTOR_ROLES: dict[str, tuple[str, ...]] = {
    "XDFG_minus1": ("C",),                # supplies C for phi of XDFG
    "beta3_Lys": ("CA", "CB"),            # D2 uses CA; alphaC criterion uses CB
    "alphaC_Glu": ("CB",),                # alphaC Cβ–Cβ criterion
    "alphaC_Glu_plus4": ("CA",),          # D1
    "XDFG": ("N", "CA", "C"),             # phi/psi of the residue before DFG; C needed for phi of DFG-Asp
    "DFG_Asp": ("N", "CA", "C"),          # phi/psi
    "DFG_Phe": ("N", "CA", "C", "CB", "CG", "CZ"),  # phi/psi, chi1 (N-CA-CB-CG), D1/D2 (CZ)
    "DFG_Gly": ("N",),                    # psi of DFG-Phe
    "DFGxxX": ("N",),                     # extended A-loop H-bond donor (6th residue counting DFG-Asp as 1)
    "XHRD": ("O",),                       # extended A-loop H-bond acceptor (residue before HRD)
}


@dataclass
class AtomTrajectory:
    """Named atoms of one run: coordinates (F, M, 3) in Å for M requested (residue, atom) entries."""

    run_id: str
    coordinates: np.ndarray  # (F, M, 3) float64, Å; NaN where status != "present"
    atoms: list[dict[str, Any]]  # one record per requested atom, in request order
    frame_index: np.ndarray
    time_ps: np.ndarray
    report: dict[str, Any] = field(default_factory=dict)

    def index(self, chain: str, resid: int, name: str, icode: str = "") -> int:
        for i, a in enumerate(self.atoms):
            if a["chain"] == chain and a["resid"] == resid and a["icode"] == icode and a["name"] == name:
                return i
        raise KeyError(f"{chain}:{resid}{icode} {name} was not requested")

    def xyz(self, chain: str, resid: int, name: str, icode: str = "") -> np.ndarray:
        """(F, 3) coordinates of one requested atom (NaN if absent)."""
        return self.coordinates[:, self.index(chain, resid, name, icode)]

    def export(self, out_dir: str | os.PathLike, stem: str | None = None) -> tuple[Path, Path]:
        import json

        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        stem = stem or self.run_id
        npz_path, json_path = out / f"{stem}.npz", out / f"{stem}.json"
        np.savez(npz_path, coordinates_angstrom=self.coordinates, frame_index=self.frame_index, time_ps=self.time_ps,
                 present=np.array([a["status"] == "present" for a in self.atoms], dtype=bool))
        rep = dict(self.report)
        rep.update(run_id=self.run_id, atoms=self.atoms,
                   export={"npz": npz_path.name, "npz_sha256": _sha256(npz_path),
                           "npz_arrays": {"coordinates_angstrom": "(frames, atoms, 3) float64, Å, NaN if absent",
                                          "present": "(atoms,) bool", "frame_index": "(frames,) int64",
                                          "time_ps": "(frames,) float64, NaN if not provided"}})
        json_path.write_text(json.dumps(rep, indent=2, allow_nan=False) + "\n")
        return npz_path, json_path


def role_request(roles: Mapping[str, ResidueSpec | Mapping[str, Any]],
                 role_atoms: Mapping[str, Sequence[str]] = KINASE_DESCRIPTOR_ROLES) -> list[tuple[ResidueSpec, list[str]]]:
    """Turn {role: residue} into an atom request, merging roles that point at the same residue. Every role in
    ``role_atoms`` must be mapped; unknown roles raise."""
    unknown = set(roles) - set(role_atoms)
    if unknown:
        raise AdapterError(f"unknown roles {sorted(unknown)}")
    missing = set(role_atoms) - set(roles)
    if missing:
        raise AdapterError(f"roles not mapped to a residue: {sorted(missing)}")
    merged: dict[tuple[str, int, str], tuple[ResidueSpec, list[str]]] = {}
    for role, spec in roles.items():
        s = spec if isinstance(spec, ResidueSpec) else ResidueSpec.from_mapping(spec)
        names = merged.setdefault(s.key, (s, []))[1]
        for n in role_atoms[role]:
            if n not in names:
                names.append(n)
    return list(merged.values())


def extract_atoms(
    universe,
    run_id: str,
    request: Sequence[tuple[ResidueSpec | Mapping[str, Any], Sequence[str]]],
    *,
    selection: str = "protein",
    on_missing: str = "mark",
    chain_attribute: str = "auto",
    altloc: str = "reject",
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    assume_angstrom_if_undeclared: bool = False,
) -> AtomTrajectory:
    """Coordinates (Å) of the requested atoms, in request order, for every frame of the slice.

    ``request`` is a sequence of (residue, atom names). Each residue key and each (residue, atom) pair may appear only
    once. See the module docstring for the missing-atom and altloc policies.
    """
    if altloc not in ALTLOC_POLICIES:
        raise AdapterError(f"unknown altloc policy {altloc!r}; expected one of {ALTLOC_POLICIES}")
    if on_missing not in ("mark", "reject"):
        raise AdapterError(f"unknown on_missing policy {on_missing!r}")
    if not run_id:
        raise AdapterError("run_id must be a non-empty string")
    if not request:
        raise AdapterError("empty atom request")

    import MDAnalysis as mda

    sel = universe.select_atoms(selection)
    if sel.n_atoms == 0:
        raise AdapterError(f"selection {selection!r} matched no atoms")
    chain_src = _chain_source(universe, sel, chain_attribute)
    residues_by_key: dict[tuple[str, int, str], Any] = {}
    for res in sel.residues:
        atoms = res.atoms & sel
        key = (_chain_of(atoms, chain_src), int(res.resid), str(getattr(res, "icode", "") or ""))
        if key in residues_by_key:
            raise AdapterError(f"duplicate residue key {key}; chain identity or insertion codes may have been lost")
        residues_by_key[key] = atoms

    specs = []
    seen = set()
    for spec, names in request:
        s = spec if isinstance(spec, ResidueSpec) else ResidueSpec.from_mapping(spec)
        if s.key in seen:
            raise AdapterError(f"residue {s.key} requested more than once; merge its atom names")
        seen.add(s.key)
        if len(set(names)) != len(list(names)) or not names:
            raise AdapterError(f"residue {s.key}: atom names must be non-empty and unique, got {list(names)}")
        specs.append((s, [str(n) for n in names]))

    records, atom_ix, decisions = [], [], []
    for s, names in specs:
        res_atoms = residues_by_key.get(s.key)
        if res_atoms is None:
            if on_missing == "reject":
                raise AdapterError(f"residue {s.key} not found in selection {selection!r}")
            for n in names:
                records.append(_record(s, n, None, "residue_not_found", None))
                atom_ix.append(-1)
            continue
        chosen, decision = _choose_residue_atoms(res_atoms, s.key, names, altloc)
        if decision is not None:
            decisions.append(decision)
        for n in names:
            a = chosen.get(n)
            if a is None and on_missing == "reject":
                raise AdapterError(f"residue {s.key} ({str(res_atoms[0].resname)}): atom {n} not found")
            records.append(_record(s, n, a, "present" if a is not None else "atom_not_found",
                                   str(res_atoms[0].resname)))
            atom_ix.append(-1 if a is None else int(a.ix))

    present = np.array([i >= 0 for i in atom_ix], dtype=bool)
    idx = np.array([i for i in atom_ix if i >= 0], dtype=np.int64)
    traj = universe.trajectory
    units = dict(getattr(traj, "units", {}) or {})
    unit_note = _check_length_unit(traj, units, assume_angstrom_if_undeclared)
    frame_identity = _check_frame_identity(traj)
    frames = list(range(len(traj)))[slice(start, stop, step)]
    if not frames:
        raise AdapterError("frame slice selected no frames")
    time_available = units.get("time") is not None
    coords = np.full((len(frames), len(records), 3), np.nan, dtype=np.float64)
    times = np.full(len(frames), np.nan, dtype=np.float64)
    read_warnings: list[str] = []
    for f_out, f_src in enumerate(frames):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ts = traj[f_src]
        read_warnings += [str(w.message) for w in caught if str(w.message) not in read_warnings]
        xyz = np.asarray(universe.atoms.positions[idx], dtype=np.float64) if idx.size else np.empty((0, 3))
        if not np.isfinite(xyz).all():
            raise AdapterError(f"non-finite atom coordinates in source frame {f_src}")
        coords[f_out, present] = xyz
        if time_available:
            times[f_out] = float(ts.time)

    fmt = _format_name(traj)
    report = {
        "schema": "topofold_inputs.atom_adapter/1",
        "engineering_note": "Engineering output; no biological or scientific validation is implied.",
        "units": {"coordinates": LENGTH_UNIT, "time": TIME_UNIT if time_available else None,
                  "source_length_unit": units.get("length"), "basis": unit_note},
        "frame_identity": frame_identity,
        "config": {"selection": selection, "chain_attribute": chain_src, "altloc": altloc, "on_missing": on_missing,
                   "frame_slice": [start, stop, step], "assume_angstrom_if_undeclared": assume_angstrom_if_undeclared},
        "altloc_decisions": decisions,
        "summary": {"atoms_requested": len(records), "present": int(present.sum()),
                    "missing": int((~present).sum()), "frames": len(frames),
                    "missing_by_status": {st: sum(1 for r in records if r["status"] == st)
                                          for st in ("residue_not_found", "atom_not_found")}},
        "provenance": {**_provenance(universe, fmt, mda.__version__), "read_warnings": read_warnings,
                       "atom_adapter_version": ATOM_ADAPTER_VERSION},
    }
    return AtomTrajectory(run_id=run_id, coordinates=coords, atoms=records,
                          frame_index=np.asarray(frames, dtype=np.int64), time_ps=times, report=report)


def _record(s: ResidueSpec, name: str, atom, status: str, resname: str | None) -> dict[str, Any]:
    return {"chain": s.chain, "resid": s.resid, "icode": s.icode, "label": s.label, "name": name, "resname": resname,
            "atom_ix": None if atom is None else int(atom.ix),
            "altloc": None if atom is None else str(getattr(atom, "altLoc", "") or ""), "status": status}


def _choose_residue_atoms(res_atoms, key, names: Sequence[str], policy: str):
    """{name: atom or None} for one residue under a residue-level altloc policy, plus the decision record."""
    by_name: dict[str, list] = {}
    for a in res_atoms:
        if a.name in names:
            by_name.setdefault(str(a.name), []).append(a)
    alt_of = lambda a: str(getattr(a, "altLoc", "") or "")  # noqa: E731
    for n, lst in by_name.items():
        alts = [alt_of(a) for a in lst]
        if len(set(alts)) != len(alts):
            raise AdapterError(f"residue {key}: atom {n} appears {len(lst)} times without distinguishing altloc "
                               f"labels ({alts})")
    labelled = {alt_of(a) for lst in by_name.values() for a in lst} - {""}
    if not labelled:
        return {n: (by_name[n][0] if n in by_name else None) for n in names}, None
    decision = {"chain": key[0], "resid": key[1], "icode": key[2], "policy": policy,
                "atoms_with_altlocs": sorted(n for n, lst in by_name.items() if any(alt_of(a) for a in lst)),
                "available": sorted(labelled)}
    if policy == "reject":
        raise AdapterError(f"residue {key}: alternate locations {sorted(labelled)} on requested atoms "
                           f"{decision['atoms_with_altlocs']}; altloc policy is 'reject'")
    if policy == "blank_then_A":
        want = "A"
    else:
        occ = {alt: np.mean([float(a.occupancy) for lst in by_name.values() for a in lst if alt_of(a) == alt])
               for alt in labelled}
        best = max(occ.values())
        winners = sorted(alt for alt, o in occ.items() if o == best)
        if len(winners) != 1:
            raise AdapterError(f"residue {key}: altlocs {winners} tie on mean occupancy {best}")
        want = winners[0]
        decision["mean_occupancies"] = {k: float(v) for k, v in sorted(occ.items())}
    chosen = {}
    for n in names:
        lst = by_name.get(n, [])
        pick = next((a for a in lst if alt_of(a) == ""), None) or next((a for a in lst if alt_of(a) == want), None)
        chosen[n] = pick
    if policy == "blank_then_A" and any(by_name.get(n) and chosen[n] is None for n in names):
        raise AdapterError(f"residue {key}: altlocs {sorted(labelled)} contain neither blank nor 'A' for some atoms")
    decision["used"] = want
    return chosen, decision


__all__ = ["ATOM_ADAPTER_VERSION", "BACKBONE", "KINASE_DESCRIPTOR_ROLES", "AtomTrajectory", "extract_atoms",
           "role_request"]
