"""Descriptive comparison of two experimental structures on an explicit residue correspondence (static, single model).

    PYTHONPATH=python .venv/bin/python -m topofold_inputs.structpair SPEC.json --out NEW_DIR [--adapter-python PATH]

SPEC.json (relative paths are resolved against the directory of SPEC.json):
    {"name": "...",
     "structures": {"A": {"label": "...", "file": "a.pdb", "chain": "A", "altloc": "blank_then_A"},
                    "B": {"label": "...", "file": "b.pdb", "chains": ["A", "B"], "primary_chain": "A", "altloc": "reject"}},
     "correspondence": {"same_numbering": "1-59,68-166"}            # identical author numbering in both chains, or
     "correspondence": {"pairs": [[1, 1], [2, 2], [3, "3A"], ...]},  # explicit A-residue -> B-residue pairs
     "on_missing": "reject",            # or "exclude": absent residues are kept as missing rows and their windows excluded
     "resname_mismatch": "reject",      # or "flag": allowed (e.g. a point mutation), flagged in every table
     "altloc_sensitivity": "highest_occupancy",   # optional second altloc policy, applied to both structures
     "modified_residues": "phospho_pdb",   # optional: admit SEP, TPO, PTR documented by MODRES records (PDB files only); default "none"
     "numeric_record": true,            # optional: also write numeric_record.json, the computed values before six-figure display rounding
     "top_n": 10}

Nothing is chosen for the user: every chain is named in the spec, with one primary chain per structure ("chain", or
"chains" plus "primary_chain" for within-entry crystal-copy variation, see docs/STRUCTPAIR_CHAINS_SPEC_2026-10-25.md); the
correspondence is never inferred (no sequence alignment, no chain ranking). Each structure is read by the existing Cα adapter (MDAnalysis), κ/τ come from the
TopoFold core via ``features``, window Cα distances from ``compare.window_distances``, local and global RMSD from
MDAnalysis QCP (``qcp_rmsd``), and exclusion causes from ``ensemble.exclusion_causes``. Outputs are descriptive only:
no p-values, populations, kinetics or functional states can be derived from two static structures.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import importlib.util
import json
import math
import os
import platform
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np

from . import chain_variation as cv
from . import numeric_record as nr
from .compare import spearman, window_distances
from .ensemble import _run, _sha256, exclusion_causes

STRUCTPAIR_VERSION = "0.5.0"   # 0.5.0: opt-in modified_residues="phospho_pdb" (SEP, TPO, PTR; PDB format); 0.4.0: opt-in numeric_record.json (full-precision window values); 0.3.1: an occupancy tie in the OPTIONAL altloc sensitivity is reported as unavailable, not fatal
SCHEMA = "topofold_inputs.structpair/3"   # /2: 0.2.0 (advisory/1 semantics); /1: 0.1.0. Older records keep their meaning
METRICS = ("d_kappa", "d_i_i2", "rmsd3", "d_tau", "d_i_i3", "rmsd4")
FOUR = ("d_tau", "d_i_i3", "rmsd4")
UNITS = {"d_kappa": "rad", "d_tau": "rad", "d_i_i2": "angstrom", "d_i_i3": "angstrom", "rmsd3": "angstrom",
         "rmsd4": "angstrom", "kappa": "rad", "tau": "rad", "d2": "angstrom", "d3": "angstrom"}
DESCRIPTION = {
    "d_kappa": "|κ_B − κ_A|, three-residue window (i−1, i, i+1), TopoFold core",
    "d_i_i2": "|d(i−1,i+1)_B − d(i−1,i+1)_A|, same three-residue window (distance baseline)",
    "rmsd3": "proper-rotation Cα RMSD of the three-residue window after optimal superposition (MDAnalysis QCP)",
    "d_tau": "circular |wrap(τ_B − τ_A)|, four-residue window (i−1 .. i+2), hinge i→i+1, TopoFold core",
    "d_i_i3": "|d(i−1,i+2)_B − d(i−1,i+2)_A|, same four-residue window (distance baseline)",
    "rmsd4": "proper-rotation Cα RMSD of the four-residue window after optimal superposition (MDAnalysis QCP)"}
COLUMNS = [("kappa_A", "rad"), ("kappa_B", "rad"), ("d2_A", "angstrom"), ("d2_B", "angstrom"), ("d_kappa", "rad"),
           ("d_i_i2", "angstrom"), ("rmsd3", "angstrom"), ("tau_A", "rad"), ("tau_B", "rad"), ("d3_A", "angstrom"),
           ("d3_B", "angstrom"), ("d_tau", "rad"), ("d_i_i3", "angstrom"), ("rmsd4", "angstrom")]
AGREEMENT_PAIRS = [("d_kappa", "d_i_i2"), ("d_kappa", "rmsd3"), ("d_i_i2", "rmsd3"),
                   ("d_tau", "d_i_i3"), ("d_tau", "rmsd4"), ("d_i_i3", "rmsd4")]
ALTLOC_POLICIES = ("reject", "blank_then_A", "highest_occupancy")
MODIFIED_POLICIES = ("none", "phospho_pdb")
LIMITS = ("Two static structures: the values are direct geometric differences between these coordinate sets. They can "
          "reflect crystal form and packing, crystallization conditions, resolution, refinement, alternate "
          "conformations and coordinate error as well as any biological difference. No ensemble, population, kinetic, "
          "functional-state or treatment-response claim follows. Agreement or disagreement between metrics does not show "
          "that one representation is better. Crystal copies of one entry are not independent replicates.")
COPY_NOTE = cv.STATEMENT


class SpecError(ValueError):
    pass


# --- spec ------------------------------------------------------------------------------------------------------------

def _residue(x) -> tuple[int, str]:
    """1, "1", "1A" or {"resid": 1, "icode": "A"} -> (resid, icode)."""
    if isinstance(x, dict):
        return int(x["resid"]), str(x.get("icode", "") or "")
    m = re.fullmatch(r"\s*(-?\d+)([A-Za-z]?)\s*", str(x))
    if not m:
        raise SpecError(f"cannot read residue {x!r}")
    return int(m.group(1)), m.group(2)


def _ranges(text: str) -> list[int]:
    out = []
    for part in str(text).split(","):
        m = re.fullmatch(r"\s*(-?\d+)\s*(?:-\s*(-?\d+)\s*)?", part)
        if not m:
            raise SpecError(f"cannot read residue range {part!r}")
        lo = int(m.group(1))
        hi = int(m.group(2)) if m.group(2) is not None else lo
        if hi < lo:
            raise SpecError(f"descending range {part!r}")
        out += list(range(lo, hi + 1))
    return out


def load_spec(path: Path) -> dict:
    spec = json.loads(Path(path).read_text())
    base = Path(path).resolve().parent
    if set(spec.get("structures", {})) != {"A", "B"}:
        raise SpecError("spec.structures must contain exactly 'A' and 'B'")
    for key, s in spec["structures"].items():
        if not s.get("file"):
            raise SpecError(f"structures.{key}.file is required")
        if "chains" in s:
            if s.get("chain"):
                raise SpecError(f"structures.{key}: give either chain or chains, not both")
            ch = s["chains"]
            if not isinstance(ch, list) or not ch or not all(isinstance(c, str) and c.strip() for c in ch):
                raise SpecError(f"structures.{key}.chains must be a non-empty list of chain names")
            if len(set(ch)) != len(ch):
                raise SpecError(f"structures.{key}.chains repeats a chain")
            prim = s.get("primary_chain")
            if prim is None and len(ch) == 1:
                prim = ch[0]
            if prim not in ch:
                raise SpecError(f"structures.{key}.primary_chain is required with more than one chain and must be one of "
                                f"chains (no primary chain is chosen for the user)")
            s["chain"], s["chains"] = prim, list(ch)
        else:
            if not s.get("chain"):
                raise SpecError(f"structures.{key}.chain is required (no default chain is chosen); or name every chain in "
                                f"structures.{key}.chains with primary_chain")
            if "primary_chain" in s:
                raise SpecError(f"structures.{key}.primary_chain needs structures.{key}.chains")
            s["chains"] = [s["chain"]]
        s.setdefault("label", key)
        s.setdefault("altloc", "reject")
        if s["altloc"] not in ALTLOC_POLICIES:
            raise SpecError(f"structures.{key}.altloc must be one of {ALTLOC_POLICIES}")
        p = Path(s["file"])
        s["path"] = str(p if p.is_absolute() else (base / p).resolve())
        if not Path(s["path"]).is_file():
            raise SpecError(f"structures.{key}.file not found: {s['path']}")
    corr = spec.get("correspondence")
    if not isinstance(corr, dict) or len(corr) != 1 or next(iter(corr)) not in ("same_numbering", "pairs"):
        raise SpecError("spec.correspondence must be {'same_numbering': 'ranges'} or {'pairs': [[a, b], ...]} "
                        "(the correspondence is never inferred)")
    if "same_numbering" in corr:
        pairs = [((r, ""), (r, "")) for r in _ranges(corr["same_numbering"])]
    else:
        pairs = [(_residue(a), _residue(b)) for a, b in corr["pairs"]]
    if len(pairs) < 4:
        raise SpecError("the correspondence needs at least 4 residues")
    for side in (0, 1):
        keys = [p[side] for p in pairs]
        if len(set(keys)) != len(keys):
            raise SpecError(f"correspondence repeats a residue of structure {'AB'[side]}")
    spec["_pairs"] = pairs
    spec.setdefault("on_missing", "reject")
    spec.setdefault("resname_mismatch", "reject")
    spec.setdefault("altloc_sensitivity", None)
    spec.setdefault("top_n", 10)
    spec.setdefault("numeric_record", False)
    spec.setdefault("modified_residues", "none")
    if spec["modified_residues"] not in MODIFIED_POLICIES:
        raise SpecError(f"modified_residues must be one of {MODIFIED_POLICIES}")
    if not isinstance(spec["numeric_record"], bool):
        raise SpecError("numeric_record must be true or false")
    try:
        spec["_variation"] = cv.parse_params(spec.get("chain_variation"))
    except cv.VariationError as e:
        raise SpecError(str(e)) from None
    if spec["on_missing"] not in ("reject", "exclude"):
        raise SpecError("on_missing must be 'reject' or 'exclude'")
    if spec["resname_mismatch"] not in ("reject", "flag"):
        raise SpecError("resname_mismatch must be 'reject' or 'flag'")
    if spec["altloc_sensitivity"] not in (None, *ALTLOC_POLICIES):
        raise SpecError(f"altloc_sensitivity must be null or one of {ALTLOC_POLICIES}")
    return spec


# --- extraction (existing components) --------------------------------------------------------------------------------

class AltlocTieRefusal(RuntimeError):
    """The adapter refused a structure because alternate locations are tied on occupancy under `highest_occupancy`.

    `ties` is the adapter's structured record (every tied residue), read from its error file, not parsed from text."""

    def __init__(self, message: str, ties: list):
        super().__init__(message)
        self.ties = ties


def extract(spec: dict, key: str, altloc: str, work: Path, adapter_python: str, tag: str, chain: str | None = None) -> dict:
    from .features import export_features, features_from_export

    s = spec["structures"][key]
    side = 0 if key == "A" else 1
    chain = chain or s["chain"]
    suffix = "" if chain == s["chain"] else f"_{chain}"
    corr = work / f"{tag}_{key}{suffix}_correspondence.json"
    corr.write_text(json.dumps([{"chain": chain, "resid": p[side][0], "icode": p[side][1], "label": f"pos{k}"}
                                for k, p in enumerate(spec["_pairs"])]))
    run_id = f"{tag}_{key}{suffix}"
    on_missing = "reject" if spec["on_missing"] == "reject" else "mark"
    err = work / f"{run_id}_error.json"
    err.unlink(missing_ok=True)
    try:
        _run(["-m", "topofold_inputs", s["path"], "--run-id", run_id, "--out", str(work), "--correspondence", str(corr),
              "--on-missing", on_missing, "--altloc", altloc, "--error-json", str(err),
              *(["--modified-residues", spec["modified_residues"]] if spec["modified_residues"] != "none" else [])], adapter_python)
    except RuntimeError as e:
        rec = None
        if err.is_file():
            try:
                rec = json.loads(err.read_text())
            except ValueError:
                rec = None
        if isinstance(rec, dict) and rec.get("kind") == "altloc_tie" and isinstance(rec.get("ties"), list) and rec["ties"]:
            raise AltlocTieRefusal(str(e), rec["ties"]) from e
        raise
    npz, js = work / f"{run_id}.npz", work / f"{run_id}.json"
    feats, rep = features_from_export(npz, js)
    export_features(feats, npz, js, work, run_id)
    a = np.load(npz)
    if a["coordinates_angstrom"].shape[0] != 1:
        raise SpecError(f"structure {key}: expected one model, found {a['coordinates_angstrom'].shape[0]}")
    fdoc = json.loads((work / f"{run_id}_features.json").read_text())
    return {"run_id": run_id, "npz": npz, "coords": a["coordinates_angstrom"][0].astype(np.float64),
            "frame_index": a["frame_index"], "kappa": feats.kappa[0], "tau": feats.tau[0],
            "kappa_valid": feats.kappa_valid[0], "tau_valid": feats.tau_valid[0],
            "kappa_reason": feats.kappa_reason, "tau_reason": feats.tau_reason, "residues": rep["residues"],
            "breaks": rep["breaks"], "altloc_decisions": rep["altloc_decisions"], "modified_block": rep.get("modified_residues"),
            "provenance": {"adapter_json_sha256": _sha256(js), "features_json_sha256": _sha256(work / f"{run_id}_features.json"),
                           "source_files": rep["provenance"]["source_files"], "adapter": rep["provenance"]["adapter"],
                           "adapter_software": rep["provenance"]["software"], "core": fdoc["provenance"]["core"],
                           "features": fdoc["provenance"]["features"], "units": rep["units"],
                           "frame_identity": rep["frame_identity"]["status"]}}


# --- geometry ---------------------------------------------------------------------------------------------------------

def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Independent numpy check of the QCP global value: SVD fit with determinant correction (proper rotation)."""
    R, cP, cQ = cv.kabsch_frame(P, Q)
    return float(np.sqrt(np.mean(np.sum(((P - cP) @ R.T + cQ - Q) ** 2, axis=1))))


def window_metrics(a: dict, b: dict, rm: dict) -> dict:
    """Per-window differences on identical coordinates; a window gets values only where κ (three-residue) or τ
    (four-residue) is valid in both structures, so every metric of a window size shares one residue coverage."""
    n = len(a["kappa"])
    d2a, d2b = window_distances(a["coords"][None], 2)[0], window_distances(b["coords"][None], 2)[0]
    d3a, d3b = window_distances(a["coords"][None], 3)[0], window_distances(b["coords"][None], 3)[0]
    m = {k: [None] * n for k in (*METRICS, "kappa_A", "kappa_B", "tau_A", "tau_B", "d2_A", "d2_B", "d3_A", "d3_B")}
    for i in range(n):
        if 1 <= i <= n - 2 and a["kappa_valid"][i] and b["kappa_valid"][i]:
            m["kappa_A"][i], m["kappa_B"][i] = float(a["kappa"][i]), float(b["kappa"][i])
            m["d2_A"][i], m["d2_B"][i] = float(d2a[i]), float(d2b[i])
            m["d_kappa"][i] = abs(float(b["kappa"][i] - a["kappa"][i]))
            m["d_i_i2"][i] = abs(float(d2b[i] - d2a[i]))
            m["rmsd3"][i] = rm["rmsd3"][i]
        if 1 <= i <= n - 3 and a["tau_valid"][i] and b["tau_valid"][i]:
            dt_ = float(b["tau"][i] - a["tau"][i])
            m["tau_A"][i], m["tau_B"][i] = float(a["tau"][i]), float(b["tau"][i])
            m["d3_A"][i], m["d3_B"][i] = float(d3a[i]), float(d3b[i])
            m["d_tau"][i] = abs(math.atan2(math.sin(dt_), math.cos(dt_)))
            m["d_i_i3"][i] = abs(float(d3b[i] - d3a[i]))
            m["rmsd4"][i] = rm["rmsd4"][i]
    return m


def compare(spec: dict, altloc: dict[str, str], work: Path, adapter_python: str, tag: str, extracted: dict | None = None) -> dict:
    a = (extracted or {}).get("A") or extract(spec, "A", altloc["A"], work, adapter_python, tag)
    b = (extracted or {}).get("B") or extract(spec, "B", altloc["B"], work, adapter_python, tag)
    out = work / f"{tag}_qcp.json"
    _run(["-m", "topofold_inputs.qcp_rmsd", str(a["npz"]), str(b["npz"]), str(out)], adapter_python)
    rm = json.loads(out.read_text())
    ok = np.isfinite(a["coords"]).all(1) & np.isfinite(b["coords"]).all(1)
    return {"A": a, "B": b, "qcp": rm, "metrics": window_metrics(a, b, rm),
            "global": {"rmsd_qcp_angstrom": rm["global"], "n_ca": rm["n_global"],
                       "rmsd_numpy_kabsch_angstrom": kabsch_rmsd(a["coords"][ok], b["coords"][ok]) if ok.sum() >= 3 else None}}


# --- within-entry crystal-copy variation (docs/STRUCTPAIR_CHAINS_SPEC_2026-10-25.md) ------------------------------------

def compute_variation(spec: dict, c: dict, work: Path, adapter_python: str, res_rows: list[dict], wrows: list[dict]) -> dict:
    """Every named chain of each entry is read with the existing adapter and compared pairwise: the cross pairings A x B and
    the pairs inside each entry. The primary pair reuses the legacy values."""
    S = spec["structures"]
    chains = {k: S[k]["chains"] for k in "AB"}
    prim = {k: S[k]["chain"] for k in "AB"}
    per = {}
    for k in "AB":
        for ch in chains[k]:
            per[(k, ch)] = c[k] if ch == prim[k] else extract(spec, k, S[k]["altloc"], work, adapter_python, "chain", chain=ch)
    for (k, ch), d in per.items():
        ref = per[(k, prim[k])]
        diff = [i for i, (x, y) in enumerate(zip(ref["residues"], d["residues"]))
                if x["status"] == "present" and y["status"] == "present" and x["resname"] != y["resname"]]
        if diff and spec["resname_mismatch"] == "reject":
            raise SpecError(f"structure {k}: chain {ch} differs from chain {prim[k]} in residue name at "
                            + ", ".join(f"{_rl(ref['residues'][i])} {ref['residues'][i]['resname']} vs {ch} "
                                        f"{d['residues'][i]['resname']}" for i in diff[:5])
                            + "; chains of one entry must agree, or set resname_mismatch to 'flag'")
    keys = cv.pair_keys(chains["A"], chains["B"])
    ent = lambda kind: "A" if kind in ("cross", "withinA") else "B"  # noqa: E731
    sides = lambda key: ((("A", key[1]), ("B", key[2])) if key[0] == "cross"  # noqa: E731
                         else ((ent(key[0]), key[1]), (ent(key[0]), key[2])))
    todo = [k for group in keys.values() for k in group if not (k[0] == "cross" and k[1:] == (prim["A"], prim["B"]))]
    qcp = {("cross", prim["A"], prim["B"]): c["qcp"]}
    if todo:
        pairs = [{"id": "|".join(k), "a": str(per[sides(k)[0]]["npz"]), "b": str(per[sides(k)[1]]["npz"])} for k in todo]
        (work / "chain_pairs.json").write_text(json.dumps({"pairs": pairs}))
        _run(["-m", "topofold_inputs.qcp_rmsd", "--batch", str(work / "chain_pairs.json"), str(work / "chain_pairs_qcp.json")],
             adapter_python)
        got = json.loads((work / "chain_pairs_qcp.json").read_text())["results"]
        qcp.update({k: got["|".join(k)] for k in todo})
    trim = spec["_variation"]["displacement_trim_fraction"]
    pair_values = {}
    for group in keys.values():
        for k in group:
            x, y = per[sides(k)[0]], per[sides(k)[1]]
            m = window_metrics(x, y, qcp[k])
            d = cv.pair_displacement(x["coords"], y["coords"], trim)
            pair_values[k] = {**{ch: m[ch] for ch in METRICS}, "disp": [None if not math.isfinite(v) else float(v) for v in d]}
    labels = [(r["A_residue"], r["B_residue"]) for r in res_rows]
    windows = [(w["window3_A"], w["window4_A"]) for w in wrows]
    rows = cv.variation_rows(pair_values, len(res_rows), labels, windows, chains["A"], chains["B"], (prim["A"], prim["B"]),
                             spec["_variation"])
    block = cv.summarize(rows, chains["A"], chains["B"], (prim["A"], prim["B"]), spec["_variation"])
    inventory = []
    for (k, ch), d in per.items():
        miss = [_rl(r) for r in d["residues"] if r["status"] != "present"]
        inventory.append({"entry": k, "chain": ch, "primary": ch == prim[k], "n_positions": len(d["residues"]),
                          "n_present": len(d["residues"]) - len(miss), "n_missing": len(miss), "missing": ";".join(miss),
                          "n_altloc_resolved": len(d["altloc_decisions"]), "n_breaks": len(d["breaks"])})
    return {"rows": rows, "block": block, "chains_rows": inventory, "n_pairs": {k: len(v) for k, v in keys.items()}}


# --- tables ------------------------------------------------------------------------------------------------------------

def _rl(r: dict) -> str:
    return f"{r['chain']}:{r['resid']}{r.get('icode', '')}"


def _mod_text(r: dict) -> str:
    m = r.get("modification")
    if not m:
        return ""
    return (f"{m['deposited_name']}; parent {m['parent']} (MODRES record, policy table); annotation: {m['annotation'] or 'none'}; record {m['record_type']}")


def residue_rows(c: dict, policy: bool = False) -> tuple[list[dict], list[int]]:
    rows, mismatch = [], []
    alt = {k: {(d["chain"], d["resid"], d.get("icode", "")): d["used"] for d in c[k]["altloc_decisions"]} for k in "AB"}
    for i, (ra, rb) in enumerate(zip(c["A"]["residues"], c["B"]["residues"])):
        both = ra["status"] == "present" and rb["status"] == "present"
        mm = both and ra["resname"] != rb["resname"]
        if mm:
            mismatch.append(i)
        rows.append({"position": i, "A_residue": _rl(ra), "A_resname": ra.get("resname") or "", "A_status": ra["status"],
                     "A_altloc_used": alt["A"].get((ra["chain"], ra["resid"], ra.get("icode", "")), ""),
                     "B_residue": _rl(rb), "B_resname": rb.get("resname") or "", "B_status": rb["status"],
                     "B_altloc_used": alt["B"].get((rb["chain"], rb["resid"], rb.get("icode", "")), ""),
                     "resname_mismatch": mm, **({"A_modification": _mod_text(ra), "B_modification": _mod_text(rb)} if policy else {})})
    return rows, mismatch


def window_rows(c: dict, res_rows: list[dict]) -> list[dict]:
    m, n = c["metrics"], len(res_rows)
    alt = {i for i, r in enumerate(res_rows) if r["A_altloc_used"] or r["B_altloc_used"]}
    mm = {i for i, r in enumerate(res_rows) if r["resname_mismatch"]}
    modified = {i for i, r in enumerate(res_rows) if r.get("A_modification") or r.get("B_modification")}
    with_mod = bool(res_rows) and "A_modification" in res_rows[0]
    rows = []
    for i in range(n):
        w3 = list(range(i - 1, i + 2)) if 1 <= i <= n - 2 else []
        w4 = list(range(i - 1, i + 3)) if 1 <= i <= n - 3 else []
        row = {"position": i, "A_residue": res_rows[i]["A_residue"], "B_residue": res_rows[i]["B_residue"],
               "resname_A": res_rows[i]["A_resname"], "resname_B": res_rows[i]["B_resname"],
               "window3_A": f"{res_rows[w3[0]]['A_residue']}..{res_rows[w3[-1]]['A_residue']}" if w3 else "",
               "window4_A": f"{res_rows[w4[0]]['A_residue']}..{res_rows[w4[-1]]['A_residue']}" if w4 else "",
               "altloc_in_window3": bool(set(w3) & alt), "altloc_in_window4": bool(set(w4) & alt),
               "mismatch_in_window3": bool(set(w3) & mm), "mismatch_in_window4": bool(set(w4) & mm),
               **({"modified_in_window3": bool(set(w3) & modified), "modified_in_window4": bool(set(w4) & modified)} if with_mod else {})}
        for k, unit in COLUMNS:
            row[f"{k}_{unit}"] = m[k][i]
        rows.append(row)
    return rows


def exclusion_rows(c: dict, n: int) -> list[dict]:
    """Windows inside the mapped range with no value, one row per structure and cause (existing records only)."""
    rows = []
    for size, ch, valid_key, lo_hi in ((3, "kappa", "kappa_valid", (1, n - 2)), (4, "tau", "tau_valid", (1, n - 3))):
        for i in range(lo_hi[0], lo_hi[1] + 1):
            if c["A"][valid_key][i] and c["B"][valid_key][i]:
                continue
            for key in "AB":
                s = c[key]
                if s[valid_key][i]:
                    continue
                for cause in exclusion_causes(s, ch, i, np.array([False])):
                    rows.append({"window_size": size, "position": i, "structure": key,
                                 "window_A": f"{_rl(c['A']['residues'][i - 1])}..{_rl(c['A']['residues'][i + size - 2])}",
                                 "cause": cause["cause"], "location": cause["location"], "source": cause["source"],
                                 "detail": cause["detail"]})
    return rows


def summarize(c: dict, wrows: list[dict], top_n: int) -> dict:
    m = c["metrics"]
    top = {}
    for k in METRICS:
        idx = sorted((i for i, v in enumerate(m[k]) if v is not None), key=lambda i: -m[k][i])[:top_n]
        top[k] = [{"position": i, "window_A": wrows[i]["window4_A" if k in FOUR else "window3_A"], "value": m[k][i],
                   "altloc": wrows[i]["altloc_in_window4" if k in FOUR else "altloc_in_window3"],
                   "mismatch": wrows[i]["mismatch_in_window4" if k in FOUR else "mismatch_in_window3"]} for i in idx]
    agree = {}
    for x, y in AGREEMENT_PAIRS:
        p = [(m[x][i], m[y][i]) for i in range(len(m[x])) if m[x][i] is not None and m[y][i] is not None]
        agree[f"{x}~{y}"] = {"spearman": spearman([u for u, _ in p], [v for _, v in p]) if len(p) >= 3 else None,
                             "n_windows": len(p)}
    return {"top": top, "agreement": agree, "n_windows": {k: sum(v is not None for v in m[k]) for k in METRICS}}


def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    return v


def _write_csv(path: Path, rows: list[dict], cols: list[str] | None = None):
    cols = cols or (list(rows[0]) if rows else [])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(v) for k, v in r.items()})


# --- figure and report ---------------------------------------------------------------------------------------------------

def figure(c: dict, wrows: list[dict], res_rows: list[dict], spec: dict, path: Path) -> str:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        from matplotlib.colors import LinearSegmentedColormap
    except ImportError:
        return "unavailable: matplotlib is not installed in this interpreter"
    ink, muted, grid, blue, band, surf = "#0b0b0b", "#52514e", "#e4e3df", "#2a78d6", "#c7dcf5", "#fcfcfb"
    m = c["metrics"]
    A, B = c["A"]["coords"], c["B"]["coords"]
    n = len(A)
    ok = np.isfinite(A).all(1) & np.isfinite(B).all(1)
    cA, cB = A[ok].mean(0), B[ok].mean(0)
    U, _, Vt = np.linalg.svd((A[ok] - cA).T @ (B[ok] - cB))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    Af = (A - cA) @ R.T + cB
    axes2 = np.linalg.svd(B[ok] - cB)[2][:2]
    a2, b2 = (Af - cB) @ axes2.T, (B - cB) @ axes2.T
    brk = {x["after_index"] for k in "AB" for x in c[k]["breaks"]}
    r3 = np.array([np.nan if v is None else v for v in m["rmsd3"]])
    fig = plt.figure(figsize=(14, 8.6), dpi=150, facecolor=surf)
    gs = fig.add_gridspec(3, 3, width_ratios=[2.1, 1, 1], wspace=0.25, hspace=0.55)
    ax0 = fig.add_subplot(gs[:, 0])
    segs, vals = [], []
    for i in range(n - 1):
        if ok[i] and ok[i + 1] and i not in brk:
            ax0.plot(*a2[i: i + 2].T, color="#9d9c96", lw=1.1, zorder=1)
            segs.append([b2[i], b2[i + 1]])
            v = [x for x in (r3[i], r3[i + 1]) if np.isfinite(x)]
            vals.append(float(np.mean(v)) if v else np.nan)
    vals = np.array(vals)
    cmap = LinearSegmentedColormap.from_list("seq", ["#cfe0f6", "#2a78d6", "#0b2f66"])
    cmap.set_bad("#9d9c96")
    lc = LineCollection(segs, cmap=cmap, array=np.ma.masked_invalid(vals), linewidths=3.0, zorder=2, capstyle="round")
    lc.set_clim(0, np.nanmax(vals) if np.isfinite(vals).any() else 1)
    ax0.add_collection(lc)
    step = max(1, n // 8)
    for i in list(range(0, n, step)) + [n - 1]:
        if ok[i]:
            ax0.annotate(res_rows[i]["B_residue"], b2[i], fontsize=7.5, color=ink, xytext=(4, 3), textcoords="offset points")
    ax0.plot([], [], color="#9d9c96", lw=1.1, label=f"A: {spec['structures']['A']['label']} (superposed)")
    ax0.legend(fontsize=7.5, frameon=False, loc="lower left", bbox_to_anchor=(0, 1.0))
    ax0.set_aspect("equal")
    ax0.set_xticks([])
    ax0.set_yticks([])
    for s in ax0.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(lc, ax=ax0, orientation="horizontal", fraction=0.035, pad=0.01, aspect=40)
    cb.set_label(f"B: {spec['structures']['B']['label']}, coloured by 3-residue local Cα RMSD (Å); gaps = breaks or "
                 "missing residues", fontsize=8, color=ink)
    cb.ax.tick_params(labelsize=7, colors=muted)
    cb.outline.set_visible(False)
    ax0.set_title("A  Cα traces after global superposition on all residues present in both\n(projection on the two "
                  "principal axes of B)", fontsize=9, loc="left", color=ink, pad=22)
    titles = {"d_kappa": "|Δκ|, 3-residue (rad)", "d_tau": "|Δτ| circular, 4-residue (rad)",
              "d_i_i2": "|Δd(i−1,i+1)| (Å)", "d_i_i3": "|Δd(i−1,i+2)| (Å)",
              "rmsd3": "local Cα RMSD, 3-residue (Å)", "rmsd4": "local Cα RMSD, 4-residue (Å)"}
    for j, k in enumerate(("d_kappa", "d_tau", "d_i_i2", "d_i_i3", "rmsd3", "rmsd4")):
        ax = fig.add_subplot(gs[j // 2, 1 + j % 2])
        y = np.array([np.nan if v is None else v for v in m[k]])
        ax.plot(np.arange(n), y, color=blue, lw=1.3)
        miss = np.where(~np.isfinite(y))[0]
        lo, hi = (1, n - 2) if k not in FOUR else (1, n - 3)
        for i in miss:
            if lo <= i <= hi:
                ax.axvspan(i - 0.5, i + 0.5, color=grid, lw=0)
        ax.set_title(f"{'BCDEFG'[j]}  {titles[k]}", fontsize=8.5, loc="left", color=ink)
        ax.tick_params(labelsize=7, colors=muted)
        ax.grid(axis="y", color=grid, lw=0.6)
        ax.set_axisbelow(True)
        ax.set_ylim(bottom=0)
        ax.set_xlim(0, n - 1)
        ticks = list(range(0, n, max(1, n // 5)))
        ax.set_xticks(ticks)
        ax.set_xticklabels([res_rows[t]["A_residue"].split(":")[1] for t in ticks])
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        if j >= 4:
            ax.set_xlabel("A residue number (grey: excluded)", fontsize=7.5, color=muted)
    fig.suptitle(f"{spec.get('name') or 'structpair'} · {spec['structures']['A']['label']} vs "
                 f"{spec['structures']['B']['label']} · static, descriptive", fontsize=10, x=0.01, ha="left", color=ink)
    fig.savefig(path, bbox_inches="tight", facecolor=surf)
    plt.close(fig)
    return "written"


def report(spec: dict, summary: dict, manifest: dict, exclusions: list[dict], res_rows: list[dict],
           var_rows: list[dict] | None = None) -> str:
    S = spec["structures"]
    g = summary["global"]
    L = [f"# Two-structure comparison: {spec.get('name') or 'unnamed'}", "",
         f"**Descriptive, static.** {LIMITS}", "",
         "## Inputs", "", "| | Label | File (SHA-256, first 12) | Chain | Altloc policy |", "|---|---|---|---|---|"]
    for k in "AB":
        shown = S[k]["chain"] if len(S[k]["chains"]) == 1 else f"{S[k]['chain']} (primary of {', '.join(S[k]['chains'])})"
        L.append(f"| {k} | {S[k]['label']} | `{Path(S[k]['path']).name}` ({manifest['inputs'][k]['sha256'][:12]}) | "
                 f"{shown} | {S[k]['altloc']} |")
    L += [""] + cv.report_lines(summary["within_entry_variation"], var_rows or [], spec["top_n"])
    n_present = sum(r["A_status"] == "present" and r["B_status"] == "present" for r in res_rows)
    mism = [r for r in res_rows if r["resname_mismatch"]]
    alts = [r for r in res_rows if r["A_altloc_used"] or r["B_altloc_used"]]
    L += ["", "## Mapping", "",
          f"- Correspondence: {manifest['correspondence']['form']} ({len(res_rows)} positions; {n_present} present in "
          f"both). Missing-residue policy: `{spec['on_missing']}`.",
          f"- Residue-name mismatches ({spec['resname_mismatch']}): "
          + (", ".join(f"{r['A_residue']} {r['A_resname']}→{r['B_residue']} {r['B_resname']}" for r in mism) or "none") + ".",
          "- Alternate locations resolved: " + (", ".join(
              f"{r['A_residue'] if r['A_altloc_used'] else r['B_residue']} ({'A' if r['A_altloc_used'] else 'B'}: "
              f"{r['A_altloc_used'] or r['B_altloc_used']})" for r in alts) or "none") + ".",
          f"- Breaks (adapter): A {[b['after_index'] for b in manifest['breaks']['A']] or 'none'}, "
          f"B {[b['after_index'] for b in manifest['breaks']['B']] or 'none'} (row index after which the chain breaks).",
          "", "## Global Cα RMSD (context only)", "",
          f"- {g['rmsd_qcp_angstrom']:.3f} Å over {g['n_ca']} Cα after optimal superposition (MDAnalysis QCP); "
          f"independent numpy SVD/Kabsch check {g['rmsd_numpy_kabsch_angstrom']:.3f} Å "
          f"(|difference| {abs(g['rmsd_qcp_angstrom'] - g['rmsd_numpy_kabsch_angstrom']):.1e} Å).",
          "- Published RMSD values often use other atom sets, residue ranges or outlier rejection; compare only with a "
          "stated definition.", "",
          f"## Largest local differences (top {spec['top_n']}; descriptive)", "",
          "| Metric | Unit | Windows (A numbering; value; * altloc, † residue-name mismatch in window) |", "|---|---|---|"]
    for k in METRICS:
        L.append(f"| {k} | {UNITS[k].replace('angstrom', 'Å')} | " + ", ".join(
            f"{t['window_A']} ({t['value']:.3f}{'*' if t['altloc'] else ''}{'†' if t['mismatch'] else ''})"
            for t in summary["top"][k]) + " |")
    L += ["", "## Agreement between equal-information metrics (Spearman over windows; descriptive)", "",
          "| Pair | Spearman | Windows |", "|---|---|---|"]
    for k, v in summary["agreement"].items():
        rho = "" if v["spearman"] is None else f"{v['spearman']:.3f}"
        L.append(f"| {k} | {rho} | {v['n_windows']} |")
    L += ["", "Three-residue metrics share one triangle, so close agreement is expected; d(i−1, i+2) is blind to the "
          "sign of τ.", "", "## Exclusions", ""]
    if exclusions:
        L += ["| Window (A) | Size | Structure | Cause | Location |", "|---|---|---|---|---|"]
        L += [f"| {e['window_A']} | {e['window_size']} | {e['structure']} | {e['cause']} | {e['location']} |"
              for e in exclusions[:60]]
        if len(exclusions) > 60:
            L.append(f"| … {len(exclusions) - 60} more rows in `exclusions.csv` | | | | |")
    else:
        L.append("No window inside the mapped range was excluded.")
    mod = summary.get("modified_residues")
    if mod:
        L += ["", "## Modified residues", "", f"Policy `{mod['policy']}`: {mod['scope']}. {mod['statement']}", ""]
        if mod["retained"]:
            L += ["| Structure | Residue | Deposited name | Parent (source) | Annotation | Record |", "|---|---|---|---|---|---|"]
            L += [f"| {r['structure']} | {r['residue']} | {r['deposited_name']} | {r['parent']} ({r['parent_source']}) | {r['annotation'] or 'none'} | {r['record_type']} |" for r in mod["retained"]]
        else:
            L.append("No modified residue was retained.")
        if mod["excluded"]:
            L += ["", "Not admitted (in scope):", "", "| Structure | Residue | Name | Reason |", "|---|---|---|---|"]
            L += [f"| {e['structure']} | {e['residue']} | {e['name']} | {e['reason']} |" for e in mod["excluded"]]
    sens = summary.get("altloc_sensitivity")
    L += ["", "## Alternate-conformation sensitivity", ""]
    if sens is None:
        L.append("Not requested.")
    elif sens.get("available") is False:
        L.append(f"**UNAVAILABLE, not zero.** The requested sensitivity (policy `{sens['policy']}`) was not computed: {sens['reason']} "
                 f"{sens['why_not_partial']} No value is reported.")
        L += ["", "| Structure | Residue | In the correspondence | Alternate locations | Occupancies |", "|---|---|---|---|---|"]
        for t in sens["ties"]:
            L.append(f"| {t['structure']} | {t['chain']}:{t['resid']}{t['icode']} | {'yes' if t['in_correspondence'] else 'no'} | "
                     f"{', '.join(a or '(blank)' for a in t['altlocs'])} | {', '.join(str(o) for o in t['occupancies'])} |")
    elif "max_abs_change" not in sens:
        L.append(sens["status"])
    else:
        L.append(f"Policy `{sens['policy']}`: largest change per metric "
                 + ", ".join(f"{k} " + ("n/a" if v is None else f"{v:.4f}") for k, v in sens["max_abs_change"].items())
                 + f"; global RMSD {sens['global_rmsd_angstrom']:.3f} Å.")
    L += ["", "## Figure", "", "![structure pair](figure.png)" if manifest["figure"] == "written" else manifest["figure"],
          "", f"Files: `windows.csv` (all metrics with units in column names), `residues.csv`, `exclusions.csv`, "
          f"`summary.json`, `manifest.json` (hashes and software)"
          + (f", `{nr.FILE_NAME}` (computed values before display rounding)" if manifest["policies"].get("numeric_record") else "")
          + f". structpair {STRUCTPAIR_VERSION}.", ""]
    return "\n".join(L)


# --- main ----------------------------------------------------------------------------------------------------------------

def modified_block(policy: str, c: dict, res_rows: list[dict]) -> dict:
    """Retained modified residues (deposited name, parent and its source, annotation) and the in-scope cases the policy did not admit."""
    retained, excluded = [], []
    for k in "AB":
        for r, row in zip(c[k]["residues"], res_rows):
            m = r.get("modification")
            if m:
                retained.append({"structure": k, "residue": _rl(r), "deposited_name": m["deposited_name"], "parent": m["parent"], "parent_source": m["parent_source"],
                                 "annotation": m["annotation"], "record_type": m["record_type"]})
        for e in (c[k].get("modified_block") or {}).get("excluded", []):
            excluded.append({"structure": k, "residue": f"{e['chain']}:{e['resid']}{e['icode']}", "name": e["name"], "record_type": e["record_type"], "reason": e["reason"]})
    return {"policy": policy, "scope": "SEP, TPO and PTR in PDB-format files with MODRES records; no other modified residue or format",
            "statement": "The deposited residue name is kept; the parent is documentation for the correspondence and is not a rename (SEP versus SER stays a visible name difference). "
                         "A modification annotation is not evidence of functional state or of any effect.",
            "retained": retained, "excluded": excluded}


def unavailable_sensitivity(policy: str, ties: list) -> dict:
    """The requested optional sensitivity could not be computed. No number is reported (in particular no zero)."""
    return {"policy": policy, "status": "unavailable: alternate locations tied on occupancy", "available": False, "values_reported": False,
            "reason": f"the `{policy}` policy cannot choose between alternate locations of equal occupancy at {len(ties)} residue(s); "
                      "the primary comparison (primary altloc policy) is complete and unaffected.",
            "why_not_partial": "A partial sensitivity would need a new rule for the tied residues (choose one or skip it), which changes the windows "
                               "around them and the global fit, so the whole optional sensitivity is unavailable. This does not mean the "
                               "structures are insensitive to the alternate-location choice.",
            "ties": ties}


def run(spec_path: Path, out: Path, adapter_python: str) -> dict:
    spec = load_spec(spec_path)
    if out.exists() and any(out.iterdir()):
        raise SpecError(f"output directory {out} is not empty; use a fresh directory")
    work = out / "work"
    work.mkdir(parents=True, exist_ok=True)
    c = compare(spec, {k: spec["structures"][k]["altloc"] for k in "AB"}, work, adapter_python, "primary")
    res_rows, mismatch = residue_rows(c, spec["modified_residues"] != "none")
    if mismatch and spec["resname_mismatch"] == "reject":
        raise SpecError("residue names differ at " + ", ".join(
            f"{res_rows[i]['A_residue']} {res_rows[i]['A_resname']} vs {res_rows[i]['B_residue']} {res_rows[i]['B_resname']}"
            for i in mismatch) + "; fix the correspondence or set resname_mismatch to 'flag'")
    wrows = window_rows(c, res_rows)
    excl = exclusion_rows(c, len(res_rows))
    summary = summarize(c, wrows, spec["top_n"])
    summary["global"] = c["global"]
    if spec["modified_residues"] != "none":
        summary["modified_residues"] = modified_block(spec["modified_residues"], c, res_rows)
    summary["altloc_sensitivity"] = None
    if spec["altloc_sensitivity"]:
        pol = spec["altloc_sensitivity"]
        if not (c["A"]["altloc_decisions"] or c["B"]["altloc_decisions"]):
            summary["altloc_sensitivity"] = {"policy": pol, "status": "not applicable: no alternate locations at the "
                                                                      "mapped Cα atoms"}
        else:
            ex, ties = {}, []
            for k in "AB":          # both structures are examined, so ties in A and in B are all reported
                try:
                    ex[k] = extract(spec, k, pol, work, adapter_python, "altloc_sens")
                except AltlocTieRefusal as e:
                    ties += [{"structure": k, **t} for t in e.ties]
            if ties:
                summary["altloc_sensitivity"] = unavailable_sensitivity(pol, ties)
        if summary["altloc_sensitivity"] is None:
            c2 = compare(spec, {"A": pol, "B": pol}, work, adapter_python, "altloc_sens", extracted=ex)
            summary["altloc_sensitivity"] = {
                "policy": pol, "status": "computed",
                "max_abs_change": {k: max((abs(c2["metrics"][k][i] - c["metrics"][k][i]) for i in range(len(res_rows))
                                           if c["metrics"][k][i] is not None and c2["metrics"][k][i] is not None),
                                          default=None) for k in METRICS},
                "global_rmsd_angstrom": c2["global"]["rmsd_qcp_angstrom"]}
    _write_csv(out / "windows.csv", wrows)
    _write_csv(out / "residues.csv", res_rows)
    _write_csv(out / "exclusions.csv", excl, ["window_size", "position", "structure", "window_A", "cause", "location",
                                              "source", "detail"])
    fig_status = figure(c, wrows, res_rows, spec, out / "figure.png")
    summary["metric_definitions"] = DESCRIPTION
    summary["units"] = {k: UNITS[k] for k in METRICS}
    summary["n_exclusion_rows"] = len(excl)
    numeric = None
    if spec["numeric_record"]:
        numeric = nr.build(c, res_rows, excl, summary, UNITS, {"structpair_version": STRUCTPAIR_VERSION},
                           {k: _sha256(spec["structures"][k]["path"]) for k in "AB"}, _sha256(spec_path))
        (out / nr.FILE_NAME).write_text(nr.dumps(numeric))
    var = None
    var_fig = None
    if any(len(spec["structures"][k]["chains"]) > 1 for k in "AB"):
        var = compute_variation(spec, c, work, adapter_python, res_rows, wrows)
        summary["within_entry_variation"] = var["block"]
        _write_csv(out / "chain_variation.csv", var["rows"])
        _write_csv(out / "chains.csv", var["chains_rows"])
        var_fig = cv.figure(var["rows"], out / "chain_variation.png")
    else:
        summary["within_entry_variation"] = cv.unavailable_block()
    (out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    outputs = (["windows.csv", "residues.csv", "exclusions.csv", "summary.json"] + (["figure.png"] if fig_status == "written" else [])
               + (["chain_variation.csv", "chains.csv"] if var else []) + (["chain_variation.png"] if var_fig == "written" else [])
               + ([nr.FILE_NAME] if numeric else []))
    manifest = {
        "schema": SCHEMA, "name": spec.get("name"), "statement": LIMITS,
        "command": " ".join([Path(sys.executable).name, "-m", "topofold_inputs.structpair", *sys.argv[1:]]),
        "spec": {"path": str(spec_path), "sha256": _sha256(spec_path)},
        "inputs": {k: {"label": spec["structures"][k]["label"], "path": spec["structures"][k]["path"],
                       "sha256": _sha256(spec["structures"][k]["path"]), "chain": spec["structures"][k]["chain"],
                       "chains": spec["structures"][k]["chains"], "altloc_policy": spec["structures"][k]["altloc"],
                       **c[k]["provenance"]} for k in "AB"},
        "correspondence": {"form": next(iter(spec["correspondence"])), "n_positions": len(res_rows),
                           "pairs": [[f"{a[0]}{a[1]}", f"{b[0]}{b[1]}"] for a, b in spec["_pairs"]]},
        "policies": {"on_missing": spec["on_missing"], "resname_mismatch": spec["resname_mismatch"],
                     "altloc_sensitivity": spec["altloc_sensitivity"],
                     "numeric_record": bool(numeric),
                     "modified_residues": spec["modified_residues"],
                     **({"modified_residues_outcome": f"{len(summary['modified_residues']['retained'])} retained, "
                                                      f"{len(summary['modified_residues']['excluded'])} excluded"}
                        if spec["modified_residues"] != "none" else {}),
                     "altloc_sensitivity_outcome": ("not requested" if summary["altloc_sensitivity"] is None
                                                    else summary["altloc_sensitivity"]["status"]),
                     "chain_variation": spec["_variation"]},
        "within_entry_variation": {"status": summary["within_entry_variation"]["status"], "advisory_semantics": cv.SEMANTICS,
                                   "statement": cv.STATEMENT, "figure": var_fig, "n_pairs": var["n_pairs"] if var else None},
        "breaks": {k: c[k]["breaks"] for k in "AB"},
        "software": {"structpair_version": STRUCTPAIR_VERSION, "module_sha256": _sha256(__file__),
                     "qcp_module_sha256": _sha256(Path(__file__).with_name("qcp_rmsd.py")),
                     "python": platform.python_version(), "numpy": np.__version__, "qcp": c["qcp"]["software"],
                     "adapter_interpreter": adapter_python},
        "figure": fig_status,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    manifest["outputs_sha256"] = {p: _sha256(out / p) for p in outputs}
    (out / "report.md").write_text(report(spec, summary, manifest, excl, res_rows, var["rows"] if var else None))
    manifest["outputs_sha256"]["report.md"] = _sha256(out / "report.md")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    return {"summary": summary, "manifest": manifest, "windows": wrows, "residues": res_rows, "exclusions": excl,
            "variation": var}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m topofold_inputs.structpair", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spec")
    ap.add_argument("--out", required=True, help="new or empty output directory")
    ap.add_argument("--adapter-python", default=os.environ.get("TOPOFOLD_ADAPTER_PYTHON"),
                    help="interpreter with MDAnalysis (default: this one, if it has MDAnalysis)")
    a = ap.parse_args(argv)
    adapter_python = a.adapter_python or (sys.executable if importlib.util.find_spec("MDAnalysis") else None)
    if adapter_python is None:
        print("no interpreter with MDAnalysis: pass --adapter-python or set TOPOFOLD_ADAPTER_PYTHON", file=sys.stderr)
        return 2
    try:
        r = run(Path(a.spec), Path(a.out), adapter_python)
    except SpecError as e:
        print(f"spec error: {e}", file=sys.stderr)
        return 2
    except RuntimeError as e:      # adapter or QCP subprocess refused the input (message from the adapter)
        print(f"input error: {str(e).strip().splitlines()[-1]}", file=sys.stderr)
        return 2
    g = r["summary"]["global"]
    print(f"{len(r['residues'])} positions; global Cα RMSD {g['rmsd_qcp_angstrom']:.3f} Å over {g['n_ca']} Cα; "
          f"{len(r['exclusions'])} exclusion rows; within-entry variation: {r['summary']['within_entry_variation']['status']}; "
          f"wrote {a.out}/report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
