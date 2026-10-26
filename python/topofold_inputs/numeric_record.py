"""Opt-in numerical record of a structpair result: the computed window values before display rounding.

The CSV files round values to six significant figures for display (`windows.csv`). With `"numeric_record": true` in the spec,
structpair also writes `numeric_record.json`, listed with its SHA-256 in `manifest.json` and checked by
`python -m topofold_inputs.structpair_verify`. It preserves the computed float64 values exactly (JSON floats round-trip), per position and
metric, with stable identifiers, units and the status of every missing value.

What it is: computational precision of the values structpair computed. What it is not: higher experimental accuracy, and it does not
recover precision lost earlier: the Cα adapter holds coordinates as float32 (about 1e-5 Å at 150 Å), and the values here inherit that.
Geometry, coordinates, tolerances and the default outputs are not changed.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

NUMERIC_SCHEMA = "topofold_inputs.structpair.numeric/1"
FILE_NAME = "numeric_record.json"
STATEMENT = ("Computed float64 values before display rounding (the CSV files round to six significant figures). This preserves "
             "computational precision only: it is not higher experimental accuracy and does not recover precision lost in the float32 "
             "coordinates held by the Cα adapter. A missing value is null with its status; it is never zero.")
THREE = ("d_kappa", "d_i_i2", "rmsd3", "kappa_A", "kappa_B", "d2_A", "d2_B")      # defined on three-residue windows; the rest on four


class NumericRecordError(ValueError):
    """The record is malformed, inconsistent with its result, or does not match the manifest."""


def build(c: dict, res_rows: list[dict], excl: list[dict], summary: dict, units: dict, software: dict, input_sha: dict, spec_sha: str) -> dict:
    n = len(res_rows)
    excluded: dict = {}
    for r in excl:
        excluded.setdefault((int(r["window_size"]), int(r["position"])), []).append(f"{r['structure']}:{r['cause']}")
    metrics = {}
    for k, vals in c["metrics"].items():
        size = 3 if k in THREE else 4
        status = []
        for i, v in enumerate(vals):
            if v is not None:
                status.append("value")
            elif (size, i) in excluded:
                status.append("excluded: " + ", ".join(sorted(set(excluded[(size, i)]))))
            else:
                status.append("not defined: outside the window range of the chain ends")
        metrics[k] = {"unit": units.get(k) or units[k.rsplit("_", 1)[0]], "window_residues": size, "values": [None if v is None else float(v) for v in vals], "status": status}
    g = c["global"]
    return {"schema": NUMERIC_SCHEMA, "statement": STATEMENT, "n_positions": n,
            "positions": [{"position": i, "A_residue": r["A_residue"], "B_residue": r["B_residue"], "A_resname": r["A_resname"], "B_resname": r["B_resname"]}
                          for i, r in enumerate(res_rows)],
            "metrics": metrics,
            "global": {"rmsd_qcp_angstrom": float(g["rmsd_qcp_angstrom"]), "n_ca": int(g["n_ca"]),
                       "rmsd_numpy_kabsch_angstrom": None if g.get("rmsd_numpy_kabsch_angstrom") is None else float(g["rmsd_numpy_kabsch_angstrom"])},
            "altloc_sensitivity": summary.get("altloc_sensitivity"),   # as in summary.json: numbers if computed, else the unavailable / not applicable record (never a zero)
            "provenance": {"structpair_version": software["structpair_version"],
                           "coordinate_precision": "float32 in the Cα adapter; values computed in float64 from them",
                           "rounded_display_files": "windows.csv, residues.csv (six significant figures)",
                           "input_sha256": input_sha, "spec_sha256": spec_sha}}


def dumps(rec: dict) -> str:
    """One line per position and per metric; JSON floats are written with repr, so reading returns the identical float64."""
    lines = ["{", f'"schema": {json.dumps(rec["schema"])},']
    for k, v in rec.items():
        if k not in ("schema", "positions", "metrics"):
            lines.append(f"{json.dumps(k)}: {json.dumps(v)},")
    lines += ['"positions": [', ",\n".join(json.dumps(p) for p in rec["positions"]), "],",
              '"metrics": {', ",\n".join(f"{json.dumps(k)}: {json.dumps(m)}" for k, m in rec["metrics"].items()), "}", "}"]
    return "\n".join(lines) + "\n"


def read(result_dir) -> dict:
    """Read and check the record of a result directory: manifest hash, schema, shape, status consistency and identifier alignment."""
    d = Path(result_dir)
    manifest = json.loads((d / "manifest.json").read_text())
    if not (manifest.get("policies") or {}).get("numeric_record"):
        raise NumericRecordError("this result was written without a numeric record (spec numeric_record is not true)")
    want = manifest["outputs_sha256"].get(FILE_NAME)
    f = d / FILE_NAME
    if want is None:
        raise NumericRecordError(f"manifest declares a numeric record but lists no hash for {FILE_NAME}")
    if not f.is_file():
        raise NumericRecordError(f"{FILE_NAME} is declared in the manifest but missing")
    if hashlib.sha256(f.read_bytes()).hexdigest() != want:
        raise NumericRecordError(f"{FILE_NAME} differs from its manifest hash")
    rec = json.loads(f.read_text())
    if rec.get("schema") != NUMERIC_SCHEMA:
        raise NumericRecordError(f"unknown numeric record schema {rec.get('schema')!r}")
    n = rec["n_positions"]
    if len(rec["positions"]) != n or any(len(m["values"]) != n or len(m["status"]) != n for m in rec["metrics"].values()):
        raise NumericRecordError("positions and metric arrays differ in length")
    for m in rec["metrics"].values():
        for v, s in zip(m["values"], m["status"]):
            if (v is None) == (s == "value"):
                raise NumericRecordError("a value and its status disagree (a missing value must carry a reason; a value has status 'value')")
    res = d / "residues.csv"
    if res.is_file():
        rows = list(csv.DictReader(open(res)))
        if len(rows) != n or any(r["position"] != str(p["position"]) or r["A_residue"] != p["A_residue"] or r["B_residue"] != p["B_residue"]
                                 for r, p in zip(rows, rec["positions"])):
            raise NumericRecordError("positions or residue identifiers do not align with residues.csv")
    return rec
