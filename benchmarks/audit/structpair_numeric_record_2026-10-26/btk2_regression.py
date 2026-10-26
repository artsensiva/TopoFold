"""Separate regression check of numeric_record.json on the original and rigidly transformed BTK-2 inputs (m13b).

    .venv/bin/python benchmarks/audit/structpair_numeric_record_2026-10-26/btk2_regression.py

Not a new evaluation and not a reclassification of KIN-PAIR-1. Uses the KIN-PAIR-1 supplement specifications of BTK-2 (`work_supplement/BTK-2/`,
read only) with `"numeric_record": true` added, runs structpair 0.4.0 on both, and checks: (1) the default outputs equal the earlier runs
byte for byte; (2) the record equals the pre-format computed values exactly; (3) what remains of the original-vs-transformed difference
(float32 coordinates in the adapter). Writes results/btk2_numeric_record_regression.json.
"""

import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
KP = ROOT / "benchmarks/audit/kin_pair_1"
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs import numeric_record as nr  # noqa: E402
from topofold_inputs import structpair as sp  # noqa: E402
from topofold_inputs import structpair_verify as V  # noqa: E402

ADAPTER = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
WS = KP / "work_supplement/BTK-2"
WORK = HERE / "work"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def run(tag, spec_src):
    d = WORK / tag
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    s = json.loads(spec_src.read_text())
    s["numeric_record"] = True
    s.pop("altloc_sensitivity", None)
    spec = d / "spec.json"
    spec.write_text(json.dumps(s))
    r = subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(spec), "--out", str(d / "out")], capture_output=True, text=True,
                       env=dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=ADAPTER))
    assert r.returncode == 0, r.stderr
    return d, spec


def main():
    out = {"id": "STRUCTPAIR-NUMERIC-RECORD-BTK2-REGRESSION-2026-10-26", "note": "engineering regression; KIN-PAIR-1 is not reclassified"}
    runs = {}
    for tag, src, ref in (("base", WS / "spec.json", KP / "results/pairs_supplement/BTK-2"), ("moved", WS / "spec_moved.json", WS / "out_moved")):
        d, spec = run(tag, src)
        rec = nr.read(d / "out")
        spec_obj = sp.load_spec(spec)
        w = d / "w"
        w.mkdir()
        c = sp.compare(spec_obj, {k: spec_obj["structures"][k]["altloc"] for k in "AB"}, w, ADAPTER, tag)
        exact = all(rec["metrics"][k]["values"] == [None if v is None else float(v) for v in vals] for k, vals in c["metrics"].items())
        csv_rows = list(csv.DictReader(open(d / "out/windows.csv")))
        runs[tag] = {"dir": d / "out", "rec": rec, "pre": c["metrics"], "csv": csv_rows}
        out[tag] = {"verification_all_ok": V.verify(d / "out")["all_ok"], "record_equals_preformat_values_exactly": exact,
                    "n_values": sum(v is not None for m in rec["metrics"].values() for v in m["values"]),
                    "default_outputs_byte_identical_to_earlier_run": {f: sha(d / "out" / f) == sha(ref / f) for f in ("windows.csv", "residues.csv", "exclusions.csv")}}
    b, m = runs["base"], runs["moved"]
    COLS = {k: f"{k}_{u}" for k, u in sp.COLUMNS}
    worst_rec = (0.0, None)
    worst_csv = (0.0, None)
    for k in b["rec"]["metrics"]:
        for j, (x, y) in enumerate(zip(b["rec"]["metrics"][k]["values"], m["rec"]["metrics"][k]["values"])):
            if x is None or y is None:
                continue
            if abs(x - y) > worst_rec[0]:
                worst_rec = (abs(x - y), (k, j))
            if k in COLS:
                cx, cy = float(b["csv"][j][COLS[k]]), float(m["csv"][j][COLS[k]])
                if abs(cx - cy) > worst_csv[0]:
                    worst_csv = (abs(cx - cy), (k, j))

    def cell(k, j):
        return {"metric": k, "position": j, "A_residue": b["rec"]["positions"][j]["A_residue"], "B_residue": b["rec"]["positions"][j]["B_residue"],
                "unit": b["rec"]["metrics"][k]["unit"],
                "base": {"computed": b["pre"][k][j], "csv": b["csv"][j][COLS[k]], "record": b["rec"]["metrics"][k]["values"][j]},
                "moved": {"computed": m["pre"][k][j], "csv": m["csv"][j][COLS[k]], "record": m["rec"]["metrics"][k]["values"][j]}}
    cells = [cell(*worst_csv[1]), cell(*worst_rec[1])]
    for cl in cells:
        cl["difference_record"] = abs(cl["base"]["record"] - cl["moved"]["record"])
        cl["difference_csv"] = abs(float(cl["base"]["csv"]) - float(cl["moved"]["csv"]))
    diag = json.loads((HERE.parent / "structpair_altloc_tie_repair_2026-10-26/results/btk2_precision_diagnosis.json").read_text())
    out.update({"worst_original_vs_transformed": {"in_the_record": {"difference": worst_rec[0], "metric": worst_rec[1][0], "position": worst_rec[1][1]},
                                                  "in_the_csv": {"difference": worst_csv[0], "metric": worst_csv[1][0], "position": worst_csv[1][1]}},
                "table_cells": cells,
                "independent_float64_from_text_max_deviation_earlier_diagnosis": diag["btk2"]["independent_float64_from_text_max_deviation"],
                "float32_spacing_at_largest_transformed_coordinate": diag["btk2"]["float32_coordinate_storage"]["moved"]["float32_spacing_at_that_magnitude_angstrom"],
                "interpretation": "the record preserves the computed values exactly; the remaining original-vs-transformed difference is the float32 coordinate step of the adapter, which the record does not remove"})
    (HERE / "results/btk2_numeric_record_regression.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
    print(json.dumps({k: out[k] for k in ("base", "moved", "worst_original_vs_transformed")}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
