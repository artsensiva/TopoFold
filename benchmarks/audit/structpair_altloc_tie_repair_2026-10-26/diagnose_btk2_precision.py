"""Diagnosis of the KIN-PAIR-1 supplement BTK-2 invariance deviation of 1.0000000000015433e-4 (tolerance 1e-4).

    .venv/bin/python benchmarks/audit/structpair_altloc_tie_repair_2026-10-26/diagnose_btk2_precision.py

Traces the value through (1) the serialized CSV of the base and moved runs, (2) the pre-format floating-point values of the same
two comparisons computed by structpair's own functions, (3) an independent float64 computation from the PDB text, (4) the decimal
arithmetic of the two serialized strings and the floating-point subtraction at the threshold, and (5) a SYNTHETIC rigid-transform
control. Reads the KIN-PAIR-1 supplement files (read only); writes results/btk2_precision_diagnosis.json. No tolerance is changed.
"""

import csv
import decimal
import json
import math
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
KP = ROOT / "benchmarks/audit/kin_pair_1"
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(KP))
from topofold_inputs import structpair as sp  # noqa: E402
import kinpair as kp  # noqa: E402

ADAPTER = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
WS = KP / "work_supplement/BTK-2"
TOL = 1e-4


def numeric_cols(row):
    return [c for c in row if c.endswith(("_rad", "_angstrom"))]


def serialized():
    base = list(csv.DictReader(open(WS / "out/windows.csv")))
    moved = list(csv.DictReader(open(WS / "out_moved/windows.csv")))
    best = (0.0, None)
    n_at_quantum = n_over_tol = n_cmp = 0
    for j, (a, b) in enumerate(zip(base, moved)):
        for c in numeric_cols(a):
            if a[c] == "" or b[c] == "":
                continue
            d = abs(float(a[c]) - float(b[c]))
            n_cmp += 1
            n_over_tol += d > TOL
            n_at_quantum += d >= 0.99e-4
            if d > best[0]:
                best = (d, (j, c, a[c], b[c]))
    j, c, sa, sb = best[1]
    dec = abs(decimal.Decimal(sa) - decimal.Decimal(sb))
    return {"n_values_compared": n_cmp, "max_serialized_deviation": best[0], "location": {"position": j, "column": c, "base_text": sa, "moved_text": sb},
            "exact_decimal_difference_of_the_two_strings": str(dec), "float_difference_repr": repr(abs(float(sa) - float(sb))),
            "float_difference_minus_tolerance": abs(float(sa) - float(sb)) - TOL, "n_values_over_tolerance": n_over_tol,
            "n_values_at_one_csv_quantum_or_more": n_at_quantum, "csv_quantum_at_this_magnitude_6_significant_figures": 10.0 ** (math.floor(math.log10(abs(float(sa)))) - 5)}


def preformat(spec_path, tag, work):
    spec = sp.load_spec(spec_path)
    c = sp.compare(spec, {k: spec["structures"][k]["altloc"] for k in "AB"}, work, ADAPTER, tag)
    return spec, c


def metric_col(k):
    return f"{k}_{sp.UNITS[k]}"


def float32_storage(c):
    X = np.concatenate([c["A"]["coords"], c["B"]["coords"]])
    X = X[np.isfinite(X).all(axis=1)]
    return {"max_abs_coordinate_angstrom": float(np.abs(X).max()), "float32_spacing_at_that_magnitude_angstrom": float(np.spacing(np.float32(np.abs(X).max())))}


def diagnose_btk2():
    s = serialized()
    loc = s["location"]
    tmp = Path(tempfile.mkdtemp(prefix="btk2_diag_", dir=HERE / "work")) if (HERE / "work").exists() else Path(tempfile.mkdtemp(prefix="btk2_diag_"))
    (tmp / "b").mkdir()
    (tmp / "m").mkdir()
    spec_b, cb = preformat(WS / "spec.json", "base", tmp / "b")
    spec_m, cm = preformat(WS / "spec_moved.json", "moved", tmp / "m")
    key = loc["column"].rsplit("_", 1)[0]
    key = next(k for k in sp.METRICS if loc["column"] == metric_col(k)) if any(loc["column"] == metric_col(k) for k in sp.METRICS) else key
    j = loc["position"]
    mb, mm = cb["metrics"], cm["metrics"]
    pre_dev = {}
    for k in mb:
        d = [abs(x - y) for x, y in zip(mb[k], mm[k]) if x is not None and y is not None]
        if d:
            pre_dev[k] = max(d)
    all_pre = max(pre_dev.values())
    at = {"metric": key, "position": j, "base_preformat": mb[key][j], "moved_preformat": mm[key][j], "preformat_difference": abs(mb[key][j] - mm[key][j]),
          "base_formatted": f"{mb[key][j]:.6g}", "moved_formatted": f"{mm[key][j]:.6g}"}
    # independent float64 from the PDB text of the same files
    files = {k: Path(spec_b["structures"][k]["path"]) for k in "AB"}
    mfiles = {k: Path(spec_m["structures"][k]["path"]) for k in "AB"}
    def windows(paths, spec):
        P = {k: kp.parse_pdb(paths[k].read_text(encoding="latin-1")) for k in "AB"}
        ra = {r.key: r for r in P["A"].residues[spec["structures"]["A"]["chain"]]}
        rb = {r.key: r for r in P["B"].residues[spec["structures"]["B"]["chain"]]}
        return kp.independent_windows([ra[a] for a, _ in [(p[0], 0) for p in spec["_pairs"]]], [rb[b] for _, b in [(0, p[1]) for p in spec["_pairs"]]])
    w0, w1 = windows(files, spec_b), windows(mfiles, spec_m)
    dev64 = max(abs(x - y) for k in ("d_kappa", "d_i_i2", "rmsd3", "d_tau", "d_i_i3", "rmsd4") for x, y in zip(w0[k], w1[k]) if x is not None and y is not None)
    return {"serialized": s, "preformat_at_the_same_value": at, "preformat_max_deviation_any_metric": all_pre, "preformat_max_deviation_by_metric": pre_dev,
            "float32_coordinate_storage": {"base": float32_storage(cb), "moved": float32_storage(cm)},
            "independent_float64_from_text_max_deviation": dev64}


def synthetic_control():
    """A SYNTHETIC helix pair, moved by exact coordinate permutations/translations: deviations before and after formatting."""
    rng = np.random.default_rng(7)
    N = 80
    t = np.arange(N) * math.radians(100)
    X = np.c_[2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(N)] + rng.normal(0, 0.05, (N, 3))
    Y = X + rng.normal(0, 0.4, X.shape)
    tmp = Path(tempfile.mkdtemp(prefix="syn_ctrl_"))

    def write(path, C, chain="A"):
        lines = [f"ATOM  {i + 1:5d}  CA  ALA {chain}{i + 1:4d}    {c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}  1.00 20.00           C" for i, c in enumerate(C)]
        Path(path).write_text("\n".join(lines) + "\nEND\n")
    write(tmp / "a.pdb", X)
    write(tmp / "b.pdb", Y)
    RA = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], float)
    RB = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], float)
    write(tmp / "am.pdb", np.round(X, 3) @ RA.T + np.array([12.5, -7.25, 3.0]))
    write(tmp / "bm.pdb", np.round(Y, 3) @ RB.T + np.array([-20.0, 5.0, 40.0]))

    def spec(a, b, name):
        s = {"name": "SYNTHETIC rigid control", "correspondence": {"same_numbering": f"1-{N}"}, "on_missing": "exclude",
             "structures": {"A": {"label": "SYNTHETIC A", "file": str(tmp / a), "chain": "A", "altloc": "reject"},
                            "B": {"label": "SYNTHETIC B", "file": str(tmp / b), "chain": "A", "altloc": "reject"}}}
        (tmp / name).write_text(json.dumps(s))
        return sp.load_spec(tmp / name)
    (tmp / "w0").mkdir()
    (tmp / "w1").mkdir()
    s0, s1 = spec("a.pdb", "b.pdb", "s0.json"), spec("am.pdb", "bm.pdb", "s1.json")
    c0 = sp.compare(s0, {"A": "reject", "B": "reject"}, tmp / "w0", ADAPTER, "base")
    c1 = sp.compare(s1, {"A": "reject", "B": "reject"}, tmp / "w1", ADAPTER, "moved")
    pre, ser, worst = 0.0, 0.0, None
    for k in c0["metrics"]:
        for x, y in zip(c0["metrics"][k], c1["metrics"][k]):
            if x is None or y is None:
                continue
            pre = max(pre, abs(x - y))
            d = abs(float(f"{x:.6g}") - float(f"{y:.6g}"))
            if d > ser:
                ser, worst = d, (k, f"{x:.6g}", f"{y:.6g}", abs(x - y))
    return {"max_preformat_deviation": pre, "max_serialized_deviation": ser, "location": worst,
            "serialized_exceeds_1e-4": ser > TOL, "serialized_minus_tolerance": ser - TOL}


def main():
    out = {"id": "BTK2-PRECISION-DIAGNOSIS-2026-10-26", "tolerance": TOL, "target": "KIN-PAIR-1 supplement BTK-2, R7 max deviation 0.0001000000000015433",
           "btk2": diagnose_btk2(), "synthetic_control": synthetic_control()}
    (HERE / "results/btk2_precision_diagnosis.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
    print(json.dumps({"serialized": out["btk2"]["serialized"]["max_serialized_deviation"], "preformat_same_value": out["btk2"]["preformat_at_the_same_value"]["preformat_difference"],
                      "preformat_max": out["btk2"]["preformat_max_deviation_any_metric"], "float64_text": out["btk2"]["independent_float64_from_text_max_deviation"],
                      "synthetic": out["synthetic_control"]}, default=str, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
