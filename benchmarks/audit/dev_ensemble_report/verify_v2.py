"""Verify the corrected (v2) demonstrations against v1 and the stated properties; writes results_v2/verification.json.

    .venv/bin/python benchmarks/audit/dev_ensemble_report/verify_v2.py
Reads git-ignored build/v2_<name>_a outputs (regenerate with the command in README.md) and the committed v1 tables.
"""

import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "python"))
import topofold  # noqa: E402
from topofold_inputs import ensemble as en  # noqa: E402
from topofold_inputs.features import core_provenance  # noqa: E402

SHIFTS = [0.3, 1.7, -2.9, 3.141592653589793, 5.0]


def table(path):
    return list(csv.DictReader(open(path)))


def main():
    out = {"note": "verification of corrected v2 demonstrations; descriptive tool, no inference", "datasets": {}}
    native = core_provenance(topofold)["native_extension"]
    for name in ("1ha8A00", "adk"):
        v1, v2 = HERE / "results" / name, HERE / "build" / f"v2_{name}_a"
        b1 = {(r["channel"], r["window"], r["group_a"], r["group_b"]): r["pooled_js"] for r in table(v1 / "tables/between.csv")}
        b2 = {(r["channel"], r["window"], r["group_a"], r["group_b"]): r["pooled_js"] for r in table(v2 / "tables/between.csv")}
        same_keys = b1.keys() == b2.keys()
        max_js = max(abs(float(b1[k]) - float(b2[k])) for k in b1) if same_keys else None
        # origin invariance of the circular diagnostic on the real members' τ series
        t2 = {(r["member"], r["window"]): r for r in table(v2 / "tables/temporal.csv") if r["channel"] == "tau"}
        inv_max, recomputed_max, n_series = 0.0, 0.0, 0
        for fz in sorted((v2 / "members").glob("*_features.npz")):
            member = fz.name[: -len("_features.npz")]
            f = np.load(fz)
            tau, valid = f["tau_rad"], f["tau_valid"]
            for w in range(tau.shape[1]):
                if not valid[:, w].all() or (member, str(w)) not in t2:
                    continue
                base, _ = en.autocorrelation(tau[:, w], True, 50)
                recomputed_max = max(recomputed_max, abs(base[1] - float(t2[(member, str(w))]["lag1"])))
                for s in SHIFTS:
                    c, _ = en.autocorrelation(np.angle(np.exp(1j * (tau[:, w] + s))), True, 50)
                    inv_max = max(inv_max, float(np.max(np.abs(c - base))))
                n_series += 1
        t1 = {(r["member"], r["window"]): float(r["lag1"]) for r in table(v1 / "tables/temporal.csv") if r["channel"] == "tau"}
        tau_change = max(abs(t1[k] - float(t2[k]["lag1"])) for k in t2 if k in t1)
        status_counts = {}
        for tname, col in (("distribution", "status"), ("replicas", "member_center_spread_status"),
                           ("between", "center_difference_status"), ("temporal", "status")):
            for r in table(v2 / f"tables/{tname}.csv"):
                key = f"{tname}:{r[col]}"
                status_counts[key] = status_counts.get(key, 0) + 1
        m = json.loads((v2 / "manifest.json").read_text())
        cores = {mm["core"]["native_extension"]["sha256"] for g in m["groups"].values() for mm in g}
        out["datasets"][name] = {
            "pooled_js_same_windows": same_keys, "pooled_js_max_abs_diff_v1_v2": max_js,
            "pooled_js_recorded_precision": "6 significant digits (CSV)",
            "tau_series_checked": n_series, "circular_autocorrelation_origin_invariance_max_abs_diff": inv_max,
            "tau_lag1_table_vs_recomputed_max_abs_diff": recomputed_max,
            "tau_lag1_max_change_v1_to_v2": tau_change,
            "status_counts": status_counts,
            "provenance": {"ensemble_version": m["software"]["ensemble_version"],
                           "module_sha256_matches_current_file": m["software"]["module_sha256"]
                           == hashlib.sha256((ROOT / "python/topofold_inputs/ensemble.py").read_bytes()).hexdigest(),
                           "native_extension_sha256_in_manifest": sorted(cores),
                           "native_extension_sha256_loaded_now": native["sha256"],
                           "thresholds": m["thresholds"], "pooling": m["pooling"]}}
    p = HERE / "results_v2" / "verification.json"
    p.parent.mkdir(exist_ok=True)
    p.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
