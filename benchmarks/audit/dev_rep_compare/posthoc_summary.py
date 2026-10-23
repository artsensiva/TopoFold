"""POST-HOC (not predeclared): distribution summaries of DEV-REP-1 primary per-window statistics.

Reads results/DEV-REP-1_results.json only and writes results/DEV-REP-1_posthoc_summary.json. Descriptive; changes no
frozen output and selects nothing.
    .venv/bin/python benchmarks/audit/dev_rep_compare/posthoc_summary.py
"""

import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
src = HERE / "results" / "DEV-REP-1_results.json"
o = json.loads(src.read_text())
out = {"posthoc": True, "source": src.name, "source_sha256": hashlib.sha256(src.read_bytes()).hexdigest(),
       "description": "per channel, over windows: fraction of the 25 between-temperature pairs above the largest of the "
                      "20 within-temperature pairs; ratio of median between to median within; pooled 320-vs-348 "
                      "distance and the widest leave-one-replica-out range", "channels": {}}
for ch, rows in o["per_window_primary"].items():
    fr = np.array([r["fraction_between_above_max_within"] for r in rows])
    ratio = np.array([r["median_between"] / r["median_within"] for r in rows])
    pooled = np.array([r["pooled"] for r in rows])
    span = np.array([r["pooled_loro_range"][1] - r["pooled_loro_range"][0] for r in rows])
    top = sorted(rows, key=lambda r: -r["fraction_between_above_max_within"])[:5]
    out["channels"][ch] = {
        "windows": len(rows),
        "fraction_between_above_max_within": {"median": float(np.median(fr)), "max": float(fr.max()),
                                              "windows_with_zero": int((fr == 0).sum())},
        "median_between_over_median_within": {"median": float(np.median(ratio)), "min": float(ratio.min()),
                                              "max": float(ratio.max()), "windows_above_1": int((ratio > 1).sum())},
        "pooled": {"median": float(np.median(pooled)), "max": float(pooled.max())},
        "pooled_loro_span_max": float(span.max()),
        "highest_fraction_windows": [{"residues": r["residues"], "fraction": r["fraction_between_above_max_within"],
                                      "pooled": r["pooled"]} for r in top],
    }
(HERE / "results" / "DEV-REP-1_posthoc_summary.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps({k: {kk: v[kk] for kk in ("fraction_between_above_max_within", "median_between_over_median_within",
                                            "pooled")} for k, v in out["channels"].items()}, indent=1))
