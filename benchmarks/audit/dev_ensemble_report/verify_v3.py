"""Verify v3 (exact exclusion causes) against v2 and the archived separate-agent gap run; writes results_v3/verification.json.

    .venv/bin/python benchmarks/audit/dev_ensemble_report/verify_v3.py
Reads git-ignored build/v3_<name>_a and build/v2_<name>_a outputs (regenerate with the commands in README.md).
"""

import csv
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REVIEW = HERE.parent / "pilot_review_2026-10-23"


def rows(p):
    return list(csv.DictReader(open(p)))


def key_set(p):
    return {(r["channel"], r["window"]) for r in rows(p)}


def npz_equal(a, b):
    x, y = np.load(a), np.load(b)
    return set(x.files) == set(y.files) and all(np.array_equal(x[k], y[k]) for k in x.files)


def main():
    out = {"note": "v3 adds exact exclusion causes; included windows, masks and pooled JS must be unchanged"}
    for name in ("1ha8A00", "adk"):
        v2, v3 = HERE / "build" / f"v2_{name}_a", HERE / "build" / f"v3_{name}_a"
        out[name] = {
            "tables_identical_to_v2": {t: (v2 / f"tables/{t}.csv").read_bytes() == (v3 / f"tables/{t}.csv").read_bytes()
                                       for t in ("distribution", "replicas", "between", "temporal")},
            "histograms_identical_to_v2": npz_equal(v2 / "distributions.npz", v3 / "distributions.npz"),
            "included_windows_identical": key_set(v2 / "tables/distribution.csv") == key_set(v3 / "tables/distribution.csv"),
            "exclusion_rows": len(rows(v3 / "tables/exclusions.csv")),
        }
    old, new = REVIEW / "gap_check", HERE / "build" / "v3_gap_1s14B00_a"
    ex = rows(new / "tables/exclusions.csv")
    out["gap_1s14B00"] = {
        "compared_with": "archived separate-agent run gap_check (tool v0.2.0, same spec)",
        "between_table_identical": (old / "tables/between.csv").read_bytes() == (new / "tables/between.csv").read_bytes(),
        "excluded_windows_identical": key_set(old / "tables/exclusions.csv") == key_set(new / "tables/exclusions.csv"),
        "excluded_windows": len(key_set(new / "tables/exclusions.csv")),
        "exclusion_rows": len(ex),
        "causes": sorted({f"{r['cause']} @ {r['location']}" for r in ex}),
        "members_with_causes": sorted({r["member"] for r in ex}),
        "source_frames_example": ex[0]["source_frames"] if ex else None,
        "version": json.loads((new / "manifest.json").read_text())["software"]["ensemble_version"],
    }
    p = HERE / "results_v3" / "verification.json"
    p.parent.mkdir(exist_ok=True)
    p.write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
