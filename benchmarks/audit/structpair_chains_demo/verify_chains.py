"""Engineering verification of structpair 0.3.0 (multiple chains, corrected advisory output) on the existing KRAS inputs and records.

    .venv/bin/python benchmarks/audit/structpair_chains_demo/verify_chains.py [--adapter-python PATH]

This checks that the implementation is consistent with earlier records; it does not check biological validity. Needs the
git-ignored local copies of the KRAS-STRUCT-1/-2 structures. Checks:
  V1  the four committed single-chain demonstrations of version 0.1.0, re-run: windows.csv, residues.csv and exclusions.csv
      byte-identical to the committed outputs, summary numbers equal;
  V2  the KRAS-STRUCT-2 chains (5US4 A, B; 5USJ A-F): per window metric, the primary value, the cross min/median/max over 12
      pairings and the within-entry maximum over 16 copy pairs equal the frozen KRAS-STRUCT-2_table.csv;
  V3  the exact K1 chains, residues and settings: displacement mean, median, primary pair, within-entry maxima, threshold and the
      flagged set equal the frozen KRAS-K1_endpoint.json.
V3 also applies the frozen K1 gate function to the flagged set and checks that its (failed) gate result equals the recorded one.
Writes results/v0_3_0/verification.json and the compact K1-regression outputs with a caption generated from the record. The
records of version 0.2.0 (schema /2, results/verification.json and results/k1_regression/) are historical and are not rewritten.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DEMO = ROOT / "benchmarks/audit/structpair_demo"
STRUCT2 = ROOT / "benchmarks/audit/kras_struct_g12d/results/KRAS-STRUCT-2_table.csv"
K1 = ROOT / "benchmarks/audit/kras_k1/results/KRAS-K1_endpoint.json"
OUT = HERE / "results" / "v0_3_0"      # 0.2.0 records stay in results/ (historical, schema /2)
METRICS = ("d_kappa", "d_i_i2", "rmsd3", "d_tau", "d_i_i3", "rmsd4")
LEGACY = {"kras1": "kras1_9O0N_vs_9O0O.json", "kras2": "kras2_5US4A_vs_5USJA.json",
          "explicit_exclude": "kras2_5US4A_vs_5USJA_explicit_exclude.json", "mismatch_flag": "demo_9O0NA_vs_5US4B_mismatch_flag.json"}
NUMERIC = ("top", "agreement", "n_windows", "global", "altloc_sensitivity", "n_exclusion_rows")
TOL = 1e-9
# V3 compares with the frozen K1 calculation, which parses the PDB text in float64, while the adapter (MDAnalysis) holds float32
# positions: the declared tolerance is 1e-4 A (float32 spacing at 60 A is 3.8e-6 A, PDB coordinates have 1e-3 A resolution)
TOL_K1 = 1e-4


def run(spec: Path, out: Path, adapter: str) -> dict:
    t0, r0 = time.time(), resource.getrusage(resource.RUSAGE_CHILDREN)
    p = subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(spec), "--out", str(out)], capture_output=True, text=True,
                       env=dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=adapter))
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    if p.returncode != 0:
        raise SystemExit(f"STOP: {spec.name}: {p.stderr.strip()}")
    size = lambda d: sum(f.stat().st_size for f in d.rglob("*") if f.is_file())  # noqa: E731
    return {"wall_s": round(time.time() - t0, 2), "cpu_s": round((r1.ru_utime + r1.ru_stime) - (r0.ru_utime + r0.ru_stime), 2),
            "output_bytes_total": size(out), "output_bytes_without_work": size(out) - size(out / "work")}


def num(x):
    return None if x in ("", None) else float(x)


def v1(adapter: str, tmp: Path) -> dict:
    res = {}
    for name, spec in LEGACY.items():
        o = tmp / f"legacy_{name}"
        rt = run(DEMO / "specs" / spec, o, adapter)
        rec = DEMO / "results" / name
        same = {f: (o / f).read_bytes() == (rec / f).read_bytes() for f in ("windows.csv", "residues.csv", "exclusions.csv")}
        new, old = json.loads((o / "summary.json").read_text()), json.loads((rec / "summary.json").read_text())
        res[name] = {"csv_byte_identical": same, "summary_numeric_equal": {k: new[k] == old[k] for k in NUMERIC},
                     "within_entry_variation": new["within_entry_variation"]["status"], **rt}
    res["all_identical"] = all(all(v["csv_byte_identical"].values()) and all(v["summary_numeric_equal"].values())
                               for k, v in res.items() if k != "all_identical")
    return res


def long_rows(path: Path) -> dict:
    out = {}
    for r in csv.DictReader(open(path)):
        out[(int(r["A_residue"].split(":")[1]), r["channel"])] = r
    return out


def v2(adapter: str, tmp: Path) -> dict:
    o = tmp / "struct2_chains"
    rt = run(HERE / "specs/kras2_all_chains_5US4_AB_vs_5USJ_AF.json", o, adapter)
    rows = long_rows(o / "chain_variation.csv")
    cells, mism, maxd = 0, [], 0.0
    pairs = {}
    for t in csv.DictReader(open(STRUCT2)):
        rid = int(t["resid"])
        for k in METRICS:
            r = rows[(rid, k)]
            for tcol, ncol in ((k, "primary"), (f"{k}_cross_min", "between_min"), (f"{k}_cross_median", "between_median"),
                               (f"{k}_cross_max", "between_max"), (f"{k}_copy_max", "within_max")):
                a, b = num(t[tcol]), num(r[ncol])
                cells += 1
                if (a is None) != (b is None) or (a is not None and abs(a - b) > TOL):
                    mism.append((rid, k, ncol))
                elif a is not None:
                    maxd = max(maxd, abs(a - b))
            pairs[(r["between_n"], r["withinA_n"], r["withinB_n"])] = pairs.get((r["between_n"], r["withinA_n"], r["withinB_n"]), 0) + 1
    return {"cells_compared": cells, "mismatches": len(mism), "first_mismatches": mism[:5], "max_abs_difference": maxd,
            "pair_counts_seen": {"/".join(k): v for k, v in pairs.items()}, "tolerance": TOL, **rt}


def coordinate_check(o: Path) -> dict:
    """Largest difference between the adapter's coordinates and a float64 parse of the PDB text, over the 8 chains and 153 positions."""
    import numpy as np
    sys.path.insert(0, str(ROOT / "benchmarks/audit/kras_k1"))
    import k1_endpoint as K
    spec = json.loads((HERE / "specs/k1_exact_chains_and_settings.json").read_text())
    worst = 0.0
    for entry, tag in (("A", "A"), ("B", "B")):
        s = spec["structures"][entry]
        ca, _ = K.parse_ca((HERE / "specs" / s["file"]).read_text())
        for ch in s["chains"]:
            run_id = f"primary_{entry}" if ch == s["primary_chain"] else f"chain_{entry}_{ch}"
            a = np.load(o / "work" / f"{run_id}.npz")["coordinates_angstrom"][0]
            res = [r for r in range(4, 165) if not 60 <= r <= 67]
            worst = max(worst, float(np.abs(a - np.stack([ca[ch][r] for r in res])).max()))
    return {"max_abs_coordinate_difference_A": worst, "adapter_value_precision": "float32 (MDAnalysis positions), exported as float64"}


def v3(adapter: str, tmp: Path) -> dict:
    o = tmp / "k1_regression"
    rt = run(HERE / "specs/k1_exact_chains_and_settings.json", o, adapter)
    coords = coordinate_check(o)
    rows = long_rows(o / "chain_variation.csv")
    rec = json.loads(K1.read_text())["per_residue"]
    m = {"between_mean": "b_mean", "between_median": "b_median", "primary": "b_primary_pair", "within_max": "w_max",
         "withinA_max": "w_gdp_entry_max", "withinB_max": "w_gnp_entry_max", "advisory_threshold": "threshold"}
    cells, mism, maxd = 0, [], 0.0
    for i, rid in enumerate(rec["resid"]):
        r = rows[(rid, "disp")]
        for ncol, kcol in m.items():
            a, b = float(r[ncol]), rec[kcol][i]
            cells += 1
            if abs(a - b) > TOL_K1:
                mism.append((rid, ncol))
            maxd = max(maxd, abs(a - b))
        cells += 1
        if (r["advisory_flag"] == "True") != rec["positive"][i]:
            mism.append((rid, "advisory_flag"))
    flagged = sorted(rid for (rid, ch), r in rows.items() if ch == "disp" and r["advisory_flag"] == "True")
    k1rec = json.loads(K1.read_text())
    recorded = k1rec["E"]
    boundary = sorted(rid for (rid, ch), r in rows.items() if ch == "disp" and r["advisory_at_boundary"] == "True")
    sys.path.insert(0, str(ROOT / "benchmarks/audit/kras_k1"))
    import k1_endpoint as K
    params = json.loads((ROOT / "benchmarks/audit/kras_k1/protocol_params.json").read_text())["endpoint"]
    gate = K.gate(flagged, len(rec["resid"]), params["gate"], params["segments"]["merge_gap"])      # the frozen K1 gate, read-only use
    keep = ("chain_variation.csv", "chains.csv", "chain_variation.png", "windows.csv", "residues.csv", "exclusions.csv", "report.md",
            "manifest.json", "summary.json")
    dest = OUT / "k1_regression"
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    for f in keep:
        shutil.copy(o / f, dest / f)
    return {"positions": len(rec["resid"]), "cells_compared": cells, "mismatches": len(mism), "first_mismatches": mism[:5],
            "max_abs_difference": maxd, "flagged_equals_recorded_E": flagged == recorded, "flagged": flagged, "recorded_E": recorded,
            "positions_at_boundary": boundary,
            "k1_gate_unchanged": {"gate_equals_recorded": gate == k1rec["gate"], "passes": gate["passes"],
                                  "segment_lengths": [x["length"] for x in gate["segments"]], "recorded_decision": k1rec["decision"]},
            "tolerance_A": TOL_K1, "tolerance_reason": "the frozen K1 code reads the PDB text in float64; the adapter holds float32 positions",
            "coordinate_check": coords,
            "note": "implementation consistency with the frozen K1 calculation, not biological validity", **rt}


def caption(rec: dict) -> str:
    k = rec["k1_regression"]
    return f"""# Figure: within-entry crystal-copy variation, K1 regression run (structpair 0.3.0)

![chain variation](chain_variation.png)

**Purpose.** To show what the new multi-chain output looks like on the chains and residues of the recorded K1 calculation, and
that it reproduces that calculation.

**What it shows.** Per residue (4–59 and 68–164), the Cα displacement between 5US4 (GDP) and 5USJ (GMPPNP) after a trimmed
fit: the min–max band and the median over the 12 chain pairings, and the maximum over the within-entry chain pairs of each
entry (dashed). Marks at the top are the residues for which the advisory rule, configured with the K1 settings (minimum difference 1 Å on the
displacement; the between-entry mean at least max(1 Å, 2 × the within-entry maximum)), is met. They are {len(k['flagged'])} residues ({', '.join(map(str, k['flagged']))}); the recorded K1 endpoint set is
{'identical' if k['flagged_equals_recorded_E'] else 'DIFFERENT'}. Across {k['cells_compared']} compared values the largest absolute
difference from the K1 record is {k['max_abs_difference']:.1e} Å (declared tolerance {k['tolerance_A']:g} Å: the adapter holds
float32 positions, which differ from a float64 parse of the PDB text by up to {k['coordinate_check']['max_abs_coordinate_difference_A']:.1e} Å).

**Main limitation.** Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of
them. The dashed lines are descriptive ranges among the selected chains, not a calibrated noise floor, a statistical
uncertainty, a nucleotide effect or an HDX replicate estimate. Agreement with the K1 record shows implementation
consistency, not biological validity. The marks are an advisory rule the user configured, not significance and not biological
separation. The frozen K1 gate applied to the marked set gives its recorded result ({'identical' if k['k1_gate_unchanged']['gate_equals_recorded'] else 'DIFFERENT'}; {'passes' if k['k1_gate_unchanged']['passes'] else 'fails'}, K1_ENDPOINT_INADEQUATE unchanged).

**Source data.** `chain_variation.csv` and `manifest.json` in this directory (inputs 5US4 and 5USJ, SHA-256 in the manifest);
`benchmarks/audit/kras_k1/results/KRAS-K1_endpoint.json`.

**Reproduce.** `.venv/bin/python benchmarks/audit/structpair_chains_demo/verify_chains.py` (about {k['cpu_s']:.0f} CPU-seconds for
this run; the figure is written by the command itself).
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--adapter-python", default=os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python"))
    a = ap.parse_args(argv)
    with tempfile.TemporaryDirectory(prefix="structpair_chains_") as t:
        tmp = Path(t)
        rec = {"id": "STRUCTPAIR-CHAINS-VERIFY-1", "date": "2026-10-25", "scope": "implementation consistency, not biological validity",
               "structpair_version": "0.3.0", "legacy_single_chain": v1(a.adapter_python, tmp), "struct2_chain_aggregates": v2(a.adapter_python, tmp),
               "k1_regression": v3(a.adapter_python, tmp)}
    OUT.mkdir(exist_ok=True)
    (OUT / "k1_regression/chain_variation_caption.md").write_text(caption(rec))
    (OUT / "verification.json").write_text(json.dumps(rec, indent=1) + "\n")
    s2, k1 = rec["struct2_chain_aggregates"], rec["k1_regression"]
    print(json.dumps({"legacy_all_identical": rec["legacy_single_chain"]["all_identical"],
                      "struct2": {k: s2[k] for k in ("cells_compared", "mismatches", "max_abs_difference", "wall_s", "cpu_s")},
                      "k1": {k: k1[k] for k in ("cells_compared", "mismatches", "max_abs_difference", "flagged_equals_recorded_E", "k1_gate_unchanged", "wall_s", "cpu_s",
                                                "output_bytes_total", "output_bytes_without_work")}}, indent=1))
    ok = (rec["legacy_single_chain"]["all_identical"] and s2["mismatches"] == 0 and k1["mismatches"] == 0 and k1["flagged_equals_recorded_E"]
          and k1["k1_gate_unchanged"]["gate_equals_recorded"] and not k1["k1_gate_unchanged"]["passes"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
