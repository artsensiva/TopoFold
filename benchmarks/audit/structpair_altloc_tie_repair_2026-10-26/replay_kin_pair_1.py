"""Engineering regression replay of the six KIN-PAIR-1 pairs with the repaired structpair (0.3.1).

    .venv/bin/python benchmarks/audit/structpair_altloc_tie_repair_2026-10-26/replay_kin_pair_1.py

Uses the ORIGINAL primary-run specifications (`kin_pair_1/work/<pair>/spec.json`, which include the optional
`altloc_sensitivity: highest_occupancy` request) and the original inputs; nothing under kin_pair_1 is modified. This is a regression
test of the repair, not a new evaluation and not a redefinition of the KIN-PAIR-1 outcome. Primary execution and optional-sensitivity
availability are reported separately. New outputs: results/replay/<pair>/ and results/replay_summary.json.
"""

import csv
import hashlib
import json
import os
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
KP = ROOT / "benchmarks/audit/kin_pair_1"
ADAPTER = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
OUT = HERE / "results/replay"
KEEP = ("windows.csv", "residues.csv", "exclusions.csv", "summary.json", "manifest.json", "report.md", "chain_variation.csv", "chains.csv")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    primary = json.loads((KP / "results/KIN-PAIR-1_results.json").read_text())
    supp = {p["pair_id"]: p for p in json.loads((KP / "results/KIN-PAIR-1_supplement.json").read_text())["pairs"]}
    shutil.rmtree(OUT, ignore_errors=True)
    OUT.mkdir(parents=True)
    rows = []
    c0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = time.time()
    for p in primary["pairs"]:
        pid = p["pair_id"]
        spec = KP / "work" / pid / "spec.json"
        s = json.loads(spec.read_text())
        assert s.get("altloc_sensitivity") == "highest_occupancy", pid           # the original request is kept
        out = HERE / "work" / f"replay_{pid}"
        shutil.rmtree(out, ignore_errors=True)
        out.parent.mkdir(exist_ok=True)
        r = subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(spec), "--out", str(out)], capture_output=True, text=True,
                           env=dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=ADAPTER))
        row = {"pair_id": pid, "frozen_primary_outcome": p["outcome"], "original_stderr": p.get("structpair", {}).get("stderr", ""),
               "primary_execution": {"returncode": r.returncode, "stdout": r.stdout.strip()[-300:], "stderr": r.stderr.strip()[-300:]}}
        if r.returncode == 0:
            ref = (KP / "results/pairs" / pid) if p["outcome"] == "PASS" else (KP / "results/pairs_supplement" / pid)
            row["reference_outputs"] = str(ref.relative_to(ROOT))
            row["primary_outputs_byte_identical_to_reference"] = {f: sha(out / f) == sha(ref / f) for f in ("windows.csv", "residues.csv", "exclusions.csv")}
            summ = json.loads((out / "summary.json").read_text())
            ref_summ = json.loads((ref / "summary.json").read_text())
            row["summary_equal_except_altloc_sensitivity"] = all(summ[k] == ref_summ[k] for k in summ if k != "altloc_sensitivity")
            sens = summ["altloc_sensitivity"]
            row["optional_sensitivity"] = {"status": sens["status"], "available": sens.get("available", sens["status"] == "computed"),
                                           "n_ties": len(sens.get("ties", [])), "ties_in_correspondence": sum(t["in_correspondence"] for t in sens.get("ties", [])),
                                           "structures_with_ties": sorted({t["structure"] for t in sens.get("ties", [])}),
                                           "values_reported": sens.get("values_reported", "max_abs_change" in sens),
                                           "max_abs_change": sens.get("max_abs_change"), "ties": sens.get("ties")}
            if p["outcome"] == "PASS":
                row["sensitivity_equals_original_primary_run"] = sens == ref_summ["altloc_sensitivity"]
            man = json.loads((out / "manifest.json").read_text())
            row["structpair_version"] = man["software"]["structpair_version"]
            row["manifest_altloc_sensitivity_outcome"] = man["policies"].get("altloc_sensitivity_outcome")
            row["manifest_hashes_verified"] = all(sha(out / f) == h for f, h in man["outputs_sha256"].items())
            dst = OUT / pid
            dst.mkdir()
            for f in KEEP:
                if (out / f).is_file():
                    shutil.copy2(out / f, dst / f)
        rows.append(row)
        print(pid, row["primary_execution"]["returncode"], row.get("optional_sensitivity", {}).get("status"), flush=True)
    c1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    summary = {"id": "KIN-PAIR-1-REPLAY-REPAIRED-STRUCTPAIR-2026-10-26", "note": "engineering regression of the repair; the frozen KIN-PAIR-1 outcome is unchanged",
               "pairs": rows,
               "counts": {"pairs": len(rows), "primary_execution_completed": sum(r["primary_execution"]["returncode"] == 0 for r in rows),
                          "primary_execution_refused": sum(r["primary_execution"]["returncode"] != 0 for r in rows),
                          "optional_sensitivity_computed": sum(r.get("optional_sensitivity", {}).get("status") == "computed" for r in rows),
                          "optional_sensitivity_not_applicable": sum(str(r.get("optional_sensitivity", {}).get("status", "")).startswith("not applicable") for r in rows),
                          "optional_sensitivity_unavailable": sum(str(r.get("optional_sensitivity", {}).get("status", "")).startswith("unavailable") for r in rows)},
               "resources": {"children_cpu_seconds": (c1.ru_utime + c1.ru_stime) - (c0.ru_utime + c0.ru_stime), "wall_seconds": time.time() - t0}}
    (HERE / "results/replay_summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary["counts"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
