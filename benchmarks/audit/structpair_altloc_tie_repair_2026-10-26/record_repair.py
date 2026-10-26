"""Failing-before / passing-after record for the structpair optional altloc-sensitivity tie repair (2026-10-26).

    .venv/bin/python benchmarks/audit/structpair_altloc_tie_repair_2026-10-26/record_repair.py [--before-commit a155ada]

Runs tests/test_structpair_altloc_ties.py (a) against the committed code, in a temporary detached git worktree of the commit before the
correction, and (b) against the working tree. Records, per test, whether it passed. Nothing in the repository is modified;
the worktree is removed afterwards.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TEST = "tests/test_structpair_altloc_ties.py"
OUT = HERE / "results/repair_before_after.json"


def run(cwd: Path, adapter: str) -> dict:
    p = subprocess.run([sys.executable, "-m", "pytest", TEST, "-q", "-p", "no:cacheprovider", "-rA"], cwd=cwd, capture_output=True, text=True,
                       env=dict(os.environ, TOPOFOLD_ADAPTER_PYTHON=adapter))
    res = {}
    for ln in p.stdout.splitlines():
        m = re.match(r"^(PASSED|FAILED|SKIPPED|ERROR) [^:]+::(\S+)", ln)
        if m:
            res[m.group(2)] = m.group(1)
    return {"passed": sorted(k for k, v in res.items() if v == "PASSED"), "failed": sorted(k for k, v in res.items() if v == "FAILED"),
            "other": sorted(k for k, v in res.items() if v not in ("PASSED", "FAILED")), "exit_code": p.returncode}


def git(*a, cwd=ROOT) -> str:
    return subprocess.run(["git", *a], cwd=cwd, capture_output=True, text=True, check=True).stdout.strip()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--before-commit", default="a155ada")
    ap.add_argument("--adapter-python", default=os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python"))
    a = ap.parse_args(argv)
    before_full = git("rev-parse", a.before_commit)
    tmp = Path(tempfile.mkdtemp(prefix="structpair_before_"))
    wt = tmp / "wt"
    try:
        git("worktree", "add", "--detach", str(wt), before_full)
        shutil.copy(ROOT / TEST, wt / TEST)                                  # the new test file, run against the old code
        before = run(wt, a.adapter_python)
        before["code"] = {"commit": before_full, "structpair_version": re.search(r'STRUCTPAIR_VERSION = "([^"]+)"', (wt / "python/topofold_inputs/structpair.py").read_text())[1]}
    finally:
        subprocess.run(["git", "worktree", "remove", "--force", str(wt)], cwd=ROOT, capture_output=True)
        shutil.rmtree(tmp, ignore_errors=True)
    after = run(ROOT, a.adapter_python)
    after["code"] = {"base_commit": git("rev-parse", "HEAD"), "working_tree": "uncommitted repaired code",
                     "structpair_version": re.search(r'STRUCTPAIR_VERSION = "([^"]+)"', (ROOT / "python/topofold_inputs/structpair.py").read_text())[1]}
    rec = {"id": "STRUCTPAIR-ALTLOC-TIE-REPAIR-2026-10-26", "test_file": TEST, "before": before, "after": after,
           "failed_before_and_passing_after": sorted(set(before["failed"]) & set(after["passed"])),
           "passing_before_and_after": sorted(set(before["passed"]) & set(after["passed"])),
           "still_failing_after": after["failed"],
           "counts": {"before": {"passed": len(before["passed"]), "failed": len(before["failed"])},
                      "after": {"passed": len(after["passed"]), "failed": len(after["failed"])}}}
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps(rec["counts"]), "failed before and passing after:", len(rec["failed_before_and_passing_after"]))
    return 0 if not rec["still_failing_after"] else 1


if __name__ == "__main__":
    sys.exit(main())
