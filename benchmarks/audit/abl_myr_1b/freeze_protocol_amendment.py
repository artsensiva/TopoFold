"""Hash-freezes the dated ABL-MYR-1b protocol amendment a1 and the independent AI review record, before the wrapper is written.

    python3 benchmarks/audit/abl_myr_1b/freeze_protocol_amendment.py
"""

import datetime as dt
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / "results/ABL-MYR-1b_protocol_amendment_a1_freeze.json"
FILES = ["benchmarks/audit/abl_myr_1b/protocol_amendment_a1.json", "benchmarks/audit/abl_myr_1b/REVIEW_2026-10-25_independent_AI.md",
         "benchmarks/audit/abl_myr_1b/protocol_1b.json", "benchmarks/audit/abl_myr_1b/PROTOCOL_1B.md",
         "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_protocol_freeze.json",
         "benchmarks/audit/abl_myr_1/run_spec.json", "benchmarks/audit/abl_myr_1/run_abl_myr.py",
         "benchmarks/audit/abl_myr_1/results/ABL-MYR-1_code_freeze.json"]


def sha(p):
    return hashlib.sha256((ROOT / p).read_bytes()).hexdigest()


def main() -> int:
    if OUT.exists():
        print("freeze exists")
        return 1
    d = ROOT / "benchmarks/audit/abl_myr_1/data_dryad"
    pf = json.loads((HERE / "results/ABL-MYR-1b_protocol_freeze.json").read_text())
    if pf["protocol_sha256"] != sha(FILES[2]) or pf["protocol_md_sha256"] != sha(FILES[3]):
        raise SystemExit("STOP: the frozen protocol differs from its freeze record")
    rec = {"id": "ABL-MYR-1b-PROTOCOL-a1", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "sha256": {p: sha(p) for p in FILES},
           "trajectories_extracted": (d / "extracted").exists() or any(d.rglob("*.xtc")),
           "review_kind": "independent AI scientific review (not human peer review)"}
    if rec["trajectories_extracted"]:
        raise SystemExit("STOP: a trajectory has been extracted; the amendment must precede it")
    OUT.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps({"frozen_utc": rec["frozen_utc"], "amendment": rec["sha256"][FILES[0]], "review": rec["sha256"][FILES[1]]}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
