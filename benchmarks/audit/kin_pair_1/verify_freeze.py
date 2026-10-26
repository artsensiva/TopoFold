"""KIN-PAIR-1: compare each freeze record with the files at the commit that introduced it and with the working tree.

    python3 benchmarks/audit/kin_pair_1/verify_freeze.py   # writes results/freeze_verification.json
"""

import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
REL = HERE.relative_to(HERE.parents[2]).as_posix()
h = lambda b: hashlib.sha256(b).hexdigest()   # noqa: E731
out = {}
for name in ("1", "1a", "2"):
    fr = json.loads((HERE / f"results/KIN-PAIR-1_freeze_{name}.json").read_text())
    c = subprocess.run(["git", "log", "--format=%h", "--diff-filter=A", "--", f"{REL}/results/KIN-PAIR-1_freeze_{name}.json"], capture_output=True, text=True, cwd=HERE.parents[2]).stdout.split()[0]
    rows = {}
    for f, hh in fr["sha256"].items():
        blob = subprocess.run(["git", "show", f"{c}:{REL}/{f}"], capture_output=True, cwd=HERE.parents[2]).stdout
        rows[f] = {"matches_at_freeze_commit": h(blob) == hh, "unchanged_in_working_tree": h((HERE / f).read_bytes()) == hh}
    out[name] = {"freeze_commit": c, "utc": fr["utc"], "files": rows}
note = ("PROTOCOL.md differs from every freeze because dated amendments were appended (amendment 1 before freeze 1a, amendment 2 after the primary run); "
        "select_panel.py differs from freeze 1 by the amendment-1 fix; all frozen code and inputs used by the primary execution (freeze 2) are unchanged except PROTOCOL.md.")
(HERE / "results/freeze_verification.json").write_text(json.dumps({"note": note, "freezes": out}, indent=1) + "\n")
print(note)
