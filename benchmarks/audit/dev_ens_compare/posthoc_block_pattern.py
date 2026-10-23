"""POST-HOC (not predeclared): tabulate which block pairs drive the DEV-ENS-1 block fractions for the top-10 windows.

Reads only results/DEV-ENS-1_results.json and writes results/DEV-ENS-1_posthoc_block_pattern.json. Descriptive;
it changes no frozen output and selects no configuration.
    .venv/bin/python benchmarks/audit/dev_ens_compare/posthoc_block_pattern.py
"""

import collections
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
src = HERE / "results" / "DEV-ENS-1_results.json"
o = json.loads(src.read_text())
out = {"posthoc": True, "source": src.name, "source_sha256": hashlib.sha256(src.read_bytes()).hexdigest(),
       "description": "for the top-10 primary windows per channel: block pairs (open i, closed j) whose between-run "
                      "distance is at or below the window's largest within-run distance, and which within-run "
                      "block pair is that largest one", "channels": {}}
pairs_w = ["open 0-1", "open 0-2", "open 1-2", "closed 0-1", "closed 0-2", "closed 1-2"]
for ch, top in o["top10_primary"].items():
    rows = {r["window"]: r for r in o["per_window_primary"][ch]}
    at_or_below, max_within_pair = collections.Counter(), collections.Counter()
    for t in top:
        r = rows[t["window"]]
        mx = r["max_block_within"]
        for k, v in enumerate(r["block_between"]):
            if v <= mx:
                at_or_below[f"open{k // 3}-closed{k % 3}"] += 1
        max_within_pair[pairs_w[r["block_within"].index(mx)]] += 1
    out["channels"][ch] = {"between_pairs_at_or_below_max_within": dict(sorted(at_or_below.items())),
                           "largest_within_pair": dict(sorted(max_within_pair.items()))}
(HERE / "results" / "DEV-ENS-1_posthoc_block_pattern.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps(out["channels"], indent=1))
