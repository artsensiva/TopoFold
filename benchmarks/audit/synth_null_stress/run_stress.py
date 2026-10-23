"""SYNTHETIC_NULLSTRESS-1: identical marginals, different temporal correlation between conditions (SYNTHETIC).

    .venv/bin/python benchmarks/audit/synth_null_stress/run_stress.py freeze
    .venv/bin/python benchmarks/audit/synth_null_stress/run_stress.py run

Spec: docs/NULL_STRESS_CHECK_SPEC.md and run_spec.json.
"""

import datetime as dt
import hashlib
import json
import math
import platform
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/synth_rep_char"))
from synthetic_generators import circular_series, linear_series  # noqa: E402
from topofold_inputs.compare import circular_hist, js_distance, linear_hist, replica_comparison  # noqa: E402
from topofold_inputs.replica_permutation import permutation_test  # noqa: E402

SPEC_JSON = HERE / "run_spec.json"
SPEC_DOC = ROOT / "docs/NULL_STRESS_CHECK_SPEC.md"
FREEZE = HERE / "results/SYNTHETIC_NULLSTRESS-1_freeze.json"
RESULTS = HERE / "results/SYNTHETIC_NULLSTRESS-1_results.json"
REPORT = HERE / "results/SYNTHETIC_NULLSTRESS-1_report.md"
CODE = ["benchmarks/audit/synth_null_stress/run_stress.py", "benchmarks/audit/synth_rep_char/synthetic_generators.py",
        "python/topofold_inputs/compare.py", "python/topofold_inputs/replica_permutation.py"]
CATS = ("consistent", "replica_dependent", "inconclusive")


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def wilson(k, n, z=1.959963984540054):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [max(0.0, c - h), min(1.0, c + h)]


def freeze():
    if FREEZE.exists():
        print("freeze exists")
        return 1
    rec = {"id": "SYNTHETIC_NULLSTRESS-1", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "spec_json_sha256": sha256(SPEC_JSON), "spec_doc_sha256": sha256(SPEC_DOC),
           "code_sha256": {c: sha256(ROOT / c) for c in CODE},
           "environment": {"python": platform.python_version(), "numpy": np.__version__, "platform": platform.platform()}}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=1))
    return 0


def draw(spec, cell, rep):
    rng = np.random.default_rng(np.random.SeedSequence([spec["seed_master"], cell["cell"], rep]))
    n, T = spec["n_rep"], spec["n_frames"]
    if cell["channel"] == "linear":
        g = lambda phi: linear_series(rng, T, spec["linear"]["mu"], spec["linear"]["sigma"], phi)  # noqa: E731
        h = lambda x: linear_hist(x, spec["linear"]["lo"], spec["linear"]["hi"], spec["bins"])[0]  # noqa: E731
    else:
        g = lambda phi: circular_series(rng, T, spec["circular"]["mu"], spec["circular"]["s"], phi)  # noqa: E731
        h = lambda x: circular_hist(x, spec["bins"], spec["circular"]["origin"])  # noqa: E731
    A = [h(g(cell["phi_a"])) for _ in range(n)]
    B = [h(g(cell["phi_b"])) for _ in range(n)]
    return A, B


def run():
    spec = json.loads(SPEC_JSON.read_text())
    frz = json.loads(FREEZE.read_text())
    bad = [c for c, v in frz["code_sha256"].items() if sha256(ROOT / c) != v]
    if frz["spec_json_sha256"] != sha256(SPEC_JSON) or frz["spec_doc_sha256"] != sha256(SPEC_DOC) or bad:
        raise SystemExit(f"STOP: spec or code differs from the freeze: {bad}")
    cap = spec["budget"]["cpu_seconds"]
    resource.setrlimit(resource.RLIMIT_CPU, (cap, cap))
    t0 = time.time()
    done, per_rep, not_run, stopped = [], [], [], False

    def cpu():
        r = resource.getrusage(resource.RUSAGE_SELF)
        return r.ru_utime + r.ru_stime

    def write(final):
        doc = {"id": "SYNTHETIC_NULLSTRESS-1", "label": "SYNTHETIC",
               "status": "complete" if final and len(done) == len(spec["cells"]) else "partial",
               "cells_planned": len(spec["cells"]), "cells_done": len(done), "cells_not_run": not_run,
               "freeze_utc": frz["frozen_utc"], "alpha": spec["alpha"], "environment": frz["environment"],
               "resources": {"cpu_s": cpu(), "wall_s": time.time() - t0, "cap_cpu_s": cap}, "cells": done}
        RESULTS.write_text(json.dumps(doc, indent=1) + "\n")
        return doc

    R = spec["repetitions"]
    for cell in spec["cells"]:
        if per_rep and cpu() + 1.5 * float(np.mean(per_rep)) * R > cap:
            not_run = [c["cell"] for c in spec["cells"] if c["cell"] >= cell["cell"]]
            stopped = True
            break
        c0 = cpu()
        rej, cats, wa, wb = 0, [], [], []
        for rep in range(R):
            A, B = draw(spec, cell, rep)
            rej += permutation_test(A, B)["p"] <= spec["alpha"]
            cats.append(replica_comparison(A, B)["category"])
            wa.append(np.median([js_distance(A[i], A[j]) for i in range(5) for j in range(i + 1, 5)]))
            wb.append(np.median([js_distance(B[i], B[j]) for i in range(5) for j in range(i + 1, 5)]))
        counts = {c: cats.count(c) for c in CATS}
        done.append({"cell": cell, "repetitions": R,
                     "permutation": {"rejections": int(rej), "frequency": rej / R, "wilson95": wilson(int(rej), R)},
                     "screen": {"counts": counts, "frequency": {c: counts[c] / R for c in CATS},
                                "wilson95": {c: wilson(counts[c], R) for c in CATS}},
                     "median_within_distance": {"A": float(np.median(wa)), "B": float(np.median(wb))}})
        per_rep.append((cpu() - c0) / R)
        write(final=False)
    doc = write(final=not stopped)
    L = ["# SYNTHETIC_NULLSTRESS-1: identical marginals, different temporal correlation (SYNTHETIC)", "",
         f"Generated from `SYNTHETIC_NULLSTRESS-1_results.json`. Status **{doc['status']}** ({doc['cells_done']} of "
         f"{doc['cells_planned']} cells). α = {doc['alpha']}; nominal permutation size 12/252 ≈ 0.048 under exact "
         "exchangeability, which these cells deliberately violate (the marginals are equal; the temporal correlation is "
         "not). 200 repetitions per cell, which gives coarse Monte Carlo precision. No real data.", "",
         f"Resources: {doc['resources']['cpu_s']:.1f} CPU-s (cap {cap}), wall {doc['resources']['wall_s']:.1f} s.", "",
         "| Channel | φ_A | φ_B | Permutation reject [95 % CI] | Screen C | Screen RD | Screen I | median within-A / within-B JS distance |",
         "|---|---|---|---|---|---|---|---|"]
    for r in doc["cells"]:
        c, p, s = r["cell"], r["permutation"], r["screen"]["frequency"]
        L.append(f"| {c['channel']} | {c['phi_a']} | {c['phi_b']} | {p['frequency']:.3f} [{p['wilson95'][0]:.3f}, "
                 f"{p['wilson95'][1]:.3f}] | {s['consistent']:.3f} | {s['replica_dependent']:.3f} | {s['inconclusive']:.3f} | "
                 f"{r['median_within_distance']['A']:.3f} / {r['median_within_distance']['B']:.3f} |")
    REPORT.write_text("\n".join(L) + "\n")
    print(json.dumps({"status": doc["status"], "cells_done": doc["cells_done"], "resources": doc["resources"]}, indent=1))
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("freeze", "run"):
        print(__doc__)
        sys.exit(2)
    sys.exit(freeze() if sys.argv[1] == "freeze" else run())
