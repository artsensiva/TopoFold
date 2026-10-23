"""SYNTHETIC_PERM-1: historical complete-separation screen vs prospective replica permutation method (SYNTHETIC).

    .venv/bin/python benchmarks/audit/synth_rep_perm/run_perm.py freeze
    .venv/bin/python benchmarks/audit/synth_rep_perm/run_perm.py run

Spec: docs/REPLICA_PERMUTATION_EVALUATION_SPEC.md and run_spec.json. Draws come from the frozen REPCHAR-1 simulate().
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
CHAR = ROOT / "benchmarks/audit/synth_rep_char"
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(CHAR))
import run_char  # noqa: E402  (frozen REPCHAR-1 runner; only simulate() and wilson() are used)
from topofold_inputs.compare import replica_comparison  # noqa: E402
from topofold_inputs.replica_permutation import permutation_test  # noqa: E402

SPEC_JSON = HERE / "run_spec.json"
SPEC_DOC = ROOT / "docs/REPLICA_PERMUTATION_EVALUATION_SPEC.md"
CHAR_SPEC = CHAR / "run_spec.json"
CHAR_RESULTS = CHAR / "results/SYNTHETIC_REPCHAR-1_results.json"
FREEZE = HERE / "results/SYNTHETIC_PERM-1_freeze.json"
RESULTS = HERE / "results/SYNTHETIC_PERM-1_results.json"
REPORT = HERE / "results/SYNTHETIC_PERM-1_report.md"
CODE = ["python/topofold_inputs/replica_permutation.py", "python/topofold_inputs/compare.py",
        "benchmarks/audit/synth_rep_char/run_char.py", "benchmarks/audit/synth_rep_char/synthetic_generators.py",
        "benchmarks/audit/synth_rep_perm/run_perm.py", "tests/test_replica_permutation.py"]
CATS = ("consistent", "replica_dependent", "inconclusive")


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def freeze():
    if FREEZE.exists():
        print("freeze exists")
        return 1
    rec = {"id": "SYNTHETIC_PERM-1", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "spec_json_sha256": sha256(SPEC_JSON), "spec_doc_sha256": sha256(SPEC_DOC),
           "repchar_spec_sha256": sha256(CHAR_SPEC), "repchar_results_sha256": sha256(CHAR_RESULTS),
           "code_sha256": {c: sha256(ROOT / c) for c in CODE},
           "environment": {"python": platform.python_version(), "numpy": np.__version__,
                           "platform": platform.platform(), "machine": platform.machine()}}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=1))
    return 0


def run_cell(cell, char_spec, alpha):
    R = char_spec["repetitions"]
    cats, rej, ps = [], [], []
    agree = {"reject&consistent": 0, "reject&replica_dependent": 0, "reject&inconclusive": 0,
             "accept&consistent": 0, "accept&replica_dependent": 0, "accept&inconclusive": 0}
    for rep in range(R):
        ha, hb, _ = run_char.simulate(cell, char_spec, rep)       # identical draws for both methods
        cat = replica_comparison(ha, hb)["category"]
        pt = permutation_test(ha, hb)
        r = pt["p"] <= alpha
        cats.append(cat)
        rej.append(r)
        ps.append(pt["p"])
        agree[("reject&" if r else "accept&") + cat] += 1
    counts = {c: cats.count(c) for c in CATS}
    n_rej = int(sum(rej))
    return {"cell": cell, "repetitions": R,
            "screen": {"counts": counts, "frequency": {c: counts[c] / R for c in CATS},
                       "wilson95": {c: run_char.wilson(counts[c], R) for c in CATS}},
            "permutation": {"rejections": n_rej, "frequency": n_rej / R, "wilson95": run_char.wilson(n_rej, R),
                            "p_median": float(np.median(ps)), "p_min": float(np.min(ps))},
            "agreement": agree}


def run():
    spec = json.loads(SPEC_JSON.read_text())
    char_spec = json.loads(CHAR_SPEC.read_text())
    frz = json.loads(FREEZE.read_text())
    bad = [c for c, h in frz["code_sha256"].items() if sha256(ROOT / c) != h]
    if (frz["spec_json_sha256"] != sha256(SPEC_JSON) or frz["spec_doc_sha256"] != sha256(SPEC_DOC)
            or frz["repchar_spec_sha256"] != sha256(CHAR_SPEC) or bad):
        raise SystemExit(f"STOP: spec or code differs from the freeze: {bad}")
    hist_cells = {r["cell"]["cell"]: r for r in json.loads(CHAR_RESULTS.read_text())["cells"]}
    cap = spec["budget"]["cpu_seconds"]
    resource.setrlimit(resource.RLIMIT_CPU, (cap, cap))
    t0 = time.time()
    done, per_rep, not_run, stopped = [], [], [], False

    def cpu():
        r = resource.getrusage(resource.RUSAGE_SELF)
        return r.ru_utime + r.ru_stime

    def write(final):
        status = "complete" if final and len(done) == len(spec["cells"]) else "partial"
        doc = {"id": "SYNTHETIC_PERM-1", "label": "SYNTHETIC", "status": status, "cells_planned": len(spec["cells"]),
               "cells_done": len(done), "cells_not_run": not_run, "freeze_utc": frz["frozen_utc"],
               "alpha": spec["alpha"], "environment": frz["environment"],
               "resources": {"cpu_s": cpu(), "wall_s": time.time() - t0, "cap_cpu_s": cap,
                             "mean_cpu_s_per_repetition": float(np.mean(per_rep)) if per_rep else None},
               "cells": done}
        RESULTS.write_text(json.dumps(doc, indent=1) + "\n")
        return doc

    for cell in spec["cells"]:
        if per_rep and cpu() + 1.5 * float(np.mean(per_rep)) * char_spec["repetitions"] > cap:
            not_run = [c["cell"] for c in spec["cells"] if c["cell"] >= cell["cell"]]
            stopped = True
            break
        c0 = cpu()
        res = run_cell({k: v for k, v in cell.items() if k != "role"}, char_spec, spec["alpha"])
        res["role"] = cell["role"]
        h = hist_cells[cell["cell"]]["counts"]
        res["screen_reproduces_repchar1"] = h == res["screen"]["counts"]
        done.append(res)
        per_rep.append((cpu() - c0) / char_spec["repetitions"])
        write(final=False)
    doc = write(final=not stopped)
    REPORT.write_text(report(doc))
    print(json.dumps({"status": doc["status"], "cells_done": doc["cells_done"], "resources": doc["resources"],
                      "screen_reproduces_repchar1_all": all(r["screen_reproduces_repchar1"] for r in done)}, indent=1))
    return 0


def desc(c):
    return ", ".join(f"{k}={c[k]}" for k in ("channel", "delta", "rho", "w_b", "D", "k", "phi") if k in c)


ROLE_TITLES = {"primary_null": "Primary identical-distribution null (S0)",
               "heterogeneous_exchangeable_null": "Heterogeneous exchangeable null (S4b, Δ = 0)",
               "single_replica_perturbation_not_ordinary_null": "Single-replica perturbation, not an ordinary null (S4a Δ = 0, S4c)",
               "condition_confounded_not_null": "Condition-confounded sampling, not a null (S5)",
               "condition_effect": "Condition effect (S1, S2, S3, S4a Δ = 1, S4b Δ = 1)"}


def report(doc):
    L = ["# SYNTHETIC_PERM-1: complete-separation screen vs replica permutation (SYNTHETIC)", "",
         f"Generated from `SYNTHETIC_PERM-1_results.json`. Status **{doc['status']}** ({doc['cells_done']} of "
         f"{doc['cells_planned']} cells). α = {doc['alpha']}; nominal permutation size 12/252 ≈ 0.048 under exact "
         "exchangeability. 200 repetitions per cell, which gives coarse Monte Carlo precision (Wilson 95 % intervals). "
         "Both methods use identical draws. Screen categories: C consistent, RD replica-dependent, I inconclusive. "
         "Per-evaluation SYNTHETIC results; not a multiple-window procedure. No real data.", "",
         f"Screen counts reproduce REPCHAR-1 in all cells: "
         f"{all(r['screen_reproduces_repchar1'] for r in doc['cells'])}. Resources: {doc['resources']['cpu_s']:.1f} CPU-s "
         f"(cap {doc['resources']['cap_cpu_s']}), wall {doc['resources']['wall_s']:.1f} s.", ""]
    for role, title in ROLE_TITLES.items():
        rows = [r for r in doc["cells"] if r["role"] == role]
        if not rows:
            continue
        L += [f"## {title}", "",
              "| Scenario | Cell | Screen C | Screen RD | Screen I | Permutation reject [95 % CI] | reject & C / RD / I | accept & C / RD |",
              "|---|---|---|---|---|---|---|---|"]
        for r in rows:
            s, p, a = r["screen"]["frequency"], r["permutation"], r["agreement"]
            L.append(f"| {r['cell']['scenario']} | {desc(r['cell'])} | {s['consistent']:.3f} | {s['replica_dependent']:.3f} | "
                     f"{s['inconclusive']:.3f} | {p['frequency']:.3f} [{p['wilson95'][0]:.3f}, {p['wilson95'][1]:.3f}] | "
                     f"{a['reject&consistent']} / {a['reject&replica_dependent']} / {a['reject&inconclusive']} | "
                     f"{a['accept&consistent']} / {a['accept&replica_dependent']} |")
        L.append("")
    return "\n".join(L)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("freeze", "run"):
        print(__doc__)
        sys.exit(2)
    sys.exit(freeze() if sys.argv[1] == "freeze" else run())
