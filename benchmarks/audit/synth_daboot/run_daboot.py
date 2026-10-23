"""SYNTHETIC_DABOOT-1: stationary-bootstrap Kuiper test of equal marginals under unequal temporal correlation.

    .venv/bin/python benchmarks/audit/synth_daboot/run_daboot.py freeze
    .venv/bin/python benchmarks/audit/synth_daboot/run_daboot.py run

Spec: docs/DEPENDENCE_AWARE_MARGINAL_COMPARISON_SPEC.md and run_spec.json. SYNTHETIC only.
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
sys.path.insert(0, str(ROOT / "benchmarks/audit/synth_null_stress"))
import run_stress  # noqa: E402  (frozen NULLSTRESS-1 runner; only draw() is used, for the identity check)
from synthetic_generators import circular_series, linear_series  # noqa: E402
from topofold_inputs.dependent_bootstrap import two_sample_test  # noqa: E402

SPEC_JSON = HERE / "run_spec.json"
SPEC_DOC = ROOT / "docs/DEPENDENCE_AWARE_MARGINAL_COMPARISON_SPEC.md"
NS_SPEC = ROOT / "benchmarks/audit/synth_null_stress/run_spec.json"
NS_RESULTS = ROOT / "benchmarks/audit/synth_null_stress/results/SYNTHETIC_NULLSTRESS-1_results.json"
FREEZE = HERE / "results/SYNTHETIC_DABOOT-1_freeze.json"
RESULTS = HERE / "results/SYNTHETIC_DABOOT-1_results.json"
REPORT = HERE / "results/SYNTHETIC_DABOOT-1_report.md"
CODE = ["python/topofold_inputs/dependent_bootstrap.py", "benchmarks/audit/synth_daboot/run_daboot.py",
        "benchmarks/audit/synth_rep_char/synthetic_generators.py", "benchmarks/audit/synth_null_stress/run_stress.py",
        "tests/test_dependent_bootstrap.py"]


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
    rec = {"id": "SYNTHETIC_DABOOT-1", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "spec_json_sha256": sha256(SPEC_JSON), "spec_doc_sha256": sha256(SPEC_DOC),
           "nullstress_spec_sha256": sha256(NS_SPEC), "code_sha256": {c: sha256(ROOT / c) for c in CODE},
           "environment": {"python": platform.python_version(), "numpy": np.__version__, "platform": platform.platform()}}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=1))
    return 0


def raw_draw(spec, cell, rep):
    """Raw series in the NULLSTRESS-1 generator order (A0..A4 then B0..B4); location shift only for positive controls."""
    rng = np.random.default_rng(np.random.SeedSequence([spec["data_seed_master"], cell["cell"], rep]))
    n, T = spec["n_rep"], spec["n_frames"]
    if cell["channel"] == "linear":
        mu, sc, gen = spec["linear"]["mu"], spec["linear"]["sigma"], linear_series
    else:
        mu, sc, gen = spec["circular"]["mu"], spec["circular"]["s"], circular_series
    A = [gen(rng, T, mu, sc, cell["phi_a"]) for _ in range(n)]
    B = [gen(rng, T, mu + cell["delta"] * sc, sc, cell["phi_b"]) for _ in range(n)]
    return A, B


def grid_for(spec, channel):
    if channel == "linear":
        lo, hi, G = spec["linear"]["grid"]
        return np.linspace(lo, hi, G)
    return np.linspace(-math.pi, math.pi, 257)[1:]


def identity_check(spec):
    """Repetition 0 of every null cell: histograms of our raw draws equal the frozen NULLSTRESS-1 draws."""
    ns = json.loads(NS_SPEC.read_text())
    ok = {}
    for cell in spec["cells"]:
        if cell["role"] != "null":
            continue
        A, B = raw_draw(spec, cell, 0)
        hA, hB = run_stress.draw(ns, ns["cells"][cell["cell"]], 0)
        if cell["channel"] == "linear":
            h = lambda x: np.histogram(x[(x >= 3) & (x <= 8)], bins=36, range=(3, 8))[0]  # noqa: E731
            ok[cell["cell"]] = all(np.array_equal(h(a)[:36], ha[:36]) for a, ha in zip(A + B, hA + hB))
        else:
            from topofold_inputs.compare import circular_hist
            ok[cell["cell"]] = all(np.array_equal(circular_hist(a, 36, -math.pi), ha) for a, ha in zip(A + B, hA + hB))
    return ok


def run():
    spec = json.loads(SPEC_JSON.read_text())
    frz = json.loads(FREEZE.read_text())
    bad = [c for c, v in frz["code_sha256"].items() if sha256(ROOT / c) != v]
    if frz["spec_json_sha256"] != sha256(SPEC_JSON) or frz["spec_doc_sha256"] != sha256(SPEC_DOC) or bad:
        raise SystemExit(f"STOP: spec or code differs from the freeze: {bad}")
    cap = spec["budget"]["cpu_seconds"]
    resource.setrlimit(resource.RLIMIT_CPU, (cap, cap))
    t0 = time.time()
    ident = identity_check(spec)
    if not all(ident.values()):
        raise SystemExit(f"STOP: null draws differ from NULLSTRESS-1: {ident}")
    ns_perm = {r["cell"]["cell"]: r["permutation"] for r in json.loads(NS_RESULTS.read_text())["cells"]}
    done, per_rep, not_run, stopped = [], [], [], False
    R = spec["repetitions"]

    def cpu():
        r = resource.getrusage(resource.RUSAGE_SELF)
        return r.ru_utime + r.ru_stime

    def write(final):
        doc = {"id": "SYNTHETIC_DABOOT-1", "label": "SYNTHETIC",
               "status": "complete" if final and len(done) == len(spec["cells"]) else "partial",
               "cells_planned": len(spec["cells"]), "cells_done": len(done), "cells_not_run": not_run,
               "freeze_utc": frz["frozen_utc"], "alpha": spec["alpha"], "B": spec["B"],
               "null_draw_identity_check": ident, "environment": frz["environment"],
               "resources": {"cpu_s": cpu(), "wall_s": time.time() - t0, "cap_cpu_s": cap}, "cells": done}
        RESULTS.write_text(json.dumps(doc, indent=1) + "\n")
        return doc

    for cell in spec["cells"]:
        if per_rep and cpu() + 1.5 * float(np.mean(per_rep)) * R > cap:
            not_run = [c["cell"] for c in spec["cells"] if c["cell"] >= cell["cell"]]
            stopped = True
            break
        c0 = cpu()
        grid = grid_for(spec, cell["channel"])
        rej, ps, la, lb = 0, [], [], []
        for rep in range(R):
            A, B = raw_draw(spec, cell, rep)
            brng = np.random.default_rng(np.random.SeedSequence([spec["bootstrap_seed_master"], cell["cell"], rep]))
            r = two_sample_test(A, B, grid, spec["B"], brng, circular=cell["channel"] == "circular")
            rej += r["p"] <= spec["alpha"]
            ps.append(r["p"])
            la.append(r["block_length_a"])
            lb.append(r["block_length_b"])
        res = {"cell": cell, "repetitions": R, "rejections": int(rej), "frequency": rej / R,
               "wilson95": wilson(int(rej), R), "p_median": float(np.median(ps)),
               "block_length_median": {"A": float(np.median(la)), "B": float(np.median(lb))}}
        if cell["cell"] in ns_perm:
            res["nullstress_permutation_frequency_same_draws"] = ns_perm[cell["cell"]]["frequency"]
        done.append(res)
        per_rep.append((cpu() - c0) / R)
        write(final=False)
    doc = write(final=not stopped)
    doc["assessment"] = assess(doc, spec)
    RESULTS.write_text(json.dumps(doc, indent=1) + "\n")
    REPORT.write_text(report(doc))
    print(json.dumps({"status": doc["status"], "cells_done": doc["cells_done"], "resources": doc["resources"],
                      "assessment": doc["assessment"]}, indent=1))
    return 0


def assess(doc, spec):
    nulls = [r for r in doc["cells"] if r["cell"]["role"] == "null"]
    pos = [r for r in doc["cells"] if r["cell"]["role"] == "positive_control"]
    fail_cells = [r["cell"]["cell"] for r in nulls if r["frequency"] > 0.10 or r["wilson95"][0] > 0.05]
    conservative = [r["cell"]["cell"] for r in nulls if r["frequency"] <= 0.01]
    d2 = [r for r in pos if r["cell"]["delta"] == 2.0]
    complete = doc["status"] == "complete"
    calib = "failure" if fail_cells else ("acceptable (practical)" if len(nulls) == 6 else "incomplete")
    sens = None
    if len(d2) == 2:
        sens = "useful" if all(r["frequency"] >= 0.80 for r in d2) else "low"
    if calib == "failure":
        decision = "calibration failure: retain descriptive comparisons; record limitation; stop"
    elif not complete or sens is None:
        decision = "incomplete run: no decision"
    elif sens == "useful":
        decision = "recommend a separately frozen real-data development evaluation"
    else:
        decision = "unresolved limitation: insufficient sensitivity"
    return {"calibration": calib, "failing_null_cells": fail_cells, "conservative_null_cells": conservative,
            "sensitivity": sens, "decision_by_predeclared_rule": decision}


def report(doc):
    a = doc["assessment"]
    L = ["# SYNTHETIC_DABOOT-1: stationary-bootstrap Kuiper test under unequal temporal correlation (SYNTHETIC)", "",
         f"Generated from `SYNTHETIC_DABOOT-1_results.json`. Status **{doc['status']}** ({doc['cells_done']} of "
         f"{doc['cells_planned']} cells). α = {doc['alpha']}, B = {doc['B']} bootstrap draws per test, 200 repetitions per "
         "cell, which gives coarse Monte Carlo precision. Null draws are identical to SYNTHETIC_NULLSTRESS-1 (checked: "
         f"{all(doc['null_draw_identity_check'].values())}). No real data.", "",
         f"Resources: {doc['resources']['cpu_s']:.1f} CPU-s (cap {doc['resources']['cap_cpu_s']}), wall "
         f"{doc['resources']['wall_s']:.1f} s.", "",
         "| Cell | Role | Channel | φ_A | φ_B | Δ | Reject [95 % CI] | median block length A / B | NULLSTRESS-1 permutation (same draws) |",
         "|---|---|---|---|---|---|---|---|---|"]
    for r in doc["cells"]:
        c = r["cell"]
        perm = r.get("nullstress_permutation_frequency_same_draws")
        L.append(f"| {c['cell']} | {c['role']} | {c['channel']} | {c['phi_a']} | {c['phi_b']} | {c['delta']} | "
                 f"{r['frequency']:.3f} [{r['wilson95'][0]:.3f}, {r['wilson95'][1]:.3f}] | "
                 f"{r['block_length_median']['A']:.1f} / {r['block_length_median']['B']:.1f} | "
                 f"{'' if perm is None else f'{perm:.3f}'} |")
    L += ["", "## Predeclared assessment", "",
          f"- Calibration: **{a['calibration']}** (failing null cells: {a['failing_null_cells'] or 'none'}; "
          f"conservative cells ≤ 0.01: {a['conservative_null_cells'] or 'none'}). Calibration failure = a null cell above "
          "0.10 or with a Wilson lower bound above 0.05. An interval containing 0.05 is not proof of validity.",
          f"- Sensitivity (both Δ = 2 controls ≥ 0.80): **{a['sensitivity']}**.",
          f"- Decision by the predeclared rule: **{a['decision_by_predeclared_rule']}**.", ""]
    return "\n".join(L)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("freeze", "run"):
        print(__doc__)
        sys.exit(2)
    sys.exit(freeze() if sys.argv[1] == "freeze" else run())
