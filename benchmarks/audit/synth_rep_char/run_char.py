"""SYNTHETIC_REPCHAR-1: sensitivity characterization of the replica-comparison rule on SYNTHETIC data.

    .venv/bin/python benchmarks/audit/synth_rep_char/run_char.py freeze
    .venv/bin/python benchmarks/audit/synth_rep_char/run_char.py run

Spec: docs/REPLICA_COMPARISON_CHARACTERIZATION_SPEC.md and run_spec.json (cell grid). No real data is read.
"""

import datetime as dt
import hashlib
import json
import math
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(HERE))
from synthetic_generators import circular_series, linear_series, mixture_series  # noqa: E402
from topofold_inputs.compare import circular_hist, linear_hist, replica_comparison  # noqa: E402

SPEC_JSON = HERE / "run_spec.json"
SPEC_DOC = ROOT / "docs" / "REPLICA_COMPARISON_CHARACTERIZATION_SPEC.md"
FREEZE = HERE / "results" / "SYNTHETIC_REPCHAR-1_freeze.json"
RESULTS = HERE / "results" / "SYNTHETIC_REPCHAR-1_results.json"
REPORT = HERE / "results" / "SYNTHETIC_REPCHAR-1_report.md"
CODE = ["python/topofold_inputs/compare.py", "benchmarks/audit/synth_rep_char/synthetic_generators.py",
        "benchmarks/audit/synth_rep_char/run_char.py", "tests/test_synthetic_generators.py"]
CATS = ("consistent", "replica_dependent", "inconclusive")


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def freeze():
    if FREEZE.exists():
        print("freeze exists")
        return 1
    rec = {"id": "SYNTHETIC_REPCHAR-1", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "spec_json_sha256": sha256(SPEC_JSON), "spec_doc_sha256": sha256(SPEC_DOC),
           "code_sha256": {c: sha256(ROOT / c) for c in CODE}}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=1))
    return 0


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return [None, None]
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return [max(0.0, c - h), min(1.0, c + h)]


def simulate(cell, spec, rep):
    """One repetition: returns (hists A, hists B, extra)."""
    fx = spec["fixed"]
    n = cell.get("n_rep", fx["n_rep"])
    T = fx["n_frames"]
    rng = np.random.default_rng(np.random.SeedSequence([spec["seed"]["master"], cell["cell"], rep]))
    sc, phi = cell["scenario"], cell.get("phi", 0.0)
    extra = {}
    if cell["channel"] == "circular":
        mu, s = fx["circular"]["mu"], fx["circular"]["s"]
        hist = lambda x: circular_hist(x, fx["bins"], fx["circular"]["origin"])  # noqa: E731
        d = cell.get("delta", 0.0)
        A = [circular_series(rng, T, mu, s, phi) for _ in range(n)]
        B = [circular_series(rng, T, mu + d * s, s, phi) for _ in range(n)]
        return [hist(x) for x in A], [hist(x) for x in B], extra
    mu, sig = fx["linear"]["mu"], fx["linear"]["sigma"]
    hist = lambda x: linear_hist(x, fx["linear"]["lo"], fx["linear"]["hi"], fx["bins"])[0]  # noqa: E731
    mA = [mu] * n
    mB = [mu] * n
    sA = [sig] * n
    sB = [sig] * n
    if sc == "S5":
        off_a, off_b = rng.normal(0, cell["k"] * sig), rng.normal(0, cell["k"] * sig)
        mA, mB = [mu + off_a] * n, [mu + off_b] * n
        extra["realized_offset_diff_sigma"] = abs(off_a - off_b) / sig
    elif sc in ("S1", "S6"):
        mB = [mu + cell["delta"] * sig] * n
    elif sc == "S2":
        sB = [sig * cell["rho"]] * n
    elif sc == "S4a":
        mA = [mu + 2 * sig] + [mu] * (n - 1)
        mB = [mu + cell["delta"] * sig] * n
    elif sc == "S4b":
        oa = rng.normal(0, cell["k"] * sig, n)
        ob = rng.normal(0, cell["k"] * sig, n)
        mA = list(mu + oa)
        mB = list(mu + cell["delta"] * sig + ob)
    if sc == "S3":
        A = [mixture_series(rng, T, mu - 2 * sig, mu + 2 * sig, sig, 0.5, cell["D"])[0] for _ in range(n)]
        B = [mixture_series(rng, T, mu - 2 * sig, mu + 2 * sig, sig, cell["w_b"], cell["D"])[0] for _ in range(n)]
    elif sc == "S4c":
        A = [mixture_series(rng, T, mu, mu + 4 * sig, sig, 0.1, cell["D"])[0]] + \
            [linear_series(rng, T, mu, sig, 0.0) for _ in range(n - 1)]
        B = [linear_series(rng, T, mu, sig, 0.0) for _ in range(n)]
    else:
        A = [linear_series(rng, T, mA[i], sA[i], phi) for i in range(n)]
        B = [linear_series(rng, T, mB[i], sB[i], phi) for i in range(n)]
    return [hist(x) for x in A], [hist(x) for x in B], extra


def run_cell(cell, spec):
    R = spec["repetitions"]
    cats, frac, ratio, pooled, s5 = [], [], [], [], {"call": [], "no_call": []}
    for rep in range(R):
        ha, hb, extra = simulate(cell, spec, rep)
        rc = replica_comparison(ha, hb)
        cats.append(rc["category"])
        frac.append(rc["fraction_between_above_max_within"])
        ratio.append(rc["median_between"] / rc["median_within"] if rc["median_within"] > 0 else float("inf"))
        pooled.append(rc["pooled"])
        if "realized_offset_diff_sigma" in extra:
            s5["call" if rc["category"] != "inconclusive" else "no_call"].append(extra["realized_offset_diff_sigma"])
    counts = {c: cats.count(c) for c in CATS}
    out = {"cell": cell, "repetitions": R,
           "counts": counts,
           "frequency": {c: counts[c] / R for c in CATS},
           "wilson95": {c: wilson(counts[c], R) for c in CATS},
           "call": {"frequency": (counts["consistent"] + counts["replica_dependent"]) / R,
                    "wilson95": wilson(counts["consistent"] + counts["replica_dependent"], R),
                    "breakdown": {"consistent": counts["consistent"], "replica_dependent": counts["replica_dependent"]}},
           "fraction_between_above_max_within": {"median": float(np.median(frac)),
                                                 "iqr": [float(np.percentile(frac, 25)), float(np.percentile(frac, 75))]},
           "median_ratio_between_within": float(np.median(ratio)), "median_pooled": float(np.median(pooled))}
    if cell["scenario"] == "S5":
        out["s5_realized_offset_diff_sigma"] = {k: (float(np.median(v)) if v else None) for k, v in s5.items()}
        out["s5_counts_by_call"] = {k: len(v) for k, v in s5.items()}
    return out


def run():
    spec = json.loads(SPEC_JSON.read_text())
    frz = json.loads(FREEZE.read_text())
    bad = [c for c, h in frz["code_sha256"].items() if sha256(ROOT / c) != h]
    if frz["spec_json_sha256"] != sha256(SPEC_JSON) or frz["spec_doc_sha256"] != sha256(SPEC_DOC) or bad:
        raise SystemExit(f"STOP: spec or code differs from the freeze: {bad}")
    cap = spec["budget"]["cpu_seconds"]
    resource.setrlimit(resource.RLIMIT_CPU, (cap, cap))
    t_wall = time.time()
    done, per_rep = [], []
    status, not_run = "complete", []

    def cpu():
        r = resource.getrusage(resource.RUSAGE_SELF)
        return r.ru_utime + r.ru_stime

    def write():
        doc = {"id": "SYNTHETIC_REPCHAR-1", "label": "SYNTHETIC", "status": status, "cells_not_run": not_run,
               "freeze_utc": frz["frozen_utc"], "spec_json_sha256": frz["spec_json_sha256"],
               "resources": {"cpu_s": cpu(), "wall_s": time.time() - t_wall, "cap_cpu_s": cap,
                             "mean_cpu_s_per_repetition": float(np.mean(per_rep)) if per_rep else None},
               "cells": done}
        RESULTS.write_text(json.dumps(doc, indent=1) + "\n")
        return doc

    for cell in spec["cells"]:
        scale = cell.get("n_rep", spec["fixed"]["n_rep"]) / spec["fixed"]["n_rep"]
        if per_rep:
            est = 1.5 * float(np.mean(per_rep)) * spec["repetitions"] * scale
            if cpu() + est > cap:
                status = "partial"
                not_run = [c["cell"] for c in spec["cells"] if c["cell"] >= cell["cell"]]
                break
        c0 = cpu()
        done.append(run_cell(cell, spec))
        per_rep.append((cpu() - c0) / spec["repetitions"] / scale)
        write()
    doc = write()
    REPORT.write_text(report(doc, spec))
    print(json.dumps({"status": status, "cells_done": len(done), "resources": doc["resources"]}, indent=1))
    return 0


def describe(cell):
    keys = [k for k in ("channel", "delta", "rho", "w_b", "D", "k", "n_rep", "phi") if k in cell]
    return ", ".join(f"{k}={cell[k]}" for k in keys)


def report(doc, spec):
    L = ["# SYNTHETIC_REPCHAR-1: sensitivity of the replica-comparison rule (SYNTHETIC data)", "",
         f"Generated from `SYNTHETIC_REPCHAR-1_results.json`. Status: **{doc['status']}**"
         + (f" (cells not run: {doc['cells_not_run']})" if doc["cells_not_run"] else "") + ". "
         f"{spec['repetitions']} repetitions per cell, which gives coarse Monte Carlo precision (Wilson 95 % "
         "intervals shown). Categories: C = consistent, RD = replica-dependent, I = inconclusive; "
         "call = C + RD, always with its breakdown. No real data were used; DEV-REP-1 is not re-scored.", "",
         f"Resources: {doc['resources']['cpu_s']:.1f} CPU-s (cap {doc['resources']['cap_cpu_s']}), wall "
         f"{doc['resources']['wall_s']:.1f} s.", ""]
    by = {}
    for r in doc["cells"]:
        by.setdefault(r["cell"]["scenario"], []).append(r)
    titles = {"S0": "S0 identical-distribution null", "S5": "S5 condition-confounded sampling",
              "S1": "S1 location shift (units of σ or s)", "S2": "S2 scale change", "S3": "S3 population shift between fixed modes",
              "S4a": "S4a one outlier replica", "S4b": "S4b random replica offsets", "S4c": "S4c rare excursion in one replica",
              "S6": "S6 replica count (Δ = 1, φ = 0.8)"}
    for sc in ("S0", "S5", "S1", "S2", "S3", "S4a", "S4b", "S4c", "S6"):
        rows = by.get(sc)
        if not rows:
            continue
        L += [f"## {titles[sc]}", "", "| Cell | C | RD | I | call [95 % CI] | median fraction (IQR) | median between/within |"
              + (" median realized offset diff, calls / no calls (σ) |" if sc == "S5" else ""),
              "|---|---|---|---|---|---|---|" + ("---|" if sc == "S5" else "")]
        for r in rows:
            f, ci = r["frequency"], r["call"]["wilson95"]
            line = (f"| {describe(r['cell'])} | {f['consistent']:.3f} | {f['replica_dependent']:.3f} | "
                    f"{f['inconclusive']:.3f} | {r['call']['frequency']:.3f} [{ci[0]:.3f}, {ci[1]:.3f}] | "
                    f"{r['fraction_between_above_max_within']['median']:.2f} "
                    f"({r['fraction_between_above_max_within']['iqr'][0]:.2f}–{r['fraction_between_above_max_within']['iqr'][1]:.2f}) | "
                    f"{r['median_ratio_between_within']:.2f} |")
            if sc == "S5":
                v = r["s5_realized_offset_diff_sigma"]
                line += f" {v['call'] if v['call'] is not None else '–'} / {v['no_call'] if v['no_call'] is not None else '–'} |"
            L.append(line)
        L.append("")
    L += ["## Reading aid (not a criterion): smallest effect with call frequency ≥ 0.5", ""]
    for sc, key in (("S1", "delta"), ("S2", "rho"), ("S3", "w_b")):
        rows = by.get(sc, [])
        groups = {}
        for r in rows:
            g = tuple((k, r["cell"][k]) for k in ("channel", "phi", "D") if k in r["cell"])
            groups.setdefault(g, []).append(r)
        for g, rs in groups.items():
            # smallest effect first: larger delta / rho, or w_b further from 0.5, means a larger effect
            order = sorted(rs, key=lambda r: abs(0.5 - r["cell"][key]) if key == "w_b" else r["cell"][key])
            hit = next((r for r in order if r["call"]["frequency"] >= 0.5), None)
            lab = ", ".join(f"{k}={v}" for k, v in g)
            L.append(f"- {sc} [{lab}]: " + (f"{key}={hit['cell'][key]} (C {hit['counts']['consistent']}, "
                                             f"RD {hit['counts']['replica_dependent']} of {hit['repetitions']})"
                                             if hit else "not reached in the grid"))
    L.append("")
    return "\n".join(L)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("freeze", "run"):
        print(__doc__)
        sys.exit(2)
    sys.exit(freeze() if sys.argv[1] == "freeze" else run())
