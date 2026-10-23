"""KRAS-STRUCT-1: descriptive comparison of GDP-bound (9O0N) and GMPPNP-bound (9O0O) KRAS crystal structures.

    .venv/bin/python benchmarks/audit/kras_struct_cmp/run_kras.py freeze
    .venv/bin/python benchmarks/audit/kras_struct_cmp/run_kras.py run

Static structures: direct geometric differences only (see run_spec.json, permitted_interpretation).
"""

import csv
import datetime as dt
import hashlib
import json
import math
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs.compare import spearman, window_distances  # noqa: E402
from topofold_inputs.features import export_features, features_from_export  # noqa: E402

ADAPTER_PY = ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python"
SPEC = HERE / "run_spec.json"
FREEZE = HERE / "results/KRAS-STRUCT-1_freeze.json"
CODE = ["benchmarks/audit/kras_struct_cmp/run_kras.py", "benchmarks/audit/kras_struct_cmp/rmsd_windows.py",
        "python/topofold_inputs/ca_adapter.py", "python/topofold_inputs/features.py",
        "python/topofold_inputs/compare.py"]
METRICS = ("d_kappa", "d_i_i2", "rmsd3", "d_tau", "d_i_i3", "rmsd4")


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def freeze():
    if FREEZE.exists():
        print("freeze exists")
        return 1
    rec = {"id": "KRAS-STRUCT-1", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "spec_sha256": sha256(SPEC), "code_sha256": {c: sha256(ROOT / c) for c in CODE}}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=1))
    return 0


def adapter(pdb, run_id, build, corr, altloc):
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"))
    p = subprocess.run([str(ADAPTER_PY), "-m", "topofold_inputs", str(ROOT / pdb), "--run-id", run_id, "--out",
                        str(build), "--correspondence", str(corr), "--on-missing", "reject", "--altloc", altloc],
                       env=env, capture_output=True, text=True)
    return p.returncode, p.stderr.strip()[-500:]


def extract(spec, build, corr, altloc, tag):
    out = {}
    for state in ("GDP", "GNP"):
        run_id = f"{tag}_{state}"
        rc, err = adapter(spec["structures"][state]["file"], run_id, build, corr, altloc)
        if rc != 0:
            return None, f"{state}: adapter rejected ({err})"
        feats, rep = features_from_export(build / f"{run_id}.npz", build / f"{run_id}.json")
        export_features(feats, build / f"{run_id}.npz", build / f"{run_id}.json", build, run_id)
        out[state] = {"rep": rep, "feats": feats, "coords": np.load(build / f"{run_id}.npz")["coordinates_angstrom"],
                      "npz": build / f"{run_id}.npz"}
    return out, "ok"


def metrics(ex, build, tag):
    a, b = ex["GDP"], ex["GNP"]
    subprocess.run([str(ADAPTER_PY), str(HERE / "rmsd_windows.py"), str(a["npz"]), str(b["npz"]),
                    str(build / f"{tag}_rmsd.json")], check=True)
    rm = json.loads((build / f"{tag}_rmsd.json").read_text())
    ka, kb = a["feats"].kappa[0], b["feats"].kappa[0]
    ta, tb = a["feats"].tau[0], b["feats"].tau[0]
    kva, kvb = a["feats"].kappa_valid[0], b["feats"].kappa_valid[0]
    tva, tvb = a["feats"].tau_valid[0], b["feats"].tau_valid[0]
    d2a, d2b = window_distances(a["coords"], 2)[0], window_distances(b["coords"], 2)[0]
    d3a, d3b = window_distances(a["coords"], 3)[0], window_distances(b["coords"], 3)[0]
    n = len(ka)
    m = {k: [None] * n for k in METRICS}
    for i in range(n):
        if kva[i] and kvb[i]:
            m["d_kappa"][i] = abs(float(kb[i] - ka[i]))
            m["d_i_i2"][i] = abs(float(d2b[i] - d2a[i]))
            m["rmsd3"][i] = rm["rmsd3"][i]
        if i < len(ta) and tva[i] and tvb[i]:
            m["d_tau"][i] = abs(math.atan2(math.sin(tb[i] - ta[i]), math.cos(tb[i] - ta[i])))
            m["d_i_i3"][i] = abs(float(d3b[i] - d3a[i]))
            m["rmsd4"][i] = rm["rmsd4"][i]
    return m, rm["global"], rm["n_global"]


def run():
    spec = json.loads(SPEC.read_text())
    frz = json.loads(FREEZE.read_text())
    bad = [c for c, v in frz["code_sha256"].items() if sha256(ROOT / c) != v]
    if frz["spec_sha256"] != sha256(SPEC) or bad:
        raise SystemExit(f"STOP: spec or code differs from the freeze: {bad}")
    for st in spec["structures"].values():
        if sha256(ROOT / st["file"]) != st["sha256"]:
            raise SystemExit(f"STOP: {st['pdb']} hash differs")
    t0 = time.time()
    build = HERE / "build"
    build.mkdir(exist_ok=True)
    corr = build / "correspondence.json"
    corr.write_text(json.dumps([{"chain": "A", "resid": r, "label": f"KRAS:{r}"} for r in range(1, 169)]))
    ex, why = extract(spec, build, corr, "blank_then_A", "primary")
    if ex is None:
        raise SystemExit(f"STOP: primary extraction failed: {why}")
    ra, rb = ex["GDP"]["rep"]["residues"], ex["GNP"]["rep"]["residues"]
    checks = {"residue_names_identical": [r["resname"] for r in ra] == [r["resname"] for r in rb],
              "all_present": all(r["status"] == "present" for r in ra + rb),
              "n_residues": len(ra), "breaks_GDP": ex["GDP"]["rep"]["breaks"], "breaks_GNP": ex["GNP"]["rep"]["breaks"],
              "altloc_decisions_GDP": ex["GDP"]["rep"]["altloc_decisions"],
              "altloc_decisions_GNP": ex["GNP"]["rep"]["altloc_decisions"]}
    if not (checks["residue_names_identical"] and checks["all_present"]):
        raise SystemExit(f"STOP: mapping check failed: {checks}")
    m, g_rmsd, n_g = metrics(ex, build, "primary")
    alt_res = {d["resid"] for d in checks["altloc_decisions_GDP"] + checks["altloc_decisions_GNP"]}
    resid = [r["resid"] for r in ra]

    def span(i, four):
        return resid[i - 1: i + (3 if four else 2)] if i >= 1 else []

    rows = []
    for i in range(len(resid)):
        row = {"index": i, "resid": resid[i], "resname": ra[i]["resname"],
               "window3": "-".join(map(str, span(i, False)[::2])) if 1 <= i <= len(resid) - 2 else "",
               "window4": f"{span(i, True)[0]}-{span(i, True)[-1]}" if 1 <= i <= len(resid) - 3 else "",
               "altloc_in_window3": bool(set(span(i, False)) & alt_res),
               "altloc_in_window4": bool(set(span(i, True)) & alt_res)}
        row.update({k: m[k][i] for k in METRICS})
        rows.append(row)
    # altloc sensitivity
    ex2, why2 = extract(spec, build, corr, "highest_occupancy", "sens")
    sens = {"status": why2}
    if ex2 is not None:
        m2, g2, _ = metrics(ex2, build, "sens")
        sens.update(max_abs_change={k: max((abs(m2[k][i] - m[k][i]) for i in range(len(resid))
                                            if m[k][i] is not None and m2[k][i] is not None), default=None)
                                    for k in METRICS}, global_rmsd=g2)
    top = {k: sorted(({"resid": r["resid"], "window": r["window4" if k in ("d_tau", "d_i_i3", "rmsd4") else "window3"],
                       "value": r[k], "altloc": r["altloc_in_window4" if k in ("d_tau", "d_i_i3", "rmsd4") else "altloc_in_window3"]}
                      for r in rows if r[k] is not None), key=lambda x: -x["value"])[:10] for k in METRICS}

    def rho(x, y):
        pairs = [(r[x], r[y]) for r in rows if r[x] is not None and r[y] is not None]
        return spearman([p[0] for p in pairs], [p[1] for p in pairs]), len(pairs)

    corr_tab = {f"{x}~{y}": dict(zip(("spearman", "n"), rho(x, y))) for x, y in
                [("d_kappa", "d_i_i2"), ("d_kappa", "rmsd3"), ("d_i_i2", "rmsd3"),
                 ("d_tau", "d_i_i3"), ("d_tau", "rmsd4"), ("d_i_i3", "rmsd4")]}
    ru = resource.getrusage(resource.RUSAGE_SELF)
    rc = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = ru.ru_utime + ru.ru_stime + rc.ru_utime + rc.ru_stime
    out = {"id": "KRAS-STRUCT-1", "freeze_utc": frz["frozen_utc"], "spec_sha256": frz["spec_sha256"],
           "checks": checks, "global_rmsd_A": g_rmsd, "global_rmsd_n_ca": n_g,
           "authors_reported_rmsd_A": 0.3, "top10": top, "rank_agreement": corr_tab, "altloc_sensitivity": sens,
           "resources": {"cpu_s": cpu, "wall_s": time.time() - t0},
           "core": json.loads((build / "primary_GDP_features.json").read_text())["provenance"]["core"]}
    (HERE / "results/KRAS-STRUCT-1_results.json").write_text(json.dumps(out, indent=1) + "\n")
    with open(HERE / "results/KRAS-STRUCT-1_table.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else ("" if v is None else v)) for k, v in r.items()})
    (HERE / "results/KRAS-STRUCT-1_report.md").write_text(report(out, spec))
    print(json.dumps({k: out[k] for k in ("checks", "global_rmsd_A", "rank_agreement", "altloc_sensitivity", "resources")},
                     indent=1, default=str))
    return 0


def report(o, spec):
    L = ["# KRAS-STRUCT-1: GDP-bound (9O0N) vs GMPPNP-bound (9O0O) wild-type KRAS, both with MRTX1133", "",
         "Generated from `KRAS-STRUCT-1_results.json` and `KRAS-STRUCT-1_table.csv`. **Static structures: direct "
         "geometric differences only.** Confounded by crystal form (P 21 21 21 vs P 2 3) and crystallization "
         "conditions, resolution (1.40 vs 1.90 Å), alternate conformations and coordinate error; both contain "
         "MRTX1133. No ensemble, population, kinetic or functional-state claim; nucleotide identity does not "
         "establish functional state.", "",
         f"- Mapping: chain A, residues 1–168 ({o['checks']['n_residues']} residues), names identical, all present; "
         f"breaks: GDP {o['checks']['breaks_GDP'] or 'none'}, GMPPNP {o['checks']['breaks_GNP'] or 'none'}.",
         f"- Global Cα RMSD after optimal superposition (MDAnalysis QCP, {o['global_rmsd_n_ca']} Cα): "
         f"{o['global_rmsd_A']:.3f} Å (the authors report 0.3 Å; atom selection may differ).",
         f"- Alternate locations: primary policy blank-then-A; sensitivity (highest occupancy): {o['altloc_sensitivity']['status']}"
         + (f"; max change per metric {json.dumps({k: (None if v is None else round(v, 4)) for k, v in o['altloc_sensitivity']['max_abs_change'].items()})}"
            if 'max_abs_change' in o['altloc_sensitivity'] else ""), "",
         "## Largest local differences (descriptive)", "",
         "| Metric | Unit | Ten largest windows (value; * = window contains an altloc residue) |", "|---|---|---|"]
    units = {"d_kappa": "rad", "d_tau": "rad", "d_i_i2": "Å", "d_i_i3": "Å", "rmsd3": "Å", "rmsd4": "Å"}
    for k in METRICS:
        L.append(f"| {k} | {units[k]} | " + ", ".join(f"{t['window']} ({t['value']:.3f}{'*' if t['altloc'] else ''})"
                                                     for t in o["top10"][k]) + " |")
    L += ["", "## Agreement between equal-information metrics (Spearman over windows; descriptive)", "",
          "| Pair | Spearman | Windows |", "|---|---|---|"]
    for k, v in o["rank_agreement"].items():
        L.append(f"| {k} | {v['spearman']:.3f} | {v['n']} |")
    L += ["", "Three-point windows: κ, d(i−1, i+1) and the 3-residue RMSD all depend on the same triangle, so close "
          "agreement is expected. Four-point windows add torsion; d(i−1, i+2) is blind to the sign of τ. Agreement or "
          "disagreement does not show that one representation is better.", "",
          f"Resources: {o['resources']['cpu_s']:.1f} CPU-s, {o['resources']['wall_s']:.1f} s wall.", ""]
    return "\n".join(L)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("freeze", "run"):
        print(__doc__)
        sys.exit(2)
    sys.exit(freeze() if sys.argv[1] == "freeze" else run())
