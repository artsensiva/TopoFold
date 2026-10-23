"""Verification of `python -m topofold_inputs.structpair` against the recorded KRAS-STRUCT-1/-2 results.

    .venv/bin/python benchmarks/audit/structpair_demo/verify.py

Runs the command (unchanged between runs) on several input specifications, each into a fresh directory under
build/<UTC stamp>/, compares with the frozen study tables and results files, and writes results/verification.json plus
compact copies of the demonstration outputs (results/<spec>/, without the work/ intermediates).
"""

import csv
import datetime as dt
import json
import os
import resource
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ADAPTER_PY = ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python"
S1 = ROOT / "benchmarks/audit/kras_struct_cmp/results"
S2 = ROOT / "benchmarks/audit/kras_struct_g12d/results"
METRICS = ("d_kappa", "d_i_i2", "rmsd3", "d_tau", "d_i_i3", "rmsd4")
COL = {"d_kappa": "d_kappa_rad", "d_i_i2": "d_i_i2_angstrom", "rmsd3": "rmsd3_angstrom", "d_tau": "d_tau_rad",
       "d_i_i3": "d_i_i3_angstrom", "rmsd4": "rmsd4_angstrom"}
TOL_REL, TOL_ABS = 1e-5, 1e-6       # recorded tables carry 6 significant figures


def cli(spec, out):
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=str(ADAPTER_PY))
    p = subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(spec), "--out", str(out)],
                       env=env, capture_output=True, text=True)
    return p.returncode, (p.stdout + p.stderr).strip()[-600:]


def rows(path):
    return list(csv.DictReader(open(path)))


def f(x):
    return None if x in ("", None) else float(x)


def close(a, b):
    return a is not None and b is not None and abs(a - b) <= max(TOL_ABS, TOL_REL * abs(b))


def by_resid(win):
    """windows.csv rows keyed by the A residue number of the window centre."""
    return {int(r["A_residue"].split(":")[1]): r for r in win}


def compare_values(new, rec, rec_col=lambda k: k):
    n_cmp, worst, fails = 0, 0.0, []
    for resid, r in rec.items():
        for k in METRICS:
            want, got = f(r[rec_col(k)]), f(new[resid][COL[k]]) if resid in new else None
            if want is None and got is None:
                continue
            n_cmp += 1
            if not close(got, want):
                fails.append({"resid": resid, "metric": k, "recorded": want, "new": got})
            else:
                worst = max(worst, abs(got - want))
    return {"n_values_compared": n_cmp, "n_mismatches": len(fails), "max_abs_difference": worst, "mismatches": fails[:10]}


def main():
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    build = HERE / "build" / stamp
    build.mkdir(parents=True)
    out = {"id": "STRUCTPAIR-VERIFY-1", "created_utc": stamp, "command": "python -m topofold_inputs.structpair",
           "module_sha256": __import__("hashlib").sha256((ROOT / "python/topofold_inputs/structpair.py").read_bytes()).hexdigest(),
           "checks": {}}
    C = out["checks"]

    # C1: KRAS-STRUCT-1 reproduction
    rc, msg = cli(HERE / "specs/kras1_9O0N_vs_9O0O.json", build / "kras1")
    s1res = json.loads((S1 / "KRAS-STRUCT-1_results.json").read_text())
    sm = json.loads((build / "kras1/summary.json").read_text())
    vals = compare_values(by_resid(rows(build / "kras1/windows.csv")),
                          {int(r["resid"]): r for r in rows(S1 / "KRAS-STRUCT-1_table.csv")})
    C["C1_struct1_reproduction"] = {
        "exit_code": rc, **vals,
        "global_rmsd": {"new": sm["global"]["rmsd_qcp_angstrom"], "recorded": s1res["global_rmsd_A"],
                        "abs_difference": abs(sm["global"]["rmsd_qcp_angstrom"] - s1res["global_rmsd_A"]),
                        "numpy_kabsch_abs_difference": abs(sm["global"]["rmsd_numpy_kabsch_angstrom"] - s1res["global_rmsd_A"])},
        "agreement_max_abs_difference": max(abs(sm["agreement"][k]["spearman"] - v["spearman"])
                                            for k, v in s1res["rank_agreement"].items()),
        "altloc_sensitivity_max_abs_difference": max(abs(sm["altloc_sensitivity"]["max_abs_change"][k] - v)
                                                     for k, v in s1res["altloc_sensitivity"]["max_abs_change"].items())}

    # C2: KRAS-STRUCT-2 primary pair
    rc, msg = cli(HERE / "specs/kras2_5US4A_vs_5USJA.json", build / "kras2")
    s2res = json.loads((S2 / "KRAS-STRUCT-2_results.json").read_text())
    s2tab = {int(r["resid"]): r for r in rows(S2 / "KRAS-STRUCT-2_table.csv")}
    sm2 = json.loads((build / "kras2/summary.json").read_text())
    new2 = by_resid(rows(build / "kras2/windows.csv"))
    C["C2_struct2_primary"] = {
        "exit_code": rc, **compare_values(new2, s2tab),
        "global_rmsd_abs_difference": abs(sm2["global"]["rmsd_qcp_angstrom"] - s2res["global_rmsd_A"]["cross"]["A-A"]),
        "altloc_sensitivity_status": sm2["altloc_sensitivity"]["status"],
        "exclusions": sorted({(e["window_A"], e["cause"]) for e in rows(build / "kras2/exclusions.csv")}),
        "recorded_exclusions": s2res["exclusions"]}

    # C3: all 12 chain pairings, different specs (chain only), unchanged command
    spec0 = json.loads((HERE / "specs/kras2_5US4A_vs_5USJA.json").read_text())
    per = {}
    for a in "AB":
        for b in "ABCDEF":
            sp = json.loads(json.dumps(spec0))
            sp["name"] = f"5US4 {a} vs 5USJ {b}"
            sp["structures"]["A"]["chain"], sp["structures"]["B"]["chain"] = a, b
            for k in "AB":
                sp["structures"][k]["file"] = str((HERE / "specs" / sp["structures"][k]["file"]).resolve())
            sp.pop("altloc_sensitivity")
            p = build / f"spec_{a}{b}.json"
            p.write_text(json.dumps(sp))
            rc, msg = cli(p, build / f"pair_{a}{b}")
            if rc != 0:
                per[f"{a}-{b}"] = {"exit_code": rc, "message": msg}
                continue
            per[f"{a}-{b}"] = {"exit_code": rc, "windows": by_resid(rows(build / f"pair_{a}{b}/windows.csv")),
                               "global": json.loads((build / f"pair_{a}{b}/summary.json").read_text())["global"]}
    agg_n, agg_fail = 0, []
    for resid, r in s2tab.items():
        for k in METRICS:
            xs = [f(v["windows"][resid][COL[k]]) for v in per.values() if "windows" in v]
            xs = [x for x in xs if x is not None]
            for stat, val in (("cross_min", min(xs) if xs else None), ("cross_median", statistics.median(xs) if xs else None),
                              ("cross_max", max(xs) if xs else None)):
                want = f(r[f"{k}_{stat}"])
                if want is None and val is None:
                    continue
                agg_n += 1
                if not close(val, want):
                    agg_fail.append({"resid": resid, "metric": k, "stat": stat, "recorded": want, "new": val})
    C["C3_struct2_all_12_pairings"] = {
        "n_pairings_run": len(per), "exit_codes": {k: v["exit_code"] for k, v in per.items()},
        "global_rmsd_max_abs_difference": max(abs(v["global"]["rmsd_qcp_angstrom"] - s2res["global_rmsd_A"]["cross"][k])
                                              for k, v in per.items() if "global" in v),
        "n_aggregate_values_compared": agg_n, "n_aggregate_mismatches": len(agg_fail), "mismatches": agg_fail[:10]}

    # C4: explicit pairs, missing residues kept and excluded
    rc, msg = cli(HERE / "specs/kras2_5US4A_vs_5USJA_explicit_exclude.json", build / "explicit_exclude")
    ex = rows(build / "explicit_exclude/exclusions.csv")
    newx = by_resid(rows(build / "explicit_exclude/windows.csv"))
    shared = {r: new2[r] for r in new2 if r in newx and all(not (59 <= x <= 68) for x in range(r - 1, r + 3))}
    missing_loc = {k: sorted({e["location"] for e in ex if e["structure"] == k and e["cause"] == "missing_coordinates"},
                             key=lambda s: int(s.split(":")[1])) for k in "AB"}
    C["C4_explicit_pairs_on_missing_exclude"] = {
        "exit_code": rc, "n_exclusion_rows": len(ex), "causes": sorted({e["cause"] for e in ex}),
        "missing_locations": missing_loc,
        "expected_missing": {"A": [f"A:{r}" for r in range(60, 65)], "B": [f"A:{r}" for r in range(61, 64)]},
        "windows_away_from_gap_identical_to_C2": compare_values(newx, {r: {k: v[COL[k]] for k in METRICS}
                                                                       for r, v in shared.items()},
                                                                rec_col=lambda k: k)}

    # C5: residue-name mismatch flagged / rejected
    rc, msg = cli(HERE / "specs/demo_9O0NA_vs_5US4B_mismatch_flag.json", build / "mismatch_flag")
    res = rows(build / "mismatch_flag/residues.csv")
    mm = [(r["A_residue"], r["A_resname"], r["B_resname"]) for r in res if r["resname_mismatch"] == "True"]
    win = rows(build / "mismatch_flag/windows.csv")
    flagged = sorted(int(r["A_residue"].split(":")[1]) for r in win if r["mismatch_in_window3"] == "True")
    sp = json.loads((HERE / "specs/demo_9O0NA_vs_5US4B_mismatch_flag.json").read_text())
    sp["resname_mismatch"] = "reject"
    for k in "AB":
        sp["structures"][k]["file"] = str((HERE / "specs" / sp["structures"][k]["file"]).resolve())
    (build / "spec_mismatch_reject.json").write_text(json.dumps(sp))
    rc2, msg2 = cli(build / "spec_mismatch_reject.json", build / "mismatch_reject")
    C["C5_resname_mismatch"] = {"flag_exit_code": rc, "mismatches": mm, "window3_centres_flagged": flagged,
                                "reject_exit_code": rc2, "reject_message": msg2.splitlines()[-1]}

    # C6/C7: refusals
    rc, msg = cli(HERE / "specs/kras1_9O0N_vs_9O0O.json", build / "kras1")
    sp = json.loads((HERE / "specs/kras2_5US4A_vs_5USJA.json").read_text())
    del sp["structures"]["B"]["chain"]
    for k in "AB":
        sp["structures"][k]["file"] = str((HERE / "specs" / sp["structures"][k]["file"]).resolve())
    (build / "spec_no_chain.json").write_text(json.dumps(sp))
    rc2, msg2 = cli(build / "spec_no_chain.json", build / "no_chain")
    sp = json.loads((HERE / "specs/kras2_5US4A_vs_5USJA.json").read_text())
    sp["structures"]["A"]["chain"] = "Z"
    for k in "AB":
        sp["structures"][k]["file"] = str((HERE / "specs" / sp["structures"][k]["file"]).resolve())
    (build / "spec_absent_chain.json").write_text(json.dumps(sp))
    rc3, msg3 = cli(build / "spec_absent_chain.json", build / "absent_chain")
    C["C6_nonempty_output_refused"] = {"exit_code": rc, "message": msg.splitlines()[-1]}
    C["C7_missing_chain_refused"] = {"exit_code": rc2, "message": msg2.splitlines()[-1]}
    C["C8_absent_chain_refused"] = {"exit_code": rc3, "message": msg3.splitlines()[-1]}

    ru, rch = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    out["resources"] = {"cpu_s": ru.ru_utime + ru.ru_stime + rch.ru_utime + rch.ru_stime}
    out["build_dir"] = str(build.relative_to(ROOT))
    (HERE / "results").mkdir(exist_ok=True)
    (HERE / "results/verification.json").write_text(json.dumps(out, indent=1) + "\n")
    for name in ("kras1", "kras2", "explicit_exclude", "mismatch_flag"):
        dst = HERE / "results" / name
        if dst.exists():
            shutil.rmtree(dst)
        shutil.copytree(build / name, dst, ignore=shutil.ignore_patterns("work"))
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ("mismatches",)} for k, v in C.items()}, indent=1,
                     default=str)[:6000])


if __name__ == "__main__":
    main()
