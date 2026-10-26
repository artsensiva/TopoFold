"""Engineering regression replay of the two affected KIN-PAIR-1 pairs (PIM1-1, PIM1-2) with phospho support (structpair 0.5.0), plus the four
unaffected pairs payload against payload. Not a new evaluation; KIN-PAIR-1 and its frozen outcome are unchanged; no tolerance is changed.

    .venv/bin/python benchmarks/audit/structpair_phospho_support_2026-10-26/replay_pim1.py

Order of work (the order matters and is recorded):
1. derived specifications are written (originals untouched) and the EXPECTED window membership and values are computed independently
   from the revised correspondence and continuity rules, then written to results/expected/ with a hash, BEFORE structpair runs on them;
2. structpair is run on the derived specs (policy phospho_pdb) and on the original specs (policy off, numeric record on, as baseline);
3. membership, values (against the independent references), unchanged windows (by complete residue identifiers), global fit and aggregates are compared.
Outputs: results/phospho_replay.json, results/derived_specs/, results/runs/ (compact).
"""

import csv
import datetime as dt
import hashlib
import json
import os
import resource
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
KP = ROOT / "benchmarks/audit/kin_pair_1"
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(KP))
sys.path.insert(0, str(ROOT / "benchmarks/audit/modified_residue_inventory_2026-10-26"))
import kinpair as kp  # noqa: E402
import inventory as inv  # noqa: E402  (parse_header only)
from topofold_inputs import numeric_record as nr  # noqa: E402

ADAPTER = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
WORK = HERE / "work"
RES = HERE / "results"
PHOSPHO = {"SEP": "SER", "TPO": "THR", "PTR": "TYR"}
TOL = 1e-4
METRICS = ("d_kappa", "d_i_i2", "rmsd3", "d_tau", "d_i_i3", "rmsd4")


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def cpu():
    a, b = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


def run(spec_path, out):
    shutil.rmtree(out, ignore_errors=True)
    r = subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(spec_path), "--out", str(out)], capture_output=True, text=True,
                       env=dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=ADAPTER))
    if r.returncode != 0:
        raise RuntimeError(f"{spec_path.name}: {r.stderr[-400:]}")
    return r.stdout.strip()


def key_of(a):
    s = str(a)
    return (int(s.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")), s[len(s.rstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")):])


def independent_admission(P, modres, chain):
    """The policy's rule, applied independently to the coordinate records: name, MODRES, N/CA/C, peptide bond to a neighbour."""
    res = P.residues[chain]
    admitted = []
    for i, r in enumerate(res):
        if r.resname not in PHOSPHO or r.ca is None or r.n is None or r.c is None:
            continue
        d = []
        if i and res[i - 1].c is not None:
            d.append(float(np.linalg.norm(res[i - 1].c - r.n)))
        if i + 1 < len(res) and res[i + 1].n is not None:
            d.append(float(np.linalg.norm(r.c - res[i + 1].n)))
        ann = modres.get((chain, r.resseq, r.icode))
        if d and min(d) <= 2.0 and ann and ann["resname"] == r.resname and ann["parent"] == PHOSPHO[r.resname]:
            r.read = True
            admitted.append(r)
    return admitted


def expected_for(spec, chain_a, chain_b, files):
    """Expected window membership and values from the (derived) correspondence and the continuity rules, independent of structpair."""
    P = {k: kp.parse_pdb(Path(files[k]).read_text(encoding="latin-1"), "blank_then_A") for k in "AB"}
    mod = {k: inv.parse_header(Path(files[k]).read_text(encoding="latin-1"))[0] for k in "AB"}
    adm = {k: independent_admission(P[k], mod[k], chain_a if k == "A" else chain_b) for k in "AB"}
    ra = {r.key: r for r in P["A"].residues[chain_a]}
    rb = {r.key: r for r in P["B"].residues[chain_b]}
    pairs = [(key_of(a), key_of(b)) for a, b in spec["correspondence"]["pairs"]]
    la, lb = [ra[a] for a, _ in pairs], [rb[b] for _, b in pairs]
    assert all(r.read for r in la + lb), "a position of the derived correspondence is not readable under the independent rule"
    w = kp.independent_windows(la, lb)
    n = len(pairs)
    ident = [f"{chain_a}:{a[0]}{a[1]}|{chain_b}:{b[0]}{b[1]}" for a, b in pairs]
    win3 = {tuple(ident[k] for k in range(j - 1, j + 2)): j for j in range(1, n - 1) if w["valid3"][j]}
    win4 = {tuple(ident[k] for k in range(j - 1, j + 3)): j for j in range(1, n - 2) if w["valid4"][j]}
    return {"n_positions": n, "positions": ident, "admitted_modified_residues": {k: [f"{r.chain}:{r.resseq}{r.icode} {r.resname}" for r in adm[k]] for k in "AB"},
            "valid3_centres": [j for j in range(n) if w["valid3"][j]], "valid4_indices": [j for j in range(n) if w["valid4"][j]],
            "excluded3": [j for j in range(1, n - 1) if not w["valid3"][j]], "excluded4": [j for j in range(1, n - 2) if not w["valid4"][j]],
            "windows3": [list(k) for k in win3], "windows4": [list(k) for k in win4], "metrics": {m: w[m] for m in METRICS}, "global_rmsd": w["global_rmsd"]}


def windows_by_identity(rec):
    """{('3'|'4', window residue-identity tuple): {metric: value}} from a numeric record (complete identifiers, not row positions)."""
    pos = [f"{p['A_residue']}|{p['B_residue']}" for p in rec["positions"]]
    out = {}
    n = len(pos)
    for j in range(n):
        if rec["metrics"]["d_kappa"]["values"][j] is not None and 1 <= j <= n - 2:
            out[("3", tuple(pos[k] for k in range(j - 1, j + 2)))] = {m: rec["metrics"][m]["values"][j] for m in ("d_kappa", "d_i_i2", "rmsd3", "kappa_A", "kappa_B", "d2_A", "d2_B")}
        if rec["metrics"]["d_tau"]["values"][j] is not None and 1 <= j <= n - 3:
            out[("4", tuple(pos[k] for k in range(j - 1, j + 3)))] = {m: rec["metrics"][m]["values"][j] for m in ("d_tau", "d_i_i3", "rmsd4", "tau_A", "tau_B", "d3_A", "d3_B")}
    return out


def main():
    c0, t0 = cpu(), dt.datetime.now(dt.timezone.utc)
    shutil.rmtree(WORK, ignore_errors=True)
    for d in ("derived_specs", "expected", "runs"):
        shutil.rmtree(RES / d, ignore_errors=True)
        (RES / d).mkdir(parents=True)
    selection = json.loads((KP / "results/KIN-PAIR-1_selection.json").read_text())
    panel = {p["pair_id"]: p for p in selection["panel"]}
    out = {"id": "PHOSPHO-SUPPORT-REPLAY-PIM1-2026-10-26", "note": "engineering regression; KIN-PAIR-1 is not reclassified; no tolerance changed; not prospective validation",
           "affected": {}, "unaffected": {}}
    affected = {"PIM1-1": 261, "PIM1-2": 261}
    # 1. derived specs and expected membership, before any structpair run on them
    expected_hashes = {}
    derived = {}
    for pid, resnum in affected.items():
        orig = json.loads((KP / "work" / pid / "spec.json").read_text())
        pairs = orig["correspondence"]["pairs"]
        flat = [(key_of(a), key_of(b)) for a, b in pairs]
        i = flat.index(((resnum - 1, ""), (resnum - 1, "")))
        assert flat[i + 1] == ((resnum + 1, ""), (resnum + 1, "")), "the omitted position is not the single gap between its neighbours"
        spec = json.loads(json.dumps(orig))
        spec["correspondence"]["pairs"] = pairs[: i + 1] + [[resnum, resnum]] + pairs[i + 1:]
        spec["modified_residues"] = "phospho_pdb"
        spec["numeric_record"] = True
        spec["name"] = orig["name"] + f" (derived regression spec: position {resnum} added, modified_residues phospho_pdb)"
        sp_path = RES / "derived_specs" / f"{pid}_derived.json"
        sp_path.write_text(json.dumps(spec, indent=1))
        files = {k: spec["structures"][k]["file"] for k in "AB"}
        exp = expected_for(spec, spec["structures"]["A"]["chain"], spec["structures"]["B"]["chain"], files)
        (RES / "expected" / f"{pid}_expected.json").write_text(json.dumps(exp) + "\n")
        expected_hashes[pid] = {"sha256": sha(RES / "expected" / f"{pid}_expected.json"), "written_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
        derived[pid] = (sp_path, exp)
    out["expected_written_before_structpair_on_derived_specs"] = expected_hashes
    # 2. runs: derived (policy on), originals with numeric record (baseline, policy off)
    for pid, (sp_path, exp) in derived.items():
        p = panel[pid]
        new_out = WORK / f"{pid}_new"
        run(sp_path, new_out)
        o = json.loads((KP / "work" / pid / "spec.json").read_text())
        o["numeric_record"] = True
        (WORK / f"{pid}_orig_numeric.json").write_text(json.dumps(o))
        old_out = WORK / f"{pid}_old"
        run(WORK / f"{pid}_orig_numeric.json", old_out)
        rn, ro = nr.read(new_out), nr.read(old_out)
        # 3a. membership against the expectation fixed beforehand
        wn = windows_by_identity(rn)
        got3 = {k[1] for k in wn if k[0] == "3"}
        got4 = {k[1] for k in wn if k[0] == "4"}
        exp3, exp4 = {tuple(x) for x in exp["windows3"]}, {tuple(x) for x in exp["windows4"]}
        ex_rows = list(csv.DictReader(open(new_out / "exclusions.csv")))
        got_ex3 = sorted({int(r["position"]) for r in ex_rows if r["window_size"] == "3"})
        got_ex4 = sorted({int(r["position"]) for r in ex_rows if r["window_size"] == "4"})
        membership = {"windows3_equal_expected": got3 == exp3, "windows4_equal_expected": got4 == exp4, "n_windows3": len(got3), "n_windows4": len(got4),
                      "excluded3_equal_expected": got_ex3 == exp["excluded3"], "excluded4_equal_expected": got_ex4 == exp["excluded4"],
                      "n_excluded_rows_by_cause": {}}
        for r in ex_rows:
            membership["n_excluded_rows_by_cause"][r["cause"]] = membership["n_excluded_rows_by_cause"].get(r["cause"], 0) + 1
        # 3b. values against the independent references
        worst = {m: max(abs(a - b) for a, b in zip(rn["metrics"][m]["values"], exp["metrics"][m]) if a is not None and b is not None) for m in METRICS}
        none_equal = all((a is None) == (b is None) for m in METRICS for a, b in zip(rn["metrics"][m]["values"], exp["metrics"][m]))
        values = {"max_abs_deviation_vs_independent": worst, "missing_pattern_equal": none_equal, "global_rmsd_workflow": rn["global"]["rmsd_qcp_angstrom"], "global_rmsd_independent": exp["global_rmsd"],
                  "global_abs_deviation": abs(rn["global"]["rmsd_qcp_angstrom"] - exp["global_rmsd"]), "tolerance": TOL,
                  "pass": all(v <= TOL for v in worst.values()) and none_equal and abs(rn["global"]["rmsd_qcp_angstrom"] - exp["global_rmsd"]) <= TOL}
        # 3c. unchanged windows by complete identifiers
        wo = windows_by_identity(ro)
        common = set(wo) & set(wn)
        unchanged = {"old_windows_with_values": len(wo), "found_in_new_with_the_same_identifiers": len(common),
                     "values_exactly_equal": sum(all(wo[k][m] == wn[k][m] for m in wo[k]) for k in common), "old_windows_missing_in_new": len(set(wo) - set(wn)),
                     "new_windows_not_in_old": len(set(wn) - set(wo)),
                     "new_windows_not_in_old_contain_modified_position": all(any("261" in x.split("|")[0] for x in k[1]) for k in set(wn) - set(wo))}
        # 3d. global fit and aggregates
        so = json.loads((old_out / "summary.json").read_text())
        sn = json.loads((new_out / "summary.json").read_text())
        pos_new = [f"{q['A_residue']}|{q['B_residue']}" for q in rn["positions"]]
        sep_idx = next(i for i, q in enumerate(rn["positions"]) if q["A_modification"] or q["B_modification"])
        tops = {}
        for m in ("d_kappa", "rmsd3", "d_tau", "rmsd4"):
            def top(rec):
                rows = [(v, i) for i, v in enumerate(rec["metrics"][m]["values"]) if v is not None]
                rows.sort(key=lambda t: (-t[0], t[1]))
                size = 3 if m in ("d_kappa", "rmsd3") else 4
                pos = [f"{q['A_residue']}|{q['B_residue']}" for q in rec["positions"]]
                return [tuple(pos[k] for k in range(i - 1, i - 1 + size)) for _, i in rows[:10]]
            to, tn = top(ro), top(rn)
            tops[m] = {"top10_old_in_new_top10": len(set(to) & set(tn)), "top10_new_not_in_old_top10": [list(t) for t in tn if t not in to]}
        agg = {"n_positions_old": ro["n_positions"], "n_positions_new": rn["n_positions"], "n_ca_old": ro["global"]["n_ca"], "n_ca_new": rn["global"]["n_ca"],
               "global_rmsd_old": ro["global"]["rmsd_qcp_angstrom"], "global_rmsd_new": rn["global"]["rmsd_qcp_angstrom"],
               "global_rmsd_change": rn["global"]["rmsd_qcp_angstrom"] - ro["global"]["rmsd_qcp_angstrom"],
               "explanation": "one more Cα (the modified residue) enters the fit, so the global RMSD is refitted over n+1 positions; local metrics of unchanged windows do not depend on the global fit",
               "n_windows_with_values_old": so["n_windows"], "n_windows_with_values_new": sn["n_windows"], "top10_changes": tops,
               "altloc_sensitivity_old": so["altloc_sensitivity"], "altloc_sensitivity_new": sn["altloc_sensitivity"],
               "modified_position_new": pos_new[sep_idx]}
        mm = sn["modified_residues"]
        out["affected"][pid] = {"derived_spec": str(sp_path.relative_to(ROOT)), "original_spec": f"benchmarks/audit/kin_pair_1/work/{pid}/spec.json (preserved)",
                                "entries": [p["A"], p["B"]], "expected": {k: exp[k] for k in ("n_positions", "admitted_modified_residues", "excluded3", "excluded4")},
                                "membership": membership, "values": values, "unchanged_windows": unchanged, "aggregates": agg,
                                "modified_residues_block": {"retained": mm["retained"], "excluded": mm["excluded"]},
                                "coverage": {"positions": [ro["n_positions"], rn["n_positions"]],
                                             "windows3_with_values": [len({k for k in wo if k[0] == "3"}), len(got3)], "windows4_with_values": [len({k for k in wo if k[0] == "4"}), len(got4)],
                                             "excluded_windows3": [len({int(r["position"]) for r in csv.DictReader(open(old_out / "exclusions.csv")) if r["window_size"] == "3"}), len(got_ex3)],
                                             "excluded_windows4": [len({int(r["position"]) for r in csv.DictReader(open(old_out / "exclusions.csv")) if r["window_size"] == "4"}), len(got_ex4)]},
                                "structpair_version": json.loads((new_out / "manifest.json").read_text())["software"]["structpair_version"],
                                "pass": all(membership[k] for k in ("windows3_equal_expected", "windows4_equal_expected", "excluded3_equal_expected", "excluded4_equal_expected")) and values["pass"]
                                and unchanged["found_in_new_with_the_same_identifiers"] == unchanged["old_windows_with_values"] == unchanged["values_exactly_equal"]}
        dst = RES / "runs" / pid
        dst.mkdir(parents=True)
        for tag, d_ in (("new", new_out), ("old_with_numeric_record", old_out)):
            (dst / tag).mkdir()
            for f in ("windows.csv", "residues.csv", "exclusions.csv", "summary.json", "manifest.json", "numeric_record.json", "report.md"):
                if (d_ / f).is_file():
                    shutil.copy2(d_ / f, dst / tag / f)
    # 4. unaffected pairs: payload against payload, policy off vs on
    for pid in ("EGFR-1", "EGFR-2", "BTK-1", "BTK-2"):
        o = json.loads((KP / "work" / pid / "spec.json").read_text())
        o["numeric_record"] = True
        (WORK / f"{pid}_off.json").write_text(json.dumps(o))
        o2 = dict(o, modified_residues="phospho_pdb")
        (WORK / f"{pid}_on.json").write_text(json.dumps(o2))
        run(WORK / f"{pid}_off.json", WORK / f"{pid}_off")
        run(WORK / f"{pid}_on.json", WORK / f"{pid}_on")
        a, b = nr.read(WORK / f"{pid}_off"), nr.read(WORK / f"{pid}_on")
        strip = lambda r: {"metrics": r["metrics"], "global": r["global"], "altloc_sensitivity": r["altloc_sensitivity"],   # noqa: E731
                           "positions": [{k: v for k, v in p.items() if not k.endswith("_modification")} for p in r["positions"]]}
        wa, wb = list(csv.DictReader(open(WORK / f"{pid}_off/windows.csv"))), list(csv.DictReader(open(WORK / f"{pid}_on/windows.csv")))
        ma, mb = (json.loads((WORK / f"{pid}_{t}/manifest.json").read_text()) for t in ("off", "on"))
        out["unaffected"][pid] = {"numeric_payload_identical": strip(a) == strip(b),
                                  "windows_csv_old_columns_identical": all({k: v for k, v in y.items() if k in x} == x for x, y in zip(wa, wb)),
                                  "residues_exclusions_identical": all((WORK / f"{pid}_off" / f).read_bytes() == (WORK / f"{pid}_on" / f).read_bytes() for f in ("exclusions.csv",)),
                                  "added_columns_only_when_policy_on": sorted(set(wb[0]) - set(wa[0])),
                                  "modified_residues_retained": len(json.loads((WORK / f"{pid}_on/summary.json").read_text())["modified_residues"]["retained"]),
                                  "expected_provenance_differences": sorted(k for k in ("policies",) if ma[k] != mb[k]) + ["spec sha256", "output hashes"]}
    out["resources"] = {"cpu_seconds": cpu() - c0, "wall_seconds": (dt.datetime.now(dt.timezone.utc) - t0).total_seconds()}
    out["all_checks_pass"] = all(v["pass"] for v in out["affected"].values()) and all(v["numeric_payload_identical"] and v["windows_csv_old_columns_identical"] and
                                                                                    v["residues_exclusions_identical"] and v["modified_residues_retained"] == 0 for v in out["unaffected"].values())
    (RES / "phospho_replay.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
    print(json.dumps({"affected": {k: {"pass": v["pass"], "coverage": v["coverage"]} for k, v in out["affected"].items()},
                      "unaffected": {k: v["numeric_payload_identical"] for k, v in out["unaffected"].items()}, "all": out["all_checks_pass"]}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
