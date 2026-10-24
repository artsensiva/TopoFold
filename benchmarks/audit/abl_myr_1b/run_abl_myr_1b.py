"""ABL-MYR-1b wrapper: the post-hoc trajectory comparison of protocol_1b.json with the dated amendment protocol_amendment_a1.json.

    PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1b/run_abl_myr_1b.py freeze-code
    PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1b/run_abl_myr_1b.py run --adapter-python <python with MDAnalysis>
    PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1b/run_abl_myr_1b.py report

The frozen ABL-MYR-1 runner (run_abl_myr.py) is imported unchanged and does all the work. This module replaces only, while
the runner runs:
  * G4 (system.gro / mdp / topology identity) by G4b: the ABL-MYR-1a and ABL-MYR-1b a2 evidence is verified by hash and by
    content; it does not claim identical production starts;
  * the required Set2 members: Set2/system.gro is not published, so system_em.gro is listed in its place only so that the
    unchanged pipeline can run; it is not used as a start and G4b does not read it;
  * G5 by the bilateral rule of the amendment: if any of the four groups keeps fewer than min_runs usable runs, the run stops;
  * the outcome wording (restricted labels, no 'reproduced') and the report.
Nothing else of the frozen runner, its thresholds or its method is changed.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import datetime as dt
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ORIG = ROOT / "benchmarks/audit/abl_myr_1"
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ORIG))
sys.path.insert(0, str(HERE))
import run_abl_myr as R0  # noqa: E402  (frozen runner, unchanged)
import write_protocol_1b as W  # noqa: E402  (protocol lineage and labels)

AMENDMENT = HERE / "protocol_amendment_a1.json"
AMEND_FREEZE = HERE / "results/ABL-MYR-1b_protocol_amendment_a1_freeze.json"
FREEZE = HERE / "results/ABL-MYR-1b_wrapper_freeze.json"
RESULTS, BLOCKED, REPORT = "ABL-MYR-1b_results.json", "ABL-MYR-1b_blocked.json", "ABL-MYR-1b_report.md"
SELF = Path(__file__).resolve()
TESTS = ["tests/test_abl_myr_1b_wrapper.py"]
AUDIT_RESULTS = ROOT / "benchmarks/audit/abl_myr_1a/results/ABL-MYR-1a_audit_results.json"
A2_RESULTS = W.A2_RESULTS
NO_SYSTEM_GRO = ("Set2",)
ALPHA_I = (518, 531)   # alphaI residues named in the source study (eLife 85216), isoform 1b numbering
STEPS = {"inconclusive": 2, "TopoFold adds reproduced residues": 3, "baselines cover TopoFold": 4, "partial": 5}
CHANNEL_ROLE = {"kappa": "TopoFold", "tau": "TopoFold", "d_i_i2": "Calpha baseline", "d_i_i3": "Calpha baseline",
                "rmsf": "Calpha baseline", "phi": "stronger-information comparator", "psi": "stronger-information comparator"}


# --- G0 lineage ------------------------------------------------------------------------------------------------------

def hash_mismatches(recorded: dict, root: Path = ROOT) -> list:
    return [p for p, h in recorded.items() if not (Path(root) / p).is_file() or R0.sha256_file(Path(root) / p) != h]


def lineage_problems() -> list:
    """Every lineage hash of the protocol freeze, the protocol itself, the amendment freeze and the original code freeze."""
    pf = json.loads(W.FREEZE.read_text())
    bad = []
    for name, ok in (("protocol_1b.json", pf["protocol_sha256"] == R0.sha256_file(W.PROTOCOL)),
                     ("PROTOCOL_1B.md", pf["protocol_md_sha256"] == R0.sha256_file(W.PROTOCOL_MD)),
                     ("write_protocol_1b.py", pf["writer_sha256"] == R0.sha256_file(W.SELF)),
                     ("a2 results", pf["a2_results_sha256"] == R0.sha256_file(W.A2_RESULTS))):
        if not ok:
            bad.append(name)
    bad += hash_mismatches({v: pf["lineage_sha256"][k] for k, v in W.LINEAGE.items()})
    bad += hash_mismatches(json.loads(AMEND_FREEZE.read_text())["sha256"])
    return bad


# --- G4b -------------------------------------------------------------------------------------------------------------

def g4b_check(audit: dict, a2: dict) -> dict:
    """Evidence replacing G4 (see protocol_1b.json). Raises GateError('G4b') unless every item holds."""
    q, d = audit.get("questions", {}), audit.get("derived", {})
    problems = []
    if q.get("Q1", {}).get("category") != "identical":
        problems.append("ABL-MYR-1a Q1 (Set1 system.gro, Myr vs noMyr) is not identical")
    for c in ("Myr", "noMyr"):
        if q.get("Q4", {}).get(c, {}).get("category") != "identical":
            problems.append(f"ABL-MYR-1a Q4 (Set2 system_em.gro is a copy of Set1's) is not identical for {c}")
    if d.get("set2_start_direct_evidence") is not True:
        problems.append("ABL-MYR-1a set2_start_direct_evidence is not true")
    crit = a2.get("decision", {}).get("criteria", {})
    if a2.get("outcome") != "comparable_for_1b" or not crit or not all(crit.values()):
        problems.append("ABL-MYR-1b a2 outcome is not comparable_for_1b with K1-K5 all true")
    if problems:
        raise R0.GateError("G4b", "comparability evidence does not hold: " + "; ".join(problems), problems)
    return {"audit_verdict_recorded": audit.get("verdict"), "audit_Q1": q["Q1"]["category"],
            "audit_Q4": {c: q["Q4"][c]["category"] for c in ("Myr", "noMyr")}, "set2_start_direct_evidence": True,
            "a2_outcome": a2["outcome"], "a2_criteria": crit,
            "note": "G4b verifies recorded evidence. It does not establish identical production starts and does not read any Set2 coordinate file"}


# --- G1 members, G5 bilateral rule ---------------------------------------------------------------------------------------

def make_wanted_members(base):
    def wanted(cfg):
        out = {a: dict(m) for a, m in base(cfg).items()}
        for arch, m in out.items():
            for s in NO_SYSTEM_GRO:
                if f"gro:{s}" in m:
                    m[f"gro:{s}"] = f"{s}/system_em.gro"       # listed so that the unchanged pipeline runs; never used as a start
        return out
    return wanted


def make_bilateral_gate(cfg: dict, base):
    """frame_gate replacement: same result per run, plus the amendment's C1: a group below min_runs stops the run."""
    order = [(s, c, r) for s, reps in cfg["sets"].items() for c in cfg["conditions"] for r in reps]   # the runner's G5 loop order
    state = {"i": 0, "usable": {}, "excluded": {}}

    def gate(times_ps, cfg_arg):
        if state["i"] >= len(order):
            raise R0.GateError("G5", "internal: more frame checks than runs")
        s, c, r = order[state["i"]]
        state["i"] += 1
        fg = base(times_ps, cfg_arg)
        key = (s, c)
        state["usable"].setdefault(key, 0)
        if fg["excluded"]:
            state["excluded"].setdefault(key, []).append({"run": f"{s}:{c}:{r}", "cause": fg["cause"]})
        else:
            state["usable"][key] += 1
        if r == cfg["sets"][s][-1] and state["usable"][key] < cfg_arg["min_runs"]:
            raise R0.GateError("G5", f"bilateral rule (amendment C1): {s} {c} has {state['usable'][key]} usable runs "
                                     f"(minimum {cfg_arg['min_runs']}); the bilateral comparison is stopped",
                               {"group": f"{s}:{c}", "usable": state["usable"][key], "min_runs": cfg_arg["min_runs"],
                                "excluded": state["excluded"].get(key, [])})
        return fg
    gate.state = state
    return gate


@contextlib.contextmanager
def patched(cfg: dict, evidence: dict):
    saved = (R0.check_g4_set, R0.wanted_members, R0.frame_gate)
    cache = {}

    def g4b(files, cfg_arg):
        if "r" not in cache:
            cache["r"] = g4b_check(evidence["audit"], evidence["a2"])
        return dict(cache["r"])
    R0.check_g4_set, R0.wanted_members = g4b, make_wanted_members(saved[1])
    R0.frame_gate = make_bilateral_gate(cfg, saved[2])
    try:
        yield
    finally:
        R0.check_g4_set, R0.wanted_members, R0.frame_gate = saved


# --- post-processing ---------------------------------------------------------------------------------------------------

def alpha_i_statement(g3: dict) -> dict:
    lo, hi = ALPHA_I
    a = g3["alignment"]
    span = a["reference_span"]
    gaps = [(x, y) for x, y in a["unobserved_reference_ranges"]]
    missing = [i for i in range(lo, hi + 1) if not (span and span[0] <= i <= span[1]) or any(x <= i <= y for x, y in gaps)]
    return {"residues_isoform_1b": f"{lo}-{hi}", "reference_span": span, "not_represented": missing,
            "represented": not missing,
            "consequence": ("all alphaI residues named in the source study are represented" if not missing else
                            "alphaI is not (fully) represented: the analysis is limited to the represented kinase-domain residues "
                            "and no claim is made about the alphaI helix or the Myr-binding region")}


def production_starts(a2: dict, retained: dict) -> dict:
    per = a2["R1_usable_coordinates"]["per_run"]
    out = {}
    for s in R0.CFG["sets"]:
        keys = [k for c in retained.get(s, {}) for k in retained[s][c]]
        out[s] = {"per_run": {k: {"ca_rmsd_to_system_em_A": per[k]["ca_rmsd_to_em_A"], "box_volume_A3": per[k].get("box_volume_A3"),
                                  "velocities_present": per[k].get("velocities_present")} for k in keys if k in per},
                  "geometry_from_a2": a2["R2_production_start_geometry"][s],
                  "note": "descriptive start-geometry diagnostic (Calpha RMSD); not a calibrated control for the window statistic E (amendment C3)"}
    return out


def _rename(x):
    if isinstance(x, dict):
        return {(k.replace("reproduced", "robust_in_both_sets") if isinstance(k, str) else k): _rename(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_rename(v) for v in x]
    return x.replace("reproduced", "robust in both sets") if isinstance(x, str) else x


def relabel(res: dict, a2: dict) -> dict:
    """Restricted outcome wording, per-set leave-one-run-out counts, alphaI statement and production starts."""
    res = dict(res)
    res["id"] = "ABL-MYR-1b"
    res["post_hoc_protocol"] = True
    o = dict(res["outcome"])
    frozen = o["label"]
    o["label"] = W.LABELS[frozen]
    o["decision_rule_step"] = STEPS[frozen]
    res["outcome"] = _rename(o)
    retained = {s: d["runs"] for s, d in res.get("sets", {}).items()}
    for s, d in res.get("sets", {}).items():
        d["leave_one_run_out_recomputations"] = sum(len(v) for v in d["runs"].values())
    res["alpha_I"] = alpha_i_statement(res["gates"]["G3"])
    res["production_starts"] = production_starts(a2, retained)
    return res


def run_wrapped(spec, cfg, data_dir, out_dir, adapter_python, evidence, work=None, log=print) -> dict:
    """The frozen pipeline under the replacements above, then the restricted wording. Raises GateError on a failed gate."""
    with patched(cfg, evidence):
        res = R0.run_pipeline(spec, cfg, data_dir, out_dir, adapter_python, work=work, log=log)
    return relabel(res, evidence["a2"])


# --- report ------------------------------------------------------------------------------------------------------------

def render_report(res: dict) -> str:
    o = res["outcome"]
    L = ["# ABL-MYR-1b results (generated from ABL-MYR-1b_results.json)", "",
         "**Post-hoc, descriptive evidence from simulated Abl kinase-domain ensembles (de Buhr and Graeter, Dryad).** The study is "
         "separate from Ayaz WT/M472I, HDX, drug response and biological function. The two thermostat sets are a consistency "
         "check, not independent replication. No causal, clinical or novelty claim is made.", "",
         f"**Outcome: {o['label']}** (decision-rule step {o['decision_rule_step']})" + (f". {o['note']}" if o.get("note") else ""), ""]
    g = res["gates"]
    archives = ", ".join("{} {} B".format(k, v["bytes"]) for k, v in g["G0"].items())
    L += ["## Gates", "", f"- G0 archives: {archives}; lineage hashes verified: "
          f"{res.get('lineage_verified')}", f"- G1: {g['G1']['members']} members extracted ({g['G1']['extracted_bytes']} B)",
          f"- G2: topology {g['G2']['topology']}, {g['G2']['n_atoms']} atoms, {g['G2']['n_residues']} residues, "
          f"φ/ψ comparator available {g['G2']['phi_psi_comparator_available']}",
          f"- G3: identity {g['G3']['alignment']['identity']:.4f}, reference span {g['G3']['alignment']['reference_span']}, "
          f"unobserved ranges {g['G3']['alignment']['unobserved_reference_ranges']}",
          f"- G4b: {next(iter(g['G4'].values()))['note']}",
          f"- G5: {len(g['G5']['excluded_runs'])} runs excluded; sets analysed {g['G5']['sets_analysed']}",
          f"- G6: pbc used {g['G6']['pbc_used']}; worst excluded-window fraction "
          f"{max(p['worst_fraction'] for p in g['G6']['passes'])}", "",
          "## αI helix", "", f"Residues {res['alpha_I']['residues_isoform_1b']} represented: {res['alpha_I']['represented']}. "
          f"{res['alpha_I']['consequence']}.", ""]
    for s, d in res.get("sets", {}).items():
        L += [f"## {s}", "", f"Runs retained: " + ", ".join(f"{c} {len(v)}" for c, v in d["runs"].items()) +
              f"; leave-one-run-out recomputations: {d['leave_one_run_out_recomputations']}.", "",
              "| Channel | Role | Windows | E > 0 | Robust (full and every leave-one-run-out) | Median E |", "|---|---|---|---|---|---|"]
        for c, v in d["channels"].items():
            m = "" if v["E_median"] is None else f"{v['E_median']:.4f}"
            L.append(f"| {c} | {CHANNEL_ROLE[c]} | {v['windows_evaluated']} | {v['windows_E_positive']} | {v['windows_exceeding']} | {m} |")
        L += ["", "Run-to-run variability (median over windows of Jensen–Shannon distances) and temporal diagnostics:", ""]
        for c, v in d["run_variability"].items():
            t = d["temporal"][c]
            L.append(f"- {c}: within Myr {v['within_Myr_js_median_of_windows']}, within noMyr {v['within_noMyr_js_median_of_windows']}, "
                     f"between {v['between_js_median_of_windows']}; lag-1 median {t['lag1_median']}, max {t['lag1_max']}; "
                     f"decorrelated within 50 frames in {t['fraction_decorrelated_within_max_lag']}")
        ps = res["production_starts"][s]
        rm = [v["ca_rmsd_to_system_em_A"] for v in ps["per_run"].values()]
        span = f"{min(rm):.3f}–{max(rm):.3f} Å" if rm else "not available for the retained runs"
        L += ["", f"Production starts (a2 results; descriptive, Calpha RMSD to system_em.gro of {len(rm)} retained runs): {span}.", ""]
    L += ["## Outcome detail", "", f"- R_T (TopoFold residues): {o['R_T']}", f"- R_B (baseline residues): {o['R_B']}",
          f"- TopoFold residues outside N2(R_B): {o['R_T_outside_N2_R_B']}",
          f"- baseline residues outside N2(R_T): {o['R_B_outside_N2_R_T']}",
          f"- φ/ψ comparator residues (reported separately): {o['R_phi_psi']}", "",
          "These are descriptive filters (E > 0 in the full set and in every leave-one-run-out recomputation; the three-residue "
          "rule). They give no error-rate control or significance, and no outcome is evidence of the superiority of a method.", ""]
    return "\n".join(L) + "\n"


# --- freeze / run / report ----------------------------------------------------------------------------------------------

def dependencies() -> list:
    return ["benchmarks/audit/abl_myr_1/run_abl_myr.py", "benchmarks/audit/abl_myr_1b/write_protocol_1b.py",
            "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_protocol_amendment_a1_freeze.json",
            "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_protocol_freeze.json",
            "benchmarks/audit/abl_myr_1b/protocol_amendment_a1_ERRATUM.md",
            "benchmarks/audit/abl_myr_1b/REVIEW_2026-10-25_independent_AI.md"] + R0.LIBRARY


def trajectories_present() -> bool:
    d = ORIG / "data_dryad"
    return (d / "extracted").exists() or any(d.rglob("*.xtc"))


def cmd_freeze_code(args) -> int:
    if FREEZE.exists():
        print("code freeze exists")
        return 1
    if trajectories_present():
        raise SystemExit("STOP: a trajectory has been extracted; the wrapper freeze must precede it")
    bad = lineage_problems()
    if bad:
        raise SystemExit(f"STOP: lineage files differ from their freeze records: {bad}")
    rec = {"id": "ABL-MYR-1b-WRAPPER", "code_frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "wrapper": {"path": str(SELF.relative_to(ROOT)), "sha256": R0.sha256_file(SELF)},
           "tests_sha256": {t: R0.sha256_file(ROOT / t) for t in TESTS},
           "dependencies_sha256": {p: R0.sha256_file(ROOT / p) for p in dependencies()},
           "amendment_sha256": R0.sha256_file(AMENDMENT), "frozen_cfg": R0.CFG, "lineage_verified": True,
           "trajectories_extracted": False, "state": "wrapper frozen before PCA.zip was opened"}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps({k: rec[k] for k in ("code_frozen_utc", "wrapper")}, indent=1))
    return 0


def verify_freeze() -> dict:
    if not FREEZE.exists():
        raise SystemExit("STOP: no wrapper freeze record; run freeze-code before opening PCA.zip")
    f = json.loads(FREEZE.read_text())
    bad = []
    if R0.sha256_file(SELF) != f["wrapper"]["sha256"]:
        bad.append("wrapper")
    bad += hash_mismatches(f["dependencies_sha256"])
    if R0.sha256_file(AMENDMENT) != f["amendment_sha256"]:
        bad.append("amendment")
    if f["frozen_cfg"] != json.loads(json.dumps(R0.CFG)):
        bad.append("runner CFG")
    if bad:
        raise SystemExit(f"STOP: files differ from the wrapper freeze: {bad}")
    R0.verify_code_freeze()                    # the frozen runner against the original code freeze
    return f


def cmd_run(args) -> int:
    f = verify_freeze()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if (out / RESULTS).exists() or (out / BLOCKED).exists():
        raise SystemExit("STOP: an ABL-MYR-1b attempt already exists in the output directory; it is preserved, not overwritten")
    bad = lineage_problems()
    if bad:
        raise SystemExit(f"STOP: lineage hashes differ: {bad}")
    spec = json.loads((ORIG / "run_spec.json").read_text())
    evidence = {"audit": json.loads(AUDIT_RESULTS.read_text()), "a2": json.loads(A2_RESULTS.read_text())}
    try:
        res = run_wrapped(spec, R0.CFG, Path(args.data_dir), out, args.adapter_python, evidence,
                          work=Path(args.work) if args.work else HERE / "work")
    except R0.GateError as e:
        rec = {"id": "ABL-MYR-1b", "status": "blocked", "gate": e.gate, "cause": R0.scrub(e.message), "detail": R0.scrub(e.detail),
               "utc": dt.datetime.now(dt.timezone.utc).isoformat(), "wrapper_sha256": f["wrapper"]["sha256"]}
        (out / BLOCKED).write_text(json.dumps(rec, indent=1, default=R0._json_default) + "\n")
        print(f"BLOCKED at {e.gate}: {e.message}")
        return 3
    res["lineage_verified"] = True
    res["wrapper_sha256"], res["amendment_sha256"] = f["wrapper"]["sha256"], f["amendment_sha256"]
    (out / RESULTS).write_text(json.dumps(res, indent=1, default=R0._json_default) + "\n")
    for s, d in res.get("sets", {}).items():
        with open(out / f"ABL-MYR-1b_windows_{s}.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["channel", "window", "E", "E_loo_min", "exceeding"], lineterminator="\n")
            w.writeheader()
            w.writerows(d["windows"])
    print(res["outcome"]["label"])
    return 0


def cmd_report(args) -> int:
    res = json.loads((Path(args.out) / RESULTS).read_text())
    (Path(args.out) / REPORT).write_text(render_report(res))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("freeze-code")
    r = sub.add_parser("run")
    r.add_argument("--adapter-python", required=True)
    r.add_argument("--data-dir", default=str(ORIG / "data_dryad"))
    r.add_argument("--out", default=str(HERE / "results"))
    r.add_argument("--work", default=None)
    p = sub.add_parser("report")
    p.add_argument("--out", default=str(HERE / "results"))
    args = ap.parse_args(argv)
    return {"freeze-code": cmd_freeze_code, "run": cmd_run, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
