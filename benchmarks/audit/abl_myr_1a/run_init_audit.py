"""ABL-MYR-1a initialization audit (post-hoc amendment, step 1). See audit_spec.json, which is hash-frozen before this code
is run on any study file.

    PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1a/run_init_audit.py freeze
    PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1a/run_init_audit.py run --adapter-python <python with MDAnalysis>
    PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1a/run_init_audit.py report

Reads only the starting files of the two model archives (system(_em).gro, production_md.mdp, topol.top, md.tpr). Opens no
trajectory and computes no descriptor. The parsers, the Myr-line rule and the archive helpers are those of the frozen
ABL-MYR-1 runner, imported unchanged (its hash is checked against the original code freeze).
"""

from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys
import time
import zipfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ORIG = ROOT / "benchmarks/audit/abl_myr_1"
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ORIG))
import run_abl_myr as R0  # noqa: E402  (frozen runner, unchanged)

SPEC = HERE / "audit_spec.json"
FREEZE = HERE / "results/ABL-MYR-1a_audit_freeze.json"
RESULTS = HERE / "results/ABL-MYR-1a_audit_results.json"
SELF = Path(__file__).resolve()
TEST = "tests/test_abl_myr_1a_audit.py"
ORIGINAL_ARTIFACTS = [f"benchmarks/audit/abl_myr_1/{p}" for p in (
    "run_spec.json", "run_abl_myr.py", "g1_member_audit.py", "results/ABL-MYR-1_freeze.json",
    "results/ABL-MYR-1_code_freeze.json", "results/ABL-MYR-1_blocked.json", "results/ABL-MYR-1_g1_member_audit.json")]
IDENTICAL_NM, SMALL_RMSD_NM = 0.001, 0.05
SETS = {"Set1": list(range(11, 21)), "Set2": list(range(21, 31))}
CONDITIONS = {"Myr": "kinase_kinked_Myr", "noMyr": "kinase_kinked_noMyr"}
INIT_KEYS = R0.INIT_KEYS


class AuditStop(Exception):
    def __init__(self, code: str, cause: str):
        self.code, self.cause = code, cause
        super().__init__(f"{code}: {cause}")


def sha256_file(p) -> str:
    return R0.sha256_file(p)


def scrub(x):
    return R0.scrub(x)


# --- comparisons (pure) ------------------------------------------------------------------------------------------------

def protein_seq_xyz(path):
    atoms = R0.protein_atoms(R0.parse_gro(path))
    return [(a[0], a[1]) for a in atoms], np.array([a[2] for a in atoms], dtype=np.float64).reshape(-1, 3)


def compare_gro_files(a_path, b_path, identical_nm=IDENTICAL_NM, small_rmsd_nm=SMALL_RMSD_NM) -> dict:
    sa, xa = protein_seq_xyz(a_path)
    sb, xb = protein_seq_xyz(b_path)
    if not sa or sa != sb:
        first = next((i for i, (p, q) in enumerate(zip(sa, sb)) if p != q), min(len(sa), len(sb)))
        return {"category": "sequence_mismatch", "n_protein_atoms": [len(sa), len(sb)], "first_difference_index": first}
    d = xa - xb
    disp = np.linalg.norm(d, axis=1)
    max_abs = round(float(np.abs(d).max()), 9)
    rmsd = round(float(np.sqrt(np.mean(disp ** 2))), 9)
    cat = "identical" if max_abs <= identical_nm else ("differs_small" if rmsd <= small_rmsd_nm else "differs_large")
    return {"category": cat, "n_protein_atoms": len(sa), "max_abs_dx_nm": max_abs, "rmsd_nm": rmsd,
            "max_atom_displacement_nm": round(float(disp.max()), 9),
            "atoms_over_0.001nm": int((disp > identical_nm).sum())}


def compare_mdp_text(a: str, b: str) -> dict:
    la, lb = R0.normalise_lines(a), R0.normalise_lines(b)
    diff = [f"Myr:{x} | noMyr:{y}" for x, y in zip(la, lb) if x != y][:10]
    if len(la) != len(lb):
        diff.append(f"line counts: {len(la)} vs {len(lb)}")
    return {"identical": la == lb, "differing_lines": diff}


def mdp_key_differences(a: str, b: str) -> dict:
    va, vb = R0.mdp_values(a), R0.mdp_values(b)
    return {k: [va.get(k), vb.get(k)] for k in sorted(set(va) | set(vb)) if va.get(k) != vb.get(k)}


def compare_top_text(myr: str, nomyr: str) -> dict:
    r = R0.compare_top(myr, nomyr)
    return {"identical_after_myr_removal": r["ok"], "differing_lines": r["differing_lines"]}


def decide(q: dict) -> dict:
    """Ordered verdict rules of audit_spec.json (stop conditions are handled before this is called)."""
    set1_pass = (q["Q1"]["category"] == "identical" and q["Q5"]["Set1"]["identical"]
                 and q["Q6"]["Set1"]["identical_after_myr_removal"])
    set2_params = q["Q5"]["Set2"]["identical"] and q["Q6"]["Set2"]["identical_after_myr_removal"]
    direct = q["Q2"]["category"] == "identical" or (
        all(q["Q4"][c]["category"] == "identical" for c in ("Myr", "noMyr")) and q["Q1"]["category"] == "identical")
    out = {"set1_original_checks_pass": bool(set1_pass), "set2_parameters_pass": bool(set2_params),
           "set2_start_direct_evidence": bool(direct)}
    if not set1_pass:
        out["verdict"] = "V4_set1_comparability_failed"
    elif direct and set2_params:
        out["verdict"] = "V1_set2_start_established_direct"
    elif q["Q2"]["category"] in ("differs_large", "sequence_mismatch") and not direct:
        out["verdict"] = "V3_set2_not_comparable"
    else:
        out["verdict"] = "V2_set2_start_not_established"
    return out


def tpr_summary(probe: dict) -> dict:
    """Q7: per set and condition from {key 'Set:cond:rep': probe record}."""
    out = {}
    for s in SETS:
        out[s] = {}
        for c in CONDITIONS:
            recs = [v for k, v in probe.items() if k.startswith(f"{s}:{c}:")]
            ok = [r for r in recs if r.get("readable")]
            out[s][c] = {"files": len(recs), "readable": len(ok), "unreadable": len(recs) - len(ok),
                         "distinct_protein_digests": len({r["protein_digest"] for r in ok}),
                         "atom_counts": sorted({r["n_atoms"] for r in ok}),
                         "positions_available_in_any": any(r.get("positions_available") for r in ok),
                         "errors": sorted({r.get("error", "")[:120] for r in recs if not r.get("readable")})[:3]}
        dm = {r["protein_digest"] for k, r in probe.items() if k.startswith(f"{s}:Myr:") and r.get("readable")}
        dn = {r["protein_digest"] for k, r in probe.items() if k.startswith(f"{s}:noMyr:") and r.get("readable")}
        out[s]["protein_digest_sets_equal_between_conditions"] = bool(dm and dm == dn)
    return out


# --- worker (MDAnalysis interpreter): tpr probe ------------------------------------------------------------------------

def _w_tpr(a: dict) -> dict:
    import warnings

    import MDAnalysis as mda
    out = {}
    for key, path in a["files"].items():
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = mda.Universe(path)
                prot = u.select_atoms("protein")
                lines = [f"{r}|{n}|{c:.6f}|{m:.4f}" for r, n, c, m in
                         zip(prot.resnames, prot.names, prot.charges, prot.masses)]
                other = collections.Counter(str(x) for x in u.select_atoms("not protein").residues.resnames)
                try:
                    _ = u.atoms.positions
                    pos = True
                except Exception:
                    pos = False
            out[key] = {"readable": True, "n_atoms": int(u.atoms.n_atoms), "n_protein_atoms": int(prot.n_atoms),
                        "protein_digest": hashlib.sha256("\n".join(lines).encode()).hexdigest(),
                        "non_protein_residues": dict(sorted(other.items())), "positions_available": pos}
        except Exception as e:  # unreadable is a recorded category, not a stop
            out[key] = {"readable": False, "error": f"{type(e).__name__}: {str(e)[:200]}"}
    return out


def probe_tprs(adapter_python: str, files: dict, work: Path) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    pj, oj = work / "audit_tpr_payload.json", work / "audit_tpr_out.json"
    pj.write_text(json.dumps({"files": files, "out": str(oj)}))
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(ROOT / "python"), os.environ.get("PYTHONPATH", "")]))
    p = subprocess.run([str(adapter_python), str(SELF), "_worker", "tpr", str(pj)], env=env, capture_output=True, text=True)
    if p.returncode != 0:
        raise AuditStop("S3", f"tpr probe worker failed: {scrub(p.stderr[-500:])}")
    return json.loads(oj.read_text())


# --- audit -------------------------------------------------------------------------------------------------------------

def members_plan() -> dict:
    """{archive: {key: suffix}} for the amended member list of audit_spec.json."""
    plan = {}
    for c, model in CONDITIONS.items():
        m = {}
        for s, reps in SETS.items():
            for f in (["system.gro"] if s == "Set1" else []) + ["system_em.gro", "topol.top", "production_md.mdp"]:
                m[f"{f}:{s}"] = f"{s}/{f}"
            for r in reps:
                m[f"tpr:{s}:{r}"] = f"{s}/{r}/md.tpr"
        plan[f"{model}.zip"] = m
    return plan


def run_audit(orig_spec: dict, data_dir: Path, extract_dir: Path, probe, spec: dict, log=print) -> dict:
    t0, cpu0 = time.time(), R0.cpu_seconds()
    budget = spec["budget"]
    res = {"id": spec["id"], "post_hoc_amendment": True, "original_frozen_attempt": "ABL-MYR-1 (stopped at G1, unchanged)"}
    try:
        res["archives"] = {}
        for a in ("kinase_kinked_Myr.zip", "kinase_kinked_noMyr.zip"):
            want, p = orig_spec["inputs"][a], Path(data_dir) / a
            if not p.is_file() or p.stat().st_size != want["bytes"] or sha256_file(p) != want["sha256"]:
                raise AuditStop("S1", f"{a} is missing or differs from the frozen size or SHA-256")
            res["archives"][a] = {"bytes": want["bytes"], "sha256": want["sha256"]}
        paths, files, total = {}, {}, 0
        for arch, keys in members_plan().items():
            with zipfile.ZipFile(Path(data_dir) / arch) as zf:
                names = zf.namelist()
                chosen = {}
                for k, suf in keys.items():
                    try:
                        chosen[k] = R0.resolve_member(names, suf)
                    except R0.GateError as e:
                        raise AuditStop("S2", e.message)
                total += sum(zf.getinfo(n).file_size for n in chosen.values())
                if total > budget["extracted_bytes_max"]:
                    raise AuditStop("S4", f"planned extraction {total} B exceeds {budget['extracted_bytes_max']} B")
                for k, member in chosen.items():
                    e = R0.safe_extract(zf, member, Path(extract_dir) / Path(arch).stem)
                    paths[(arch, k)] = e["path"]
                    files[f"{arch}:{k}"] = {"member": member, "bytes": e["bytes"], "sha256": e["sha256"]}
        res["extracted"] = {"files": files, "bytes": sum(f["bytes"] for f in files.values())}
        P = lambda c, k: paths[(f"{CONDITIONS[c]}.zip", k)]  # noqa: E731
        q = {}
        try:
            q["Q1"] = compare_gro_files(P("Myr", "system.gro:Set1"), P("noMyr", "system.gro:Set1"))
            q["Q2"] = compare_gro_files(P("Myr", "system_em.gro:Set2"), P("noMyr", "system_em.gro:Set2"))
            q["Q3"] = compare_gro_files(P("Myr", "system_em.gro:Set1"), P("noMyr", "system_em.gro:Set1"))
            q["Q4"] = {c: compare_gro_files(P(c, "system_em.gro:Set1"), P(c, "system_em.gro:Set2")) for c in CONDITIONS}
            q["Q3_pre_to_post_set1"] = {c: compare_gro_files(P(c, "system.gro:Set1"), P(c, "system_em.gro:Set1"))
                                        for c in CONDITIONS}
            txt = lambda c, k: Path(P(c, k)).read_text()  # noqa: E731
            q["Q5"] = {s: compare_mdp_text(txt("Myr", f"production_md.mdp:{s}"), txt("noMyr", f"production_md.mdp:{s}"))
                       for s in SETS}
            q["Q5_set1_vs_set2_key_differences"] = {
                c: mdp_key_differences(txt(c, "production_md.mdp:Set1"), txt(c, "production_md.mdp:Set2")) for c in CONDITIONS}
            q["Q6"] = {s: compare_top_text(txt("Myr", f"topol.top:{s}"), txt("noMyr", f"topol.top:{s}")) for s in SETS}
            q["Q8"] = {s: {c: {k: R0.mdp_values(txt(c, f"production_md.mdp:{s}")).get(k) for k in INIT_KEYS}
                           for c in CONDITIONS} for s in SETS}
        except (ValueError, IndexError) as e:
            raise AuditStop("S3", f"cannot parse a starting file: {type(e).__name__}: {e}")
        tpr_files = {f"{s}:{c}:{r}": P(c, f"tpr:{s}:{r}") for s, reps in SETS.items() for c in CONDITIONS for r in reps}
        probe_out = probe(tpr_files)
        q["Q7"] = tpr_summary(probe_out)
        q["Q7_per_file"] = probe_out
        res["questions"] = q
        res["derived"] = decide(q)
        res["verdict"] = res["derived"]["verdict"]
        if R0.cpu_seconds() - cpu0 > budget["cpu_seconds_max"]:
            raise AuditStop("S4", "CPU budget exceeded")
    except AuditStop as e:
        res["verdict"] = "audit_incomplete"
        res["stop"] = {"code": e.code, "cause": e.cause}
    res["resources"] = {"cpu_s": R0.cpu_seconds() - cpu0, "wall_s": time.time() - t0, "download_bytes": 0}
    return res


def render_report(res: dict) -> str:
    L = ["# ABL-MYR-1a initialization audit (generated from ABL-MYR-1a_audit_results.json)", "",
         "Post-hoc amendment; not the original frozen test. Engineering audit of starting files only: no trajectory, no "
         "descriptor, no condition difference. Not evidence about Ayaz WT/M472I, HDX, drug response or biological function.", "",
         f"**Verdict: {res['verdict']}**" + (f" (stop {res['stop']['code']}: {res['stop']['cause']})" if "stop" in res else ""), ""]
    if "questions" not in res:
        return "\n".join(L) + "\n"
    q, d = res["questions"], res["derived"]

    def row(name, r):
        if r["category"] == "sequence_mismatch":
            return f"| {name} | sequence_mismatch | | | |"
        return (f"| {name} | {r['category']} | {r['n_protein_atoms']} | {r['max_abs_dx_nm']} | {r['rmsd_nm']} |")
    L += ["| Comparison | Category | Protein atoms | max abs dx (nm) | RMSD (nm) |", "|---|---|---|---|---|",
          row("Q1 Set1 system.gro, Myr vs noMyr", q["Q1"]), row("Q2 Set2 system_em.gro, Myr vs noMyr", q["Q2"]),
          row("Q3 Set1 system_em.gro, Myr vs noMyr (calibration)", q["Q3"]),
          row("Q4 system_em.gro Set1 vs Set2, Myr", q["Q4"]["Myr"]), row("Q4 system_em.gro Set1 vs Set2, noMyr", q["Q4"]["noMyr"]),
          row("Set1 Myr: system.gro vs system_em.gro", q["Q3_pre_to_post_set1"]["Myr"]),
          row("Set1 noMyr: system.gro vs system_em.gro", q["Q3_pre_to_post_set1"]["noMyr"]), ""]
    for s in SETS:
        L.append(f"- {s}: production_md.mdp identical between conditions: {q['Q5'][s]['identical']}; topol.top identical after "
                 f"removing Myr lines: {q['Q6'][s]['identical_after_myr_removal']}")
        for ln in q["Q5"][s]["differing_lines"] + q["Q6"][s]["differing_lines"]:
            L.append(f"  - differing: `{ln}`")
    L += ["", "Run inputs (md.tpr, MDAnalysis topology reader only):", ""]
    for s in SETS:
        for c in CONDITIONS:
            r = q["Q7"][s][c]
            L.append(f"- {s} {c}: {r['readable']} of {r['files']} readable; distinct protein digests {r['distinct_protein_digests']}; "
                     f"positions available: {r['positions_available_in_any']}")
    L += ["", "Initialization stated in production_md.mdp (stated, not verified):", ""]
    for s in SETS:
        for c in CONDITIONS:
            L.append(f"- {s} {c}: " + ", ".join(f"{k}={v}" for k, v in q["Q8"][s][c].items() if v is not None))
    L += ["", f"Derived: Set1 original checks pass = {d['set1_original_checks_pass']}; Set2 parameters pass = "
          f"{d['set2_parameters_pass']}; Set2 start direct evidence = {d['set2_start_direct_evidence']}.", ""]
    return "\n".join(L) + "\n"


# --- freeze / run / report ---------------------------------------------------------------------------------------------

def _hashes(paths) -> dict:
    return {p: sha256_file(ROOT / p) for p in paths}


def cmd_freeze(args) -> int:
    if FREEZE.exists():
        print("freeze exists")
        return 1
    data = ORIG / "data_dryad"
    opened = [str(d.relative_to(ROOT)) for d in (data / "extracted_audit", data / "extracted", ORIG / "work") if d.exists()]
    if opened:
        raise SystemExit(f"STOP: {opened} exist; the audit freeze must precede any extraction")
    rec = {"id": json.loads(SPEC.read_text())["id"], "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "spec_sha256": sha256_file(SPEC), "audit_code": {"path": str(SELF.relative_to(ROOT)), "sha256": sha256_file(SELF)},
           "test": {TEST: sha256_file(ROOT / TEST)},
           "frozen_runner_library_sha256": _hashes(["benchmarks/audit/abl_myr_1/run_abl_myr.py"]),
           "original_artifacts_sha256": _hashes(ORIGINAL_ARTIFACTS),
           "study_members_extracted_before_freeze": [],
           "contents_viewed_before_freeze": "none of Set1 or Set2 (only zip member names and the dataset README)"}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps({k: rec[k] for k in ("frozen_utc", "spec_sha256", "audit_code")}, indent=1))
    return 0


def verify_freeze() -> dict:
    if not FREEZE.exists():
        raise SystemExit("STOP: no audit freeze record")
    f = json.loads(FREEZE.read_text())
    bad = []
    if sha256_file(SPEC) != f["spec_sha256"]:
        bad.append("audit_spec.json")
    if sha256_file(SELF) != f["audit_code"]["sha256"]:
        bad.append("run_init_audit.py")
    bad += [p for p, h in {**f["frozen_runner_library_sha256"], **f["original_artifacts_sha256"]}.items()
            if sha256_file(ROOT / p) != h]
    orig_cf = json.loads((ORIG / "results/ABL-MYR-1_code_freeze.json").read_text())
    if sha256_file(ORIG / "run_abl_myr.py") != orig_cf["runner"]["sha256"]:
        bad.append("run_abl_myr.py differs from the original code freeze")
    if bad:
        raise SystemExit(f"STOP: frozen files differ: {bad}")
    return f


def cmd_run(args) -> int:
    f = verify_freeze()
    orig_spec = json.loads((ORIG / "run_spec.json").read_text())
    spec = json.loads(SPEC.read_text())
    data = ORIG / "data_dryad"
    work = Path(args.work)
    res = run_audit(orig_spec, data, data / "extracted_audit",
                    lambda files: probe_tprs(args.adapter_python, files, work), spec)
    res["audit_spec_sha256"], res["audit_code_sha256"] = f["spec_sha256"], f["audit_code"]["sha256"]
    res["original_artifacts_unchanged"] = True
    Path(args.out).mkdir(parents=True, exist_ok=True)
    (Path(args.out) / RESULTS.name).write_text(json.dumps(scrub(res), indent=1) + "\n")
    print(res["verdict"], json.dumps(res["resources"]))
    return 0


def cmd_report(args) -> int:
    res = json.loads((Path(args.out) / RESULTS.name).read_text())
    (Path(args.out) / "ABL-MYR-1a_audit_report.md").write_text(render_report(res))
    return 0


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "_worker":
        a = json.loads(Path(argv[2]).read_text())
        Path(a["out"]).write_text(json.dumps(_w_tpr(a)))
        return 0
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("freeze")
    r = sub.add_parser("run")
    r.add_argument("--adapter-python", required=True)
    r.add_argument("--out", default=str(HERE / "results"))
    r.add_argument("--work", default=str(HERE / "work"))
    p = sub.add_parser("report")
    p.add_argument("--out", default=str(HERE / "results"))
    args = ap.parse_args(argv)
    return {"freeze": cmd_freeze, "run": cmd_run, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
