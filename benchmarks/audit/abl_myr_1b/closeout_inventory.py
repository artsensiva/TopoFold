"""ABL-MYR-1b closeout: provenance inventory of the executed comparison (nothing is rerun or recomputed).

    PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1b/closeout_inventory.py [--copy-dir ~/Downloads]

Hashes (SHA-256, full) what is present locally and records for each artifact whether it is committed to Git, present only
locally (git-ignored) or absent. It checks:
  * the Dryad archives against the source digests in debuhr_abl_intake/intake_record.json;
  * the extracted archive members against the G1 member digests recorded in results/ABL-MYR-1b_results.json;
  * every (path, SHA-256) pair in the ABL-MYR-1 / 1a / 1b freeze records against the files as they now are, and the
    frozen wrapper's own lineage check (run_abl_myr_1b.lineage_problems);
  * optional further copies of the archives (--copy-dir), stating whether they sit on the same device.
No archive is opened, no trajectory is read and no comparison number is recomputed.
Writes results/ABL-MYR-1b_closeout_inventory.json and results/ABL-MYR-1b_closeout_SHA256SUMS.txt.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DRYAD = ROOT / "benchmarks/audit/abl_myr_1/data_dryad"
INTAKE = ROOT / "benchmarks/audit/debuhr_abl_intake/intake_record.json"
RESULTS = HERE / "results/ABL-MYR-1b_results.json"
OUT = HERE / "results/ABL-MYR-1b_closeout_inventory.json"
SUMS = HERE / "results/ABL-MYR-1b_closeout_SHA256SUMS.txt"
ADAPTER_PY = ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python"
CORE_PY = ROOT / ".venv/bin/python"
GMX_DEFAULT = Path.home() / "tools/gromacs-a2/env/bin/gmx"
USED_ARCHIVES = ("PCA.zip", "kinase_kinked_Myr.zip", "kinase_kinked_noMyr.zip")
FREEZES = ["benchmarks/audit/abl_myr_1/results/ABL-MYR-1_freeze.json",
           "benchmarks/audit/abl_myr_1/results/ABL-MYR-1_code_freeze.json",
           "benchmarks/audit/abl_myr_1a/results/ABL-MYR-1a_audit_freeze.json",
           "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_freeze.json",
           "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_a1_freeze.json",
           "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_a2_freeze.json",
           "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_protocol_freeze.json",
           "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_protocol_amendment_a1_freeze.json",
           "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_wrapper_freeze.json"]
COMMITTED_DIRS = ["benchmarks/audit/abl_myr_1", "benchmarks/audit/abl_myr_1a", "benchmarks/audit/abl_myr_1b",
                  "benchmarks/audit/debuhr_abl_intake", "benchmarks/audit/abl_prereg/reference"]
COMMITTED_TESTS = r"^tests/test_abl_(myr|prereg)"
LOCAL_DIRS = {  # git-ignored derived data: role, how it can be regenerated
    "benchmarks/audit/abl_myr_1/data_dryad/extracted":
        ("members extracted by the ABL-MYR-1b comparison (G1)", "run_abl_myr_1b.py run (re-extracts from the 3 archives)"),
    "benchmarks/audit/abl_myr_1b/work":
        ("comparison intermediates (Cα coordinates, φ/ψ, RMSF, per-run features, worker I/O)", "run_abl_myr_1b.py run"),
    "benchmarks/audit/abl_myr_1/data_dryad/extracted_audit":
        ("members extracted by the ABL-MYR-1a initialisation audit", "run_init_audit.py run"),
    "benchmarks/audit/abl_myr_1/data_dryad/extracted_1b":
        ("members extracted by the 1b check, attempt 1 (S3 stop)", "run_comp_check.py run"),
    "benchmarks/audit/abl_myr_1/data_dryad/extracted_1b_a1":
        ("members extracted by check amendment a1", "run_comp_check_a1.py run"),
    "benchmarks/audit/abl_myr_1/data_dryad/extracted_1b_a2":
        ("members extracted by check amendment a2 (read with GROMACS)", "run_comp_check_a2.py run"),
    "benchmarks/audit/abl_myr_1a/work": ("ABL-MYR-1a audit intermediates", "run_init_audit.py run"),
}
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def sha256_file(p: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def rel(p: Path, root: Path | None = None) -> str:
    return str(Path(os.path.abspath(p)).relative_to(root or ROOT))   # abspath: keep venv symlinks inside the repository


def tracked_files(root: Path = ROOT, others: bool = False) -> set:
    """Tracked files; with others=True the untracked, not ignored files instead (those a normal commit would add)."""
    cmd = ["git", "ls-files", "-z"] + (["--others", "--exclude-standard"] if others else [])
    out = subprocess.run(cmd, cwd=root, capture_output=True, check=True).stdout.decode()
    return {p for p in out.split("\0") if p}


def hash_pairs(obj, root: Path = ROOT) -> list:
    """Every (repo path, SHA-256) pair in a freeze record: {path: hex} entries and {'path': .., 'sha256': ..} objects."""
    out = []
    if isinstance(obj, dict):
        if isinstance(obj.get("path"), str) and isinstance(obj.get("sha256"), str) and HEX64.match(obj["sha256"]) \
                and "/" in obj["path"]:
            out.append((obj["path"], obj["sha256"]))
        for k, v in obj.items():
            if isinstance(v, str) and HEX64.match(v) and "/" in k and (root / k).suffix:
                out.append((k, v))
            else:
                out += hash_pairs(v, root)
    elif isinstance(obj, list):
        for v in obj:
            out += hash_pairs(v, root)
    return out


def check_freezes(freezes: list, root: Path = ROOT) -> dict:
    per, bad = {}, []
    for fz in freezes:
        pairs = sorted(set(hash_pairs(json.loads((root / fz).read_text()), root)))
        rows = []
        for p, h in pairs:
            f = root / p
            got = sha256_file(f) if f.is_file() else None
            rows.append({"path": p, "recorded": h, "matches": got == h, "present": f.is_file()})
            if got != h:
                bad.append({"freeze": fz, "path": p, "recorded": h, "now": got})
        per[fz] = {"sha256": sha256_file(root / fz), "pairs_checked": len(rows),
                   "pairs_matching": sum(r["matches"] for r in rows)}
    return {"records": per, "pairs_checked": sum(v["pairs_checked"] for v in per.values()),
            "pairs_matching": sum(v["pairs_matching"] for v in per.values()), "mismatches": bad}


def archives(intake: dict, dryad: Path, copy_dir: Path | None) -> list:
    rows = []
    for f in intake["raw_inventory"]["files"]:
        p = dryad / f["path"]
        row = {"name": f["path"], "bytes_source": f["bytes"], "sha256_source": f["sha256_source"],
               "used_by_abl_myr_1b_comparison": f["path"] in USED_ARCHIVES, "local_path": rel(p) if p.exists() else None}
        if p.is_file():
            row.update(bytes_local=p.stat().st_size, sha256_local=sha256_file(p))
            row["matches_source"] = row["sha256_local"] == f["sha256_source"] and row["bytes_local"] == f["bytes"]
        else:
            row.update(bytes_local=None, sha256_local=None, matches_source=None)
        if copy_dir is not None and (copy_dir / f["path"]).is_file():
            c = copy_dir / f["path"]
            row["further_copy"] = {"location": "~/" + str(c.resolve().relative_to(Path.home())) if c.resolve().is_relative_to(Path.home()) else str(c),
                                   "sha256": sha256_file(c), "same_device_as_local": p.is_file() and os.stat(c).st_dev == os.stat(p).st_dev}
        rows.append(row)
    return rows


def members(results: dict, ext_root: Path) -> dict:
    """G1 member digests (recorded during the run) against the extracted files now on disk."""
    files = results["gates"]["G1"]["files"]
    rows, bad, missing = [], [], []
    for key, f in sorted(files.items()):
        archive = key.split(":")[0]
        p = ext_root / archive.removesuffix(".zip") / f["member"]
        cand = [p, ext_root / f["member"]]
        q = next((c for c in cand if c.is_file()), None)
        now = sha256_file(q) if q else None
        rows.append({"key": key, "member": f["member"], "bytes": f["bytes"], "sha256_g1": f["sha256"],
                     "local": rel(q) if q else None, "matches": now == f["sha256"]})
        (missing if q is None else bad if now != f["sha256"] else []).append(key)
    return {"recorded_in_g1": len(files), "present": len(files) - len(missing), "matching": sum(r["matches"] for r in rows),
            "missing": missing, "mismatching": bad, "files": rows}


def reconcile_extracted(mem: dict, ext_tree: dict) -> dict:
    """Why the extracted directory's file count and bytes differ from G1's member list (duplicates, reader caches)."""
    seen = {}
    for f in mem["files"]:
        if f["local"]:
            seen.setdefault(f["local"], []).append(f["key"])
    other = [e for e in ext_tree["entries"] if e["path"] not in seen]
    kinds = {}
    for e in other:
        name = e["path"].rsplit("/", 1)[-1]
        kind = "MDAnalysis XTC offset cache" if name.startswith(".") and "_offsets." in name else "other"
        kinds.setdefault(kind, {"files": 0, "bytes": 0})
        kinds[kind]["files"] += 1
        kinds[kind]["bytes"] += e["bytes"]
    return {"g1_entries": len(mem["files"]), "unique_member_files": len(seen),
            "files_listed_under_several_g1_keys": {p: k for p, k in seen.items() if len(k) > 1},
            "files_not_in_g1": kinds}


def tree(d: Path) -> dict:
    files = sorted(p for p in d.rglob("*") if p.is_file()) if d.is_dir() else []
    rows = [{"path": rel(p), "bytes": p.stat().st_size, "sha256": sha256_file(p)} for p in files]
    digest = hashlib.sha256("".join(f"{r['sha256']}  {r['path']}\n" for r in rows).encode()).hexdigest()
    return {"present": d.is_dir(), "files": len(rows), "bytes": sum(r["bytes"] for r in rows),
            "manifest_sha256": digest if rows else None, "entries": rows}


def committed(tracked: set, untracked: set = frozenset(), root: Path = ROOT, exclude: tuple = ()) -> list:
    """Lineage files in Git, plus untracked files in the same places that are not ignored (added by the closeout commit)."""
    keep = lambda p: any(p.startswith(d + "/") for d in COMMITTED_DIRS) or re.match(COMMITTED_TESTS, p)  # noqa: E731
    paths = sorted(p for p in tracked | set(untracked) if keep(p) and p not in exclude and "__pycache__" not in p)
    return [{"path": p, "bytes": (root / p).stat().st_size, "sha256": sha256_file(root / p), "tracked_at_inventory": p in tracked}
            for p in paths if (root / p).is_file()]


def wrapper_lineage() -> list:
    """The frozen wrapper's own G0 lineage check (protocol freeze lineage by name, protocol, amendment freeze)."""
    sys.path.insert(0, str(HERE))
    import run_abl_myr_1b as W1b
    return W1b.lineage_problems()


def _run(cmd: list) -> str | None:
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=120, check=True).stdout
    except (OSError, subprocess.SubprocessError):
        return None


def python_env(py: Path) -> dict:
    if not py.exists():
        return {"present": False}
    code = ("import sys,json,importlib.metadata as m;"
            "print(json.dumps({'python':sys.version.split()[0],"
            "'distributions':sorted(f\"{d.metadata['Name']}=={d.version}\" for d in m.distributions())}))")
    out = _run([str(py), "-I", "-c", code])
    return {"present": True, "path": rel(py), **(json.loads(out) if out else {"error": "could not query"})}


def environment(gmx: Path) -> dict:
    so = sorted((ROOT / ".venv/lib").glob("python*/site-packages/topofold/topofold*.so"))
    ver = _run([str(gmx), "--version"]) if gmx.exists() else None
    gv = {k: m.group(1).strip() for k in ("GROMACS version", "Precision", "MPI library")
          if ver and (m := re.search(rf"^{k}:\s*(.+)$", ver, re.M))}
    head = _run(["git", "-C", str(ROOT), "rev-parse", "HEAD"])
    crates = _run(["git", "-C", str(ROOT), "log", "-1", "--format=%H %cI", "--", "crates", "Cargo.lock"])
    return {
        "os": platform.platform(), "machine": platform.machine(),
        "git_head_at_inventory": head.strip() if head else None,
        "core_python": python_env(CORE_PY), "adapter_python": python_env(ADAPTER_PY),
        "compiled_core": [{"path": rel(p), "sha256": sha256_file(p),
                           "mtime_utc": dt.datetime.fromtimestamp(p.stat().st_mtime, dt.timezone.utc).isoformat()} for p in so],
        "last_commit_touching_rust_sources": crates.strip() if crates else None,
        "compiled_core_note": "built locally with maturin and git-ignored; that it was built from the commit above is not "
                              "recorded by the build and is not verified here",
        "gromacs": {"present": gmx.exists(), "location": "outside the repository, under the user's home (TOPOFOLD_GMX)",
                    **gv, "launcher_sha256": sha256_file(gmx) if gmx.exists() else None,
                    "used_by": "check amendment a2 only (run-input reading); the comparison itself does not call GROMACS",
                    "acquisition_record": "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_gromacs_acquisition.json"},
    }


SOURCE_NOTE = ("Dryad API file downloads need a bearer token and the web download passes a bot challenge; the owner downloaded "
               "the archives in a browser (D57). Recovery therefore needs the same manual step.")

COMMANDS = [
    "PY=benchmarks/audit/wpd_local_retrieval/build/venv/bin/python   # Python with MDAnalysis 2.10.0",
    "sha256sum -c benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_closeout_SHA256SUMS.txt   # committed files and archives",
    "PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1b/run_abl_myr_1b.py run --adapter-python $PY   # comparison",
    "PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1b/run_abl_myr_1b.py report",
    "python3 benchmarks/audit/abl_myr_1b/summarize_results_1b.py",
]

RECOVERY = [
    "1. Check out scientific-rebuild at or after the closeout commit; build the core (maturin develop --release) and the "
    "MDAnalysis environment (benchmarks/audit/wpd_local_retrieval).",
    "2. Download PCA.zip, kinase_kinked_Myr.zip and kinase_kinked_noMyr.zip from doi:10.5061/dryad.9cnp5hqnx (version 5) "
    "in a browser into benchmarks/audit/abl_myr_1/data_dryad/.",
    "3. Verify them: sha256sum -c results/ABL-MYR-1b_closeout_SHA256SUMS.txt (the archive lines must match; they equal the "
    "Dryad source digests).",
    "4. To inspect the recorded result, nothing needs to be rerun: results/ABL-MYR-1b_results.json, the window CSVs, the "
    "report and the residue summary are committed.",
    "5. To regenerate the intermediates for method development, run the comparison into a new output directory "
    "(--out and --work pointing outside results/ and work/), so that the recorded result files are not overwritten; "
    "a fresh output is a new run, not the recorded one, and byte-identity with the recorded intermediates has not been tested.",
    "6. The a2 run-input check additionally needs GROMACS 2024.5 (results/ABL-MYR-1b_gromacs_acquisition.json gives the "
    "installer, package digests and commands).",
]


def build(copy_dir: Path | None, gmx: Path) -> dict:
    intake = json.loads(INTAKE.read_text())
    res = json.loads(RESULTS.read_text())
    tracked = tracked_files()
    arch = archives(intake, DRYAD, copy_dir)
    mem = members(res, DRYAD / "extracted")
    fz = check_freezes(FREEZES)
    local = {d: {"role": r, "regenerate_with": how, "git_ignored": d not in tracked and
                 subprocess.run(["git", "check-ignore", "-q", d], cwd=ROOT).returncode == 0, **tree(ROOT / d)}
             for d, (r, how) in LOCAL_DIRS.items()}
    com = committed(tracked, tracked_files(others=True), exclude=(rel(OUT), rel(SUMS)))
    lineage = wrapper_lineage()
    used = [a for a in arch if a["used_by_abl_myr_1b_comparison"]]
    copies = [a for a in used if "further_copy" in a]
    gaps = []
    if any(a["sha256_local"] is None for a in used):
        gaps.append("an archive used by the comparison is not present locally")
    if not copies or any(c["further_copy"]["same_device_as_local"] for c in copies):
        gaps.append("no independent backup of the three archives is verified: any further copy found is on the same device "
                    "as the working copy; the recoverable source is the Dryad record (manual browser download)")
    gaps.append("the git-ignored derived data (extracted members, work/ intermediates) exist only on this machine; they are "
                "hashed here and regenerable from the archives with the committed code, but no regeneration was run")
    gaps.append("the compiled core (.so) and both Python environments are local builds; their versions and hashes are "
                "recorded, not stored")
    gaps.append("GROMACS 2024.5 is installed outside the repository; only its acquisition record is committed")
    gaps.append("Set2/system.gro is not published in the Dryad record (D58); it cannot be recovered from the source")
    if lineage:
        gaps.append(f"the wrapper lineage check reports {lineage}")
    if mem["missing"] or mem["mismatching"]:
        gaps.append(f"extracted members missing {len(mem['missing'])}, mismatching {len(mem['mismatching'])} against G1")
    return {
        "id": "ABL-MYR-1b-CLOSEOUT",
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "scope": "provenance inventory only; no comparison, threshold or frozen file was rerun or changed",
        "source": {"dataset_doi": "10.5061/dryad.9cnp5hqnx", "url": intake["source"]["url"],
                   "version": intake["source"]["version_or_release"], "licence": intake["licence"]["data_licence_as_stated"],
                   "study_doi": "10.7554/eLife.85216", "access_note": SOURCE_NOTE},
        "result_files": {p: sha256_file(HERE / "results" / p) for p in
                         ("ABL-MYR-1b_results.json", "ABL-MYR-1b_report.md", "ABL-MYR-1b_residue_summary.json",
                          "ABL-MYR-1b_windows_Set1.csv", "ABL-MYR-1b_windows_Set2.csv", "ABL-MYR-1b_wrapper_freeze.json")},
        "result_records": {"outcome_label": res["outcome"].get("label"), "status": res["status"],
                           "wrapper_sha256": res["wrapper_sha256"], "amendment_sha256": res["amendment_sha256"],
                           "g0_archives": res["gates"]["G0"]},
        "archives": arch,
        "g1_members": mem,
        "extracted_vs_g1": reconcile_extracted(mem, local["benchmarks/audit/abl_myr_1/data_dryad/extracted"]),
        "freeze_checks": fz,
        "wrapper_lineage_problems": lineage,
        "committed": {"files": len(com), "bytes": sum(c["bytes"] for c in com), "entries": com},
        "local_only": local,
        "environment": environment(gmx),
        "commands": COMMANDS,
        "recovery_steps": RECOVERY,
        "gaps": gaps,
    }


def sums_text(inv: dict) -> str:
    lines = [f"{c['sha256']}  {c['path']}" for c in inv["committed"]["entries"]]
    lines += [f"{a['sha256_source']}  benchmarks/audit/abl_myr_1/data_dryad/{a['name']}" for a in inv["archives"]
              if a["used_by_abl_myr_1b_comparison"]]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copy-dir", default=None, help="a directory that may hold further copies of the archives")
    ap.add_argument("--gmx", default=os.environ.get("TOPOFOLD_GMX", str(GMX_DEFAULT)))
    args = ap.parse_args(argv)
    inv = build(Path(args.copy_dir).expanduser() if args.copy_dir else None, Path(args.gmx))
    OUT.write_text(json.dumps(inv, indent=1, ensure_ascii=False) + "\n")
    SUMS.write_text(sums_text(inv))
    a = [x for x in inv["archives"] if x["used_by_abl_myr_1b_comparison"]]
    print(json.dumps({"archives_used_matching_source": sum(bool(x["matches_source"]) for x in a),
                      "g1_members": {k: inv["g1_members"][k] for k in ("recorded_in_g1", "present", "matching")},
                      "freeze_pairs": [inv["freeze_checks"]["pairs_checked"], inv["freeze_checks"]["pairs_matching"]],
                      "wrapper_lineage_problems": inv["wrapper_lineage_problems"],
                      "committed_files": inv["committed"]["files"], "gaps": len(inv["gaps"])}, indent=1))
    ok = not (inv["freeze_checks"]["mismatches"] or inv["g1_members"]["mismatching"] or inv["wrapper_lineage_problems"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
