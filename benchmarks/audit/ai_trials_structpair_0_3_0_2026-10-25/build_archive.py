"""Build the compact archive of the two AI trials of structpair 0.3.0 (commit 19b5268) and hash-verify it.

    .venv/bin/python benchmarks/audit/ai_trials_structpair_0_3_0_2026-10-25/build_archive.py [--verify-only]

The two trial directories are read only. Their environments, clones, build output and `work/` folders are not copied or hashed.
`inventory.json` records the SHA-256, size and original path of every other file of each trial (including all screenshots,
results and the corrupted copies); `archive/` holds a compact copy of the reports, specs, manifests, summaries, evidence records
and the screenshots of the viewer finding. `--verify-only` re-hashes the archive copies and the originals against the inventory.
"""

import hashlib
import json
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRIALS = {
    "trial1": {
        "agent": "GPT-6 family, Codex desktop (exact variant not exposed)",
        "root": "/home/artem/Documents/Codex/2026-10-25/files-pasted-by-the-user-you/outputs/topofold-independent-trial",
        "hash_dirs": ["evidence", "screenshots", "results", "specs", "logs", "integrity", "inputs", "kras_struct_cmp", "kras_struct_g12d"],
        "hash_files": ["REPORT.md", "run.py", "audit_inputs.py", "interpret_results.py", "verify_manifests.py", "finalize_report.py"],
        "copy": ["REPORT.md", "specs", "evidence/preregistration.md", "evidence/independent-interpretation.md", "evidence/input-audit.json",
                 "evidence/manifest-verification.json", "evidence/frozen-evidence-verification.json", "evidence/tamper-test.json",
                 "evidence/timing.json", "evidence/departures.md", "evidence/artifact-sha256.json", "screenshots/CAPTIONS.md",
                 "logs/COMMANDS.md", "results/abl/manifest.json", "results/abl/summary.json", "results/kras2/manifest.json",
                 "results/kras2/summary.json"],
    },
    "trial2": {
        "agent": "Claude Opus 5.5 (claude-opus-5-5), Claude desktop app",
        "root": "/home/artem/.config/Claude/scratch-workspaces/e9fd6719-de25-4a0b-8f80-7ee2955de949/78f8b148-a578-47e6-b26f-1b088f117afe/scratch-2026-10-25-04801e/trial",
        "hash_dirs": ["notes", "screenshots", "outputs", "specs", "logs", "tamper", "scripts", "inputs"],
        "hash_files": ["REPORT.md"],
        "copy": ["REPORT.md", "specs", "notes", "screenshots/CAPTIONS.md", "screenshots/24_kras_focus_desync_top.jpg",
                 "screenshots/25_kras_focus_desync_settled.jpg", "screenshots/26_kras_focus_resync.jpg",
                 "screenshots/27_repro_focus39.jpg", "screenshots/28_repro_desync.jpg", "logs", "scripts",
                 "outputs/abl_1OPJ_AB_vs_1OPK_A/manifest.json", "outputs/abl_1OPJ_AB_vs_1OPK_A/summary.json",
                 "outputs/kras2/manifest.json", "outputs/kras2/summary.json"],
    },
}


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def files_of(root: Path, spec: dict):
    out = []
    for d in spec["hash_dirs"]:
        for p in sorted((root / d).rglob("*")):
            if p.is_file() and "work" not in p.relative_to(root).parts and "__pycache__" not in p.parts:
                out.append(p)
    out += [root / f for f in spec["hash_files"] if (root / f).is_file()]
    return sorted(set(out))


def build():
    inv, arch = {}, HERE / "archive"
    for name, spec in TRIALS.items():
        root = Path(spec["root"])
        rows = [{"path": str(p.relative_to(root)), "bytes": p.stat().st_size, "sha256": sha(p)} for p in files_of(root, spec)]
        inv[name] = {"agent": spec["agent"], "root": spec["root"], "n_files": len(rows), "files": rows}
        for rel in spec["copy"]:
            src = root / rel
            for f in ([src] if src.is_file() else sorted(q for q in src.rglob("*") if q.is_file())):
                dst = arch / name / f.relative_to(root)
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(f, dst)
    (HERE / "inventory.json").write_text(json.dumps(inv, indent=1) + "\n")


def verify() -> int:
    inv = json.loads((HERE / "inventory.json").read_text())
    res = {}
    for name, t in inv.items():
        bad_orig, n_arch, bad_arch = [], 0, []
        by = {r["path"]: r["sha256"] for r in t["files"]}
        for r in t["files"]:
            p = Path(t["root"]) / r["path"]
            if not p.is_file() or sha(p) != r["sha256"]:
                bad_orig.append(r["path"])
        for f in sorted((HERE / "archive" / name).rglob("*")):
            if f.is_file():
                rel = str(f.relative_to(HERE / "archive" / name))
                n_arch += 1
                if rel in by and sha(f) != by[rel]:
                    bad_arch.append(rel)
                elif rel not in by:
                    bad_arch.append(rel + " (not in inventory)")
        res[name] = {"originals_checked": len(t["files"]), "originals_changed_or_missing": bad_orig,
                     "archive_files_checked": n_arch, "archive_mismatch": bad_arch}
    print(json.dumps(res, indent=1))
    return 0 if all(not v["originals_changed_or_missing"] and not v["archive_mismatch"] for v in res.values()) else 1


if __name__ == "__main__":
    if "--verify-only" not in sys.argv:
        build()
    sys.exit(verify())
