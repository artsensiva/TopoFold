"""Re-check a structpair result against its manifest.json (standard library only).

    PYTHONPATH=python python -m topofold_inputs.structpair_verify RESULT_DIR [--json OUT.json]

Recomputes the SHA-256 of every file in `outputs_sha256`, of both input structures and of the spec, and compares each with the
manifest. A path in the manifest is used as recorded; a relative one is read from the current directory. Prints one line per
file and the counts. Exit code 0 only if every listed file exists and matches. manifest.json does not hash itself, so a
changed manifest is not detected by this check; compare its hash with a copy kept elsewhere if that matters.
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verify(result_dir) -> dict:
    d = Path(result_dir)
    m = json.loads((d / "manifest.json").read_text())
    todo = [("output", name, d / name, h) for name, h in m["outputs_sha256"].items()]
    todo += [("input " + k, v["path"], Path(v["path"]), v["sha256"]) for k, v in m["inputs"].items()]
    todo.append(("spec", m["spec"]["path"], Path(m["spec"]["path"]), m["spec"]["sha256"]))
    checks = []
    for kind, name, path, expected in todo:
        if not path.is_file():
            checks.append({"kind": kind, "path": str(name), "status": "missing", "expected": expected, "actual": None})
            continue
        actual = _sha256(path)
        checks.append({"kind": kind, "path": str(name), "status": "ok" if actual == expected else "mismatch",
                       "expected": expected, "actual": actual})
    n_ok = sum(c["status"] == "ok" for c in checks)
    return {"result": str(d), "schema": m.get("schema"), "checks": checks, "n_checked": len(checks), "n_ok": n_ok,
            "n_mismatch": sum(c["status"] == "mismatch" for c in checks), "n_missing": sum(c["status"] == "missing" for c in checks),
            "all_ok": n_ok == len(checks)}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("result_dir")
    ap.add_argument("--json", help="write the full record to this file")
    a = ap.parse_args(argv)
    try:
        r = verify(a.result_dir)
    except (OSError, ValueError, KeyError) as e:
        print(f"cannot read the manifest of {a.result_dir}: {e}", file=sys.stderr)
        return 2
    for c in r["checks"]:
        print(f"{c['status']:8s} {c['kind']:8s} {c['path']}")
    print(f"{r['n_ok']} of {r['n_checked']} files match the manifest ({r['n_mismatch']} mismatch, {r['n_missing']} missing); schema {r['schema']}")
    if a.json:
        Path(a.json).write_text(json.dumps(r, indent=1) + "\n")
    return 0 if r["all_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
