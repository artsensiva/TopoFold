"""ABL-MYR-1b comparability check, amendment a1 (see check_amendment_a1.json).

    <python with MDAnalysis> benchmarks/audit/abl_myr_1b/run_comp_check_a1.py freeze
    <python with MDAnalysis> benchmarks/audit/abl_myr_1b/run_comp_check_a1.py run
    <python with MDAnalysis> benchmarks/audit/abl_myr_1b/run_comp_check_a1.py report

The frozen check (run_comp_check.py) stopped at S3 because its index parser rejected repeated group names, which GROMACS
index files may contain. This module imports the frozen check unchanged and replaces only parse_ndx and tc_group_check
while run_check runs; every other rule, constant and criterion of check_spec.json is applied as frozen.
"""

from __future__ import annotations

import argparse
import contextlib
import datetime as dt
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_comp_check as B  # noqa: E402  (frozen check, unchanged)

AMENDMENT = HERE / "check_amendment_a1.json"
FREEZE = HERE / "results/ABL-MYR-1b_check_a1_freeze.json"
RESULTS = "ABL-MYR-1b_check_a1_results.json"
CA_NPZ = "ABL-MYR-1b_production_start_ca_a1.npz"
SELF = Path(__file__).resolve()
TEST = "tests/test_abl_myr_1b_check_a1.py"
ATTEMPT1 = ["benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_results_attempt1.json",
            "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_report_attempt1.md",
            "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_production_start_ca_attempt1.npz"]
AMBIGUOUS = None
_FROZEN_TC_GROUP_CHECK = B.tc_group_check     # bound once, so the patch cannot recurse


def parse_ndx(text: str) -> dict:
    """Index groups; a repeated name keeps its atom list if every occurrence is identical, else it is ambiguous (None)."""
    occ, order, cur = {}, [], None
    for ln in text.splitlines():
        s = ln.strip()
        m = re.match(r"^\[\s*(.+?)\s*\]$", s)
        if m:
            cur = [m.group(1), []]
            order.append(cur)
        elif s and cur is not None:
            cur[1] += [int(x) - 1 for x in s.split()]
        elif s:
            raise ValueError("index entries before the first group")
    for name, idx in order:
        occ.setdefault(name, []).append(np.array(idx, dtype=np.int64))
    out = {}
    for name, arrs in occ.items():
        out[name] = arrs[0] if all(np.array_equal(arrs[0], a) for a in arrs[1:]) else AMBIGUOUS
    return out


def tc_group_check(groups, tc_names, atom_resnames, atom_resindex, protein_resnames, condition):
    amb = [g for g in tc_names if g in groups and groups[g] is AMBIGUOUS]
    if amb:
        return {"tc_grps": tc_names, "groups": {}, "ok": False,
                "problems": [f"ambiguous group {g!r}: repeated in index.ndx with different atoms" for g in amb]}
    clean = {k: v for k, v in groups.items() if v is not AMBIGUOUS}
    return _FROZEN_TC_GROUP_CHECK(clean, tc_names, atom_resnames, atom_resindex, protein_resnames, condition)


@contextlib.contextmanager
def amended():
    saved = (B.parse_ndx, B.tc_group_check)
    B.parse_ndx, B.tc_group_check = parse_ndx, tc_group_check
    try:
        yield
    finally:
        B.parse_ndx, B.tc_group_check = saved


def run_check(*args, **kw):
    with amended():
        return B.run_check(*args, **kw)


def cmd_freeze(args) -> int:
    if FREEZE.exists():
        print("freeze exists")
        return 1
    if (B.ORIG / "data_dryad/extracted_1b_a1").exists():
        raise SystemExit("STOP: extracted_1b_a1 exists; the a1 freeze must precede its run")
    B.verify_freeze()
    rec = {"id": "ABL-MYR-1b-CHECK-a1", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "amendment_sha256": B.R0.sha256_file(AMENDMENT), "code": {"path": str(SELF.relative_to(B.ROOT)), "sha256": B.R0.sha256_file(SELF)},
           "test": {TEST: B._h(TEST)}, "frozen_check": {"spec_sha256": B.R0.sha256_file(B.SPEC), "code_sha256": B.R0.sha256_file(B.SELF),
                                                      "freeze_record_sha256": B.R0.sha256_file(B.FREEZE)},
           "attempt1_sha256": {p: B._h(p) for p in ATTEMPT1}}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps({k: rec[k] for k in ("frozen_utc", "amendment_sha256", "code")}, indent=1))
    return 0


def verify() -> dict:
    B.verify_freeze()
    f = json.loads(FREEZE.read_text())
    bad = [n for n, ok in (("amendment", B.R0.sha256_file(AMENDMENT) == f["amendment_sha256"]),
                           ("a1 code", B.R0.sha256_file(SELF) == f["code"]["sha256"]),
                           ("frozen check freeze record", B.R0.sha256_file(B.FREEZE) == f["frozen_check"]["freeze_record_sha256"]))
           if not ok]
    bad += [p for p, h in f["attempt1_sha256"].items() if B._h(p) != h]
    if bad:
        raise SystemExit(f"STOP: files differ from the a1 freeze: {bad}")
    return f


def cmd_run(args) -> int:
    f = verify()
    spec = json.loads(B.SPEC.read_text())
    res, ca = run_check(json.loads((B.ORIG / "run_spec.json").read_text()), B.ORIG / "data_dryad",
                        B.ORIG / "data_dryad/extracted_1b_a1", spec)
    res["amendment"] = {"id": "ABL-MYR-1b-CHECK-a1", "amendment_sha256": f["amendment_sha256"], "code_sha256": f["code"]["sha256"],
                        "frozen_check_spec_sha256": f["frozen_check"]["spec_sha256"],
                        "frozen_check_code_sha256": f["frozen_check"]["code_sha256"],
                        "original_attempt": "results/ABL-MYR-1b_check_results_attempt1.json (S3 stop)"}
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    np.savez(out / CA_NPZ, **{k.replace(":", "_"): v for k, v in ca.items() if v is not None})
    res["production_start_ca_npz"] = {"file": CA_NPZ, "sha256": B.R0.sha256_file(out / CA_NPZ), "units": "Å", "atoms": "protein Cα, whole"}
    (out / RESULTS).write_text(json.dumps(B.R0.scrub(res), indent=1, default=B.R0._json_default) + "\n")
    print(res["outcome"], json.dumps(res["resources"]))
    return 0


def render_report(res: dict) -> str:
    head = ("> Amendment a1 of the frozen check (check_amendment_a1.json): index files with repeated group names are accepted. "
            "The original attempt stopped at S3 and is kept as results/ABL-MYR-1b_check_results_attempt1.json.\n\n")
    return head + B.render_report(res)


def cmd_report(args) -> int:
    res = json.loads((Path(args.out) / RESULTS).read_text())
    (Path(args.out) / "ABL-MYR-1b_check_a1_report.md").write_text(render_report(res))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("freeze")
    for name in ("run", "report"):
        p = sub.add_parser(name)
        p.add_argument("--out", default=str(HERE / "results"))
    args = ap.parse_args(argv)
    return {"freeze": cmd_freeze, "run": cmd_run, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
