"""ABL-MYR-1b comparability check, amendment a2 (see check_amendment_a2.json): independent reading of the stored box, positions
and velocities of the 40 production md.tpr with GROMACS, cross-checked against the earlier reader.

    TOPOFOLD_GMX=<path to gmx> <python with MDAnalysis> benchmarks/audit/abl_myr_1b/run_comp_check_a2.py freeze
    TOPOFOLD_GMX=<path to gmx> <python with MDAnalysis> benchmarks/audit/abl_myr_1b/run_comp_check_a2.py run
    <python with MDAnalysis> benchmarks/audit/abl_myr_1b/run_comp_check_a2.py report

The frozen check (run_comp_check.py) and amendment a1 (index groups) are imported unchanged. This module only supplies a
different md.tpr loader, whose record carries the GROMACS box, positions and velocities, and records the GROMACS
provenance, the reader agreement and descriptive run-input settings. One bounded attempt: no further workaround.
"""

from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import run_comp_check as B  # noqa: E402  (frozen check, unchanged)
import run_comp_check_a1 as A1  # noqa: E402  (amendment a1, unchanged)

AMENDMENT = HERE / "check_amendment_a2.json"
FREEZE = HERE / "results/ABL-MYR-1b_check_a2_freeze.json"
RESULTS = "ABL-MYR-1b_check_a2_results.json"
CA_NPZ = "ABL-MYR-1b_production_start_ca_a2.npz"
SELF = Path(__file__).resolve()
TEST = "tests/test_abl_myr_1b_check_a2.py"
ACQUISITION = "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_gromacs_acquisition.json"
A1_RECORDS = ["benchmarks/audit/abl_myr_1b/check_amendment_a1.json", "benchmarks/audit/abl_myr_1b/run_comp_check_a1.py",
              "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_a1_freeze.json",
              "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_a1_results.json",
              "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_a1_report.md",
              "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_production_start_ca_a1.npz",
              "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_tpr_header_diagnostic.json"]
POS_TOL_NM, VEL_TOL_NM_PS, BOX_TOL_NM = 1e-3, 1e-3, 1e-4      # 1e-3 nm = 0.01 Å
DUMP_TIMEOUT_S = 600
IR_KEYS = ("integrator", "dt", "nsteps", "tcoupl", "nsttcouple", "pcoupl", "pcoupltype", "tau-p", "constraints", "continuation",
           "ld-seed", "refcoord-scaling", "comm-mode", "cutoff-scheme", "coulombtype", "vdw-modifier")


class GmxUnavailable(Exception):
    pass


# --- GROMACS ------------------------------------------------------------------------------------------------------------

def gmx_executable(env=None):
    env = os.environ if env is None else env
    p = env.get("TOPOFOLD_GMX")
    if p:
        return p if Path(p).is_file() else None
    return shutil.which("gmx")


def gmx_provenance(gmx: str) -> dict:
    """Version, precision, package record and executable hashes, with the install prefix scrubbed to <env>."""
    p = subprocess.run([gmx, "--version"], capture_output=True, text=True, timeout=120)
    if p.returncode != 0:
        raise GmxUnavailable(f"gmx --version failed: {p.stderr[-300:]}")
    txt = p.stdout + p.stderr
    field = lambda k: next((ln.split(":", 1)[1].strip() for ln in txt.splitlines() if ln.startswith(k)), None)  # noqa: E731
    exe, prefix = field("Executable"), field("Data prefix")
    out = {"version": field("GROMACS version"), "precision": field("Precision"), "mpi": field("MPI library"),
           "launcher_sha256": B.R0.sha256_file(gmx), "executable_sha256": B.R0.sha256_file(exe) if exe and Path(exe).is_file() else None}
    pk = sorted(Path(prefix or ".").glob("conda-meta/gromacs-*.json")) if prefix else []
    if pk:
        d = json.loads(pk[0].read_text())
        out["package"] = {k: d.get(k) for k in ("name", "version", "build", "channel", "sha256", "url")}
    return json.loads(json.dumps(out).replace(prefix, "<env>")) if prefix else out


def parse_dump(lines) -> dict:
    """Header flags, natoms, box (nm), x (nm), v (nm/ps), selected inputrec values and molecule blocks from `gmx dump -s`."""
    out = {"flags": {}, "natoms": None, "box": None, "x": None, "v": None, "ir": {}, "ref_t": None, "tau_t": None,
           "nrdf": None, "molblocks": []}
    mode, buf, expect, rows = None, [], 0, None
    mb_cur = None
    for ln in lines:
        s = ln.rstrip("\n")
        if mode in ("box", "x", "v"):
            m = re.match(r"^\s+(box|x|v)\[\s*(\d+)\]=\{\s*([^}]*)\}", s)
            if m and m.group(1) == mode:
                if int(m.group(2)) != len(buf):
                    raise ValueError(f"{mode} row index {m.group(2)} out of order at {len(buf)}")
                buf.append([float(t) for t in m.group(3).replace(",", " ").split()])
                if len(buf) == expect:
                    out[mode] = np.array(buf, dtype=np.float64)
                    mode, buf = None, []
                continue
            raise ValueError(f"truncated {mode} block: {len(buf)} of {expect} rows")
        m = re.match(r"^(box|x|v) \((\d+)x(\d+)\):", s)
        if m:
            mode, expect, buf = m.group(1), int(m.group(2)), []
            if int(m.group(3)) != 3:
                raise ValueError(f"{m.group(1)} block is not N x 3")
            continue
        m = re.match(r"^\s+(bBox|bX|bV|bF)\s+=\s+(\w+(?: \w+)?)\s*$", s)
        if m:
            out["flags"][m.group(1)] = m.group(2) == "present"
            continue
        m = re.match(r"^\s+natoms\s+=\s+(\d+)\s*$", s)
        if m:
            out["natoms"] = int(m.group(1))
            continue
        m = re.match(r"^   ([A-Za-z][\w-]*)\s+=\s+(.*?)\s*$", s)
        if m and m.group(1) in IR_KEYS and m.group(1) not in out["ir"]:
            out["ir"][m.group(1)] = m.group(2)
            continue
        m = re.match(r"^   (ref-t|tau-t|nrdf):\s+(.*?)\s*$", s)
        if m:
            key = {"ref-t": "ref_t", "tau-t": "tau_t", "nrdf": "nrdf"}[m.group(1)]
            if out[key] is None:
                out[key] = m.group(2)
            continue
        m = re.match(r'^\s+moltype\s+=\s+(\d+)\s+"([^"]*)"', s)
        if m:
            mb_cur = {"moltype": int(m.group(1)), "name": m.group(2)}
            out["molblocks"].append(mb_cur)
            continue
        m = re.match(r"^\s+#molecules\s+=\s+(\d+)", s)
        if m and mb_cur is not None and "molecules" not in mb_cur:
            mb_cur["molecules"] = int(m.group(1))
    if mode is not None:
        raise ValueError(f"truncated {mode} block")
    if out["box"] is None or out["x"] is None or out["natoms"] is None:
        raise ValueError("dump has no box, x block or natoms")
    return out


def read_gmx(path, gmx: str, runner=subprocess.run) -> dict:
    p = runner([gmx, "dump", "-quiet", "-s", str(path)], capture_output=True, text=True, timeout=DUMP_TIMEOUT_S)
    if p.returncode != 0:
        raise ValueError(f"gmx dump exit code {p.returncode}: {p.stderr.strip()[-200:]}")
    d = parse_dump(p.stdout.splitlines())
    d["returncode"] = p.returncode
    return d


def mda_box(path):
    """The box read by MDAnalysis's own tpr routine (box 3x3 in nm), or None if unavailable."""
    try:
        from MDAnalysis.topology.tpr import utils as tu
        data = tu.TPXUnpacker(Path(path).read_bytes())
        th = tu.read_tpxheader(data)
        return np.asarray(tu.extract_box_info(data, th.fver)[0], dtype=np.float64)
    except Exception:
        return None


def box_valid(box) -> tuple:
    b = np.asarray(box, float)
    if b.shape != (3, 3) or not np.isfinite(b).all():
        return False, "box is not a finite 3x3"
    if abs(b[0, 1]) > 1e-9 or abs(b[0, 2]) > 1e-9 or abs(b[1, 2]) > 1e-9:
        return False, "box is not lower-triangular"
    if not (b[0, 0] > 0 and b[1, 1] > 0 and b[2, 2] > 0):
        return False, "box diagonal is not positive"
    return True, None


def box_dims_A(box_nm) -> np.ndarray:
    from MDAnalysis.lib.mdamath import triclinic_box
    b = np.asarray(box_nm, float) * 10.0
    return np.asarray(triclinic_box(b[0], b[1], b[2]), dtype=np.float64)


def agreement(rec: dict, g: dict, mbox) -> dict:
    """Validation (d): GROMACS vs the earlier reader, both in nm and nm/ps. MDAnalysis's TPRReader returns tpr coordinates in nm
    numerically (its units attribute says nm and nm/ps) although the frozen a1 loader labelled them Å; the comparison here is in
    nm, and a wrong unit assumption would show up as a disagreement on every run."""
    out = {"checks": {}, "problems": []}
    n = rec["n_atoms"]
    out["checks"]["natoms_equal"] = bool(g["natoms"] == n == len(g["x"]) and (g["v"] is None or len(g["v"]) == n))
    out["checks"]["header_box_and_x_present"] = bool(g["flags"].get("bBox") and g["flags"].get("bX"))
    ok_box, why = box_valid(g["box"])
    out["checks"]["gmx_box_valid"] = ok_box
    if not ok_box:
        out["problems"].append(why)
    if out["checks"]["natoms_equal"]:
        dx = float(np.max(np.abs(np.asarray(rec["positions"]) - g["x"])))
        out["max_abs_dx_nm"] = round(dx, 7)
        out["checks"]["positions_agree"] = dx <= POS_TOL_NM
        has_v = (g["v"] is not None, rec["velocities"] is not None)
        out["velocities"] = {"gmx": has_v[0], "earlier_reader": has_v[1]}
        if all(has_v):
            dv = float(np.max(np.abs(np.asarray(rec["velocities"]) - g["v"])))
            out["max_abs_dv_nm_ps"] = round(dv, 7)
            out["checks"]["velocities_agree"] = dv <= VEL_TOL_NM_PS
        else:
            out["checks"]["velocities_agree"] = not any(has_v)
    else:
        out["checks"]["positions_agree"] = out["checks"]["velocities_agree"] = False
    if mbox is not None:
        db = float(np.max(np.abs(np.asarray(mbox) - g["box"])))
        out["box_vs_mdanalysis_extract_box_info_nm"] = round(db, 9)
        out["checks"]["box_agrees_with_extract_box_info"] = db <= BOX_TOL_NM
    else:
        out["box_vs_mdanalysis_extract_box_info_nm"] = None
    out["ok"] = all(out["checks"].values())
    for k, v in out["checks"].items():
        if not v:
            out["problems"].append(k)
    return out


def unwrap_by_bonds(positions_A, dims, bonds, heavy_idx):
    """Heavy-atom coordinates made whole by walking the protein bond graph with minimum-image steps (box from GROMACS).

    Returns (coordinates of the heavy atoms in the largest connected component, in index order, number of components)."""
    from collections import deque

    from MDAnalysis.lib.distances import minimize_vectors
    pos = np.asarray(positions_A, dtype=np.float64)
    adj = collections.defaultdict(list)
    for i, j in np.asarray(bonds, dtype=np.int64):
        adj[int(i)].append(int(j))
        adj[int(j)].append(int(i))
    box = np.asarray(dims, dtype=np.float32)
    whole, comp, comps = {}, {}, 0
    nodes = sorted(adj)
    for root in nodes:
        if root in whole:
            continue
        whole[root], comp[root] = pos[root].copy(), comps
        q = deque([root])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if v in whole:
                    continue
                d = minimize_vectors((pos[v] - pos[u])[None].astype(np.float32), box)[0]
                whole[v], comp[v] = whole[u] + d, comps
                q.append(v)
        comps += 1
    sizes = collections.Counter(comp.values())
    main = max(sizes, key=sizes.get) if sizes else None
    keep = [] if main is None else [int(i) for i in heavy_idx if comp.get(int(i)) == main]
    return np.array([whole[i] for i in keep], dtype=np.float64).reshape(-1, 3), comps


def heavy_unwrapped(rec, positions_A, dims, heavy_idx):
    xyz, _ = unwrap_by_bonds(positions_A, dims, rec["protein_bonds"], heavy_idx)
    return xyz


SIDE: dict = {}


def load_tpr_gmx(path, gmx: str, mda_loader=None, gmx_reader=None, box_reader=None, heavy_fn=None) -> dict:
    """Record for the frozen validations with the GROMACS box, positions and velocities; raises on any disagreement."""
    mda_loader, gmx_reader = mda_loader or B.load_tpr, gmx_reader or (lambda p: read_gmx(p, gmx))
    box_reader, heavy_fn = box_reader or mda_box, heavy_fn or heavy_unwrapped
    rec = mda_loader(path)
    g = gmx_reader(path)
    ag = agreement(rec, g, box_reader(path))
    key = side_key(path)
    SIDE[key] = {"agreement": ag, "ir": g["ir"], "ref_t": g["ref_t"], "tau_t": g["tau_t"], "nrdf": g["nrdf"],
                 "molblocks": g["molblocks"], "flags": g["flags"], "box_nm": np.asarray(g["box"]).round(6).tolist(),
                 "box_volume_nm3": round(float(abs(np.linalg.det(np.asarray(g["box"])))), 6)}
    if not ag["ok"]:
        raise ValueError("GROMACS and the earlier reader disagree or the box is invalid: " + ", ".join(ag["problems"]))
    pos_A, dims = g["x"] * 10.0, box_dims_A(g["box"])
    names, resn = np.asarray(rec["names"]), np.asarray(rec["resnames"])
    prot = np.isin(resn, sorted(rec["protein_resnames"]))
    ca_idx = np.flatnonzero(np.isin(resn, sorted(B.R0.AMINO)) & (names == "CA"))
    heavy_idx = np.flatnonzero(prot & ~np.char.startswith(names.astype(str), "H"))
    out = dict(rec)
    out.update(positions=pos_A, dims=dims, velocities=None if g["v"] is None else g["v"] * 10.0, ca_positions=pos_A[ca_idx],
               ca_resnames=resn[ca_idx].tolist(), heavy_whole=heavy_fn(rec, pos_A, dims, heavy_idx),
               protein_vel=None if g["v"] is None else g["v"][prot] * 10.0)       # Å and Å/ps, from GROMACS (nm * 10)
    return out


def side_key(path) -> str:
    parts = Path(path).parts
    return f"{parts[-3]}:{'noMyr' if 'noMyr' in parts[-4] else 'Myr'}:{parts[-2]}"


def side_summary(side: dict) -> dict:
    runs = list(side.values())
    num = lambda k: [r["agreement"][k] for r in runs if r["agreement"].get(k) is not None]  # noqa: E731
    out = {"runs": len(runs), "agreement_ok": sum(r["agreement"]["ok"] for r in runs),
           "max_abs_dx_nm": {"max": max(num("max_abs_dx_nm"), default=None)},
           "max_abs_dv_nm_ps": {"max": max(num("max_abs_dv_nm_ps"), default=None)},
           "box_vs_extract_box_info_nm": {"max": max(num("box_vs_mdanalysis_extract_box_info_nm"), default=None)},
           "velocities_in_both": sum(bool(r["agreement"].get("velocities", {}).get("gmx")) and
                                     bool(r["agreement"].get("velocities", {}).get("earlier_reader")) for r in runs)}
    desc = {}
    for s in B.SETS:
        for c in B.CONDITIONS:
            rs = [v for k, v in side.items() if k.startswith(f"{s}:{c}:")]
            if not rs:
                continue
            distinct = lambda f: sorted({json.dumps(f(r), sort_keys=True) for r in rs})  # noqa: E731
            desc[f"{s}:{c}"] = {"runs": len(rs), "ir": {k: distinct(lambda r, k=k: r["ir"].get(k)) for k in IR_KEYS},
                                "ref_t": distinct(lambda r: r["ref_t"]), "tau_t": distinct(lambda r: r["tau_t"]),
                                "nrdf": distinct(lambda r: r["nrdf"]), "molblocks": distinct(lambda r: r["molblocks"]),
                                "box_volume_nm3": {"min": min(r["box_volume_nm3"] for r in rs), "max": max(r["box_volume_nm3"] for r in rs)},
                                "box_nm_distinct": len(distinct(lambda r: r["box_nm"]))}
    out["run_input_settings"] = desc
    diffs = {}
    for s in B.SETS:
        a, b = desc.get(f"{s}:Myr"), desc.get(f"{s}:noMyr")
        if a and b:
            diffs[s] = {k: {"Myr": a["ir"][k], "noMyr": b["ir"][k]} for k in IR_KEYS if a["ir"][k] != b["ir"][k]}
            for k in ("ref_t", "tau_t", "nrdf", "molblocks"):
                if a[k] != b[k]:
                    diffs[s][k] = {"Myr": a[k], "noMyr": b[k]}
    out["myr_vs_nomyr_run_input_differences_descriptive"] = diffs
    return out


# --- run ----------------------------------------------------------------------------------------------------------------

def run_check(orig_spec, data_dir, extract_dir, spec, gmx, loader=None, log=print):
    SIDE.clear()
    loader = loader or (lambda p: load_tpr_gmx(p, gmx))
    with A1.amended():
        res, ca = B.run_check(orig_spec, data_dir, extract_dir, spec, loader=loader, log=log)
    res["a2"] = side_summary(SIDE) if SIDE else {"runs": 0}
    res["a2"]["per_run"] = {k: {"agreement": v["agreement"], "box_nm": v["box_nm"], "box_volume_nm3": v["box_volume_nm3"]}
                            for k, v in SIDE.items()}
    return res, ca


def _h(p):
    return B._h(p)


def cmd_freeze(args) -> int:
    if FREEZE.exists():
        print("freeze exists")
        return 1
    if (B.ORIG / "data_dryad/extracted_1b_a2").exists():
        raise SystemExit("STOP: extracted_1b_a2 exists; the a2 freeze must precede its run")
    A1.verify()
    gmx = gmx_executable()
    if not gmx:
        raise SystemExit("STOP: TOPOFOLD_GMX is not set to a gmx executable")
    rec = {"id": "ABL-MYR-1b-CHECK-a2", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "amendment_sha256": B.R0.sha256_file(AMENDMENT), "code": {"path": str(SELF.relative_to(B.ROOT)), "sha256": B.R0.sha256_file(SELF)},
           "test": {TEST: _h(TEST)}, "acquisition_record_sha256": _h(ACQUISITION),
           "gromacs_at_freeze": gmx_provenance(gmx),
           "frozen_check": {"spec_sha256": B.R0.sha256_file(B.SPEC), "code_sha256": B.R0.sha256_file(B.SELF)},
           "a1_records_sha256": {p: _h(p) for p in A1_RECORDS}, "extracted_before_freeze": []}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps({k: rec[k] for k in ("frozen_utc", "amendment_sha256", "code")}, indent=1))
    return 0


def verify() -> dict:
    A1.verify()
    f = json.loads(FREEZE.read_text())
    bad = [n for n, ok in (("amendment", B.R0.sha256_file(AMENDMENT) == f["amendment_sha256"]),
                           ("a2 code", B.R0.sha256_file(SELF) == f["code"]["sha256"]),
                           ("acquisition record", _h(ACQUISITION) == f["acquisition_record_sha256"])) if not ok]
    bad += [p for p, h in f["a1_records_sha256"].items() if _h(p) != h]
    if bad:
        raise SystemExit(f"STOP: files differ from the a2 freeze: {bad}")
    return f


def cmd_run(args) -> int:
    f = verify()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    gmx = gmx_executable()
    if not gmx:
        res = {"id": "ABL-MYR-1b-CHECK-a2", "outcome": "check_incomplete",
               "stop": {"code": "S5", "cause": "TOPOFOLD_GMX is not set to a gmx executable and none is on PATH"}}
    else:
        try:
            prov = gmx_provenance(gmx)
        except (GmxUnavailable, subprocess.SubprocessError, OSError) as e:
            res = {"id": "ABL-MYR-1b-CHECK-a2", "outcome": "check_incomplete", "stop": {"code": "S5", "cause": f"GROMACS cannot be run: {e}"}}
        else:
            spec = json.loads(B.SPEC.read_text())
            res, ca = run_check(json.loads((B.ORIG / "run_spec.json").read_text()), B.ORIG / "data_dryad",
                                B.ORIG / "data_dryad/extracted_1b_a2", spec, gmx)
            res["a2"]["gromacs"] = prov
            res["a2"]["gromacs_matches_freeze"] = prov == f["gromacs_at_freeze"]
            np.savez(out / CA_NPZ, **{k.replace(":", "_"): v for k, v in ca.items() if v is not None})
            res["production_start_ca_npz"] = {"file": CA_NPZ, "sha256": B.R0.sha256_file(out / CA_NPZ), "units": "Å",
                                              "atoms": "protein Cα, whole"}
    res["amendment"] = {"id": "ABL-MYR-1b-CHECK-a2", "amendment_sha256": f["amendment_sha256"], "code_sha256": f["code"]["sha256"],
                        "frozen_check_spec_sha256": f["frozen_check"]["spec_sha256"],
                        "frozen_check_code_sha256": f["frozen_check"]["code_sha256"]}
    (out / RESULTS).write_text(json.dumps(B.R0.scrub(res), indent=1, default=B.R0._json_default) + "\n")
    print(res["outcome"], json.dumps(res.get("resources")))
    return 0


def render_report(res: dict) -> str:
    head = ("> Amendment a2 of the ABL-MYR-1b check (check_amendment_a2.json): the stored box, positions and velocities of the 40 "
            "production md.tpr are read independently with GROMACS and compared with the earlier reader. One bounded attempt; "
            "criteria K1-K5 and their thresholds are those of the frozen check. 'not_comparable' means criteria not met, not "
            "that the systems are physically incomparable; 'comparable_for_1b' means criteria passed, not that no confound exists.\n\n")
    body = B.render_report(res)
    a2 = res.get("a2")
    if a2 and a2.get("runs"):
        g = a2.get("gromacs", {})
        extra = ["", "## a2: independent reader", "",
                 f"- GROMACS {g.get('version')} ({g.get('precision')}); matches the freeze: {a2.get('gromacs_matches_freeze')}",
                 f"- reader agreement (validation d): {a2['agreement_ok']} of {a2['runs']} run inputs; max |dx| "
                 f"{a2['max_abs_dx_nm']['max']} nm (tolerance {POS_TOL_NM}), max |dv| {a2['max_abs_dv_nm_ps']['max']} nm/ps "
                 f"(tolerance {VEL_TOL_NM_PS}); velocities in both readers in {a2['velocities_in_both']} runs; GROMACS box vs "
                 f"MDAnalysis extract_box_info max {a2['box_vs_extract_box_info_nm']['max']} nm (tolerance {BOX_TOL_NM})",
                 "- descriptive Myr vs noMyr differences read from the run inputs: " +
                 json.dumps(a2["myr_vs_nomyr_run_input_differences_descriptive"], sort_keys=True), ""]
        body += "\n".join(extra) + "\n"
    return head + body


def cmd_report(args) -> int:
    res = json.loads((Path(args.out) / RESULTS).read_text())
    (Path(args.out) / "ABL-MYR-1b_check_a2_report.md").write_text(render_report(res))
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
