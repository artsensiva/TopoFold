"""DEV-CA-1: run the Cα pipeline once on the local AdK open-start trajectory under run_spec.json.

Run from the repository root with the core interpreter:
    .venv/bin/python benchmarks/audit/dev_ca_pipeline/run_dev.py
Writes bulky arrays to build/ (git-ignored) and summaries to results/.
"""

import hashlib
import json
import os
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
ADAPTER_PY = ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    spec_path = HERE / "run_spec.json"
    spec = json.loads(spec_path.read_text())
    run_id = spec["id"]
    build = HERE / "build"
    build.mkdir(exist_ok=True)

    hashes = {}
    for rel, expected in spec["dataset"]["files"].items():
        got = sha256(ROOT / rel)
        hashes[rel] = {"expected": expected, "observed": got, "match": got == expected}
    if not all(h["match"] for h in hashes.values()):
        print("STOP: source hashes differ from the G0 record", json.dumps(hashes, indent=1))
        return 2

    fs, pol = spec["frame_subset"], spec["policies"]
    gro, xtc = list(spec["dataset"]["files"])
    cmd = ["/usr/bin/time", "-v", sys.executable, "-m", "topofold_inputs.pipeline", str(ROOT / gro), str(ROOT / xtc),
           "--run-id", run_id, "--out", str(build), "--adapter-python", str(ADAPTER_PY),
           "--selection", spec["selection"], "--chain-attribute", pol["chain_attribute"], "--altloc", pol["altloc"],
           "--pbc", pol["pbc"], "--max-ca-ca", str(pol["max_ca_ca_angstrom"]),
           "--min-segment", str(pol["min_segment_angstrom"]), "--collinear-sin-tol", str(pol["collinear_sin_tol"]),
           "--start", str(fs["start"]), "--step", str(fs["step"])]
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"))
    cap = spec["resource_cap"]["cpu_seconds"]
    t0 = time.time()
    proc = subprocess.run(cmd, env=env, capture_output=True, text=True,
                          preexec_fn=lambda: resource.setrlimit(resource.RLIMIT_CPU, (cap, cap)))
    wall = time.time() - t0
    (build / "pipeline_stdout.txt").write_text(proc.stdout)
    (build / "pipeline_stderr.txt").write_text(proc.stderr)

    def tv(label):
        m = re.search(re.escape(label) + r":\s*([0-9.:]+)", proc.stderr)
        return m.group(1) if m else None

    resources = {"wall_seconds": round(wall, 2), "user_cpu_seconds": tv("User time (seconds)"),
                 "system_cpu_seconds": tv("System time (seconds)"),
                 "max_rss_kib_largest_process": tv("Maximum resident set size (kbytes)"),
                 "cpu_cap_seconds": cap, "exit_status": proc.returncode,
                 "measurement": "GNU time -v over the pipeline process tree (adapter runs as a child process)"}
    summary = {"run_id": run_id, "spec_sha256": sha256(spec_path), "source_hashes": hashes, "resources": resources}
    if proc.returncode != 0:
        summary["status"] = "failed"
        (HERE / "results" / f"{run_id}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(proc.stderr[-2000:])
        return proc.returncode

    a = json.loads((build / f"{run_id}.json").read_text())
    f = json.loads((build / f"{run_id}_features.json").read_text())
    exp = spec["dataset"]["expected_from_g0"]
    res = a["residues"]
    frames = a["frames"]["frame_index"]
    times = a["frames"]["time_ps"]
    checks = {
        "source_hashes_match_g0": all(h["match"] for h in hashes.values()),
        "residue_count": len(res) == exp["n_ca"],
        "resids_1_to_214_consecutive": [r["resid"] for r in res] == list(range(1, exp["n_ca"] + 1)),
        "all_present": all(r["status"] == "present" for r in res),
        "frames_expected": frames == list(range(0, exp["n_frames"], fs["step"]))
                           and len(frames) == fs["expected_frames"],
        "times_equal_dt_times_index": all(abs(t - exp["nominal_dt_ps"] * i) < 1e-3 for t, i in zip(times, frames)),
        "units_converted_from_nm": a["units"]["source_length_unit"] == "nm" and a["units"]["coordinates"] == "angstrom",
        "frame_identity_not_verifiable": a["frame_identity"]["status"] == "not_verifiable",
        "no_nonfinite_valid_features": f["non_finite_valid_values"] == 0,
        "no_breaks": len(a["breaks"]) == 0,
        "single_segment": a["segments_half_open"] == [[0, exp["n_ca"]]],
        "kappa_valid_full": f["counts"]["kappa"]["valid"] == fs["expected_frames"] * (exp["n_ca"] - 2),
        "tau_valid_full": f["counts"]["tau"]["valid"] == fs["expected_frames"] * (exp["n_ca"] - 3),
        "frames_within_cap": len(frames) <= spec["resource_cap"]["max_frames"],
    }
    # independent NumPy reference on the exported coordinates (textbook turning angle and IUPAC dihedral)
    import numpy as np
    c = np.load(build / f"{run_id}.npz")["coordinates_angstrom"].astype(np.float32).astype(np.float64)
    fz = np.load(build / f"{run_id}_features.npz")
    u, v = c[:, 1:-1] - c[:, :-2], c[:, 2:] - c[:, 1:-1]
    k_ref = np.arccos(np.clip((u * v).sum(-1) / np.linalg.norm(u, axis=-1) / np.linalg.norm(v, axis=-1), -1, 1))
    b1, b2, b3 = c[:, 1:-2] - c[:, :-3], c[:, 2:-1] - c[:, 1:-2], c[:, 3:] - c[:, 2:-1]
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    t_ref = np.arctan2(np.linalg.norm(b2, axis=-1) * (b1 * n2).sum(-1), (n1 * n2).sum(-1))
    kv, tv = fz["kappa_valid"][:, 1:-1], fz["tau_valid"][:, 1:-1]
    dk = np.abs(fz["kappa_rad"][:, 1:-1] - k_ref)[kv]
    dt = np.abs(np.angle(np.exp(1j * (fz["tau_rad"][:, 1:-1] - t_ref))))[tv]
    dist = np.linalg.norm(np.diff(c, axis=1), axis=-1)
    reference = {"max_abs_kappa_diff_rad": float(dk.max()), "max_circular_tau_diff_rad": float(dt.max()),
                 "cells_compared": {"kappa": int(kv.sum()), "tau": int(tv.sum())},
                 "ca_ca_distance_angstrom": {"min": float(dist.min()), "max": float(dist.max())}}
    checks["independent_reference_agrees_1e-9"] = reference["max_abs_kappa_diff_rad"] < 1e-9 and \
        reference["max_circular_tau_diff_rad"] < 1e-9
    summary.update({
        "status": "completed",
        "independent_reference": reference,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "mapping": {"residues": len(res), "segments": a["segments_half_open"], "breaks": a["breaks"],
                    "chain_attribute": a["config"]["chain_attribute"], "chains": sorted({r["chain"] for r in res}),
                    "frame_identity": a["frame_identity"], "units": a["units"], "pbc": a["pbc"],
                    "precision": a["precision"], "read_warnings": a["provenance"]["read_warnings"]},
        "features": {"counts": f["counts"], "non_finite_valid_values": f["non_finite_valid_values"],
                     "diagnostics": f["diagnostics"], "core": f["provenance"]["core"]},
        "outputs_sha256": {p.name: sha256(p) for p in sorted(build.glob(f"{run_id}*"))},
        "evidence_level": "development evidence (engineering); no biological, novelty, superiority or state claim",
    })
    (HERE / "results" / f"{run_id}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for name in (f"{run_id}_report.md", f"{run_id}_features.json", f"{run_id}.json"):
        (HERE / "results" / name).write_bytes((build / name).read_bytes())
    print(json.dumps({"checks": checks, "resources": resources}, indent=1))
    return 0 if summary["all_checks_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
