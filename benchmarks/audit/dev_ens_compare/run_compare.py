"""DEV-ENS-1: descriptive per-window comparison of the AdK open-start and closed-start trajectories.

    .venv/bin/python benchmarks/audit/dev_ens_compare/run_compare.py freeze   # record spec + code hashes
    .venv/bin/python benchmarks/audit/dev_ens_compare/run_compare.py run      # execute under the freeze

Bulky arrays go to build/ (git-ignored); machine-readable results and the generated report go to results/.
Development evidence only: no test, no state assignment, no biological interpretation.
"""

import datetime as dt
import hashlib
import json
import os
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs.compare import (  # noqa: E402
    circular_hist, compare_channel, equal_duration_blocks, equalize, linear_hist, spearman, window_distances,
)

ADAPTER_PY = ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python"
SPEC = HERE / "run_spec.json"
FREEZE = HERE / "results" / "DEV-ENS-1_freeze.json"
CODE = ["python/topofold_inputs/compare.py", "python/topofold_inputs/features.py",
        "python/topofold_inputs/ca_adapter.py", "python/topofold_inputs/pipeline.py",
        "python/topofold_inputs/__main__.py", "benchmarks/audit/dev_ens_compare/run_compare.py",
        "tests/test_ensemble_compare.py"]


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def freeze():
    if FREEZE.exists():
        print(f"freeze already exists: {FREEZE}")
        return 1
    rec = {"id": "DEV-ENS-1", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "spec_sha256": sha256(SPEC), "code_sha256": {c: sha256(ROOT / c) for c in CODE}}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=1))
    return 0


def run_pipeline(label, run, spec, build, cpu_left):
    obs, pol = spec["observation_period"]["source_frame_slice"], spec["pipeline_policies"]
    run_id = f"DEV-ENS-1_{label}"
    cmd = ["/usr/bin/time", "-v", sys.executable, "-m", "topofold_inputs.pipeline", str(ROOT / run["topology"]),
           str(ROOT / run["trajectory"]), "--run-id", run_id, "--out", str(build), "--adapter-python", str(ADAPTER_PY),
           "--selection", pol["selection"], "--chain-attribute", pol["chain_attribute"], "--altloc", pol["altloc"],
           "--pbc", pol["pbc"], "--max-ca-ca", str(pol["max_ca_ca_angstrom"]),
           "--min-segment", str(pol["min_segment_angstrom"]), "--collinear-sin-tol", str(pol["collinear_sin_tol"]),
           "--start", str(obs["start"]), "--stop", str(obs["stop"]), "--step", str(obs["step"])]
    cap = max(int(cpu_left), 1)
    proc = subprocess.run(cmd, env=dict(os.environ, PYTHONPATH=str(ROOT / "python")), capture_output=True, text=True,
                          preexec_fn=lambda: resource.setrlimit(resource.RLIMIT_CPU, (cap, cap)))
    (build / f"{run_id}_stdout.txt").write_text(proc.stdout)
    (build / f"{run_id}_stderr.txt").write_text(proc.stderr)

    def tv(label_):
        m = re.search(re.escape(label_) + r":\s*([0-9.:]+)", proc.stderr)
        return float(m.group(1)) if m and ":" not in m.group(1) else (m.group(1) if m else None)

    res = {"exit_status": proc.returncode, "user_cpu_s": tv("User time (seconds)"),
           "system_cpu_s": tv("System time (seconds)"), "max_rss_kib": tv("Maximum resident set size (kbytes)"),
           "cpu_cap_s": cap}
    if proc.returncode != 0:
        raise SystemExit(f"pipeline failed for {label}: {proc.stderr[-1500:]}")
    return run_id, res


def load(build, run_id):
    rep = json.loads((build / f"{run_id}.json").read_text())
    a = np.load(build / f"{run_id}.npz")
    f = np.load(build / f"{run_id}_features.npz")
    return rep, a, f


def make_hist(channel, cfg):
    bins = cfg["bins"]
    if channel["type"] == "circular":
        origin = channel["origin"] + cfg["tau_origin_shift_bins"] * 2 * np.pi / bins
        return lambda x: circular_hist(x, bins, origin)
    return lambda x: linear_hist(x, channel["lo"], channel["hi"], bins)[0]


def summarize(rows):
    use = [r for r in rows if not r["excluded"]]
    fb = np.array([r["full_between"] for r in use])
    frac = np.array([r["fraction_block_between_above_max_within"] for r in use])
    return {"windows_compared": len(use), "windows_excluded": len(rows) - len(use),
            "full_between": {"min": float(fb.min()), "median": float(np.median(fb)), "max": float(fb.max())},
            "windows_all_block_between_above_max_within": int((frac == 1.0).sum()),
            "windows_no_block_between_above_max_within": int((frac == 0.0).sum())}


def main_run():
    spec = json.loads(SPEC.read_text())
    frz = json.loads(FREEZE.read_text())
    mismatch = {c: (h, sha256(ROOT / c)) for c, h in frz["code_sha256"].items() if sha256(ROOT / c) != h}
    if frz["spec_sha256"] != sha256(SPEC) or mismatch:
        raise SystemExit(f"STOP: spec or code differs from the freeze: {mismatch}")
    build = HERE / "build"
    build.mkdir(exist_ok=True)
    t_wall = time.time()

    src = {}
    for label, run in spec["dataset"]["runs"].items():
        for k in ("topology", "trajectory"):
            got = sha256(ROOT / run[k])
            src[run[k]] = {"expected": run[f"{k}_sha256"], "observed": got, "match": got == run[f"{k}_sha256"]}
    if not all(v["match"] for v in src.values()):
        raise SystemExit(f"STOP: source hashes differ: {src}")

    cap = spec["resource_cap"]["cpu_seconds"]
    used, resources, loaded = 0.0, {}, {}
    for label, run in spec["dataset"]["runs"].items():
        run_id, res = run_pipeline(label, run, spec, build, cap - used)
        used += (res["user_cpu_s"] or 0) + (res["system_cpu_s"] or 0)
        resources[label] = res
        loaded[label] = load(build, run_id)

    # --- observation-period and mapping checks ---------------------------------------------------------------------
    (ra, aa, fa), (rb, ab, fb) = loaded["open"], loaded["closed"]
    ta, tb = aa["time_ps"], ab["time_ps"]
    obs = spec["observation_period"]
    ident = lambda r: [(x["chain"], x["resid"], x["icode"], x["resname"]) for x in r["residues"]]  # noqa: E731
    checks = {
        "source_hashes_match": True,
        "samples_per_run": int(len(ta)) == obs["expected_samples_per_run"] == int(len(tb)),
        "total_within_cap": int(len(ta) + len(tb)) <= spec["resource_cap"]["max_total_frames"],
        "timestamps_identical": bool(np.array_equal(ta, tb)),
        "timestamps_800ps_grid": bool(np.allclose(ta, 800.0 * np.arange(len(ta)))),
        "within_common_interval": bool(ta.min() >= 0 and ta.max() <= 335400.0),
        "residue_identity_identical": ident(ra) == ident(rb),
        "no_breaks_either_run": not ra["breaks"] and not rb["breaks"],
        "native_extension_resolved": all(json.loads((build / f"DEV-ENS-1_{k}_features.json").read_text())
                                         ["provenance"]["core"]["native_extension"]["status"] == "resolved"
                                         for k in ("open", "closed")),
    }
    if not all(checks.values()):
        raise SystemExit(f"STOP: pre-comparison checks failed: {checks}")

    # --- blocks ----------------------------------------------------------------------------------------------------
    nb = spec["blocks"]["n_blocks"]
    raw_a = equal_duration_blocks(ta, 0.0, 335400.0, nb)
    raw_b = equal_duration_blocks(tb, 0.0, 335400.0, nb)
    eq, n_block = equalize(raw_a + raw_b)
    blocks_a, blocks_b = eq[:nb], eq[nb:]
    full = np.arange(len(ta))
    block_record = {"raw_sizes": {"open": [len(b) for b in raw_a], "closed": [len(b) for b in raw_b]},
                    "equalized_size": n_block,
                    "time_ranges_ps": [[float(ta[b[0]]), float(ta[b[-1]])] for b in blocks_a],
                    "dropped_samples": int(sum(len(b) for b in raw_a + raw_b) - 2 * nb * n_block)}

    # --- channels ----------------------------------------------------------------------------------------------------
    Xa, Xb = aa["coordinates_angstrom"], ab["coordinates_angstrom"]
    data = {
        "kappa": (fa["kappa_rad"], fb["kappa_rad"], fa["kappa_valid"], fb["kappa_valid"]),
        "d_i_i2": (window_distances(Xa, 2), window_distances(Xb, 2), None, None),
        "tau": (fa["tau_rad"], fb["tau_rad"], fa["tau_valid"], fb["tau_valid"]),
        "d_i_i3": (window_distances(Xa, 3), window_distances(Xb, 3), None, None),
    }
    # matched validity per pair: feature valid and distance finite, in both channels of the pair
    valid = {}
    for p, q in spec["pairs"]:
        va = data[p][2] & np.isfinite(data[q][0])
        vb = data[p][3] & np.isfinite(data[q][1])
        valid[p] = valid[q] = (va, vb)

    def window_label(ch, w):
        res = ra["residues"]
        span = (w - 1, w + 1) if ch in ("kappa", "d_i_i2") else (w - 1, w + 2)
        return f"{res[span[0]]['resid']}-{res[span[1]]['resid']}" if span[0] >= 0 and span[1] < len(res) else None

    configs = [dict(name="primary", **spec["configurations"]["primary"])] + spec["configurations"]["sensitivity"]
    results = {}
    for cfg in configs:
        per = {}
        for ch, (xa, xb, _, _) in data.items():
            hist = make_hist(spec["channels"][ch], cfg)
            va, vb = valid[ch]
            rows = compare_channel(np.nan_to_num(xa), np.nan_to_num(xb), va, vb, blocks_a, blocks_b, hist)
            for r in rows:
                r["residues"] = window_label(ch, r["window"])
            oor = None
            if spec["channels"][ch]["type"] == "linear":
                c = spec["channels"][ch]
                both = np.concatenate([xa[va], xb[vb]])
                oor = {"underflow": int((both < c["lo"]).sum()), "overflow": int((both > c["hi"]).sum()),
                       "observations": int(both.size)}
            per[ch] = {"rows": rows, "summary": summarize(rows), "out_of_range": oor}
        results[cfg["name"]] = per

    # --- descriptive rank stability and sensitivity -------------------------------------------------------------------
    prim = results["primary"]

    def vec(per, ch, key="full_between", k=None):
        rows = [r for r in per[ch]["rows"] if not r["excluded"]]
        return np.array([r[key][k] if k is not None else r[key] for r in rows]), [r["window"] for r in rows]

    stability = {}
    for ch in data:
        m = [vec(prim, ch, "matched_block_between", k)[0] for k in range(nb)]
        stability[ch] = {"spearman_matched_block_pairs": {f"{i}-{j}": spearman(m[i], m[j])
                                                           for i in range(nb) for j in range(i + 1, nb)}}
    paired = {f"{p}~{q}": spearman(vec(prim, p)[0], vec(prim, q)[0]) for p, q in spec["pairs"]}
    top = {}
    for ch in data:
        v, w = vec(prim, ch)
        order = np.argsort(-v, kind="mergesort")[:10]
        top[ch] = [{"window": int(w[i]), "residues": window_label(ch, int(w[i])), "full_between": float(v[i]),
                    "fraction_block_between_above_max_within":
                        next(r for r in prim[ch]["rows"] if r["window"] == w[i])
                        ["fraction_block_between_above_max_within"]} for i in order]
    sensitivity = {}
    for cfg in configs[1:]:
        sensitivity[cfg["name"]] = {}
        for ch in data:
            v0, w0 = vec(prim, ch)
            v1, w1 = vec(results[cfg["name"]], ch)
            assert w0 == w1
            t0 = set(np.array(w0)[np.argsort(-v0, kind="mergesort")[:10]])
            t1 = set(np.array(w1)[np.argsort(-v1, kind="mergesort")[:10]])
            sensitivity[cfg["name"]][ch] = {"spearman_vs_primary": spearman(v0, v1), "top10_overlap": len(t0 & t1),
                                            "summary": results[cfg["name"]][ch]["summary"]}

    self_ru = resource.getrusage(resource.RUSAGE_SELF)
    resources["comparison_stage"] = {"user_cpu_s": self_ru.ru_utime, "system_cpu_s": self_ru.ru_stime,
                                     "max_rss_kib": self_ru.ru_maxrss}
    total_cpu = sum((r.get("user_cpu_s") or 0) + (r.get("system_cpu_s") or 0) for r in resources.values())
    resources["total_cpu_s"] = total_cpu
    resources["wall_s"] = time.time() - t_wall
    resources["within_cpu_cap"] = total_cpu <= cap

    out = {"id": "DEV-ENS-1", "evidence_level": "descriptive development comparison; no test, state or biological claim",
           "spec_sha256": frz["spec_sha256"], "freeze_utc": frz["frozen_utc"], "source_hashes": src,
           "checks": checks, "observation": {"samples_per_run": int(len(ta)), "time_first_ps": float(ta[0]),
                                             "time_last_ps": float(ta[-1]), "sample_interval_ps": 800.0},
           "blocks": block_record, "summaries": {c: {ch: results[c][ch]["summary"] for ch in data} for c in results},
           "out_of_range_primary": {ch: prim[ch]["out_of_range"] for ch in data},
           "rank_stability_primary": stability, "paired_channel_spearman_primary": paired,
           "top10_primary": top, "sensitivity": sensitivity, "resources": resources,
           "per_window_primary": {ch: prim[ch]["rows"] for ch in data},
           "outputs_sha256": {p.name: sha256(p) for p in sorted(build.glob("DEV-ENS-1_*"))}}
    res_path = HERE / "results" / "DEV-ENS-1_results.json"
    res_path.write_text(json.dumps(out, indent=1) + "\n")
    (HERE / "results" / "DEV-ENS-1_report.md").write_text(report(out, spec, sha256(res_path)))
    print(json.dumps({"checks": checks, "blocks": block_record, "resources": resources}, indent=1))
    return 0


def report(o, spec, res_hash):
    L = [f"# DEV-ENS-1: descriptive per-window comparison (AdK open-start vs closed-start)", "",
         "Generated from `DEV-ENS-1_results.json` (SHA-256 " + res_hash[:16] + "…). Development evidence only: "
         "one run per starting condition, so differences mix the starting structure with relaxation from it. "
         "No test, state assignment, population or biological interpretation. Histogram Jensen–Shannon distance "
         "(base 2) inspired by established per-feature comparison; not a reproduction of PENSA.", "",
         "## Observation period and blocks", "",
         f"- Samples: {o['observation']['samples_per_run']} per run, {o['observation']['time_first_ps']:.0f}–"
         f"{o['observation']['time_last_ps']:.0f} ps every {o['observation']['sample_interval_ps']:.0f} ps; "
         f"timestamps identical in both runs: {o['checks']['timestamps_identical']}",
         f"- Blocks: raw sizes {o['blocks']['raw_sizes']}; equalized to {o['blocks']['equalized_size']} samples "
         f"(dropped {o['blocks']['dropped_samples']}); time ranges (ps) {o['blocks']['time_ranges_ps']}", "",
         "## Primary configuration (36 bins)", "",
         "| Channel | Windows | Excluded | Full-run JS distance min / median / max | Windows with all 9 between-run "
         "block distances above the largest within-run one | … with none |", "|---|---|---|---|---|---|"]
    for ch, s in o["summaries"]["primary"].items():
        fb = s["full_between"]
        L.append(f"| {ch} | {s['windows_compared']} | {s['windows_excluded']} | {fb['min']:.3f} / {fb['median']:.3f}"
                 f" / {fb['max']:.3f} | {s['windows_all_block_between_above_max_within']} | "
                 f"{s['windows_no_block_between_above_max_within']} |")
    L += ["", "Block fractions are descriptive (blocks are not independent replicates). Full-run and block distances "
          "use different sample sizes and are not compared in magnitude.", "",
          "Out-of-range observations (primary): " + ", ".join(
              f"{ch} {v['underflow']} under / {v['overflow']} over of {v['observations']}"
              for ch, v in o["out_of_range_primary"].items() if v), "",
          "## Rank stability across time (Spearman between matched block-pair rankings, primary)", ""]
    for ch, s in o["rank_stability_primary"].items():
        L.append(f"- {ch}: " + ", ".join(f"blocks {k}: {v:.3f}" for k, v in s["spearman_matched_block_pairs"].items()))
    L += ["", "Paired channels on the same windows (Spearman of full-run distances): " + ", ".join(
        f"{k} {v:.3f}" for k, v in o["paired_channel_spearman_primary"].items()), "",
          "κ and the i→i+2 distance are near-deterministic functions of each other for near-constant bond lengths; the "
          "i→i+3 distance cannot distinguish τ from −τ. Larger distances or different rankings do not show that one "
          "representation is better.", "", "## Largest full-run distances (primary, descriptive)", ""]
    for ch, rows in o["top10_primary"].items():
        L.append(f"- {ch}: " + ", ".join(f"{r['residues']} ({r['full_between']:.3f}; "
                                          f"{r['fraction_block_between_above_max_within']:.2f})" for r in rows))
    L += ["", "Each entry: window residues (full-run JS distance; fraction of between-run block distances above the "
          "largest within-run block distance).", "", "## Sensitivity (predeclared; not used for selection)", "",
          "| Configuration | Channel | Spearman vs primary | Top-10 overlap | Full-run median |", "|---|---|---|---|---|"]
    for cfg, chs in o["sensitivity"].items():
        for ch, s in chs.items():
            L.append(f"| {cfg} | {ch} | {s['spearman_vs_primary']:.3f} | {s['top10_overlap']}/10 | "
                     f"{s['summary']['full_between']['median']:.3f} |")
    r = o["resources"]
    L += ["", "## Resources", "",
          f"- Total CPU {r['total_cpu_s']:.1f} s (cap {spec['resource_cap']['cpu_seconds']} s); wall {r['wall_s']:.1f} s; "
          f"max RSS per stage (KiB): open {r['open']['max_rss_kib']}, closed {r['closed']['max_rss_kib']}, "
          f"comparison {r['comparison_stage']['max_rss_kib']}", ""]
    return "\n".join(L)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("freeze", "run"):
        print(__doc__)
        sys.exit(2)
    sys.exit(freeze() if sys.argv[1] == "freeze" else main_run())
