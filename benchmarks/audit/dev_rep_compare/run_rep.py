"""DEV-REP-1: replicated descriptive comparison of mdCATH 1ha8A00 at 320 K vs 348 K (5 replicas each).

    .venv/bin/python benchmarks/audit/dev_rep_compare/run_rep.py freeze   # record spec + code hashes
    .venv/bin/python benchmarks/audit/dev_rep_compare/run_rep.py run      # execute under the freeze

Stage 1 (HDF5 bridge + adapter) runs in the WP-D audit venv; features and comparison run here. Arrays go to build/
(git-ignored); results/ holds machine-readable results and the generated report. Development evidence only.
"""

import collections
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
    circular_hist, linear_hist, replica_comparison, spearman, window_distances,
)
from topofold_inputs.features import export_features, features_from_export  # noqa: E402

ADAPTER_PY = ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python"
SPEC = HERE / "run_spec.json"
FREEZE = HERE / "results" / "DEV-REP-1_freeze.json"
CODE = ["python/topofold_inputs/compare.py", "python/topofold_inputs/features.py",
        "python/topofold_inputs/ca_adapter.py", "python/topofold_inputs/mdcath_bridge.py",
        "benchmarks/audit/dev_rep_compare/run_rep.py", "tests/test_ensemble_compare.py",
        "tests/test_mdcath_bridge.py"]
CATS = ("consistent", "replica_dependent", "inconclusive")


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
    rec = {"id": "DEV-REP-1", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "spec_sha256": sha256(SPEC), "code_sha256": {c: sha256(ROOT / c) for c in CODE}}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps(rec, indent=1))
    return 0


def bridge(spec, T, r, build, cpu_left):
    s = spec["sampling"]["used_frames"]
    run_id = f"DEV-REP-1_{T}_{r}"
    cmd = ["/usr/bin/time", "-v", str(ADAPTER_PY), "-m", "topofold_inputs.mdcath_bridge", str(ROOT / spec["dataset"]["file"]),
           spec["dataset"]["domain"], T, r, "--run-id", run_id, "--out", str(build), "--start", str(s["start"]),
           "--stop", str(s["stop"]), "--step", str(s["step"]), "--pbc", "make_whole_sequential",
           "--altloc", spec["pipeline_policies"]["altloc"], "--max-ca-ca", str(spec["pipeline_policies"]["max_ca_ca_angstrom"])]
    cap = max(int(cpu_left), 1)
    p = subprocess.run(cmd, env=dict(os.environ, PYTHONPATH=str(ROOT / "python")), capture_output=True, text=True,
                       preexec_fn=lambda: resource.setrlimit(resource.RLIMIT_CPU, (cap, cap)))
    (build / f"{run_id}_stderr.txt").write_text(p.stderr)
    if p.returncode != 0:
        raise SystemExit(f"STOP: bridge failed for {T}/{r}: {p.stderr[-1500:]}")

    def tv(label):
        m = re.search(re.escape(label) + r":\s*([0-9.]+)\n", p.stderr)
        return float(m.group(1)) if m else 0.0
    return run_id, {"user_cpu_s": tv("User time (seconds)"), "system_cpu_s": tv("System time (seconds)"),
                    "max_rss_kib": tv("Maximum resident set size (kbytes)")}


def hist_fn(ch, cfg):
    b = cfg["bins"]
    if ch["type"] == "circular":
        o = ch["origin"] + cfg["tau_origin_shift_bins"] * 2 * np.pi / b
        return lambda x: circular_hist(x, b, o)
    return lambda x: linear_hist(x, ch["lo"], ch["hi"], b)[0]


def run():
    spec = json.loads(SPEC.read_text())
    frz = json.loads(FREEZE.read_text())
    bad = {c: h for c, h in frz["code_sha256"].items() if sha256(ROOT / c) != h}
    if frz["spec_sha256"] != sha256(SPEC) or bad:
        raise SystemExit(f"STOP: spec or code differs from the freeze: {sorted(bad)}")
    build = HERE / "build"
    build.mkdir(exist_ok=True)
    t0 = time.time()
    src_hash = sha256(ROOT / spec["dataset"]["file"])
    if src_hash != spec["dataset"]["sha256"]:
        raise SystemExit("STOP: source hash differs from the WP-D record")

    cap = spec["resource_budget"]["cpu_seconds"]
    used, res_bridge, runs = 0.0, {}, {}
    for T, reps in spec["dataset"]["conditions"].items():
        for r in reps:
            run_id, ru = bridge(spec, T, r, build, cap - used)
            used += ru["user_cpu_s"] + ru["system_cpu_s"]
            res_bridge[f"{T}/{r}"] = ru
            pol = spec["pipeline_policies"]
            feats, rep = features_from_export(build / f"{run_id}.npz", build / f"{run_id}.json",
                                              min_segment_angstrom=pol["min_segment_angstrom"],
                                              collinear_sin_tol=pol["collinear_sin_tol"])
            export_features(feats, build / f"{run_id}.npz", build / f"{run_id}.json", build, run_id)
            runs[(T, r)] = {"rep": rep, "a": np.load(build / f"{run_id}.npz"), "f": feats,
                            "fdoc": json.loads((build / f"{run_id}_features.json").read_text())}

    # --- checks -------------------------------------------------------------------------------------------------------
    u = spec["sampling"]["used_frames"]
    expect_frames = list(range(u["start"], u["stop"], u["step"]))
    ident0 = [(x["chain"], x["resid"], x["icode"], x["resname"]) for x in runs[("320", "0")]["rep"]["residues"]]
    checks = {
        "source_hash_matches": True,
        "frames_per_replica": all(len(v["a"]["frame_index"]) == u["per_replica"] for v in runs.values()),
        "source_frame_indices_identical": all(v["a"]["frame_index"].tolist() == expect_frames for v in runs.values()),
        "inferred_time_labelled": all("INFERRED" in v["rep"]["frames_inferred_time"]["basis"] for v in runs.values()),
        "residue_identity_identical": all([(x["chain"], x["resid"], x["icode"], x["resname"])
                                           for x in v["rep"]["residues"]] == ident0 for v in runs.values()),
        "mapping_check_no_mismatch": all(v["rep"]["provenance"]["bridge"]["mapping_check"]["mismatches"] == 0
                                         for v in runs.values()),
        "no_breaks": all(not v["rep"]["breaks"] for v in runs.values()),
        "native_extension_resolved": all(v["fdoc"]["provenance"]["core"]["native_extension"]["status"] == "resolved"
                                         for v in runs.values()),
        "total_frames_within_budget": sum(len(v["a"]["frame_index"]) for v in runs.values())
                                      <= spec["resource_budget"]["max_total_frames"],
    }
    if not all(checks.values()):
        raise SystemExit(f"STOP: checks failed: {checks}")

    # --- channels and matched validity ---------------------------------------------------------------------------------
    def channel_data(v):
        X = v["a"]["coordinates_angstrom"]
        f = v["f"]
        return {"kappa": (f.kappa, f.kappa_valid), "d_i_i2": (window_distances(X, 2), None),
                "tau": (f.tau, f.tau_valid), "d_i_i3": (window_distances(X, 3), None)}

    data = {k: channel_data(v) for k, v in runs.items()}
    valid_windows, excluded = {}, {}
    for p, q in spec["pairs"]:
        n_w = data[("320", "0")][p][0].shape[1]
        ok = np.ones(n_w, bool)
        for k in data:
            ok &= data[k][p][1].all(axis=0) & np.isfinite(data[k][q][0]).all(axis=0)
        valid_windows[p] = valid_windows[q] = np.flatnonzero(ok)
        excluded[p] = excluded[q] = {"windows": int((~ok).sum()), "reason": "chain-end windows (feature undefined) or "
                                     "any invalid frame in any replica for either channel of the pair"}

    res_list = runs[("320", "0")]["rep"]["residues"]

    def label(ch, w):
        a, b = (w - 1, w + 1) if ch in ("kappa", "d_i_i2") else (w - 1, w + 2)
        return f"{res_list[a]['resid']}-{res_list[b]['resid']}"

    configs = [dict(name="primary", **spec["configurations"]["primary"])] + spec["configurations"]["sensitivity"]
    out_cfg, oor = {}, {}
    for cfg in configs:
        per = {}
        for ch, chspec in spec["channels"].items():
            h = hist_fn(chspec, cfg)
            rows = []
            for w in valid_windows[ch]:
                ha = [h(data[("320", r)][ch][0][:, w]) for r in spec["dataset"]["conditions"]["320"]]
                hb = [h(data[("348", r)][ch][0][:, w]) for r in spec["dataset"]["conditions"]["348"]]
                rc = replica_comparison(ha, hb)
                rc.update(window=int(w), residues=label(ch, int(w)))
                if cfg["name"] != "primary":
                    rc.pop("pairs")
                    rc.pop("loro")
                rows.append(rc)
            per[ch] = rows
            if cfg["name"] == "primary" and chspec["type"] == "linear":
                allv = np.concatenate([data[k][ch][0][:, valid_windows[ch]].ravel() for k in data])
                oor[ch] = {"underflow": int((allv < chspec["lo"]).sum()), "overflow": int((allv > chspec["hi"]).sum()),
                           "observations": int(allv.size)}
        out_cfg[cfg["name"]] = per

    prim = out_cfg["primary"]
    counts = {c: {ch: dict(collections.Counter(r["category"] for r in rows)) for ch, rows in per.items()}
              for c, per in out_cfg.items()}
    cross = {}
    for p, q in spec["pairs"]:
        cq = {r["window"]: r["category"] for r in prim[q]}
        cross[f"{p}~{q}"] = dict(collections.Counter(f"{r['category']}|{cq[r['window']]}" for r in prim[p]))
    pooled_spearman = {f"{p}~{q}": spearman([r["pooled"] for r in prim[p]], [r["pooled"] for r in prim[q]])
                       for p, q in spec["pairs"]}
    sens = {}
    for c in out_cfg:
        if c == "primary":
            continue
        sens[c] = {ch: {"same_category_windows": sum(a["category"] == b["category"] for a, b in zip(prim[ch], out_cfg[c][ch])),
                        "windows": len(prim[ch]),
                        "spearman_pooled_vs_primary": spearman([a["pooled"] for a in prim[ch]],
                                                               [b["pooled"] for b in out_cfg[c][ch]])}
                   for ch in prim}
    listing = {ch: {cat: [{"residues": r["residues"], "pooled": r["pooled"], "pooled_loro_range": r["pooled_loro_range"],
                           "fraction": r["fraction_between_above_max_within"],
                           "separating_when_left_out": r["separating_when_left_out"]}
                          for r in sorted(rows, key=lambda x: -x["pooled"]) if r["category"] == cat]
                    for cat in ("consistent", "replica_dependent")} for ch, rows in prim.items()}

    ru = resource.getrusage(resource.RUSAGE_SELF)
    total = used + ru.ru_utime + ru.ru_stime
    resources = {"bridge_stage_per_replica": res_bridge, "bridge_cpu_s": used,
                 "features_and_comparison": {"user_cpu_s": ru.ru_utime, "system_cpu_s": ru.ru_stime,
                                             "max_rss_kib": ru.ru_maxrss},
                 "total_cpu_s": total, "wall_s": time.time() - t0, "within_cpu_budget": total <= cap,
                 "max_rss_kib_bridge": max(v["max_rss_kib"] for v in res_bridge.values())}
    core = runs[("320", "0")]["fdoc"]["provenance"]["core"]
    out = {"id": "DEV-REP-1", "evidence_level": "replicated descriptive development comparison; no test or biological claim",
           "spec_sha256": frz["spec_sha256"], "freeze_utc": frz["frozen_utc"], "source_sha256": src_hash,
           "checks": checks, "frames": {"source_indices": [expect_frames[0], expect_frames[-1]], "per_replica": len(expect_frames),
                                        "inferred_time_ns": [float(expect_frames[0]), float(expect_frames[-1])],
                                        "inferred_time_basis": runs[("320", "0")]["rep"]["frames_inferred_time"]["basis"]},
           "pbc": {k[0] + "/" + k[1]: v["rep"]["pbc"] for k, v in runs.items()},
           "core_provenance": {"wrapper": core["wrapper"], "native_extension": core["native_extension"]},
           "excluded_windows": excluded, "category_counts": counts, "paired_channel_crosstab_primary": cross,
           "pooled_distance_spearman_primary": pooled_spearman, "sensitivity": sens, "windows_by_category_primary": listing,
           "out_of_range_primary": oor, "resources": resources, "per_window_primary": prim,
           "outputs_sha256": {p.name: sha256(p) for p in sorted(build.glob("DEV-REP-1_*")) if p.suffix in (".npz", ".json")}}
    rp = HERE / "results" / "DEV-REP-1_results.json"
    rp.write_text(json.dumps(out, indent=1) + "\n")
    (HERE / "results" / "DEV-REP-1_report.md").write_text(report(out, spec, sha256(rp)))
    print(json.dumps({"checks": checks, "category_counts": counts["primary"], "resources": {
        k: resources[k] for k in ("total_cpu_s", "wall_s", "within_cpu_budget", "max_rss_kib_bridge")}}, indent=1))
    return 0


def report(o, spec, h):
    L = ["# DEV-REP-1: replicated descriptive comparison, mdCATH 1ha8A00, 320 K vs 348 K", "",
         f"Generated from `DEV-REP-1_results.json` (SHA-256 {h[:16]}…). Development evidence only. Five replicas per "
         "temperature, all started from one equilibrated configuration, so they measure reproducibility given that start. "
         "Pairwise replica distances share replicas and are dependent; categories are descriptive, not tests. "
         "320 K vs 348 K is a simulation-condition contrast with no experimental endpoint.", "",
         "## Inputs", "",
         f"- Source SHA-256 {o['source_sha256'][:16]}… (matches the WP-D record). Source frames "
         f"{o['frames']['source_indices'][0]}–{o['frames']['source_indices'][1]} ({o['frames']['per_replica']} per replica; "
         f"frames 0–42 excluded as predeclared). Time: {o['frames']['inferred_time_basis']}",
         f"- Checks: " + ", ".join(f"{k} {v}" for k, v in o["checks"].items()),
         f"- Core: wrapper SHA-256 {o['core_provenance']['wrapper']['sha256'][:16]}…, native extension "
         f"{o['core_provenance']['native_extension']['status']} (SHA-256 {str(o['core_provenance']['native_extension']['sha256'])[:16]}…)",
         f"- Excluded windows per channel: " + ", ".join(f"{k} {v['windows']}" for k, v in o["excluded_windows"].items())
         + " (chain-end windows where the feature is undefined)", "",
         "## Categories (primary, 36 bins)", "",
         "| Channel | Consistent | Replica-dependent | Inconclusive |", "|---|---|---|---|"]
    for ch, c in o["category_counts"]["primary"].items():
        L.append(f"| {ch} | {c.get('consistent', 0)} | {c.get('replica_dependent', 0)} | {c.get('inconclusive', 0)} |")
    L += ["", "*Consistent*: every between-temperature replica-pair distance exceeds every within-temperature one. "
          "*Replica-dependent*: that holds only after leaving out one specific replica. *Inconclusive*: neither; this is "
          "not evidence of no difference.", "",
          "Paired channels on the same windows (category κ|d(i,i+2) and τ|d(i,i+3), primary):", ""]
    for k, v in o["paired_channel_crosstab_primary"].items():
        L.append(f"- {k}: " + ", ".join(f"{a} {n}" for a, n in sorted(v.items())))
    L += ["", "Spearman of pooled 320-vs-348 distances between paired channels: " + ", ".join(
        f"{k} {v:.3f}" for k, v in o["pooled_distance_spearman_primary"].items()), "",
          "κ and d(i,i+2) are near-deterministic functions of each other; d(i,i+3) cannot distinguish τ from −τ. "
          "Agreement or disagreement does not show that one representation is better.", "",
          "## Windows by category (primary; pooled JS distance with its leave-one-replica-out range)", ""]
    for ch, cats in o["windows_by_category_primary"].items():
        for cat, rows in cats.items():
            if rows:
                L.append(f"- {ch}, {cat}: " + ", ".join(
                    f"{r['residues']} ({r['pooled']:.3f} [{r['pooled_loro_range'][0]:.3f}–{r['pooled_loro_range'][1]:.3f}]"
                    + (f"; left out {','.join(r['separating_when_left_out'])}" if r["separating_when_left_out"] else "") + ")"
                    for r in rows))
    L += ["", "## Sensitivity (predeclared; not used for selection)", "",
          "| Configuration | Channel | Same category as primary | Spearman of pooled distance vs primary |",
          "|---|---|---|---|"]
    for c, chs in o["sensitivity"].items():
        for ch, s in chs.items():
            L.append(f"| {c} | {ch} | {s['same_category_windows']}/{s['windows']} | {s['spearman_pooled_vs_primary']:.3f} |")
    r = o["resources"]
    L += ["", "## Resources", "",
          f"- Total CPU {r['total_cpu_s']:.1f} s (budget {spec['resource_budget']['cpu_seconds']} s); wall {r['wall_s']:.1f} s; "
          f"max RSS: bridge stage {r['max_rss_kib_bridge']:.0f} KiB, features and comparison "
          f"{r['features_and_comparison']['max_rss_kib']} KiB", ""]
    return "\n".join(L)


if __name__ == "__main__":
    if len(sys.argv) != 2 or sys.argv[1] not in ("freeze", "run"):
        print(__doc__)
        sys.exit(2)
    sys.exit(freeze() if sys.argv[1] == "freeze" else run())
