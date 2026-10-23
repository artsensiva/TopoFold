"""Descriptive comparison of named groups of trajectories: residue-mapped tables, distributions, replica variability and
temporal diagnostics. No p-values, significance labels or uncertainty intervals are produced (see README).

    PYTHONPATH=python .venv/bin/python -m topofold_inputs.ensemble SPEC.json --out OUTDIR [--adapter-python PATH]

SPEC.json:
    {"name": "...",
     "groups": {"<group>": [<member>, ...], ...},
     "options": {"bins": 36, "max_lag": 50, "pbc": "none"}}          # optional
    member (MDAnalysis files): {"type": "mdanalysis", "topology": "...", "trajectories": ["..."],
                                "start": 0, "stop": null, "step": 1, "selection": "protein", "altloc": "reject"}
    member (mdCATH HDF5):      {"type": "mdcath", "h5": "...", "domain": "...", "temperature": "320", "replica": "0",
                                "start": 43, "stop": 430}

Each member is read through the existing adapter (MDAnalysis) or mdCATH bridge, then κ/τ come from the TopoFold core
via ``features``. Every member must have the same residue identity (chain, resid, icode, resname) in the same order.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import json
import math
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np

from .compare import circular_hist, js_distance, linear_hist, window_distances

PKG_ROOT = Path(__file__).resolve().parents[1]
CHANNELS = {
    "kappa": {"kind": "linear", "lo": 0.0, "hi": math.pi, "unit": "rad", "span": (1, 1)},
    "tau": {"kind": "circular", "unit": "rad", "span": (1, 2)},
    "d_i_i2": {"kind": "linear", "lo": 3.0, "hi": 8.0, "unit": "angstrom", "span": (1, 1)},
    "d_i_i3": {"kind": "linear", "lo": 3.0, "hi": 12.0, "unit": "angstrom", "span": (1, 2)},
}
DEFAULTS = {"bins": 36, "max_lag": 50, "pbc": "none"}


# --- extraction (existing components) -------------------------------------------------------------------------------

def _sha256(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _run(cmd, adapter_python):
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(PKG_ROOT), os.environ.get("PYTHONPATH", "")]))
    p = subprocess.run([adapter_python, *cmd], env=env, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"extraction failed: {' '.join(cmd)}\n{p.stderr[-2000:]}")


def extract_member(member: dict, run_id: str, work: Path, adapter_python: str, options: dict) -> dict:
    """Adapter or bridge in the adapter interpreter, then κ/τ here. Returns arrays and identity."""
    from .features import export_features, features_from_export

    def slice_args(m):
        out = []
        for k in ("start", "stop", "step"):
            if m.get(k) is not None:
                out += [f"--{k}", str(m[k])]
        return out

    if member["type"] == "mdanalysis":
        cmd = ["-m", "topofold_inputs", member["topology"], *member.get("trajectories", []), "--run-id", run_id,
               "--out", str(work), "--selection", member.get("selection", "protein"),
               "--altloc", member.get("altloc", "reject"), "--pbc", member.get("pbc", options["pbc"]), *slice_args(member)]
    elif member["type"] == "mdcath":
        cmd = ["-m", "topofold_inputs.mdcath_bridge", member["h5"], member["domain"], str(member["temperature"]),
               str(member["replica"]), "--run-id", run_id, "--out", str(work),
               "--pbc", member.get("pbc", options["pbc"]), *slice_args(member)]
    else:
        raise ValueError(f"unknown member type {member['type']!r}")
    _run(cmd, adapter_python)
    npz, js = work / f"{run_id}.npz", work / f"{run_id}.json"
    feats, rep = features_from_export(npz, js)
    export_features(feats, npz, js, work, run_id)
    a = np.load(npz)
    fdoc = json.loads((work / f"{run_id}_features.json").read_text())
    return {"run_id": run_id, "coords": a["coordinates_angstrom"], "present": a["present"],
            "frame_index": a["frame_index"], "kappa": feats.kappa, "kappa_valid": feats.kappa_valid,
            "tau": feats.tau, "tau_valid": feats.tau_valid, "residues": rep["residues"],
            "provenance": {"member_spec": member, "source_files": rep["provenance"]["source_files"],
                           "adapter_json_sha256": _sha256(js),
                           "features_json_sha256": _sha256(work / f"{run_id}_features.json"),
                           "core": fdoc["provenance"]["core"], "frame_identity": rep["frame_identity"]["status"],
                           "units": rep["units"]}}


# --- descriptive statistics ------------------------------------------------------------------------------------------

def _circ_mean(t):
    return float(np.angle(np.mean(np.exp(1j * t))))


def _circ_sd(t):
    R = abs(np.mean(np.exp(1j * t)))
    return float(math.sqrt(-2.0 * math.log(R))) if R > 0 else float("inf")


def _wrap(x):
    return float(math.pi - ((math.pi - x) % (2 * math.pi)))


def _autocorr(x, max_lag):
    x = np.asarray(x, dtype=np.float64) - np.mean(x)
    v = float(np.dot(x, x))
    n = len(x)
    if v == 0:
        return np.full(min(max_lag, n - 1) + 1, np.nan)
    return np.array([np.dot(x[: n - k], x[k:]) / v for k in range(min(max_lag, n - 1) + 1)])


def temporal_diagnostics(x, circular: bool, max_lag: int) -> dict:
    """Lag-1 autocorrelation and first lag with autocorrelation below 1/e (descriptive). Circular: mean of the
    autocorrelation functions of cos θ and sin θ."""
    ac = (0.5 * (_autocorr(np.cos(x), max_lag) + _autocorr(np.sin(x), max_lag))) if circular else _autocorr(x, max_lag)
    below = np.flatnonzero(ac[1:] < math.exp(-1))
    return {"lag1": None if np.isnan(ac[1]) else float(ac[1]),
            "decorrelation_lag": int(below[0] + 1) if below.size else None, "max_lag": int(len(ac) - 1)}


def _hist(x, ch, bins):
    spec = CHANNELS[ch]
    if spec["kind"] == "circular":
        return circular_hist(x, bins, -math.pi)
    return linear_hist(x, spec["lo"], spec["hi"], bins)[0]


def channel_arrays(m: dict) -> dict:
    """Per channel: (frames, windows) values and validity (distances valid where finite and the window lies in one
    adapter segment, taken from the matching κ/τ validity)."""
    d2, d3 = window_distances(m["coords"], 2), window_distances(m["coords"], 3)
    return {"kappa": (m["kappa"], m["kappa_valid"]), "tau": (m["tau"], m["tau_valid"]),
            "d_i_i2": (d2, m["kappa_valid"] & np.isfinite(d2)), "d_i_i3": (d3, m["tau_valid"] & np.isfinite(d3))}


def window_label(residues, ch, w):
    lo, hi = CHANNELS[ch]["span"]
    a, b = w - lo, w + hi
    if a < 0 or b >= len(residues):
        return None
    ra, rb = residues[a], residues[b]
    return f"{ra['chain']}:{ra['resid']}{ra['icode'] or ''}-{rb['resid']}{rb['icode'] or ''}"


def analyze(groups: dict[str, list[dict]], options: dict) -> dict:
    """All tables from extracted members. Pure function of the arrays (deterministic)."""
    names = list(groups)
    members = [m for g in names for m in groups[g]]
    ident = [(r["chain"], r["resid"], r["icode"], r["resname"]) for r in members[0]["residues"]]
    for m in members[1:]:
        if [(r["chain"], r["resid"], r["icode"], r["resname"]) for r in m["residues"]] != ident:
            raise ValueError(f"residue identity of {m['run_id']} differs from {members[0]['run_id']}; "
                             "groups must share one residue mapping")
    residues = members[0]["residues"]
    bins, max_lag = options["bins"], options["max_lag"]
    arrays = {m["run_id"]: channel_arrays(m) for m in members}
    rows = {"distribution": [], "replicas": [], "between": [], "temporal": [], "exclusions": []}
    hists = {}
    for ch, spec in CHANNELS.items():
        circ = spec["kind"] == "circular"
        n_win = arrays[members[0]["run_id"]][ch][0].shape[1]
        for w in range(n_win):
            lab = window_label(residues, ch, w)
            if lab is None:
                continue
            bad = {m["run_id"]: int((~arrays[m["run_id"]][ch][1][:, w]).sum()) for m in members}
            if any(bad.values()):
                rows["exclusions"].append({"channel": ch, "window": w, "residues": lab,
                                           "reason": "invalid frames (missing residue, break or degenerate geometry)",
                                           "invalid_frames_by_member": json.dumps({k: v for k, v in bad.items() if v})})
                continue
            per_group = {}
            for g in names:
                vals = [arrays[m["run_id"]][ch][0][:, w] for m in groups[g]]
                mh = [_hist(v, ch, bins) for v in vals]
                pooled = np.concatenate(vals)
                hists[(ch, w, g)] = (np.sum(mh, axis=0), np.array(mh))
                centers = [(_circ_mean(v) if circ else float(np.mean(v))) for v in vals]
                d = {"channel": ch, "window": w, "residues": lab, "group": g, "n_members": len(vals),
                     "n_frames": int(pooled.size)}
                if circ:
                    d.update(center=_circ_mean(pooled), spread=_circ_sd(pooled), q05=None, q50=None, q95=None)
                else:
                    q = np.quantile(pooled, [0.05, 0.5, 0.95])
                    d.update(center=float(np.mean(pooled)), spread=float(np.std(pooled, ddof=1)),
                             q05=float(q[0]), q50=float(q[1]), q95=float(q[2]))
                rows["distribution"].append(d)
                wj = [js_distance(mh[i], mh[j]) for i in range(len(mh)) for j in range(i + 1, len(mh))]
                rows["replicas"].append({
                    "channel": ch, "window": w, "residues": lab, "group": g, "n_members": len(vals),
                    "member_center_spread": (None if len(vals) < 2 else
                                             (_circ_sd(np.array(centers)) if circ else float(np.std(centers, ddof=1)))),
                    "within_member_js_median": float(np.median(wj)) if wj else None,
                    "within_member_js_max": float(np.max(wj)) if wj else None})
                for m, v in zip(groups[g], vals):
                    t = temporal_diagnostics(v, circ, max_lag)
                    rows["temporal"].append({"channel": ch, "window": w, "residues": lab, "group": g,
                                             "member": m["run_id"], "n_frames": int(v.size), **t})
                per_group[g] = (mh, d["center"])
            for i, ga in enumerate(names):
                for gb in names[i + 1:]:
                    (ma, ca), (mb, cb) = per_group[ga], per_group[gb]
                    bt = [js_distance(x, y) for x in ma for y in mb]
                    wa = [js_distance(ma[p], ma[q]) for p in range(len(ma)) for q in range(p + 1, len(ma))]
                    wb = [js_distance(mb[p], mb[q]) for p in range(len(mb)) for q in range(p + 1, len(mb))]
                    rows["between"].append({
                        "channel": ch, "window": w, "residues": lab, "group_a": ga, "group_b": gb,
                        "pooled_js": js_distance(np.sum(ma, axis=0), np.sum(mb, axis=0)),
                        "between_member_js_median": float(np.median(bt)),
                        "within_member_js_median_a": float(np.median(wa)) if wa else None,
                        "within_member_js_median_b": float(np.median(wb)) if wb else None,
                        "center_difference": _wrap(ca - cb) if circ else ca - cb})
    return {"rows": rows, "hists": hists, "residues": residues}


# --- output ----------------------------------------------------------------------------------------------------------

def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return "inf" if math.isinf(v) else f"{v:.6g}"
    return v


def write_outputs(res: dict, groups: dict, spec: dict, options: dict, out: Path, spec_path: Path | None) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    tables = out / "tables"
    tables.mkdir(exist_ok=True)
    for name, rows in res["rows"].items():
        path = tables / f"{name}.csv"
        cols = list(rows[0].keys()) if rows else ["channel"]
        with open(path, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
            wr.writeheader()
            for r in rows:
                wr.writerow({k: _fmt(v) for k, v in r.items()})
    arrs = {}
    for (ch, w, g), (pooled, per) in sorted(res["hists"].items()):
        arrs[f"{ch}__w{w}__{g}__pooled"] = pooled
        arrs[f"{ch}__w{w}__{g}__members"] = per
    for ch, s in CHANNELS.items():
        arrs[f"{ch}__bin_edges"] = (np.linspace(-math.pi, math.pi, options["bins"] + 1) if s["kind"] == "circular"
                                    else np.linspace(s["lo"], s["hi"], options["bins"] + 1))
    np.savez(out / "distributions.npz", **arrs)
    manifest = {
        "schema": "topofold_inputs.ensemble/1", "name": spec.get("name"),
        "statement": "Descriptive comparison only: no p-values, significance labels or uncertainty intervals.",
        "spec_sha256": _sha256(spec_path) if spec_path else None, "options": options,
        "channels": {k: {kk: vv for kk, vv in v.items() if kk != "span"} for k, v in CHANNELS.items()},
        "groups": {g: [{"run_id": m["run_id"], "frames": int(len(m["frame_index"])),
                        "frame_index_first_last": [int(m["frame_index"][0]), int(m["frame_index"][-1])],
                        **m["provenance"]} for m in ms] for g, ms in groups.items()},
        "outputs_sha256": {p.name: _sha256(p) for p in sorted(tables.glob("*.csv"))}
        | {"distributions.npz": _sha256(out / "distributions.npz")},
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "module_sha256": _sha256(__file__)},
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    (out / "report.md").write_text(report(res, groups, manifest))
    return manifest


def report(res: dict, groups: dict, manifest: dict) -> str:
    rows = res["rows"]
    names = list(groups)
    L = [f"# Descriptive ensemble comparison: {manifest.get('name') or 'unnamed'}", "",
         "Generated from the tables in `tables/`. **Descriptive only:** no p-values, significance labels or "
         "uncertainty intervals. Differences between groups can reflect geometry, dynamics, sampling or starting "
         "conditions; replica variability and temporal diagnostics are shown so that readers can judge this. Marginal "
         "per-window distributions can miss joint changes.", "", "## Inputs", ""]
    for g, ms in manifest["groups"].items():
        L.append(f"- **{g}**: {len(ms)} member(s): " + "; ".join(
            f"{m['run_id']} ({m['frames']} frames, source frames {m['frame_index_first_last'][0]}–"
            f"{m['frame_index_first_last'][1]}, frame identity {m['frame_identity']})" for m in ms))
    core = next(iter(manifest["groups"].values()))[0]["core"]
    L += [f"- Core: native extension {core['native_extension']['status']} "
          f"(SHA-256 {str(core['native_extension']['sha256'])[:16]}…); analysis options {manifest['options']} "
          "(per-member reading options, e.g. frame slice and PBC policy, are in `manifest.json`)", "",
          "## Exclusions", "",
          f"{len(rows['exclusions'])} channel–window combinations excluded (invalid frames in at least one member); "
          "listed in `tables/exclusions.csv`.", ""]
    for i, ga in enumerate(names):
        for gb in names[i + 1:]:
            L += [f"## {ga} vs {gb}: largest pooled JS distances per channel (descriptive)", "",
                  "| Channel | Window | Pooled JS | Between-member JS (median) | Within-member JS median "
                  f"{ga} / {gb} | Center difference |", "|---|---|---|---|---|---|"]
            for ch in CHANNELS:
                sel = sorted((r for r in rows["between"] if r["channel"] == ch and r["group_a"] == ga
                              and r["group_b"] == gb), key=lambda r: -r["pooled_js"])[:5]
                for r in sel:
                    wa, wb = r["within_member_js_median_a"], r["within_member_js_median_b"]
                    L.append(f"| {ch} | {r['residues']} | {r['pooled_js']:.3f} | {r['between_member_js_median']:.3f} | "
                             f"{'n/a' if wa is None else f'{wa:.3f}'} / {'n/a' if wb is None else f'{wb:.3f}'} | "
                             f"{r['center_difference']:+.3f} |")
            L.append("")
    L += ["## Temporal diagnostics (median over windows; descriptive)", "",
          "| Channel | Group | Lag-1 autocorrelation | Decorrelation lag (frames; first lag < 1/e) |", "|---|---|---|---|"]
    for ch in CHANNELS:
        for g in names:
            t = [r for r in rows["temporal"] if r["channel"] == ch and r["group"] == g]
            l1 = [r["lag1"] for r in t if r["lag1"] is not None]
            dl = [r["decorrelation_lag"] for r in t if r["decorrelation_lag"] is not None]
            beyond = sum(1 for r in t if r["decorrelation_lag"] is None)
            L.append(f"| {ch} | {g} | {np.median(l1):.3f} | {np.median(dl) if dl else 'n/a'}"
                     f"{f' ({beyond} member-windows beyond max lag)' if beyond else ''} |" if l1 else
                     f"| {ch} | {g} | n/a | n/a |")
    L += ["", "Single-member groups have no replica variability (shown as n/a).", ""]
    return "\n".join(L)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m topofold_inputs.ensemble", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spec")
    ap.add_argument("--out", required=True)
    ap.add_argument("--adapter-python", default=os.environ.get("TOPOFOLD_ADAPTER_PYTHON"))
    a = ap.parse_args(argv)
    spec_path = Path(a.spec)
    spec = json.loads(spec_path.read_text())
    options = {**DEFAULTS, **spec.get("options", {})}
    adapter_python = a.adapter_python or (sys.executable if importlib.util.find_spec("MDAnalysis") else None)
    if adapter_python is None:
        print("no interpreter with MDAnalysis: pass --adapter-python", file=sys.stderr)
        return 2
    out = Path(a.out)
    work = out / "members"
    work.mkdir(parents=True, exist_ok=True)
    groups = {}
    for g, members in spec["groups"].items():
        groups[g] = [extract_member(m, f"{g}_{k}", work, adapter_python, options) for k, m in enumerate(members)]
    res = analyze(groups, options)
    write_outputs(res, groups, spec, options, out, spec_path)
    n = {k: len(v) for k, v in res["rows"].items()}
    print(f"groups {list(groups)}; rows {n}; wrote {out}/report.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
