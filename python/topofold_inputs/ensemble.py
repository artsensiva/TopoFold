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

ENSEMBLE_VERSION = "0.2.0"
R_MIN = 1e-6            # numerical conditioning threshold for a circular mean direction (not a biological criterion)
CONSTANT_VAR = 1e-20    # sample variance (squared units) at or below which a series is treated as constant
UNAVAILABLE_R = f"unavailable: resultant length below {R_MIN:g} (mean direction numerically ill-conditioned)"


def circular_summary(t) -> dict:
    """Mean direction, resultant length R = |mean(e^{iθ})| and circular SD sqrt(−2 ln R).

    R is clipped to [0, 1] so rounding can never give R > 1. When R < R_MIN the mean direction is numerically
    ill-conditioned and is reported as unavailable (as is the spread, which is then very large). R_MIN is a numerical
    rule only; it says nothing about whether a distribution is biologically meaningful.
    """
    m = np.mean(np.exp(1j * np.asarray(t, dtype=np.float64)))
    R = float(min(max(abs(m), 0.0), 1.0))
    if R < R_MIN:
        return {"center": None, "spread": None, "resultant_length": R, "status": UNAVAILABLE_R}
    return {"center": float(np.angle(m)), "spread": float(math.sqrt(max(-2.0 * math.log(R), 0.0))),
            "resultant_length": R, "status": "ok"}


def _wrap(x):
    return float(math.pi - ((math.pi - x) % (2 * math.pi)))


def autocorrelation(x, circular: bool, max_lag: int):
    """Normalized autocorrelation C(0..L) with L = min(max_lag, n − 1), or (None, reason).

    Linear: C(k) = sum_t u_t u_{t+k} / sum_t u_t², u = x − mean(x).
    Circular: with z_t = e^{iθ_t} and u_t = z_t − mean(z), C(k) = Re sum_t conj(u_t) u_{t+k} / sum_t |u_t|².
    A common rotation θ → θ + α multiplies every u_t by e^{iα}, so C(k) is independent of the angular origin.
    Both use the full-series denominator (C(0) = 1, |C(k)| <= 1); they are descriptive, not estimates with stated
    error.
    """
    x = np.asarray(x, dtype=np.float64)
    n = len(x)
    if n < 2:
        return None, "unavailable: fewer than 2 frames"
    u = np.exp(1j * x) if circular else x.astype(np.complex128)
    u = u - u.mean()
    v = float(np.sum(np.abs(u) ** 2))
    if v / n <= CONSTANT_VAR:
        return None, "unavailable: constant series"
    L = min(max_lag, n - 1)
    return np.array([float(np.real(np.sum(np.conj(u[: n - k]) * u[k:]))) / v for k in range(L + 1)]), "ok"


def temporal_diagnostics(x, circular: bool, max_lag: int) -> dict:
    """Lag-1 autocorrelation and the first lag with C(k) < 1/e (descriptive), with explicit statuses.

    decorrelation_status is "reached", "not reached within L lags" (defined series whose C stays >= 1/e up to the
    largest available lag L), or the reason the diagnostic is unavailable.
    """
    ac, why = autocorrelation(x, circular, max_lag)
    if ac is None:
        return {"lag1": None, "decorrelation_lag": None, "lags_available": 0, "status": why,
                "decorrelation_status": why}
    L = len(ac) - 1
    below = np.flatnonzero(ac[1:] < math.exp(-1))
    return {"lag1": float(ac[1]), "decorrelation_lag": int(below[0] + 1) if below.size else None, "lags_available": L,
            "status": "ok", "decorrelation_status": "reached" if below.size else f"not reached within {L} lags"}


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


def _member_spread(centers, circ: bool) -> dict:
    if len(centers) < 2:
        return {"member_center_spread": None, "member_center_spread_status": "unavailable: single member"}
    if any(c is None for c in centers):
        return {"member_center_spread": None,
                "member_center_spread_status": "unavailable: a member's circular mean is ill-conditioned"}
    if circ:
        cs = circular_summary(np.array(centers))
        return {"member_center_spread": cs["spread"], "member_center_spread_status": cs["status"]}
    return {"member_center_spread": float(np.std(centers, ddof=1)), "member_center_spread_status": "ok"}


def _center_difference(ca, cb, sa, sb, ga, gb, circ: bool) -> dict:
    if ca is None or cb is None:
        bad = ga if ca is None else gb
        return {"center_difference": None,
                "center_difference_status": f"unavailable: center of {bad} {UNAVAILABLE_R.split(': ', 1)[1]}"}
    return {"center_difference": _wrap(ca - cb) if circ else ca - cb, "center_difference_status": "ok"}


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
                msum = [circular_summary(v) if circ else None for v in vals]
                centers = [c["center"] for c in msum] if circ else [float(np.mean(v)) for v in vals]
                d = {"channel": ch, "window": w, "residues": lab, "group": g, "n_members": len(vals),
                     "n_frames": int(pooled.size), "member_frames": json.dumps([int(v.size) for v in vals])}
                if circ:
                    cs = circular_summary(pooled)
                    d.update(center=cs["center"], spread=cs["spread"], resultant_length=cs["resultant_length"],
                             q05=None, q50=None, q95=None, status=cs["status"])
                else:
                    q = np.quantile(pooled, [0.05, 0.5, 0.95])
                    d.update(center=float(np.mean(pooled)),
                             spread=float(np.std(pooled, ddof=1)) if pooled.size > 1 else None,
                             resultant_length=None, q05=float(q[0]), q50=float(q[1]), q95=float(q[2]),
                             status="ok" if pooled.size > 1 else "spread unavailable: fewer than 2 frames")
                rows["distribution"].append(d)
                wj = [js_distance(mh[i], mh[j]) for i in range(len(mh)) for j in range(i + 1, len(mh))]
                rows["replicas"].append({
                    "channel": ch, "window": w, "residues": lab, "group": g, "n_members": len(vals),
                    **_member_spread(centers, circ),
                    "within_member_js_median": float(np.median(wj)) if wj else None,
                    "within_member_js_max": float(np.max(wj)) if wj else None})
                for m, v in zip(groups[g], vals):
                    t = temporal_diagnostics(v, circ, max_lag)
                    rows["temporal"].append({"channel": ch, "window": w, "residues": lab, "group": g,
                                             "member": m["run_id"], "n_frames": int(v.size), **t})
                per_group[g] = (mh, d["center"], d["status"])
            for i, ga in enumerate(names):
                for gb in names[i + 1:]:
                    (ma, ca, sa), (mb, cb, sb) = per_group[ga], per_group[gb]
                    bt = [js_distance(x, y) for x in ma for y in mb]
                    wa = [js_distance(ma[p], ma[q]) for p in range(len(ma)) for q in range(p + 1, len(ma))]
                    wb = [js_distance(mb[p], mb[q]) for p in range(len(mb)) for q in range(p + 1, len(mb))]
                    rows["between"].append({
                        "channel": ch, "window": w, "residues": lab, "group_a": ga, "group_b": gb,
                        "pooled_js": js_distance(np.sum(ma, axis=0), np.sum(mb, axis=0)),
                        "between_member_js_median": float(np.median(bt)),
                        "within_member_js_median_a": float(np.median(wa)) if wa else None,
                        "within_member_js_median_b": float(np.median(wb)) if wb else None,
                        **_center_difference(ca, cb, sa, sb, ga, gb, circ)})
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
                     "ensemble_version": ENSEMBLE_VERSION, "module_sha256": _sha256(__file__)},
        "thresholds": {"circular_resultant_length_min": R_MIN, "constant_series_variance": CONSTANT_VAR},
        "pooling": "member frames concatenated: members weighted by retained frame counts",
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
         "per-window distributions can miss joint changes.", "",
         "**Pooling:** each group's pooled distribution concatenates its members' retained frames, so members are "
         "weighted by their frame counts (`member_frames` in `tables/distribution.csv`).", "",
         f"**Circular summaries:** the resultant length R is reported; mean directions (and differences involving them) "
         f"are shown as unavailable when R < {R_MIN:g}, a numerical conditioning rule, not a biological criterion.", "",
         "## Inputs", ""]
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
                    cd = r["center_difference"]
                    L.append(f"| {ch} | {r['residues']} | {r['pooled_js']:.3f} | {r['between_member_js_median']:.3f} | "
                             f"{'n/a' if wa is None else f'{wa:.3f}'} / {'n/a' if wb is None else f'{wb:.3f}'} | "
                             f"{'n/a' if cd is None else f'{cd:+.3f}'} |")
            L.append("")
    L += ["## Temporal diagnostics (median over member-windows with defined values; descriptive)", "",
          "Autocorrelation C(k): linear channels use the centered series; the circular channel uses centered unit "
          "vectors, C(k) = Re Σ conj(u_t) u_{t+k} / Σ |u_t|², which does not depend on the angular origin. The "
          "decorrelation lag is the first lag with C(k) < 1/e.", "",
          "| Channel | Group | Lag-1 autocorrelation | Decorrelation lag (frames) | Not reached within available lags | "
          "Unavailable (reason: count) |", "|---|---|---|---|---|---|"]
    for ch in CHANNELS:
        for g in names:
            t = [r for r in rows["temporal"] if r["channel"] == ch and r["group"] == g]
            l1 = [r["lag1"] for r in t if r["lag1"] is not None]
            dl = [r["decorrelation_lag"] for r in t if r["decorrelation_lag"] is not None]
            nr = sum(1 for r in t if r["decorrelation_status"].startswith("not reached"))
            un = {}
            for r in t:
                if r["status"] != "ok":
                    un[r["status"]] = un.get(r["status"], 0) + 1
            L.append(f"| {ch} | {g} | {f'{np.median(l1):.3f}' if l1 else 'n/a'} | "
                     f"{np.median(dl) if dl else 'n/a'} | {nr} | "
                     f"{'; '.join(f'{k}: {v}' for k, v in sorted(un.items())) or 'none'} |")
    nun = sum(1 for r in rows["distribution"] if r["status"] != "ok")
    L += ["", f"Distribution rows with an unavailable summary: {nun} (reasons in the `status` column). "
          "Single-member groups have no replica variability (shown as n/a).", ""]
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
