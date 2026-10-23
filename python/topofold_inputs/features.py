"""Residue-attributed Cα curvature (κ) and torsion (τ) from adapter output, computed by the TopoFold core.

Values come from ``topofold.compute_invariants`` (the existing Rust implementation); this module only decides
*where* the core may be applied and records why every other cell is undefined. It introduces no estimator.

Definitions (core conventions, radians):
* κ at residue i: turning angle at vertex i between the Cα(i−1)→Cα(i) and Cα(i)→Cα(i+1) directions, in [0, π].
* τ at hinge i (bond i→i+1): signed dihedral over Cα(i−1), Cα(i), Cα(i+1), Cα(i+2), in (−π, π].

A feature is computed only when all its atoms are present and lie in one adapter segment (no chain boundary,
numbering discontinuity, missing residue or over-long Cα–Cα distance between them), no segment it uses is shorter
than ``min_segment_angstrom``, and no vertex whose osculating plane it needs has |sin κ| below ``collinear_sin_tol``.
The core rejects a whole trace for such geometry, so runs are split there; nothing is interpolated or invented.

The core binding takes float32 input. Adapter coordinates are float32 values from MDAnalysis (possibly shifted by
whole box vectors under ``pbc=make_whole_sequential``); the float32 rounding applied here is recorded.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import platform
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

import numpy as np

FEATURES_VERSION = "0.1.1"

REASONS = {
    0: "valid",
    1: "missing_coordinates",
    2: "terminal_or_break",
    3: "degenerate_segment",
    4: "near_collinear",
    5: "run_too_short_for_core",
    6: "core_error",
    7: "non_finite_core_output",
}
VALID, MISSING, TERMINAL, DEGENERATE, COLLINEAR, SHORT, CORE_ERROR, NONFINITE = range(8)


@dataclass
class CAFeatures:
    kappa: np.ndarray  # (F, N) radians, NaN where invalid
    kappa_reason: np.ndarray  # (F, N) uint8, REASONS
    tau: np.ndarray  # (F, N-1) radians, hinge i = bond i -> i+1, NaN where invalid
    tau_reason: np.ndarray  # (F, N-1) uint8
    kappa_atoms: np.ndarray  # (N, 3) residue indices (i-1, i, i+1), -1 outside the chain
    tau_atoms: np.ndarray  # (N-1, 4) residue indices (i-1, i, i+1, i+2), -1 outside the chain
    config: dict[str, Any] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)

    @property
    def kappa_valid(self) -> np.ndarray:
        return self.kappa_reason == VALID

    @property
    def tau_valid(self) -> np.ndarray:
        return self.tau_reason == VALID


def compute_ca_features(
    coordinates: np.ndarray,
    present: np.ndarray,
    segments: Sequence[Sequence[int]],
    *,
    min_segment_angstrom: float = 1e-6,
    collinear_sin_tol: float = 1e-6,
    core=None,
) -> CAFeatures:
    """κ and τ for every frame, computed by the core only inside adapter segments on non-degenerate runs."""
    X = np.asarray(coordinates, dtype=np.float64)
    if X.ndim != 3 or X.shape[2] != 3:
        raise ValueError(f"coordinates must be (frames, residues, 3); got {X.shape}")
    F, N, _ = X.shape
    present = np.asarray(present, dtype=bool)
    if present.shape != (N,):
        raise ValueError("present must have one entry per residue")
    _check_segments(segments, present, N)
    if np.isnan(X[:, present]).any() or not np.isfinite(X[:, present]).all():
        raise ValueError("non-finite coordinates for present residues")
    if not (min_segment_angstrom > 0 and collinear_sin_tol > 0):
        raise ValueError("tolerances must be positive")
    if core is None:
        import topofold as core

    idx = np.arange(N)
    kappa_atoms = np.stack([idx - 1, idx, idx + 1], axis=1)
    kappa_atoms[kappa_atoms >= N] = -1
    hinge = np.arange(max(N - 1, 0))
    tau_atoms = np.stack([hinge - 1, hinge, hinge + 1, hinge + 2], axis=1) if N > 1 else np.zeros((0, 4), int)
    tau_atoms[tau_atoms >= N] = -1

    kappa = np.full((F, N), np.nan)
    tau = np.full((F, max(N - 1, 0)), np.nan)
    k_reason = np.full((F, N), TERMINAL, dtype=np.uint8)
    t_reason = np.full((F, max(N - 1, 0)), TERMINAL, dtype=np.uint8)
    k_reason[:, ~present] = MISSING
    for h in range(N - 1):
        if any(a >= 0 and not present[a] for a in tau_atoms[h]):
            t_reason[:, h] = MISSING  # atoms beyond the chain ends stay TERMINAL

    Xf32 = X.astype(np.float32)
    rounding = float(np.nanmax(np.abs(Xf32.astype(np.float64) - X))) if present.any() else 0.0
    core_errors: list[dict[str, Any]] = []
    for f in range(F):
        P = Xf32[f].astype(np.float64)  # exactly what the core receives
        for s, e in segments:
            _segment_features(core, Xf32[f], P, s, e, kappa[f], tau[f], k_reason[f], t_reason[f],
                              min_segment_angstrom, collinear_sin_tol, core_errors, f)

    return CAFeatures(
        kappa=kappa, kappa_reason=k_reason, tau=tau, tau_reason=t_reason,
        kappa_atoms=kappa_atoms, tau_atoms=tau_atoms,
        config={"min_segment_angstrom": min_segment_angstrom, "collinear_sin_tol": collinear_sin_tol,
                "core_input_dtype": "float32"},
        diagnostics={"max_float32_rounding_angstrom": rounding, "core_errors": core_errors[:20],
                     "n_core_errors": len(core_errors)},
    )


def _segment_features(core, Xf32, P, s, e, kappa, tau, k_reason, t_reason, min_seg, sin_tol, core_errors, frame):
    n = e - s
    if n < 3:
        return
    d = np.diff(P[s:e], axis=0)
    length = np.linalg.norm(d, axis=1)
    deg_seg = length < min_seg  # segment k joins local atoms k and k+1
    sin_k = np.zeros(n)
    with np.errstate(invalid="ignore", divide="ignore"):
        u = d / length[:, None]
        sin_k[1:-1] = np.linalg.norm(np.cross(u[:-1], u[1:]), axis=1)
    collinear = np.zeros(n, dtype=bool)
    collinear[1:-1] = (sin_k[1:-1] < sin_tol) & ~deg_seg[:-1] & ~deg_seg[1:]

    # κ reasons for internal vertices of the segment (endpoints stay TERMINAL)
    for v in range(1, n - 1):
        g = s + v
        if deg_seg[v - 1] or deg_seg[v]:
            k_reason[g] = DEGENERATE
        elif collinear[v]:
            k_reason[g] = COLLINEAR
        else:
            k_reason[g] = SHORT  # overwritten below when a run of >= 4 atoms covers it
    # τ reasons for hinges whose four atoms lie in the segment
    for v in range(1, n - 2):
        g = s + v
        if deg_seg[v - 1] or deg_seg[v] or deg_seg[v + 1]:
            t_reason[g] = DEGENERATE
        elif collinear[v] or collinear[v + 1]:
            t_reason[g] = COLLINEAR
        else:
            t_reason[g] = SHORT

    # split into runs: between atoms of a degenerate segment, and at (sharing) each collinear vertex
    runs, a = [], 0
    for k in range(n - 1):
        if deg_seg[k]:
            runs.append((a, k + 1))
            a = k + 1
        elif k + 1 < n - 1 and collinear[k + 1]:
            runs.append((a, k + 2))
            a = k + 1
    runs.append((a, n))
    for a, b in runs:
        if b - a < 4:
            continue
        try:
            out = core.compute_invariants(np.ascontiguousarray(Xf32[s + a: s + b]))
        except Exception as exc:  # the core rejected this run: record, never substitute
            for v in range(a + 1, b - 1):
                k_reason[s + v] = CORE_ERROR
            for v in range(a + 1, b - 2):
                t_reason[s + v] = CORE_ERROR
            core_errors.append({"frame": frame, "residues": [s + a, s + b], "error": str(exc)})
            continue
        kap, ta = np.asarray(out[0], dtype=np.float64), np.asarray(out[1], dtype=np.float64)
        if kap.shape != (b - a - 2,) or ta.shape != (b - a - 3,):
            raise RuntimeError(f"unexpected core output shapes {kap.shape}, {ta.shape} for run of {b - a} atoms")
        for j, val in enumerate(kap):
            g = s + a + 1 + j
            if math.isfinite(val):
                kappa[g], k_reason[g] = val, VALID
            else:
                k_reason[g] = NONFINITE
        for j, val in enumerate(ta):
            g = s + a + 1 + j
            if math.isfinite(val):
                tau[g], t_reason[g] = val, VALID
            else:
                t_reason[g] = NONFINITE


def _check_segments(segments, present, N):
    covered = np.zeros(N, dtype=bool)
    last = -1
    for s, e in segments:
        if not (0 <= s < e <= N) or s < last:
            raise ValueError(f"invalid or overlapping segment [{s}, {e})")
        if not present[s:e].all():
            raise ValueError(f"segment [{s}, {e}) contains residues without coordinates")
        covered[s:e] = True
        last = e
    if (covered != present).any():
        raise ValueError("segments must cover exactly the present residues")


# --- export --------------------------------------------------------------------------------------------------------

def features_from_export(npz_path: str | Path, json_path: str | Path, **kwargs) -> tuple[CAFeatures, dict[str, Any]]:
    """Load adapter output (checking its recorded npz hash) and compute features."""
    rep = json.loads(Path(json_path).read_text())
    if rep.get("export", {}).get("npz_sha256") != _sha256(npz_path):
        raise ValueError("adapter npz does not match the SHA-256 recorded in its JSON report")
    data = np.load(npz_path)
    feats = compute_ca_features(data["coordinates_angstrom"], data["present"], rep["segments_half_open"], **kwargs)
    return feats, rep


def export_features(feats: CAFeatures, adapter_npz: str | Path, adapter_json: str | Path, out_dir: str | Path,
                    stem: str, core=None) -> tuple[Path, Path, Path]:
    """Write ``<stem>_features.npz``, ``<stem>_features.json`` and ``<stem>_report.md``."""
    if core is None:
        import topofold as core
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rep = json.loads(Path(adapter_json).read_text())
    data = np.load(adapter_npz)
    npz = out / f"{stem}_features.npz"
    np.savez(npz, kappa_rad=feats.kappa, kappa_valid=feats.kappa_valid, kappa_reason=feats.kappa_reason,
             kappa_atoms=feats.kappa_atoms, tau_rad=feats.tau, tau_valid=feats.tau_valid,
             tau_reason=feats.tau_reason, tau_atoms=feats.tau_atoms, frame_index=data["frame_index"],
             time_ps=data["time_ps"], chain=data["chain"], resid=data["resid"], icode=data["icode"],
             label=data["label"])
    counts = {"kappa": _reason_counts(feats.kappa_reason), "tau": _reason_counts(feats.tau_reason)}
    nonfinite_valid = int((~np.isfinite(feats.kappa[feats.kappa_valid])).sum()
                          + (~np.isfinite(feats.tau[feats.tau_valid])).sum())
    doc = {
        "schema": "topofold_inputs.features/2",
        "engineering_note": "Engineering output: geometric descriptors only; no state, pocket or allostery inference.",
        "run_id": rep["run_id"],
        "units": {"kappa": "radian, [0, pi]", "tau": "radian, (-pi, pi]", "time": rep["units"].get("time")},
        "definitions": {
            "kappa[f, i]": "turning angle at residue i from Ca(i-1), Ca(i), Ca(i+1); atoms in kappa_atoms[i]",
            "tau[f, i]": "dihedral Ca(i-1), Ca(i), Ca(i+1), Ca(i+2) about bond i->i+1; atoms in tau_atoms[i]",
        },
        "reason_codes": {str(k): v for k, v in REASONS.items()},
        "shape": {"frames": int(feats.kappa.shape[0]), "residues": int(feats.kappa.shape[1]),
                  "hinges": int(feats.tau.shape[1])},
        "counts": counts,
        "non_finite_valid_values": nonfinite_valid,
        "config": feats.config,
        "diagnostics": feats.diagnostics,
        "inputs": {"adapter_npz": Path(adapter_npz).name, "adapter_npz_sha256": _sha256(adapter_npz),
                   "adapter_json": Path(adapter_json).name, "adapter_json_sha256": _sha256(adapter_json),
                   "source_files": rep["provenance"]["source_files"]},
        "provenance": {
            "features": {"version": FEATURES_VERSION, "source_sha256": _sha256(__file__)},
            "core": core_provenance(core),
            "software": {"python": platform.python_version(), "numpy": np.__version__},
            "created_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        },
    }
    npz_hash = _sha256(npz)
    doc["export"] = {"npz": npz.name, "npz_sha256": npz_hash}
    js = out / f"{stem}_features.json"
    js.write_text(json.dumps(doc, indent=2, allow_nan=False) + "\n")
    md = out / f"{stem}_report.md"
    md.write_text(_inspection_report(feats, rep, doc))
    return npz, js, md


def _reason_counts(reason: np.ndarray) -> dict[str, int]:
    return {REASONS[k]: int((reason == k).sum()) for k in REASONS}


def core_provenance(core) -> dict[str, Any]:
    """Identify the Python wrapper and, separately, the native extension that defines ``compute_invariants``.

    The native module is taken from the function's own ``__module__`` and accepted only if Python loaded it with an
    extension-module loader from a file with an extension-module suffix. Nothing is inferred from file names; if this
    identification fails, the native entry says ``unresolved`` and carries no hash.
    """
    import importlib.machinery
    import sys

    wrapper_file = getattr(core, "__file__", None)
    wrapper = {"module": getattr(core, "__name__", None), "path": wrapper_file,
               "sha256": _sha256(wrapper_file) if wrapper_file and Path(wrapper_file).is_file() else None}
    fn = getattr(core, "compute_invariants", None)
    mod_name = getattr(fn, "__module__", None)
    mod = sys.modules.get(mod_name) if mod_name else None
    native: dict[str, Any] = {"module": mod_name, "path": None, "sha256": None}
    spec = getattr(mod, "__spec__", None)
    path = getattr(mod, "__file__", None)
    if fn is None:
        native.update(status="unresolved", reason="core has no compute_invariants")
    elif mod is None:
        native.update(status="unresolved", reason=f"defining module {mod_name!r} is not loaded")
    elif not isinstance(getattr(spec, "loader", None), importlib.machinery.ExtensionFileLoader):
        native.update(status="unresolved", reason=f"defining module {mod_name!r} is not a loaded extension module")
    elif not (path and Path(path).is_file()
              and any(str(path).endswith(sfx) for sfx in importlib.machinery.EXTENSION_SUFFIXES)):
        native.update(status="unresolved", reason=f"extension module file {path!r} not found or unexpected suffix")
    else:
        native.update(status="resolved", path=str(Path(path).resolve()), sha256=_sha256(path),
                      loader="ExtensionFileLoader")
    return {"version": _core_version(core), "function": "compute_invariants", "wrapper": wrapper,
            "native_extension": native}


def _core_version(core) -> str | None:
    v = getattr(core, "__version__", None)
    if v:
        return str(v)
    try:
        from importlib.metadata import version
        return version("topofold")
    except Exception:
        return None


def _inspection_report(feats: CAFeatures, rep: dict[str, Any], doc: dict[str, Any]) -> str:
    s = rep["summary"]
    br: dict[str, int] = {}
    for b in rep["breaks"]:
        for r in b["reasons"]:
            br[r] = br.get(r, 0) + 1
    kv, tv = feats.kappa[feats.kappa_valid], feats.tau[feats.tau_valid]
    lines = [
        f"# Cα feature inspection report: {rep['run_id']}",
        "",
        "Generated from the exported results files. Engineering output: descriptive geometry only; no state, "
        "pocket, allostery or biological inference.",
        "",
        "## Input",
        f"- Source files: " + ", ".join(f"`{Path(f['path']).name}` (SHA-256 {str(f['sha256'])[:16]}…)"
                                        for f in rep["provenance"]["source_files"]),
        f"- Reader: {rep['provenance']['reader']['class']}; units: coordinates {rep['units']['coordinates']} "
        f"({rep['units']['basis']}); frame identity: {rep['frame_identity']['status']}",
        f"- Frames: {s['frames']} (source frame indices {rep['frames']['frame_index'][0]}–"
        f"{rep['frames']['frame_index'][-1]}); residues: {s['residues']} ({s['present']} with Cα, "
        f"{s['missing']} without)",
        f"- Breaks: {s['breaks']} " + (str(dict(sorted(br.items()))) if br else "") + f"; segments: {s['segments']}",
        f"- Altloc decisions: {len(rep['altloc_decisions'])}; PBC policy: {rep['pbc']['policy']} "
        f"(frames changed: {rep['pbc']['frames_changed']})",
        "",
        "## Feature validity (cells = frames × positions)",
        "",
        "| Reason | κ | τ |",
        "|---|---|---|",
    ]
    for k, name in REASONS.items():
        lines.append(f"| {name} | {doc['counts']['kappa'][name]} | {doc['counts']['tau'][name]} |")
    lines += [
        "",
        f"- Non-finite values among valid cells: {doc['non_finite_valid_values']}",
        f"- Core errors: {feats.diagnostics['n_core_errors']}; max float32 rounding of core input: "
        f"{feats.diagnostics['max_float32_rounding_angstrom']:.3g} Å",
        f"- κ over valid cells (rad): min {kv.min():.4f}, median {np.median(kv):.4f}, max {kv.max():.4f}"
        if kv.size else "- κ: no valid cells",
        f"- τ over valid cells (rad): min {tv.min():.4f}, max {tv.max():.4f}; circular mean resultant length "
        f"{abs(np.mean(np.exp(1j * tv))):.4f}" if tv.size else "- τ: no valid cells",
        "",
        "## Positions never valid in any frame",
        "",
    ]
    res = rep["residues"]
    never_k = [i for i in range(feats.kappa.shape[1]) if not feats.kappa_valid[:, i].any()]
    never_t = [i for i in range(feats.tau.shape[1]) if not feats.tau_valid[:, i].any()]

    def name(i):
        r = res[i]
        return f"{r['chain']}:{r['resid']}{r['icode'] or ''}"

    lines.append(f"- κ ({len(never_k)}): " + (", ".join(
        f"{name(i)} [{REASONS[int(feats.kappa_reason[0, i])]}]" for i in never_k[:30]) or "none")
        + (" …" if len(never_k) > 30 else ""))
    lines.append(f"- τ hinges ({len(never_t)}): " + (", ".join(
        f"{name(i)}→{name(i + 1)} [{REASONS[int(feats.tau_reason[0, i])]}]" for i in never_t[:30]) or "none")
        + (" …" if len(never_t) > 30 else ""))
    lines += ["", f"Features file: `{doc['export']['npz']}` (SHA-256 {doc['export']['npz_sha256'][:16]}…).", ""]
    return "\n".join(lines)


def _sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()
