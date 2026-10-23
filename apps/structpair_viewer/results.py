"""Read-only loader for `topofold_inputs.structpair` output directories (no Streamlit, no TopoFold core, no MDAnalysis).

Everything here only reads files. Values are taken from windows.csv, residues.csv, exclusions.csv, summary.json and
manifest.json exactly as written; empty CSV cells stay missing (NaN / None) and are never replaced by zero. The only
computation is for display: a rigid Cα superposition of the two input structures for the molecular view.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

SCHEMA = "topofold_inputs.structpair/1"
REQUIRED = ("manifest.json", "summary.json", "windows.csv", "residues.csv", "exclusions.csv")
# essential scientific outputs: each must have a valid SHA-256 entry in manifest["outputs_sha256"] and match it.
# Optional outputs (report.md, figure.png) are verified when listed; listed-but-absent ones are noted, not fatal.
# The manifest does not hash itself.
HASHED_REQUIRED = ("summary.json", "windows.csv", "residues.csv", "exclusions.csv")
METRICS = {  # key -> (windows.csv column, unit, window size, short label)
    "d_kappa": ("d_kappa_rad", "rad", 3, "|Δκ| (curvature)"),
    "d_i_i2": ("d_i_i2_angstrom", "Å", 3, "|Δd(i−1,i+1)| (distance baseline)"),
    "rmsd3": ("rmsd3_angstrom", "Å", 3, "local Cα RMSD, 3 residues"),
    "d_tau": ("d_tau_rad", "rad", 4, "|Δτ| (torsion, circular)"),
    "d_i_i3": ("d_i_i3_angstrom", "Å", 4, "|Δd(i−1,i+2)| (distance baseline)"),
    "rmsd4": ("rmsd4_angstrom", "Å", 4, "local Cα RMSD, 4 residues"),
}
WINDOW_COLS = ["position", "A_residue", "B_residue", "resname_A", "resname_B", "window3_A", "window4_A",
               "altloc_in_window3", "altloc_in_window4", "mismatch_in_window3", "mismatch_in_window4"]
CAUSE_TEXT = {
    "numbering_discontinuity": "The two adjacent mapped residues are not consecutive in author numbering (a gap in "
                               "the correspondence or unobserved residues), so no continuous window spans them.",
    "chain_boundary": "The window would cross from one chain to another.",
    "missing_CA": "A residue next to this junction has no Cα coordinates in this structure.",
    "missing_coordinates": "A residue in the window has no Cα coordinates in this structure (unobserved, or absent "
                           "from the file).",
    "ca_ca_distance": "Consecutive Cα atoms are farther apart than the adapter's break threshold (chain break or "
                      "coordinate error).",
    "terminal_or_break": "The window touches a chain end or break, but no adapter record names which.",
    "degenerate_segment": "Two consecutive Cα positions coincide, so curvature or torsion is undefined.",
    "near_collinear": "Three consecutive Cα atoms are nearly collinear, so the torsion is numerically undefined.",
    "run_too_short_for_core": "The continuous run of residues is too short for the core calculation.",
    "core_error": "The TopoFold core reported an error for this window.",
    "non_finite_core_output": "The core returned a non-finite value for this window.",
    "distance_not_finite": "A Cα–Cα distance in the window is not finite.",
}


class ResultError(ValueError):
    """The directory is not a readable, intact structpair result."""


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@dataclass
class StructpairResult:
    path: Path
    manifest: dict
    summary: dict
    windows: pd.DataFrame
    residues: pd.DataFrame
    exclusions: pd.DataFrame
    integrity: dict
    report_md: str | None = None
    figure: Path | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def name(self) -> str:
        return self.manifest.get("name") or self.path.name

    def label(self, key: str) -> str:
        return self.manifest["inputs"][key].get("label", key)


def _bool_col(s: pd.Series) -> pd.Series:
    return s.map({"True": True, "False": False, True: True, False: False})


def load_result(path: str | Path) -> StructpairResult:
    """Load one structpair output directory; raise ResultError with a plain message if it is not valid."""
    p = Path(path)
    if not p.is_dir():
        raise ResultError(f"{p} is not a directory")
    missing = [f for f in REQUIRED if not (p / f).is_file()]
    if missing:
        raise ResultError(f"{p} is not a structpair result: missing {', '.join(missing)}")
    try:
        manifest = json.loads((p / "manifest.json").read_text())
        summary = json.loads((p / "summary.json").read_text())
    except json.JSONDecodeError as e:
        raise ResultError(f"{p}: unreadable JSON ({e})") from e
    if manifest.get("schema") != SCHEMA:
        raise ResultError(f"{p}: manifest schema {manifest.get('schema')!r}, expected {SCHEMA!r}")
    # integrity: essential outputs need a valid hash entry; every recorded output must match its hash
    recorded = manifest.get("outputs_sha256")
    if not isinstance(recorded, dict):
        raise ResultError(f"{p}: manifest has no outputs_sha256 table")
    no_entry = [f for f in HASHED_REQUIRED
                if not (isinstance(recorded.get(f), str) and re.fullmatch(r"[0-9a-f]{64}", recorded[f]))]
    if no_entry:
        raise ResultError(f"{p}: manifest lacks a valid SHA-256 entry for {', '.join(no_entry)}; "
                          "these outputs cannot be verified")
    checked = {f: (sha256(p / f) == h) if (p / f).is_file() else None for f, h in recorded.items()}
    bad = [f for f, ok in checked.items() if ok is False]
    absent = [f for f, ok in checked.items() if ok is None]
    if bad or any(f in absent for f in HASHED_REQUIRED):
        raise ResultError(f"{p}: files differ from the manifest hashes: {', '.join(bad + absent)}")
    # keep_default_na=False + explicit "" handling: empty cells -> NaN, never 0; residue labels stay strings
    read = lambda f: pd.read_csv(p / f, dtype=str, keep_default_na=False)  # noqa: E731
    w, r, x = read("windows.csv"), read("residues.csv"), read("exclusions.csv")
    need_w = WINDOW_COLS + [c for c, *_ in METRICS.values()]
    lost = [c for c in need_w if c not in w.columns]
    if lost:
        raise ResultError(f"{p}: windows.csv lacks columns {lost}")
    for c in w.columns:
        if c.endswith(("_rad", "_angstrom")):
            w[c] = pd.to_numeric(w[c].replace("", np.nan), errors="raise").astype(float)
    for c in ("altloc_in_window3", "altloc_in_window4", "mismatch_in_window3", "mismatch_in_window4"):
        w[c] = _bool_col(w[c])
    w["position"] = w["position"].astype(int)
    r["position"] = r["position"].astype(int)
    r["resname_mismatch"] = _bool_col(r["resname_mismatch"])
    if len(r) != len(w) or not (w["position"].to_numpy() == r["position"].to_numpy()).all():
        raise ResultError(f"{p}: windows.csv and residues.csv do not share one position index")
    if len(x):
        x["position"] = x["position"].astype(int)
        x["window_size"] = x["window_size"].astype(int)
    x["explanation"] = [CAUSE_TEXT.get(c, "No standard explanation; see the detail column.") for c in x["cause"]]
    notes = []
    if absent:
        notes.append(f"Optional outputs listed in the manifest but absent: {', '.join(absent)}")
    return StructpairResult(path=p, manifest=manifest, summary=summary, windows=w, residues=r, exclusions=x,
                            integrity=checked, report_md=(p / "report.md").read_text() if (p / "report.md").is_file() else None,
                            figure=(p / "figure.png") if (p / "figure.png").is_file() else None, notes=notes)


def discover(root: str | Path) -> dict[str, Path]:
    """Structpair result directories directly under root (or root itself), keyed by directory name."""
    root = Path(root)
    if (root / "manifest.json").is_file():
        return {root.name: root}
    return {d.name: d for d in sorted(root.iterdir()) if (d / "manifest.json").is_file()} if root.is_dir() else {}


def snapshot(paths) -> dict[str, str]:
    """SHA-256 of every file under the given directories (used to show that viewing changes nothing)."""
    out = {}
    for d in paths:
        for f in sorted(Path(d).rglob("*")):
            if f.is_file():
                out[str(f)] = sha256(f)
    return out


# --- per-position views ---------------------------------------------------------------------------------------------

def metric_frame(res: StructpairResult, metric: str) -> pd.DataFrame:
    """One metric along the chain: value (NaN where unavailable), status and residue labels of both structures."""
    col, unit, size, _ = METRICS[metric]
    w = res.windows
    n = len(w)
    lo, hi = 1, n - (2 if size == 3 else 3)
    in_range = (w["position"] >= lo) & (w["position"] <= hi)
    excluded = set(res.exclusions.loc[res.exclusions["window_size"] == size, "position"]) if len(res.exclusions) else set()
    status = np.where(~in_range, "outside mapped range",
                      np.where(w[col].notna(), "value", np.where(w["position"].isin(excluded), "excluded",
                                                                 "unavailable")))
    win = w[f"window{size}_A"]
    df = pd.DataFrame({"position": w["position"], "A_residue": w["A_residue"], "B_residue": w["B_residue"],
                       "window_A": win, "value": w[col], "status": status,
                       "altloc_in_window": w[f"altloc_in_window{size}"],
                       "mismatch_in_window": w[f"mismatch_in_window{size}"]})
    # contiguous runs of available values, so a line is never drawn across a gap
    seg, k, prev = [], 0, False
    for ok in df["value"].notna():
        if ok and not prev:
            k += 1
        seg.append(k if ok else None)
        prev = ok
    df["segment"] = seg
    df.attrs.update(metric=metric, unit=unit, window_size=size, column=col)
    return df


def window_members(res: StructpairResult, position: int, size: int) -> list[int]:
    """Positions (rows of residues.csv) covered by the window indexed at `position` (i−1 .. i+size−2)."""
    n = len(res.residues)
    return [i for i in range(position - 1, position + size - 1) if 0 <= i < n]


def exclusion_explanations(res: StructpairResult, position: int, size: int) -> pd.DataFrame:
    x = res.exclusions
    if not len(x):
        return x
    return x[(x["position"] == position) & (x["window_size"] == size)]


# --- molecular view (display only) ---------------------------------------------------------------------------------

def _label(chain: str, resid: int, icode: str) -> str:
    return f"{chain}:{resid}{icode}"


def read_ca(pdb: Path, chain: str, altloc_policy: str, expected_sha256: str | None = None) -> dict[str, np.ndarray]:
    """Cα coordinates of one chain from a PDB file (first model), keyed by structpair residue label 'C:123A'.

    Alternate locations follow the policy recorded in the manifest (blank_then_A or highest_occupancy; 'reject'
    uses the single location, which structpair guaranteed). With expected_sha256, the bytes parsed are the bytes
    hashed, and a mismatch raises ResultError."""
    data = Path(pdb).read_bytes()
    if expected_sha256 is not None and hashlib.sha256(data).hexdigest() != expected_sha256:
        raise ResultError(f"{Path(pdb).name}: SHA-256 differs from the manifest; molecular view disabled for this "
                          "structure")
    cand: dict[str, list] = {}
    for line in data.decode().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")) or line[12:16].strip() != "CA" or line[21] != chain:
            continue
        key = _label(chain, int(line[22:26]), line[26].strip())
        occ = float(line[54:60]) if line[54:60].strip() else 1.0
        cand.setdefault(key, []).append((line[16], occ, np.array([float(line[30:38]), float(line[38:46]),
                                                                  float(line[46:54])])))
    out = {}
    for key, alts in cand.items():
        if altloc_policy == "highest_occupancy":
            pick = max(alts, key=lambda a: a[1])
        else:
            pick = next((a for a in alts if a[0] == " "), None) or next((a for a in alts if a[0] == "A"), alts[0])
        out[key] = pick[2]
    return out


def locate_input(res: StructpairResult, key: str, repo_root: Path) -> tuple[Path | None, str]:
    """Find the input structure recorded in the manifest and confirm its SHA-256; return (path, status message)."""
    inp = res.manifest["inputs"][key]
    rec = Path(inp["path"])
    cands = [rec]
    parts = rec.parts
    if "benchmarks" in parts:   # same repository checked out elsewhere
        cands.append(repo_root.joinpath(*parts[parts.index("benchmarks"):]))
    for c in cands:
        if c.is_file():
            if sha256(c) == inp["sha256"]:
                return c, f"{c.name}: SHA-256 matches the manifest"
            return None, f"{c}: SHA-256 differs from the manifest; molecular view disabled for this structure"
    return None, (f"input file not found ({rec.name}); run benchmarks/audit/structpair_demo/fetch_inputs.py or place "
                  "the file at the recorded path")


def display_coordinates(res: StructpairResult, coords: dict[str, dict[str, np.ndarray]]):
    """Per-position Cα arrays for A (superposed on B, display only) and B; NaN where a residue is missing.

    The superposition uses every position present in both structures (proper rotation, SVD with determinant
    correction). It is a display aid, not a new result."""
    r = res.residues
    A = np.array([coords["A"].get(lab, [np.nan] * 3) for lab in r["A_residue"]], float)
    B = np.array([coords["B"].get(lab, [np.nan] * 3) for lab in r["B_residue"]], float)
    A[r["A_status"].to_numpy() != "present"] = np.nan
    B[r["B_status"].to_numpy() != "present"] = np.nan
    ok = np.isfinite(A).all(1) & np.isfinite(B).all(1)
    if ok.sum() >= 3:
        cA, cB = A[ok].mean(0), B[ok].mean(0)
        U, _, Vt = np.linalg.svd((A[ok] - cA).T @ (B[ok] - cB))
        d = math.copysign(1.0, np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
        A = (A - cA) @ R.T + cB
    return A, B
