"""ABL-MYR-1 runner: Myr vs noMyr kinase-domain ensembles from the public de Buhr & Gräter Dryad set.

    PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1/run_abl_myr.py freeze-code
    PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1/run_abl_myr.py run --adapter-python <python with MDAnalysis>
    PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1/run_abl_myr.py report

Implements run_spec.json (frozen, SHA-256 in results/ABL-MYR-1_freeze.json) exactly: gates G0–G6, then the comparison.
A failed gate stops the run with the gate, the observed cause and nothing else computed. Engineering code only; the
synthetic tests are tests/test_abl_myr_runner*.py. The orchestrator needs the TopoFold core (κ/τ); the worker
subcommand (`_worker`) runs in the MDAnalysis interpreter, because no one environment has both.

Interpretations fixed here, before any study file is opened (the spec text is ambiguous on them):
* Archive members are located by their path suffix (e.g. ``PCA/traj_kinked_Myr_11.xtc``, ``Set1/system.gro``); exactly one
  match must exist, otherwise G1 fails. Members are extracted with validated relative paths only.
* G2 "atom names are backbone": allowed names are BACKBONE_NAMES below (heavy backbone atoms, terminal oxygens and
  hydrogens bonded to backbone atoms). Any other name fails G2.
* G4 "topol.top differ only by the Myr molecule": lines mentioning MYR/GLYM are removed from both files; the remainder
  (including [ molecules ] counts of every other molecule) must be identical. "Protein atoms" of system.gro are residues
  named as standard amino acids or CHARMM protonation variants.
* G5: non-increasing time stamps or an interval outside 2 ns ± 0.01 ns in the analysis frames stops the run; only a run
  with fewer than 450 analysis frames is excluded.
* G6: ≤ 5 % excluded windows per channel and set. If the first pass (pbc 'none') exceeds it and every excluded-window
  cause of the κ/τ/distance channels is a Cα–Cα break, one second pass is run for those channels with the adapter's
  only make-whole policy, ``make_whole_sequential`` (the spec's "make_whole").
  φ/ψ always use stored coordinates (a residue split across the box gives invalid φ/ψ, counted as excluded).
* Outcome step 2 ("inconclusive") also applies if a set is not analysed, because reproduction needs both sets.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import re
import resource
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path, PurePosixPath

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs.compare import circular_hist, js_distance  # noqa: E402
from topofold_inputs.ensemble import analyze, extract_member  # noqa: E402
from topofold_inputs.kinase_roles import Numbering, Reference, RoleMappingError, map_kinase_roles  # noqa: E402

SPEC = HERE / "run_spec.json"
FREEZE = HERE / "results/ABL-MYR-1_freeze.json"
CODE_FREEZE = HERE / "results/ABL-MYR-1_code_freeze.json"
SELF = Path(__file__).resolve()
LIBRARY = [f"python/topofold_inputs/{m}.py" for m in
           ("ensemble", "features", "compare", "ca_adapter", "atom_adapter", "kinase_geometry", "kinase_roles",
            "__main__")]
TESTS = ["tests/test_abl_myr_runner.py", "tests/test_abl_myr_runner_mda.py"]

CFG = {
    "conditions": {"Myr": "kinase_kinked_Myr", "noMyr": "kinase_kinked_noMyr"},
    "sets": {"Set1": list(range(11, 21)), "Set2": list(range(21, 31))},
    "t_min_ps": 1_000_000.0, "dt_ps": 2000.0, "dt_tol_ps": 10.0,
    "min_frames": 450, "min_runs": 8, "max_excluded_fraction": 0.05,
    "gro_tol_nm": 0.001, "bins": 36, "max_lag": 50,
    "anchors": {"beta3_Lys": 290, "alphaC_Glu": 305, "HRD_H": 380, "DFG_Asp": 400},
    "reference": {"accession": "P00519", "isoform": "P00519-2",
                  "fasta": "benchmarks/audit/abl_prereg/reference/UniProt_P00519-2.fasta"},
}
CHANNELS_E = ("kappa", "tau", "d_i_i2", "d_i_i3")      # existing ensemble channels
CHANNELS_ALL = CHANNELS_E + ("rmsf", "phi", "psi")
TOPOFOLD = ("kappa", "tau")
BASELINES = ("d_i_i2", "d_i_i3", "rmsf")
SPANS = {"kappa": (1, 1), "d_i_i2": (1, 1), "tau": (1, 2), "d_i_i3": (1, 2), "rmsf": (0, 0), "phi": (0, 0),
         "psi": (0, 0)}
BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT", "OT1", "OT2", "OC1", "OC2", "O1", "O2", "H", "HN", "HA", "HA1", "HA2",
                  "H1", "H2", "H3", "HT1", "HT2", "HT3"}
AMINO = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
         "THR", "TRP", "TYR", "VAL", "HSD", "HSE", "HSP", "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN"}
MYR_RE = re.compile(r"MYR|GLYM", re.I)
INIT_KEYS = ("integrator", "tcoupl", "ref_t", "tau_t", "pcoupl", "dt", "nsteps", "gen_vel", "gen_seed",
             "continuation", "constraints", "nstxout_compressed", "nstxtcout", "compressed_x_grps")


class GateError(Exception):
    def __init__(self, gate: str, message: str, detail=None):
        self.gate, self.message, self.detail = gate, message, detail
        super().__init__(f"{gate}: {message}")


def sha256_file(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cpu_seconds() -> float:
    a, b = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


# --- G0 / G1: archives and members -------------------------------------------------------------------------------------

def check_g0(data_dir: Path, spec: dict) -> dict:
    rec = {}
    for name in ("PCA.zip", "kinase_kinked_Myr.zip", "kinase_kinked_noMyr.zip"):
        want, p = spec["inputs"][name], Path(data_dir) / name
        if not p.is_file():
            raise GateError("G0", f"{name} is not present in the data directory", {"archive": name})
        size = p.stat().st_size
        if size != want["bytes"]:
            raise GateError("G0", f"{name} has {size} bytes, expected {want['bytes']}", {"archive": name})
        digest = sha256_file(p)
        if digest != want["sha256"]:
            raise GateError("G0", f"{name} SHA-256 {digest} differs from the frozen {want['sha256']}", {"archive": name})
        rec[name] = {"bytes": size, "sha256": digest}
    return rec


def resolve_member(names, suffix: str) -> str:
    hits = [n for n in names if (n == suffix or n.endswith("/" + suffix)) and not n.endswith("/")
            and not n.startswith("__MACOSX/")]
    if len(hits) != 1:
        raise GateError("G1", f"member {suffix!r}: {len(hits)} matches in the archive (exactly one required)",
                        {"suffix": suffix, "hits": hits[:5]})
    return hits[0]


def wanted_members(cfg: dict) -> dict:
    """{archive: {key: suffix}} for the frozen member list."""
    pca = {"ref": "PCA/KD_noIa_bb_ref.pdb", "avg": "PCA/average.pdb", "steps": "PCA/steps"}
    models = {}
    for cond, model in cfg["conditions"].items():
        m = {}
        for s, reps in cfg["sets"].items():
            for k, f in (("gro", "system.gro"), ("gro_em", "system_em.gro"), ("top", "topol.top"),
                         ("mdp_prod", "production_md.mdp"), ("mdp_nvt", "nvt.mdp"), ("mdp_npt", "npt.mdp")):
                m[f"{k}:{s}"] = f"{s}/{f}"
            for r in reps:
                m[f"tpr:{s}:{r}"] = f"{s}/{r}/md.tpr"
                pca[f"xtc:{s}:{cond}:{r}"] = f"PCA/traj_kinked_{cond}_{r}.xtc"
        models[f"{model}.zip"] = m
    return {"PCA.zip": pca, **models}


def select_members(archive_names: dict, cfg: dict) -> dict:
    out = {}
    for arch, keys in wanted_members(cfg).items():
        out[arch] = {k: resolve_member(archive_names[arch], suf) for k, suf in keys.items()}
    return out


def safe_extract(zf: zipfile.ZipFile, name: str, dest_root: Path) -> dict:
    parts = PurePosixPath(name).parts
    if not parts or name.startswith("/") or ".." in parts or "\\" in name:
        raise GateError("G1", f"unsafe member path {name!r}")
    dest = Path(dest_root).joinpath(*parts)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with zf.open(name) as src, open(dest, "wb") as out:      # CRC is verified when the member is read to the end
        shutil.copyfileobj(src, out, 1 << 20)
    return {"member": name, "path": str(dest), "bytes": dest.stat().st_size, "sha256": sha256_file(dest)}


# --- G4: shared start ---------------------------------------------------------------------------------------------------

def parse_gro(path) -> list:
    lines = Path(path).read_text().splitlines()
    n = int(lines[1])
    out = []
    for ln in lines[2:2 + n]:
        out.append((ln[5:10].strip(), ln[10:15].strip(), (float(ln[20:28]), float(ln[28:36]), float(ln[36:44]))))
    if len(out) != n:
        raise GateError("G4", f"{path}: {len(out)} atom lines, header says {n}")
    return out


def protein_atoms(atoms) -> list:
    return [a for a in atoms if a[0] in AMINO]


def compare_gro(a_path, b_path, tol_nm: float) -> dict:
    a, b = protein_atoms(parse_gro(a_path)), protein_atoms(parse_gro(b_path))
    if len(a) != len(b):
        return {"ok": False, "cause": f"protein atom counts differ ({len(a)} vs {len(b)})"}
    if [x[:2] for x in a] != [x[:2] for x in b]:
        i = next(i for i, (x, y) in enumerate(zip(a, b)) if x[:2] != y[:2])
        return {"ok": False, "cause": f"protein residue/atom names differ at atom {i + 1}: {a[i][:2]} vs {b[i][:2]}"}
    # rounded to 1e-9 nm so that a one-digit .gro difference equal to the tolerance is not failed by subtraction noise
    d = round(float(np.max(np.abs(np.array([x[2] for x in a]) - np.array([y[2] for y in b])))), 9) if a else 0.0
    return {"ok": d <= tol_nm, "max_abs_dx_nm": d, "n_protein_atoms": len(a),
            "cause": None if d <= tol_nm else f"max |Δx| {d:.4f} nm > {tol_nm} nm"}


def normalise_lines(text: str) -> list:
    out = []
    for ln in text.splitlines():
        ln = ln.split(";", 1)[0]
        ln = re.sub(r"\s+", "", ln)
        if ln:
            out.append(ln)
    return out


def mdp_values(text: str) -> dict:
    d = {}
    for ln in text.splitlines():
        ln = ln.split(";", 1)[0].strip()
        if "=" in ln:
            k, v = ln.split("=", 1)
            d[k.strip().lower().replace("-", "_")] = v.strip()
    return d


def compare_top(myr_text: str, nomyr_text: str) -> dict:
    def squeeze(t):
        return [re.sub(r"\s+", " ", ln.split(";", 1)[0]).strip() for ln in t.splitlines()
                if ln.split(";", 1)[0].strip() and not MYR_RE.search(ln.split(";", 1)[0])]
    a, b = squeeze(myr_text), squeeze(nomyr_text)
    if a == b:
        return {"ok": True, "differing_lines": []}
    diff = [f"Myr:{x} | noMyr:{y}" for x, y in zip(a, b) if x != y][:10]
    if len(a) != len(b):
        diff.append(f"line counts after removing Myr lines: {len(a)} vs {len(b)}")
    return {"ok": False, "differing_lines": diff}


def check_g4_set(files: dict, cfg: dict) -> dict:
    """files: {cond: {'gro','top','mdp_prod'} -> path} for one set. Raises GateError on any mismatch."""
    a, b = files["Myr"], files["noMyr"]
    gro = compare_gro(a["gro"], b["gro"], cfg["gro_tol_nm"])
    if not gro["ok"]:
        raise GateError("G4", f"system.gro protein atoms differ: {gro['cause']}", gro)
    ma, mb = Path(a["mdp_prod"]).read_text(), Path(b["mdp_prod"]).read_text()
    if normalise_lines(ma) != normalise_lines(mb):
        la, lb = normalise_lines(ma), normalise_lines(mb)
        d = [f"Myr:{x} | noMyr:{y}" for x, y in zip(la, lb) if x != y][:10]
        raise GateError("G4", "production_md.mdp differs between Myr and noMyr", {"differing_lines": d, "n_lines": [len(la), len(lb)]})
    top = compare_top(Path(a["top"]).read_text(), Path(b["top"]).read_text())
    if not top["ok"]:
        raise GateError("G4", "topol.top differs beyond the Myr molecule and its parameters", top)
    vals = mdp_values(ma)
    return {"gro": gro, "mdp_production_identical": True, "topology": top,
            "initialisation_recorded": {k: vals.get(k) for k in INIT_KEYS}}


# --- G2 / G5 -----------------------------------------------------------------------------------------------------------

def choose_topology(n_ref: int, n_avg: int, n_xtc: int) -> str:
    if n_ref == n_xtc:
        return "ref"
    if n_avg == n_xtc:
        return "avg"
    raise GateError("G2", f"trajectory atom count {n_xtc} matches neither KD_noIa_bb_ref.pdb ({n_ref}) nor average.pdb "
                          f"({n_avg})")


def check_g2(topo: dict, xtc_facts: dict) -> dict:
    counts = {k: v["n_atoms"] for k, v in xtc_facts.items()}
    if len(set(counts.values())) != 1:
        raise GateError("G2", "trajectories have different atom counts", {"counts": sorted(set(counts.values()))})
    n = next(iter(counts.values()))
    if topo["n_atoms"] != n:
        raise GateError("G2", f"topology has {topo['n_atoms']} atoms, trajectories {n}")
    bad = sorted(set(topo["names"]) - BACKBONE_NAMES)
    if bad:
        raise GateError("G2", f"atom names are not backbone: {bad}")
    missing = [r for r in topo["residues"] if not all(r[4:7])]
    return {"n_atoms": n, "n_trajectories": len(counts), "n_residues": len(topo["residues"]),
            "phi_psi_comparator_available": not missing, "residues_missing_N_CA_C": [r[:3] for r in missing][:10]}


def frame_gate(times_ps, cfg: dict) -> dict:
    t = np.asarray(times_ps, dtype=np.float64)
    if len(t) < 2 or not np.all(np.diff(t) > 0):
        raise GateError("G5", "time stamps are not strictly increasing", {"n_frames": int(len(t))})
    ok = t >= cfg["t_min_ps"] - 1e-3
    start = int(np.argmax(ok)) if ok.any() else len(t)
    n = len(t) - start
    if n >= 2:
        dts = np.diff(t[start:])
        if np.any(np.abs(dts - cfg["dt_ps"]) > cfg["dt_tol_ps"]):
            raise GateError("G5", "analysis-frame interval is outside 2 ns ± 0.01 ns",
                            {"min_dt_ps": float(dts.min()), "max_dt_ps": float(dts.max())})
    excluded = n < cfg["min_frames"]
    return {"start_index": start, "n_frames_file": int(len(t)), "n_analysis": int(n),
            "first_time_ns": float(t[0] / 1000.0), "last_time_ns": float(t[-1] / 1000.0), "excluded": bool(excluded),
            "cause": f"{n} analysis frames < {cfg['min_frames']}" if excluded else None}


# --- statistics --------------------------------------------------------------------------------------------------------

def e_stat(a, b, dist):
    """E = median cross-condition pair distance - max(median within-a, median within-b); None if a group has < 2 runs."""
    wa = [dist(a[i], a[j]) for i in range(len(a)) for j in range(i + 1, len(a))]
    wb = [dist(b[i], b[j]) for i in range(len(b)) for j in range(i + 1, len(b))]
    if not wa or not wb or not a or not b:
        return None
    cross = [dist(x, y) for x in a for y in b]
    return float(np.median(cross) - max(np.median(wa), np.median(wb)))


def rows_to_e(rows_between: list) -> dict:
    out = {}
    for r in rows_between:
        wa, wb = r["within_member_js_median_a"], r["within_member_js_median_b"]
        if wa is None or wb is None:
            continue
        out[(r["channel"], r["window"])] = float(r["between_member_js_median"] - max(wa, wb))
    return out


def scalar_e(values: dict, ga: list, gb: list) -> dict:
    """values: {run_id: (n,) array}; per-residue E with |difference| as the run-pair distance."""
    n = len(values[ga[0]])
    out = {}
    for i in range(n):
        xa, xb = [values[r][i] for r in ga], [values[r][i] for r in gb]
        if not np.all(np.isfinite(xa + xb)):
            continue
        e = e_stat(xa, xb, lambda p, q: abs(p - q))
        if e is not None:
            out[i] = e
    return out


def hist_e(hists: dict, ga: list, gb: list) -> dict:
    """hists: {run_id: {residue index: counts}}; JS distance between run histograms."""
    idx = set.intersection(*(set(hists[r]) for r in ga + gb))
    return {i: e_stat([hists[r][i] for r in ga], [hists[r][i] for r in gb], js_distance) for i in sorted(idx)}


def window_residues(channel: str, w: int) -> range:
    lo, hi = SPANS[channel]
    return range(w - lo, w + hi + 1)


def classify_outcome(reproduced: dict, n_residues: int, sets_analysed: bool = True) -> dict:
    """reproduced: {channel: set of window indices reproduced in both sets}. Ordered steps of the frozen spec."""
    def res(chs):
        s = set()
        for c in chs:
            for w in reproduced.get(c, ()):
                s.update(i for i in window_residues(c, w) if 0 <= i < n_residues)
        return s
    r_t, r_b = res(TOPOFOLD), res(BASELINES)
    n2 = set()
    for i in r_b:
        n2.update(range(i - 2, i + 3))
    only_t, only_b = sorted(r_t - n2), sorted(r_b - {j for k in r_t for j in range(k - 2, k + 3)})
    per = {c: res((c,)) for c in CHANNELS_ALL}
    jac = {f"{a}|{b}": (len(per[a] & per[b]) / len(per[a] | per[b]) if per[a] | per[b] else None)
           for i, a in enumerate(CHANNELS_ALL) for b in CHANNELS_ALL[i + 1:]}
    out = {"jaccard_reproduced_residues": jac, "R_T": sorted(r_t), "R_B": sorted(r_b), "R_T_outside_N2_R_B": only_t, "R_B_outside_N2_R_T": only_b,
           "R_phi_psi": sorted(res(("phi", "psi"))),
           "reproduced_window_counts": {c: len(reproduced.get(c, ())) for c in CHANNELS_ALL}}
    if not sets_analysed or not any(reproduced.get(c) for c in CHANNELS_ALL):
        out["label"] = "inconclusive"
        out["note"] = "a set was not analysed" if not sets_analysed else "no channel has a reproduced window"
    elif len(only_t) >= 3:
        out["label"] = "TopoFold adds reproduced residues"
    elif r_t <= n2:
        out["label"] = "baselines cover TopoFold"
    else:
        out["label"] = "partial"
    return out


def g6_fractions(rows: dict, extra_valid: dict, extra_excluded: dict) -> dict:
    """Excluded-window fraction per channel from analyze() rows plus φ/ψ counts."""
    valid = {c: set() for c in CHANNELS_E}
    for r in rows["distribution"]:
        valid[r["channel"]].add(r["window"])
    excl = {c: set() for c in CHANNELS_E}
    for r in rows["exclusions"]:
        excl[r["channel"]].add(r["window"])
    out = {}
    for c in CHANNELS_E:
        n = len(valid[c]) + len(excl[c])
        out[c] = {"windows": n, "excluded": len(excl[c]), "fraction": len(excl[c]) / n if n else 0.0}
    for c, n in extra_valid.items():
        e = extra_excluded.get(c, 0)
        out[c] = {"windows": n + e, "excluded": e, "fraction": e / (n + e) if n + e else 0.0}
    return out


def temporal_summary(rows_temporal: list) -> dict:
    out = {}
    for c in CHANNELS_E:
        rs = [r for r in rows_temporal if r["channel"] == c]
        lag1 = [r["lag1"] for r in rs if r["lag1"] is not None]
        dec = [r["decorrelation_lag"] for r in rs if r["decorrelation_lag"] is not None]
        out[c] = {"series": len(rs), "lag1_median": float(np.median(lag1)) if lag1 else None,
                  "lag1_max": float(np.max(lag1)) if lag1 else None,
                  "fraction_decorrelated_within_max_lag": len(dec) / len(rs) if rs else None,
                  "decorrelation_lag_frames_median": float(np.median(dec)) if dec else None}
    return out


# --- worker (MDAnalysis interpreter) -----------------------------------------------------------------------------------

def _json_default(o):
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))


def _w_xtcs(a: dict) -> dict:
    import warnings

    from MDAnalysis.coordinates.XTC import XTCReader
    out = {}
    for key, path in a["files"].items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r = XTCReader(path)
            times = [float(ts.time) for ts in r]
            out[key] = {"n_atoms": int(r.n_atoms), "n_frames": len(times), "times_ps": times}
            r.close()
    return out


def _w_count(a: dict) -> dict:
    import warnings

    from topofold_inputs.ca_adapter import load_universe
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return {"n_atoms": int(load_universe(a["topology"]).atoms.n_atoms)}


def _w_topo(a: dict) -> dict:
    import warnings

    from topofold_inputs.ca_adapter import _chain_of, _chain_source, load_universe
    from topofold_inputs.kinase_roles import observed_chain
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        u = load_universe(a["topology"])
        sel = u.select_atoms("protein")
        src = _chain_source(u, sel, "auto")
        chains = sorted({_chain_of(r.atoms & sel, src) for r in sel.residues})
        out = {"n_atoms": int(u.atoms.n_atoms), "names": sorted(set(str(x) for x in u.atoms.names)),
               "chains": chains, "n_protein_atoms": int(sel.n_atoms)}
        if len(chains) != 1:
            out["residues"] = []
            return out
        chain = chains[0]
        obs = observed_chain(u, chain)
        names_of = {(int(r.resid), str(getattr(r, "icode", "") or "")): set(str(x) for x in (r.atoms & sel).names)
                    for r in sel.residues}
        res = []
        for o in obs:
            nm = names_of.get((o.resid, o.icode), set())      # a modified residue outside 'protein': flagged, not a crash
            res.append([o.resid, o.icode, o.resname, o.bonded_to_previous, "N" in nm, "CA" in nm, "C" in nm])
        out["residues"], out["chain"] = res, chain
        ref = Reference.from_fasta(Path(a["fasta"]).read_text(), a["accession"], a["isoform"], "UniProt")
        try:
            m = map_kinase_roles(obs, chain, ref, Numbering("reference", 0, "isoform 1b numbering, declared"),
                                 a["anchors"])
            out["mapping"] = {"accepted": True, "report": m.report}
        except RoleMappingError as e:
            out["mapping"] = {"accepted": False, "problems": e.problems, "report": e.report}
    return out


def _w_phipsi(a: dict) -> dict:
    import warnings

    from topofold_inputs.atom_adapter import extract_atoms
    from topofold_inputs.ca_adapter import ResidueSpec, load_universe
    from topofold_inputs.kinase_geometry import dihedral
    res, chain = a["residues"], a["chain"]
    n = len(res)
    bonded = [bool(r[3]) if r[3] is not None else False for r in res]
    arrays, stats = {}, {}
    for run_id, spec in a["runs"].items():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = load_universe(a["topology"], spec["xtc"])
            req = [(ResidueSpec(chain, r[0], r[1]), ["N", "CA", "C"]) for r in res]
            tr = extract_atoms(u, run_id, req, on_missing="mark", altloc="reject", start=spec["start"])
        X = tr.coordinates.reshape(len(tr.frame_index), n, 3, 3)      # (F, residue, N/CA/C, xyz)
        N_, CA, C = X[:, :, 0], X[:, :, 1], X[:, :, 2]
        phi = np.full((len(X), n), np.nan)
        psi = np.full((len(X), n), np.nan)
        lim = a["bond_max"]

        def ok(p, q):
            return np.linalg.norm(p - q, axis=-1) <= lim
        for i in range(1, n):
            if bonded[i]:
                good = ok(C[:, i - 1], N_[:, i]) & ok(N_[:, i], CA[:, i]) & ok(CA[:, i], C[:, i])
                v = dihedral(C[:, i - 1], N_[:, i], CA[:, i], C[:, i])
                phi[:, i] = np.where(good, v, np.nan)
        for i in range(n - 1):
            if bonded[i + 1]:
                good = ok(N_[:, i], CA[:, i]) & ok(CA[:, i], C[:, i]) & ok(C[:, i], N_[:, i + 1])
                v = dihedral(N_[:, i], CA[:, i], C[:, i], N_[:, i + 1])
                psi[:, i] = np.where(good, v, np.nan)
        arrays[f"{run_id}__phi"], arrays[f"{run_id}__psi"] = phi, psi
        stats[run_id] = {"frames": int(len(X)), "frame_index_first": int(tr.frame_index[0])}
    np.savez(a["npz"], **arrays)
    return {"stats": stats}


def _w_rmsf(a: dict) -> dict:
    from MDAnalysis.analysis.align import rotation_matrix
    src = np.load(a["npz_in"])
    out = {}
    for run_id in src.files:
        X = np.asarray(src[run_id], dtype=np.float64)
        ok = np.isfinite(X).all(axis=(0, 2))
        Y = X[:, ok]
        ref = Y[0] - Y[0].mean(axis=0)
        al = np.empty_like(Y)
        for f in range(len(Y)):
            m = Y[f] - Y[f].mean(axis=0)
            R, _ = rotation_matrix(m, ref)
            al[f] = m @ R.T
        r = np.sqrt(((al - al.mean(axis=0)) ** 2).sum(axis=-1).mean(axis=0))
        full = np.full(X.shape[1], np.nan)
        full[ok] = r
        out[run_id] = full
    np.savez(a["npz_out"], **out)
    return {"runs": len(out)}


WORKERS = {"xtcs": _w_xtcs, "count": _w_count, "topo": _w_topo, "phipsi": _w_phipsi, "rmsf": _w_rmsf}


def call_worker(adapter_python: str, cmd: str, payload: dict, work: Path) -> dict:
    work.mkdir(parents=True, exist_ok=True)
    k = len(list(work.glob(f"worker_{cmd}_*.json")))
    pj, oj = work / f"worker_{cmd}_{k}.json", work / f"worker_{cmd}_{k}.out.json"
    pj.write_text(json.dumps({**payload, "out": str(oj)}, default=_json_default))
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(ROOT / "python"), os.environ.get("PYTHONPATH", "")]))
    p = subprocess.run([adapter_python, str(SELF), "_worker", cmd, str(pj)], env=env, capture_output=True, text=True)
    if p.returncode != 0:
        raise GateError("worker", f"worker {cmd} failed: {p.stderr[-1500:]}")
    return json.loads(oj.read_text())


# --- pipeline ----------------------------------------------------------------------------------------------------------

def run_pipeline(spec: dict, cfg: dict, data_dir: Path, out_dir: Path, adapter_python: str, work: Path | None = None,
                 log=print) -> dict:
    """All gates and the comparison. Returns the results dict; raises GateError on a failed gate."""
    t0, cpu0 = time.time(), cpu_seconds()
    data_dir, out_dir = Path(data_dir), Path(out_dir)
    work = Path(work) if work else HERE / "work"
    budget = spec["budget"]
    res: dict = {"id": "ABL-MYR-1", "gates": {}}

    def budget_check(stage):
        used = cpu_seconds() - cpu0
        if used > budget["cpu_hours_max"] * 3600:
            raise GateError("budget", f"CPU budget exceeded after {stage}: {used / 3600:.2f} h")

    res["gates"]["G0"] = check_g0(data_dir, spec)
    log("G0 passed")
    # G1: members
    names = {}
    for arch in res["gates"]["G0"]:
        with zipfile.ZipFile(data_dir / arch) as zf:
            names[arch] = zf.namelist()
    selected = select_members(names, cfg)
    ext_root = data_dir / "extracted"
    extracted = {}
    for arch, keys in selected.items():
        with zipfile.ZipFile(data_dir / arch) as zf:
            for k, member in keys.items():
                extracted[(arch, k)] = safe_extract(zf, member, ext_root / Path(arch).stem)
    total = sum(e["bytes"] for e in extracted.values()) + sum(v["bytes"] for v in res["gates"]["G0"].values())
    if total > budget["storage_bytes_max"]:
        raise GateError("budget", f"storage {total} B exceeds {budget['storage_bytes_max']} B")
    res["gates"]["G1"] = {"members": len(extracted), "extracted_bytes": sum(e["bytes"] for e in extracted.values()),
                          "files": {f"{a}:{k}": {"member": e["member"], "bytes": e["bytes"], "sha256": e["sha256"]}
                                    for (a, k), e in extracted.items()}}
    log(f"G1 passed ({len(extracted)} members)")
    path = lambda arch, key: extracted[(arch, key)]["path"]  # noqa: E731
    pca = "PCA.zip"
    # G2: atom counts, topology choice, backbone names
    xtc_paths = {f"{s}:{c}:{r}": path(pca, f"xtc:{s}:{c}:{r}") for s, reps in cfg["sets"].items()
                 for c in cfg["conditions"] for r in reps}
    facts = call_worker(adapter_python, "xtcs", {"files": xtc_paths}, work)
    n_xtc = {v["n_atoms"] for v in facts.values()}
    if len(n_xtc) != 1:
        raise GateError("G2", "trajectories have different atom counts", {"counts": sorted(n_xtc)})
    n_ref = call_worker(adapter_python, "count", {"topology": path(pca, "ref")}, work)["n_atoms"]
    n_avg = call_worker(adapter_python, "count", {"topology": path(pca, "avg")}, work)["n_atoms"]
    which = choose_topology(n_ref, n_avg, n_xtc.pop())
    topology = path(pca, which)
    fasta = ROOT / cfg["reference"]["fasta"]
    topo = call_worker(adapter_python, "topo", {"topology": topology, "fasta": str(fasta), "anchors": cfg["anchors"],
                                                "accession": cfg["reference"]["accession"],
                                                "isoform": cfg["reference"]["isoform"]}, work)
    if len(topo["chains"]) != 1:
        raise GateError("G2", f"topology has {len(topo['chains'])} protein chains ({topo['chains']}); one required")
    res["gates"]["G2"] = {"topology": which, "topology_sha256": extracted[(pca, which)]["sha256"],
                          **check_g2(topo, facts)}
    log("G2 passed")
    # G3: mapping
    mp = topo["mapping"]
    if not mp["accepted"]:
        raise GateError("G3", "residue-role mapping refused: " + "; ".join(mp["problems"]), mp["report"])
    rep = mp["report"]
    res["gates"]["G3"] = {"roles": rep["roles"], "alignment": {k: rep["alignment"][k] for k in (
        "aligned_pairs", "identity", "reference_span", "mutations", "insertions", "unaligned_terminal",
        "unobserved_reference_ranges", "chain_breaks_before", "author_numbering_diagnostic")}}
    log("G3 passed")
    # G4: shared start, per set
    g4 = {}
    for s in cfg["sets"]:
        files = {c: {"gro": path(f"{m}.zip", f"gro:{s}"), "top": path(f"{m}.zip", f"top:{s}"),
                     "mdp_prod": path(f"{m}.zip", f"mdp_prod:{s}")} for c, m in cfg["conditions"].items()}
        g4[s] = check_g4_set(files, cfg)
        g4[s]["md_tpr_sha256"] = {f"{c}:{r}": extracted[(f"{m}.zip", f"tpr:{s}:{r}")]["sha256"]
                                  for c, m in cfg["conditions"].items() for r in cfg["sets"][s]}
    res["gates"]["G4"] = g4
    log("G4 passed")
    # G5: frames and run exclusions
    runs, g5, excluded_runs = {}, {}, []
    for s, reps in cfg["sets"].items():
        for c in cfg["conditions"]:
            for r in reps:
                key = f"{s}:{c}:{r}"
                fg = frame_gate(facts[key]["times_ps"], cfg)
                g5[key] = fg
                if fg["excluded"]:
                    excluded_runs.append({"run": key, "cause": fg["cause"]})
                else:
                    runs[key] = {"set": s, "cond": c, "rep": r, "xtc": xtc_paths[key], "start": fg["start_index"]}
    analysed = {s: all(sum(1 for k in runs if runs[k]["set"] == s and runs[k]["cond"] == c) >= cfg["min_runs"]
                       for c in cfg["conditions"]) for s in cfg["sets"]}
    res["gates"]["G5"] = {"runs": g5, "excluded_runs": excluded_runs, "sets_analysed": analysed}
    log(f"G5 passed ({len(runs)} runs kept, {len(excluded_runs)} excluded; sets analysed {analysed})")
    if not any(analysed.values()):
        res["outcome"] = classify_outcome({}, len(topo["residues"]), sets_analysed=False)
        res["status"] = "complete"
        return _finish(res, out_dir, t0, cpu0)
    budget_check("G5")

    # extraction with the pbc rule (G6)
    residues_ref = [(r[0], r[1], r[2]) for r in topo["residues"]]
    n_res = len(residues_ref)
    options = {"bins": cfg["bins"], "max_lag": cfg["max_lag"], "pbc": "none"}
    phipsi = call_worker(adapter_python, "phipsi", {
        "topology": topology, "chain": topo["chain"], "residues": topo["residues"], "bond_max": 2.0,
        "runs": {k: {"xtc": v["xtc"], "start": v["start"]} for k, v in runs.items()},
        "npz": str(work / "phipsi.npz")}, work) if res["gates"]["G2"]["phi_psi_comparator_available"] else None
    pp = np.load(work / "phipsi.npz") if phipsi else None

    def extract(pbc):
        members = {}
        for k, v in runs.items():
            member = {"type": "mdanalysis", "topology": topology, "trajectories": [v["xtc"]], "start": v["start"],
                      "selection": "protein", "altloc": "reject", "pbc": pbc}
            m = extract_member(member, k.replace(":", "_"), work / "extract" / pbc, str(adapter_python), options)
            if [(r["resid"], r["icode"], r["resname"]) for r in m["residues"]] != residues_ref:
                raise GateError("G3", f"residues of {k} differ from the mapped topology residues")
            members[k] = m
        return members

    def per_set(members, s):
        g = {c: [k for k in members if runs[k]["set"] == s and runs[k]["cond"] == c] for c in cfg["conditions"]}
        return g

    def phipsi_hists(g):
        """{channel: (hists {run: {residue: counts}}, n_excluded_residues, n_valid_residues)}; empty if unavailable."""
        if pp is None:
            return {}
        out, ids = {}, [k for c in g for k in g[c]]
        for ch in ("phi", "psi"):
            hs, nbad = {}, 0
            for i in range(n_res):
                cols = {k: pp[f"{k}__{ch}"][:, i] for k in ids}
                if any(not np.isfinite(v).all() for v in cols.values()):
                    nbad += 1
                    continue
                for k, v in cols.items():
                    hs.setdefault(k, {})[i] = circular_hist(v, cfg["bins"], -math.pi)
            out[ch] = (hs, nbad, n_res - nbad)
        return out

    set_names = [s for s in cfg["sets"] if analysed[s]]
    passes = []
    for pbc in ("none", "make_whole_sequential"):
        members = extract(pbc)
        budget_check(f"extraction ({pbc})")
        fr, causes, full_by_set, ph_by_set = {}, set(), {}, {}
        for s in set_names:
            g = per_set(members, s)
            full_by_set[s] = analyze({c: [members[k] for k in g[c]] for c in g}, options)
            ph_by_set[s] = phipsi_hists(g)
            ph = ph_by_set[s]
            fr[s] = g6_fractions(full_by_set[s]["rows"], {c: ph[c][2] for c in ph}, {c: ph[c][1] for c in ph})
            causes |= {r["cause"] for r in full_by_set[s]["rows"]["exclusions"]}
        worst = max((v["fraction"] for d in fr.values() for v in d.values()), default=0.0)
        passes.append({"pbc": pbc, "g6_fractions": fr, "worst_fraction": worst, "exclusion_causes": sorted(causes)})
        if worst <= cfg["max_excluded_fraction"]:
            break
        if not (pbc == "none" and causes and causes <= {"ca_ca_distance"}):
            raise GateError("G6", f"more than {cfg['max_excluded_fraction']:.0%} of windows are excluded (worst "
                                  f"{worst:.3f}, pbc {pbc})", {"passes": passes})
    res["gates"]["G6"] = {"passes": passes, "pbc_used": passes[-1]["pbc"]}
    log(f"G6 passed (pbc {passes[-1]['pbc']})")

    # RMSF (MDAnalysis QCP in the worker) from the extracted Cα coordinates
    np.savez(work / "ca_coords.npz", **{k.replace(":", "_"): m["coords"] for k, m in members.items()})
    call_worker(adapter_python, "rmsf", {"npz_in": str(work / "ca_coords.npz"), "npz_out": str(work / "rmsf.npz")}, work)
    rm = np.load(work / "rmsf.npz")
    rmsf = {k: rm[k.replace(":", "_")] for k in members}

    res["sets"], reproduced_by_set = {}, {}
    for s in set_names:
        g = per_set(members, s)
        ga, gb = g["Myr"], g["noMyr"]
        full, ph = full_by_set[s], ph_by_set[s]
        e_full = rows_to_e(full["rows"]["between"])
        e_full.update({("rmsf", i): e for i, e in scalar_e(rmsf, ga, gb).items()})
        for ch in ph:
            e_full.update({(ch, i): e for i, e in hist_e(ph[ch][0], ga, gb).items() if e is not None})
        loo = []
        for drop in ga + gb:
            sub = {"Myr": [k for k in ga if k != drop], "noMyr": [k for k in gb if k != drop]}
            e = rows_to_e(analyze({c: [members[k] for k in sub[c]] for c in sub}, options)["rows"]["between"])
            e.update({("rmsf", i): v for i, v in scalar_e(rmsf, sub["Myr"], sub["noMyr"]).items()})
            for ch in ph:
                e.update({(ch, i): v for i, v in hist_e(ph[ch][0], sub["Myr"], sub["noMyr"]).items() if v is not None})
            loo.append(e)
        exceed = {key: bool(e > 0 and all(key in lo and lo[key] > 0 for lo in loo)) for key, e in e_full.items()}
        reproduced_by_set[s] = {c: {w for (ch, w), v in exceed.items() if ch == c and v} for c in CHANNELS_ALL}
        per_channel = {}
        for c in CHANNELS_ALL:
            ev = [e for (ch, w), e in e_full.items() if ch == c]
            per_channel[c] = {"windows_evaluated": len(ev), "windows_E_positive": int(sum(e > 0 for e in ev)),
                              "windows_exceeding": int(sum(1 for (ch, w), v in exceed.items() if ch == c and v)),
                              "E_median": float(np.median(ev)) if ev else None,
                              "E_max": float(np.max(ev)) if ev else None}
        within = {}
        for c in CHANNELS_E:
            wa = [r["within_member_js_median"] for r in full["rows"]["replicas"] if r["channel"] == c and r["group"] == "Myr"
                  and r["within_member_js_median"] is not None]
            wb = [r["within_member_js_median"] for r in full["rows"]["replicas"] if r["channel"] == c and r["group"] == "noMyr"
                  and r["within_member_js_median"] is not None]
            bt = [r["between_member_js_median"] for r in full["rows"]["between"] if r["channel"] == c]
            within[c] = {"within_Myr_js_median_of_windows": float(np.median(wa)) if wa else None,
                         "within_noMyr_js_median_of_windows": float(np.median(wb)) if wb else None,
                         "between_js_median_of_windows": float(np.median(bt)) if bt else None}
        res["sets"][s] = {"runs": {c: g[c] for c in g}, "channels": per_channel, "run_variability": within,
                          "temporal": temporal_summary(full["rows"]["temporal"]),
                          "windows": [{"channel": ch, "window": w, "E": e, "E_loo_min": float(min(lo[(ch, w)] for lo in loo)),
                                       "exceeding": exceed[(ch, w)]} for (ch, w), e in sorted(e_full.items())],
                          "exclusion_rows": len(full["rows"]["exclusions"]),
                          "exclusion_causes": sorted({r["cause"] for r in full["rows"]["exclusions"]})}
    both = len(reproduced_by_set) == len(cfg["sets"])
    repro = {c: set.intersection(*(reproduced_by_set[s][c] for s in reproduced_by_set)) if both else set()
             for c in CHANNELS_ALL}
    res["outcome"] = classify_outcome(repro, n_res, sets_analysed=both)
    res["outcome"]["reproduced_windows"] = {c: sorted(v) for c, v in repro.items()}
    res["residues"] = [{"index": i, "resid": r[0], "icode": r[1], "resname": r[2]} for i, r in enumerate(residues_ref)]
    res["status"] = "complete"
    return _finish(res, out_dir, t0, cpu0)


def _finish(res, out_dir, t0, cpu0):
    res["resources"] = {"cpu_s": cpu_seconds() - cpu0, "wall_s": time.time() - t0}
    return res


# --- freeze, run, report -----------------------------------------------------------------------------------------------

def _need_freeze():
    frz = json.loads(FREEZE.read_text())
    if sha256_file(SPEC) != frz["spec_sha256"]:
        raise SystemExit("STOP: run_spec.json differs from the freeze record")
    return frz


def cmd_freeze_code(args) -> int:
    if CODE_FREEZE.exists():
        print("code freeze exists")
        return 1
    frz = _need_freeze()
    data = HERE / "data_dryad"
    opened = [str(d.relative_to(HERE)) for d in (data / "extracted", HERE / "work") if d.exists()]
    if opened:
        raise SystemExit(f"STOP: {opened} exist, so the runner has already opened study files; the code freeze must precede it")
    # Presence is not opening: names, sizes and mtimes only (no hashing, no reading of contents).
    present = [{"name": q.name, "bytes": q.stat().st_size,
                "mtime_utc": dt.datetime.fromtimestamp(q.stat().st_mtime, dt.timezone.utc).isoformat()}
               for q in sorted(data.glob("*"))] if data.exists() else []
    rec = {"id": "ABL-MYR-1", "code_frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "spec_sha256": frz["spec_sha256"], "spec_freeze_record_sha256": sha256_file(FREEZE),
           "runner": {"path": "benchmarks/audit/abl_myr_1/run_abl_myr.py", "sha256": sha256_file(SELF)},
           "tests_sha256": {t: sha256_file(ROOT / t) for t in TESTS},
           "library_sha256": {m: sha256_file(ROOT / m) for m in LIBRARY},
           "reference_fasta_sha256": sha256_file(ROOT / CFG["reference"]["fasta"]),
           "input_files_present_at_code_freeze": present,
           "input_files_note": "listed by name, size and mtime only; not hashed, read or opened by any TopoFold code",
           "frozen_cfg": CFG}
    CODE_FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps({k: rec[k] for k in ("code_frozen_utc", "runner", "tests_sha256")}, indent=1))
    return 0


def verify_code_freeze() -> dict:
    if not CODE_FREEZE.exists():
        raise SystemExit("STOP: no code freeze record; run freeze-code before opening any study archive")
    _need_freeze()
    cf = json.loads(CODE_FREEZE.read_text())
    bad = []
    if sha256_file(SELF) != cf["runner"]["sha256"]:
        bad.append("runner")
    bad += [m for m, h in cf["library_sha256"].items() if sha256_file(ROOT / m) != h]
    if sha256_file(ROOT / CFG["reference"]["fasta"]) != cf["reference_fasta_sha256"]:
        bad.append("reference fasta")
    if cf["frozen_cfg"] != json.loads(json.dumps(CFG)):
        bad.append("CFG")
    if bad:
        raise SystemExit(f"STOP: code differs from the code freeze: {bad}")
    return cf


def scrub(x):
    """Remove machine-specific absolute paths from text recorded in the repository."""
    if isinstance(x, str):
        return x.replace(str(ROOT), "<repo>").replace(str(Path.home()), "~")
    if isinstance(x, dict):
        return {k: scrub(v) for k, v in x.items()}
    if isinstance(x, list):
        return [scrub(v) for v in x]
    return x


def cmd_run(args) -> int:
    cf = verify_code_freeze()
    spec = json.loads(SPEC.read_text())
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    data_dir = Path(args.data_dir)
    try:
        res = run_pipeline(spec, CFG, data_dir, out_dir, args.adapter_python, work=Path(args.work) if args.work else None)
    except GateError as e:
        blocked = {"id": "ABL-MYR-1", "status": "blocked", "gate": e.gate, "cause": scrub(e.message), "detail": scrub(e.detail),
                   "utc": dt.datetime.now(dt.timezone.utc).isoformat(), "spec_sha256": cf["spec_sha256"],
                   "runner_sha256": cf["runner"]["sha256"]}
        (out_dir / "ABL-MYR-1_blocked.json").write_text(json.dumps(blocked, indent=1, default=_json_default) + "\n")
        print(f"BLOCKED at {e.gate}: {e.message}")
        return 3
    res["spec_sha256"], res["runner_sha256"] = cf["spec_sha256"], cf["runner"]["sha256"]
    (out_dir / "ABL-MYR-1_results.json").write_text(json.dumps(res, indent=1, default=_json_default) + "\n")
    for s, d in res.get("sets", {}).items():
        with open(out_dir / f"ABL-MYR-1_windows_{s}.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["channel", "window", "E", "E_loo_min", "exceeding"], lineterminator="\n")
            w.writeheader()
            w.writerows(d["windows"])
    print(json.dumps(res["outcome"], default=_json_default))
    return 0


def render_report(res: dict) -> str:
    L = ["# ABL-MYR-1 results (generated from ABL-MYR-1_results.json)", "",
         f"Status: **{res['status']}**. Outcome: **{res['outcome']['label']}**"
         + (f" ({res['outcome']['note']})" if res["outcome"].get("note") else ""), "",
         "Descriptive evidence from simulated Abl ensembles; separate from Ayaz WT/M472I; not biological validation.", ""]
    for s, d in res.get("sets", {}).items():
        L += [f"## {s}", "", "| channel | windows | E > 0 | exceeding (full and all leave-one-run-out) | median E |", "|---|---|---|---|---|"]
        for c, v in d["channels"].items():
            m = "" if v["E_median"] is None else f"{v['E_median']:.4f}"
            L.append(f"| {c} | {v['windows_evaluated']} | {v['windows_E_positive']} | {v['windows_exceeding']} | {m} |")
        L += ["", "Run variability (median over windows of Jensen–Shannon distances):", ""]
        for c, v in d["run_variability"].items():
            L.append(f"- {c}: within Myr {v['within_Myr_js_median_of_windows']}, within noMyr "
                     f"{v['within_noMyr_js_median_of_windows']}, between {v['between_js_median_of_windows']}")
        L += ["", "Temporal diagnostics:", ""]
        for c, v in d["temporal"].items():
            L.append(f"- {c}: lag-1 median {v['lag1_median']}, max {v['lag1_max']}; decorrelated within the lag limit in "
                     f"{v['fraction_decorrelated_within_max_lag']}")
        L.append("")
    o = res["outcome"]
    L += ["## Outcome", "", f"- R_T (TopoFold residues): {o['R_T']}", f"- R_B (baseline residues): {o['R_B']}",
          f"- TopoFold residues outside N2(R_B): {o['R_T_outside_N2_R_B']}",
          f"- baseline residues outside N2(R_T): {o['R_B_outside_N2_R_T']}", f"- φ/ψ comparator residues: {o['R_phi_psi']}", ""]
    return "\n".join(L) + "\n"


def cmd_report(args) -> int:
    p = Path(args.out) / "ABL-MYR-1_results.json"
    (Path(args.out) / "ABL-MYR-1_report.md").write_text(render_report(json.loads(p.read_text())))
    return 0


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "_worker":
        a = json.loads(Path(argv[2]).read_text())
        out = WORKERS[argv[1]](a)
        Path(a["out"]).write_text(json.dumps(out, default=_json_default))
        return 0
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("freeze-code")
    r = sub.add_parser("run")
    r.add_argument("--adapter-python", required=True, help="interpreter with MDAnalysis")
    r.add_argument("--data-dir", default=str(HERE / "data_dryad"))
    r.add_argument("--out", default=str(HERE / "results"))
    r.add_argument("--work", default=None)
    p = sub.add_parser("report")
    p.add_argument("--out", default=str(HERE / "results"))
    args = ap.parse_args(argv)
    return {"freeze-code": cmd_freeze_code, "run": cmd_run, "report": cmd_report}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
