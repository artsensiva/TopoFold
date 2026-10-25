"""KRAS K1 step S1: the crystal-only endpoint (frozen implementation). No MD coordinate is read.

    .venv/bin/python benchmarks/audit/kras_k1/k1_endpoint.py freeze    # record hashes; verifies the input identities
    .venv/bin/python benchmarks/audit/kras_k1/k1_endpoint.py run       # compute the endpoint once, apply the gate

Implements section 4.1 and 4.2 of PROTOCOL.md with the parameters of protocol_params.json: per-residue Calpha displacement
between 5US4 (GDP) and 5USJ (GMPPNP) chains, the positive rule, the segments and the adequacy gate. `run` refuses to start
unless the protocol, the parameter file, this file and the tests match the freeze record, refuses to overwrite an existing
endpoint record, and computes the endpoint once. Distances are in Å, float64.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import itertools
import json
import math
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PARAMS = HERE / "protocol_params.json"
PROTOCOL = HERE / "PROTOCOL.md"
FREEZE = HERE / "results/KRAS-K1_freeze.json"
RECORD = HERE / "results/KRAS-K1_endpoint.json"
FEAS = ROOT / "benchmarks/audit/next_app_feasibility/results/feasibility_metadata.json"
PROV = HERE / "results/k1_provenance_checks.json"
SPEC = ROOT / "benchmarks/audit/kras_struct_g12d/run_spec.json"
SELF = Path(__file__).resolve()
TESTS = ["tests/test_kras_k1_endpoint.py", "tests/test_kras_k1_protocol.py"]
FROZEN_STRUCT2_SET = frozenset(range(1, 60)) | frozenset(range(68, 167))     # KRAS-STRUCT-2 mapping (158 residues)
TOL = 1e-9                                                                    # Å, inclusive comparisons
ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
       "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while b := f.read(1 << 20):
            h.update(b)
    return h.hexdigest()


# --- parsing ---------------------------------------------------------------------------------------------------------

def parse_ca(text: str) -> tuple:
    """({chain: {resid: xyz}}, {chain: {resid: resname}}) from the ATOM records' Calpha atoms.

    Raises on alternate locations, insertion codes, a duplicate Calpha or zero occupancy (none is expected in 5US4/5USJ);
    HETATM records are ignored.
    """
    ca, names = {}, {}
    for ln in text.splitlines():
        if not ln.startswith("ATOM  ") or ln[12:16].strip() != "CA":
            continue
        chain, resid = ln[21], int(ln[22:26])
        if ln[16] != " ":
            raise ValueError(f"alternate location at {chain}{resid}")
        if ln[26] != " ":
            raise ValueError(f"insertion code at {chain}{resid}")
        if float(ln[54:60]) <= 0.0:
            raise ValueError(f"zero occupancy at {chain}{resid}")
        if resid in ca.setdefault(chain, {}):
            raise ValueError(f"duplicate Calpha at {chain}{resid}")
        ca[chain][resid] = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])], dtype=np.float64)
        names.setdefault(chain, {})[resid] = ln[17:20].strip()
    return ca, names


# --- geometry --------------------------------------------------------------------------------------------------------

def kabsch(moving: np.ndarray, reference: np.ndarray) -> tuple:
    """Proper rotation R (det +1) and translation t with moving @ R.T + t closest to reference (least squares)."""
    cm, cr = moving.mean(axis=0), reference.mean(axis=0)
    u, _, vt = np.linalg.svd((moving - cm).T @ (reference - cr))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d if d != 0 else 1.0]) @ u.T
    return r, cr - r @ cm


def trim_indices(disp: np.ndarray, resids: list, m: int) -> list:
    """Indices of the m most displaced residues; equal displacements drop the lower residue number first."""
    order = sorted(range(len(resids)), key=lambda i: (-float(disp[i]), resids[i]))
    return order[:m]


def pair_displacement(reference: np.ndarray, moving: np.ndarray, resids: list, trim_fraction: float) -> np.ndarray:
    """Per-residue displacement (Å) of `moving` from `reference` after the trimmed fit of the protocol (4.1)."""
    n = len(resids)
    r, t = kabsch(moving, reference)
    first = np.linalg.norm(reference - (moving @ r.T + t), axis=1)
    drop = set(trim_indices(first, resids, int(math.floor(trim_fraction * n + 1e-9))))
    keep = [i for i in range(n) if i not in drop]
    r2, t2 = kabsch(moving[keep], reference[keep])
    return np.linalg.norm(reference - (moving @ r2.T + t2), axis=1)


# --- endpoint --------------------------------------------------------------------------------------------------------

def segments(residues, merge_gap: int) -> list:
    """Groups of positives: a residue joins the previous group if its number exceeds the previous one by at most gap + 1."""
    out = []
    for r in sorted(set(residues)):
        if out and r - out[-1][-1] - 1 <= merge_gap:
            out[-1].append(r)
        else:
            out.append([r])
    return out


def compute_positive(b: np.ndarray, w: np.ndarray, min_disp: float, factor: float) -> tuple:
    """(threshold, positive mask, boundary mask): positive if b >= max(min_disp, factor * w), inclusive within TOL."""
    thr = np.maximum(min_disp, factor * w)
    return thr, b >= thr - TOL, np.abs(b - thr) <= TOL


def gate(positives: list, n_scored: int, g: dict, merge_gap: int) -> dict:
    max_pos = int(math.floor(g["max_positive_fraction"] * n_scored + 1e-9))
    segs = segments(positives, merge_gap)
    long_segs = [s for s in segs if len(s) >= g["min_segment_length"]]
    checks = {"n_positive_at_least_min": len(positives) >= g["min_positive"],
              "n_positive_at_most_max": len(positives) <= max_pos,
              "segments_of_min_length_at_least_min": len(long_segs) >= g["min_segments"]}
    return {"n_positive": len(positives), "min_positive": g["min_positive"], "max_positive": max_pos,
            "segments": [{"residues": s, "length": len(s)} for s in segs], "segments_meeting_min_length": len(long_segs),
            "min_segments": g["min_segments"], "min_segment_length": g["min_segment_length"], "checks": checks,
            "passes": all(checks.values())}


def check_identity(ca: dict, names: dict, chains: list, topo_seq: str, expected_observed=FROZEN_STRUCT2_SET,
                   topo_range=(4, 164), expected_n=153) -> list:
    """Scored residues after the identity checks of the protocol; raises ValueError on any failed check."""
    if sorted(ca) != sorted(chains):
        raise ValueError(f"chains found {sorted(ca)} differ from the frozen chains {sorted(chains)}")
    observed = set.intersection(*(set(ca[c]) for c in chains))
    within = {r for r in observed if 1 <= r <= 166}
    if within != set(expected_observed):
        raise ValueError("the residues observed in all chains differ from the frozen KRAS-STRUCT-2 set")
    scored = sorted(r for r in observed if topo_range[0] <= r <= topo_range[1])
    if expected_n is not None and len(scored) != expected_n:
        raise ValueError(f"{len(scored)} scored residues, expected {expected_n}")
    for c in chains:
        if names[c].get(12) != "ASP":
            raise ValueError(f"chain {c}: residue 12 is {names[c].get(12)}, not ASP (G12D)")
    bad = [(c, r) for c in chains for r in scored if r != 12 and ONE.get(names[c][r]) != topo_seq[r - topo_range[0]]]
    if bad:
        raise ValueError(f"residue names differ from the MD topology sequence: {bad[:5]}")
    return scored


def endpoint(coords: dict, scored: list, p: dict) -> dict:
    """coords: {('GDP', chain): (n, 3), ('GNP', chain): (n, 3)} in the order of `scored`. Returns the arrays and the sets."""
    tf = p["pair_fit"]["trim_fraction"]
    gdp = sorted(k for k in coords if k[0] == "GDP")
    gnp = sorted(k for k in coords if k[0] == "GNP")
    between = {(a[1], b[1]): pair_displacement(coords[a], coords[b], scored, tf) for a in gdp for b in gnp}
    within = {f"{s}:{x[1]}-{y[1]}": pair_displacement(coords[x], coords[y], scored, tf)
              for s, grp in (("GDP", gdp), ("GNP", gnp)) for x, y in itertools.combinations(grp, 2)}
    stack = np.stack([between[k] for k in sorted(between)])
    b_mean, b_med = stack.mean(axis=0), np.median(stack, axis=0)
    b_primary = between[(gdp[0][1], gnp[0][1])]
    w = np.max(np.stack([within[k] for k in sorted(within)]), axis=0)
    out = {"n_between_pairs": len(between), "n_within_pairs": len(within), "b_mean": b_mean, "b_median": b_med,
           "b_primary_pair": b_primary, "w_max": w,
           "w_gdp_entry": np.max(np.stack([v for k, v in within.items() if k.startswith("GDP")]), axis=0),
           "w_gnp_entry": np.max(np.stack([v for k, v in within.items() if k.startswith("GNP")]), axis=0)}
    md, wf = p["min_displacement_A"], p["within_state_factor"]
    thr, pos, bnd = compute_positive(b_mean, w, md, wf)
    out.update(threshold=thr, positive=pos, boundary=bnd)
    sens = {}
    for name, arr in (("primary_pair_only", b_primary), ("median_over_pairs", b_med)):
        _, ps, bd = compute_positive(arr, w, md, wf)
        sens[name] = {"positive": ps, "boundary": bd}
    out["sensitivity"] = sens
    return out


def analyse(text_gdp: str, text_gnp: str, p: dict, topo_seq: str, **kw) -> dict:
    """Pure pipeline from the two PDB texts to the endpoint, the sets, the gate and the sensitivity analyses."""
    e = p["endpoint"]
    ca_a, nm_a = parse_ca(text_gdp)
    ca_b, nm_b = parse_ca(text_gnp)
    chains_a, chains_b = e["structures"]["GDP"]["chains"], e["structures"]["GNP"]["chains"]
    # identity is checked over all eight chains together: observed in all, residue 12, sequence
    ca = {("GDP", c): ca_a[c] for c in ca_a}
    ca.update({("GNP", c): ca_b[c] for c in ca_b})
    names = {("GDP", c): nm_a[c] for c in nm_a}
    names.update({("GNP", c): nm_b[c] for c in nm_b})
    keys = sorted(ca)
    flat_ca = {f"{k[0]}:{k[1]}": ca[k] for k in keys}
    flat_nm = {f"{k[0]}:{k[1]}": names[k] for k in keys}
    expected = sorted([f"GDP:{c}" for c in chains_a] + [f"GNP:{c}" for c in chains_b])
    scored = check_identity(flat_ca, flat_nm, expected, topo_seq, **kw)
    coords = {k: np.stack([ca[k][r] for r in scored]) for k in keys}
    res = endpoint(coords, scored, e)
    merge_gap = e["segments"]["merge_gap"]
    pos = [r for r, f in zip(scored, res["positive"]) if f]
    res.update(scored=scored, E=pos, gate=gate(pos, len(scored), e["gate"], merge_gap), observed={
        f"{k[0]}:{k[1]}": sorted(ca[k]) for k in keys})
    res["sensitivity_sets"] = {}
    for name, d in res["sensitivity"].items():
        s = [r for r, f in zip(scored, d["positive"]) if f]
        res["sensitivity_sets"][name] = {"E": s, "n": len(s), "segments": segments(s, merge_gap),
                                         "jaccard_with_E": (len(set(s) & set(pos)) / len(set(s) | set(pos))) if (s or pos) else None}
    return res


# --- freeze and run --------------------------------------------------------------------------------------------------

def input_files(spec: dict) -> dict:
    out = {}
    for key, pdb in (("GDP", "5US4"), ("GNP", "5USJ")):
        s = spec["structures"][key]
        path = ROOT / s["file"]
        if not path.is_file():
            raise SystemExit(f"STOP: input file missing: {s['file']}")
        out[pdb] = {"path": s["file"], "sha256": sha256_file(path), "spec_sha256": s["sha256"],
                    "matches_spec": sha256_file(path) == s["sha256"], "bytes": path.stat().st_size}
        if not out[pdb]["matches_spec"]:
            raise SystemExit(f"STOP: {pdb} differs from the hash recorded in KRAS-STRUCT-2")
    return out


def _tests_hash() -> dict:
    return {t: sha256_file(ROOT / t) for t in TESTS}


def cmd_freeze(args) -> int:
    if FREEZE.exists() or RECORD.exists():
        raise SystemExit("STOP: a freeze or an endpoint record already exists")
    spec = json.loads(SPEC.read_text())
    rec = {"id": "KRAS-K1-S0-freeze", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "protocol": {"path": "benchmarks/audit/kras_k1/PROTOCOL.md", "sha256": sha256_file(PROTOCOL)},
           "params": {"path": "benchmarks/audit/kras_k1/protocol_params.json", "sha256": sha256_file(PARAMS)},
           "endpoint_code": {"path": "benchmarks/audit/kras_k1/k1_endpoint.py", "sha256": sha256_file(SELF)},
           "tests_sha256": _tests_hash(),
           "evidence_sha256": {str(p.relative_to(ROOT)): sha256_file(p) for p in (PROV, FEAS, SPEC)},
           "inputs": input_files(spec),
           "state": "frozen before any crystal endpoint was computed; no MD trajectory downloaded, opened or scored",
           "endpoint_record_exists": False}
    FREEZE.write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps({"frozen_utc": rec["frozen_utc"], "inputs": {k: v["matches_spec"] for k, v in rec["inputs"].items()}}))
    return 0


def verify_freeze() -> dict:
    f = json.loads(FREEZE.read_text())
    bad = [n for n, ok in (("protocol", sha256_file(PROTOCOL) == f["protocol"]["sha256"]),
                           ("params", sha256_file(PARAMS) == f["params"]["sha256"]),
                           ("endpoint code", sha256_file(SELF) == f["endpoint_code"]["sha256"]),
                           ("tests", _tests_hash() == f["tests_sha256"]),
                           ("evidence", all(sha256_file(ROOT / k) == v for k, v in f["evidence_sha256"].items()))) if not ok]
    if bad:
        raise SystemExit(f"STOP: files differ from the freeze record: {bad}")
    return f


def git_state() -> dict:
    """The commit that introduced the freeze record, HEAD, and any uncommitted change to the frozen files."""
    def git(*a):
        return subprocess.run(["git", "-C", str(ROOT), *a], capture_output=True, text=True, check=True).stdout.strip()
    frozen = [str(p.relative_to(ROOT)) for p in (FREEZE, PROTOCOL, PARAMS, SELF)] + TESTS
    return {"frozen_commit": git("log", "-1", "--format=%H", "--", str(FREEZE.relative_to(ROOT))),
            "head": git("rev-parse", "HEAD"), "uncommitted_frozen_files": git("status", "--porcelain", "--", *frozen)}


def cmd_run(args) -> int:
    f = verify_freeze()
    g = git_state()
    if not g["frozen_commit"] or g["uncommitted_frozen_files"]:
        raise SystemExit(f"STOP: the freeze record must be committed and the frozen files unmodified: {g}")
    if RECORD.exists():
        raise SystemExit("STOP: the endpoint has already been computed once; it is not recomputed")
    t0, c0 = time.time(), time.process_time()
    p = json.loads(PARAMS.read_text())
    spec = json.loads(SPEC.read_text())
    inputs = input_files(spec)
    if inputs != f["inputs"]:
        raise SystemExit("STOP: input identities differ from the freeze record")
    topo_seq = json.loads(FEAS.read_text())["kras_topology_check"]["sequence"]
    res = analyse((ROOT / inputs["5US4"]["path"]).read_text(), (ROOT / inputs["5USJ"]["path"]).read_text(), p, topo_seq)
    sc = res["scored"]
    arr = lambda k: [float(x) for x in res[k]]  # noqa: E731
    context = {"median_b": float(np.median(res["b_mean"])), "median_w": float(np.median(res["w_max"])),
               "max_b": float(res["b_mean"].max()), "max_w": float(res["w_max"].max()),
               "n_b_at_least_min_displacement": int((res["b_mean"] >= p["endpoint"]["min_displacement_A"] - TOL).sum()),
               "n_b_at_least_factor_times_w": int((res["b_mean"] >= p["endpoint"]["within_state_factor"] * res["w_max"] - TOL).sum()),
               "w_exceeds_b_residues": [r for r, x, y in zip(sc, res["w_max"], res["b_mean"]) if x > y],
               "positives_with_threshold_set_by_within_state": [
                   r for r, x, y in zip(sc, res["w_max"], res["positive"])
                   if y and p["endpoint"]["within_state_factor"] * x > p["endpoint"]["min_displacement_A"]]}
    context["w_exceeds_b_count"] = len(context["w_exceeds_b_residues"])
    boundary = [r for r, x in zip(sc, res["boundary"]) if x]
    missing = {k: [r for r in range(1, 170) if r not in set(v)] for k, v in res["observed"].items()}
    record = {"id": "KRAS-K1-S1-endpoint", "computed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "freeze_record_sha256": sha256_file(FREEZE), "frozen_commit": g["frozen_commit"], "head_at_run": g["head"], "frozen": {k: f[k] for k in ("protocol", "params", "endpoint_code", "tests_sha256")},
              "evidence_class": p["endpoint"]["evidence_class"],
              "inputs": inputs, "units": "A (angstrom), float64",
              "n_scored": len(sc), "scored_residues": sc, "residues_not_observed_in_all_chains_1_to_169": missing,
              "n_between_pairs": res["n_between_pairs"], "n_within_pairs": res["n_within_pairs"],
              "per_residue": {"resid": sc, "b_mean": arr("b_mean"), "b_median": arr("b_median"), "b_primary_pair": arr("b_primary_pair"),
                              "w_max": arr("w_max"), "w_gdp_entry_max": arr("w_gdp_entry"), "w_gnp_entry_max": arr("w_gnp_entry"),
                              "threshold": arr("threshold"), "positive": [bool(x) for x in res["positive"]]},
              "E": res["E"], "n_E": len(res["E"]), "boundary_residues": boundary, "gate": res["gate"],
              "context": context, "sensitivity": res["sensitivity_sets"],
              "decision": "K1_ENDPOINT_ADEQUATE" if res["gate"]["passes"] else "K1_ENDPOINT_INADEQUATE",
              "decision_meaning": ("the gate passed: the proxy contains a sufficiently distributed comparison target; this does not "
                                   "show that the endpoint can discriminate methods; owner approval is required before any MD download")
              if res["gate"]["passes"] else "the gate failed: K1 is closed with this reason; no alternative chain, threshold, endpoint or system",
              "resources": {"cpu_s": time.process_time() - c0, "wall_s": time.time() - t0,
                            "max_rss_kb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
              "software": {"python": sys.version.split()[0], "numpy": np.__version__},
              "md_opened": False}
    with open(RECORD, "x") as fh:
        fh.write(json.dumps(record, indent=1) + "\n")
    print(json.dumps({"n_scored": len(sc), "n_E": len(res["E"]), "segments": [s["length"] for s in res["gate"]["segments"]],
                      "gate": res["gate"]["passes"], "decision": record["decision"], "cpu_s": round(record["resources"]["cpu_s"], 3)}))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=("freeze", "run"))
    args = ap.parse_args(argv)
    return {"freeze": cmd_freeze, "run": cmd_run}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
