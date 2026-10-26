"""KIN-PAIR-1 execution: structpair on the selected panel, checked against independent references (PROTOCOL.md section 7).

    .venv/bin/python benchmarks/audit/kin_pair_1/run_study.py freeze   # hashes the panel, protocol and this code
    .venv/bin/python benchmarks/audit/kin_pair_1/run_study.py run      # refuses to run if any hash differs from the freeze

Needs data_pdb/ (written by select_panel.py; SHA-256 verified against the selection record), network access to UniProt for
the three FASTA files (hash-verified), and an interpreter with MDAnalysis for the structpair adapter (TOPOFOLD_ADAPTER_PYTHON).
All numbers go to results/KIN-PAIR-1_results.json; make_report.py reads only that file.
"""

import csv
import datetime as dt
import hashlib
import json
import os
import resource
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import kinpair as kp  # noqa: E402
import select_panel as sel  # noqa: E402

FREEZE2 = HERE / "results/KIN-PAIR-1_freeze_2.json"
RESULTS = HERE / "results/KIN-PAIR-1_results.json"
WORK = HERE / "work"
PAIRS_DIR = HERE / "results/pairs"
DATA_DIR = HERE / "data_pdb"
ADAPTER_PY = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
FROZEN = ["protocol.json", "PROTOCOL.md", "kinpair.py", "select_panel.py", "run_study.py", "results/KIN-PAIR-1_selection.json",
          "results/KIN-PAIR-1_freeze_1.json", "results/KIN-PAIR-1_freeze_1a.json"]
COLS = {"d_kappa": "d_kappa_rad", "d_i_i2": "d_i_i2_angstrom", "rmsd3": "rmsd3_angstrom", "d_tau": "d_tau_rad",
        "d_i_i3": "d_i_i3_angstrom", "rmsd4": "rmsd4_angstrom"}
ROT_A = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], float)      # cyclic permutation of the axes: exact proper rotation
ROT_B = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], float)     # 90 degrees about z: exact proper rotation
T_A, T_B = np.array([12.5, -7.25, 3.0]), np.array([-20.0, 5.0, 40.0])


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def fnum(x):
    return None if x in ("", None) else float(x)


def cpu():
    a, b = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    return a.ru_utime + a.ru_stime + b.ru_utime + b.ru_stime


def label(r):
    return f"{r.chain}:{r.resseq}{r.icode}"


# --- one pair ---------------------------------------------------------------------------------------------------------

def chains_for(P, entity_chains, primary, ref_seq, accession, pairs_keys):
    """Primary chain plus, by the second-chain rule of protocol.json, at most one more chain of the same entity."""
    q = {c: sum(r.read for r in P.residues[c]) for c in entity_chains if c in P.residues}
    order = sorted(q, key=lambda c: (-q[c], c))
    rec = {"candidates": order, "read_counts": q, "second_chain": None, "second_chain_accepted": False, "reason": "single chain in the entity"}
    rest = [c for c in order if c != primary]
    if not rest:
        return [primary], rec
    second = rest[0]
    rec["second_chain"] = second
    m2 = kp.map_chain(P, second, ref_seq, accession)
    if m2["problems"]:
        rec["reason"] = "second chain fails the mapping rules: " + "; ".join(m2["problems"])[:200]
        return [primary], rec
    res2 = {r.key: r for r in P.residues[second]}
    agree = disagree = 0
    for key, uidx in pairs_keys:
        r = res2.get(key)
        if r is not None and r.read and key in m2["index"]:
            agree += m2["index"][key] == uidx
            disagree += m2["index"][key] != uidx
    rec.update(n_same_mapping=agree, n_different_mapping=disagree)
    if disagree or agree < 180:
        rec["reason"] = f"second chain shares the mapping at {agree} positions with {disagree} different (needs 100 % and at least 180)"
        return [primary], rec
    rec.update(second_chain_accepted=True, reason="accepted")
    return [primary, second], rec


def key_str(r):
    return f"{r.resseq}{r.icode}" if r.icode else r.resseq


def transform(text, R, t):
    out = []
    for line in text.splitlines():
        if line[:6] in ("ATOM  ", "HETATM"):
            v = R @ np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) + t
            f = "".join(f"{x:8.3f}" for x in v)
            if len(f) != 24:
                raise ValueError("transformed coordinate does not fit the PDB field")
            line = line[:30] + f + line[54:]
        out.append(line)
    return "\n".join(out) + "\n"


def run_structpair(spec_path, out_dir):
    env = dict(os.environ, PYTHONPATH=str(ROOT / "python"), TOPOFOLD_ADAPTER_PYTHON=ADAPTER_PY)
    c0 = cpu()
    t0 = time.time()
    p = subprocess.run([sys.executable, "-m", "topofold_inputs.structpair", str(spec_path), "--out", str(out_dir)], env=env,
                       capture_output=True, text=True)
    return {"returncode": p.returncode, "stdout": p.stdout.strip()[-400:], "stderr": p.stderr.strip()[-600:], "cpu_seconds": cpu() - c0,
            "wall_seconds": time.time() - t0}


def read_windows(d):
    return list(csv.DictReader(open(Path(d) / "windows.csv")))


def window_dev(rows, w):
    """Worst deviation per metric between a structpair windows.csv and the independent windows; plus validity differences."""
    worst, nval, bad_valid = {}, {}, {}
    for k, col in COLS.items():
        mx, n, bad = 0.0, 0, []
        for j, r in enumerate(rows):
            got, want = fnum(r[col]), w[k][j]
            if (got is None) != (want is None):
                bad.append(j)
            elif got is not None:
                mx = max(mx, abs(got - want))
                n += 1
        worst[k], nval[k], bad_valid[k] = mx, n, bad
    return worst, nval, bad_valid


def disp_check(pair_dir, sets, pairs, tol):
    """Independent trimmed displacement for every chain pairing against chain_variation.csv (channel disp)."""
    f = Path(pair_dir) / "chain_variation.csv"
    if not f.is_file():
        return {"status": "not applicable: one chain per entry (structpair writes no displacement channel)"}
    rows = [r for r in csv.DictReader(open(f)) if r["channel"] == "disp"]
    byj = {int(r["position"]): r for r in rows}
    n = len(pairs)
    per_pair = {}
    for ca, cb in [(a, b) for a in sets["A"] for b in sets["B"]]:
        idx = [j for j in range(n) if sets["A"][ca][j] is not None and sets["B"][cb][j] is not None]
        if len(idx) < 4:
            continue
        XA = np.array([sets["A"][ca][j] for j in idx])
        XB = np.array([sets["B"][cb][j] for j in idx])
        d, _ = kp.trimmed_displacement(XA, XB, 0.2)
        per_pair[(ca, cb)] = dict(zip(idx, d))
    primary = (sets["primary"]["A"], sets["primary"]["B"])
    worst = 0.0
    n_cmp = 0
    bad = []
    for j in range(n):
        vals = sorted(v[j] for v in per_pair.values() if j in v)
        row = byj.get(j)
        if row is None:
            bad.append(j)
            continue
        want = {"between_n": len(vals)}
        if vals:
            m = len(vals)
            want.update(between_min=vals[0], between_max=vals[-1], between_mean=sum(vals) / m,
                        between_median=vals[m // 2] if m % 2 else (vals[m // 2 - 1] + vals[m // 2]) / 2,
                        primary=per_pair[primary].get(j))
        if int(row["between_n"]) != want["between_n"]:
            bad.append(j)
            continue
        for k, v in want.items():
            if k == "between_n" or v is None:
                continue
            got = fnum(row[k])
            if got is None:
                bad.append(j)
                continue
            worst = max(worst, abs(got - v))
            n_cmp += 1
    return {"status": "compared", "n_chain_pairings": len(per_pair), "n_values_compared": n_cmp, "max_abs_deviation": worst,
            "positions_with_count_or_missing_mismatch": bad, "pass": worst <= tol and not bad}


def analyze_pair(pair, protocol, panel_targets, texts_cache, net):
    tol = protocol["tolerances"]["numerical_reference_abs"]
    pid = pair["pair_id"]
    gene = pair["gene"]
    acc = pair["accession"]
    out = {"pair_id": pid, "gene": gene, "accession": acc, "entries": {"A": pair["A"], "B": pair["B"]}, "meta": {"A": pair["meta_A"], "B": pair["meta_B"]},
           "checks": {}, "outcome": None}
    ref_seq, ref_sha = net.uniprot_fasta(acc)
    out["uniprot_fasta_sha256"] = ref_sha
    out["uniprot_fasta_matches_selection_record"] = panel_targets.get(gene, {}).get("uniprot_fasta_sha256") == ref_sha
    if acc not in texts_cache:
        texts_cache[acc] = (ref_seq, ref_sha)
    P, text, mp = {}, {}, {}
    for k in "AB":
        f = DATA_DIR / f"{pair[k]}.pdb"
        text[k] = f.read_text(encoding="latin-1")
        out.setdefault("inputs", {})[k] = {"file": f"data_pdb/{pair[k]}.pdb", "sha256": sha(f), "matches_selection_record": sha(f) == pair[f"sha256_{k}"]}
        if sha(f) != pair[f"sha256_{k}"]:
            out["outcome"] = "INPUT_OR_MAPPING_FAILURE"
            out["failure"] = f"{pair[k]}: file hash differs from the selection record"
            return out
        P[k] = kp.parse_pdb(text[k], "blank_then_A")
        mp[k] = kp.map_chain(P[k], pair[f"chain_{k}"], ref_seq, acc)
    # R1 mapping
    probs = [f"{k}: " + "; ".join(mp[k]["problems"]) for k in "AB" if mp[k]["problems"]]
    corr = kp.correspondence(P["A"], mp["A"], pair["chain_A"], P["B"], mp["B"], pair["chain_B"])
    out["mapping"] = {k: {"chain": pair[f"chain_{k}"], **{x: mp[k][x] for x in ("n_residues", "n_aligned", "n_unaligned", "identity", "coverage", "mismatches",
                                                                                 "dbref_segments", "n_checked_against_dbref", "n_dbref_disagreements", "alignment_score")},
                      "author_minus_uniprot_offsets": sorted({key[0] - idx for key, idx in mp[k]["index"].items() if key[1] == ""})[:12],
                      "modified_residues_not_read": sorted({r.resname for r in P[k].residues[pair[f"chain_{k}"]] if not r.read and r.ca is not None}),
                      "altloc_ca_in_chain": sum(r.altloc_ca for r in P[k].residues[pair[f"chain_{k}"]])} for k in "AB"}
    out["correspondence"] = {"n_positions": corr["n_shared"], "uniprot_span": corr["span"], "n_span_positions_not_shared": corr["n_span_not_shared"]}
    ra = {r.key: r for r in P["A"].residues[pair["chain_A"]]}
    rb = {r.key: r for r in P["B"].residues[pair["chain_B"]]}
    reasons = {"A": {"not in the deposited chain (unobserved or unaligned)": 0, "present but not read by the adapter": 0},
               "B": {"not in the deposited chain (unobserved or unaligned)": 0, "present but not read by the adapter": 0}}
    shared_idx = {i for _, _, i in corr["pairs"]}
    if corr["span"]:
        for i in range(corr["span"][0], corr["span"][1] + 1):
            if i in shared_idx:
                continue
            for k, m, rr in (("A", mp["A"], ra), ("B", mp["B"], rb)):
                keys = [key for key, idx in m["index"].items() if idx == i]
                if not keys:
                    reasons[k]["not in the deposited chain (unobserved or unaligned)"] += 1
                elif not rr[keys[0]].read:
                    reasons[k]["present but not read by the adapter"] += 1
    out["correspondence"]["span_positions_not_shared_by_reason"] = reasons
    # R8 corruption control
    out["checks"]["R8_shifted_mapping"] = {k: {"identity_true": kp.shifted_identity(P[k], pair[f"chain_{k}"], ref_seq, mp[k], 0),
                                                "identity_shifted_by_1": kp.shifted_identity(P[k], pair[f"chain_{k}"], ref_seq, mp[k], 1)} for k in "AB"}
    out["checks"]["R8_shifted_mapping"]["pass"] = all(v["identity_shifted_by_1"] < protocol["tolerances"]["shifted_mapping_identity_must_be_below"]
                                                      for v in (out["checks"]["R8_shifted_mapping"]["A"], out["checks"]["R8_shifted_mapping"]["B"]))
    if probs or corr["n_shared"] < protocol["quality"]["min_shared_positions"]:
        out["outcome"] = "INPUT_OR_MAPPING_FAILURE"
        out["failure"] = "; ".join(probs) or f"only {corr['n_shared']} shared positions"
        return out
    out["checks"]["R1_mapping"] = {"pass": True, "identity": {k: mp[k]["identity"] for k in "AB"}, "dbref_disagreements": {k: mp[k]["n_dbref_disagreements"] for k in "AB"}}
    # chains for the structpair spec
    keysA = [(a, i) for a, _, i in corr["pairs"]]
    keysB = [(b, i) for _, b, i in corr["pairs"]]
    chA, recA = chains_for(P["A"], pair["meta_A"]["target_chains"], pair["chain_A"], ref_seq, acc, keysA)
    chB, recB = chains_for(P["B"], pair["meta_B"]["target_chains"], pair["chain_B"], ref_seq, acc, keysB)
    out["chains"] = {"A": {"used": chA, **recA}, "B": {"used": chB, **recB}}
    pairs = [[key_str(ra[a]), key_str(rb[b])] for a, b, _ in corr["pairs"]]
    wdir = WORK / pid
    shutil.rmtree(wdir, ignore_errors=True)
    wdir.mkdir(parents=True)

    def spec_for(files, outname):
        s = {"name": f"KIN-PAIR-1 {pid}: {pair['A']} vs {pair['B']} ({gene})",
             "structures": {k: ({"label": f"{pair[k]} chain {c[0]}", "file": str(Path(files[k]).resolve()), "chain": c[0], "altloc": "blank_then_A"} if len(c) == 1 else
                                {"label": f"{pair[k]} chains {','.join(c)}", "file": str(Path(files[k]).resolve()), "chains": c, "primary_chain": c[0], "altloc": "blank_then_A"})
                            for k, c in (("A", chA), ("B", chB))},
             "correspondence": {"pairs": pairs}, "on_missing": "exclude", "resname_mismatch": "flag", "altloc_sensitivity": "highest_occupancy"}
        if len(chA) > 1 or len(chB) > 1:
            s["chain_variation"] = {"displacement_trim_fraction": protocol["structpair"]["displacement_trim_fraction"], "advisory": {"min_difference": {}}}
        (wdir / outname).write_text(json.dumps(s, indent=1))
        return wdir / outname

    base = run_structpair(spec_for({k: DATA_DIR / f"{pair[k]}.pdb" for k in "AB"}, "spec.json"), wdir / "out")
    out["structpair"] = {k: base[k] for k in ("returncode", "stdout", "cpu_seconds", "wall_seconds")}
    if base["returncode"] != 0:
        out["structpair"]["stderr"] = base["stderr"]
        out["outcome"] = "IMPLEMENTATION_FAILURE"
        out["failure"] = "structpair failed on inputs that passed the mapping checks"
        return out
    d = wdir / "out"
    # independent windows on the primary pair
    w = kp.independent_windows([ra[a] for a, _, _ in corr["pairs"]], [rb[b] for _, b, _ in corr["pairs"]])
    rows = read_windows(d)
    n = len(corr["pairs"])
    assert len(rows) == n, (len(rows), n)
    worst, nval, badv = window_dev(rows, w)
    out["checks"]["R2_window_set"] = {"pass": not any(badv.values()), "n_positions": n, "n_valid3_independent": int(sum(w["valid3"])), "n_valid4_independent": int(sum(w["valid4"])),
                                      "validity_mismatch_positions": {k: v for k, v in badv.items() if v}}
    out["checks"]["R3_kappa_tau_values"] = {"pass": all(worst[k] <= tol for k in ("d_kappa", "d_i_i2", "d_tau", "d_i_i3")), "max_abs_deviation": {k: worst[k] for k in ("d_kappa", "d_i_i2", "d_tau", "d_i_i3")},
                                            "n_values": {k: nval[k] for k in ("d_kappa", "d_i_i2", "d_tau", "d_i_i3")}}
    out["checks"]["R4_local_rmsd"] = {"pass": all(worst[k] <= tol for k in ("rmsd3", "rmsd4")), "max_abs_deviation": {k: worst[k] for k in ("rmsd3", "rmsd4")},
                                      "n_values": {k: nval[k] for k in ("rmsd3", "rmsd4")}}
    summ = json.loads((d / "summary.json").read_text())
    gdev = abs(summ["global"]["rmsd_qcp_angstrom"] - w["global_rmsd"])
    out["checks"]["R5_global_rmsd"] = {"pass": gdev <= tol and summ["global"]["n_ca"] == n, "structpair": summ["global"]["rmsd_qcp_angstrom"], "independent": w["global_rmsd"],
                                       "abs_deviation": gdev, "n_ca_structpair": summ["global"]["n_ca"], "n_positions": n}
    # R9 exclusions
    exp3 = {j for j in range(1, n - 1) if not w["valid3"][j]}
    exp4 = {j for j in range(1, n - 2) if not w["valid4"][j]}
    exrows = list(csv.DictReader(open(d / "exclusions.csv")))
    got3 = {int(r["position"]) for r in exrows if r["window_size"] == "3"}
    got4 = {int(r["position"]) for r in exrows if r["window_size"] == "4"}
    causes = {}
    for r in exrows:
        causes[r["cause"]] = causes.get(r["cause"], 0) + 1
    out["checks"]["R9_exclusions"] = {"pass": got3 == exp3 and got4 == exp4, "n_rows": len(exrows), "n_excluded_windows_independent": {"3": len(exp3), "4": len(exp4)},
                                      "n_excluded_windows_structpair": {"3": len(got3), "4": len(got4)}, "only_in_structpair": {"3": sorted(got3 - exp3), "4": sorted(got4 - exp4)},
                                      "only_independent": {"3": sorted(exp3 - got3), "4": sorted(exp4 - got4)}, "rows_by_cause": causes}
    # R6 displacement
    sets = {"A": {c: [(r.ca if (r := {x.key: x for x in P["A"].residues[c]}.get(a)) is not None and r.read else None) for a, _, _ in corr["pairs"]] for c in chA},
            "B": {c: [(r.ca if (r := {x.key: x for x in P["B"].residues[c]}.get(b)) is not None and r.read else None) for _, b, _ in corr["pairs"]] for c in chB},
            "primary": {"A": chA[0], "B": chB[0]}}
    out["checks"]["R6_displacement"] = disp_check(d, sets, corr["pairs"], tol)
    XA = np.array([ra[a].ca for a, _, _ in corr["pairs"]])
    XB = np.array([rb[b].ca for _, b, _ in corr["pairs"]])
    disp_trim, disp_all = kp.trimmed_displacement(XA, XB, protocol["structpair"]["displacement_trim_fraction"])
    # R7 invariance
    files = {}
    for k, R, t in (("A", ROT_A, T_A), ("B", ROT_B, T_B)):
        f = wdir / f"moved_{pair[k]}.pdb"
        f.write_text(transform(text[k], R, t))
        files[k] = f
    inv = run_structpair(spec_for(files, "spec_moved.json"), wdir / "out_moved")
    if inv["returncode"] != 0:
        out["checks"]["R7_rigid_motion"] = {"pass": False, "reason": "structpair failed on the moved copies: " + inv["stderr"][-200:]}
    else:
        a, b = read_windows(d), read_windows(wdir / "out_moved")
        mx, nm, patt = 0.0, 0, 0
        for ra_, rb_ in zip(a, b):
            for c in ra_:
                if c.endswith(("_rad", "_angstrom")):
                    x, y = fnum(ra_[c]), fnum(rb_[c])
                    if (x is None) != (y is None):
                        patt += 1
                    elif x is not None:
                        mx, nm = max(mx, abs(x - y)), nm + 1
        s2 = json.loads((wdir / "out_moved/summary.json").read_text())
        gd = abs(s2["global"]["rmsd_qcp_angstrom"] - summ["global"]["rmsd_qcp_angstrom"])
        cv = 0.0
        if (d / "chain_variation.csv").is_file():
            ca_, cb_ = (list(csv.DictReader(open(x / "chain_variation.csv"))) for x in (d, wdir / "out_moved"))
            for r1, r2 in zip(ca_, cb_):
                for c in ("between_min", "between_median", "between_mean", "between_max", "primary"):
                    x, y = fnum(r1.get(c)), fnum(r2.get(c))
                    if x is not None and y is not None and r1["channel"] == "disp":
                        cv = max(cv, abs(x - y))
        out["checks"]["R7_rigid_motion"] = {"pass": mx <= protocol["tolerances"]["invariance_abs"] and gd <= protocol["tolerances"]["invariance_abs"] and not patt and cv <= protocol["tolerances"]["invariance_abs"],
                                            "max_abs_deviation_window_columns": mx, "n_values": nm, "validity_pattern_differences": patt, "global_rmsd_deviation": gd,
                                            "max_abs_deviation_displacement_channel": cv, "rotations": "A: cyclic axis permutation, B: 90 degrees about z; translations (12.5,-7.25,3) and (-20,5,40)"}
    # classification and descriptive summary
    out["classification"] = kp.classify_pair(w, protocol["classification"]["top_k"], protocol["classification"]["min_valid_windows"])
    lab = lambda j: f"{label(ra[corr['pairs'][j][0]])}/{label(rb[corr['pairs'][j][1]])}"   # noqa: E731
    def tops(vals, k=5):
        return [{"position": j, "A_B_residue": lab(j), "uniprot_index": corr["pairs"][j][2], "value": vals[j]} for j in kp.top_k(vals, k)]
    valid3 = [j for j in range(n) if w["rmsd3"][j] is not None]
    dj = [float(disp_trim[j]) for j in range(n)]
    out["descriptive"] = {"n_positions": n, "global_rmsd_all_positions": w["global_rmsd"],
                          "displacement_trimmed_fit": {"max": float(disp_trim.max()), "median": float(np.median(disp_trim)), "argmax": tops(dj, 1)[0],
                                                       "max_all_position_fit": float(disp_all.max())},
                          "local_shape": {"max_rmsd3": max(w["rmsd3"][j] for j in valid3), "max_rmsd4": max(v for v in w["rmsd4"] if v is not None),
                                          "max_d_kappa": max(v for v in w["d_kappa"] if v is not None), "max_d_tau": max(v for v in w["d_tau"] if v is not None)},
                          "top_windows": {k: tops(w[k]) for k in COLS},
                          "spearman_displacement_vs_rmsd3": kp.spearman([dj[j] for j in valid3], [w["rmsd3"][j] for j in valid3]),
                          "top_displacement_positions": tops(dj, 5)}
    out["_series"] = {"w": {k: w[k] for k in COLS}, "disp": dj}
    # chain sensitivity (descriptive)
    cvf = d / "chain_variation.csv"
    if cvf.is_file():
        cvrows = list(csv.DictReader(open(cvf)))
        sens = {}
        for ch in ("d_kappa", "rmsd3", "d_tau", "rmsd4", "disp"):
            rr = [r for r in cvrows if r["channel"] == ch and r["primary"] != ""]
            if rr:
                prim = [float(r["primary"]) for r in rr]
                mean = [float(r["between_mean"]) for r in rr]
                rng = max(float(r["between_max"]) - float(r["between_min"]) for r in rr)
                tp = lambda v: set(sorted(range(len(v)), key=lambda i: (-v[i], i))[:10])   # noqa: E731
                sens[ch] = {"n_positions": len(rr), "max_range_across_chain_pairings": rng, "top10_overlap_primary_vs_mean_of_pairings": len(tp(prim) & tp(mean))}
        out["chain_sensitivity"] = sens
    else:
        out["chain_sensitivity"] = {"status": "not available: one chain per entry in both"}
    out["altloc_sensitivity_structpair"] = summ.get("altloc_sensitivity")
    out["exclusions_compact"] = [{k: r[k] for k in ("window_size", "position", "structure", "cause", "location")} for r in exrows][:60]
    passed = [out["checks"][k]["pass"] for k in ("R1_mapping", "R2_window_set", "R3_kappa_tau_values", "R4_local_rmsd", "R5_global_rmsd", "R6_displacement", "R7_rigid_motion", "R8_shifted_mapping", "R9_exclusions")
              if k in out["checks"] and "pass" in out["checks"][k]]
    out["outcome"] = "PASS" if all(passed) else "IMPLEMENTATION_FAILURE"
    # compact copy
    dst = PAIRS_DIR / pid
    shutil.rmtree(dst, ignore_errors=True)
    dst.mkdir(parents=True)
    for f in ("windows.csv", "residues.csv", "exclusions.csv", "summary.json", "manifest.json", "report.md", "chain_variation.csv", "chains.csv"):
        if (d / f).is_file():
            shutil.copy2(d / f, dst / f)
    (dst / "spec.json").write_text((wdir / "spec.json").read_text())
    out["compact_outputs"] = {f.name: sha(f) for f in sorted(dst.iterdir())}
    return out


def study_outcome(pairs):
    ex = [p for p in pairs if p["outcome"] in ("PASS", "IMPLEMENTATION_FAILURE")]
    if not pairs or not ex:
        return "NOT EXECUTABLE"
    if any(p["outcome"] == "IMPLEMENTATION_FAILURE" for p in ex):
        return "IMPLEMENTATION FAILURE"
    if len(ex) >= 4 and len({p["gene"] for p in ex}) >= 2:
        return "SUPPORTED ON THIS PANEL"
    return "PARTIAL"


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "freeze":
        if FREEZE2.exists():
            print("STOP: freeze record exists")
            return 1
        FREEZE2.write_text(json.dumps({"study": "KIN-PAIR-1", "stage": "analysis code and panel frozen before any structpair execution",
                                       "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), "sha256": {f: sha(HERE / f) for f in FROZEN},
                                       "results_present_at_freeze": sorted(p.name for p in (HERE / "results").glob("KIN-PAIR-1_results*"))}, indent=1) + "\n")
        print("frozen", FREEZE2)
        return 0
    if cmd != "run":
        print(__doc__)
        return 2
    fr = json.loads(FREEZE2.read_text())
    for f, h in fr["sha256"].items():
        if sha(HERE / f) != h:
            print(f"STOP: {f} differs from the freeze")
            return 1
    protocol = json.loads((HERE / "protocol.json").read_text())
    selection = json.loads((HERE / "results/KIN-PAIR-1_selection.json").read_text())
    net = sel.Net()
    t0, c0 = time.time(), cpu()
    pairs, cache = [], {}
    for pair in selection["panel"]:
        print("pair", pair["pair_id"], pair["A"], pair["B"], flush=True)
        r = analyze_pair(pair, protocol, selection["targets"], cache, net)
        print("  ->", r["outcome"], flush=True)
        pairs.append(r)
    series = {p["pair_id"]: p.pop("_series", None) for p in pairs}
    (HERE / "results/KIN-PAIR-1_series.json").write_text(json.dumps(series) + "\n")
    res = {"study": "KIN-PAIR-1", "protocol_version": protocol["protocol_version"], "freeze_utc": fr["utc"], "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "structpair_version": "0.3.0", "pairs": pairs, "study_outcome": study_outcome(pairs),
           "counts": {"pairs_selected": len(selection["panel"]), "pairs_executed": sum(p["outcome"] in ("PASS", "IMPLEMENTATION_FAILURE") for p in pairs),
                      "pairs_pass": sum(p["outcome"] == "PASS" for p in pairs), "input_or_mapping_failures": sum(p["outcome"] == "INPUT_OR_MAPPING_FAILURE" for p in pairs),
                      "implementation_failures": sum(p["outcome"] == "IMPLEMENTATION_FAILURE" for p in pairs), "targets": sorted({p["gene"] for p in pairs})},
           "resources": {"execution_cpu_seconds": cpu() - c0, "execution_wall_seconds": time.time() - t0, "selection_files_downloaded": selection["files_downloaded"],
                         "selection_bytes_downloaded": selection["bytes_downloaded"],
                         "pdb_bytes_on_disk": sum(f.stat().st_size for f in DATA_DIR.glob("*.pdb")),
                         "installation_cost": "not part of this study (existing environment)"},
           "software": {"python": sys.version.split()[0], "numpy": np.__version__, "adapter_python": ADAPTER_PY}}
    RESULTS.write_text(json.dumps(res, indent=1, default=str) + "\n")
    print("outcome:", res["study_outcome"], res["counts"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
