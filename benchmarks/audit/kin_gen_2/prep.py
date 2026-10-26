"""KIN-GEN-2 preparation (development only): download, map, label, fit the D1D2 groups, choose k, run the correctness gates.

    .venv/bin/python benchmarks/audit/kin_gen_2/prep.py          # needs results/KIN-GEN-2_freeze_1.json; network: RCSB, UniProt
Output: results/dev_model.json (development entries, groups, k per representation, leave-one-kinase-out table, primary baseline, gates, preparation cost).
No evaluation kinase or structure is touched here. Protocol: PROTOCOL.md.
"""

import hashlib
import json
import os
import resource
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "benchmarks/audit/kin_pair_1"))
import kingen2 as g  # noqa: E402
import select_panel as sel  # noqa: E402  (RCSB search only)

DATA = HERE / "data_pdb"
CACHE = HERE / "uniprot_cache"
ADAPTER = os.environ.get("TOPOFOLD_ADAPTER_PYTHON") or str(ROOT / "benchmarks/audit/wpd_local_retrieval/build/venv/bin/python")
RCSB_CRIT = {"entry_criteria": {"method": "X-RAY DIFFRACTION", "max_resolution_angstrom": 2.5}}


def download(entry: str, counters: dict):
    f = DATA / f"{entry}.pdb"
    if not f.is_file():
        DATA.mkdir(exist_ok=True)
        try:
            raw = urllib.request.urlopen(f"https://files.rcsb.org/download/{entry}.pdb", timeout=90).read()
        except Exception:      # noqa: BLE001
            return None
        f.write_bytes(raw)
        counters["bytes"] += len(raw)
        counters["files"] += 1
    return f.read_bytes()


def candidate_chains(info, chains):
    out = []
    for ch in sorted(chains):
        rl = chains[ch]
        if len(rl) < 30:
            continue
        m, ident, n = g.pre_map(info["sequence"], rl)
        if ident >= 0.90 and n >= 30:
            out.append((ch, m))
    return out


def process_kinase(name: str, acc: str, entries: list, counters: dict, log: list):
    """Entries -> accepted records of one kinase. Anchors by the protocol rules; one chain per entry."""
    info = g.uniprot_info(acc, CACHE)
    rec = {"kinase": name, "accession": acc, "problems": info["problems"], "entries_considered": len(entries), "info": {k: v for k, v in info.items() if k not in ("sequence", "domain_sequence")}}
    if info["problems"]:
        log.append({"kinase": name, "status": "excluded", "reason": "; ".join(info["problems"])})
        return rec, []
    parsed = {}
    for e in entries:
        raw = download(e, counters)
        if raw is None:
            log.append({"kinase": name, "entry": e, "status": "excluded", "reason": "no PDB-format file at RCSB"})
            continue
        parsed[e] = (g.sha256(raw), g.parse_chains(raw.decode("latin-1")))
    mapped, cand = [], {}
    for e, (h, chains) in parsed.items():
        cc = candidate_chains(info, chains)
        cand[e] = cc
        if cc:
            mapped.append(cc[0][1])
    glu = g.alphac_glu(info, mapped)
    rec["alphaC_glu"] = glu
    if glu["selected"] is None:
        log.append({"kinase": name, "status": "excluded", "reason": glu["reason"]})
        return rec, []
    out = []
    for e, (h, chains) in parsed.items():
        if not cand[e]:
            log.append({"kinase": name, "entry": e, "status": "excluded", "reason": "no chain maps to the reference (identity >= 0.90, >= 30 aligned)"})
            continue
        first_reason = None
        for ch, _ in cand[e]:
            r = g.analyse_chain(info, glu["selected"], chains[ch])
            first_reason = first_reason or f"chain {ch}: {r['reason']}"
            if r["ok"]:
                out.append({"id": e, "kinase": name, "accession": acc, "chain": ch, "sha256": h, **{k: r[k] for k in ("window_ca", "window_ids", "window_resnames", "D1", "D2", "anchor_residues", "identity")}})
                break
        else:
            log.append({"kinase": name, "entry": e, "status": "excluded", "reason": first_reason})
    return rec, out


def main():
    t0, c0 = time.time(), resource.getrusage(resource.RUSAGE_SELF)
    bad = g.check_freeze(HERE / "results/KIN-GEN-2_freeze_1.json")
    if bad:
        print("STOP: frozen file differs:", bad)
        return 1
    proto = json.loads((HERE / "protocol.json").read_text())
    amend = json.loads((HERE / "amendment_1.json").read_text())
    dev = proto["development"]
    net, counters, log = sel.Net(), {"bytes": 0, "files": 0}, []
    entries_all, kin_rec = [], {}
    for name, acc in dev["kinases"].items():
        ids, total = net.search(acc, RCSB_CRIT)
        take = ids[: dev["max_entries_per_kinase"]]
        for c, k in (("1IEP", "ABL1"), ("1KV2", "MAPK14"), ("1UWH", "BRAF")):
            if k == name and c not in take:
                take.append(c)
        rec, out = process_kinase(name, acc, take, counters, log)
        rec["search_total"] = total
        kin_rec[name] = rec
        entries_all += out
        print(name, len(take), "->", len(out), flush=True)
    pts = np.array([[e["D1"], e["D2"]] for e in entries_all])
    groups = g.fit_groups(pts)
    for e, lab in zip(entries_all, groups["labels"]):
        e["group"] = lab
    gates = {}
    # G6 group sizes and G5 controls
    gates["G6"] = {"pass": groups["sizes"].get("out", 0) >= 8 and groups["sizes"].get("inter", 0) >= 3, "sizes": groups["sizes"], "centroids": groups["centroids"], "ranges": groups["ranges"]}
    ctrl = {c: next((e["group"] for e in entries_all if e["id"] == c), "not testable") for c in amend["controls"]["out"]}
    testable = {c: v for c, v in ctrl.items() if v != "not testable"}
    pka = [e["group"] for e in entries_all if e["accession"] == "P17612"]
    pka_in = sum(x == "in" for x in pka) / len(pka) if pka else 0.0
    gates["G5"] = {"pass": len(testable) >= 2 and all(v == "out" for v in testable.values()) and pka_in >= 0.9, "out_controls_assigned": ctrl, "PRKACA_entries": len(pka), "PRKACA_fraction_in": pka_in, "rule": "amendments 1 and 2"}
    # G1 core kappa/tau against the numpy reference over every development window
    dk = dt = 0.0
    for e in entries_all:
        kap, tau, rk, rt = g.core_kappa_tau(np.array(e["window_ca"]))
        dk, dt = max(dk, float(np.abs(kap - rk).max())), max(dt, float(np.abs(np.arctan2(np.sin(tau - rt), np.cos(tau - rt))).max()))
    gates["G1"] = {"pass": max(dk, dt) <= 1e-4, "max_kappa_deviation": dk, "max_tau_deviation_circular": dt, "n_windows": len(entries_all)}
    # G2 rigid-motion invariance
    rng = np.random.default_rng(1)
    worst = 0.0
    for e in entries_all[:60]:
        ca = np.array(e["window_ca"])
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        Q *= np.sign(np.linalg.det(Q))
        moved = ca @ Q.T + rng.normal(size=3) * 30
        a, b = g.raw_features(ca), g.raw_features(moved)
        w64 = max(max(float(np.abs(a[r] - b[r]).max()) for r in ("B1", "B3")), abs(g.rmsd(ca, moved)))
        wT = float(np.abs(a["T"] - b["T"]).max())
        worst = (max(worst[0], w64), max(worst[1], wT)) if isinstance(worst, tuple) else (w64, wT)
    gates["G2"] = {"pass": worst[0] <= amend["G2"]["float64_representations_B1_B2_B3_max"] and worst[1] <= amend["G2"]["T_max_float32_core"], "max_deviation_float64_representations": worst[0],
                   "max_deviation_T": worst[1], "limits": amend["G2"], "n_windows": min(60, len(entries_all))}
    # G4 role residues, shifted control
    bad = [e["id"] for e in entries_all if e["window_resnames"][3:6] != ["ASP", "PHE", "GLY"]]
    shifted = sum(1 for e in entries_all if e["window_resnames"][4:7] == ["ASP", "PHE", "GLY"])
    gates["G4"] = {"pass": not bad and shifted == 0, "windows_without_DFG_at_the_anchor": bad, "windows_where_a_shift_by_one_still_reads_DFG": shifted}
    # G3 and G8 against MDAnalysis (adapter interpreter)
    sample = entries_all[:: max(1, len(entries_all) // 20)][:20]
    pairs = [(int(i), int(j)) for i, j in zip(rng.integers(0, len(entries_all), 50), rng.integers(0, len(entries_all), 50))]
    job = {"pairs": [[entries_all[i]["window_ca"], entries_all[j]["window_ca"]] for i, j in pairs],
           "d1d2": [{"file": str(DATA / f"{e['id']}.pdb"), "chain": e["chain"], "anchors": e["anchor_residues"], "D1": e["D1"], "D2": e["D2"]} for e in sample]}
    (HERE / "work").mkdir(exist_ok=True)
    (HERE / "work/gate_job.json").write_text(json.dumps(job))
    r = subprocess.run([ADAPTER, str(HERE / "gate_mda.py"), str(HERE / "work/gate_job.json"), str(HERE / "work/gate_out.json")], capture_output=True, text=True)
    if r.returncode == 0:
        go = json.loads((HERE / "work/gate_out.json").read_text())
        mine = [g.rmsd(entries_all[i]["window_ca"], entries_all[j]["window_ca"]) for i, j in pairs]
        d3 = float(np.max(np.abs(np.array(mine) - np.array(go["rmsd"]))))
        gates["G3"] = {"pass": d3 <= 1e-4, "max_deviation": d3, "n_pairs": len(pairs)}
        gates["G8"] = {"pass": go["max_d1d2_deviation"] <= 1e-3, "max_deviation": go["max_d1d2_deviation"], "n_entries": len(sample)}
    else:
        gates["G3"] = gates["G8"] = {"pass": False, "reason": r.stderr[-300:]}
    # development model: in / out only; k by leave-one-kinase-out
    used = [e for e in entries_all if e["group"] in ("in", "out")]
    model = g.Model(used, [e["group"] for e in used])
    table = g.loko(model, [e["kinase"] for e in used])
    k_sel = {rep: sorted(table[rep], key=lambda k: (-table[rep][k]["ba"], k))[0] for rep in g.REPS}
    base = max(("B1", "B2", "B3"), key=lambda r: (table[r][k_sel[r]]["ba"], -("B1", "B2", "B3").index(r)))
    cost = {"download_bytes": counters["bytes"], "files": counters["files"], "wall_seconds": time.time() - t0,
            "cpu_seconds": resource.getrusage(resource.RUSAGE_SELF).ru_utime + resource.getrusage(resource.RUSAGE_SELF).ru_stime - c0.ru_utime - c0.ru_stime}
    per_kin = {k: {"in": sum(e["kinase"] == k and e["group"] == "in" for e in entries_all), "out": sum(e["kinase"] == k and e["group"] == "out" for e in entries_all),
                   "inter": sum(e["kinase"] == k and e["group"] == "inter" for e in entries_all)} for k in dev["kinases"]}
    out = {"id": "KIN-GEN-2-DEVELOPMENT-MODEL", "route": "B (independent D1D2 groups; not KinCoRe labels)", "kinases": kin_rec, "per_kinase_group_counts": per_kin, "entries": entries_all,
           "groups": {k: groups[k] for k in ("centroids", "sizes", "ranges")}, "gates": gates, "loko": {rep: {str(k): v for k, v in t.items()} for rep, t in table.items()},
           "k_selected": k_sel, "primary_baseline": base, "development_ids": [e["id"] for e in entries_all], "development_kinases": list(dev["kinases"].values()),
           "exclusions": log, "preparation_cost": cost, "all_gates_pass": all(v["pass"] for v in gates.values())}
    (HERE / "results/dev_model.json").write_text(json.dumps(out, indent=1, default=str) + "\n")
    print(json.dumps({"n_dev": len(entries_all), "groups": groups["sizes"], "gates": {k: v["pass"] for k, v in gates.items()}, "k": k_sel, "baseline": base,
                      "loko_ba": {r: round(table[r][k_sel[r]]["ba"], 3) for r in g.REPS}, "cost": cost}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
