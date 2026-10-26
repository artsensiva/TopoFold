"""KIN-GEN-2 deciding evaluation (protocol frozen: PROTOCOL.md, amendments, freeze 2).

    .venv/bin/python benchmarks/audit/kin_gen_2/evaluate.py freeze   # hashes code, development model and protocol files
    .venv/bin/python benchmarks/audit/kin_gen_2/evaluate.py run      # refuses to run if anything differs from the freeze; network: RCSB, UniProt

Selects evaluation kinases and entries deterministically, labels them with the development D1D2 groups (never refitted), scores every representation with the
development-only k-NN procedure, and applies the outcome rules. Evidence scope: structural-label recovery only. Output: results/eval_results.json.
"""

import datetime as dt
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import kingen2 as g  # noqa: E402
import prep  # noqa: E402  (process_kinase, download, RCSB search settings)
import select_panel as sel  # noqa: E402

FREEZE2 = HERE / "results/KIN-GEN-2_freeze_2.json"
FROZEN = ["kingen2.py", "prep.py", "evaluate.py", "gate_mda.py", "protocol.json", "amendment_1.json", "amendment_2.json", "results/dev_model.json", "results/KIN-GEN-2_freeze_1.json",
          "results/reference_source_check.json"]
PANEL_TARGETS = {"P00533", "P11309", "Q06187"}       # KIN-PAIR-1 panel (EGFR, PIM1, BTK)


def cpu():
    r = resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime + r.ru_stime


def domain_identity(a: str, b: str) -> float:
    """Kinase-domain identity: identical aligned positions of a local alignment over the length of the shorter domain."""
    score, pairs = g._align(a, b, np.zeros(len(b), dtype=bool))
    return sum(a[i] == b[j] for i, j in pairs) / min(len(a), len(b))


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "freeze":
        import hashlib
        if FREEZE2.exists():
            print("STOP: freeze record exists")
            return 1
        FREEZE2.write_text(json.dumps({"study": "KIN-GEN-2", "stage": "freeze 2: analysis code, development model and selection code before any evaluation structure is selected or scored",
                                       "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                                       "sha256": {f: hashlib.sha256((HERE / f).read_bytes()).hexdigest() for f in FROZEN},
                                       "evaluation_files_present": sorted(p.name for p in (HERE / "data_pdb").glob("*.pdb"))}, indent=1) + "\n")
        print("frozen")
        return 0
    if cmd != "run":
        print(__doc__)
        return 2
    import hashlib
    fr = json.loads(FREEZE2.read_text())
    bad = [f for f, h in fr["sha256"].items() if hashlib.sha256((HERE / f).read_bytes()).hexdigest() != h]
    if bad or g.check_freeze(HERE / "results/KIN-GEN-2_freeze_1.json"):
        print("STOP: frozen files differ:", bad, g.check_freeze(HERE / "results/KIN-GEN-2_freeze_1.json"))
        return 1
    t0, c0 = time.time(), cpu()
    proto = json.loads((HERE / "protocol.json").read_text())
    ev = proto["evaluation"]
    dm = json.loads((HERE / "results/dev_model.json").read_text())
    if not dm["all_gates_pass"]:
        print("STOP: a development gate failed")
        return 1
    used = [e for e in dm["entries"] if e["group"] in ("in", "out")]
    model = g.Model(used, [e["group"] for e in used], dm["k_selected"])
    groups = {"centroids": dm["groups"]["centroids"]}
    dev_acc, dev_ids = set(dm["development_kinases"]), set(dm["development_ids"])
    expo = json.loads((ROOT / "benchmarks/audit/kin_pair_1/results/exposure_history.json").read_text())
    expo_acc, expo_ids = set(expo["exposed_uniprot_accessions"]), {k.upper() for k in expo["entries"]}
    cands = json.loads((ROOT / "benchmarks/audit/kin_pair_1/results/candidate_targets_verified.json").read_text())
    net, counters, log = sel.Net(), {"bytes": 0, "files": 0}, []
    ranking = []
    for gene, v in cands.items():
        acc = v["accession"]
        if not v["verified"] or acc in expo_acc or acc in PANEL_TARGETS or acc in dev_acc:
            ranking.append({"gene": gene, "accession": acc, "status": "excluded", "reason": "exposed, KIN-PAIR-1 panel target or development kinase"})
            continue
        ids, total = net.search(acc, prep.RCSB_CRIT)
        ranking.append({"gene": gene, "accession": acc, "status": "candidate", "search_total": total, "ids": ids})
    pool = sorted((r for r in ranking if r["status"] == "candidate"), key=lambda r: (-r["search_total"], r["accession"]))
    entries, kin_rec, processed = [], {}, []
    dev_info = {n: g.uniprot_info(a, prep.CACHE) for n, a in zip(dm["kinases"], dm["development_kinases"])}
    for r in pool:
        n_out = sum(e["group"] == "out" for e in entries)
        k_out = len({e["kinase"] for e in entries if e["group"] == "out"})
        if len(processed) >= ev["max_kinases"] or (n_out >= 15 and k_out >= 3) or counters["bytes"] > ev["max_download_mb"] * 1e6:
            r["status"] = "not processed (stopping rule)"
            continue
        ids = [i for i in r["ids"] if i.upper() not in expo_ids and i not in dev_ids][: ev["max_entries_per_kinase"]]
        rec, out = prep.process_kinase(r["gene"], r["accession"], ids, counters, log)
        for e in out:
            e["group"] = g.assign_group([e["D1"], e["D2"]], groups)
        rec["search_total"] = r["search_total"]
        r["status"] = "processed" if out else "processed, no usable entry"
        r["n_usable"] = len(out)
        kin_rec[r["gene"]] = rec
        processed.append(r["gene"])
        entries += out
        print(r["gene"], len(ids), "->", len(out), flush=True)
    # family relation to the development kinases
    fam = {}
    for gene in kin_rec:
        info = g.uniprot_info(kin_rec[gene]["accession"], prep.CACHE)
        if info["problems"]:
            continue
        idn = {n: domain_identity(info["domain_sequence"], di["domain_sequence"]) for n, di in dev_info.items() if not di["problems"]}
        fam[gene] = {"max_identity": max(idn.values()), "nearest_development_kinase": max(idn, key=idn.get), "identities": idn, "family_flag": max(idn.values()) >= ev["family_flag_identity"]}
    # G7 disjointness
    g7 = {"pass": not ({e["accession"] for e in entries} & dev_acc) and not ({e["id"] for e in entries} & dev_ids) and not ({e["id"] for e in entries} & expo_ids), "n_evaluation_entries": len(entries)}
    assert not (set(model.ids) & {e["id"] for e in entries})            # the model holds development ids only
    res = {"id": "KIN-GEN-2-EVALUATION", "route": "B (D1D2 groups fitted on development; not KinCoRe labels)", "evidence_scope": "structural-label recovery only", "ranking": [{k: v for k, v in r.items() if k != "ids"} for r in ranking],
           "processed_kinases": processed, "kinases": kin_rec, "family_relation": fam, "exclusions": log, "G7": g7,
           "entries": [{k: v for k, v in e.items()} for e in entries]}
    scored = [e for e in entries if e["group"] in ("in", "out")]
    cnt = {k: {"in": sum(e["kinase"] == k and e["group"] == "in" for e in entries), "out": sum(e["kinase"] == k and e["group"] == "out" for e in entries), "inter": sum(e["kinase"] == k and e["group"] == "inter" for e in entries)} for k in processed}
    res["counts_per_kinase"] = cnt
    n_out = sum(e["group"] == "out" for e in scored)
    k_out = sum(1 for c in cnt.values() if c["out"] >= 1)
    k_in5 = sum(1 for c in cnt.values() if c["in"] >= 5)
    insufficient = (not g7["pass"]) or n_out < 10 or k_out < 3 or k_in5 < 3
    res["eligibility"] = {"n_scored": len(scored), "n_out": n_out, "n_in": len(scored) - n_out, "n_inter_excluded_from_scoring": sum(e["group"] == "inter" for e in entries), "kinases_with_out": k_out, "kinases_with_5_in": k_in5, "insufficient": insufficient}
    ps = {}
    if scored:
        for rep in g.REPS:
            k = dm["k_selected"][rep]
            rows = []
            for e in scored:
                s = model.score(rep, np.array(e["window_ca"]), k)
                rows.append({"id": e["id"], "kinase": e["kinase"], "y": 1 if e["group"] == "out" else 0, "score": s, "pred": 1 if s >= 0.5 else 0})
            wb = g.weighted_ba(rows)
            y, pr = [r["y"] for r in rows], [r["pred"] for r in rows]
            pooled = g.balanced_accuracy(y, pr)
            per = {}
            for kin in processed:
                rr = [r for r in rows if r["kinase"] == kin]
                per[kin] = {"n": len(rr), "n_out": sum(r["y"] for r in rr), "tpr": (float(np.mean([r["pred"] == 1 for r in rr if r["y"] == 1])) if any(r["y"] == 1 for r in rr) else None),
                            "tnr": (float(np.mean([r["pred"] == 0 for r in rr if r["y"] == 0])) if any(r["y"] == 0 for r in rr) else None)}
            dist = [r for r in rows if not fam.get(r["kinase"], {"family_flag": True})["family_flag"]]
            ps[rep] = {"k": k, "weighted": wb, "pooled_ba": pooled[0], "pooled_tpr": pooled[1], "pooled_tnr": pooled[2], "auroc": g.auroc(y, [r["score"] for r in rows]), "per_kinase": per,
                       "distant_kinases_only": g.weighted_ba(dist) if dist else None, "n_distant_entries": len(dist), "rows": rows}
    res["scores"] = ps
    base = dm["primary_baseline"]
    outcome, detail = "insufficient evidence", {}
    if ps and not insufficient:
        d = ps["T"]["weighted"]["ba_w"] - ps[base]["weighted"]["ba_w"]
        both = [k for k, v in ps["T"]["per_kinase"].items() if v["tpr"] is not None and v["tnr"] is not None]
        wins = sum(1 for k in both if (ps["T"]["per_kinase"][k]["tpr"] + ps["T"]["per_kinase"][k]["tnr"]) >= (ps[base]["per_kinase"][k]["tpr"] + ps[base]["per_kinase"][k]["tnr"]))
        detail = {"delta_ba_w_T_minus_primary_baseline": d, "primary_baseline": base, "kinases_with_both_classes": both, "T_at_least_baseline_in": wins}
        if abs(d) <= 0.05:
            outcome = "comparable"
        elif d >= 0.10 and ps["T"]["weighted"]["ba_w"] >= 0.70 and len(both) >= 3 and wins >= (2 / 3) * len(both):
            outcome = "advantage"
        elif d <= -0.10:
            outcome = "baseline advantage"
        else:
            outcome = "unclear"
    best_ca = max(ps[r]["weighted"]["ba_w"] for r in g.REPS) if ps else float("nan")
    res["outcome"] = {"label": outcome, "detail": detail, "best_ca_only_ba_w": best_ca, "ca_only_limitation_supported": bool(ps) and best_ca < 0.80,
                      "full_atom_reference_note": "the full-atom D1/D2 rule defines the labels, so its agreement with them is 1 by construction; it is not a performance result"}
    res["resources"] = {"evaluation_download_bytes": counters["bytes"], "files": counters["files"], "cpu_seconds": cpu() - c0, "wall_seconds": time.time() - t0, "limits": {"download_mb": 100, "cpu_hours": 2.0},
                        "preparation_recorded_separately": dm["preparation_cost"]}
    (HERE / "results/eval_results.json").write_text(json.dumps(res, indent=1, default=str) + "\n")
    print(json.dumps({"processed": processed, "eligibility": res["eligibility"], "outcome": res["outcome"]["label"], "ba_w": {r: round(ps[r]["weighted"]["ba_w"], 3) for r in ps}}, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
