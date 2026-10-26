"""KIN-PAIR-1 supplementary analysis (protocol amendment 2, post-hoc): the 4 pairs on which the optional highest_occupancy
sensitivity policy refused an occupancy tie. Same driver and checks as the primary run, spec without `altloc_sensitivity`;
plus an independent altloc-B sensitivity of the window metrics (computed here, outside structpair).

    .venv/bin/python benchmarks/audit/kin_pair_1/supplement_altloc_ties.py
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import kinpair as kp  # noqa: E402
import run_study as rs  # noqa: E402
import select_panel as sel  # noqa: E402

OUT = HERE / "results/KIN-PAIR-1_supplement.json"


def choose_b(recs, policy):
    """Alternate location B where the Cα has one, else blank, else A."""
    if not recs:
        return None
    for want in ("B", "", "A"):
        for alt, _, xyz in recs:
            if alt == want:
                return xyz
    return recs[0][2]


def main():
    primary = json.loads((HERE / "results/KIN-PAIR-1_results.json").read_text())
    selection = json.loads((HERE / "results/KIN-PAIR-1_selection.json").read_text())
    protocol = json.loads((HERE / "protocol.json").read_text())
    todo = [p["pair_id"] for p in primary["pairs"] if p["outcome"] == "IMPLEMENTATION_FAILURE" and "tie on occupancy" in p.get("structpair", {}).get("stderr", "")]
    rs.WORK = HERE / "work_supplement"
    rs.PAIRS_DIR = HERE / "results/pairs_supplement"
    orig = rs.run_structpair

    def without_sensitivity(spec_path, out_dir):
        s = json.loads(Path(spec_path).read_text())
        s.pop("altloc_sensitivity", None)
        Path(spec_path).write_text(json.dumps(s, indent=1))
        return orig(spec_path, out_dir)
    rs.run_structpair = without_sensitivity
    net = sel.Net()
    res, t0, c0 = [], time.time(), rs.cpu()
    for pair in selection["panel"]:
        if pair["pair_id"] not in todo:
            continue
        print("pair", pair["pair_id"], flush=True)
        r = rs.analyze_pair(pair, protocol, selection["targets"], {}, net)
        series = r.pop("_series", None)
        # independent altloc-B sensitivity
        ref_seq, _ = net.uniprot_fasta(pair["accession"])
        P = {k: kp.parse_pdb((rs.DATA_DIR / f"{pair[k]}.pdb").read_text(encoding="latin-1"), "blank_then_A") for k in "AB"}
        mp = {k: kp.map_chain(P[k], pair[f"chain_{k}"], ref_seq, pair["accession"]) for k in "AB"}
        corr = kp.correspondence(P["A"], mp["A"], pair["chain_A"], P["B"], mp["B"], pair["chain_B"])
        keep = kp._choose
        kp._choose = choose_b
        try:
            Q = {k: kp.parse_pdb((rs.DATA_DIR / f"{pair[k]}.pdb").read_text(encoding="latin-1"), "blank_then_A") for k in "AB"}
        finally:
            kp._choose = keep
        ra, rb = ({x.key: x for x in P[k].residues[pair[f"chain_{k}"]]} for k in "AB")
        qa, qb = ({x.key: x for x in Q[k].residues[pair[f"chain_{k}"]]} for k in "AB")
        w0 = kp.independent_windows([ra[a] for a, _, _ in corr["pairs"]], [rb[b] for _, b, _ in corr["pairs"]])
        w1 = kp.independent_windows([qa[a] for a, _, _ in corr["pairs"]], [qb[b] for _, b, _ in corr["pairs"]])
        moved = [j for j, (a, b, _) in enumerate(corr["pairs"]) if ra[a].altloc_ca or rb[b].altloc_ca]
        sens = {"n_positions_with_altloc_ca": len(moved), "positions": moved[:20],
                "residues": [f"{ra[corr['pairs'][j][0]].chain}:{ra[corr['pairs'][j][0]].resseq}" for j in moved[:20]],
                "validity_sets_equal": w0["valid3"] == w1["valid3"] and w0["valid4"] == w1["valid4"], "max_abs_change": {}, "top10_overlap": {}}
        for k in rs.COLS:
            d = [abs(x - y) for x, y in zip(w0[k], w1[k]) if x is not None and y is not None]
            sens["max_abs_change"][k] = max(d) if d else None
            sens["top10_overlap"][k] = len(set(kp.top_k(w0[k], 10)) & set(kp.top_k(w1[k], 10)))
        sens["global_rmsd_change"] = abs(w0["global_rmsd"] - w1["global_rmsd"])
        r["independent_altloc_B_sensitivity"] = sens
        res.append(r)
        print("  ->", r["outcome"], flush=True)
    out = {"study": "KIN-PAIR-1", "analysis": "supplementary (post-hoc, protocol amendment 2)", "pairs": res,
           "cpu_seconds": rs.cpu() - c0, "wall_seconds": time.time() - t0,
           "note": "reported separately from the frozen primary outcome, which is not changed by this analysis"}
    OUT.write_text(json.dumps(out, indent=1, default=str) + "\n")
    print([(p["pair_id"], p["outcome"]) for p in res])
    return 0


if __name__ == "__main__":
    sys.exit(main())
