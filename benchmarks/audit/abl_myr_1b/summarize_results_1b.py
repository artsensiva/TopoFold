"""Post-run reader for results/ABL-MYR-1b_results.json: residue indices to isoform 1b residue numbers and contiguous segments.

    python3 benchmarks/audit/abl_myr_1b/summarize_results_1b.py

Reads only the machine-readable results (nothing is recomputed and no trajectory is opened). The topology's author numbering
equals the isoform 1b reference numbering (G3 offset 0), which is checked here. Not part of the frozen wrapper.
"""

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / "results/ABL-MYR-1b_results.json"
OUT = HERE / "results/ABL-MYR-1b_residue_summary.json"


def segments(nums):
    out, cur = [], []
    for n in sorted(set(nums)):
        if cur and n == cur[-1] + 1:
            cur.append(n)
        else:
            if cur:
                out.append(cur)
            cur = [n]
    if cur:
        out.append(cur)
    return [[s[0], s[-1]] if len(s) > 1 else [s[0]] for s in out]


def summarize(res: dict) -> dict:
    offs = res["gates"]["G3"]["alignment"]["author_numbering_diagnostic"]["reference_index_minus_author_resid"]
    if set(offs) != {"0"}:
        raise ValueError(f"author numbering is not the reference numbering: {offs}")
    num = {r["index"]: r["resid"] for r in res["residues"]}
    o = res["outcome"]
    conv = lambda idx: [num[i] for i in idx]  # noqa: E731
    out = {"numbering": "isoform 1b (UniProt P00519-2); topology resid equals the reference index",
           "outcome_label": o["label"], "decision_rule_step": o["decision_rule_step"],
           "topology_residue_range": [min(num.values()), max(num.values())], "n_residues": len(num),
           "sets": {s: {"retained_runs": {c: len(v) for c, v in d["runs"].items()},
                        "leave_one_run_out_recomputations": d["leave_one_run_out_recomputations"],
                        "windows_robust_in_set": {c: v["windows_exceeding"] for c, v in d["channels"].items()}}
                    for s, d in res["sets"].items()},
           "windows_robust_in_both_sets": {c: len(v) for c, v in o["robust_in_both_sets_windows"].items()}}
    for key in ("R_T", "R_B", "R_T_outside_N2_R_B", "R_B_outside_N2_R_T", "R_phi_psi"):
        nums = conv(o[key])
        out[key] = {"residues": nums, "segments": segments(nums), "count": len(nums)}
    both = sorted(set(o["R_T"]) & set(o["R_B"]))
    out["R_T_and_R_B_overlap"] = {"residues": conv(both), "count": len(both)}
    out["alpha_I"] = res["alpha_I"]
    return out


def main() -> int:
    s = summarize(json.loads(RES.read_text()))
    OUT.write_text(json.dumps(s, indent=1) + "\n")
    print(json.dumps({k: s[k] for k in ("outcome_label", "R_T", "R_B", "R_T_outside_N2_R_B", "R_B_outside_N2_R_T")}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
