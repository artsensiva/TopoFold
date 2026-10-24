"""ABL-MYR-1b trajectory-comparison protocol (post-hoc; DRAFT for independent scientific review). No trajectory is opened here.

    python3 benchmarks/audit/abl_myr_1b/write_protocol_1b.py write     # builds protocol_1b.json and PROTOCOL_1B.md
    python3 benchmarks/audit/abl_myr_1b/write_protocol_1b.py freeze    # hashes them into results/ABL-MYR-1b_protocol_freeze.json

Every measured number in the protocol is read from results/ABL-MYR-1b_check_a2_results.json; the text around them is fixed
here. The comparison method, thresholds and baselines are those of the original ABL-MYR-1 specification, restated and
referenced by hash; only the items listed under "differences_from_original" change.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
A2_RESULTS = HERE / "results/ABL-MYR-1b_check_a2_results.json"
PROTOCOL = HERE / "protocol_1b.json"
PROTOCOL_MD = HERE / "PROTOCOL_1B.md"
FREEZE = HERE / "results/ABL-MYR-1b_protocol_freeze.json"
SELF = Path(__file__).resolve()
LINEAGE = {
    "ABL-MYR-1 original specification": "benchmarks/audit/abl_myr_1/run_spec.json",
    "ABL-MYR-1 specification freeze": "benchmarks/audit/abl_myr_1/results/ABL-MYR-1_freeze.json",
    "ABL-MYR-1 runner (frozen)": "benchmarks/audit/abl_myr_1/run_abl_myr.py",
    "ABL-MYR-1 code freeze": "benchmarks/audit/abl_myr_1/results/ABL-MYR-1_code_freeze.json",
    "ABL-MYR-1 blocked G1 record": "benchmarks/audit/abl_myr_1/results/ABL-MYR-1_blocked.json",
    "ABL-MYR-1 member audit": "benchmarks/audit/abl_myr_1/results/ABL-MYR-1_g1_member_audit.json",
    "ABL-MYR-1a audit freeze": "benchmarks/audit/abl_myr_1a/results/ABL-MYR-1a_audit_freeze.json",
    "ABL-MYR-1a audit results (verdict V4)": "benchmarks/audit/abl_myr_1a/results/ABL-MYR-1a_audit_results.json",
    "ABL-MYR-1b check spec": "benchmarks/audit/abl_myr_1b/check_spec.json",
    "ABL-MYR-1b check freeze": "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_freeze.json",
    "ABL-MYR-1b check attempt 1 (S3 stop)": "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_results_attempt1.json",
    "ABL-MYR-1b amendment a1": "benchmarks/audit/abl_myr_1b/check_amendment_a1.json",
    "ABL-MYR-1b a1 results (K1 not established)": "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_a1_results.json",
    "ABL-MYR-1b amendment a2": "benchmarks/audit/abl_myr_1b/check_amendment_a2.json",
    "ABL-MYR-1b a2 freeze": "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_a2_freeze.json",
    "ABL-MYR-1b a2 results (comparable_for_1b)": "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_check_a2_results.json",
    "GROMACS acquisition record": "benchmarks/audit/abl_myr_1b/results/ABL-MYR-1b_gromacs_acquisition.json",
}
LABELS = {
    "blocked": "blocked (post-hoc protocol)",
    "inconclusive": "inconclusive (post-hoc, descriptive)",
    "TopoFold adds reproduced residues": "TopoFold flags residues outside the baseline neighbourhoods, consistent across both thermostat sets (post-hoc, descriptive)",
    "baselines cover TopoFold": "baselines cover TopoFold (post-hoc, descriptive)",
    "partial": "partial (post-hoc, descriptive)",
}


def sha256(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def measured(res: dict) -> dict:
    """Every measured number of the protocol, read from the a2 results."""
    r1, r2, a2 = res["R1_usable_coordinates"], res["R2_production_start_geometry"], res["a2"]
    per = r1["per_run"].values()
    ca = [v["ca_rmsd_to_em_A"] for v in per]
    bonds = [(v["bonds"]["min_A"], v["bonds"]["max_A"]) for v in per]
    out = {"criteria": res["decision"]["criteria"], "outcome": res["outcome"], "runs_validated": r1["runs_ok"], "runs": r1["runs"],
           "production_start_ca_rmsd_to_system_em_A": {"min": min(ca), "max": max(ca), "limit": 1.0},
           "protein_bond_length_range_A": {"min": min(b[0] for b in bonds), "max": max(b[1] for b in bonds)},
           "reader_agreement": {"runs_ok": a2["agreement_ok"], "runs": a2["runs"], "max_abs_dx_nm": a2["max_abs_dx_nm"]["max"],
                                "max_abs_dv_nm_ps": a2["max_abs_dv_nm_ps"]["max"],
                                "box_vs_extract_box_info_nm": a2["box_vs_extract_box_info_nm"]["max"],
                                "velocities_in_both_readers": a2["velocities_in_both"]},
           "gromacs": {"version": a2["gromacs"]["version"], "build": a2["gromacs"]["package"]["build"]}, "sets": {}}
    for s in ("Set1", "Set2"):
        d = {}
        for k in ("within_Myr", "within_noMyr", "between_Myr_noMyr"):
            v = r2[s][k]
            d[k] = {"ca_rmsd_A": {x: v["ca_rmsd_A"][x] for x in ("n", "min", "median", "max")},
                    "heavy_atom_rmsd_median_A": v["heavy_rmsd_A"]["median"], "pairs_below_0.01A": v["pairs_below_0.01A"]}
        d["start_contrast_E_A"] = (d["between_Myr_noMyr"]["ca_rmsd_A"]["median"]
                                   - max(d["within_Myr"]["ca_rmsd_A"]["median"], d["within_noMyr"]["ca_rmsd_A"]["median"]))
        comp = res["R5_composition"][s]["comparison"]
        d["composition_differences"] = comp["differences"]
        d["net_charge_e"] = comp["net_charge_e"]
        d["box_volume_nm3"] = {c: a2["run_input_settings"][f"{s}:{c}"]["box_volume_nm3"] for c in ("Myr", "noMyr")}
        d["run_input_thermostat"] = {c: {"tcoupl": a2["run_input_settings"][f"{s}:{c}"]["ir"]["tcoupl"],
                                         "ref_t": a2["run_input_settings"][f"{s}:{c}"]["ref_t"],
                                         "tau_t": a2["run_input_settings"][f"{s}:{c}"]["tau_t"]} for c in ("Myr", "noMyr")}
        out["sets"][s] = d
    out["cross_set_ca_rmsd_median_A"] = {c: r2["cross_set"][c]["ca_rmsd_A"]["median"] for c in ("Myr", "noMyr")}
    return out


def build(res: dict) -> dict:
    return {
        "id": "ABL-MYR-1b", "date": "2026-10-25", "post_hoc": True,
        "status": "PROTOCOL FROZEN FOR INDEPENDENT SCIENTIFIC REVIEW. No trajectory has been opened or scored. Nothing in this file is a result",
        "kind": "bounded, descriptive comparison of simulated Abl kinase-domain ensembles with and without myristoyl (Myr) from the public de Buhr and Graeter data (Dryad 10.5061/dryad.9cnp5hqnx, CC0)",
        "why_post_hoc": "The original ABL-MYR-1 specification could not run (G1: Set2/system.gro is not published) and its G4 identity rules cannot hold for any Myr vs noMyr pair (ABL-MYR-1a, verdict V4). This protocol was written after the starting files and the 40 production run inputs had been read, but before any trajectory frame was opened. It keeps the original gates, baselines, statistic, thresholds and budget, changes only what is listed under differences_from_original, and restricts its interpretation. The original frozen attempt and every earlier amendment are preserved unchanged",
        "separate_from": "Ayaz et al. WT/M472I, HDX-MS, imatinib resistance, drug response, clinical relevance and biological function. No result of this study is evidence about any of them",
        "lineage_sha256": {k: sha256(ROOT / v) for k, v in LINEAGE.items()},
        "question": "At which kinase-domain windows do the Myr and noMyr simulated ensembles differ by more than replicate runs of the same condition differ, for TopoFold kappa/tau and for equal-information C-alpha baselines, and are those windows found in both thermostat sets? Descriptive only",
        "measured_comparability_evidence": measured(res),
        "what_the_evidence_supports": [
            "both conditions derive from identical pre-minimization protein coordinates in Set1 (ABL-MYR-1a Q1) and the Set2 post-minimization structures are copies of Set1's (ABL-MYR-1a Q4)",
            "the 40 production run inputs were read by two independent readers that agree within the stated tolerances, and pass the frozen K1-K5 criteria",
            "the thermostat targets and coupling settings, tc-group content and composition differ between conditions only in ways that follow from the presence of Myr (group names, MYRP, a small difference in water count)"],
        "what_the_evidence_does_not_support": [
            "that the production starts are identical: they are distinct, with the measured Calpha RMSDs above",
            "that no confound exists: criteria K1-K5 passed, and other differences were not examined (for example the equilibration trajectories, which are not published, and the random velocity seeds)",
            "that the two thermostat sets are independent replicates of the system: they share the post-minimization structures and differ in the thermostat and in random seeds",
            "any causal statement that Myr produces an observed difference: Myr presence co-varies with solvent count and with the equilibration of a different system",
            "anything about the simulation force field or the model's relation to the real protein"],
        "inputs": {"archives": ["PCA.zip", "kinase_kinked_Myr.zip", "kinase_kinked_noMyr.zip"], "hashes": "sizes and SHA-256 of ABL-MYR-1/run_spec.json (verified again at G0)",
                   "trajectories": "PCA.zip: traj_kinked_{Myr,noMyr}_{11..30}.xtc, kinase-domain backbone, one frame per 2 ns; topology KD_noIa_bb_ref.pdb (average.pdb only under the original G2 rule)",
                   "set_definition": "Set1 = replicates 11-20 (v-rescale), Set2 = replicates 21-30 (Nose-Hoover); sets are analysed separately and never pooled; whole runs are the replicate unit"},
        "differences_from_original": {
            "G1_members": "Set2/system.gro is removed from the required members (it is not published); every other member is required as before",
            "G4_replaced_by_G4b": "The original G4 (system.gro identity, production_md.mdp identity, topol.top identity) is replaced by G4b: verify, by hash, that the ABL-MYR-1a audit results show Q1 identical, Q4 identical for both conditions and set2_start_direct_evidence true, and that the ABL-MYR-1b a2 results show outcome comparable_for_1b with all of K1-K5 true; verify the three archive hashes again. G4b does not claim identical production starts",
            "outcome_wording": "the frozen outcome steps and thresholds are unchanged, but the labels are replaced by the restricted labels below, and the word 'reproduced' is not used",
            "added_descriptive_reports": ["per-run production-start descriptors from the a2 results (C-alpha RMSD to system_em.gro, box volume) listed beside the results, never used in the statistic", "the contrast statistic E applied to the production-start C-alpha traces, reported as a negative control (measured values above)", "the sizes of the temporal and run-to-run variability per channel as in the original"]},
        "unchanged_from_original": {"gates": "G0, G2, G3, G5, G6 exactly as in ABL-MYR-1/run_spec.json (including the t >= 1000 ns analysis window, 2 ns interval, 450-frame minimum, 8-run minimum per group, 5 % excluded-window limit, anchors beta3_Lys 290, alphaC_Glu 305, HRD_H 380, DFG_Asp 400 in isoform 1b numbering)",
                                    "method": "topofold_inputs.ensemble analysis with bins 36 and max_lag 50 frames; channels kappa, tau (TopoFold), d(i,i+2), d(i,i+3) and per-residue C-alpha RMSF (equal-information baselines), phi/psi (stronger-information comparator, reported separately); statistic E (cross-condition run-pair distance minus the larger within-condition median); a window is robust if E > 0 in the full set and in all 20 leave-one-run-out recomputations; no p-values, confidence intervals or significance labels",
                                    "thresholds": "outcome steps and the >= 3 residue rule are those of ABL-MYR-1 and are not changed",
                                    "budget": {"download_bytes_max": 1200000000, "storage_bytes_max": 4000000000, "cpu_hours_max": 2.0},
                                    "stop_rule": "the first failed gate stops the run and is reported; no threshold is relaxed afterwards"},
        "outcomes_ordered_with_restricted_labels": [
            {"step": 1, "label": LABELS["blocked"], "rule": "any gate G0-G6 (with G4b in place of G4) fails"},
            {"step": 2, "label": LABELS["inconclusive"], "rule": "no channel has a window that is robust in both sets"},
            {"step": 3, "label": LABELS["TopoFold adds reproduced residues"], "rule": "|R_T \\ N2(R_B)| >= 3, with R_c built from windows robust in both sets"},
            {"step": 4, "label": LABELS["baselines cover TopoFold"], "rule": "R_T is a subset of N2(R_B)"},
            {"step": 5, "label": LABELS["partial"], "rule": "every remaining case"}],
        "interpretation_restrictions": [
            "descriptive evidence about differences between two simulated ensembles of one modelled kinase domain, beyond run-to-run variability, and about which of the compared methods localizes them",
            "'consistent across both thermostat sets' is not independent replication: see what_the_evidence_does_not_support",
            "production starts differ between runs by about the measured C-alpha RMSDs, within and between conditions alike; any difference in dynamics smaller than what such starting variation could produce must not be interpreted",
            "the topology is named KD_noIa_bb_ref: if it does not contain the alphaI helix, the region where Myr binds the C-lobe is not analysed; this is checked at G2 and stated in the report",
            "the analysis window starts at 1000 ns as in the source study; the equilibration and the first microsecond are not examined",
            "no causal claim about Myr, no claim about the real protein, drug response, clinical relevance or any other study's data"],
        "implementation_plan": {"summary": "a thin wrapper that imports the frozen ABL-MYR-1 runner unchanged and replaces exactly G4 (by G4b), the required Set2 member list, and the outcome labels; everything else is the frozen code",
                                "required_before_trajectories": ["independent scientific review of this protocol recorded in the repository", "the wrapper written, tested on SYNTHETIC archives, and hash-frozen (code freeze) before PCA.zip is opened", "verification that the frozen runner, its code freeze and every lineage file above are unchanged"]},
        "review_gate": "No trajectory (PCA.zip) is opened until an independent scientific review has been recorded and the wrapper is code-frozen",
        "questions_for_independent_reviewers": [
            "Is the evidence in the comparability section sufficient for the restricted interpretation, or is a further control needed before opening trajectories?",
            "Are the two thermostat sets usefully distinct for a robustness statement, given identical post-minimization structures?",
            "Does a measured production-start contrast of this size (see the negative-control statistic) threaten the Myr vs noMyr comparison at the stated window?",
            "Is carrying over the original thresholds and the >= 3 residue rule appropriate for a post-hoc protocol, and are the restricted labels adequate?",
            "Is the exclusion of the first 1000 ns and the use of backbone-only 2 ns frames adequate for kappa/tau and the baselines?",
            "Is the possible absence of the alphaI helix from the topology disqualifying for the question asked?"],
        "stop_conditions": "G0-G6 as in the original (G4b for G4); any failure stops the run and is reported; no workaround amendment without an owner decision",
        "not_in_scope": ["Ayaz WT/M472I or HDX", "free energies or kinetics", "causal Myr claims", "a second comparison with different thresholds"]}


def render_md(p: dict) -> str:
    m = p["measured_comparability_evidence"]
    L = ["# ABL-MYR-1b protocol (post-hoc; frozen for independent scientific review)", "",
         f"**{p['status']}.**", "",
         "Generated by `write_protocol_1b.py` from `protocol_1b.json`; every measured number is read from "
         "`results/ABL-MYR-1b_check_a2_results.json`.", "", "## Scope", "", p["kind"] + ".", "",
         "**Separate from:** " + p["separate_from"] + ".", "", "## Why post-hoc", "", p["why_post_hoc"], "", "## Question", "", p["question"], "",
         "## Measured comparability evidence (ABL-MYR-1b check, amendment a2)", "",
         f"- Criteria K1–K5: {', '.join(f'{k} {v}' for k, v in m['criteria'].items())} (outcome `{m['outcome']}`: criteria passed, not proof of absence of confounds).",
         f"- Run inputs validated: {m['runs_validated']} of {m['runs']}. GROMACS {m['gromacs']['version']} ({m['gromacs']['build']}) and the earlier reader agree on all "
         f"{m['reader_agreement']['runs_ok']} files: max |Δx| {m['reader_agreement']['max_abs_dx_nm']} nm, max |Δv| {m['reader_agreement']['max_abs_dv_nm_ps']} nm/ps, box vs "
         f"MDAnalysis `extract_box_info` {m['reader_agreement']['box_vs_extract_box_info_nm']} nm.",
         f"- Production-start Cα RMSD to `system_em.gro` (after superposition): {m['production_start_ca_rmsd_to_system_em_A']['min']}–"
         f"{m['production_start_ca_rmsd_to_system_em_A']['max']} Å (frozen limit {m['production_start_ca_rmsd_to_system_em_A']['limit']} Å); protein bond lengths "
         f"{m['protein_bond_length_range_A']['min']}–{m['protein_bond_length_range_A']['max']} Å.", "",
         "Production-start Cα RMSD after optimal superposition (Å; min / median / max over run pairs; heavy-atom median):", "",
         "| Set | Pairs | Min | Median | Max | Heavy-atom median |", "|---|---|---|---|---|---|"]
    for s, d in m["sets"].items():
        for k, name in (("within_Myr", "within Myr"), ("within_noMyr", "within noMyr"), ("between_Myr_noMyr", "between Myr and noMyr")):
            c = d[k]["ca_rmsd_A"]
            L.append(f"| {s}, {name} | {c['n']} | {c['min']:.3f} | {c['median']:.3f} | {c['max']:.3f} | {d[k]['heavy_atom_rmsd_median_A']:.3f} |")
    L += ["", "Start-contrast statistic E applied to the production-start Cα traces (the analysis statistic's form, as a negative control): " +
          ", ".join(f"{s} {d['start_contrast_E_A']:+.3f} Å" for s, d in m["sets"].items()) + ". No pair of production starts is closer than 0.01 Å "
          f"(pairs below: {sum(d[k]['pairs_below_0.01A'] for d in m['sets'].values() for k in ('within_Myr', 'within_noMyr', 'between_Myr_noMyr'))}). "
          "Set1 vs Set2 Cα medians: " + ", ".join(f"{c} {v:.3f} Å" for c, v in m["cross_set_ca_rmsd_median_A"].items()) + ".", "",
          "Composition and settings read from the run inputs:", ""]
    for s, d in m["sets"].items():
        th = d["run_input_thermostat"]
        L.append(f"- {s}: composition differences between conditions {d['composition_differences']}; net charge {d['net_charge_e']}; thermostat Myr {th['Myr']['tcoupl']} "
                 f"ref_t {th['Myr']['ref_t']}, noMyr {th['noMyr']['tcoupl']} ref_t {th['noMyr']['ref_t']}; box volumes (nm³) {d['box_volume_nm3']}.")
    L += ["", "### What the evidence supports", ""] + [f"- {x}" for x in p["what_the_evidence_supports"]]
    L += ["", "### What it does not support", ""] + [f"- {x}" for x in p["what_the_evidence_does_not_support"]]
    L += ["", "## Inputs", ""] + [f"- {k}: {v}" for k, v in p["inputs"].items()]
    L += ["", "## Differences from the original ABL-MYR-1 specification", ""]
    for k, v in p["differences_from_original"].items():
        L.append(f"- **{k}:** " + (v if isinstance(v, str) else "; ".join(v)))
    L += ["", "## Unchanged from the original", ""]
    for k, v in p["unchanged_from_original"].items():
        L.append(f"- **{k}:** " + (v if isinstance(v, str) else json.dumps(v)))
    L += ["", "## Outcomes (ordered; restricted labels)", "", "| Step | Label | Rule |", "|---|---|---|"]
    L += [f"| {o['step']} | {o['label']} | {o['rule']} |" for o in p["outcomes_ordered_with_restricted_labels"]]
    L += ["", "## Interpretation restrictions", ""] + [f"- {x}" for x in p["interpretation_restrictions"]]
    L += ["", "## Implementation plan and review gate", "", p["implementation_plan"]["summary"] + ".", ""]
    L += [f"- {x}" for x in p["implementation_plan"]["required_before_trajectories"]] + ["", "**" + p["review_gate"] + ".**", "",
                                                                                           "## Questions for independent reviewers", ""]
    L += [f"{i}. {q}" for i, q in enumerate(p["questions_for_independent_reviewers"], 1)]
    L += ["", "## Lineage (SHA-256 recorded in `protocol_1b.json`)", ""] + [f"- {k}: `{v}`" for k, v in LINEAGE.items()] + [""]
    return "\n".join(L) + "\n"


def cmd_write() -> int:
    res = json.loads(A2_RESULTS.read_text())
    if res.get("outcome") != "comparable_for_1b":
        raise SystemExit("STOP: the a2 outcome is not comparable_for_1b; no ABL-MYR-1b protocol is written")
    p = build(res)
    PROTOCOL.write_text(json.dumps(p, indent=1) + "\n")
    PROTOCOL_MD.write_text(render_md(p))
    print("written", PROTOCOL.name, PROTOCOL_MD.name)
    return 0


def cmd_freeze() -> int:
    if FREEZE.exists():
        print("freeze exists")
        return 1
    rec = {"id": "ABL-MYR-1b-PROTOCOL", "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "protocol_sha256": sha256(PROTOCOL), "protocol_md_sha256": sha256(PROTOCOL_MD), "writer_sha256": sha256(SELF),
           "a2_results_sha256": sha256(A2_RESULTS), "lineage_sha256": {k: sha256(ROOT / v) for k, v in LINEAGE.items()},
           "trajectories_opened": False, "xtc_extracted": False,
           "state": "frozen for independent scientific review; no trajectory opened or scored"}
    FREEZE.write_text(json.dumps(rec, indent=2) + "\n")
    print(json.dumps({k: rec[k] for k in ("frozen_utc", "protocol_sha256", "protocol_md_sha256")}, indent=1))
    return 0


def main(argv=None) -> int:
    a = (sys.argv[1:] if argv is None else argv) or ["write"]
    return {"write": cmd_write, "freeze": cmd_freeze}[a[0]]()


if __name__ == "__main__":
    sys.exit(main())
