"""SYNTHETIC tests of the ABL-MYR-1b wrapper (benchmarks/audit/abl_myr_1b/run_abl_myr_1b.py). Engineering checks only: archives
are generated here and nothing is read from Dryad. The pure tests run everywhere; the end-to-end tests need the TopoFold core
and TOPOFOLD_ADAPTER_PYTHON (an MDAnalysis interpreter) and are skipped otherwise."""

import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_myr_1"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_myr_1b"))
import run_abl_myr as R0  # noqa: E402
import run_abl_myr_1b as B  # noqa: E402
import test_abl_myr_runner as T  # noqa: E402
import write_protocol_1b as W  # noqa: E402

AUDIT = json.loads(B.AUDIT_RESULTS.read_text())
A2 = json.loads(B.A2_RESULTS.read_text())
EVIDENCE = {"audit": AUDIT, "a2": A2}
ADAPTER = T.ADAPTER
e2e = T.e2e


# --- G0 lineage and the freeze guards -----------------------------------------------------------------------------------

def test_real_lineage_is_verified_and_tampering_is_detected(tmp_path):
    assert B.lineage_problems() == []
    f = tmp_path / "x.txt"
    f.write_text("a")
    rec = {"x.txt": R0.sha256_file(f), "missing.txt": "0" * 64}
    assert B.hash_mismatches(rec, tmp_path) == ["missing.txt"]
    f.write_text("b")
    assert B.hash_mismatches(rec, tmp_path) == ["x.txt", "missing.txt"]


def test_amendment_freeze_precedes_the_wrapper_and_records_the_review():
    f = json.loads(B.AMEND_FREEZE.read_text())
    assert f["review_kind"].startswith("independent AI scientific review") and f["trajectories_extracted"] is False
    for k in ("benchmarks/audit/abl_myr_1b/protocol_amendment_a1.json", "benchmarks/audit/abl_myr_1b/REVIEW_2026-10-25_independent_AI.md"):
        assert k in f["sha256"]
    assert f["frozen_utc"] < json.loads(B.FREEZE.read_text())["code_frozen_utc"] if B.FREEZE.exists() else True
    assert "not human peer review" in (ROOT / "benchmarks/audit/abl_myr_1b/REVIEW_2026-10-25_independent_AI.md").read_text()


def test_run_and_freeze_refuse_without_prerequisites(tmp_path, monkeypatch):
    monkeypatch.setattr(B, "FREEZE", tmp_path / "none.json")
    with pytest.raises(SystemExit, match="no wrapper freeze"):
        B.verify_freeze()
    monkeypatch.setattr(B, "trajectories_present", lambda: True)
    with pytest.raises(SystemExit, match="must precede"):
        B.cmd_freeze_code(None)


def test_amendment_content_matches_the_wrapper_constants():
    a = json.loads(B.AMENDMENT.read_text())
    assert "fewer than 8 usable runs" in a["changes"]["C1_bilateral_G5"] and "every retained run" in a["changes"]["C2_leave_one_run_out_count"]
    assert f"{B.ALPHA_I[0]}-{B.ALPHA_I[1]}" in a["changes"]["C4_alphaI_rule"]
    assert "not independent replication" in a["changes"]["C3_wording"]["thermostat_sets"]
    assert "not a calibrated negative control" in a["changes"]["C3_wording"]["start_contrast"]
    assert a["unchanged"]["budget"]["cpu_hours_max"] == 2.0 and R0.CFG["min_runs"] == 8 and R0.CFG["min_frames"] == 450


# --- G4b ----------------------------------------------------------------------------------------------------------------

def test_g4b_accepts_the_real_evidence_and_names_each_failure():
    ok = B.g4b_check(AUDIT, A2)
    assert ok["set2_start_direct_evidence"] and ok["a2_outcome"] == "comparable_for_1b" and "does not establish identical production starts" in ok["note"]
    bad = copy.deepcopy(AUDIT)
    bad["questions"]["Q1"]["category"] = "differs_small"
    with pytest.raises(R0.GateError, match="Q1") as e:
        B.g4b_check(bad, A2)
    assert e.value.gate == "G4b"
    bad = copy.deepcopy(AUDIT)
    bad["questions"]["Q4"]["noMyr"]["category"] = "differs_large"
    with pytest.raises(R0.GateError, match="Q4.*noMyr"):
        B.g4b_check(bad, A2)
    bad = copy.deepcopy(AUDIT)
    bad["derived"]["set2_start_direct_evidence"] = False
    with pytest.raises(R0.GateError, match="set2_start_direct_evidence"):
        B.g4b_check(bad, A2)
    bad = copy.deepcopy(A2)
    bad["decision"]["criteria"]["K3_tc_groups"] = False
    with pytest.raises(R0.GateError, match="K1-K5"):
        B.g4b_check(AUDIT, bad)
    bad = copy.deepcopy(A2)
    bad["outcome"] = "not_comparable"
    with pytest.raises(R0.GateError):
        B.g4b_check(AUDIT, bad)
    with pytest.raises(R0.GateError):
        B.g4b_check({}, {})


# --- G1 members and the bilateral G5 rule ----------------------------------------------------------------------------------

def test_wanted_members_drop_set2_system_gro_and_keep_everything_else():
    base = R0.wanted_members
    new = B.make_wanted_members(base)(R0.CFG)
    old = base(R0.CFG)
    for arch in ("kinase_kinked_Myr.zip", "kinase_kinked_noMyr.zip"):
        assert old[arch]["gro:Set2"] == "Set2/system.gro" and new[arch]["gro:Set2"] == "Set2/system_em.gro"
        assert new[arch]["gro:Set1"] == "Set1/system.gro"
        assert {k: v for k, v in new[arch].items() if k != "gro:Set2"} == {k: v for k, v in old[arch].items() if k != "gro:Set2"}
    assert new["PCA.zip"] == old["PCA.zip"]
    assert not any(v.endswith("Set2/system.gro") for m in new.values() for v in m.values())


def times(n, t0=1.0e6):
    return [t0 + 2000.0 * i for i in range(n)]


def run_gate(cfg, usable):
    """Feed the runner's loop order through the bilateral gate; ``usable`` {(set, cond): n usable runs}."""
    gate = B.make_bilateral_gate(cfg, R0.frame_gate)
    out = []
    for s, reps in cfg["sets"].items():
        for c in cfg["conditions"]:
            for k, r in enumerate(reps):
                n = cfg["min_frames"] + 5 if k < usable[(s, c)] else 3
                out.append(gate(times(n), cfg))
    return gate, out


SMALL = {**R0.CFG, "sets": {"Set1": [11, 12, 13, 14, 15], "Set2": [21, 22, 23, 24, 25]}, "min_runs": 4, "min_frames": 30}


def full_usable(**over):
    u = {(s, c): 5 for s in SMALL["sets"] for c in SMALL["conditions"]}
    u.update({tuple(k.split(":")): v for k, v in over.items()})
    return u


def test_bilateral_gate_passes_when_every_group_has_enough_runs():
    gate, out = run_gate(SMALL, full_usable(**{"Set2:noMyr": 4}))
    assert len(out) == 20 and sum(o["excluded"] for o in out) == 1 and gate.state["usable"][("Set2", "noMyr")] == 4
    assert gate.state["excluded"][("Set2", "noMyr")][0]["run"] == "Set2:noMyr:25"


@pytest.mark.parametrize("group", ["Set1:Myr", "Set1:noMyr", "Set2:Myr", "Set2:noMyr"])
def test_bilateral_gate_stops_when_any_one_group_falls_below_the_minimum(group):
    with pytest.raises(R0.GateError, match="bilateral rule") as e:
        run_gate(SMALL, full_usable(**{group: 3}))
    assert e.value.gate == "G5" and e.value.detail["group"] == group and e.value.detail["usable"] == 3
    assert len(e.value.detail["excluded"]) == 2 and "analysis frames" in e.value.detail["excluded"][0]["cause"]


def test_frozen_runner_alone_would_continue_with_one_set():
    """The behaviour the amendment replaces: with Set2 noMyr below the minimum the unmodified rule marks the set not analysed."""
    cfg = SMALL
    runs = {f"{s}:{c}:{r}": 1 for s, reps in cfg["sets"].items() for c in cfg["conditions"] for r in reps if not (s == "Set2" and c == "noMyr" and r > 22)}
    analysed = {s: all(sum(1 for k in runs if k.startswith(f"{s}:{c}:")) >= cfg["min_runs"] for c in cfg["conditions"]) for s in cfg["sets"]}
    assert analysed == {"Set1": True, "Set2": False}


def test_patched_context_restores_the_frozen_runner():
    before = (R0.check_g4_set, R0.wanted_members, R0.frame_gate)
    with B.patched(SMALL, EVIDENCE):
        assert R0.check_g4_set is not before[0] and R0.frame_gate is not before[2]
        assert R0.check_g4_set({}, SMALL)["set2_start_direct_evidence"]
    assert (R0.check_g4_set, R0.wanted_members, R0.frame_gate) == before
    with pytest.raises(R0.GateError):
        with B.patched(SMALL, {"audit": {}, "a2": A2}):
            R0.check_g4_set({}, SMALL)
    assert (R0.check_g4_set, R0.wanted_members, R0.frame_gate) == before


# --- alphaI, production starts, wording --------------------------------------------------------------------------------------

def g3(span, gaps=()):
    return {"alignment": {"reference_span": span, "unobserved_reference_ranges": [list(g) for g in gaps]}}


def test_alpha_i_statement():
    inside = B.alpha_i_statement(g3([250, 540]))
    assert inside["represented"] and inside["not_represented"] == [] and "all alphaI residues" in inside["consequence"]
    cut = B.alpha_i_statement(g3([250, 524]))
    assert not cut["represented"] and cut["not_represented"] == list(range(525, 532)) and "no claim is made about the alphaI helix" in cut["consequence"]
    gap = B.alpha_i_statement(g3([250, 540], [[520, 522]]))
    assert gap["not_represented"] == [520, 521, 522]
    assert not B.alpha_i_statement(g3([250, 500]))["represented"] and not B.alpha_i_statement(g3(None))["represented"]


def synthetic_result(label="TopoFold adds reproduced residues"):
    chans = {c: {"windows_evaluated": 10, "windows_E_positive": 6, "windows_exceeding": 2, "E_median": 0.01} for c in B.CHANNEL_ROLE}
    var = {c: {"within_Myr_js_median_of_windows": 0.1, "within_noMyr_js_median_of_windows": 0.1, "between_js_median_of_windows": 0.11}
           for c in R0.CHANNELS_E}
    temp = {c: {"lag1_median": 0.5, "lag1_max": 0.9, "fraction_decorrelated_within_max_lag": 0.8} for c in R0.CHANNELS_E}
    sets = {s: {"runs": {"Myr": ["a", "b", "c"], "noMyr": ["d", "e", "f"]}, "channels": chans, "run_variability": var, "temporal": temp,
                "windows": []} for s in R0.CFG["sets"]}
    g = {"G0": {"PCA.zip": {"bytes": 1, "sha256": "x"}}, "G1": {"members": 3, "extracted_bytes": 9},
         "G2": {"topology": "ref", "n_atoms": 100, "n_residues": 25, "phi_psi_comparator_available": True},
         "G3": g3([250, 524], [[300, 302]]) | {"alignment": {"identity": 1.0, "reference_span": [250, 524], "unobserved_reference_ranges": [[300, 302]]}},
         "G4": {"Set1": {"note": "G4b note"}}, "G5": {"excluded_runs": [], "sets_analysed": {"Set1": True, "Set2": True}},
         "G6": {"pbc_used": "none", "passes": [{"worst_fraction": 0.01}]}}
    o = {"label": label, "R_T": [1], "R_B": [2], "R_T_outside_N2_R_B": [1], "R_B_outside_N2_R_T": [2], "R_phi_psi": [],
         "reproduced_window_counts": {"kappa": 1}, "jaccard_reproduced_residues": {"a|b": 0.5}, "reproduced_windows": {"kappa": [3]},
         "note": "no channel has a reproduced window"}
    return {"id": "ABL-MYR-1", "gates": g, "sets": sets, "outcome": o, "status": "complete"}


@pytest.mark.parametrize("frozen,step", [("inconclusive", 2), ("TopoFold adds reproduced residues", 3), ("baselines cover TopoFold", 4), ("partial", 5)])
def test_relabel_uses_restricted_wording_without_reproduced(frozen, step):
    res = B.relabel(synthetic_result(frozen), A2)
    o = res["outcome"]
    assert o["label"] == W.LABELS[frozen] and o["decision_rule_step"] == step and "post-hoc" in o["label"]
    assert "reproduced" not in json.dumps(o) and "reproduced" not in o["label"] and "independent" not in o["label"]
    assert "robust_in_both_sets_windows" in o and o["note"] == "no channel has a robust in both sets window"
    assert res["id"] == "ABL-MYR-1b" and res["post_hoc_protocol"] is True
    assert {s: d["leave_one_run_out_recomputations"] for s, d in res["sets"].items()} == {"Set1": 6, "Set2": 6}
    assert res["alpha_I"]["represented"] is False and set(res["production_starts"]) == {"Set1", "Set2"}


def test_production_starts_come_from_the_a2_results_and_the_retained_runs_only():
    keep = {"Set1": {"Myr": ["Set1:Myr:11", "Set1:Myr:12"], "noMyr": ["Set1:noMyr:11"]}, "Set2": {"Myr": [], "noMyr": []}}
    ps = B.production_starts(A2, keep)
    assert set(ps["Set1"]["per_run"]) == {"Set1:Myr:11", "Set1:Myr:12", "Set1:noMyr:11"} and ps["Set2"]["per_run"] == {}
    assert ps["Set1"]["per_run"]["Set1:Myr:11"]["ca_rmsd_to_system_em_A"] == A2["R1_usable_coordinates"]["per_run"]["Set1:Myr:11"]["ca_rmsd_to_em_A"]
    assert ps["Set1"]["geometry_from_a2"] == A2["R2_production_start_geometry"]["Set1"] and "not a calibrated control" in ps["Set1"]["note"]


def test_render_report_states_the_restrictions_and_every_channel_role():
    md = B.render_report(B.relabel(synthetic_result("partial"), A2) | {"lineage_verified": True})
    for needle in ("Post-hoc, descriptive evidence", "consistency check, not independent replication", "No causal, clinical or novelty claim",
                   "stronger-information comparator", "Calpha baseline", "TopoFold", "leave-one-run-out recomputations: 6", "αI helix",
                   "no claim is made about the alphaI helix", "no error-rate control", "Production starts"):
        assert needle in md, needle
    assert "reproduced" not in md.lower().replace("reproduced residues", "")


# --- end to end on SYNTHETIC archives (core + MDAnalysis interpreter) ---------------------------------------------------------

def build_synth(root, early=()):
    """SYNTHETIC archives; runs named in ``early`` get time stamps that end before 1000 ns (excluded at G5)."""
    code = ("import sys, json; sys.path.insert(0, sys.argv[1]); import test_abl_myr_runner as T\n"
            "early = set(json.loads(sys.argv[3])); orig = T.write_xtc\n"
            "T.write_xtc = lambda path, pdb, X, t0_ps=1.0e6, dt_ps=2000.0, times=None: orig(path, pdb, X, t0_ps=0.0 if path.name in early else t0_ps, dt_ps=dt_ps, times=times)\n"
            "spec, data = T.build_dataset(sys.argv[2]); print(json.dumps(spec))")
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(ROOT / "python"), os.environ.get("PYTHONPATH", "")]))
    p = subprocess.run([ADAPTER, "-c", code, str(Path(__file__).parent), str(root), json.dumps(sorted(early))], env=env,
                       capture_output=True, text=True, check=True)
    return json.loads(p.stdout.strip().splitlines()[-1]), Path(root) / "data"


def run_e2e(tmp_path, early=(), evidence=EVIDENCE):
    spec, data = build_synth(tmp_path, early)
    return B.run_wrapped(spec, T.SMALL, data, tmp_path / "out", ADAPTER, evidence, work=tmp_path / "work", log=lambda *_: None)


@e2e
def test_end_to_end_wrapper_on_synthetic_archives(tmp_path):
    res = run_e2e(tmp_path)
    assert res["status"] == "complete" and res["id"] == "ABL-MYR-1b"
    assert res["outcome"]["label"] in set(W.LABELS.values()) and "reproduced" not in json.dumps(res["outcome"])
    assert all(d["leave_one_run_out_recomputations"] == 8 for d in res["sets"].values())          # 4 + 4 retained runs per set
    assert res["gates"]["G4"]["Set2"]["set2_start_direct_evidence"] and "does not read any Set2" in res["gates"]["G4"]["Set2"]["note"]
    assert all("system_em.gro" in v["member"] for k, v in res["gates"]["G1"]["files"].items() if k.endswith("gro:Set2"))
    assert res["alpha_I"]["represented"] is False and res["alpha_I"]["not_represented"] == list(range(518, 532))
    assert set(res["production_starts"]) == {"Set1", "Set2"}
    md = B.render_report(res | {"lineage_verified": True})
    assert "consistency check, not independent replication" in md and "αI helix" in md


@e2e
def test_end_to_end_bilateral_stop_when_one_group_loses_runs(tmp_path):
    early = {f"traj_kinked_noMyr_{r}.xtc" for r in (23, 24)}                                     # Set2 noMyr: 2 of 4 usable (< 3)
    with pytest.raises(R0.GateError, match="bilateral rule") as e:
        run_e2e(tmp_path, early=early)
    assert e.value.gate == "G5" and e.value.detail["group"] == "Set2:noMyr" and e.value.detail["usable"] == 2
    assert {x["run"] for x in e.value.detail["excluded"]} == {"Set2:noMyr:23", "Set2:noMyr:24"}
    assert "0 analysis frames" in e.value.detail["excluded"][0]["cause"]
    assert not (tmp_path / "data/extracted/PCA").exists() or not list((tmp_path / "data/extracted/PCA").rglob("*.xlsx"))


@e2e
def test_end_to_end_g4b_stops_when_the_evidence_does_not_hold(tmp_path):
    bad = copy.deepcopy(AUDIT)
    bad["questions"]["Q1"]["category"] = "differs_small"
    with pytest.raises(R0.GateError, match="Q1") as e:
        run_e2e(tmp_path, evidence={"audit": bad, "a2": A2})
    assert e.value.gate == "G4b"


@e2e
def test_end_to_end_archive_hash_failure_stops_at_g0(tmp_path):
    spec, data = build_synth(tmp_path)
    spec["inputs"]["PCA.zip"]["sha256"] = "0" * 64
    with pytest.raises(R0.GateError, match="SHA-256") as e:
        B.run_wrapped(spec, T.SMALL, data, tmp_path / "out", ADAPTER, EVIDENCE, work=tmp_path / "work", log=lambda *_: None)
    assert e.value.gate == "G0"
