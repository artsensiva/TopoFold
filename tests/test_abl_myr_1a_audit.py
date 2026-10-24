"""SYNTHETIC tests of the ABL-MYR-1a initialization audit (benchmarks/audit/abl_myr_1a/run_init_audit.py). Engineering checks
only: the archives are generated here, and nothing is read from Dryad. The tpr worker test runs only where MDAnalysis is
installed."""

import hashlib
import json
import sys
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_myr_1a"))
sys.path.insert(0, str(ROOT / "benchmarks/audit/abl_myr_1"))
import run_init_audit as A  # noqa: E402
from test_abl_myr_runner import MDP, base_backbone, gro_text, residues, top_text  # noqa: E402

SPEC = json.loads((ROOT / "benchmarks/audit/abl_myr_1a/audit_spec.json").read_text())
BUDGET = {"budget": {"extracted_bytes_max": 10 ** 9, "cpu_seconds_max": 600}}
NO_TPR = lambda files: {k: {"readable": False, "error": "SYNTHETIC unreadable"} for k in files}  # noqa: E731


def build(tmp_path, *, pre=None, em=None, mdp=None, sol=None, drop=(), bad_gro=False, same_em_as_set1=False):
    """SYNTHETIC model archives. ``pre``/``em``: {(set, cond): shift in nm}; ``mdp``: {(set, cond): text};
    ``sol``: {(set, cond): water count}; ``drop``: member suffixes to leave out."""
    pre, em, mdp, sol = pre or {}, em or {}, mdp or {}, sol or {}
    res, X0 = residues(), base_backbone()
    stage = tmp_path / "stage"
    for c, model in A.CONDITIONS.items():
        for s, reps in A.SETS.items():
            d = stage / model / s
            d.mkdir(parents=True, exist_ok=True)
            e_key = ("Set1", c) if same_em_as_set1 and s == "Set2" else (s, c)
            if s == "Set1":
                (d / "system.gro").write_text(gro_text(res, X0, c == "Myr", pre.get((s, c), 0.0)))
            (d / "system_em.gro").write_text(gro_text(res, X0, c == "Myr", em.get(e_key, 0.0)))
            (d / "production_md.mdp").write_text(mdp.get((s, c), MDP))
            (d / "topol.top").write_text(top_text(c == "Myr", sol.get((s, c), 5000)))
            for r in reps:
                (d / str(r)).mkdir(exist_ok=True)
                (d / str(r) / "md.tpr").write_bytes(b"SYNTHETIC tpr")
    if bad_gro:
        (stage / "kinase_kinked_Myr/Set2/system_em.gro").write_text("SYNTHETIC\nnot-a-number\n")
    data = tmp_path / "data"
    data.mkdir()
    orig = {"inputs": {}}
    for model in A.CONDITIONS.values():
        with zipfile.ZipFile(data / f"{model}.zip", "w") as zf:
            for p in sorted((stage / model).rglob("*")):
                rel = str(p.relative_to(stage))
                if p.is_file() and not any(rel.endswith(x) for x in drop):
                    zf.write(p, f"{model}/{rel[len(model) + 1:]}")
        b = (data / f"{model}.zip").read_bytes()
        orig["inputs"][f"{model}.zip"] = {"bytes": len(b), "sha256": hashlib.sha256(b).hexdigest()}
    return orig, data


def audit(tmp_path, spec=BUDGET, **kw):
    orig, data = build(tmp_path, **kw)
    return A.run_audit(orig, data, tmp_path / "ext", NO_TPR, {"id": "X", **spec}, log=lambda *_: None)


def test_spec_constants_and_rules_match_the_code():
    assert (A.IDENTICAL_NM, A.SMALL_RMSD_NM) == (SPEC["constants"]["identical_nm"], SPEC["constants"]["small_rmsd_nm"])
    assert A.IDENTICAL_NM == A.R0.CFG["gro_tol_nm"]
    assert [r["verdict"] for r in SPEC["verdict_rules_ordered"]] == [
        "audit_incomplete", "V4_set1_comparability_failed", "V1_set2_start_established_direct", "V3_set2_not_comparable",
        "V2_set2_start_not_established"]
    assert set(SPEC["what_each_file_can_and_cannot_establish"]) >= {"Set1/system.gro", "<rep>/md.tpr (all 40 runs)"}
    assert A.SETS == {"Set1": list(range(11, 21)), "Set2": list(range(21, 31))}
    plan = A.members_plan()["kinase_kinked_Myr.zip"]
    assert "system.gro:Set2" not in plan and "system.gro:Set1" in plan and "system_em.gro:Set2" in plan
    assert sum(k.startswith("tpr:") for k in plan) == 20


def test_coordinate_categories_and_boundaries(tmp_path):
    res, X0 = residues(), base_backbone()

    def cmp(shift):
        (tmp_path / "a.gro").write_text(gro_text(res, X0, False, 0.0))
        (tmp_path / "b.gro").write_text(gro_text(res, X0, False, shift))
        return A.compare_gro_files(tmp_path / "a.gro", tmp_path / "b.gro")
    assert cmp(0.0)["category"] == "identical" and cmp(0.0)["max_abs_dx_nm"] == 0.0
    assert cmp(0.001)["category"] == "identical"                                      # exactly one .gro digit: inside the tolerance
    r = cmp(0.002)
    assert r["category"] == "differs_small" and r["atoms_over_0.001nm"] == r["n_protein_atoms"]
    assert cmp(0.0288)["category"] == "differs_small"                                  # RMSD 0.0499
    assert cmp(0.0292)["category"] == "differs_large"                                  # RMSD 0.0506
    g = (tmp_path / "b.gro").read_text().splitlines()
    (tmp_path / "c.gro").write_text("\n".join([g[0], str(int(g[1]) - 1)] + g[3:]) + "\n")
    assert A.compare_gro_files(tmp_path / "a.gro", tmp_path / "c.gro")["category"] == "sequence_mismatch"
    h = ["SYN", "1", "    1SOL     OW    1   0.000   0.000   0.000", "   1.0 1.0 1.0"]
    (tmp_path / "w.gro").write_text("\n".join(h) + "\n")
    assert A.compare_gro_files(tmp_path / "w.gro", tmp_path / "w.gro")["category"] == "sequence_mismatch"   # no protein atoms


def test_mdp_and_top_comparisons_name_their_differences():
    assert A.compare_mdp_text(MDP, MDP + "; comment\n\n")["identical"]
    r = A.compare_mdp_text(MDP, MDP.replace("gen_vel = no", "gen_vel = yes"))
    assert not r["identical"] and any("gen_vel" in x for x in r["differing_lines"])
    assert A.mdp_key_differences(MDP, MDP.replace("ref_t = 300", "ref_t = 310")) == {"ref_t": ["300", "310"]}
    assert A.compare_top_text(top_text(True), top_text(False))["identical_after_myr_removal"]
    t = A.compare_top_text(top_text(True), top_text(False, sol=4990))
    assert not t["identical_after_myr_removal"] and any("SOL" in x for x in t["differing_lines"])


def q(q1="identical", q2="identical", q4="identical", mdp1=True, top1=True, mdp2=True, top2=True):
    c = lambda x: {"category": x}  # noqa: E731
    return {"Q1": c(q1), "Q2": c(q2), "Q4": {"Myr": c(q4), "noMyr": c(q4)},
            "Q5": {"Set1": {"identical": mdp1}, "Set2": {"identical": mdp2}},
            "Q6": {"Set1": {"identical_after_myr_removal": top1}, "Set2": {"identical_after_myr_removal": top2}}}


@pytest.mark.parametrize("kw,verdict", [
    (dict(), "V1_set2_start_established_direct"),
    (dict(q2="differs_small", q4="identical"), "V1_set2_start_established_direct"),                  # copy of Set1 structures
    (dict(q2="differs_small", q4="differs_small"), "V2_set2_start_not_established"),
    (dict(q2="differs_large", q4="differs_small"), "V3_set2_not_comparable"),
    (dict(q2="sequence_mismatch", q4="sequence_mismatch"), "V3_set2_not_comparable"),
    (dict(q2="differs_large", q4="identical"), "V1_set2_start_established_direct"),                  # direct evidence wins
    (dict(q1="differs_small"), "V4_set1_comparability_failed"),
    (dict(q1="differs_small", q2="identical"), "V4_set1_comparability_failed"),
    (dict(mdp1=False), "V4_set1_comparability_failed"),
    (dict(top1=False), "V4_set1_comparability_failed"),
    (dict(mdp2=False), "V2_set2_start_not_established"),                                              # direct evidence, parameters differ
    (dict(top2=False, q2="differs_small", q4="differs_small"), "V2_set2_start_not_established"),
    (dict(q2="identical", q4="differs_small", q1="identical", top2=False), "V2_set2_start_not_established"),
])
def test_decision_table(kw, verdict):
    assert A.decide(q(**kw))["verdict"] == verdict


def test_q4_route_requires_q1_identical():
    d = A.decide(q(q1="differs_small", q2="differs_small", q4="identical"))
    assert d["verdict"] == "V4_set1_comparability_failed" and d["set2_start_direct_evidence"] is False


def test_tpr_summary_from_probe_records():
    probe = {}
    for s, reps in A.SETS.items():
        for c in A.CONDITIONS:
            for r in reps:
                probe[f"{s}:{c}:{r}"] = {"readable": True, "n_atoms": 10, "protein_digest": "d1", "positions_available": False}
    probe["Set2:Myr:21"] = {"readable": False, "error": "boom"}
    t = A.tpr_summary(probe)
    assert t["Set1"]["Myr"]["readable"] == 10 and t["Set1"]["Myr"]["positions_available_in_any"] is False
    assert t["Set2"]["Myr"]["readable"] == 9 and t["Set2"]["Myr"]["unreadable"] == 1 and t["Set2"]["Myr"]["errors"] == ["boom"]
    assert t["Set1"]["protein_digest_sets_equal_between_conditions"] is True


# --- run_audit on SYNTHETIC archives: one test per verdict and stop ------------------------------------------------------

def test_audit_v1_when_set2_structures_are_identical(tmp_path):
    r = audit(tmp_path)
    assert r["verdict"] == "V1_set2_start_established_direct" and r["questions"]["Q2"]["category"] == "identical"
    assert r["extracted"]["bytes"] > 0 and r["resources"]["download_bytes"] == 0
    assert r["questions"]["Q7"]["Set2"]["Myr"]["unreadable"] == 10
    assert set(r["questions"]["Q8"]["Set1"]["Myr"]) >= {"gen_vel", "continuation"}


def test_audit_v1_through_the_cross_set_copy_route(tmp_path):
    em = {("Set1", "noMyr"): 0.01}
    r = audit(tmp_path, em=em, same_em_as_set1=True)
    assert r["questions"]["Q2"]["category"] == "differs_small" and r["questions"]["Q3"]["category"] == "differs_small"
    assert all(r["questions"]["Q4"][c]["category"] == "identical" for c in ("Myr", "noMyr"))
    assert r["verdict"] == "V1_set2_start_established_direct"


def test_audit_v2_when_set2_differs_slightly_and_is_not_a_copy(tmp_path):
    r = audit(tmp_path, em={("Set1", "noMyr"): 0.01, ("Set2", "noMyr"): 0.02})
    assert r["verdict"] == "V2_set2_start_not_established" and r["derived"]["set2_start_direct_evidence"] is False


def test_audit_v3_when_set2_differs_greatly(tmp_path):
    r = audit(tmp_path, em={("Set2", "noMyr"): 0.5})
    assert r["verdict"] == "V3_set2_not_comparable" and r["questions"]["Q2"]["category"] == "differs_large"


def test_audit_v4_when_set1_start_differs_or_set1_topology_differs(tmp_path):
    r = audit(tmp_path / "a", pre={("Set1", "noMyr"): 0.002})
    assert r["verdict"] == "V4_set1_comparability_failed" and r["questions"]["Q1"]["category"] == "differs_small"
    r = audit(tmp_path / "b", sol={("Set1", "noMyr"): 4990})
    assert r["verdict"] == "V4_set1_comparability_failed" and not r["questions"]["Q6"]["Set1"]["identical_after_myr_removal"]
    assert r["questions"]["Q6"]["Set2"]["identical_after_myr_removal"]


def test_audit_stops_on_each_stop_condition(tmp_path):
    orig, data = build(tmp_path / "s1")
    orig["inputs"]["kinase_kinked_Myr.zip"]["sha256"] = "0" * 64
    r = A.run_audit(orig, data, tmp_path / "x1", NO_TPR, {"id": "X", **BUDGET})
    assert r["verdict"] == "audit_incomplete" and r["stop"]["code"] == "S1" and "extracted" not in r
    r = audit(tmp_path / "s2", drop=("Set2/topol.top",))
    assert r["stop"]["code"] == "S2" and "0 matches" in r["stop"]["cause"]
    r = audit(tmp_path / "s3", bad_gro=True)
    assert r["stop"]["code"] == "S3" and "cannot parse" in r["stop"]["cause"]
    r = audit(tmp_path / "s4", spec={"budget": {"extracted_bytes_max": 100, "cpu_seconds_max": 600}})
    assert r["stop"]["code"] == "S4"
    assert "verdict" in r and r["verdict"] == "audit_incomplete"


def test_report_is_generated_from_the_results_json(tmp_path):
    r = audit(tmp_path)
    md = A.render_report(r)
    assert "V1_set2_start_established_direct" in md and "Q2 Set2 system_em.gro" in md and "not evidence about" in md.lower()
    assert "audit_incomplete" in A.render_report({"verdict": "audit_incomplete", "stop": {"code": "S1", "cause": "x"}})


def test_original_artifacts_are_preserved_when_a_freeze_exists():
    orig_cf = json.loads((ROOT / "benchmarks/audit/abl_myr_1/results/ABL-MYR-1_code_freeze.json").read_text())
    assert A.sha256_file(ROOT / "benchmarks/audit/abl_myr_1/run_abl_myr.py") == orig_cf["runner"]["sha256"]
    frz = ROOT / "benchmarks/audit/abl_myr_1a/results/ABL-MYR-1a_audit_freeze.json"
    if not frz.exists():
        pytest.skip("audit freeze not yet written")
    rec = json.loads(frz.read_text())
    for p, h in rec["original_artifacts_sha256"].items():
        assert A.sha256_file(ROOT / p) == h, p


def test_tpr_worker_records_unreadable_files_without_stopping(tmp_path):
    pytest.importorskip("MDAnalysis")
    bad = tmp_path / "bad.tpr"
    bad.write_bytes(b"not a tpr")
    out = A._w_tpr({"files": {"Set1:Myr:11": str(bad), "Set1:Myr:12": str(tmp_path / "missing.tpr")}})
    assert all(v["readable"] is False and v["error"] for v in out.values())
