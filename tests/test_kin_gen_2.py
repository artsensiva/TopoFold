"""KIN-GEN-2 library: numerical references, scoring procedure, group naming, leakage controls, mapping control (benchmarks/audit/kin_gen_2).
SYNTHETIC windows with known truth; the real-structure checks use the local KIN-PAIR-1 file and the cached UniProt record and skip without them."""

import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
KG = ROOT / "benchmarks/audit/kin_gen_2"
sys.path.insert(0, str(KG))
import kingen2 as g  # noqa: E402


def helix_window(shift=0.0, noise=0.0, seed=0, radius=2.3, rise=1.5):
    t = np.arange(g.N_WIN) * math.radians(100)
    X = np.c_[radius * np.cos(t), radius * np.sin(t), rise * np.arange(g.N_WIN)] + shift
    return X + np.random.default_rng(seed).normal(0, noise, X.shape)


def test_balanced_accuracy_weighted_ba_and_auroc():
    assert g.balanced_accuracy([1, 1, 0, 0], [1, 0, 0, 0]) == (0.75, 0.5, 1.0)
    rows = [{"kinase": "A", "y": 1, "pred": 1}, {"kinase": "A", "y": 1, "pred": 1}, {"kinase": "B", "y": 1, "pred": 0},        # TPR: A 1.0, B 0.0 -> 0.5 (pooled would be 2/3)
            {"kinase": "A", "y": 0, "pred": 0}, {"kinase": "C", "y": 0, "pred": 0}, {"kinase": "C", "y": 0, "pred": 1}]       # TNR: A 1.0, C 0.5 -> 0.75
    w = g.weighted_ba(rows)
    assert w["tpr"] == 0.5 and w["tnr"] == 0.75 and w["ba_w"] == 0.625 and w["n_kinases_tpr"] == 2 and w["n_kinases_tnr"] == 2
    assert g.auroc([1, 1, 0, 0], [0.9, 0.4, 0.4, 0.1]) == pytest.approx(0.875) and math.isnan(g.auroc([1, 1], [0.1, 0.2]))


def test_group_fit_names_by_rule_and_assignment_uses_the_fitted_centroids_only():
    rng = np.random.default_rng(3)
    pts = np.vstack([rng.normal([6, 14.5], 0.4, (80, 2)), rng.normal([13.5, 11], 0.5, (20, 2)), rng.normal([7.5, 8.5], 0.4, (8, 2))])
    grp = g.fit_groups(pts)
    assert grp["sizes"] == {"in": 80, "out": 20, "inter": 8}
    assert g.assign_group([13.0, 11.5], grp) == "out" and g.assign_group([6.2, 14.0], grp) == "in" and g.assign_group([7.4, 8.4], grp) == "inter"
    assert grp["centroids"]["out"][0] > grp["centroids"]["inter"][0] or grp["centroids"]["out"][0] > grp["centroids"]["in"][0]
    assert max(grp["centroids"], key=lambda k: grp["centroids"][k][0]) == "out"


def test_core_kappa_tau_equal_the_numpy_reference_and_all_representations_are_rigid_motion_invariant():
    rng = np.random.default_rng(5)
    for s in range(5):
        ca = helix_window(noise=0.4, seed=s)
        kap, tau, rk, rt = g.core_kappa_tau(ca)
        assert np.abs(kap - rk).max() < 1e-4 and np.abs(np.arctan2(np.sin(tau - rt), np.cos(tau - rt))).max() < 1e-4
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        Q *= np.sign(np.linalg.det(Q))
        moved = ca @ Q.T + rng.normal(size=3) * 40
        a, b = g.raw_features(ca), g.raw_features(moved)
        assert np.abs(a["T"] - b["T"]).max() < 1e-4                                       # the core receives float32
        assert max(np.abs(a["B1"] - b["B1"]).max(), np.abs(a["B3"] - b["B3"]).max()) < 1e-9 and g.rmsd(ca, moved) < 1e-9
    assert len(g.raw_features(helix_window())["T"]) == 19 and len(g.raw_features(helix_window())["B1"]) == 13 and len(g.raw_features(helix_window())["B3"]) == 36


def test_a_mirror_image_changes_tau_but_not_distances():
    ca = helix_window(noise=0.4, seed=2)
    mirror = ca * np.array([-1.0, 1.0, 1.0])
    a, b = g.raw_features(ca), g.raw_features(mirror)
    assert np.abs(a["B1"] - b["B1"]).max() < 1e-9 and np.abs(a["T"][:7] - b["T"][:7]).max() < 1e-4        # kappa and distances are chirality blind
    assert np.abs(a["T"][7:13] - b["T"][7:13]).max() > 0.1                                                 # sin(tau) flips: the T representation sees handedness


def _toy_model():
    dev, y, kin = [], [], []
    for i in range(6):
        dev.append({"id": f"A{i}", "window_ca": helix_window(noise=0.15, seed=10 + i).tolist()}); y.append("out"); kin.append("KA")      # helical windows, label out
    for i in range(6):
        dev.append({"id": f"B{i}", "window_ca": helix_window(noise=0.15, seed=20 + i, radius=0.4, rise=3.6).tolist()}); y.append("in"); kin.append("KB")   # extended windows, label in
    for i in range(6):
        dev.append({"id": f"C{i}", "window_ca": helix_window(noise=0.15, seed=30 + i, radius=0.4, rise=3.6).tolist()}); y.append("in"); kin.append("KC")
    return g.Model(dev, y), kin


def test_leave_one_kinase_out_never_uses_the_held_out_kinase_as_a_neighbour():
    m, kin = _toy_model()
    ca = m.dev_ca[0]
    assert all(m.score(r, ca, 3) == 1.0 for r in g.REPS)                                                    # with its own kinase present, a helical window is out
    held = np.array(kin) == "KA"
    assert all(m.score(r, ca, 3, exclude=held) == 0.0 for r in g.REPS)                                      # without it, only extended 'in' neighbours remain
    table = g.loko(m, kin, k_grid=(1, 3))
    assert all(table[r][1]["tpr"] == 0.0 for r in g.REPS)                                                    # the held-out out-kinase is mispredicted: no leakage through itself


def test_model_holds_development_entries_only_and_scaling_comes_from_them():
    m, _ = _toy_model()
    assert set(m.ids) == {f"{p}{i}" for p in "ABC" for i in range(6)}
    assert np.allclose(m.F["T"].mean(0), 0, atol=1e-9) and np.allclose(m.F["B1"].std(0), 1, atol=1e-6)


def test_vote_threshold_is_fixed_and_k_nearest_selection_is_deterministic():
    m, _ = _toy_model()
    ca = helix_window(noise=0.15, seed=99)
    assert m.score("T", ca, 5) == m.score("T", ca, 5) and 0.0 <= m.score("B2", ca, 7) <= 1.0


PIM1 = ROOT / "benchmarks/audit/kin_pair_1/data_pdb/8R0Q.pdb"
CACHE = KG / "uniprot_cache/P11309.json"


@pytest.mark.skipif(not (PIM1.is_file() and CACHE.is_file()), reason="local PIM1 structure or cached UniProt record not available")
def test_real_structure_window_d1_d2_and_the_shifted_mapping_control():
    info = g.uniprot_info("P11309", KG / "uniprot_cache")
    assert info["problems"] == [] and info["beta3_lys"] == 67 and info["hrd_hits"] == [165] and info["dfg_hits"] == [186]
    chains = g.parse_chains(PIM1.read_text(encoding="latin-1"))
    r = g.analyse_chain(info, 89, chains["A"])
    assert r["ok"] and r["window_ids"][3] == "A:186" and r["window_resnames"][3:6] == ["ASP", "PHE", "GLY"] and len(r["window_ca"]) == 9
    res = {(x["resid"]): x for x in chains["A"]}
    d1 = float(np.linalg.norm(g.xyz(res[93], "CA") - g.xyz(res[187], "CZ")))
    d2 = float(np.linalg.norm(g.xyz(res[67], "CA") - g.xyz(res[187], "CZ")))
    assert r["D1"] == pytest.approx(d1, abs=1e-9) and r["D2"] == pytest.approx(d2, abs=1e-9)                    # independent recomputation from the atoms
    shifted = dict(info, beta3_lys=68)
    assert not g.analyse_chain(shifted, 89, chains["A"])["ok"]                                                 # a mapping shifted by one residue is refused by the helper


def test_freeze_check_accepts_appended_amendments_only():
    sh = subprocess.run(["git", "rev-parse", "--is-shallow-repository"], capture_output=True, text=True, cwd=ROOT).stdout.strip()
    fr = KG / "results/KIN-GEN-2_freeze_1.json"
    if sh == "true" or not fr.is_file():
        pytest.skip("needs the full git history")
    assert g.check_freeze(fr) == []
