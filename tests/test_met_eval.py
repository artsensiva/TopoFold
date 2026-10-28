"""MET exploratory evaluation: synthetic tests of joins, endpoint arithmetic, eligibility, identical rows, weighting, buffering,
fold-local preprocessing, failure handling, resource accounting and the decision rules. SYNTHETIC data only; no scientific result."""
import copy
import importlib.util
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
D = ROOT / "benchmarks/audit/met_eval_2026-10-28"
spec = importlib.util.spec_from_file_location("met_evaluate", D / "evaluate.py")
ev = importlib.util.module_from_spec(spec)
sys.modules["met_evaluate"] = ev
spec.loader.exec_module(ev)
pytest.importorskip("sklearn")
pytest.importorskip("scipy")

HEADER = "variants\tposition\twildtype\tmutation\ttype\tROSACE_effects\tROSACE_tests\tROSACE_alt\tROSACE.GOF\tROSACE.LOF\tinhibitor\n"


def tsv(tmp_path, rows):
    p = tmp_path / "assay.tsv"
    p.write_text(HEADER + "".join(f"v\t{pos}\tA\t{mut}\t{typ}\t{eff}\t0\tFALSE\tFALSE\t{'TRUE' if lof else 'FALSE'}\t{inh}\n" for inh, pos, mut, typ, eff, lof in rows))
    return p


# ------------------------------------------------------------------------------------------------ endpoint, eligibility, joins


def test_endpoint_arithmetic_and_eligibility(tmp_path):
    rows = [("Crizo", 1, "C", "missense", 1.0, False), ("DMSO", 1, "C", "missense", 0.25, False),          # gamma = 0.75
            ("Crizo", 2, "C", "missense", 3.0, False), ("DMSO", 2, "C", "missense", 1.0, True),            # DMSO LOF: excluded
            ("Crizo", 3, "C", "missense", 3.0, False),                                                      # no DMSO row: excluded
            ("Crizo", 4, "C", "synonymous", 3.0, False), ("DMSO", 4, "C", "synonymous", 1.0, False),       # not missense: excluded
            ("Crizo", 5, "C", "missense", "NA", False), ("DMSO", 5, "C", "missense", 1.0, False),          # non-numeric: excluded
            ("Tiv", 6, "C", "missense", 3.0, False), ("DMSO", 6, "C", "missense", 1.0, False),             # other inhibitor: ignored
            ("Crizo", 7, "D", "missense", -2.0, False), ("DMSO", 7, "D", "missense", -0.5, False)]
    a = ev.read_assay(tsv(tmp_path, rows), ["Crizo"])
    g = ev.gamma_table(a, "Crizo")
    assert set(g) == {(1, "C"), (7, "D")} and g[(1, "C")][0] == pytest.approx(0.75) and g[(7, "D")][0] == pytest.approx(-1.5)


def test_duplicate_assay_rows_raise(tmp_path):
    rows = [("Crizo", 1, "C", "missense", 1.0, False), ("Crizo", 1, "C", "missense", 2.0, False)]
    with pytest.raises(ValueError):
        ev.read_assay(tsv(tmp_path, rows), ["Crizo"])


def synthetic_lookup(n_pos=20, gly_pos=5, bad_structure=(7,), no_esm=((3, "D"),)):
    struct, mutation, esm, gamma = {}, {}, {}, {}
    rng = np.random.default_rng(0)
    for p in range(1, n_pos + 1):
        row = {"T_reason": "ok", "L_reason": "no_side_chain" if p == gly_pos else "ok", "L_min_heavy_dist_A": "5.0", "L_min_sidechain_dist_A": "" if p == gly_pos else "7.0", "L_n_contact_atoms_4p5": "2", "L_ca_min_dist_A": "6.0"}
        for c in ev.T_COLS + ev.C_COLS:
            row[c] = f"{rng.normal():.4f}"
        if p in bad_structure:
            row["T_reason"] = "missing_residue:1"
            for c in ev.T_COLS + ev.C_COLS:
                row[c] = ""
        struct[("Crizo", p)] = row
    for wt in "ACD":
        for mut in "ACD":
            if wt != mut:
                mutation[(wt, mut)] = {c: str(rng.normal()) for c in ev.M_COLS}
    for p in range(1, n_pos + 1):
        for mut in "CD":
            gamma[(p, mut)] = (float(rng.normal()), "A")
            if (p, mut) not in no_esm:
                esm[(p, mut)] = float(rng.normal())
    return struct, mutation, esm, gamma


def test_join_keeps_identical_rows_for_every_feature_set_and_counts_exclusions_explicitly():
    struct, mutation, esm, gamma = synthetic_lookup()
    a = ev.build_drug_arrays("Crizo", gamma, struct, mutation, esm)
    c = a["counts"]
    assert c["eligible"] == 40 and c["no_structure"] == 2 and c["no_esm"] == 1 and c["used"] == 37
    assert ev.feature_matrix(a, ["M"]).shape[0] == ev.feature_matrix(a, ["M", "L", "C"]).shape[0] == ev.feature_matrix(a, ["M", "L", "ESM", "T"]).shape[0] == 37
    assert 7 not in a["pos"] and (3, "D") not in set(zip(a["pos"], a["mut"]))
    assert a["row_id"][0] == f"Crizo:{a['pos'][0]}{a['mut'][0]}"


def test_glycine_distance_replacement_is_deterministic_and_flagged():
    struct, mutation, esm, gamma = synthetic_lookup()
    a = ev.build_drug_arrays("Crizo", gamma, struct, mutation, esm)
    L = a["X"]["L"]
    i = int(np.where(a["pos"] == 5)[0][0])
    assert L[i, 1] == L[i, 0] == 5.0 and L[i, 4] == 1.0
    j = int(np.where(a["pos"] == 6)[0][0])
    assert L[j, 1] == 7.0 and L[j, 4] == 0.0


def test_esm_reader_uses_only_two_columns_and_rejects_non_invariant_scores(tmp_path):
    p = tmp_path / "esm.csv"
    p.write_text("pos_mut,score,key,mean,pos,resistance\n1K,-8.2,Crizo,0.5,1,0\n1K,-8.2,Camp,9.9,1,1\n2R,-1.0,Crizo,0.1,2,0\n")
    assert ev.read_esm(p) == {(1, "K"): -8.2, (2, "R"): -1.0}
    p.write_text("pos_mut,score,key,mean,pos,resistance\n1K,-8.2,Crizo,0.5,1,0\n1K,-7.0,Camp,9.9,1,1\n")
    with pytest.raises(ValueError):
        ev.read_esm(p)


def test_no_feature_block_name_refers_to_an_outcome():
    names = [c for b in ev.BLOCK_COLS.values() for c in b]
    assert not [n for n in names if any(t in n.lower() for t in ("gamma", "effect", "fitness", "resist", "score", "rosace"))]
    assert len(ev.M_COLS) == 14 and len(ev.L_COLS) == 5 and len(ev.T_COLS) == 19 and len(ev.C_COLS) == 36


# ------------------------------------------------------------------------------------------------ folds, buffer, weighting


def test_blocks_cover_positions_once_and_training_respects_the_buffer():
    f = ev.block_assignment(287, 10)
    assert len(f) == 287 and sorted(np.bincount(f).tolist()) == [28, 28, 28, 29, 29, 29, 29, 29, 29, 29]
    assert all(f[i] <= f[i + 1] for i in range(286))
    pos = np.tile(np.arange(1, 288), 3)  # three "drugs" sharing positions: a position is never split across folds
    for k in range(10):
        train, test = ev.train_test_masks(pos, f, k, 8)
        held = np.where(f == k)[0] + 1
        assert (f[pos[test] - 1] == k).all() and not (train & test).any()
        assert np.abs(pos[train][:, None] - held[None]).min() > 8
        assert set(pos[train]) | set(range(held.min() - 8, held.max() + 9)) >= set(range(1, 288))


def test_position_weighting_is_not_row_weighting():
    pos, e2 = np.array([1, 1, 1, 2]), np.array([1.0, 1.0, 1.0, 9.0])
    assert ev.position_weighted_mse(pos, e2) == pytest.approx(5.0) and e2.mean() == pytest.approx(3.0)
    assert ev.position_mean_sq(pos, e2) == {1: 1.0, 2: 9.0}


# ------------------------------------------------------------------------------------------------ fold-local preprocessing, alpha


def synth_arrays(n_pos=60, n_sub=10, seed=1, signal=None, noise=0.3):
    rng = np.random.default_rng(seed)
    pos = np.repeat(np.arange(1, n_pos + 1), n_sub)
    n = len(pos)
    X = {"M": rng.normal(size=(n, 14)), "L": rng.normal(size=(n, 5)), "T": rng.normal(size=(n, 19)), "C": rng.normal(size=(n, 36)), "ESM": rng.normal(size=(n, 1))}
    # position-level structural blocks are constant within a position (as in the exporter)
    for b in ("L", "T", "C"):
        base = rng.normal(size=(n_pos, X[b].shape[1]))
        X[b] = base[pos - 1]
    y = rng.normal(scale=noise, size=n)
    if signal == "L":
        y = y + 2.0 * X["L"][:, 0]
    return {"pos": pos, "mut": ["C"] * n, "y": y, "sd": np.full(n, 0.4), "X": X, "row_id": [f"S:{p}C{i}" for i, p in enumerate(pos)], "counts": {}}


def test_standardization_uses_training_statistics_only():
    Xtr = np.array([[1.0, 5.0], [3.0, 5.0], [5.0, 5.0]])
    Xte = np.array([[1000.0, 5.0]])
    A, B = ev.standardize(Xtr, Xte)
    assert A.mean(0)[0] == pytest.approx(0) and A.std(0)[0] == pytest.approx(1) and (A[:, 1] == 0).all() and B[0, 1] == 0
    A2, _ = ev.standardize(Xtr, np.array([[-1e9, 5.0]]))
    np.testing.assert_allclose(A, A2)


def test_predictions_do_not_depend_on_held_out_or_buffer_outcomes():
    a = synth_arrays()
    f = ev.block_assignment(60, 10)
    p0, info, test = ev.run_unit(a, f, 4)
    assert info["flag"] == "ok"
    b = copy.deepcopy(a)
    held = np.where(f == 4)[0] + 1
    zone = (a["pos"] >= held.min() - 8) & (a["pos"] <= held.max() + 8)  # held-out block and its buffer
    b["y"] = np.where(zone, b["y"] + 1000.0, b["y"])
    p1, _, _ = ev.run_unit(b, f, 4)
    for k in p0:
        np.testing.assert_array_equal(p0[k], p1[k])
    c = copy.deepcopy(a)
    far = ~zone
    c["y"] = np.where(far, c["y"] + 50.0, c["y"])  # training outcomes do matter
    p2, _, _ = ev.run_unit(c, f, 4)
    assert not np.allclose(p0[("F0", "ridge")], p2[("F0", "ridge")])


def test_alpha_selection_default_tie_rule_and_position_weighted_loss(monkeypatch):
    a = synth_arrays(n_pos=6, n_sub=2)
    train = np.ones(len(a["pos"]), bool)
    alpha, losses, used = ev.select_alpha(ev.feature_matrix(a, ["M"]), a["y"], a["pos"], train)
    assert (alpha, losses, used) == (ev.DESIGN["default_alpha_if_no_inner_fold"], None, 0)  # buffer leaves no usable inner training rows
    b = synth_arrays(n_pos=60, n_sub=5)
    X0 = np.zeros((len(b["pos"]), 3))  # uninformative features: all alphas tie -> the larger alpha
    alpha, losses, used = ev.select_alpha(X0, b["y"], b["pos"], np.ones(len(b["pos"]), bool))
    assert used == ev.DESIGN["inner_blocks"] and alpha == max(ev.DESIGN["alphas"])
    calls = []
    real = ev.position_weighted_mse
    monkeypatch.setattr(ev, "position_weighted_mse", lambda p, e: calls.append(len(p)) or real(p, e))
    ev.select_alpha(ev.feature_matrix(b, ["M"]), b["y"], b["pos"], np.ones(len(b["pos"]), bool))
    assert len(calls) == ev.DESIGN["inner_blocks"] * len(ev.DESIGN["alphas"])


def test_inner_folds_respect_position_grouping_and_buffer():
    pos = np.repeat(np.arange(1, 61), 3)
    for g in ev.inner_blocks(pos, 5):
        far = (pos < g.min() - 8) | (pos > g.max() + 8)
        assert not (np.isin(pos, g) & far).any() and np.isin(pos, g).sum() == 3 * len(g)


# ------------------------------------------------------------------------------------------------ failure handling


def test_empty_insufficient_and_untrainable_cells_are_flagged_not_dropped():
    a = synth_arrays()
    f = ev.block_assignment(60, 10)
    empty = copy.deepcopy(a)
    keep = f[empty["pos"] - 1] != 3
    for k in ("pos", "y", "sd"):
        empty[k] = empty[k][keep]
    empty["X"] = {b: v[keep] for b, v in empty["X"].items()}
    empty["mut"], empty["row_id"] = [m for m, kp in zip(a["mut"], keep) if kp], [r for r, kp in zip(a["row_id"], keep) if kp]
    p, info, test = ev.run_unit(empty, f, 3)
    assert info["flag"] == "empty_test" and p == {}
    few = synth_arrays(n_pos=60, n_sub=1)  # 1 row per position -> 60 rows; training for a block drops 6+16 -> 38 rows < 50
    p, info, test = ev.run_unit(few, f, 4)
    assert info["flag"] == "insufficient_training_fallback_mean"
    vals = {tuple(np.round(v, 9)) for v in p.values()}
    assert len(vals) == 1 and len(p) == len(ev.DESIGN["feature_sets"]) * 2  # the same fallback for every feature set and model
    tiny = synth_arrays(n_pos=12, n_sub=3)
    p, info, test = ev.run_unit(tiny, ev.block_assignment(12, 10), 0)
    assert info["flag"] in ("no_training_rows", "insufficient_training_fallback_mean", "ok")


def test_run_units_resumes_without_redoing_and_stops_at_the_budget(tmp_path, monkeypatch):
    data = {"S1": synth_arrays(seed=2), "S2": synth_arrays(seed=3)}
    f = ev.block_assignment(60, 10)
    prog = tmp_path / "progress.jsonl"
    calls = []
    real = ev.run_unit
    monkeypatch.setattr(ev, "run_unit", lambda a, ff, k: calls.append(k) or real(a, ff, k))
    done, status = ev.run_units(data, f, prog, cpu_cap=0.0, rss_cap=1e9)
    assert status == "budget exhausted: partial" and done == {} and calls == []
    monkeypatch.setattr(ev, "DESIGN", dict(ev.DESIGN, outer_blocks=10))
    done, status = ev.run_units({"S1": data["S1"]}, f, prog, cpu_cap=1e9, rss_cap=1e9)
    assert status == "complete" and len(done) == 10 and len(calls) == 10
    done2, status2 = ev.run_units({"S1": data["S1"]}, f, prog, cpu_cap=1e9, rss_cap=1e9)
    assert status2 == "complete" and len(calls) == 10 and len(done2) == 10  # nothing recomputed
    assert len(prog.read_text().strip().splitlines()) == 10


def test_cpu_accounting_includes_children_and_threads_are_bounded():
    import os
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        assert os.environ[k] == "1"
    t0 = ev.cpu_seconds()
    subprocess.run([sys.executable, "-c", "import time\nt=time.process_time()\nwhile time.process_time()-t<0.4: sum(range(1000))"], check=True)
    assert ev.cpu_seconds() - t0 >= 0.3


# ------------------------------------------------------------------------------------------------ bootstrap, rules, end to end


def test_bootstrap_is_reproducible_jointly_resampled_and_shift_equivariant():
    rng = np.random.default_rng(0)
    Eb = rng.uniform(0.5, 1.5, size=(6, 60))
    cfg = dict(ev.DESIGN["bootstrap"], resamples=300)
    d1 = ev.bootstrap_deltas(Eb - 0.25, Eb, 60, cfg)
    d2 = ev.bootstrap_deltas(Eb - 0.25, Eb, 60, cfg)
    np.testing.assert_array_equal(d1, d2)
    np.testing.assert_allclose(d1, -0.25)
    np.testing.assert_allclose(ev.bootstrap_deltas(Eb, Eb, 60, cfg), 0.0)


def test_decision_rules_truth_table():
    c = ev.classify
    assert c(-0.1, -0.02, 9, 6) == "criterion met" and c(-0.1, -0.02, 8, 5) == "criterion met"
    assert c(-0.1, 0.01, 9, 6) == "inconclusive" and c(-0.1, -0.02, 7, 6) == "inconclusive" and c(-0.1, -0.02, 9, 4) == "inconclusive"
    assert c(0.02, 0.1, 9, 6) == "criterion not met" and c(-0.1, -0.02, 5, 6) == "criterion not met"
    assert ev.family_outcome("criterion met", "criterion met") == "criterion met in both families"
    assert "one family only" in ev.family_outcome("criterion met", "inconclusive") and "either family" in ev.family_outcome("inconclusive", "criterion not met")

    def res(rel3, cls3, cls4):
        return {"C3": {f: {"rel": rel3[i], "cls": cls3[i]} for i, f in enumerate(("ridge", "gbm"))}, "C4": {f: {"rel": 0.1, "cls": cls4[i]} for i, f in enumerate(("ridge", "gbm"))}}
    met = "criterion met"
    assert ev.engineering_outcome(res([0.06, 0.06], [met, met], [met, met])).startswith("continue; angular")
    assert ev.engineering_outcome(res([0.06, 0.06], [met, met], [met, "inconclusive"])).startswith("continue with geometry")
    assert ev.engineering_outcome(res([0.049, 0.06], [met, met], [met, met])).startswith("do not continue")  # below the 5% threshold in one family
    assert ev.engineering_outcome(res([0.06, 0.06], [met, "inconclusive"], [met, met])).startswith("do not continue")  # one family disagrees


def test_end_to_end_recovers_a_planted_signal_and_finds_none_in_noise(monkeypatch):
    monkeypatch.setattr(ev, "DESIGN", dict(ev.DESIGN, bootstrap=dict(ev.DESIGN["bootstrap"], resamples=300)))
    f = ev.block_assignment(60, 10)
    data = {f"S{i}": synth_arrays(seed=10 + i, signal="L") for i in range(6)}
    preds = {d: {} for d in data}
    for d, a in data.items():
        for k in range(10):
            p, info, test = ev.run_unit(a, f, k)
            assert info["flag"] == "ok"
            for key, v in p.items():
                preds[d].setdefault(key, np.full(len(a["y"]), np.nan))[test] = v
    res = ev.evaluate_predictions(data, preds, f, 60)
    c1, c3, c4 = res["comparisons"]["C1"], res["comparisons"]["C3"], res["comparisons"]["C4"]
    assert c1["families"]["ridge"]["classification"] == "criterion met" and c1["families"]["ridge"]["relative_mse_reduction"] > 0.3
    assert c3["families"]["ridge"]["classification"] != "criterion met"  # T is pure noise
    assert res["engineering_outcome"].startswith("do not continue")
    assert math.isfinite(c4["families"]["gbm"]["delta"])


# ------------------------------------------------------------------------------------------------ freeze


def test_freeze_record_matches_once_it_exists():
    fz = D / "freeze_eval.json"
    if not fz.is_file():
        pytest.skip("not frozen yet")
    rec = json.loads(fz.read_text())
    for f, h in rec["files"].items():
        if f.endswith("test_met_eval.py"):
            continue  # checked by `evaluate.py run`
        assert ev.sha256_file(D / f) == h, f
