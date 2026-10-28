"""Consistency of the committed MET evaluation records: stored classifications and the engineering outcome follow from the stored numbers by the frozen rules."""
import gzip
import importlib.util
import json
import sys
from pathlib import Path

D = Path(__file__).resolve().parent.parent / "benchmarks/audit/met_eval_2026-10-28"
spec = importlib.util.spec_from_file_location("met_evaluate_r", D / "evaluate.py")
ev = importlib.util.module_from_spec(spec)
sys.modules["met_evaluate_r"] = ev
spec.loader.exec_module(ev)
R = json.loads((D / "results/results.json").read_text())


def test_run_is_complete_with_every_cell_present_and_the_freeze_recorded():
    assert R["status"] == "complete" and len(R["cells"]) == 60 and not R["non_ok_cells"]
    assert all(c["flag"] == "ok" for c in R["cells"]) and R["resources"]["child_processes_started"] == 0
    assert R["resources"]["cpu_seconds_total_process_and_children"] < R["resources"]["cpu_cap_seconds"] and R["resources"]["max_rss_mb"] < R["resources"]["rss_cap_mb"]
    assert json.loads((D / "freeze_eval.json").read_text())["files"]["evaluate.py"] == R["freeze"]["files"]["evaluate.py"]


def test_stored_classifications_follow_from_the_stored_numbers():
    for cid, c in R["evaluation"]["comparisons"].items():
        for m, f in c["families"].items():
            assert f["classification"] == ev.classify(f["delta"], f["interval_95_descriptive"][1], f["blocks_improving"], f["drugs_improving"]), (cid, m)
            assert f["relative_mse_reduction"] == (f["mse_b"] - f["mse_a"]) / f["mse_b"] or abs(f["relative_mse_reduction"] - (f["mse_b"] - f["mse_a"]) / f["mse_b"]) < 1e-12
            assert sum(x < 0 for x in f["block_deltas"]) == f["blocks_improving"]
        assert c["outcome"] == ev.family_outcome(c["families"]["ridge"]["classification"], c["families"]["gbm"]["classification"])
    res = {c: {m: {"rel": R["evaluation"]["comparisons"][c]["families"][m]["relative_mse_reduction"], "cls": R["evaluation"]["comparisons"][c]["families"][m]["classification"]} for m in ("ridge", "gbm")} for c in ("C3", "C4")}
    assert R["evaluation"]["engineering_outcome"] == ev.engineering_outcome(res)


def test_predictions_file_has_every_row_set_model_once_and_no_outcome_values():
    with gzip.open(D / "results/oof_predictions.csv.gz", "rt") as f:
        header = f.readline().strip().split(",")
        n = sum(1 for _ in f)
    assert header == ["drug", "row_id", "outer_fold", "feature_set", "model", "prediction"]
    assert n == sum(c["used"] for c in R["row_counts"].values()) * len(ev.DESIGN["feature_sets"]) * len(ev.DESIGN["models"])
