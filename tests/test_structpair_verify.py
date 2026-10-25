"""structpair_verify: the documented manifest-check procedure (D77)."""

import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs import structpair_verify as V  # noqa: E402

KRAS2 = ROOT / "benchmarks/audit/structpair_demo/results/kras2"


def _copy(tmp_path):
    d = tmp_path / "r"
    shutil.copytree(KRAS2, d)
    return d


def test_shipped_result_matches_its_manifest():
    r = V.verify(KRAS2)
    n_outputs = len(json.loads((KRAS2 / "manifest.json").read_text())["outputs_sha256"])
    assert r["n_checked"] == n_outputs + 3 and r["all_ok"] is True        # outputs + two inputs + spec


def test_one_altered_value_is_a_mismatch_naming_the_file(tmp_path):
    d = _copy(tmp_path)
    text = (d / "windows.csv").read_text()
    (d / "windows.csv").write_text(text.replace("1.", "2.", 1))
    r = V.verify(d)
    bad = [c["path"] for c in r["checks"] if c["status"] == "mismatch"]
    assert bad == ["windows.csv"] and r["n_mismatch"] == 1 and not r["all_ok"]


def test_a_missing_file_is_not_a_pass_and_the_command_exit_codes(tmp_path):
    d = _copy(tmp_path)
    ok = subprocess.run([sys.executable, "-m", "topofold_inputs.structpair_verify", str(d), "--json", str(tmp_path / "o.json")],
                        capture_output=True, text=True, cwd=ROOT, env={"PYTHONPATH": str(ROOT / "python"), "PATH": ""})
    assert ok.returncode == 0 and "files match the manifest (0 mismatch, 0 missing)" in ok.stdout
    assert json.loads((tmp_path / "o.json").read_text())["all_ok"] is True
    (d / "residues.csv").unlink()
    bad = subprocess.run([sys.executable, "-m", "topofold_inputs.structpair_verify", str(d)], capture_output=True, text=True,
                         cwd=ROOT, env={"PYTHONPATH": str(ROOT / "python"), "PATH": ""})
    assert bad.returncode == 1 and "missing  output   residues.csv" in bad.stdout
    assert subprocess.run([sys.executable, "-m", "topofold_inputs.structpair_verify", str(tmp_path / "none")],
                          capture_output=True, text=True, cwd=ROOT, env={"PYTHONPATH": str(ROOT / "python"), "PATH": ""}).returncode == 2
