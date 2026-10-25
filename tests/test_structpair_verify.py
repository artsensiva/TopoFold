"""structpair_verify: the documented manifest-check procedure (D77). Uses a SYNTHETIC result built in a temporary directory,
because the shipped manifests record absolute input paths of the machine that made them."""

import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from topofold_inputs import structpair_verify as V  # noqa: E402

ENV = {"PYTHONPATH": str(ROOT / "python"), "PATH": ""}


def _h(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def _synthetic(tmp_path):
    d = tmp_path / "SYNTHETIC_result"
    d.mkdir()
    (tmp_path / "SYNTHETIC_A.pdb").write_text("END\n")
    (tmp_path / "SYNTHETIC_B.pdb").write_text("TER\nEND\n")
    (tmp_path / "SYNTHETIC_spec.json").write_text("{}\n")
    for name in ("windows.csv", "residues.csv", "summary.json"):
        (d / name).write_text(f"SYNTHETIC {name}\n1.5,2.5\n")
    manifest = {"schema": "topofold_inputs.structpair/3",
                "outputs_sha256": {n: _h(d / n) for n in ("windows.csv", "residues.csv", "summary.json")},
                "inputs": {k: {"path": str(tmp_path / f"SYNTHETIC_{k}.pdb"), "sha256": _h(tmp_path / f"SYNTHETIC_{k}.pdb")} for k in "AB"},
                "spec": {"path": str(tmp_path / "SYNTHETIC_spec.json"), "sha256": _h(tmp_path / "SYNTHETIC_spec.json")}}
    (d / "manifest.json").write_text(json.dumps(manifest))
    return d


def test_every_listed_file_is_checked_and_matches(tmp_path):
    r = V.verify(_synthetic(tmp_path))
    assert r["n_checked"] == 6 and r["n_ok"] == 6 and r["all_ok"] is True        # 3 outputs + 2 inputs + the spec


def test_one_altered_value_is_a_mismatch_naming_the_file(tmp_path):
    d = _synthetic(tmp_path)
    (d / "windows.csv").write_text((d / "windows.csv").read_text().replace("1.5", "2.5"))
    r = V.verify(d)
    assert [c["path"] for c in r["checks"] if c["status"] == "mismatch"] == ["windows.csv"]
    assert r["n_mismatch"] == 1 and r["n_ok"] == 5 and not r["all_ok"]


def test_a_missing_file_is_not_a_pass_and_the_command_exit_codes(tmp_path):
    d = _synthetic(tmp_path)
    run = lambda *a: subprocess.run([sys.executable, "-m", "topofold_inputs.structpair_verify", *map(str, a)],   # noqa: E731
                                    capture_output=True, text=True, cwd=ROOT, env=ENV)
    ok = run(d, "--json", tmp_path / "o.json")
    assert ok.returncode == 0 and "6 of 6 files match the manifest (0 mismatch, 0 missing)" in ok.stdout
    assert json.loads((tmp_path / "o.json").read_text())["all_ok"] is True
    (d / "residues.csv").unlink()
    (tmp_path / "SYNTHETIC_A.pdb").unlink()                                      # a moved input is a failure of the check
    bad = run(d)
    assert bad.returncode == 1 and "missing  output   residues.csv" in bad.stdout and "missing  input A" in bad.stdout
    assert "4 of 6 files match the manifest (0 mismatch, 2 missing)" in bad.stdout
    assert run(tmp_path / "none").returncode == 2
