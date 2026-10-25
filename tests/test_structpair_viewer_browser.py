"""Real-browser regression checks of the structpair viewer (headless Chrome, DevTools protocol).

Skipped when google-chrome, websocket-client or streamlit is missing. Each case starts the viewer on a free local port.

- Folder switch (AI tester 2, D77): select comparison kras2, a non-default focus, change the Results folder to a location whose
  result fails its manifest check, change it back. After the page settles, the sidebar, details, table row, chart marker and 3D
  markers must show the same focus. On commit 19b5268 the sidebar showed position 0 and the rest position 39.
- Ordinary switching (AI tester 1, D69): the table selection follows the focus through sort, filter and comparison switches.
"""

import json
import shutil
import socket
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
APP = ROOT / "apps/structpair_viewer"
RESULTS = ROOT / "benchmarks/audit/structpair_demo/results"
ABL = ROOT / "benchmarks/audit/structpair_abl_demo_2026-10-25/results/run"

pytestmark = pytest.mark.skipif(
    shutil.which("google-chrome") is None or shutil.which("streamlit") is None and not (Path(sys.executable).parent / "streamlit").exists(),
    reason="needs google-chrome and streamlit")
try:
    import websocket  # noqa: F401
except ImportError:                                           # pragma: no cover
    pytestmark = pytest.mark.skip(reason="needs websocket-client")


def _free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture
def viewer(tmp_path):
    procs = []

    def start(results: Path) -> str:
        port = _free_port()
        exe = Path(sys.executable).parent / "streamlit"
        p = subprocess.Popen([str(exe) if exe.exists() else "streamlit", "run", str(APP / "app.py"), "--server.address", "localhost",
                              "--server.port", str(port), "--server.headless", "true", "--browser.gatherUsageStats", "false",
                              "--", "--results", str(results)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=ROOT)
        procs.append(p)
        url = f"http://localhost:{port}"
        for _ in range(100):
            try:
                urllib.request.urlopen(url + "/_stcore/health", timeout=2)
                return url
            except OSError:
                time.sleep(0.2)
        raise RuntimeError("viewer did not start")

    yield start
    for p in procs:
        p.terminate()


def _run(script, *args, timeout=240):
    return subprocess.run([sys.executable, str(APP / script), *map(str, args)], capture_output=True, text=True, timeout=timeout)


def test_focus_survives_a_results_folder_change_to_a_failing_result_and_back(tmp_path, viewer):
    outputs, tamper = tmp_path / "outputs", tmp_path / "tamper"
    shutil.copytree(RESULTS / "kras2", outputs / "kras2")
    shutil.copytree(ABL, outputs / "abl_1OPJ_AB_vs_1OPK_A")
    shutil.copytree(ABL, tamper / "abl_tampered")
    w = tamper / "abl_tampered/windows.csv"
    lines = w.read_text().splitlines(keepends=True)
    lines[2] = lines[2].replace("0.", "0.9", 1) if "0." in lines[2] else lines[2] + "x"
    w.write_text("".join(lines))                                         # one altered value; the manifest is not updated
    out = tmp_path / "check"
    r = _run("browser_check_folder.py", viewer(outputs), out, outputs, tamper)
    rec = json.loads((out / "folder_switch_check.json").read_text())
    assert rec["steps"][2]["error_banner"].startswith("This result cannot be displayed"), rec["steps"][2]   # the folder really failed
    assert r.returncode == 0 and rec["all_agree"], rec["agree_after_return"]
    assert rec["steps"][3]["focus"].split(" ")[0] == "39"


def test_focus_survives_a_results_folder_change_to_another_valid_result_and_back(tmp_path, viewer):
    first, second = tmp_path / "first", tmp_path / "second"
    shutil.copytree(RESULTS / "kras2", first / "kras2")
    shutil.copytree(ABL, second / "abl_1OPJ_AB_vs_1OPK_A")
    out = tmp_path / "check"
    r = _run("browser_check_folder.py", viewer(first), out, first, second)
    rec = json.loads((out / "folder_switch_check.json").read_text())
    assert r.returncode == 0 and rec["all_agree"], rec["agree_after_return"]
    assert rec["steps"][2]["comparison"].startswith("abl_1OPJ")                # the other location really was displayed


def test_table_selection_follows_focus_through_ordinary_switching(tmp_path, viewer):
    shown = tmp_path / "shown"
    shutil.copytree(RESULTS / "kras1", shown / "kras1")
    shutil.copytree(RESULTS / "kras2", shown / "kras2")
    out = tmp_path / "check"
    r = _run("browser_check.py", viewer(shown), out)
    rec = json.loads((out / "browser_check.json").read_text())
    assert r.returncode == 0 and rec["all_agree"] and len(rec["steps"]) == 7
