"""Real-browser check of table selection vs focus in the structpair viewer (headless Chrome via DevTools protocol).

    .venv/bin/python apps/structpair_viewer/browser_check.py URL OUT_DIR

Replays the sequence reported by AI tester 1 (load kras2; metrics d_kappa + rmsd3; sort by d_kappa descending; rows
"without values"; select the second row; rows "all") and records, after each step, the focus shown in the sidebar and
the rows the table marks as selected (from the grid's accessibility tree). Needs google-chrome and websocket-client.
Exit code 0 when the table selection and the focus agree at every step, 1 otherwise.
"""

import base64
import json
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

import websocket

URL, OUT = sys.argv[1], Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
PORT = 9334
SIZE = 3


class Browser:
    def __init__(self):
        self.proc = subprocess.Popen(
            ["google-chrome", "--headless=new", "--disable-gpu", "--no-sandbox", "--hide-scrollbars",
             f"--remote-debugging-port={PORT}", f"--remote-allow-origins=http://127.0.0.1:{PORT}",
             f"--user-data-dir={tempfile.mkdtemp()}", "about:blank"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for _ in range(100):
            try:
                tabs = json.load(urllib.request.urlopen(f"http://127.0.0.1:{PORT}/json"))
                break
            except OSError:
                time.sleep(0.2)
        page = [t for t in tabs if t["type"] == "page"][0]
        self.ws = websocket.create_connection(page["webSocketDebuggerUrl"], timeout=60, suppress_origin=True)
        self.n = 0
        self.call("Emulation.setDeviceMetricsOverride", width=1600, height=2400, deviceScaleFactor=1, mobile=False)

    def call(self, method, **params):
        self.n += 1
        self.ws.send(json.dumps({"id": self.n, "method": method, "params": params}))
        while True:
            m = json.loads(self.ws.recv())
            if m.get("id") == self.n:
                if "error" in m:
                    raise RuntimeError(m["error"])
                return m.get("result", {})

    def js(self, expr):
        r = self.call("Runtime.evaluate", expression=expr, returnByValue=True, awaitPromise=True)
        return r.get("result", {}).get("value")

    def wait_idle(self, t=2.5):
        time.sleep(t)
        for _ in range(60):
            if not self.js("!!document.querySelector('[data-testid=stStatusWidget] [data-testid=stStatusWidgetRunningIcon]')"):
                break
            time.sleep(0.25)
        time.sleep(0.8)

    def click_xy(self, x, y):
        for kind in ("mousePressed", "mouseReleased"):
            self.call("Input.dispatchMouseEvent", type=kind, x=x, y=y, button="left", clickCount=1)

    def click_text(self, selector, text):
        box = self.js(f"""(() => {{ const e = [...document.querySelectorAll({json.dumps(selector)})]
            .find(b => b.innerText.trim() === {json.dumps(text)}); if (!e) return null;
            e.scrollIntoView({{block: 'center'}}); const r = e.getBoundingClientRect();
            return [r.x + r.width / 2, r.y + r.height / 2]; }})()""")
        if box is None:
            raise RuntimeError(f"no {selector} with text {text!r}")
        self.click_xy(*box)

    def choose(self, label, option_prefix):
        box = self.js(f"""(() => {{ const w = [...document.querySelectorAll('[data-testid=stSelectbox]')]
            .find(e => e.innerText.startsWith({json.dumps(label)})); w.scrollIntoView({{block: 'center'}});
            const r = w.querySelector('button[aria-label=Open]').getBoundingClientRect();
            return [r.x + r.width / 2, r.y + r.height / 2]; }})()""")
        self.click_xy(*box)
        time.sleep(0.8)
        box = self.js(f"""(() => {{ const o = [...document.querySelectorAll('[role=option]')]
            .find(e => e.innerText.trim().startsWith({json.dumps(option_prefix)})); if (!o) return null;
            o.scrollIntoView({{block: 'center'}}); const r = o.getBoundingClientRect();
            return [r.x + r.width / 2, r.y + r.height / 2]; }})()""")
        if box is None:
            raise RuntimeError(f"option {option_prefix!r} not found for {label!r}")
        self.click_xy(*box)

    def state(self):
        return self.js("""(() => {
            const focusBox = [...document.querySelectorAll('[data-testid=stSelectbox]')]
                .find(e => e.innerText.startsWith('Focus window'));
            const focus = focusBox ? focusBox.querySelector('input[role=combobox]').value : null;
            const details = ([...document.querySelectorAll('p')].map(e => e.innerText)
                .find(t => t.startsWith('Focus window at position')) || null);
            const tables = [...document.querySelectorAll('[data-testid=stDataFrame] table')];
            let selected = [];
            for (const t of tables) {
                const head = [...t.querySelectorAll('th')].map(h => h.innerText.trim());
                if (!head.includes('window3_A')) continue;
                const pi = head.indexOf('position');
                for (const tr of t.querySelectorAll('tr[aria-selected=true]'))
                    selected.push(tr.querySelectorAll('td, th')[pi]?.innerText.trim());
            }
            return {focus, details, selected_positions: selected};
        })()""")

    def selected_rows_full_scan(self):
        """Scroll the (virtualised) window table from top to bottom and collect every row marked selected."""
        return self.js("""(async () => {
            const d = [...document.querySelectorAll('[data-testid=stDataFrame]')]
                .filter(d => [...d.querySelectorAll('th')].some(h => h.innerText.trim() === 'window3_A'))[0];
            if (!d) return null;
            const sc = d.querySelector('.dvn-scroller'); const found = new Set(); let seen = 0;
            const pc = [...d.querySelectorAll('th')].find(h => h.innerText.trim() === 'position').getAttribute('aria-colindex');
            for (let top = 0; top <= sc.scrollHeight + 200; top += 150) {
                sc.scrollTop = top; await new Promise(r => setTimeout(r, 350));
                for (const tr of d.querySelectorAll('tr[aria-rowindex]')) {
                    const cell = tr.querySelector(`td[aria-colindex="${pc}"]`); if (!cell) continue;
                    seen = Math.max(seen, +tr.getAttribute('aria-rowindex') - 1);
                    if (tr.getAttribute('aria-selected') === 'true') found.add(cell.innerText.trim());
                }
            }
            sc.scrollTop = 0; await new Promise(r => setTimeout(r, 150));
            return {selected: [...found], rows_seen: seen};
        })()""")

    def grid_row_marker(self, row):
        """Screen point of the selection checkbox of data row `row` (0-based) in the window table."""
        return self.js(f"""(() => {{ const d = [...document.querySelectorAll('[data-testid=stDataFrame]')]
            .filter(d => d.querySelector('table') && [...d.querySelectorAll('th')].some(h => h.innerText.trim() === 'window3_A'))[0];
            d.scrollIntoView({{block: 'start'}}); const r = d.getBoundingClientRect();
            return [r.x + 14, r.y + 35 * ({row} + 1) + 17]; }})()""")

    def shot(self, name):
        data = self.call("Page.captureScreenshot", format="png")["data"]
        (OUT / f"{name}.png").write_bytes(base64.b64decode(data))

    def close(self):
        self.proc.terminate()


def main():
    b = Browser()
    log, ok = [], True
    try:
        b.call("Page.navigate", url=URL)
        b.wait_idle(8)
        steps = [("load kras2", lambda: b.choose("Comparison", "kras2")),
                 ("sort by d_kappa descending", lambda: b.choose("Sort by", "d_kappa")),
                 ("rows: without values", lambda: b.click_text("button", "without values")),
                 ("select second row", lambda: b.click_xy(*b.grid_row_marker(1))),
                 ("rows: all", lambda: b.click_text("button", "all")),
                 ("sort by position", lambda: b.choose("Sort by", "position")),
                 ("switch to kras1", lambda: b.choose("Comparison", "kras1"))]
        # whether the focus row is among the displayed rows after each step (known from the data: positions 79 and
        # 84 have values, so "without values" hides them; 58 has no three-residue value in kras2)
        expected_displayed = {"rows: without values": False}
        for i, (name, act) in enumerate(steps):
            act()
            b.wait_idle()
            s = b.state()
            scan = b.selected_rows_full_scan()
            s["selected_positions"], s["table_rows_scanned"] = scan["selected"], scan["rows_seen"]
            s["focus_in_displayed_rows"] = expected_displayed.get(name, True)
            pos = s["focus"].split(" ")[0] if s["focus"] else None
            want = [pos] if s["focus_in_displayed_rows"] else []
            agree = s["selected_positions"] == want and s["details"] == f"Focus window at position {pos} ({SIZE} residues)"
            ok &= agree
            log.append({"step": i + 1, "action": name, **s, "selection_agrees_with_focus": agree})
            b.shot(f"{i + 1:02d}-{name.replace(' ', '_').replace(':', '')}")
    finally:
        b.close()
    (OUT / "browser_check.json").write_text(json.dumps({"url": URL, "steps": log, "all_agree": ok}, indent=1) + "\n")
    print(json.dumps(log, indent=1))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
