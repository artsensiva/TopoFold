"""Real-browser check of the focus after the Results folder is changed and changed back (headless Chrome, DevTools protocol).

    .venv/bin/python apps/structpair_viewer/browser_check_folder.py URL OUT_DIR FOLDER_A FOLDER_B [--expect-sync]

URL is a viewer started with `-- --results FOLDER_A`; FOLDER_A holds kras1 and kras2, FOLDER_B is another valid location
(or one whose result fails its manifest check). Replays the sequence reported by AI tester 2: select comparison kras2, focus a
non-default position, change the Results folder to FOLDER_B, change it back to FOLDER_A. After each step, once the page has
settled, records the sidebar focus, the details heading, the focus marker of the plot and the table rows marked selected.
Exit code 0 when all four agree with the focus remembered for kras2 after the return, 1 otherwise.
"""

import json
import sys
import time
from pathlib import Path

URL, OUT, FOLDER_A, FOLDER_B = sys.argv[1], Path(sys.argv[2]), sys.argv[3], sys.argv[4]
sys.argv = sys.argv[:3]                          # browser_check reads URL and OUT at import
sys.path.insert(0, str(Path(__file__).parent))
import browser_check as bc                       # noqa: E402

FOCUS = "39"


def set_folder(b, text):
    box = b.js("""(() => { const i = document.querySelector('[data-testid=stTextInput] input');
        i.scrollIntoView({block: 'center'}); i.focus(); i.select(); const r = i.getBoundingClientRect();
        return [r.x + r.width / 2, r.y + r.height / 2]; })()""")
    b.click_xy(*box)
    b.js("document.querySelector('[data-testid=stTextInput] input').select()")
    b.call("Input.insertText", text=text)
    for kind in ("keyDown", "keyUp"):
        b.call("Input.dispatchKeyEvent", type=kind, key="Enter", code="Enter", windowsVirtualKeyCode=13, text="\r")


def choose_typed(b, label, text, option_prefix):
    """Select an option of a long selectbox by typing into it (the list is virtualised)."""
    box = b.js(f"""(() => {{ const w = [...document.querySelectorAll('[data-testid=stSelectbox]')]
        .find(e => e.innerText.startsWith({json.dumps(label)})); w.scrollIntoView({{block: 'center'}});
        const r = w.querySelector('input[role=combobox]').getBoundingClientRect();
        return [r.x + r.width / 2, r.y + r.height / 2]; }})()""")
    b.click_xy(*box)
    time.sleep(0.5)
    b.call("Input.insertText", text=text)
    time.sleep(0.8)
    box = b.js(f"""(() => {{ const o = [...document.querySelectorAll('[role=option]')]
        .find(e => e.innerText.trim().startsWith({json.dumps(option_prefix)})); if (!o) return null;
        const r = o.getBoundingClientRect(); return [r.x + r.width / 2, r.y + r.height / 2]; }})()""")
    if box is None:
        raise RuntimeError(f"option {option_prefix!r} not found for {label!r}")
    b.click_xy(*box)


def plot_marker(b):
    """Focus markers as drawn: the orange rule of each metric chart (x of the #eb6834 mark within the chart, pixels) and the
    residues of the orange 3D markers (text held by the Plotly figure)."""
    return b.js("""(() => {
        const rules = [];
        for (const c of document.querySelectorAll('[data-testid=stVegaLiteChart]')) {
            const marks = [...c.querySelectorAll('svg path, svg line')].filter(e => (e.getAttribute('stroke') || '').toLowerCase() === '#eb6834');
            const cb = c.getBoundingClientRect();
            rules.push(marks.map(e => Math.round(e.getBoundingClientRect().x - cb.x)));
        }
        const g = document.querySelector('[data-testid=stPlotlyChart] .js-plotly-plot');
        const orange = g && g.data ? g.data.filter(t => t.name && t.name.startsWith('focus window')).map(t => [t.name, t.text]) : null;
        return {chart_rule_x: rules, plotly_focus_traces: orange};
    })()""")


def full_state(b):
    s = b.state()
    scan = b.selected_rows_full_scan()
    s["selected_positions"] = scan["selected"] if scan else None
    s["plot_marker"] = plot_marker(b)
    s["comparison"] = b.js("""(() => { const w = [...document.querySelectorAll('[data-testid=stSelectbox]')]
        .find(e => e.innerText.startsWith('Comparison')); return w ? w.querySelector('input[role=combobox]').value : null; })()""")
    s["folder"] = b.js("document.querySelector('[data-testid=stTextInput] input').value")
    s["error_banner"] = b.js("(document.querySelector('[data-testid=stAlert]') || {}).innerText || null")
    return s


def main():
    b = bc.Browser()
    log = []
    try:
        b.call("Page.navigate", url=URL)
        b.wait_idle(8)
        steps = [("select comparison kras2", lambda: b.choose("Comparison", "kras2")),
                 (f"select focus {FOCUS}", lambda: choose_typed(b, "Focus window", f"{FOCUS} ·", f"{FOCUS} ·")),
                 ("results folder -> B", lambda: set_folder(b, FOLDER_B)),
                 ("results folder -> A (back)", lambda: set_folder(b, FOLDER_A))]
        for i, (name, act) in enumerate(steps):
            act()
            b.wait_idle(4)
            s = full_state(b)
            log.append({"step": i + 1, "action": name, **s})
            b.shot(f"{i + 1:02d}-{name.replace(' ', '_').replace('>', 'to').replace('(', '').replace(')', '')}")
            time.sleep(0.5)
    finally:
        b.close()
    final = log[-1]
    pos = final["focus"].split(" ")[0] if final["focus"] else None
    synced = log[1]                                                      # focus just selected, page settled
    agree = {"sidebar_is_remembered_focus": pos == FOCUS,
             "details_agree": final["details"] == f"Focus window at position {FOCUS} (3 residues)",
             "table_row_agrees": final["selected_positions"] == [FOCUS],
             "chart_rule_unchanged_from_step_2": final["plot_marker"]["chart_rule_x"]
             == synced["plot_marker"]["chart_rule_x"] and bool(synced["plot_marker"]["chart_rule_x"]),
             "3d_markers_unchanged_from_step_2": final["plot_marker"]["plotly_focus_traces"]
             == synced["plot_marker"]["plotly_focus_traces"]}
    ok = all(agree.values())
    (OUT / "folder_switch_check.json").write_text(json.dumps(
        {"url": URL, "folders": [FOLDER_A, FOLDER_B], "focus_set": FOCUS, "steps": log, "agree_after_return": agree,
         "all_agree": ok}, indent=1) + "\n")
    print(json.dumps({"agree_after_return": agree, "last": final}, indent=1))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
