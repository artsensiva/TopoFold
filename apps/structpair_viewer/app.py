"""TopoFold Structural Comparison: local, read-only viewer for `topofold_inputs.structpair` results.

    streamlit run apps/structpair_viewer/app.py --server.address localhost [-- --results DIR]

DIR is one structpair output directory or a folder of them (default: the KRAS demonstrations in
benchmarks/audit/structpair_demo/results). Nothing is recomputed: values come from the result files, which are only
read. Needs streamlit, pandas, altair, plotly and numpy; not the compiled TopoFold core or MDAnalysis.
"""

import argparse
import functools
import json
import sys
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
from results import (METRICS, ResultError, discover, display_coordinates, exclusion_explanations,  # noqa: E402
                     table_key_and_selection,
                     load_result, locate_input, metric_frame, read_ca, window_members)

REPO = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS = REPO / "benchmarks/audit/structpair_demo/results"
BLUE, GREY, HIGHLIGHT, BAND = "#2a78d6", "#9d9c96", "#eb6834", "#e4e3df"


def cli_results_dir() -> Path:
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--results", default=str(DEFAULT_RESULTS))
    a, _ = ap.parse_known_args(sys.argv[1:])
    return Path(a.results)


@st.cache_data(max_entries=32)
def cached_ca(path: str, sha: str, chain: str, altloc: str):
    """Parsed Cα coordinates, keyed by the file's verified SHA-256 (a changed file is a new key, never a stale hit)."""
    return read_ca(path, chain, altloc, expected_sha256=sha)


def structure_coords(res):
    """Locate and hash both input structures on every rerun (not cached); a missing or changed file disables the view."""
    coords, status = {}, {}
    for k in "AB":
        f, msg = locate_input(res, k, REPO)
        status[k] = msg
        if f is not None:
            inp = res.manifest["inputs"][k]
            try:
                coords[k] = cached_ca(str(f), inp["sha256"], inp["chain"], inp["altloc_policy"])
            except ResultError as e:
                status[k] = str(e)
    return coords, status


def load_name(path: Path) -> str:
    try:
        return json.loads((path / "manifest.json").read_text()).get("name") or path.name
    except (OSError, ValueError):
        return f"{path.name} (unreadable manifest)"


def focus_from_table(key: str, positions: list[int]):
    """Selecting a table row moves the focus window (runs before the rerun, so the focus widget can be updated).

    The row index is translated to a position through the rows displayed when the click happened. Clearing the
    selection keeps the focus; the table is then rebuilt marking the focus row again."""
    rows = st.session_state[key].selection.rows
    if rows:
        st.session_state.focus = positions[rows[0]]


st.set_page_config(page_title="TopoFold Structural Comparison", page_icon=":material/compare:", layout="wide")
st.title("TopoFold Structural Comparison")
st.caption("Local, read-only viewer for two-structure comparisons made with `python -m topofold_inputs.structpair`. "
           "Static and descriptive: direct geometric differences between two coordinate sets.")

# --- choose a result -------------------------------------------------------------------------------------------------
# True only if the previous run drew the focus widget. A run that stopped earlier (no result found, or one that failed its
# manifest check) drops the widget from the browser while its value survives on the server, so the next run must push it again.
focus_widget_was_drawn = st.session_state.pop("_focus_drawn", False)
if "results_root" not in st.session_state:
    st.session_state.results_root = str(cli_results_dir())
with st.sidebar:
    st.header("Results")
    root_text = st.text_input("Results folder", key="results_root",
                              help="One structpair output directory, or a folder containing several. After a change the first comparison "
                                   "of the new folder is shown; the focus of each comparison is remembered.")
    found = discover(root_text)
    if not found:
        st.error(f"No structpair results (manifest.json) found in `{root_text}`.")
        st.stop()
    order = sorted(found, key=lambda k: (not k.startswith(("kras1", "kras2")), k))   # the two KRAS studies first
    names = {k: f"{k} — {load_name(found[k])}" for k in order}
    choice = st.selectbox("Comparison", order, key="result_choice", format_func=names.get)

path = found[choice]
try:
    res = load_result(path)          # re-read and hash-verified on every rerun; deliberately not cached
except ResultError as e:
    st.error(f"This result cannot be displayed: {e}")
    st.stop()

with st.sidebar:
    metric_keys = st.pills("Metrics", list(METRICS), selection_mode="multi", default=["d_kappa", "rmsd3"],
                           format_func=lambda k: f"{k} ({METRICS[k][1]})", key="metrics")
    size = st.segmented_control("Highlight window size", [3, 4], default=3, key="hl_size",
                                format_func=lambda s: f"{s} residues")
    positions = res.windows["position"].tolist()
    lab = {p: f"{p} · A {a} ↔ B {b}" for p, a, b in zip(res.windows["position"], res.windows["A_residue"],
                                                        res.windows["B_residue"])}
    # focus is remembered per comparison: a position index means a different residue pair in another comparison
    remembered = st.session_state.setdefault("_focus_by_result", {})
    previous = st.session_state.get("_focus_result")
    if previous != str(path):
        if previous is not None and "focus" in st.session_state:
            remembered[previous] = st.session_state.focus
        st.session_state.focus = remembered.get(str(path), positions[len(positions) // 2])
        st.session_state["_focus_result"] = str(path)
    elif not focus_widget_was_drawn:
        st.session_state.focus = remembered.get(str(path), positions[len(positions) // 2])
    if st.session_state.get("focus") not in positions:
        st.session_state.focus = positions[len(positions) // 2]
    focus = st.selectbox("Focus window (position · residue in A ↔ in B)", positions, key="focus",
                         format_func=lambda p: lab[p])
    remembered[str(path)] = focus
    st.session_state["_focus_drawn"] = True

st.subheader(res.name)
S = res.summary
with st.container(horizontal=True):
    for k in "AB":
        inp = res.manifest["inputs"][k]
        chs = inp.get("chains") or [inp["chain"]]
        chain_text = f"chain {inp['chain']}" if len(chs) == 1 else f"chains {', '.join(chs)} (primary {inp['chain']})"
        st.metric(f"Structure {k}", inp["label"], f"{chain_text} · altloc {inp['altloc_policy']}",
                  delta_color="off", delta_arrow="off", border=True)
    g = S["global"]
    st.metric("Global Cα RMSD (context only)", f"{g['rmsd_qcp_angstrom']:.3f} Å", f"{g['n_ca']} Cα, MDAnalysis QCP",
              delta_color="off", delta_arrow="off", border=True)
    st.metric("Excluded window rows", str(len(res.exclusions)), "see Exclusions tab", delta_color="off",
              delta_arrow="off", border=True)

members = window_members(res, focus, size or 3)
member_rows = res.residues.iloc[members]

# --- sequence plots and molecular view -------------------------------------------------------------------------------
left, right = st.columns([1.15, 1.0])
with left:
    st.markdown("**Along the chain** — one plot per metric (units differ, so axes are never shared). Grey bands: "
                "windows inside the mapped range without a value; see Exclusions.")
    if not metric_keys:
        st.info("Select at least one metric in the sidebar.")
    for mk in metric_keys:
        df = metric_frame(res, mk)
        unit, wsize, label = METRICS[mk][1], METRICS[mk][2], METRICS[mk][3]
        df = df.assign(x0=df.position - 0.5, x1=df.position + 0.5, focus=df.position == focus,
                       has_value=df["value"].notna(), seg=df["segment"].fillna(-1).astype(int))
        df = df.astype({c: object for c in ("A_residue", "B_residue", "window_A", "status")})
        ymax = float(df["value"].max()) if df["has_value"].any() else 1.0
        xs = alt.Scale(domain=[-0.5, len(res.windows) - 0.5], nice=False)
        x = alt.X("position:Q", title="window position (hover for residues in A and B)", scale=xs)
        y = alt.Y("value:Q", title=f"{unit}", scale=alt.Scale(domain=[0, ymax * 1.05], nice=False))
        tip = [alt.Tooltip("window_A:N", title="window (A numbering)"), alt.Tooltip("A_residue:N", title="centre in A"),
               alt.Tooltip("B_residue:N", title="centre in B"), alt.Tooltip("value:Q", title=f"value ({unit})", format=".4f"),
               alt.Tooltip("status:N"), alt.Tooltip("altloc_in_window:N", title="altloc in window"),
               alt.Tooltip("mismatch_in_window:N", title="residue-name mismatch")]
        base = alt.Chart(df)
        layers = [base.transform_filter("datum.status == 'excluded' || datum.status == 'unavailable'")
                      .mark_rect(color=BAND).encode(x=alt.X("x0:Q", scale=xs, title="window position (hover for residues in A and B)"), x2="x1:Q", tooltip=tip),
                  base.transform_filter("datum.has_value").mark_line(color=BLUE, strokeWidth=1.5)
                      .encode(x=x, y=y, detail="seg:N"),
                  base.transform_filter("datum.has_value").mark_circle(size=50, opacity=0).encode(x=x, y=y, tooltip=tip),
                  base.transform_filter("datum.focus").mark_rule(color=HIGHLIGHT, strokeWidth=2).encode(x=x)]
        st.altair_chart(alt.layer(*layers).properties(height=180, title=f"{mk}: {label}, {wsize}-residue windows"),
                        key=f"chart_{mk}")

with right:
    st.markdown("**Molecular view** — Cα traces; A superposed on B for display only (all residues present in both). "
                f"Orange: residues of the focus window ({size or 3} residues). The 3D view needs WebGL in the browser; if it stays "
                "empty, the plots, the focus table and the details below show the same window.")
    coords, cstatus = structure_coords(res)
    if set(coords) != {"A", "B"}:
        st.warning("Molecular view unavailable. " + " · ".join(f"{k}: {m}" for k, m in cstatus.items()))
    else:
        A, B = display_coordinates(res, coords)
        brk = {b["after_index"] for k in "AB" for b in res.manifest["breaks"][k]}

        def with_breaks(X):
            rows = []
            for i in range(len(X)):
                rows.append(X[i])
                if i in brk:
                    rows.append([np.nan] * 3)
            return np.array(rows)

        def hover(col):
            out = []
            for i, r in res.residues.iterrows():
                out.append(f"pos {r['position']}: {r[col]} {r[col.replace('residue', 'resname')]}")
                if r["position"] in brk:
                    out.append("")
            return out

        fig = go.Figure()
        for key, X, color, width in (("A", A, GREY, 3), ("B", B, BLUE, 5)):
            fig.add_trace(go.Scatter3d(x=with_breaks(X)[:, 0], y=with_breaks(X)[:, 1], z=with_breaks(X)[:, 2],
                                       mode="lines", line=dict(color=color, width=width),
                                       text=hover(f"{key}_residue"), hoverinfo="text",
                                       name=f"{key}: {res.label(key)}"))
        for key, X, sym in (("A", A, "circle-open"), ("B", B, "circle")):
            sel = [i for i in members if np.isfinite(X[i]).all()]
            fig.add_trace(go.Scatter3d(x=X[sel, 0], y=X[sel, 1], z=X[sel, 2], mode="markers",
                                       marker=dict(size=6, color=HIGHLIGHT, symbol=sym),
                                       text=[f"{res.residues.iloc[i][f'{key}_residue']}" for i in sel], hoverinfo="text",
                                       name=f"focus window in {key}"))
        fig.update_layout(height=460, margin=dict(l=0, r=0, t=0, b=0), legend=dict(orientation="h", y=-0.02, yanchor="top"),
                          scene=dict(xaxis_visible=False, yaxis_visible=False, zaxis_visible=False, aspectmode="data"))
        st.plotly_chart(fig, key="mol")
        st.caption(" · ".join(cstatus.values()))
    with st.container(border=True):
        st.markdown(f"**Focus window at position {focus}** ({size or 3} residues)")
        st.dataframe(member_rows[["position", "A_residue", "A_resname", "A_status", "B_residue", "B_resname",
                                  "B_status"]], hide_index=True, key="focus_rows")
        ex = exclusion_explanations(res, focus, size or 3)
        if len(ex):
            for _, e in ex.iterrows():
                st.warning(f"Excluded in structure {e['structure']}: **{e['cause']}** at {e['location'] or 'n/a'}. "
                           f"{e['explanation']} Record: {e['detail']} ({e['source']}).{' ' + e['detail_reading'] if e['detail_reading'] else ''}")
        else:
            vals = {mk: res.windows.loc[res.windows.position == focus, METRICS[mk][0]].iloc[0]
                    for mk in METRICS if METRICS[mk][2] == (size or 3)}
            st.markdown(" · ".join(f"{mk} = {'unavailable' if pd.isna(v) else f'{v:.4f}'} {METRICS[mk][1]}"
                                   for mk, v in vals.items()))

# --- tables ----------------------------------------------------------------------------------------------------------
VAR = S.get("within_entry_variation")
if VAR and VAR.get("status") in ("available", "partial"):
    st.warning(f"Within-entry crystal-copy variation ({VAR['status']}, see its tab; exclusions and provenance are unchanged): "
               + VAR["statement"])
elif VAR:
    st.caption("Within-entry crystal-copy variation: unavailable (one chain per entry).")
MOD = S.get("modified_residues")
if isinstance(MOD, dict):
    kept = "; ".join(f"{r['structure']} {r['residue']} {r['deposited_name']} (parent {r['parent']}, {r['annotation'] or 'no annotation text'})" for r in MOD["retained"]) or "none"
    st.info(f"**Modified residues (policy `{MOD['policy']}`).** Retained: {kept}. Deposited names are not renamed to the parent, and SEP versus SER stays a visible "
            "name difference in the tables; an annotation is not evidence of functional state.")
    if MOD["excluded"]:
        st.warning("**Modified-residue names excluded (not admitted by the policy):** " + "; ".join(f"{e['structure']} {e['residue']} {e['name']}: {e['reason']}" for e in MOD["excluded"]))
ALT = S.get("altloc_sensitivity")
if isinstance(ALT, dict) and ALT.get("available") is False:
    st.warning("**The requested alternate-conformation sensitivity is UNAVAILABLE (not zero).** " + ALT["reason"] + " " + ALT["why_not_partial"]
               + " Tied residues: " + "; ".join(f"{t['structure']} {t['chain']}:{t['resid']}{t['icode']} (altlocs "
                                                 f"{', '.join(a or 'blank' for a in t['altlocs'])}, occupancies "
                                                 f"{', '.join(str(o) for o in t['occupancies'])}"
                                                 f"{'' if t['in_correspondence'] else ', outside the correspondence'})" for t in ALT["ties"]) + ".")
t_win, t_var, t_exc, t_res, t_prov, t_rep = st.tabs(["Windows", "Crystal-copy variation", "Exclusions",
                                                     "Residue correspondence", "Provenance", "Report"])
with t_win:
    with st.container(horizontal=True):
        show = st.segmented_control("Rows", ["all", "with values", "without values", "altloc in window",
                                             "name mismatch in window"], default="all", key="row_filter")
        top_metric = st.selectbox("Sort by", ["position"] + list(METRICS), key="sort_by",
                                  format_func=lambda k: k if k == "position" else f"{k} ({METRICS[k][1]}, descending)")
    w = res.windows.copy()
    mcols = [c for c, *_ in METRICS.values()]
    if show == "with values":
        w = w[w[mcols].notna().any(axis=1)]
    elif show == "without values":
        w = w[w[mcols].isna().all(axis=1)]
    elif show == "altloc in window":
        w = w[w["altloc_in_window3"] | w["altloc_in_window4"]]
    elif show == "name mismatch in window":
        w = w[w["mismatch_in_window3"] | w["mismatch_in_window4"]]
    if top_metric != "position":
        w = w.sort_values(METRICS[top_metric][0], ascending=False, na_position="last")
    st.caption("Empty cells are unavailable values (outside the window range or excluded); they are not zero. "
               "Select a row to focus it in the plots and the molecular view.")
    numcfg = {c: st.column_config.NumberColumn(c, format="%.4f") for c in w.columns if c.endswith(("_rad", "_angstrom"))}
    shown = w["position"].tolist()
    tkey, tsel = table_key_and_selection(str(path), shown, focus)
    if focus not in shown:
        st.caption(f"The focus window (position {focus}) is not among the displayed rows, so no row is marked.")
    st.dataframe(w, hide_index=True, column_config=numcfg, on_select=functools.partial(focus_from_table, tkey, shown),
                 selection_mode="single-row", selection_default=tsel, key=tkey)
with t_var:
    if res.variation is None:
        st.info("Not available for this result: " + ((VAR or {}).get("reason") or "no chain_variation.csv was written")
                + ". Name every chain of an entry in `structures.X.chains` (with `primary_chain`) to obtain it.")
    else:
        st.warning(VAR["statement"] if VAR else "Descriptive only.")
        adv = VAR["parameters"]["advisory"]
        if res.advisory_semantics == "advisory/1":
            st.info("Record of schema /2 (structpair 0.2.0), read with its historical meaning: the flag column `separated` used a "
                    "default minimum difference of 0 and a tolerance-widened comparison, so it can be true where every value is zero. "
                    "Read it only where a minimum difference was set explicitly; results of schema /3 use the advisory rule below.")
            st.markdown(f"Chains A {VAR['chains']['A']}, B {VAR['chains']['B']}; pairs {VAR['n_pairs']}. Recorded flag rule: the between-entry "
                        f"**{adv['statistic']}** at least max(floor, {adv['within_factor']:g} × the within-entry maximum), floors "
                        f"{adv['min_difference'] or 'none (0)'}. Empty cells are unavailable values, never zero; `*_n` counts the chain pairs.")
        else:
            floors = adv["min_difference"]
            rule = ("**No advisory rule is configured, so no flag is assessed**; the numeric values are shown. Set "
                    "`chain_variation.advisory.min_difference` to a finite positive number for a channel to obtain a flag for it."
                    if not floors else
                    "Advisory rule (not significance, not biological separation), configured for "
                    + ", ".join(f"{c} (≥ {v:g})" for c, v in floors.items())
                    + f": the between-entry **{adv['statistic']}** is at least max(minimum difference, {adv['within_factor']:g} × the "
                    "within-entry maximum); equality meets the rule; `advisory_at_boundary` marks values within 1e-9 of the threshold.")
            st.markdown(f"Chains A {VAR['chains']['A']}, B {VAR['chains']['B']}; pairs {VAR['n_pairs']}. {rule} Empty cells are "
                        "unavailable values, never zero; the counts (`*_n`) say how many chain pairs contributed.")
        chan = st.selectbox("Channel", list(res.variation["channel"].unique()), key="var_channel")
        d = res.variation[res.variation["channel"] == chan]
        st.dataframe(d.drop(columns=["unit"]), hide_index=True, key="var_table")
        if res.variation_figure:
            st.image(str(res.variation_figure), caption="chain_variation.png as written by structpair (displacement channel)")
        if res.chains is not None:
            st.markdown("**Chains taking part** (missing residues and alternate-location decisions per chain):")
            st.dataframe(res.chains, hide_index=True, key="var_chains")
with t_exc:
    if not len(res.exclusions):
        st.success("No window inside the mapped range was excluded.")
    else:
        st.markdown("One row per excluded window, structure and cause, from `exclusions.csv` (adapter break records, "
                    "residue status and feature reason codes). The explanation column describes the cause type; the "
                    "detail column is the record itself. In a break record, 'output residues' are row indices (the position column), "
                    "not residue numbers; the reading column translates them using residues.csv.")
        st.dataframe(res.exclusions[["window_size", "position", "window_A", "structure", "cause", "location",
                                     "explanation", "detail", "detail_reading", "source"]], hide_index=True, key="exc_table")
with t_res:
    mm = res.residues[res.residues["resname_mismatch"]]
    miss = res.residues[(res.residues["A_status"] != "present") | (res.residues["B_status"] != "present")]
    st.markdown(f"Correspondence: **{res.manifest['correspondence']['form']}**, {len(res.residues)} positions. "
                f"Missing in A or B: {len(miss)}. Residue-name mismatches: {len(mm)}. Labels are `chain:number"
                "[insertion code]` exactly as in the input files.")
    st.dataframe(res.residues, hide_index=True, key="res_table")
with t_prov:
    m = res.manifest
    st.markdown(f"**Statement recorded with the result:** {m['statement']}")
    st.table({"Result directory": str(res.path), "Created (UTC)": m["created_utc"], "Spec SHA-256": m["spec"]["sha256"],
              "Command": m["command"], "structpair version": m["software"]["structpair_version"],
              "structpair module SHA-256": m["software"]["module_sha256"], "Python / NumPy": f"{m['software']['python']} / {m['software']['numpy']}",
              "QCP": f"{m['software']['qcp']['function']} (MDAnalysis {m['software']['qcp']['MDAnalysis']})",
              "Policies": ", ".join(f"{k}={v}" for k, v in m["policies"].items())}, border="horizontal")
    st.dataframe(pd.DataFrame([{"structure": k, "label": v["label"], "file": Path(v["path"]).name, "chain": v["chain"],
                                "chains": ", ".join(v.get("chains") or [v["chain"]]), "altloc policy": v["altloc_policy"], "SHA-256": v["sha256"],
                                "core native extension SHA-256": (v.get("core", {}).get("native_extension") or {}).get("sha256"),
                                "MDAnalysis": v.get("adapter_software", {}).get("MDAnalysis")}
                               for k, v in m["inputs"].items()]), hide_index=True, key="prov_inputs")
    st.dataframe(pd.DataFrame([{"output": f, "matches manifest SHA-256": ok} for f, ok in res.integrity.items()]),
                 hide_index=True, key="prov_outputs")
    for n in res.notes:
        st.info(n)
with t_rep:
    if res.report_md:
        st.markdown("\n".join(l for l in res.report_md.splitlines() if not l.startswith("![")))
    if res.figure:
        st.image(str(res.figure), caption="figure.png as written by structpair")
