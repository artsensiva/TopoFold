#!/usr/bin/env python3
"""
TopoFold Interactive 3D Web Dashboard
=====================================

Real-time, unsupervised cryptic pocket discovery and 3D molecular viewer
powered by TopoFold SE(3)-invariant differential geometry.

Usage:
    streamlit run apps/streamlit_app.py
"""

import os
import sys
import json
import tempfile
import numpy as np

try:
    import streamlit as st
except ImportError:
    st = None

try:
    import plotly.graph_objects as go
except ImportError:
    go = None

import topofold as tf

BENCHMARK_DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "benchmarks", "data")


def extract_pdb_residue_ids(pdb_str: str):
    """
    Extracts authentic PDB residue sequence numbers from CA records.
    """
    res_ids = []
    for line in pdb_str.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if line[12:16].strip() == "CA":
                try:
                    rid = int(line[22:26].strip())
                    res_ids.append(rid)
                except ValueError:
                    pass
    return res_ids


def render_3dmol_viewer(pdb_str: str, highlight_range=None, width=700, height=480):
    """
    Renders an interactive 3D molecular structure using 3Dmol.js embedded via HTML.
    Layered representation:
      - Scaffold: Semi-transparent cyan cartoon
      - Cryptic Pocket: Vivid orange cartoon + stick side-chains (radius 0.22)
    """
    hl_start, hl_end = highlight_range if highlight_range else (-1, -1)
    escaped_pdb = json.dumps(pdb_str)

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="utf-8">
      <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
      <style>
        body {{
          margin: 0;
          padding: 0;
          overflow: hidden;
          background-color: #0e1117;
          font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
        }}
        #container {{
          width: {width}px;
          height: {height}px;
          position: relative;
        }}
        #legend {{
          position: absolute;
          top: 10px;
          left: 10px;
          background: rgba(14, 17, 23, 0.85);
          padding: 8px 12px;
          border-radius: 6px;
          font-size: 12px;
          color: #e0e0e0;
          border: 1px solid #30363d;
          z-index: 10;
        }}
        .legend-item {{
          display: flex;
          align-items: center;
          margin-bottom: 4px;
        }}
        .legend-color {{
          width: 12px;
          height: 12px;
          border-radius: 3px;
          margin-right: 8px;
        }}
      </style>
    </head>
    <body>
      <div id="container">
        <div id="legend">
          <div class="legend-item">
            <div class="legend-color" style="background: #4a90e2;"></div>
            <span>Rigid Scaffold (SE(3) Invariant)</span>
          </div>
          <div class="legend-item">
            <div class="legend-color" style="background: #e67e22;"></div>
            <span>Discovered Pocket (Res {hl_start}..{hl_end})</span>
          </div>
        </div>
      </div>
      <script>
        let viewer = $3Dmol.createViewer(document.getElementById('container'), {{
          backgroundColor: '#0e1117'
        }});
        let pdbData = {escaped_pdb};
        viewer.addModel(pdbData, "pdb");
        
        // Base structure style: Semi-transparent cyan cartoon
        viewer.setStyle({{}}, {{cartoon: {{color: '#4a90e2', opacity: 0.80}}}});
        
        // Highlight discovered cryptic pocket: Cartoon + stick
        let hlStart = {hl_start};
        let hlEnd = {hl_end};
        if (hlStart > 0 && hlEnd >= hlStart) {{
          let pocketResis = [];
          for (let r = hlStart; r <= hlEnd; r++) {{
            pocketResis.push(r);
          }}
          viewer.addStyle({{resi: pocketResis}}, {{
            cartoon: {{color: '#e67e22', opacity: 1.0}},
            stick: {{colorscheme: 'orangeCarbon', radius: 0.22}}
          }});
        }}
        
        viewer.zoomTo();
        viewer.render();
        viewer.spin('y', 0.4);
      </script>
    </body>
    </html>
    """
    return html_content


def load_dataset(preset_name: str, uploaded_pdb=None, uploaded_dcd=None):
    """
    Loads dataset and persists it into session state.
    """
    if preset_name == "Human Abl1 Kinase DFG Flip (2GQG / 1IEP)":
        pdb_path = os.path.join(BENCHMARK_DATA_DIR, "abl_reference.pdb")
        dcd_path = os.path.join(BENCHMARK_DATA_DIR, "abl_dfg_trajectory.dcd")
        if os.path.exists(pdb_path) and os.path.exists(dcd_path):
            with open(pdb_path, "r") as f:
                ref_pdb_str = f.read()
            traj_data = tf.read_dcd(dcd_path)
            res_ids = extract_pdb_residue_ids(ref_pdb_str)
            return {
                "name": preset_name,
                "traj": traj_data,
                "pdb": ref_pdb_str,
                "res_ids": res_ids,
                "desc": f"Human c-Abl1 Kinase Domain (1,500 frames, {traj_data.shape[1]} residues, PDB 225..498)",
            }

    elif preset_name == "BPTI Catalytic P1 Loop Flip (Shaw et al. Science 2010)":
        pdb_path = os.path.join(BENCHMARK_DATA_DIR, "5PTI.pdb")
        dcd_path = os.path.join(BENCHMARK_DATA_DIR, "bpti_equilibrium.dcd")
        if os.path.exists(pdb_path) and os.path.exists(dcd_path):
            with open(pdb_path, "r") as f:
                ref_pdb_str = f.read()
            traj_data = tf.read_dcd(dcd_path)
            res_ids = extract_pdb_residue_ids(ref_pdb_str)
            return {
                "name": preset_name,
                "traj": traj_data,
                "pdb": ref_pdb_str,
                "res_ids": res_ids,
                "desc": f"Bovine Pancreatic Trypsin Inhibitor (2,500 frames, {traj_data.shape[1]} residues, PDB 1..58)",
            }

    elif preset_name == "Synthetic Bistable Gating Ensemble (60 residues)":
        dcd_path = os.path.join(BENCHMARK_DATA_DIR, "bistable_ensemble.dcd")
        pdb_path = os.path.join(BENCHMARK_DATA_DIR, "state_a.pdb")
        if os.path.exists(pdb_path) and os.path.exists(dcd_path):
            with open(pdb_path, "r") as f:
                ref_pdb_str = f.read()
            traj_data = tf.read_dcd(dcd_path)
            res_ids = extract_pdb_residue_ids(ref_pdb_str)
            return {
                "name": preset_name,
                "traj": traj_data,
                "pdb": ref_pdb_str,
                "res_ids": res_ids,
                "desc": f"Synthetic Bistable Loop (2,000 frames, {traj_data.shape[1]} residues)",
            }

    elif uploaded_pdb and uploaded_dcd:
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_pdb:
            tmp_pdb.write(uploaded_pdb.read())
            tmp_pdb_path = tmp_pdb.name
        with tempfile.NamedTemporaryFile(suffix=".dcd", delete=False) as tmp_dcd:
            tmp_dcd.write(uploaded_dcd.read())
            tmp_dcd_path = tmp_dcd.name

        with open(tmp_pdb_path, "r") as f:
            ref_pdb_str = f.read()
        traj_data = tf.read_dcd(tmp_dcd_path)
        res_ids = extract_pdb_residue_ids(ref_pdb_str)
        return {
            "name": "Custom Upload",
            "traj": traj_data,
            "pdb": ref_pdb_str,
            "res_ids": res_ids,
            "desc": f"Custom Upload ({traj_data.shape[0]} frames, {traj_data.shape[1]} residues)",
        }

    return None


def run_dashboard():
    if st is None:
        print("Streamlit is not installed in the environment. Install via `pip install streamlit plotly`.")
        return

    st.set_page_config(
        page_title="TopoFold: Autonomous Cryptic Pocket Explorer",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    st.title("🧬 TopoFold: Autonomous Cryptic Pocket & Allosteric Explorer")
    st.markdown(
        "*SE(3)-Invariant Discrete Differential Geometry & Vantage-Point Metric Search for Protein Dynamics*"
    )
    st.markdown("---")

    # Initialize Session State
    if "current_dataset" not in st.session_state:
        st.session_state["current_dataset"] = None
    if "scan_results" not in st.session_state:
        st.session_state["scan_results"] = None
    if "selected_pocket_idx" not in st.session_state:
        st.session_state["selected_pocket_idx"] = 0

    # Sidebar: Data Ingestion & Configuration
    st.sidebar.header("📁 Trajectory Ingestion")
    preset_choice = st.sidebar.selectbox(
        "Choose Benchmark Dataset:",
        [
            "Human Abl1 Kinase DFG Flip (2GQG / 1IEP)",
            "BPTI Catalytic P1 Loop Flip (Shaw et al. Science 2010)",
            "Synthetic Bistable Gating Ensemble (60 residues)",
            "Custom Upload (PDB + DCD)",
        ],
        index=0,
    )

    uploaded_pdb = None
    uploaded_dcd = None
    if preset_choice == "Custom Upload (PDB + DCD)":
        uploaded_pdb = st.sidebar.file_uploader("Upload Structure (PDB)", type=["pdb"])
        uploaded_dcd = st.sidebar.file_uploader("Upload Trajectory (DCD)", type=["dcd"])

    # Load dataset if changed or not loaded yet
    dataset_needs_load = (
        st.session_state["current_dataset"] is None
        or st.session_state["current_dataset"]["name"] != preset_choice
    )

    if dataset_needs_load:
        loaded = load_dataset(preset_choice, uploaded_pdb, uploaded_dcd)
        if loaded is not None:
            st.session_state["current_dataset"] = loaded
            st.session_state["scan_results"] = None  # Invalidate previous scan
            st.session_state["selected_pocket_idx"] = 0

    curr_ds = st.session_state["current_dataset"]

    if curr_ds is not None:
        st.sidebar.success(f"✓ {curr_ds['desc']}")
    else:
        st.sidebar.info("Awaiting trajectory upload...")

    st.sidebar.header("⚙️ Geometry Engine Parameters")
    window_size = st.sidebar.slider("Sliding Window Size (W residues):", min_value=4, max_value=16, value=8)
    bc_threshold = st.sidebar.slider("Sarle's Bimodality Threshold (BC):", min_value=0.50, max_value=0.98, value=0.60, step=0.02)
    use_cbeta = st.sidebar.checkbox("Enable C-beta Ribbon Geometry (Side-Chain Gating)", value=True)

    if curr_ds is None:
        st.info("👈 Select a pre-packaged benchmark dataset or upload a PDB+DCD trajectory to begin.")
        return

    traj_data = curr_ds["traj"]
    ref_pdb_str = curr_ds["pdb"]
    res_ids = curr_ds["res_ids"]
    n_residues = traj_data.shape[1]
    if len(res_ids) != n_residues:
        res_ids = list(range(1, n_residues + 1))

    # Action Button & Status
    col_btn, col_stats = st.columns([1, 3])
    with col_btn:
        run_scan = st.button("🚀 Run Autonomous Pocket Scan", type="primary")

    if run_scan:
        with st.spinner("Analyzing intrinsic curve invariants & Sarle's bimodality moments in Rust core..."):
            scores, bc_tau, bc_kappa = tf.compute_bimodality_profile(traj_data, window_size=window_size)
            candidates = tf.scan_cryptic_pockets(traj_data, window_size=window_size, bc_threshold=bc_threshold)
            st.session_state["scan_results"] = {
                "scores": scores,
                "bc_tau": bc_tau,
                "bc_kappa": bc_kappa,
                "candidates": candidates,
                "window_size": window_size,
            }
            st.session_state["selected_pocket_idx"] = 0

    scan_res = st.session_state["scan_results"]

    if scan_res is not None:
        scores = scan_res["scores"]
        bc_tau = scan_res["bc_tau"]
        bc_kappa = scan_res["bc_kappa"]
        candidates = scan_res["candidates"]
        w_size = scan_res["window_size"]

        with col_stats:
            st.success(
                f"**Scan Complete**: Evaluated {len(scores)} sliding windows across {traj_data.shape[0]:,} frames. "
                f"Discovered **{len(candidates)}** candidate cryptic pocket / functional loop clusters."
            )

        # Layout: Left column = Profile Plot & Table, Right column = 3D Viewer & Selection
        col_plot, col_view = st.columns([1.15, 1.0])

        with col_plot:
            st.subheader("📈 Sequence-Wide Bimodality Profile")
            mid_offset = w_size // 2
            x_res_pdb = [res_ids[min(i + mid_offset, n_residues - 1)] for i in range(len(scores))]

            if go is not None:
                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=x_res_pdb, y=scores, mode="lines",
                    name="Composite max(BC_tau, BC_kappa)",
                    line=dict(color="#2ecc71", width=3)
                ))
                fig.add_trace(go.Scatter(
                    x=x_res_pdb, y=bc_tau, mode="lines",
                    name="Discrete Torsion (BC_tau)",
                    line=dict(color="#3498db", width=1.5, dash="dash")
                ))
                fig.add_trace(go.Scatter(
                    x=x_res_pdb, y=bc_kappa, mode="lines",
                    name="Discrete Curvature (BC_kappa)",
                    line=dict(color="#9b59b6", width=1.5, dash="dot")
                ))
                fig.add_hline(
                    y=0.555, line_dash="dash", line_color="#e74c3c",
                    annotation_text="Unimodal Threshold (BC=0.555)", annotation_position="bottom right"
                )
                fig.update_layout(
                    xaxis_title="PDB Residue Sequence Number",
                    yaxis_title="Sarle's Bimodality Coefficient (BC)",
                    yaxis_range=[-0.02, 1.05],
                    template="plotly_dark",
                    margin=dict(l=40, r=20, t=30, b=40),
                    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
                )
                st.plotly_chart(fig)
            else:
                st.line_chart(scores)

            st.subheader("🏆 Discovered Cryptic Pocket Clusters")
            if candidates:
                cand_data = []
                for r, (s_res, e_res, sc) in enumerate(candidates, 1):
                    pdb_start = res_ids[min(s_res, n_residues - 1)]
                    pdb_end = res_ids[min(e_res, n_residues - 1)]
                    cand_data.append({
                        "Rank": f"#{r}",
                        "PDB Residues": f"Res {pdb_start} .. {pdb_end}",
                        "Length": e_res - s_res + 1,
                        "Peak Bimodality (BC)": f"{sc:.4f}",
                        "0-Based Index": f"[{s_res} .. {e_res}]",
                    })
                st.dataframe(cand_data)
            else:
                st.info("No clusters detected above current threshold. Lower threshold to discover subtle transitions.")

        with col_view:
            st.subheader("🔬 3D Molecular Structure & Pocket Highlight")
            selected_pdb_range = (res_ids[0], res_ids[min(10, n_residues - 1)])

            if candidates:
                cand_options = []
                for i, (s_res, e_res, sc) in enumerate(candidates):
                    pdb_start = res_ids[min(s_res, n_residues - 1)]
                    pdb_end = res_ids[min(e_res, n_residues - 1)]
                    cand_options.append(f"Rank #{i+1}: Res {pdb_start}..{pdb_end} (BC = {sc:.4f})")

                selected_idx = st.selectbox(
                    "Select Detected Pocket to Highlight:",
                    range(len(candidates)),
                    format_func=lambda i: cand_options[i],
                    key="pocket_selectbox",
                )
                st.session_state["selected_pocket_idx"] = selected_idx
                sel_cand = candidates[selected_idx]
                selected_pdb_range = (
                    res_ids[min(sel_cand[0], n_residues - 1)],
                    res_ids[min(sel_cand[1], n_residues - 1)],
                )

            st.markdown(
                f"**Active Highlight**: PDB Residues `{selected_pdb_range[0]}..{selected_pdb_range[1]}` "
                f"*(Orange Licorice & Sticks)* on scaffold *(Cyan Cartoon)*."
            )

            html_viewer = render_3dmol_viewer(
                ref_pdb_str,
                highlight_range=selected_pdb_range,
                width=650,
                height=480,
            )
            st.components.v1.html(html_viewer, height=490)
    else:
        st.info("💡 Click **'Run Autonomous Pocket Scan'** above to detect cryptic pockets and mobile functional loops.")


if __name__ == "__main__":
    if st is None:
        print("TopoFold Interactive 3D Dashboard")
        print("Please install Streamlit to run GUI: pip install streamlit plotly")
    else:
        run_dashboard()
