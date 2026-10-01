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


def render_3dmol_viewer(pdb_str, highlight_range=None, width=700, height=450):
    """
    Renders an interactive 3D molecular structure using 3Dmol.js embedded via HTML.
    """
    hl_start, hl_end = highlight_range if highlight_range else (-1, -1)

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
      <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    </head>
    <body style="margin: 0; padding: 0; overflow: hidden; background-color: #0e1117;">
      <div id="container" style="width: {width}px; height: {height}px; position: relative;"></div>
      <script>
        let viewer = $3Dmol.createViewer(document.getElementById('container'), {{
          backgroundColor: '#0e1117'
        }});
        let pdbData = `{pdb_str}`;
        viewer.addModel(pdbData, "pdb");
        
        // Base structure style: Semi-transparent cyan cartoon
        viewer.setStyle({{}}, {{cartoon: {{color: '#4a90e2', opacity: 0.85}}}});
        
        // Highlight discovered cryptic pocket
        if ({hl_start} > 0 && {hl_end} >= {hl_start}) {{
          viewer.setStyle({{resi: Array.from({{length: {hl_end} - {hl_start} + 1}}, (_, i) => i + {hl_start})}},
            {{cartoon: {{color: '#e67e22', opacity: 1.0}}, stick: {{colorscheme: 'orangeCarbon', radius: 0.18}}}}
          );
        }}
        
        viewer.zoomTo();
        viewer.render();
        viewer.spin('y', 0.5);
      </script>
    </body>
    </html>
    """
    return html_content


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

    # Sidebar: Data Ingestion & Configuration
    st.sidebar.header("📁 Trajectory Ingestion")
    preset = st.sidebar.selectbox(
        "Choose Pre-Packaged Benchmark Dataset:",
        [
            "BPTI Catalytic P1 Loop Flip (Shaw et al. Science 2010)",
            "Synthetic Bistable Gating Ensemble (60 residues)",
            "Custom Upload (PDB + DCD)",
        ],
    )

    window_size = st.sidebar.slider("Sliding Window Size (W residues):", min_value=4, max_value=16, value=8)
    bc_threshold = st.sidebar.slider("Sarle's Bimodality Threshold (BC):", min_value=0.50, max_value=0.95, value=0.60, step=0.05)
    use_cbeta = st.sidebar.checkbox("Enable C-beta Ribbon Geometry (Side-Chain Gating)", value=True)

    pdb_file = None
    traj_data = None
    ref_pdb_str = ""

    if preset == "BPTI Catalytic P1 Loop Flip (Shaw et al. Science 2010)":
        pdb_path = os.path.join(BENCHMARK_DATA_DIR, "5PTI.pdb")
        dcd_path = os.path.join(BENCHMARK_DATA_DIR, "bpti_equilibrium.dcd")
        if os.path.exists(pdb_path) and os.path.exists(dcd_path):
            with open(pdb_path, "r") as f:
                ref_pdb_str = f.read()
            traj_data = tf.read_dcd(dcd_path)
            st.sidebar.success("Loaded BPTI MD Trajectory (2,500 frames, 58 residues)")
        else:
            st.sidebar.error("BPTI trajectory files not found in benchmarks/data.")

    elif preset == "Synthetic Bistable Gating Ensemble (60 residues)":
        dcd_path = os.path.join(BENCHMARK_DATA_DIR, "bistable_ensemble.dcd")
        pdb_path = os.path.join(BENCHMARK_DATA_DIR, "state_a.pdb")
        if os.path.exists(pdb_path) and os.path.exists(dcd_path):
            with open(pdb_path, "r") as f:
                ref_pdb_str = f.read()
            traj_data = tf.read_dcd(dcd_path)
            st.sidebar.success("Loaded Synthetic Ensemble (2,000 frames, 60 residues)")
        else:
            st.sidebar.warning("Run benchmarks/generate_bistable_trajectory.py to generate dataset.")

    else:
        uploaded_pdb = st.sidebar.file_uploader("Upload Structure (PDB)", type=["pdb"])
        uploaded_dcd = st.sidebar.file_uploader("Upload Trajectory (DCD)", type=["dcd"])
        if uploaded_pdb and uploaded_dcd:
            with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp_pdb:
                tmp_pdb.write(uploaded_pdb.read())
                tmp_pdb_path = tmp_pdb.name
            with tempfile.NamedTemporaryFile(suffix=".dcd", delete=False) as tmp_dcd:
                tmp_dcd.write(uploaded_dcd.read())
                tmp_dcd_path = tmp_dcd.name

            with open(tmp_pdb_path, "r") as f:
                ref_pdb_str = f.read()
            traj_data = tf.read_dcd(tmp_dcd_path)
            st.sidebar.success(f"Ingested Trajectory: {traj_data.shape[0]} frames, {traj_data.shape[1]} residues.")

    if traj_data is None:
        st.info("👈 Select a pre-packaged benchmark dataset or upload a PDB+DCD trajectory to begin.")
        return

    # Action Button
    col_btn, col_stats = st.columns([1, 3])
    with col_btn:
        run_scan = st.button("🚀 Run Autonomous Pocket Scan", type="primary")

    if run_scan or "scan_results" in st.session_state:
        if run_scan:
            with st.spinner("Analyzing intrinsic curve invariants & Sarle's bimodality moments..."):
                scores, bc_tau, bc_kappa = tf.compute_bimodality_profile(traj_data, window_size=window_size)
                candidates = tf.scan_cryptic_pockets(traj_data, window_size=window_size, bc_threshold=bc_threshold)
                st.session_state["scan_results"] = (scores, bc_tau, bc_kappa, candidates)
        else:
            scores, bc_tau, bc_kappa, candidates = st.session_state["scan_results"]

        with col_stats:
            st.success(
                f"**Scan Complete**: Evaluated {len(scores)} sliding windows across {traj_data.shape[0]:,} frames. "
                f"Discovered **{len(candidates)}** candidate cryptic pocket / functional loop clusters."
            )

        # Layout: Left column = Profile Plot, Right column = 3D Viewer & Selection
        col_plot, col_view = st.columns([1.1, 1.0])

        with col_plot:
            st.subheader("📈 Sequence-Wide Bimodality Profile")
            x_res = np.arange(len(scores)) + (window_size // 2) + 1

            if go is not None:
                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=x_res, y=scores, mode="lines",
                    name="Composite Score max(BC_tau, BC_kappa)",
                    line=dict(color="#27ae60", width=3)
                ))
                fig.add_trace(go.Scatter(
                    x=x_res, y=bc_tau, mode="lines",
                    name="Discrete Torsion (BC_tau)",
                    line=dict(color="#2980b9", width=1.5, dash="dash")
                ))
                fig.add_trace(go.Scatter(
                    x=x_res, y=bc_kappa, mode="lines",
                    name="Discrete Curvature (BC_kappa)",
                    line=dict(color="#8e44ad", width=1.5, dash="dot")
                ))
                fig.add_hline(
                    y=0.555, line_dash="dash", line_color="#e74c3c",
                    annotation_text="Unimodal Threshold (BC=0.555)", annotation_position="bottom right"
                )
                fig.update_layout(
                    xaxis_title="Residue Index",
                    yaxis_title="Sarle's Bimodality Coefficient",
                    yaxis_range=[-0.02, 1.05],
                    template="plotly_dark",
                    margin=dict(l=40, r=20, t=30, b=40),
                    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
                )
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.line_chart(scores)

            st.subheader("🏆 Discovered Cryptic Pocket Clusters")
            if candidates:
                cand_data = []
                for r, (s_res, e_res, sc) in enumerate(candidates, 1):
                    cand_data.append({
                        "Rank": f"#{r}",
                        "PDB Residues": f"Res {s_res + 1} .. {e_res + 1}",
                        "Residue Count": e_res - s_res + 1,
                        "Peak Bimodality": f"{sc:.4f}",
                    })
                st.dataframe(cand_data, use_container_width=True)
            else:
                st.info("No clusters detected above current threshold. Lower threshold to see subtle transitions.")

        with col_view:
            st.subheader("🔬 3D Molecular Structure & Pocket Highlight")
            selected_range = (10, 18)
            if candidates:
                selected_cand_idx = st.selectbox(
                    "Highlight Detected Pocket on 3D Scaffold:",
                    range(len(candidates)),
                    format_func=lambda i: f"Rank #{i+1}: Residues {candidates[i][0]+1}..{candidates[i][1]+1} (Peak BC = {candidates[i][2]:.4f})"
                )
                sel_cand = candidates[selected_cand_idx]
                selected_range = (sel_cand[0] + 1, sel_cand[1] + 1)

            st.markdown(
                f"**Active Highlight**: PDB Residues `{selected_range[0]}..{selected_range[1]}` "
                f"*(Orange Licorice)* amidst scaffold *(Cyan Cartoon)*."
            )

            html_viewer = render_3dmol_viewer(ref_pdb_str, highlight_range=selected_range, width=650, height=440)
            st.components.v1.html(html_viewer, height=450)


if __name__ == "__main__":
    if st is None:
        print("TopoFold Interactive 3D Dashboard")
        print("Please install Streamlit to run GUI: pip install streamlit plotly")
    else:
        run_dashboard()
