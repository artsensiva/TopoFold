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
import html
import json
import tempfile
from pathlib import Path
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

REPO_ROOT = Path(__file__).resolve().parent.parent
BENCHMARK_DATA_DIR = REPO_ROOT / "benchmarks" / "data"


def extract_pdb_residue_ids(pdb_str: str) -> list[int]:
    """
    Extracts authentic PDB residue sequence numbers from CA (protein) or C4' (RNA) records.
    """
    res_ids = []
    seen = set()
    for line in pdb_str.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_name = line[12:16].strip()
            if atom_name in ("CA", "C4'"):
                try:
                    rid = int(line[22:26].strip())
                    if rid not in seen:
                        seen.add(rid)
                        res_ids.append(rid)
                except ValueError:
                    pass
    return res_ids


def render_3dmol_viewer(
    pdb_str: str,
    highlight_range: tuple[int, int] | None = None,
    res_count: int = 0,
    height: int = 520,
    is_rna: bool = False,
    is_idp: bool = False,
) -> str:
    """
    Renders an interactive 3D molecular structure using 3Dmol.js embedded via HTML.

    Features:
      - Robust multi-CDN loader with sequential fallbacks (cdnjs -> Pitt -> jsdelivr -> 3Dmol.org).
      - HTML-escaped raw PDB text embedding via hidden textarea (immune to JS quoting/newline issues).
      - Scaffold: semi-transparent light grey cartoon / RNA ribbon / disordered wireframe + trace.
      - Discovered Pocket / Hinge / Motif: vivid orange/red (#FF5722) cartoon + sticks / spheres.
      - Responsive ResizeObserver & staggered redraws (guarantees no blank grey canvas on load/resize).
      - Explicit dimensions (width: 100%, height: 520px) with WebGL initialization error boundary.
      - Status text indicator: 'Status: WebGL rendered N residues'.
    """
    hl_start, hl_end = highlight_range if highlight_range else (-1, -1)
    escaped_pdb_html = html.escape(pdb_str)
    pocket_label = f"Res {hl_start}..{hl_end}" if (hl_start > 0 and hl_end >= hl_start) else "None"
    if is_idp:
        scaffold_label = "Disordered Wireframe/Trace"
        highlight_label = "Transient Motif"
    elif is_rna:
        scaffold_label = "RNA Ribbon"
        highlight_label = "Switching Hinge"
    else:
        scaffold_label = "Cartoon"
        highlight_label = "Pocket"

    html_content = f"""<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <style>
    * {{
      box-sizing: border-box;
    }}
    body {{
      margin: 0;
      padding: 0;
      overflow: hidden;
      background-color: #0e1117;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
    }}
    #viewer-container {{
      width: 100%;
      height: {height}px;
      min-height: {height}px;
      position: relative;
      background-color: #111318;
      border-radius: 8px;
      overflow: hidden;
      border: 1px solid #30363d;
      box-shadow: 0 4px 12px rgba(0, 0, 0, 0.35);
    }}
    #viewport {{
      width: 100%;
      height: 100%;
      position: absolute;
      top: 0;
      left: 0;
    }}
    #legend {{
      position: absolute;
      top: 12px;
      left: 12px;
      background: rgba(14, 17, 23, 0.88);
      backdrop-filter: blur(4px);
      padding: 8px 14px;
      border-radius: 6px;
      font-size: 12px;
      color: #e0e0e0;
      border: 1px solid #30363d;
      z-index: 10;
      pointer-events: none;
      user-select: none;
    }}
    .legend-item {{
      display: flex;
      align-items: center;
      margin-bottom: 5px;
    }}
    .legend-item:last-child {{
      margin-bottom: 0;
    }}
    .legend-color {{
      width: 14px;
      height: 14px;
      border-radius: 3px;
      margin-right: 8px;
      flex-shrink: 0;
    }}
    #status-bar {{
      margin-top: 8px;
      font-size: 12px;
      color: #8b949e;
      font-family: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Monaco, Consolas, monospace;
      padding: 0 2px;
      display: flex;
      justify-content: space-between;
      align-items: center;
    }}
    .status-badge {{
      display: inline-block;
      width: 8px;
      height: 8px;
      border-radius: 50%;
      background-color: #2ea043;
      margin-right: 6px;
    }}
    #loading-overlay {{
      position: absolute;
      top: 0;
      left: 0;
      width: 100%;
      height: 100%;
      display: flex;
      flex-direction: column;
      align-items: center;
      justify-content: center;
      background-color: #111318;
      z-index: 20;
      transition: opacity 0.3s ease;
    }}
    .spinner {{
      width: 32px;
      height: 32px;
      border: 3px solid #21262d;
      border-top: 3px solid #58a6ff;
      border-radius: 50%;
      animation: spin 0.8s linear infinite;
      margin-bottom: 12px;
    }}
    @keyframes spin {{
      0% {{ transform: rotate(0deg); }}
      100% {{ transform: rotate(360deg); }}
    }}
  </style>
  <!-- Primary CDN: Cloudflare cdnjs -->
  <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"></script>
</head>
<body>
  <!-- Raw PDB text embedded inside hidden textarea to avoid JS quoting issues -->
  <textarea id="pdb-data" style="display: none;">{escaped_pdb_html}</textarea>

  <div id="viewer-container">
    <div id="loading-overlay">
      <div class="spinner"></div>
      <div id="loading-text" style="color: #8b949e; font-size: 13px; font-family: ui-monospace, SFMono-Regular, monospace;">Initializing 3Dmol Canvas...</div>
    </div>
    <div id="viewport"></div>
    <div id="legend">
      <div class="legend-item">
        <div class="legend-color" style="background: #b0bec5; opacity: 0.6; border: 1px solid #78909c;"></div>
        <span>Scaffold ({scaffold_label}, opacity: 0.6)</span>
      </div>
      <div class="legend-item">
        <div class="legend-color" style="background: #FF5722; border: 1px solid #e64a19;"></div>
        <span>Discovered {highlight_label} ({pocket_label})</span>
      </div>
    </div>
  </div>

  <div id="status-bar">
    <span id="status-text"><span class="status-badge" id="status-dot"></span>Status: Initializing WebGL viewer...</span>
    <span style="color: #6e7681;">3Dmol.js • Drag to rotate • Scroll to zoom</span>
  </div>

  <script>
    const FALLBACK_CDNS = [
      "https://3dmol.csb.pitt.edu/build/3Dmol-min.js",
      "https://cdn.jsdelivr.net/npm/3dmol@2.4.2/build/3Dmol-min.js",
      "https://3Dmol.org/build/3Dmol-min.js"
    ];

    function setStatus(msg, isError) {{
      const textEl = document.getElementById('status-text');
      const dotEl = document.getElementById('status-dot');
      if (textEl) {{
        textEl.innerHTML = '<span class="status-badge" style="background-color: ' + (isError ? '#cf222e' : '#2ea043') + ';"></span>' + msg;
      }}
    }}

    function hideLoader() {{
      const loader = document.getElementById('loading-overlay');
      if (loader) {{
        loader.style.opacity = '0';
        setTimeout(() => {{ loader.style.display = 'none'; }}, 300);
      }}
    }}

    function initViewer() {{
      try {{
        const container = document.getElementById('viewport');
        if (!container) return;

        if (typeof $3Dmol === 'undefined') {{
          setStatus("Status: Error - 3Dmol library unavailable", true);
          const loadText = document.getElementById('loading-text');
          if (loadText) loadText.innerHTML = '<span style="color: #cf222e;">Error: 3Dmol library failed to load</span>';
          return;
        }}

        // Initialize WebGL Viewer with explicit container
        const viewer = $3Dmol.createViewer(container, {{
          backgroundColor: '#111318'
        }});

        if (!viewer) {{
          setStatus("Status: Error - WebGL context creation failed", true);
          const loadText = document.getElementById('loading-text');
          if (loadText) loadText.innerHTML = '<span style="color: #cf222e;">Error: WebGL context creation failed</span>';
          return;
        }}

        // Ingest raw PDB text from hidden textarea
        const pdbData = document.getElementById('pdb-data').value;
        if (!pdbData || pdbData.trim().length === 0) {{
          setStatus("Status: Warning - No PDB coordinates found", true);
          const loadText = document.getElementById('loading-text');
          if (loadText) loadText.innerHTML = '<span style="color: #d29922;">Warning: No PDB coordinates found</span>';
          return;
        }}

        viewer.addModel(pdbData, "pdb");

        const isIdp = {"true" if is_idp else "false"};
        const isRna = {"true" if is_rna else "false"};

        // 1. Macromolecule-aware scaffold representation
        if (isIdp) {{
          // Intrinsically Disordered Protein: wireframe and backbone trace (no false cartoon tags)
          viewer.setStyle({{}}, {{
            trace: {{ color: '#9E9E9E', radius: 0.2 }},
            line: {{ color: '#BDBDBD', linewidth: 1.5 }}
          }});
        }} else if (isRna) {{
          // RNA Riboswitch: formatted specifically for nucleic P-C4' chains
          viewer.setStyle({{}}, {{
            cartoon: {{ color: 'spectrum', style: 'oval', ribbon: true }},
            stick: {{ radius: 0.12 }}
          }});
        }} else {{
          // Standard Globular Protein: semi-transparent light grey cartoon
          viewer.setStyle({{}}, {{
            cartoon: {{ color: 'lightgrey', opacity: 0.6 }}
          }});
        }}

        // 2. Discovered Pocket / Motif / Switching Hinge
        const hlStart = {hl_start};
        const hlEnd = {hl_end};
        if (hlStart > 0 && hlEnd >= hlStart) {{
          const pocketResis = [];
          for (let r = hlStart; r <= hlEnd; r++) {{
            pocketResis.push(r);
          }}
          if (isIdp) {{
            // Highlighted nucleated NACore motif
            viewer.addStyle(
              {{ resi: pocketResis }},
              {{
                cartoon: {{ color: '#FF5722' }},
                stick: {{ color: '#FF5722', radius: 0.25 }},
                sphere: {{ color: '#FF5722', radius: 0.4 }}
              }}
            );
          }} else if (isRna) {{
            // Highlighted P1 switching stem
            viewer.addStyle(
              {{ resi: pocketResis }},
              {{
                cartoon: {{ color: '#FF5722', style: 'oval', ribbon: true }},
                stick: {{ color: '#FF5722', radius: 0.25 }},
                sphere: {{ color: '#FF5722', radius: 0.4 }}
              }}
            );
          }} else {{
            // Highlighted standard pocket / loop
            viewer.addStyle(
              {{ resi: pocketResis }},
              {{
                cartoon: {{ color: '#FF5722', opacity: 1.0 }},
                stick: {{ color: '#FF5722', radius: 0.2 }}
              }}
            );
          }}
        }}

        viewer.zoomTo();
        viewer.render();
        viewer.spin('y', 0.3);
        hideLoader();

        // Staggered resize and observer to prevent empty canvas issues on layout load
        function doResize() {{
          if (viewer) {{
            viewer.resize();
            viewer.render();
          }}
        }}

        window.addEventListener('resize', doResize);
        if (window.ResizeObserver) {{
          const ro = new ResizeObserver(() => {{
            doResize();
          }});
          ro.observe(container);
        }}

        // Staggered resize triggers for iframe/tabs settling
        [50, 150, 350, 700, 1200].forEach(delay => {{
          setTimeout(doResize, delay);
        }});

        // Count unique rendered residues
        const atoms = viewer.selectedAtoms({{}});
        const resSet = new Set();
        for (let i = 0; i < atoms.length; i++) {{
          if (atoms[i].resi !== undefined) {{
            resSet.add(atoms[i].resi);
          }}
        }}
        const nRes = resSet.size > 0 ? resSet.size : {res_count};
        setStatus("Status: WebGL rendered " + nRes + " residues", false);

      }} catch (err) {{
        console.error("3Dmol WebGL Initialization Error:", err);
        setStatus("Status: WebGL error - " + err.message, true);
        const loadText = document.getElementById('loading-text');
        if (loadText) loadText.innerHTML = '<span style="color: #cf222e;">WebGL Error: ' + err.message + '</span>';
      }}
    }}

    // Check if primary CDN loaded; if not, try fallbacks sequentially
    if (typeof $3Dmol !== 'undefined') {{
      initViewer();
    }} else {{
      let cdnIdx = 0;
      function tryNextCDN() {{
        if (cdnIdx >= FALLBACK_CDNS.length) {{
          setStatus("Status: Error - Failed to load 3Dmol.js from all CDNs", true);
          const loadText = document.getElementById('loading-text');
          if (loadText) loadText.innerHTML = '<span style="color: #cf222e;">Error: Failed to load 3Dmol.js</span>';
          return;
        }}
        const nextUrl = FALLBACK_CDNS[cdnIdx++];
        setStatus("Status: Loading 3Dmol.js from backup CDN (" + cdnIdx + "/" + FALLBACK_CDNS.length + ")...", false);
        const loadText = document.getElementById('loading-text');
        if (loadText) loadText.innerText = 'Loading 3Dmol.js from backup CDN (' + cdnIdx + '/' + FALLBACK_CDNS.length + ')...';
        const script = document.createElement('script');
        script.src = nextUrl;
        script.onload = function() {{
          initViewer();
        }};
        script.onerror = function() {{
          tryNextCDN();
        }};
        document.head.appendChild(script);
      }}
      tryNextCDN();
    }}
  </script>
</body>
</html>
"""
    return html_content


def load_dataset(preset_name: str, uploaded_pdb=None, uploaded_dcd=None) -> dict | None:
    """
    Loads dataset and prepares coordinate array and PDB reference string using universal absolute paths.
    """
    if preset_name == "Human Abl1 Kinase DFG Flip (2GQG / 1IEP)":
        pdb_path = BENCHMARK_DATA_DIR / "abl_reference.pdb"
        dcd_path = BENCHMARK_DATA_DIR / "abl_dfg_trajectory.dcd"
        if pdb_path.is_file() and dcd_path.is_file():
            with open(pdb_path, "r") as f:
                ref_pdb_str = f.read()
            traj_data = tf.read_dcd(str(dcd_path))
            res_ids = extract_pdb_residue_ids(ref_pdb_str)
            return {
                "name": preset_name,
                "traj": traj_data,
                "pdb": ref_pdb_str,
                "res_ids": res_ids,
                "desc": f"Human c-Abl1 Kinase Domain (1,500 frames, {traj_data.shape[1]} residues, PDB 225..498)",
            }

    elif preset_name == "Human Alpha-Synuclein IDP (Parkinson's Disease - 140 Residues)":
        pdb_path = BENCHMARK_DATA_DIR / "alphasynuclein_reference.pdb"
        dcd_path = BENCHMARK_DATA_DIR / "alphasynuclein_ensemble.dcd"
        if pdb_path.is_file() and dcd_path.is_file():
            with open(pdb_path, "r") as f:
                ref_pdb_str = f.read()
            traj_data = tf.read_dcd(str(dcd_path))
            res_ids = extract_pdb_residue_ids(ref_pdb_str)
            return {
                "name": preset_name,
                "traj": traj_data,
                "pdb": ref_pdb_str,
                "res_ids": res_ids,
                "desc": f"Human Alpha-Synuclein IDP (2,000 frames, {traj_data.shape[1]} residues, PDB 1..140)",
                "is_idp": True,
            }

    elif preset_name == "Adenine Riboswitch (RNA 3D Dynamics - PDB 1Y26)":
        pdb_path = BENCHMARK_DATA_DIR / "1Y26.pdb"
        dcd_path = BENCHMARK_DATA_DIR / "rna_riboswitch_trajectory.dcd"
        base_path = BENCHMARK_DATA_DIR / "rna_riboswitch_base.npy"
        if pdb_path.is_file() and dcd_path.is_file():
            with open(pdb_path, "r") as f:
                ref_pdb_str = f.read()
            traj_data = tf.read_dcd(str(dcd_path))
            res_ids = extract_pdb_residue_ids(ref_pdb_str)
            base_traj = np.load(str(base_path)) if base_path.is_file() else None
            return {
                "name": preset_name,
                "traj": traj_data,
                "base_traj": base_traj,
                "pdb": ref_pdb_str,
                "res_ids": res_ids,
                "desc": f"Adenine Riboswitch Aptamer Domain (1,000 frames, {traj_data.shape[1]} nucleotides, PDB 13..83)",
                "is_rna": True,
            }

    elif preset_name == "BPTI Catalytic P1 Loop Flip (Shaw et al. Science 2010)":
        pdb_path = BENCHMARK_DATA_DIR / "5PTI.pdb"
        dcd_path = BENCHMARK_DATA_DIR / "bpti_equilibrium.dcd"
        if pdb_path.is_file() and dcd_path.is_file():
            with open(pdb_path, "r") as f:
                ref_pdb_str = f.read()
            traj_data = tf.read_dcd(str(dcd_path))
            res_ids = extract_pdb_residue_ids(ref_pdb_str)
            return {
                "name": preset_name,
                "traj": traj_data,
                "pdb": ref_pdb_str,
                "res_ids": res_ids,
                "desc": f"Bovine Pancreatic Trypsin Inhibitor (2,500 frames, {traj_data.shape[1]} residues, PDB 1..58)",
            }

    elif preset_name == "Synthetic Bistable Gating Ensemble (60 residues)":
        dcd_path = BENCHMARK_DATA_DIR / "bistable_ensemble.dcd"
        pdb_path = BENCHMARK_DATA_DIR / "state_a.pdb"
        if pdb_path.is_file() and dcd_path.is_file():
            with open(pdb_path, "r") as f:
                ref_pdb_str = f.read()
            traj_data = tf.read_dcd(str(dcd_path))
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
    if "current_pdb_text" not in st.session_state:
        st.session_state["current_pdb_text"] = ""
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
            "Human Alpha-Synuclein IDP (Parkinson's Disease - 140 Residues)",
            "Adenine Riboswitch (RNA 3D Dynamics - PDB 1Y26)",
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

    # Detect dataset load or change
    dataset_needs_load = (
        st.session_state["current_dataset"] is None
        or st.session_state["current_dataset"]["name"] != preset_choice
    )

    if dataset_needs_load:
        loaded = load_dataset(preset_choice, uploaded_pdb, uploaded_dcd)
        if loaded is not None:
            st.session_state["current_dataset"] = loaded
            st.session_state["current_pdb_text"] = loaded["pdb"]
            st.session_state["scan_results"] = None  # Invalidate scan on dataset switch
            st.session_state["selected_pocket_idx"] = 0

    curr_ds = st.session_state["current_dataset"]

    if curr_ds is not None:
        st.sidebar.success(f"✓ {curr_ds['desc']}")
        # Ensure session state has the raw PDB text cached
        if not st.session_state.get("current_pdb_text"):
            st.session_state["current_pdb_text"] = curr_ds["pdb"]
    else:
        st.sidebar.info("Awaiting trajectory upload...")

    is_idp = bool(curr_ds and curr_ds.get("is_idp", False))
    is_rna = bool(curr_ds and curr_ds.get("is_rna", False))

    if is_idp:
        st.sidebar.header("⚙️ IDP Spectral Parameters")
        window_radius = st.sidebar.slider(
            "Sliding Window Half-Width (R residues):",
            min_value=2,
            max_value=8,
            value=4,
            help="Subchain radius R around residue i (total window 2*R + 1 residues).",
        )
        z_threshold = st.sidebar.slider(
            "Topological Z-Score Threshold (Z_c):",
            min_value=0.5,
            max_value=3.0,
            value=1.0,
            step=0.1,
            help="Sequence-normalized Z-score threshold for pre-nucleation motifs.",
        )
    elif is_rna:
        st.sidebar.header("⚙️ RNA Ribbon Parameters")
        window_size = st.sidebar.slider(
            "Sliding Window Size (W nucleotides):",
            min_value=2,
            max_value=8,
            value=3,
            help="Sliding window width for RNA curvature and base dihedral moments.",
        )
        bc_threshold = st.sidebar.slider(
            "Sarle's Bimodality Threshold (BC):",
            min_value=0.50,
            max_value=0.98,
            value=0.70,
            step=0.02,
            help="Bimodality cutoff for identifying bistable switching hinges.",
        )
        st.sidebar.caption("Ribonucleic Ribbon: Phosphorus backbone (P) + Glycosidic Base (C1'→N9/N1)")
    else:
        st.sidebar.header("⚙️ Geometry Engine Parameters")
        window_size = st.sidebar.slider("Sliding Window Size (W residues):", min_value=4, max_value=16, value=8)
        bc_threshold = st.sidebar.slider("Sarle's Bimodality Threshold (BC):", min_value=0.50, max_value=0.98, value=0.60, step=0.02)
        use_cbeta = st.sidebar.checkbox("Enable C-beta Ribbon Geometry (Side-Chain Gating)", value=True)

    if curr_ds is None:
        st.info("👈 Select a pre-packaged benchmark dataset or upload a PDB+DCD trajectory to begin.")
        return

    traj_data = curr_ds["traj"]
    ref_pdb_str = st.session_state.get("current_pdb_text") or curr_ds["pdb"]
    res_ids = curr_ds["res_ids"]
    n_residues = traj_data.shape[1]
    if len(res_ids) != n_residues:
        res_ids = list(range(1, n_residues + 1))

    # Action Button & Status Header
    col_btn, col_stats = st.columns([1, 3])
    with col_btn:
        if is_idp:
            button_label = "⚡ Run Spectral IDP Scan"
        elif is_rna:
            button_label = "⚡ Run Ribonucleic Hinge Scan"
        else:
            button_label = "🚀 Run Autonomous Pocket Scan"
        run_scan = st.button(button_label, type="primary")

    if run_scan:
        if is_idp:
            with st.spinner("Computing Spectral Topological Density S_topo and detecting transient motifs..."):
                s_topo = tf.compute_idp_topological_density(traj_data, window_radius=window_radius)
                mean_dens, var_dens, std_dens, z_scores = tf.compute_idp_density_profile(
                    traj_data, window_radius=window_radius
                )
                motifs = tf.detect_transient_motifs(
                    traj_data, window_size=window_radius * 2, z_threshold=z_threshold
                )
                st.session_state["scan_results"] = {
                    "is_idp": True,
                    "is_rna": False,
                    "s_topo": s_topo,
                    "mean_dens": mean_dens,
                    "var_dens": var_dens,
                    "std_dens": std_dens,
                    "z_scores": z_scores,
                    "motifs": motifs,
                    "window_radius": window_radius,
                    "z_threshold": z_threshold,
                }
                st.session_state["selected_pocket_idx"] = 0
        elif is_rna:
            with st.spinner("Analyzing RNA ribbon invariants & Sarle's bimodality moments in Rust core..."):
                base_traj = curr_ds.get("base_traj")
                scores, bc_theta, bc_kappa, bc_tau = tf.compute_rna_bimodality_profile(
                    traj_data, base_ensemble=base_traj, window_size=window_size
                )
                candidates = tf.scan_rna_switching_hinges(
                    traj_data, base_ensemble=base_traj, window_size=window_size, threshold=bc_threshold
                )
                st.session_state["scan_results"] = {
                    "is_idp": False,
                    "is_rna": True,
                    "scores": scores,
                    "bc_theta": bc_theta,
                    "bc_kappa": bc_kappa,
                    "bc_tau": bc_tau,
                    "candidates": candidates,
                    "window_size": window_size,
                }
                st.session_state["selected_pocket_idx"] = 0
        else:
            with st.spinner("Analyzing intrinsic curve invariants & Sarle's bimodality moments in Rust core..."):
                scores, bc_tau, bc_kappa = tf.compute_bimodality_profile(traj_data, window_size=window_size)
                candidates = tf.scan_cryptic_pockets(traj_data, window_size=window_size, bc_threshold=bc_threshold)
                st.session_state["scan_results"] = {
                    "is_idp": False,
                    "is_rna": False,
                    "scores": scores,
                    "bc_tau": bc_tau,

                    "bc_kappa": bc_kappa,
                    "candidates": candidates,
                    "window_size": window_size,
                }
                st.session_state["selected_pocket_idx"] = 0

    scan_res = st.session_state["scan_results"]

    if scan_res is not None:
        if scan_res.get("is_idp", False):
            s_topo = scan_res["s_topo"]
            z_scores = scan_res["z_scores"]
            motifs = scan_res["motifs"]
            z_th = scan_res["z_threshold"]

            with col_stats:
                st.success(
                    f"**IDP Analysis Complete**: Evaluated {traj_data.shape[0]:,} disordered conformations across {n_residues} residues. "
                    f"Discovered **{len(motifs)}** transient pre-nucleation motif(s) (Peak Z = {float(np.max(z_scores)):.2f})."
                )

            col_plot, col_view = st.columns([1.15, 1.0])

            with col_plot:
                st.subheader("📈 Spectral Topological Density Profile")
                x_res_pdb = res_ids

                if go is not None:
                    fig = go.Figure()
                    # Add shaded region for NAC domain (61..95)
                    fig.add_vrect(
                        x0=61, x1=95,
                        fillcolor="rgba(255, 87, 34, 0.15)",
                        layer="below", line_width=1,
                        line_color="rgba(255, 87, 34, 0.4)",
                        annotation_text="NAC Domain (61–95)",
                        annotation_position="top left",
                    )
                    # Add trace for S_topo
                    fig.add_trace(go.Scatter(
                        x=x_res_pdb, y=s_topo, mode="lines",
                        name="Spectral Density S_topo(i)",
                        line=dict(color="#00E676", width=2.5),
                        yaxis="y1",
                    ))
                    # Add trace for Z-scores
                    fig.add_trace(go.Scatter(
                        x=x_res_pdb, y=z_scores, mode="lines",
                        name="Topological Z-Score Z(i)",
                        line=dict(color="#FF9800", width=2.0, dash="dash"),
                        yaxis="y2",
                    ))
                    # Z threshold line on secondary y-axis
                    fig.add_hline(
                        y=z_th, line_dash="dot", line_color="#E91E63",
                        annotation_text=f"Z Cutoff ({z_th:.1f})",
                        annotation_position="bottom right",
                        yref="y2",
                    )
                    fig.update_layout(
                        xaxis_title="PDB Residue Sequence Number",
                        yaxis=dict(
                            title="Spectral Density S_topo",
                            title_font=dict(color="#00E676"),
                            tickfont=dict(color="#00E676"),
                        ),
                        yaxis2=dict(
                            title="Topological Z-Score",
                            title_font=dict(color="#FF9800"),
                            tickfont=dict(color="#FF9800"),
                            overlaying="y",
                            side="right",
                        ),
                        template="plotly_dark",
                        margin=dict(l=40, r=40, t=30, b=40),
                        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
                    )
                    st.plotly_chart(fig)
                else:
                    st.line_chart(s_topo)

                st.subheader("🏆 Detected Transient Pre-Nucleation Motifs")
                if motifs:
                    motif_data = []
                    for r, m in enumerate(motifs, 1):
                        s_res, e_res, p_res, mean_c, peak_c, z_sc = m
                        pdb_start = res_ids[min(s_res, n_residues - 1)]
                        pdb_end = res_ids[min(e_res, n_residues - 1)]
                        pdb_peak = res_ids[min(p_res, n_residues - 1)]
                        motif_data.append({
                            "Rank": f"#{r}",
                            "Residue Range": f"Res {pdb_start} .. {pdb_end}",
                            "Peak Residue": f"Res {pdb_peak}",
                            "Mean S_topo": f"{mean_c:.4f}",
                            "Peak S_topo": f"{peak_c:.4f}",
                            "Peak Z-Score": f"{z_sc:.2f}",
                        })
                    st.dataframe(motif_data)
                else:
                    st.info("No transient motifs detected above current Z threshold.")

            with col_view:
                st.subheader("🔬 3D Molecular Structure & Nucleation Focus")
                if motifs:
                    top_m = motifs[0]
                    def_start = res_ids[min(top_m[0], n_residues - 1)]
                    def_end = res_ids[min(top_m[1], n_residues - 1)]
                else:
                    def_start, def_end = (66, 78) if n_residues >= 78 else (res_ids[0], res_ids[-1])

                selected_pdb_range = (def_start, def_end)

                if motifs:
                    motif_options = []
                    for i, m in enumerate(motifs):
                        s_res, e_res, p_res, _, _, z_sc = m
                        p_start = res_ids[min(s_res, n_residues - 1)]
                        p_end = res_ids[min(e_res, n_residues - 1)]
                        motif_options.append(f"Motif #{i+1}: Res {p_start}..{p_end} (Z = {z_sc:.2f})")

                    sel_idx = st.selectbox(
                        "Select Detected Motif to Highlight:",
                        range(len(motifs)),
                        format_func=lambda i: motif_options[i],
                        key="idp_motif_selectbox",
                    )
                    sel_m = motifs[sel_idx]
                    selected_pdb_range = (
                        res_ids[min(sel_m[0], n_residues - 1)],
                        res_ids[min(sel_m[1], n_residues - 1)],
                    )

                st.markdown(
                    f"**Active Highlight**: Transiently nucleated segment PDB Residues `{selected_pdb_range[0]}..{selected_pdb_range[1]}` "
                    f"*(Vivid Orange/Red Cartoon & Sticks)* with disordered tails *(Disordered Wireframe/Trace)*."
                )

                html_viewer = render_3dmol_viewer(
                    ref_pdb_str,
                    highlight_range=selected_pdb_range,
                    res_count=len(res_ids),
                    height=520,
                    is_idp=True,
                )
                st.components.v1.html(html_viewer, height=565)

        elif scan_res.get("is_rna", False):
            scores = scan_res["scores"]
            bc_theta = scan_res["bc_theta"]
            bc_kappa = scan_res["bc_kappa"]
            bc_tau = scan_res["bc_tau"]
            candidates = scan_res["candidates"]
            w_size = scan_res["window_size"]

            with col_stats:
                st.success(
                    f"**RNA Analysis Complete**: Evaluated {len(scores)} sliding windows across {traj_data.shape[0]:,} frames. "
                    f"Discovered **{len(candidates)}** candidate conformational switching hinge(s) (Peak BC = {float(np.max(scores)):.4f})."
                )

            col_plot, col_view = st.columns([1.15, 1.0])

            with col_plot:
                st.subheader("📈 Sequence-Wide Ribonucleic Bimodality Profile")
                mid_offset = 1
                x_res_pdb = [res_ids[min(i + mid_offset, n_residues - 1)] for i in range(len(scores))]

                if go is not None:
                    fig = go.Figure()
                    # Add shaded regions for P1 stems: 13..21 and 74..82 if present
                    if 13 in res_ids and 21 in res_ids:
                        fig.add_vrect(
                            x0=13, x1=21,
                            fillcolor="rgba(255, 87, 34, 0.12)",
                            layer="below", line_width=1,
                            line_color="rgba(255, 87, 34, 0.3)",
                            annotation_text="P1 5' Stem (13–21)",
                            annotation_position="top left",
                        )
                    if 74 in res_ids and 82 in res_ids:
                        fig.add_vrect(
                            x0=74, x1=82,
                            fillcolor="rgba(255, 87, 34, 0.22)",
                            layer="below", line_width=1,
                            line_color="rgba(255, 87, 34, 0.5)",
                            annotation_text="P1 Switching Stem (74–82)",
                            annotation_position="top right",
                        )
                    fig.add_trace(go.Scatter(
                        x=x_res_pdb, y=scores, mode="lines",
                        name="Composite max(BC_θ, BC_κ, BC_τ)",
                        line=dict(color="#00E676", width=3)
                    ))
                    fig.add_trace(go.Scatter(
                        x=x_res_pdb, y=bc_theta, mode="lines",
                        name="Base Ribbon Dihedral (BC_θ)",
                        line=dict(color="#FF9800", width=2.0, dash="dash")
                    ))
                    fig.add_trace(go.Scatter(
                        x=x_res_pdb, y=bc_kappa, mode="lines",
                        name="Phosphorus Curvature (BC_κ)",
                        line=dict(color="#AB47BC", width=1.5, dash="dot")
                    ))
                    fig.add_trace(go.Scatter(
                        x=x_res_pdb, y=bc_tau, mode="lines",
                        name="Phosphorus Torsion (BC_τ)",
                        line=dict(color="#29B6F6", width=1.5, dash="dashdot")
                    ))
                    fig.add_hline(
                        y=0.555, line_dash="dash", line_color="#e74c3c",
                        annotation_text="Unimodal Threshold (BC=0.555)", annotation_position="bottom right"
                    )
                    fig.update_layout(
                        xaxis_title="PDB Nucleotide Sequence Number",
                        yaxis_title="Sarle's Bimodality Coefficient (BC)",
                        yaxis_range=[-0.02, 1.05],
                        template="plotly_dark",
                        margin=dict(l=40, r=20, t=30, b=40),
                        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
                    )
                    st.plotly_chart(fig)
                else:
                    st.line_chart(scores)

                st.subheader("🏆 Discovered RNA Switching Hinges")
                if candidates:
                    cand_data = []
                    for r, c in enumerate(candidates, 1):
                        s_res = c["start"]
                        e_res = c["end"]
                        sc = c["score"]
                        pdb_start = res_ids[min(s_res, n_residues - 1)]
                        pdb_end = res_ids[min(e_res, n_residues - 1)]
                        cand_data.append({
                            "Rank": f"#{r}",
                            "PDB Nucleotides": f"Res {pdb_start} .. {pdb_end}",
                            "Length": e_res - s_res + 1,
                            "Peak BC": f"{sc:.4f}",
                            "BC_θ (Base)": f"{c['bc_theta']:.4f}",
                            "BC_κ (Curv)": f"{c['bc_kappa']:.4f}",
                            "BC_τ (Tors)": f"{c['bc_tau']:.4f}",
                        })
                    st.dataframe(cand_data)
                else:
                    st.info("No switching hinges detected above current threshold.")

            with col_view:
                st.subheader("🔬 3D Molecular Structure & Switching Hinge")
                selected_pdb_range = (74, 82) if (74 in res_ids and 82 in res_ids) else (res_ids[0], res_ids[-1])

                if candidates:
                    cand_options = []
                    for i, c in enumerate(candidates):
                        s_res = c["start"]
                        e_res = c["end"]
                        sc = c["score"]
                        pdb_start = res_ids[min(s_res, n_residues - 1)]
                        pdb_end = res_ids[min(e_res, n_residues - 1)]
                        cand_options.append(f"Rank #{i+1}: Res {pdb_start}..{pdb_end} (BC = {sc:.4f})")

                    selected_idx = st.selectbox(
                        "Select Switching Hinge to Highlight:",
                        range(len(candidates)),
                        format_func=lambda i: cand_options[i],
                        key="rna_hinge_selectbox",
                    )
                    st.session_state["selected_pocket_idx"] = selected_idx
                    sel_cand = candidates[selected_idx]
                    selected_pdb_range = (
                        res_ids[min(sel_cand["start"], n_residues - 1)],
                        res_ids[min(sel_cand["end"], n_residues - 1)],
                    )

                st.markdown(
                    f"**Active Highlight**: RNA Switching Hinge PDB Nucleotides `{selected_pdb_range[0]}..{selected_pdb_range[1]}` "
                    f"*(Vivid Orange/Red Cartoon Ribbon & Sticks)* on scaffold *(Semi-Transparent Light Grey)*."
                )

                html_viewer = render_3dmol_viewer(
                    ref_pdb_str,
                    highlight_range=selected_pdb_range,
                    res_count=len(res_ids),
                    height=520,
                    is_rna=True,
                )
                st.components.v1.html(html_viewer, height=565)

        else:
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
                    f"*(Orange/Red Cartoon & Sticks)* on scaffold *(Light Grey Cartoon)*."
                )

                html_viewer = render_3dmol_viewer(
                    ref_pdb_str,
                    highlight_range=selected_pdb_range,
                    res_count=len(res_ids),
                    height=520,
                )
                st.components.v1.html(html_viewer, height=565)

    else:
        # Initial Structure Preview before scan execution
        col_preview_info, col_preview_view = st.columns([1, 1.25])
        with col_preview_info:
            if is_idp:
                st.info(
                    "💡 Click **'⚡ Run Spectral IDP Scan'** above to compute the ensemble Spectral Topological Density and identify pre-nucleation motifs."
                )
                st.markdown(f"**Loaded Scaffold**: {curr_ds['desc']}")
                st.markdown(
                    f"- **Trajectory Frames**: `{traj_data.shape[0]:,}`\n"
                    f"- **Residue Count**: `{n_residues}` (PDB `{res_ids[0]}..{res_ids[-1]}`)\n"
                    f"- **Ensemble Characteristic**: Intrinsically Disordered Protein (IDP) ensemble with large Cartesian variance (RMSD > 18 Å)\n"
                    f"- **Biophysical Target**: Transient nucleation within hydrophobic NACore (residues 61..95)"
                )
            elif is_rna:
                st.info(
                    "💡 Click **'⚡ Run Ribonucleic Hinge Scan'** above to detect RNA conformational switching hinges and bistable regulatory elements."
                )
                st.markdown(f"**Loaded Scaffold**: {curr_ds['desc']}")
                st.markdown(
                    f"- **Trajectory Frames**: `{traj_data.shape[0]:,}`\n"
                    f"- **Nucleotide Count**: `{n_residues}` (PDB `{res_ids[0]}..{res_ids[-1]}`)\n"
                    f"- **Ribonucleic Invariants**: Phosphorus backbone $\\kappa_P, \\tau_P, \\text{{Wr}}_P$, and Glycosidic Ribbon $\\theta_{{\\text{{base}}}}$\n"
                    f"- **Biophysical Target**: Autonomous detection of P1 switching terminator hinge (residues 74..82)"
                )
            else:
                st.info(
                    "💡 Click **'🚀 Run Autonomous Pocket Scan'** above to detect cryptic pockets and mobile functional loops."
                )
                st.markdown(f"**Loaded Scaffold**: {curr_ds['desc']}")
                st.markdown(
                    f"- **Trajectory Frames**: `{traj_data.shape[0]:,}`\n"
                    f"- **Residue Count**: `{n_residues}` (PDB `{res_ids[0]}..{res_ids[-1]}`)\n"
                    f"- **Differential Invariants**: Curvature $\\kappa$, Torsion $\\tau$, Side-Chain Ribbon $\\theta_\\beta$"
                )
        with col_preview_view:
            st.subheader("🔬 3D Molecular Structure (Scaffold Preview)")
            if is_idp and n_residues >= 78:
                preview_hl = (66, 78)
            elif is_rna and 74 in res_ids and 82 in res_ids:
                preview_hl = (74, 82)
            else:
                preview_hl = None
            html_viewer = render_3dmol_viewer(
                ref_pdb_str,
                highlight_range=preview_hl,
                res_count=len(res_ids),
                height=520,
                is_rna=is_rna,
                is_idp=is_idp,
            )
            st.components.v1.html(html_viewer, height=565)


if __name__ == "__main__":
    if st is None:
        print("TopoFold Interactive 3D Dashboard")
        print("Please install Streamlit to run GUI: pip install streamlit plotly")
    else:
        run_dashboard()
