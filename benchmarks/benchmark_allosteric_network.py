#!/usr/bin/env python3
"""
TopoFold Benchmark: Intrinsic Allosteric Communication Networks (Abl1 Kinase)
=============================================================================

Computes the sequence-wide, SE(3)-invariant Allosteric Communication Network
matrix on the Human c-Abl1 Kinase Domain (PDB 2GQG / 1IEP, 274 residues, 1,500 frames)
via Mutual Information of discrete curve invariants (curvature, torsion, ribbon dihedrals).

Key Advantages over Traditional Cross-Correlation:
  - Coordinate Invariance: Independent of global Cartesian rotation or Kabsch alignment.
  - Non-Linear Coupling: Mutual Information via Gaussian Copula captures non-linear
    conformational correlations missed by linear Pearson/Cartesian covariance.
  - High Performance: Fully parallelized Rust Rayon core computes 37,000+ residue pairs
    across thousands of frames in milliseconds.

Usage:
    python benchmarks/benchmark_allosteric_network.py
"""

import os
import sys
import time
import argparse
import numpy as np

os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib_topofold"
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import topofold as tf

BENCHMARK_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "assets")


def run_benchmark():
    parser = argparse.ArgumentParser(description="TopoFold Allosteric Network Benchmark: Abl1 Kinase")
    parser.add_argument("--trajectory", default=os.path.join(BENCHMARK_DATA_DIR, "abl_dfg_trajectory.dcd"))
    parser.add_argument("--reference-pdb", default=os.path.join(BENCHMARK_DATA_DIR, "abl_reference.pdb"))
    parser.add_argument("--output-figure", default=os.path.join(ASSETS_DIR, "abl_allosteric_network_matrix.png"))
    parser.add_argument("--regularizer-eps", type=float, default=1e-6)
    args = parser.parse_args()

    os.makedirs(ASSETS_DIR, exist_ok=True)

    if not os.path.exists(args.trajectory):
        print(f"Error: Trajectory {args.trajectory} not found. Run benchmark_abl_kinase_dfg.py first.")
        sys.exit(1)

    print("=" * 78)
    print(" TOPOFOLD BENCHMARK: INTRINSIC ALLOSTERIC NETWORKS (ABL1 KINASE)")
    print("=" * 78)

    # Ingest Trajectory
    print(f"Ingesting binary DCD trajectory: {args.trajectory} ...")
    traj = tf.read_dcd(args.trajectory)
    n_frames, n_residues, _ = traj.shape
    print(f"Trajectory Shape: {n_frames} frames across {n_residues} residues")

    # Offset to authentic PDB residue numbers (Abl1 kinase catalytic domain starts at PDB Met225)
    pdb_offset = 225

    # Compute Allosteric Network in Rust Core
    print(f"\nComputing Intrinsic Allosteric Communication Matrix (Rayon parallel core)...")
    t0 = time.perf_counter()
    net = tf.compute_allosteric_network(traj, regularizer_eps=args.regularizer_eps)
    elapsed = time.perf_counter() - t0
    total_pairs = (n_residues * (n_residues - 1)) // 2

    print(f"✓ Allosteric Network Computed in {elapsed*1000:.2f} ms")
    print(f"  -> Total Residue Pairs: {total_pairs:,}")
    print(f"  -> Throughput:          {elapsed / total_pairs * 1e6:.3f} µs/pair ({elapsed / n_frames * 1e6:.1f} µs/frame)")

    # Compute Allosteric Centrality (Sum of non-local correlations |j - i| >= 4)
    centrality = np.zeros(n_residues)
    for i in range(n_residues):
        for j in range(n_residues):
            if abs(i - j) >= 4:
                centrality[i] += net[i, j]

    # Identify Top Communicating Hubs
    hub_indices = np.argsort(centrality)[::-1]
    print("\n--- Top Allosteric Communication Hubs in Abl1 Kinase ---")
    print(f"{'Rank':^6}|{'0-Idx':^8}|{'PDB Residue':^14}|{'Centrality (Σ r_MI)':^22}")
    print("-" * 52)
    for r in range(10):
        idx = hub_indices[r]
        pdb_res = idx + pdb_offset
        print(f"{r+1:^6}|{idx:^8}|{f'Res {pdb_res}':^14}|{centrality[idx]:^22.4f}")

    # Top Non-Local Correlated Pairs (|j - i| >= 8)
    non_local_pairs = []
    for i in range(n_residues):
        for j in range(i + 8, n_residues):
            non_local_pairs.append((i, j, net[i, j]))
    non_local_pairs.sort(key=lambda x: x[2], reverse=True)

    print("\n--- Top Non-Local Allosteric Coupling Pathways (|j - i| >= 8) ---")
    print(f"{'Rank':^6}|{'Pair (PDB)':^24}|{'Coupling r_MI':^16}")
    print("-" * 50)
    for r in range(8):
        i, j, val = non_local_pairs[r]
        print(f"{r+1:^6}|{f'Res {i+pdb_offset} ↔ Res {j+pdb_offset}':^24}|{val:^16.4f}")

    # Key Functional Motifs in Abl1
    # P-loop: residues 248..256 (0-idx 23..31)
    # alphaC-helix: residues 280..295 (0-idx 55..70)
    # Catalytic Loop: residues 314..325 (0-idx 89..100)
    # DFG motif: residues 380..385 (0-idx 155..160)
    # Activation loop: residues 386..405 (0-idx 161..180)
    motifs = [
        ("P-loop (248..256)", 23, 31, "#3498db"),
        ("αC-helix (280..295)", 55, 70, "#2ecc71"),
        ("Cat-loop (314..325)", 89, 100, "#9b59b6"),
        ("DFG / A-loop (380..400)", 155, 175, "#e67e22"),
    ]

    # Generate 300 DPI Publication Figure
    fig = plt.figure(figsize=(19, 5.8), dpi=300)
    fig.patch.set_facecolor("#0e1117")

    # -------------------------------------------------------------
    # Panel A: Sequence-Wide Allosteric Correlation Heatmap
    # -------------------------------------------------------------
    ax_a = fig.add_subplot(1, 3, 1)
    ax_a.set_facecolor("#161b22")

    # Mask diagonal for visual contrast in non-local pathways
    net_display = net.copy()
    np.fill_diagonal(net_display, 0.0)

    pdb_ticks = np.arange(0, n_residues, 40)
    pdb_labels = [str(t + pdb_offset) for t in pdb_ticks]

    im = ax_a.imshow(
        net_display,
        cmap="magma",
        origin="lower",
        norm=mcolors.PowerNorm(gamma=0.5, vmin=0.0, vmax=1.0),
    )
    ax_a.set_title(
        "A. Sequence-Wide Allosteric Network Matrix\nGeneralized Correlation r_MI (SE(3)-Invariant)",
        fontsize=11, fontweight="bold", color="#f0f6fc", pad=12
    )
    ax_a.set_xlabel("PDB Residue Sequence Number", fontsize=9, color="#c9d1d9")
    ax_a.set_ylabel("PDB Residue Sequence Number", fontsize=9, color="#c9d1d9")
    ax_a.set_xticks(pdb_ticks)
    ax_a.set_xticklabels(pdb_labels, fontsize=8, color="#8b949e")
    ax_a.set_yticks(pdb_ticks)
    ax_a.set_yticklabels(pdb_labels, fontsize=8, color="#8b949e")

    # Annotate functional boxes
    for name, s, e, color in motifs:
        rect = plt.Rectangle((s, s), e - s, e - s, fill=False, edgecolor=color, linewidth=1.5, linestyle="--")
        ax_a.add_patch(rect)

    cbar = fig.colorbar(im, ax=ax_a, fraction=0.046, pad=0.04)
    cbar.set_label("Coupling Strength r_MI", color="#c9d1d9", fontsize=9)
    cbar.ax.tick_params(colors="#8b949e", labelsize=8)

    # -------------------------------------------------------------
    # Panel B: Residue-Resolved Communication Centrality Profile
    # -------------------------------------------------------------
    ax_b = fig.add_subplot(1, 3, 2)
    ax_b.set_facecolor("#161b22")

    res_pdb_axis = np.arange(n_residues) + pdb_offset
    ax_b.plot(res_pdb_axis, centrality, color="#2ecc71", linewidth=2.0, label="Centrality Σ r_MI (|j-i|≥4)")
    ax_b.fill_between(res_pdb_axis, 0, centrality, color="#2ecc71", alpha=0.15)

    # Highlight motifs
    for name, s, e, color in motifs:
        ax_b.axvspan(s + pdb_offset, e + pdb_offset, color=color, alpha=0.18, label=name)

    ax_b.set_title(
        "B. Residue Allosteric Communication Centrality\nIdentification of Kinase Allosteric Drivers",
        fontsize=11, fontweight="bold", color="#f0f6fc", pad=12
    )
    ax_b.set_xlabel("PDB Residue Sequence Number", fontsize=9, color="#c9d1d9")
    ax_b.set_ylabel("Communication Centrality (Σ r_MI)", fontsize=9, color="#c9d1d9")
    ax_b.tick_params(colors="#8b949e", labelsize=8)
    ax_b.grid(True, linestyle=":", alpha=0.25, color="#8b949e")
    ax_b.legend(loc="upper right", fontsize=7.5, facecolor="#161b22", edgecolor="#30363d", labelcolor="#f0f6fc")

    # -------------------------------------------------------------
    # Panel C: Catalytic Triad Sub-Network Coupling
    # -------------------------------------------------------------
    ax_c = fig.add_subplot(1, 3, 3)
    ax_c.set_facecolor("#161b22")

    # Zoom in on P-loop (23..35) vs DFG / Activation Loop (150..175)
    sub_matrix = net[150:175, 20:45]
    im_c = ax_c.imshow(
        sub_matrix,
        cmap="viridis",
        origin="lower",
        aspect="auto",
        vmin=0.0, vmax=0.4,
    )
    ax_c.set_title(
        "C. Allosteric Cross-Talk: DFG Switch ↔ P-Loop\nLong-Range Mechanical Coupling Sub-Matrix",
        fontsize=11, fontweight="bold", color="#f0f6fc", pad=12
    )
    ax_c.set_xlabel("P-Loop / Hinge Residues (PDB 245..270)", fontsize=9, color="#c9d1d9")
    ax_c.set_ylabel("DFG / A-Loop Residues (PDB 375..400)", fontsize=9, color="#c9d1d9")

    p_ticks = np.arange(0, 25, 5)
    p_labels = [str(t + 20 + pdb_offset) for t in p_ticks]
    d_ticks = np.arange(0, 25, 5)
    d_labels = [str(t + 150 + pdb_offset) for t in d_ticks]

    ax_c.set_xticks(p_ticks)
    ax_c.set_xticklabels(p_labels, fontsize=8, color="#8b949e")
    ax_c.set_yticks(d_ticks)
    ax_c.set_yticklabels(d_labels, fontsize=8, color="#8b949e")

    cbar_c = fig.colorbar(im_c, ax=ax_c, fraction=0.046, pad=0.04)
    cbar_c.set_label("Coupling r_MI", color="#c9d1d9", fontsize=9)
    cbar_c.ax.tick_params(colors="#8b949e", labelsize=8)

    plt.tight_layout()
    plt.savefig(args.output_figure, dpi=300, facecolor=fig.get_facecolor(), edgecolor="none")
    plt.close()

    print(f"\n✓ Generated 300 DPI publication figure: {args.output_figure}")
    print("=" * 78)


if __name__ == "__main__":
    run_benchmark()
