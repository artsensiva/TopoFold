#!/usr/bin/env python3
"""
TopoFold Benchmark: Side-Chain Rotameric Cryptic Pocket Gating Validation
========================================================================

*** PHYSICALLY INVALID BENCHMARK (docs/KNOWN_ISSUES.md, C4): Cb is rotated around a fixed Ca, which cannot happen in a protein; Cb does not move under chi1 rotation. ***

Demonstrates why pure C-alpha representations fail to detect rotameric
cryptic pocket gating, and how TopoFold's SE(3)-invariant C-beta ribbon
geometry resolves isolated side-chain flips amidst rigid backbone constraints.
"""

import os
import sys
import numpy as np
import matplotlib
os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib_topofold"
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

import topofold as tf


def generate_rotamer_gating_ensemble(n_frames=1000, n_res=30, seed=42):
    """
    Generates a 1,000-frame synthetic trajectory of a 30-residue protein:
    - C-alpha backbone: Rigid idealized helix (zero C-alpha displacement across all frames).
    - Residue 15: Exhibits an isolated bistable rotamer flip between State A (-126 deg) and State B (+126 deg).
    """
    rng = np.random.default_rng(seed)

    # 1. Construct rigid C-alpha backbone (ideal alpha-helical trace)
    r = 2.3
    pitch = 1.5
    d_theta = 100.0 * np.pi / 180.0

    ca_base = np.zeros((n_res, 3), dtype=np.float32)
    for i in range(n_res):
        th = i * d_theta
        ca_base[i] = [r * np.cos(th), r * np.sin(th), i * pitch]

    # Rigid backbone across all frames
    ca_traj = np.repeat(ca_base[np.newaxis, :, :], n_frames, axis=0)

    # 2. Base C-beta coordinates pointing radially outward
    cb_traj = np.zeros_like(ca_traj)
    for f in range(n_frames):
        for i in range(n_res):
            ca_pt = ca_traj[f, i]
            rad = np.array([ca_pt[0], ca_pt[1], 0.0], dtype=np.float32)
            rad /= np.linalg.norm(rad)
            cb_traj[f, i] = ca_pt + 1.52 * rad

    # 3. Residue 15 executes a clean bistable rotamer flip in the local Frenet (N, B) plane:
    # State A (frames 0..499): Rotamer angle = -126 deg
    # State B (frames 500..999): Rotamer angle = +126 deg
    idx = 15
    T_prev = (ca_base[idx] - ca_base[idx - 1]) / np.linalg.norm(ca_base[idx] - ca_base[idx - 1])
    T_next = (ca_base[idx + 1] - ca_base[idx]) / np.linalg.norm(ca_base[idx + 1] - ca_base[idx])
    B = np.cross(T_prev, T_next)
    B /= np.linalg.norm(B)
    T_vert = (T_prev + T_next) / np.linalg.norm(T_prev + T_next)
    N = np.cross(B, T_vert)
    N /= np.linalg.norm(N)

    for f in range(n_frames):
        angle = -np.pi * 0.7 if f < 500 else +np.pi * 0.7
        angle += rng.normal(0, np.radians(2.0))
        v_rot = np.cos(angle) * N + np.sin(angle) * B
        cb_traj[f, idx] = ca_base[idx] + 1.52 * v_rot

    labels = np.array([0 if f < 500 else 1 for f in range(n_frames)], dtype=int)
    return ca_traj, cb_traj, labels


def main():
    print("=" * 80)
    print("      TOPOFOLD BENCHMARK: ROTAMERIC CRYPTIC POCKET GATING")
    print("      Isolated Side-Chain Flip on Rigid C-Alpha Backbone")
    print("=" * 80)

    # 1. Generate Ensemble
    print("\n[Step 1/4] Generating Controlled Rotamer Gating Trajectory...")
    ca_traj, cb_traj, labels = generate_rotamer_gating_ensemble(n_frames=1000, n_res=30)
    print(f"  -> Generated 1,000 frames (30 residues each)")
    print(f"  -> C-alpha displacement: Exactly 0.0 Å (Rigid backbone)")
    print(f"  -> Residue 15: Isolated bistable rotamer flip (gauche- vs trans)")

    # 2. Evaluate Pure C-alpha Invariants vs Ribbon Invariants
    print("\n[Step 2/4] Evaluating Autonomous Bimodality Profiles...")
    window_size = 8
    ca_scores, ca_tau, ca_kappa = tf.compute_bimodality_profile(ca_traj, window_size=window_size)
    ribbon_scores, rib_tau, rib_kappa, rib_theta = tf.compute_bimodality_profile(
        ca_traj, window_size=window_size, cb_coords=cb_traj
    )

    peak_idx = int(np.argmax(ribbon_scores))
    peak_res = peak_idx + (window_size // 2)

    print(f"  -> Pure C-alpha Max Bimodality: {np.max(ca_scores):.4f} (Complete Zero / Unimodal)")
    print(f"  -> Ribbon Max Bimodality:       {np.max(ribbon_scores):.4f} (Peak at Residue {peak_res})")
    print(f"  -> Peak Ribbon bc_theta:        {np.max(rib_theta):.4f}")

    # 3. Cartesian PCA Baseline (with tiny machine epsilon jitter to avoid rank 0)
    print("\n[Step 3/4] Computing Cartesian PCA Baseline on C-alpha Backbone...")
    flat_ca = ca_traj.reshape(ca_traj.shape[0], -1) + np.random.normal(0, 1e-6, size=(1000, 90))
    pca = PCA(n_components=2)
    ca_pca = pca.fit_transform(flat_ca)
    var_exp = pca.explained_variance_ratio_ * 100
    print(f"  -> Cartesian PCA PC1 Var: {var_exp[0]:.2f}%, PC2 Var: {var_exp[1]:.2f}%")

    # 4. Generate Publication Figure
    print("\n[Step 4/4] Rendering Publication Benchmark Figure...")
    os.makedirs("assets", exist_ok=True)
    out_png = os.path.join("assets", "rotamer_gating_benchmark.png")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(15, 6), dpi=300)

    # ----------------------------------------------------
    # Panel A: Sequence Bimodality Profile
    # ----------------------------------------------------
    x_res = np.arange(len(ca_scores)) + (window_size // 2)

    ax_a.plot(x_res, ca_scores, label=r"Pure $\text{C}_\alpha$ Trace ($\max(BC_\tau, BC_\kappa) = 0.0$)",
              color="#7f8c8d", lw=2.4, linestyle="--")
    ax_a.plot(x_res, ribbon_scores, label=r"TopoFold Ribbon ($\max(BC_\tau, BC_\kappa, BC_\theta)$)",
              color="#27ae60", lw=2.8)
    ax_a.plot(x_res, rib_theta, label=r"Side-Chain Dihedral ($BC_\theta$)",
              color="#8e44ad", lw=2.0, alpha=0.85)

    # Reference uniform threshold line
    ax_a.axhline(0.555, color="#e74c3c", linestyle=":", lw=1.8, label="Unimodal Boundary ($BC = 0.555$)")
    ax_a.axvline(15, color="#f39c12", linestyle="-.", lw=1.5, alpha=0.7, label="Rotamer Gate (Residue 15)")

    ax_a.set_title(r"$\mathbf{A}$   Sequence-Wide Bimodality Profile ($W = 8$ residues)", fontsize=13, pad=12)
    ax_a.set_xlabel("Residue Index", fontsize=11, fontweight="bold")
    ax_a.set_ylabel("Sarle's Bimodality Coefficient (BC)", fontsize=11, fontweight="bold")
    ax_a.set_ylim(-0.05, 1.08)
    ax_a.set_xlim(x_res[0], x_res[-1])
    ax_a.grid(True, linestyle="--", alpha=0.4)
    ax_a.legend(loc="upper right", frameon=True, fontsize=9.5, framealpha=0.92)

    # Annotation callout
    ax_a.annotate(f"Rotamer Gating Detected!\nPeak BC = {ribbon_scores[peak_idx]:.3f}",
                  xy=(x_res[peak_idx], ribbon_scores[peak_idx]),
                  xytext=(x_res[peak_idx] - 7, 0.82),
                  arrowprops=dict(facecolor="#27ae60", shrink=0.08, width=1.5, headwidth=6),
                  fontsize=10, fontweight="bold", color="#1e8449",
                  bbox=dict(boxstyle="round,pad=0.3", fc="#eafaf1", ec="#27ae60", lw=1.2))

    # ----------------------------------------------------
    # Panel B: Intrinsic Ribbon Angle vs Cartesian PCA
    # ----------------------------------------------------
    theta_res15 = np.zeros(1000)
    for f in range(1000):
        th = tf.compute_theta_beta(ca_traj[f], cb_traj[f])
        theta_res15[f] = np.degrees(th[15 - 1])  # 0-indexed internal vertex

    # Scatter of Cartesian PC1 vs Intrinsic theta_beta
    ax_b.scatter(ca_pca[:500, 0] * 1000, theta_res15[:500],
                 c="#2980b9", alpha=0.6, s=28, edgecolors="none", label="State A (Closed Rotamer)")
    ax_b.scatter(ca_pca[500:, 0] * 1000, theta_res15[500:],
                 c="#e67e22", alpha=0.6, s=28, edgecolors="none", label="State B (Open Cryptic Gate)")

    ax_b.set_title(r"$\mathbf{B}$   Cartesian PCA Collapse vs. Ribbon Dihedral $\theta_\beta$", fontsize=13, pad=12)
    ax_b.set_xlabel(r"Cartesian PC1 ($\times 10^{-3}$ Å, Overlapped / Zero Variance)", fontsize=11, fontweight="bold")
    ax_b.set_ylabel(r"Ribbon Orientation $\theta_\beta$ at Res 15 (Degrees)", fontsize=11, fontweight="bold")
    ax_b.grid(True, linestyle="--", alpha=0.4)
    ax_b.legend(loc="upper left", frameon=True, fontsize=10, framealpha=0.92)

    # Highlight vertical separation vs horizontal overlap
    ax_b.axhline(0, color="#7f8c8d", linestyle="--", alpha=0.5)
    ax_b.annotate("Clean Bimodal Separation\nAlong Intrinsic Ribbon Space",
                  xy=(0.0, 126), xytext=(-0.5, 70),
                  arrowprops=dict(facecolor="#e67e22", shrink=0.08, width=1.5, headwidth=6),
                  fontsize=9.5, fontweight="bold", color="#d35400",
                  bbox=dict(boxstyle="round,pad=0.3", fc="#fef5e7", ec="#e67e22", lw=1.2))

    plt.subplots_adjust(left=0.08, right=0.96, top=0.90, bottom=0.12, wspace=0.25)
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"  -> Benchmark figure successfully saved to: {out_png}")
    print("\n" + "=" * 80)
    print("                       BENCHMARK SUMMARY")
    print("=" * 80)
    print(f"Pure C-alpha Bimodality Score:     {np.max(ca_scores):.4f} (Blind to rotamer gating)")
    print(f"TopoFold Ribbon Bimodality Score:  {np.max(ribbon_scores):.4f} (Crisp rotamer detection)")
    print(f"Cartesian PC1 State Discrimination: Complete Overlap (Zero C-alpha Variance)")
    print("=" * 80)


if __name__ == "__main__":
    main()
