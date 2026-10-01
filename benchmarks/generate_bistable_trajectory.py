#!/usr/bin/env python3
"""
TopoFold Benchmark Trajectory Generator (Phase 4).

Generates a synthetic 2,000-frame molecular dynamics trajectory of a 60-residue
C-alpha backbone exhibiting a clean bistable allosteric transition in a core functional
loop (residues 25..35), surrounded by large-amplitude stochastic Brownian motions in the
terminal arms (residues 0..24 and 36..59).

Scientific Model:
-----------------
1. Backbone geometry is generated under strict physical peptide chain constraints:
   - C_alpha - C_alpha bond lengths: 3.80 +/- 0.03 Å.
   - Non-zero discrete turning angles (kappa in [0.7, 2.2] rad).
2. Bistable transition in residues 25..35:
   - Frames 0..999 (State A): Compact helical loop turn (closed cryptic pocket).
     kappa_loop ~ 1.55 +/- 0.03 rad, tau_loop ~ 0.85 +/- 0.04 rad.
   - Frames 1000..1999 (State B): Extended hairpin loop (open cryptic pocket).
     kappa_loop ~ 0.95 +/- 0.03 rad, tau_loop ~ 2.75 +/- 0.04 rad.
3. Terminal arms (0..24 and 36..59):
   - Undergo continuous thermal stochastic fluctuations.
   - Creates massive global Cartesian variance (> 80% of total variance),
     which severely confounds Cartesian PCA despite optimal superposition.

Outputs:
--------
- benchmarks/data/bistable_ensemble.dcd
- benchmarks/data/state_a.pdb
- benchmarks/data/state_b.pdb
"""

import os
import struct
import numpy as np


def integrate_discrete_curve(lengths: np.ndarray, kappas: np.ndarray, torsions: np.ndarray) -> np.ndarray:
    """
    Reconstructs 3D polygonal curve coordinates from intrinsic differential invariants
    using the discrete Frenet-Serret frame integration.
    """
    n_atoms = len(lengths) + 1
    coords = np.zeros((n_atoms, 3), dtype=np.float64)

    # Initial Frenet-Serret frame: Tangent T0, Normal N0, Binormal B0
    T = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    N = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    B = np.array([0.0, 0.0, 1.0], dtype=np.float64)

    coords[0] = [0.0, 0.0, 0.0]
    coords[1] = coords[0] + lengths[0] * T

    for i in range(len(kappas)):
        kap = kappas[i]
        tau = torsions[i] if i < len(torsions) else 0.0

        # Rotate tangent in (T, N) plane by turning angle kappa
        T_new = np.cos(kap) * T + np.sin(kap) * N
        norm_t = np.linalg.norm(T_new)
        if norm_t > 0:
            T_new /= norm_t

        # Normal vector before torsion rotation
        N_prime = -np.sin(kap) * T + np.cos(kap) * N
        norm_np = np.linalg.norm(N_prime)
        if norm_np > 0:
            N_prime /= norm_np
        B_prime = np.cross(T_new, N_prime)

        # Rotate osculating plane around T_new by dihedral angle tau
        N_new = np.cos(tau) * N_prime + np.sin(tau) * B_prime
        norm_nn = np.linalg.norm(N_new)
        if norm_nn > 0:
            N_new /= norm_nn
        B_new = np.cross(T_new, N_new)

        T, N, B = T_new, N_new, B_new
        coords[i + 2] = coords[i + 1] + lengths[i + 1] * T

    return coords.astype(np.float32)


def write_dcd(filename: str, trajectory: np.ndarray, step_interval: int = 100, timestep: float = 0.0416):
    """
    Writes a trajectory of shape (F, N, 3) to standard binary DCD format using Fortran records.
    """
    n_frames, n_atoms, _ = trajectory.shape

    with open(filename, "wb") as f:
        # 1. Header Record (84 bytes)
        rec1 = bytearray(84)
        rec1[0:4] = b"CORD"
        struct.pack_into("<i", rec1, 4, n_frames)
        struct.pack_into("<i", rec1, 8, 0)                  # start_step
        struct.pack_into("<i", rec1, 12, step_interval)      # step_interval
        struct.pack_into("<i", rec1, 16, n_frames * step_interval) # n_steps
        struct.pack_into("<i", rec1, 36, 0)                  # has_unit_cell = False
        struct.pack_into("<f", rec1, 40, timestep)           # integration timestep
        struct.pack_into("<i", rec1, 80, 24)                 # charmm_version
        f.write(struct.pack("<I", 84) + rec1 + struct.pack("<I", 84))

        # 2. Title Record
        title_text = b"TopoFold Bistable Allosteric Benchmark Trajectory (Phase 4)".ljust(80, b" ")
        rec2 = struct.pack("<i", 1) + title_text
        f.write(struct.pack("<I", len(rec2)) + rec2 + struct.pack("<I", len(rec2)))

        # 3. Atom Count Record
        rec3 = struct.pack("<i", n_atoms)
        f.write(struct.pack("<I", 4) + rec3 + struct.pack("<I", 4))

        # 4. Coordinate Frames
        for frame_idx in range(n_frames):
            frame = trajectory[frame_idx]
            for axis in range(3):
                axis_coords = [float(frame[i, axis]) for i in range(n_atoms)]
                coord_bytes = struct.pack(f"<{n_atoms}f", *axis_coords)
                rec_len = len(coord_bytes)
                f.write(struct.pack("<I", rec_len) + coord_bytes + struct.pack("<I", rec_len))


def write_pdb(filename: str, coords: np.ndarray, chain_id: str = "A"):
    """
    Writes a 3D coordinate array of shape (N, 3) to a standard PDB file.
    """
    amino_acids = ["GLY", "ALA", "VAL", "LEU", "ILE", "PHE", "PRO", "TYR", "TRP", "SER",
                   "THR", "CYS", "MET", "ASN", "GLN", "LYS", "ARG", "HIS", "ASP", "GLU"]

    with open(filename, "w") as f:
        f.write(f"HEADER    TOPOFOLD BENCHMARK REFERENCE STRUCTURE\n")
        for i in range(len(coords)):
            res_name = amino_acids[i % len(amino_acids)]
            x, y, z = coords[i]
            line = (
                f"ATOM  {i+1:5d}  CA  {res_name} {chain_id}{i+1:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 10.00           C\n"
            )
            f.write(line)
        f.write("END\n")


def generate_bistable_ensemble(
    n_frames: int = 2000,
    n_atoms: int = 60,
    loop_start: int = 25,
    loop_end: int = 35,
    random_seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Generates the bistable conformational ensemble.

    Parameters:
    -----------
    n_frames : int, default=2000
        Number of frames (1000 in State A, 1000 in State B).
    n_atoms : int, default=60
        Chain length.
    loop_start : int, default=25
        Start residue of functional cryptic loop.
    loop_end : int, default=35
        End residue of functional cryptic loop.
    """
    np.random.seed(random_seed)
    trajectory = np.zeros((n_frames, n_atoms, 3), dtype=np.float32)
    labels = np.zeros(n_frames, dtype=int)  # 0: State A, 1: State B

    half = n_frames // 2

    # Loop invariant index ranges:
    # Curvatures kappa_i: residue index i corresponds to array index i-1.
    # Residues loop_start..loop_end (inclusive) correspond to:
    # Internal vertices in loop: loop_start..loop_end -> kappa indices: loop_start - 1 .. loop_end - 1
    # Internal segments in loop: loop_start..loop_end-1 -> tau indices: loop_start - 1 .. loop_end - 2
    kap_loop_slice = slice(loop_start - 1, loop_end)
    tau_loop_slice = slice(loop_start - 1, loop_end - 1)

    print(f"Generating {n_frames} frames (60 residues each)...")
    print(f"  - Frames 0..{half - 1}: State A (Closed Cryptic Pocket, Helical Loop Turn)")
    print(f"  - Frames {half}..{n_frames - 1}: State B (Open Cryptic Pocket, Extended Hairpin)")

    for f in range(n_frames):
        is_state_b = f >= half
        labels[f] = 1 if is_state_b else 0

        # Physical peptide bond lengths: 3.80 +/- 0.02 Å
        lengths = 3.80 + np.random.normal(0.0, 0.015, n_atoms - 1)

        # Baseline curvature across the chain: alpha-beta conformational spectrum
        kappas = 1.45 + np.random.normal(0.0, 0.08, n_atoms - 2)
        kappas = np.clip(kappas, 0.85, 2.10)

        # Terminal stochastic fluctuations: random walk in torsions
        # Terminal residues 0..24:
        torsions = np.random.uniform(-np.pi, np.pi, n_atoms - 3)

        # Correlated thermal breathing in terminal tails
        tail1_phase = np.sin(f * 0.02) * 0.8
        tail2_phase = np.cos(f * 0.03) * 0.8
        torsions[:loop_start - 1] += tail1_phase
        torsions[loop_end - 1:] += tail2_phase

        # Wrap torsions into (-pi, pi]
        torsions = (torsions + np.pi) % (2 * np.pi) - np.pi

        # Apply bistable state parameters to the core functional loop (25..35)
        if not is_state_b:
            # State A: Closed / compact helical turn
            kappas[kap_loop_slice] = 1.55 + np.random.normal(0.0, 0.025, kap_loop_slice.stop - kap_loop_slice.start)
            torsions[tau_loop_slice] = 0.85 + np.random.normal(0.0, 0.035, tau_loop_slice.stop - tau_loop_slice.start)
        else:
            # State B: Open / extended cryptic pocket turn
            kappas[kap_loop_slice] = 0.95 + np.random.normal(0.0, 0.025, kap_loop_slice.stop - kap_loop_slice.start)
            torsions[tau_loop_slice] = 2.75 + np.random.normal(0.0, 0.035, tau_loop_slice.stop - tau_loop_slice.start)

        coords = integrate_discrete_curve(lengths, kappas, torsions)

        # Add stochastic global rigid-body translation and orientation drift
        # to replicate unconstrained MD simulation coordinates
        drift_translation = np.random.normal(0.0, 5.0, (1, 3))
        coords += drift_translation.astype(np.float32)

        trajectory[f] = coords

    state_a_ref = trajectory[0].copy()
    state_b_ref = trajectory[half].copy()

    return trajectory, state_a_ref, state_b_ref


def main():
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    os.makedirs(data_dir, exist_ok=True)

    dcd_path = os.path.join(data_dir, "bistable_ensemble.dcd")
    state_a_pdb = os.path.join(data_dir, "state_a.pdb")
    state_b_pdb = os.path.join(data_dir, "state_b.pdb")

    # Generate 2,000 frames
    traj, ref_a, ref_b = generate_bistable_ensemble(n_frames=2000, n_atoms=60)

    # Save binary DCD
    print(f"\nWriting binary DCD trajectory to {dcd_path}...")
    write_dcd(dcd_path, traj)
    dcd_size_mb = os.path.getsize(dcd_path) / (1024 * 1024)
    print(f"  -> File size: {dcd_size_mb:.2f} MB")

    # Save reference PDBs
    print(f"Writing State A reference PDB to {state_a_pdb}...")
    write_pdb(state_a_pdb, ref_a)

    print(f"Writing State B reference PDB to {state_b_pdb}...")
    write_pdb(state_b_pdb, ref_b)

    print("\nDataset Generation Complete!")
    print(f"  - Total frames: {traj.shape[0]}")
    print(f"  - Residues per frame: {traj.shape[1]}")
    print(f"  - State A frames: 0..999 (Closed Cryptic Pocket)")
    print(f"  - State B frames: 1000..1999 (Open Cryptic Pocket)")


if __name__ == "__main__":
    main()
