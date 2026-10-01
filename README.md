# TopoFold Engine

High-performance, SE(3)-invariant discrete differential geometry engine for protein conformational trajectories and metric indexing.

## Motivation & Scientific Foundations

In structural biology and molecular dynamics (MD), conventional conformational analysis relies on Cartesian coordinates and linear Principal Component Analysis (PCA). This paradigm suffers from two fundamental limitations:
1. **Rotational & Translational Confounding:** Cartesian coordinates require global structural alignment (Kabsch/RMSD superpositions), which is expensive ($O(M \cdot N)$ across $10^5$ frames) and prone to domain-motion artifacts.
2. **Loss of Non-Linear Transitions:** Cryptic pockets and allosteric state transitions often involve localized torsional hinging. Linear PCA smears these transitions into overlapping clouds in low-dimensional projections.

**TopoFold** models the protein backbone ($\text{C}_\alpha$ trace) as an intrinsic 3D polygonal space curve using **discrete differential geometry**:
- **Segment Lengths ($l_i$):** Distance between consecutive $\text{C}_\alpha$ atoms (${\approx}3.80\text{ \AA}$).
- **Discrete Curvature ($\kappa_i$):** Turning angle between consecutive segment vectors in $[0, \pi]$.
- **Discrete Torsion ($\tau_i$):** Signed dihedral angle between consecutive osculating planes in $(-\pi, \pi]$.

By the Fundamental Theorem of Discrete Space Curves, the sequence $(l_i, \kappa_i, \tau_i)$ uniquely determines the 3D backbone conformation up to global rigid-body motion ($\text{SE}(3)$) without requiring any coordinate superposition.

## Project Structure

```text
TopoFold/
├── crates/
│   ├── topofold-core/     # Foundational mathematical kernel and invariants
│   ├── topofold-io/       # (Phase 2) Trajectory streaming (XTC, DCD, PDB)
│   ├── topofold-index/    # (Phase 2) Metric indexing (VP-tree, Cover Tree)
│   └── topofold-python/   # (Phase 3) PyO3 Python bindings (import topofold)
└── tests/
```

## Guarantees

- **Strict Memory Safety:** 100% safe Rust (`#![forbid(unsafe_code)]`).
- **Mathematical Invariance:** Curvature and torsion are invariant under $\text{SE}(3)$ (rotations $R \in \text{SO}(3)$ and translations $t \in \mathbb{R}^3$) to floating-point precision ($< 10^{-6}$).
- **Chiral Sensitivity:** Correctly distinguishes mirror reflections ($O(3) \setminus \text{SO}(3)$), preserving the handedness of $\alpha$-helices and $\beta$-sheets.

## Building & Testing

```bash
# Build the workspace
cargo build --release

# Run unit tests and SE(3) invariance validation
cargo test -p topofold-core
```

## License

Licensed under either of:
- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.
