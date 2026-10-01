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
│   ├── topofold-core/     # Foundational mathematical kernel and SE(3) invariants
│   ├── topofold-index/    # Vantage-Point Tree metric indexing and Fréchet search
│   ├── topofold-io/       # High-throughput trajectory streaming (DCD, Multi-Model PDB)
│   └── topofold-python/   # PyO3 Python bindings (import topofold)
└── tests/                 # Integration and Python test suites
```

## Guarantees

- **Strict Memory Safety:** 100% safe Rust (`#![forbid(unsafe_code)]` across all Rust core/index/io crates).
- **Mathematical Invariance:** Curvature and torsion are invariant under $\text{SE}(3)$ (rotations $R \in \text{SO}(3)$ and translations $t \in \mathbb{R}^3$) down to machine precision ($< 10^{-12}$).
- **Chiral Sensitivity:** Correctly distinguishes mirror reflections ($O(3) \setminus \text{SO}(3)$), preserving the handedness of $\alpha$-helices and $\beta$-sheets.

## Building & Testing

### Rust Workspace

```bash
# Build the workspace in release mode
cargo build --release

# Run all unit, integration, and doc tests across all crates
cargo test --workspace

# Run Clippy linter
cargo clippy --workspace --all-targets
```

### Python Bindings

```bash
# Build and install editable extension via maturin
maturin develop

# Run Python test suite
pytest tests/test_python_bindings.py
```

## Quick Start (Python)

```python
import numpy as np
import topofold as tf

# 1. Parse PDB C-alpha coordinates
coords = tf.read_pdb("protein.pdb", chain="A")

# 2. Compute SE(3)-invariant geometric primitives
kappa, tau, writhe = tf.compute_invariants(coords)

# 3. Stream trajectory (DCD or Multi-Model PDB) directly into metric index
index = tf.ConformationalIndex.from_dcd(
    "simulation.dcd", 
    ca_indices=[0, 4, 8, 12, ...],
    window_size=10, 
    coarse_radius=0.5
)

# 4. Search for k-nearest conformers (global backbone or localized loop/cryptic pocket)
hits = index.query(coords, k=10)
loop_hits = index.query_subcurve(start_res=120, end_res=135, query_coords=coords, k=5)
```

## Benchmark: TopoFold vs. Cartesian PCA

In molecular dynamics ensembles, localized allosteric changes (such as cryptic pocket openings) are frequently obscured in linear Cartesian PCA due to high-amplitude thermal fluctuations in flexible terminal tails.

Using a synthetic bistable 60-residue trajectory ($N = 2,000$ frames) exhibiting a 10-residue loop transition (Closed State A $\leftrightarrow$ Open State B) amidst high Brownian terminal noise:

| Evaluation Criterion                  | Cartesian PCA (Kabsch Aligned) | TopoFold Intrinsic Subcurve Index |
|:--------------------------------------|:-------------------------------|:-----------------------------------|
| **Mathematical Formulation**          | Global Euclidean $\mathbb{R}^{3N}$ SVD | $SE(3)$-Invariant $(\kappa, \tau)$ Fréchet  |
| **Silhouette Score ($S$)**            | **0.3374** (Severe overlap)    | **0.9850** (Near-ideal clustering) |
| **Separation Quality**                | Overlapping / Obscured         | Complete Bimodal Separation        |
| **Sensitivity to Terminal Noise**     | **Severe** ($> 80\%$ variance) | **Zero** (Strictly Localized)      |
| **Coordinate Superposition Needed**   | **Yes** ($O(M \cdot N)$ Kabsch)| **None** (Intrinsic Geometry)     |
| **Sub-second Pocket Search**          | **No** (Recomputes Covariance) | **Yes** (6.9 µs per frame)         |
| **Memory Ingestion**                  | Full Cartesian RAM buffer      | Zero-Copy Trajectory Stream (DCD)  |

![TopoFold vs Cartesian PCA Benchmark](assets/benchmark_pca_vs_topofold.png)

To reproduce the benchmark:
```bash
python benchmarks/run_benchmark.py
# Or open the interactive Jupyter notebook:
# jupyter notebook benchmarks/benchmark_pca_vs_topofold.ipynb
```

## License

Licensed under either of:
- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.
