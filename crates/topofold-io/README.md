# topofold-io

High-throughput molecular dynamics trajectory streaming I/O for TopoFold.

## Supported Formats

- **DCD Trajectory Reader & Writer:** Pure binary, zero-copy streaming for CHARMM, NAMD, and OpenMM trajectories.
- **Multi-Model PDB Streaming:** Memory-efficient chunked parsing for conformational ensembles (NMR models, MD snapshots).
- **C-alpha Atom Filtering:** Directly streams and isolates $C_\alpha$ backbones from full-atom trajectory files.
