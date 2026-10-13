# topofold-io

Molecular dynamics trajectory I/O for TopoFold.

## Supported formats

- **DCD reader and writer** (CHARMM, NAMD, OpenMM). The Rust reader iterates frame by frame
  (`next_frame`), and each frame is copied into owned buffers. It is not zero-copy or
  memory-mapped. The Python `read_dcd` binding collects all frames into one array.
- **Multi-model PDB.** `read_pdb_trajectory` parses all models into memory and returns them as
  a `Vec`.
- **Cα selection.** Cα (and optionally Cβ) atoms are selected from full-atom frames.

Not supported: XTC/TRR and other GROMACS formats; chain-break detection (breaks are treated as
bonds; see `docs/KNOWN_ISSUES.md`, m1). Throughput has not been benchmarked against other
readers.
