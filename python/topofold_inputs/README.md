# topofold_inputs — Cα input adapter (engineering, optional)

Converts an MDAnalysis `Universe` into Cα coordinates (Å, float64) with explicit residue identity, breaks and
provenance. It reads nothing itself; supported formats are exactly MDAnalysis's. MDAnalysis is **not** a dependency of
the core `topofold` package; this package is imported from `python/` when needed.

**Output** (`CATrajectory`, exported as `<run_id>.npz` + `<run_id>.json`):
- `coordinates_angstrom` (frames, residues, 3), NaN for residues without a Cα;
- frame indices, time in ps (NaN when the format has no time), run id;
- per residue: chain (and whether it came from chainID or segid), segid, resid, insertion code, resname, caller label,
  MDAnalysis residue/atom index, altloc used, status (`present`, `no_CA_atom`, `residue_not_found`);
- breaks (`chain_boundary`, `missing_CA`, `numbering_discontinuity`, `ca_ca_distance` with frames) and the
  half-open `segments` between them;
- source units and reader precision, configuration, altloc decisions, PBC statistics, MDAnalysis warnings, software
  versions, SHA-256 of the source files, the adapter source and the exported npz.

**Policies.** Correspondence to external numbering only as supplied (`--correspondence`, list of
`{chain, resid, icode, label}`); duplicate residue keys, unresolved alternate locations and non-finite coordinates
raise `AdapterError`; `--altloc reject|blank_then_A|highest_occupancy` (default `reject`); `--pbc none|make_whole_sequential`
(default `none`; orthorhombic boxes only; never across a break). GRO input loses chain IDs and insertion codes; such
cases are rejected, not guessed.

**Frame identity (fixed-topology contract).** Atom identity comes from the topology and is assumed identical, in the
same order, in every frame. Multi-model PDB identity records are checked per model and differences raise; for
coordinate-only formats (XTC, DCD, TRR) and in-memory readers this cannot be checked and is reported as
`not_verifiable`. Output order (topology or correspondence) is a separate choice. Readers with unit conversion
disabled are rejected unless already in Å; readers declaring no unit need `assume_angstrom_if_undeclared=True`.

**Continuity.** Adjacency is taken from numbering only, conservatively: 52 → 52A → 52B → 53 is continuous; skipped,
reversed or unusual insertion codes, 52 → 53A, chain/segment changes and all other jumps are breaks. Distance never
establishes continuity; it only adds breaks (Cα–Cα above `--max-ca-ca`).

**Features** (`features.py`, needs the `topofold` extension, not MDAnalysis): κ per residue and τ per hinge from the
existing core `compute_invariants`, computed only within segments and split at degenerate segments
(`--min-segment`) and near-collinear vertices (`--collinear-sin-tol`). Every cell is valid or carries a reason
(`missing_coordinates`, `terminal_or_break`, `degenerate_segment`, `near_collinear`, `run_too_short_for_core`,
`core_error`, `non_finite_core_output`); invalid cells are NaN. The core receives float32 input.

**Comparison** (`compare.py`, NumPy only): histogram Jensen–Shannon distance (base 2, no pseudocounts), linear bins
with explicit under/overflow, circular bins with a stated origin, κ/τ window-matched Cα distances, equal-duration
blocks. Descriptive only; used by `benchmarks/audit/dev_ens_compare/`.

**mdCATH bridge** (`mdcath_bridge.py`, needs h5py and MDAnalysis): one HDF5 replica → in-memory Universe from the
file's own topology. It checks the mapping atom by atom, requires `coords` in Å, preserves source frame indices,
labels time as INFERRED (the files have no timestamps) and records hashes.

**Replica-level comparison** (development and SYNTHETIC evidence only; neither validates the BC scan):
- `compare.replica_comparison` is the complete-separation screen, with categories consistent, replica-dependent and
  inconclusive.
  - REPCHAR-1: no calls under the identical-distribution null; low sensitivity under temporal correlation; calls
    condition-confounded sampling.
- `replica_permutation.permutation_test` is an exact replica-label permutation: 252 allocations for 5 vs 5, α chosen
  by the caller.
  - PERM-1: near-nominal under identical distributions and more sensitive than the screen.
  - NULLSTRESS-1: it also rejects when only temporal correlation differs, so a rejection means the replica
    histograms differ (geometry and/or dynamics).
- Not yet evaluated on real data, and not validated for isolating marginal-geometry differences (the real-data design is deferred).
- `dependent_bootstrap.py` (stationary bootstrap + Kuiper, replicas resampled separately) is an audit prototype. It
  failed SYNTHETIC calibration under unequal temporal correlation (DABOOT-1) and is not recommended for inference.

**Descriptive group comparison** (`ensemble.py`):
`python -m topofold_inputs.ensemble SPEC.json --out DIR [--adapter-python PATH]`.
- **Input:** named groups of members, either MDAnalysis files or mdCATH HDF5 replicas, sharing one residue mapping.
- **Outputs:**
  - residue-mapped CSV tables: per-group distributions, replica variability, pairwise group differences,
    per-member temporal diagnostics, and exclusions with reasons;
  - `distributions.npz`, `report.md` and a hashed `manifest.json`.
- **Scope:** descriptive only, with no p-values, significance labels or uncertainty intervals. Demonstration in
  `benchmarks/audit/dev_ensemble_report/`.
- **v0.2.0 details:**
  - **Circular temporal diagnostic:** C(k) on centered unit vectors, which is independent of the angular origin.
  - **Circular summaries:** report the resultant length; centers with R < 10⁻⁶ are reported as unavailable (a
    numerical rule).
  - **Degenerate series:** single-frame, constant or short series carry explicit reasons.
  - **Pooling:** pooled distributions weight members by their retained frame counts.

**Run** (repository root):

```bash
# whole pipeline: input -> checked coordinates -> kappa/tau -> NPZ + JSON + inspection report
PYTHONPATH=python .venv/bin/python -m topofold_inputs.pipeline TOP [TRAJ ...] --run-id ID --out OUTDIR \
    --adapter-python benchmarks/audit/wpd_local_retrieval/build/venv/bin/python [adapter options]
# SYNTHETIC example input for the command above
PYTHONPATH=python benchmarks/audit/wpd_local_retrieval/build/venv/bin/python python/examples/ca_adapter_example.py --out OUTDIR
# tests
PYTHONPATH=python benchmarks/audit/wpd_local_retrieval/build/venv/bin/python -m pytest tests/test_ca_adapter.py
.venv/bin/python -m pytest tests/test_ca_features.py tests/test_ensemble_compare.py
```

**Not included, and what would activate it.**
- N/C atoms and φ/ψ: a selected comparison that needs the peptide-plane or φ/ψ target (e.g. the WP-D φ/ψ task).
- Side-chain features: an application whose state definition depends on side chains (e.g. kinase DFG χ).
- Streaming: a permitted run whose arrays exceed the memory ceiling (EB-1 R4, 8 GiB).
- New readers: a permitted dataset in a format MDAnalysis cannot read.
- Indexing or optimization: a measured bottleneck in a real workflow.

Engineering component only: tested on SYNTHETIC fixtures, one local PDB read and the DEV-CA-1 development run
(`benchmarks/audit/dev_ca_pipeline/`); no biological validation is implied.
