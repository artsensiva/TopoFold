# EB-1 Engineering Benchmark: Frozen Protocol (2026-10-12)

**Status:** FROZEN before any EB-1 timing was collected. SHA-256 in `benchmarks/audit/eb1/results/eb1_protocol_freeze.json`.
Supersedes the proposal `ENGINEERING_BENCHMARK_EB1_PROPOSAL.md`, which is preserved unchanged.

**Scope:** engineering only. It does not re-estimate retrieval quality or reopen method selection; WP-D results and the
2026-10-12 addendum stand. It uses existing development data only: no held-out file is read, nothing is downloaded,
and no production code is changed.

## 1. Operational vs research-only paths

| Path | Stages | Decides? |
|---|---|---|
| **Operational A/B2 retrieval** | S1 ingestion → S2 reconstruction → S3 features → S4 exhaustive top-k search | **Yes** |
| Research-only | R1 φ/ψ and peptide-plane targets; R2 tie-aware NR metrics; R3 domain bootstrap | No (reported only) |

### Operational stage definitions

**S1. Ingestion.** Stored file → Cα coordinates as stored (float32, Å), residues in ascending resid order. Readers:
- AdK XTC: MDAnalysis 2.10.0. The XTC holds the full solvated system, so every atom is decoded.
- mdCATH HDF5: h5py 3.16.0. Selected-atom chunked reads with the frozen `mdcath_reader_r1`, which reads N, CA, C and O,
  the only existing HDF5 path.

**S2. Reconstruction.** Sequential minimum image along the Cα chain (WP-B `adk_lib.make_whole`), the burn-in slice
where it applies, then float32 rounding of the whole Cα coordinates, held in float64.

**S3. Features.**
- A: κ (turning angles) and τ (dihedrals) per frame, from the Rust binding `compute_invariants` on float32 Cα.
- B2: window Cα coordinates; a slice, no computation.

**S4. Exhaustive top-k search** (k = 10). For each window, distances from queries to all database frames, then
`argpartition` top-k.
- A: `legacy_lib.d_aligned` (NumPy).
- B2: `wpd_lib.horn_rmsd` (NumPy).
- Batches of 128 queries (the WP-D setting).

### Retrieval modes and work units

- **Window:** W consecutive residues, start s = 0 … L − W, so n_win = L − W + 1. W ∈ {5, 8}.
- **Single-query, all windows:** one query frame against n_D database frames in every window.
  Comparisons C = n_win · n_D. This is the latency case.
- **Batch:** n_Q query frames against n_D database frames in every window. Comparisons C = n_Q · n_win · n_D.
- The search cost model is linear in C. Ingestion and feature extraction are modeled separately, linear in decoded
  atom-frames (S1) and in frames × residues (S2, S3).

## 2. Workloads

### Measured (development data)

| ID | Data | L | Database | Queries | n_win (W5 / W8) |
|---|---|---|---|---|---|
| D1 | AdK `open_start` (WP-B raw XTC, SHA-256 bound to WP-B G0) | 214 | frames 1238–2252 (n_D = 1015) | from frames 0–1012 | 210 / 207 |
| D2 | mdCATH schema sample `1ha8A00` (development; SHA-256 verified), 320 K | 51 | replica 1 after 10 % burn-in (n_D = 405) | replica 0 after burn-in (up to 405) | 47 / 44 |

### Declared target (projection only)

| ID | Definition |
|---|---|
| **P2** | L = 1,000 residues; W = 8 (n_win = 993); database n_D = 10,000 frames; batch n_Q = 1,000 frames; single-query = 1 frame |

P2 ingestion is projected for two storage layouts:
- HDF5 protein-only: 16,000 atoms per frame, 11,000 frames.
- XTC full system: 160,000 atoms per frame, assuming 10× solvent, 11,000 frames.

P2 is an engineering stress assumption, not a scientific requirement.

## 3. Requirements (absolute)

| ID | Requirement | Component |
|---|---|---|
| R1 | P2 batch search: S3 + S4 CPU time ≤ **96 CPU-hours** (single-process time summed; no parallel credit assumed) | Engine |
| R2 | P2 single-query, all windows: query S3 + S4 wall time ≤ **10 s**, single process | Engine |
| R3 | P2 loading of all 11,000 frames (S1 + S2) ≤ **2 h** wall, single process, per layout | Loader |
| R4 | Absolute peak RSS of each operational stage at P2 ≤ **8 GB** | Stage owner (S1/S2 → loader; S3/S4 → engine) |

Absolute peak RSS is measured for every cell. Baseline-adjusted RSS (minus an import-only process) is supplementary.

## 4. Correctness references, one per processing boundary (mandatory)

| Boundary | Compared object | Reference | Pass criterion |
|---|---|---|---|
| C0 stored | Operational Cα read (D1: MDAnalysis CA selection; D2: chunked selected-atom h5py read) | Same stored file read independently. D1: MDAnalysis full-system positions indexed at the recorded Cα atom indices. D2: whole `coords` dataset read, then indexed. *Both are same-library consistency checks; no second XTC decoder is available offline* | Bit-identical float32 |
| C1 reconstructed | Operational Cα make-whole output (float64, before rounding) | D1: WP-B derived `open_start_ca_whole.npy` (independent pipeline, hash-recorded). D2: frozen harness `make_whole_r1` on N–CA–C, Cα column | max \|Δ\| ≤ 1e-6 Å |
| C2 features | Binding κ/τ (float32 input) | NumPy turning-angle/dihedral on the identical float32-rounded coordinates (float64) | ≤ 1e-9 rad (τ circular) |
| C3 search | Operational batched top-10 and distances, on 20 seeded queries × 5 seeded windows per dataset × W × method | Per-pair recomputation: A with the frozen `d_aligned` on single pairs; B2 with Kabsch SVD per pair | A distance \|Δ\| ≤ 1e-10; B2 \|Δ\| ≤ 1e-6 Å |

**Comparison rules:**
- Stored coordinates are never compared with reconstructed ones.
- At every boundary the record covers: atoms (Cα; plus N, C, O for the research targets), order (ascending resid),
  frames, units (Å), dtype (float32 stored, float64 compute) and PBC treatment (sequential minimum image along the
  Cα chain; along N–CA–C for the research path).

**Top-k and near ties.** Let d_(k) and d_(k+1) be the 10th and 11th smallest distances.
- The boundary is a **near tie** if d_(k+1) − d_(k) ≤ 1e-9 · (1 + d_(k)).
- Non-tie queries need identical top-10 index sets.
- Near-tie queries are compared by tie-aware membership vectors (WP-D rule), which must agree to ≤ 1e-12. Their count
  is reported.

Any mandatory failure stops the benchmark, and the outcome is **insufficient evidence**.

## 5. Minimum measurement matrix (complete before any decision)

Every cell runs in a **fresh subprocess, single process, BLAS threads = 1**. Each cell records wall time, CPU time,
absolute peak RSS (`ru_maxrss`), bytes read (`/proc/self/io`) and MemAvailable at start.

| Block | Cells | Sizes | Repeats |
|---|---|---|---|
| M0 | Import-only baseline | — | 3 |
| M1 S1 ingestion | D1 frames {563, 1126, 2253}; D2 trajectories {1, 2, 4} (320/0, 320/1, 348/0, 348/1) | 6 | 3 (the first is reported separately as first access) |
| M2 S2 reconstruction | Same sizes as M1 | 6 | 3 |
| M3 S3 features | D1 frames {563, 1126, 2253}; D2 frames {1, 2, 4} trajectories' frames | 6 | 3 (plus 1 warm-up excluded) |
| M4 S4 single-query, all windows | D1, D2 × W {5, 8} × {A, B2} × n_D fraction {¼, ½, 1} | 24 | 5 query frames each |
| M5 S4 batch | D1, D2 × W {5, 8} × {A, B2} × n_Q {32, 64, 128}; n_D full; windows: D1 every 10th start, D2 all | 24 | 1 |
| M6 research-only | D1 W8, n_Q = 128, every 10th window: R1 targets, R2 metrics; R3 bootstrap on SYNTHETIC 30-domain × 7-method arrays, 2,000 resamples | 3 | 1 |
| M7 correctness | C0–C3 as in §4 | — | 1 |
| M8 parallel efficiency (descriptive) | D1 W8 batch with 4 workers vs 1 | 1 | 1. **Only if MemAvailable ≥ 2 GB**; never decisive |

## 6. Projection and support rules

**Time.** Fit t = c0 + c1 · x by least squares over the measured sizes of a block. The work variable x is:
- S1: decoded atom-frames;
- S2, S3: frames × residues;
- S4: comparisons C.

A projection is **supported** only if:
1. there are ≥ 3 sizes;
2. c1 > 0;
3. every relative residual is ≤ 0.25.

The extrapolation factor x_P2 / x_max is reported.

Single-query (R2) uses the M4 fit for the dataset and W = 8. Batch (R1) uses the M5 fit. The P2 projection uses the
larger of the D1- and D2-based projections (conservative).

**Memory.** Fit absolute peak RSS = m0 + m1 · b, where b is the stage's resident array bytes from analytic shapes,
with the same support rules. Projected peak RSS at P2 = m0 + m1 · b_P2.

**Unsupported projections.** If a projection needed by R1–R4 is unsupported, the outcome is **insufficient
evidence**. The cause is stated.

## 7. Decision rule (applied in order)

1. A mandatory correctness failure, a missing matrix cell (M0–M7), or an unsupported required projection →
   **insufficient evidence**.
2. **Loader flagged** if R3 is violated for either layout, or R4 is violated by S1 or S2.
3. **Engine flagged** if R1 or R2 is violated, or R4 is violated by S3 or S4.
4. Outcome:
   - neither flagged → **neither**;
   - one flagged → that one;
   - both flagged → **both**, ordered by the larger ratio of projected value to limit first.

Projections for implementations that do not exist are never reported as measured savings. The outcome names which
component to build, not how much it would save.

## 8. Resource limits and stop rules

- **Budget:** ≤ 2 aggregate CPU-hours, including every subprocess. The runner tracks self and children CPU time.
- **CPU stop:** if the cumulative total reaches 1.8 CPU-h, the remaining cells are not run, and missing cells yield
  insufficient evidence.
- **Memory:** each cell has `RLIMIT_AS` = 6 GB.
  - Before each cell, if MemAvailable < 500 MB, wait up to 10 min, then skip the cell (it becomes missing).
  - A cell exceeding its limit is recorded as failed, which also means insufficient evidence.
- **No silent changes:** the matrix, sizes and rules above are not changed after observing performance. Any amendment
  is recorded separately, before any affected measurement, with its reason.

## 9. Environment disclosure

- Shared desktop machine: 6 cores, 15 GB RAM.
- At freeze, other user applications held most memory: about 1.1 GB available and swap full.
- Timings and RSS are therefore conditional on contention, and page caches cannot be dropped.
- Cache state is recorded per cell (first access vs repeated access).
