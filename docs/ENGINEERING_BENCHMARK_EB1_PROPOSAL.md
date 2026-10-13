# EB-1: Bounded Engineering Benchmark (proposal, 2026-10-12)

**Status:** proposal for review. Nothing in it has been run, and no data will be downloaded.

**Purpose:** decide, from measurements, whether the next engineering effort should go to the **loader**
(ingestion and reconstruction) or to the **engine** (feature extraction and distance/search) for the two retained
candidates:
- **A** (residue-aligned κ–τ), for angular local-conformation retrieval;
- **B2** (window Cα RMSD), for Cartesian peptide-plane retrieval.

The legacy κ–τ Fréchet index and its writhe cascade are **not** targets. They were deprioritized as an engineering
decision (LE-1 E1.3).

## 1. Data (existing development data only)

| Workload | Source | Role |
|---|---|---|
| D1 | AdK `open_start` and `closed_start` (WP-B raw XTC + GRO, SHA-256 bound to the WP-B G0 record) | Development; long trajectories, one 214-residue protein |
| D2 | mdCATH schema sample `1ha8A00` (HDF5, SHA-256 verified) | Development; HDF5 path, small domain, four trajectories |

**Not used:**
- The 30 held-out mdCATH files. They remain held-out evidence; timing on them is not needed for this decision.
- Any new download.

**Scaling series.** Each workload is measured at frame fractions {0.25, 0.5, 1.0} and window sizes W ∈ {5, 8}, so
that time and memory can be fitted against frames × windows.

**Projections.** These are labelled as projections, not measurements, and are made to two declared target workloads:
- **P1:** the completed 30-domain held-out job, re-expressed per stage;
- **P2:** an engineering stress target of 1,000 residues × 10,000 frames per trajectory.

P2 is a stated engineering assumption, not a biological requirement.

## 2. Stages, measured separately

| Stage | What is timed | Implementations compared |
|---|---|---|
| S1 Ingestion | File → required atoms (N, CA, C, O) as stored float32, including decoding | MDAnalysis XTC reader (D1); h5py chunked selected-atom reads (D2); h5py whole-dataset read (D2, context) |
| S2 Reconstruction | Periodic make-whole, burn-in slice, float32 rounding | NumPy sequential minimum image (current harness) |
| S3 Feature extraction | κ/τ, Cα pair vectors, φ/ψ and peptide-plane atoms (targets) | Rust binding `compute_invariants` per frame; vectorized NumPy recomputation |
| S4 Distance/search | Exhaustive top-k (k = 10) over all windows, query vs database | A, B2 (and M1, B1 for context) in vectorized NumPy at batch sizes {64, 128, 256}, with the best batch reported per implementation; the existing Rust legacy `query_subcurve` path reported as context only (different distance) |
| S5 Aggregation | Per-window metrics, domain summaries | Current harness |

**What is reported per stage:**
- wall time (median of ≥ 5 steady-state repeats) and CPU time;
- bytes read (`/proc/self/io`);
- peak RSS (isolated subprocess per stage, minus an import-only baseline);
- tracemalloc peak for Python/NumPy allocations.

Setup (imports, index fit, cache warm-up) is reported separately. Runs are single-threaded and at 4 threads, with
BLAS threads pinned and recorded.

**File cache.** First-read and warm-cache timings are both reported. Caches cannot be dropped without
administrative rights, so "first read" is the first access in a session and is labelled as such.

## 3. Fairness and correctness checks (mandatory; a failure invalidates that stage's timing)

**Fairness:**
- Compare implementations of the same quantity, on the same machine, dtype, thread count and batch policy.
- Do not attribute Rust vs Python differences to the representation.

**Correctness, per stage:**

| Stage | Check |
|---|---|
| S1 | Arrays from each reader are bit-identical to the reference arrays: WP-B derived Cα arrays for D1; the schema-check arrays for D2 |
| S2 | Equal to the frozen harness functions; no backbone bond > 2.0 Å after reconstruction |
| S3 | Binding vs NumPy ≤ 1e-9 rad on identical float32-rounded inputs |
| S4 | Distances agree with the frozen harness ≤ 1e-12 (same formula) or ≤ 1e-6 Å (Horn vs Kabsch); tie-aware top-10 memberships identical across implementations |

## 4. Decision rule (fixed before measurement)

Let **share(Sx)** be a stage's fraction of end-to-end wall time, and **T_P** and **M_P** the projected time and peak
memory for target workload P.

| Order | Condition | Next work |
|---|---|---|
| 1 | Any correctness check fails | Fix or document the failing stage first; no performance decision |
| 2 | share(S1 + S2) ≥ 0.5 on either development workload, or the S1/S2 projected peak memory for P2 exceeds 8 GB | **Loader work:** streaming or chunked readers, reconstruction without full-trajectory copies |
| 3 | share(S3 + S4) ≥ 0.5 and T_P2 > 24 h on one 6-core machine | **Engine work**, limited to A and B2: compiled batched exhaustive top-k, possibly reusing the Rust core |
| 4 | Both 2 and 3 hold | The stage with the larger projected absolute saving on P2 first; both estimates recorded |
| 5 | Neither holds | Neither. Keep the NumPy harness and spend effort on validation (e.g., a second held-out dataset) |

The 0.5 share, 8 GB and 24 h thresholds are **engineering choices**, stated so that the decision is not made after
seeing the numbers. Projections carry the uncertainty of the fitted scaling and are reported with the fit residuals.

## 5. Budget and outputs

- **Budget:** ≤ 2 h CPU; no downloads; no production edits.
- **Outputs:**
  - a frozen EB-1 protocol;
  - a timing harness under `benchmarks/audit/eb1/`;
  - `results/eb1_*.json` with versions and hashes;
  - a report.

## 6. What EB-1 cannot decide

- Retrieval quality. That is fixed by the held-out results and is not re-estimated.
- Whether the loader generalizes to other formats beyond XTC and the mdCATH HDF5 layout.
- Performance of a not-yet-written implementation. Only existing implementations are measured; any proposed engine
  or loader work needs its own correctness and performance acceptance test.
