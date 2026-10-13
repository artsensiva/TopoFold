# EB-1 Amendment 1: Corrected Accounting and Targeted Confirmation (frozen 2026-10-12)

**Amends:** `EB1_ENGINEERING_BENCHMARK_PROTOCOL.md`. That file is unchanged, and so are the original EB-1
measurements, analysis, report and outcome (*insufficient evidence*). Recorded after the original results were seen.
- The amended design responds to the identified specification error and to an independent review.
- Thresholds R1–R4 are **unchanged**.
- Frozen by SHA-256 in `benchmarks/audit/eb1/results/eb1a1_freeze.json` before any amended measurement.

## §1. What required correction

| # | Problem in the original | Correction |
|---|---|---|
| 1 | S4 memory variable used a fixed batch of 128, ignored the actual n_D in batch cells, ignored method-specific temporaries and retained outputs, and pooled A with B2 and single-query with batch | Per method **and** per mode (single-query vs batch), the predictor is P = resident inputs + retained top-k outputs + transient bytes (§3), with b = min(n_Q, 128). Never pooled across methods or modes |
| 2 | The S3 memory cell also ran the NumPy reference; the S2/S3 cells loaded the full cached array and sliced it, so the backing allocation stayed resident; S4 cells loaded arrays the method does not use | Operational cells load **owned copies** of only the frames they use (`owned_subset`). The S3 cell never runs a reference. References run in a separate correctness subprocess. S2 and S3 are **re-measured** |
| 3 | B2 was charged for κ/τ extraction (S3) | Method-specific paths: A = S1 → S2 → S3 → S4_A; B2 = S1 → S2 → S4_B2 |
| 4 | No held-out confirmation of the fitted models | Five confirmation workloads per method, excluded from fitting (§4), each with a prediction-error criterion |
| 5 | Incidental harness work was not separated | S4 timing is decomposed into per-call slice copies (incidental harness overhead), distance kernel, top-k and the remaining loop. This instrumentation performs the identical computation (tested). Both total and decomposition are reported |

## §2. Reused measurements (with justification)

- **S1 ingestion** time and peak RSS are reused from the original run. The ingestion process did not use cached
  slices or references, so it is unaffected by issues 1–3. Only its memory predictor is corrected:
  - D1 (XTC): the output Cα array (n × L × 12 B); the MDAnalysis topology and frame buffer are constant.
  - D2 (HDF5): the largest per-trajectory N/CA/C/O read (F × L × (9 + 3) × 8 B), plus all retained Cα outputs,
    plus one chunk buffer (64 × 4L × 12 B).
- **Correctness C0–C3** from the original run are reused for the original paths; new paths get C1′–C3′ (§5).
- Everything else is re-measured.

## §3. Memory and work accounting

**Resident inputs.** These are owned copies of the frames used.
- A: n_frames × (2L − 5) × 8 B (κ, τ).
- B2: n_frames × 3L × 8 B (Cα).
- S2: n × L × (12 + 24 + 12 + 24) B (stored, whole, float32 copy, rounded).
- S3: n × L × (24 + 12) + n × (2L − 5) × 8 B.

**Retained outputs.** n_Q × n_win × k × 8 B.

**Transient bytes per kernel call.** α_m · b · n_D + β_m · n_D, with α and β fixed by the pre-freeze SYNTHETIC
tracemalloc calibration `results/eb1a1_temp_calibration.json` (W = 8; maximum calibration error is recorded there).

**Work.**
- Comparisons C = n_Q · n_win · n_D.
- Kernel calls = n_win · ⌈n_Q / b⌉.
- Call database rows = calls · n_D.

## §4. Targeted design (W = 8)

**Data.**
- D1: AdK `open_start` (database frames 1238–2252, n_D ≤ 1015; queries from frames 0–1012).
- D1+: the same database plus all 1,678 `closed_start` frames, n_D = 2,693. These are **real** frames of the same
  214-residue protein, read through the operational S1/S2/S3 path from the WP-B raw XTC, which is SHA-256 bound to
  G0. Engineering use only: no retrieval-quality claim.
- D2: mdCATH schema sample (development data).
- No constructed repetition is used.

**Fit set** (per method; D1):

| Sweep | n_Q | n_D | n_win |
|---|---|---|---|
| Batch, n_Q | {32, 64, 128, 192, 256} (below, at and above b = 128) | 1015 | 20 |
| Batch, n_D | {128, 256} | {250, 500} | 20 |
| Batch, n_win | 64 | 1015 | {5, 80} |
| Single-query (n_Q = 1, batch 1) | 1 | {250, 500, 1015} | {20, 207} |

**Confirmation set** (excluded from every fit):

| ID | Data | Mode | n_Q | n_D | n_win |
|---|---|---|---|---|---|
| CW1 | D1 | batch | 512 | 1015 | 207 |
| CW2 | D1+ | batch | 256 | 2693 | 40 |
| CW3 | D1+ | single | 1 | 2693 | 207 |
| CW4 | D2 | batch | 256 | 405 | 44 |
| CW5 | D2 | single | 1 | 405 | 44 |

**Repeats.** Each S4 cell is timed 3 times when its single run is under 5 s, else once. Each runs in a fresh
subprocess (BLAS threads 1, `RLIMIT_AS` = 6 GB).

**S2 and S3 re-measurement** (S3 only for A): D1 n ∈ {563, 1126, 2253}; D2 n ∈ {405, 810, 1620}; 3 repeats each.

## §5. Correctness (separate subprocess; mandatory)

| Check | What it compares | Criterion |
|---|---|---|
| C1′ | `closed_start` operational reconstruction vs WP-B derived `closed_start_ca_whole.npy` | ≤ 1e-6 Å |
| C2′ | Binding vs NumPy features on the closed frames | ≤ 1e-9 rad |
| C3′ | D1+ search: 20 queries × 3 windows per method, instrumented search vs per-pair reference (A: `d_aligned`; B2: Kabsch) | Original tolerances and near-tie rule |
| C4′ | Instrumented search outputs identical to `eb1_lib.search` on CW4 | Identical arrays |

## §6. Models, support and confirmation

**Time models, per method:**
- Batch: t = c0 + c1·C + c2·(call database rows).
- Single: t = c0 + c1·C + c2·n_win.
- If c2 < 0, the model is refit without the c2 term (pre-declared).

**Memory model, per method and mode:** absolute peak RSS = m0 + m1·P.

**S2/S3 time and memory:** linear in work and in P respectively (original §6).

**A model is supported only if all of:**
1. ≥ 3 distinct fit sizes;
2. the leading coefficient is > 0;
3. every fit relative residual is ≤ 0.25;
4. every confirmation workload of that method/mode is predicted within **25 %** relative error.

The prediction error on each confirmation workload is always reported. P2 extrapolation factors are reported, and
P2 projections remain extrapolations beyond the largest confirmed size.

## §7. Amended requirements and decision

P2, the limits and the §7 logic of the original protocol are unchanged. Method-specific composition:

| Requirement | A | B2 |
|---|---|---|
| R1 batch | S3 (11,000 frames, CPU) + S4_A batch (P2, CPU) | S4_B2 batch (P2, CPU) |
| R2 single-query | S3 (1 frame) + S4_A single (P2) | S4_B2 single (P2) |
| R4 S4 memory | per method, batch and single | per method, batch and single |

**P2 memory predictor.**
- Resident inputs: 11,000 frames at L = 1,000.
- Retained outputs: n_Q = 1,000 (batch) or 1 (single), n_win = 993.
- Transient term: at b = min(n_Q, 128) and n_D = 10,000.

R3 and the S1/S2 memory requirements use the S1 reuse in §2 and the S2 re-measurement.

**Outcome:** insufficient evidence if any mandatory check fails, any planned cell is missing, or any required
projection is unsupported (including a failed confirmation); otherwise loader / engine / both / neither by original
§7.

## §8. Budget and stop rules

- **CPU:** aggregate including the original run ≤ 2 CPU-h. The amendment runner stops at a cumulative 1.8 CPU-h,
  counting the original 0.072 h.
- **Memory:** before each cell, MemAvailable ≥ 500 MB (wait up to 10 min, then skip → missing); `RLIMIT_AS` = 6 GB.
- No desktop process is touched.
- No change after observing performance.
