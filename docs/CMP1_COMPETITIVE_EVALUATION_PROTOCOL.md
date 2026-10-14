# CMP-1 Bounded Competitive Evaluation of Window Cα RMSD Retrieval: Protocol (frozen 2026-10-13)

**Supersedes:** `docs/B2_COMPETITIVE_EVALUATION_PROTOCOL_DRAFT.md`, which is kept as the reviewed draft.
**Starting checkpoint:** `f55def1`, plus the uncommitted BQ-3 closeout.
**Scope:** harness only. No production integration, held-out scoring, raw-data download, commits or pushes.

**Disclosure:** this protocol was written after **adapter-development checks**
(`benchmarks/audit/cmp1/results/cmp1_adapter_dev_checks.json`). Those checks use queries, sizes and seeds separate
from the deciding workloads below; they verified library semantics and changed the F32 tolerance (§4). No deciding
measurement had been taken when this protocol was frozen.

**Claim boundary:** CMP-1 compares four implementations on five workloads on one machine. It is not a field-wide
best-tool claim. If an existing library is preferable for this task, recommending it is a valid outcome; it is not a
failure of the broader TopoFold mission.

## §1. Task and output contract

**Inputs:** one query frame, n_D database frames, an explicit window list, W = 8 and k = 10.

**Distance:** for every window and database frame, the minimal Cα RMSD over **proper rotations**, after independently
removing each window's centroid. Each window is fitted on its own; no global or group fit is used. Unweighted. Search
is exhaustive.

**Output per query** (retained):
- top-k database positions, int64, of shape (n_win, k);
- the top-k distances, float64 Å, of shape (n_win, k).

The top-k order is ascending (distance, position), computed by one shared `topk_lex` for every implementation, with
its cost included.

**Per-window distance vectors** (n_D) are transient in timed runs; each is discarded after its top-k. Full
(n_win, n_D) distance arrays are retained only in verification processes, outside timing.

**Difference from BQ-3:**
- BQ-3 retained only top-k indices.
- CMP-1 retains indices **and** distances.
- Output identity hashes cover both: a timed run's hash must equal the gate's hash for the same implementation and
  query.

**Memory accounting (per fresh process):**
- lifetime peak RSS;
- RSS after imports and warm-up (`rss_after_setup`);
- RSS after loading and after preparation + query 1;
- the exact bytes of prepared structures (V2b: centred windows + G_B; MDTraj: per-window `xyz` + cached traces;
  MDAnalysis: centred float64 windows);
- retained output bytes.

The memory criterion uses **incremental** peak RSS (peak − `rss_after_setup`), so that library import footprints do
not dominate; absolute peaks are reported too.

## §2. Implementations (pinned; semantics verified in installed source)

| ID | Implementation | Precision | Preparation (equal opportunity) |
|---|---|---|---|
| T-V2b | `b2opt_lib` V2b (unchanged, hash-checked) | float64 | centred windows + G_B |
| T-V0 | EB-1 harness kernel (internal baseline) | float64 | none |
| E-MDT | MDTraj **1.11.2** `md.rmsd(..., precentered=True, parallel=False)` | float32 coordinates and distances; Theobald QCP; polynomial coefficients in float32, quartic solved in double (`theobald_rmsd.cpp`) | full trajectory in nm float32 → `atom_slice` per window → `center_coordinates()` (cached traces) |
| E-MDA | MDAnalysis **2.10.0** `lib.qcprot.CalcRMSDRotationalMatrix` on centred float64 windows, one call per pair | float64 inputs, QCP in double; Newton (evalprec 1e-14, ≤ 50 iterations), no eigensolver fallback | centred float64 windows (`qcprot` requires centred inputs) |

**Semantics verified in the pinned MDTraj source and by adapter-development checks:**
- `precentered` is honoured only when `atom_indices` is None. With a subset it is ignored with a warning, so
  pre-slicing is required.
- The non-precentered `atom_indices=None` path centres the input in place.
- Pre-sliced precentered windows agree with subset RMSDs to about 1e-6 nm.

**MDAnalysis:** the local-window comparator is the independent per-pair QCP fit above. Globally fitted group RMSDs
(`analysis.rms.RMSD` group selections) are **not** used.

**CPPTRAJ is excluded from this run.**
- It needs AmberTools (conda/source), which is outside the granted authorization (MDTraj and its dependencies only).
- Its 2-D RMS tooling would need separate verification of per-window semantics.
- Excluding it means CMP-1 does not cover every strong C++ implementation. This limitation is stated in the result.

**Environment:** isolated venv `benchmarks/audit/cmp1/build/venv` (Python 3.14.4; NumPy 2.5.3; SciPy 1.18.1), installed
offline from a hashed wheelhouse. Commands, wheel hashes, the MDTraj PyPI digest match and inspected source hashes are
in `results/cmp1_environment.json`. Existing audit environments are untouched.

## §3. Canonical input, units and immutability

**Canonical physical input:** one float64 array in Å (`cmp1_data.load`):
- Real workloads use the S2 rounded arrays, verified float32-representable.
- SYNTHETIC_tiled rows (BQ-2 builder) are recast once through float32, a documented change.

The array is passed **read-only**; any in-place write raises. Each process also checks its input hash before and
after.

**Conversions:**
- MDTraj receives float32 nm: (x · 0.1) cast to float32. Its distances are converted back to float64 Å by × 10.
  Representation error is recorded separately: the reference recomputed on the converted representation versus the
  canonical reference.
- All other implementations receive float64 Å.

## §4. Accuracy tracks (frozen; identical requirements within each track)

**Reference:** float64 batched SVD Kabsch (`ref_svd`) on the canonical input, computed outside timing.

| Track | Distance requirement | Ranking requirement | Basis |
|---|---|---|---|
| **F64** (strict, secondary) | every pair \|d − d_ref\| ≤ 1e-6 Å | equivalent with τ = 2e-6 Å (BQ-1 §1.2) | as BQ-1–BQ-3 |
| **F32** (common, **primary**) | every pair \|d − d_ref\| ≤ 1e-3 Å | equivalent with τ = 2e-3 Å | one tenth of input coordinate quantization. AdK XTC precision is 1000 (0.001 nm = 0.01 Å). mdCATH `1ha8A00` coordinates are float32 Å on a 0.01 Å grid. Both were checked locally |

**F32 history (disclosed):** the first F32 proposal used an a-priori msd bound, 64·u·(G_q + G_j)/W. Adapter
development showed this bound does not model float32 polynomial-coefficient cancellation: MDTraj reached about 152×
u-scaled error while staying below 2.5e-4 Å. The bound was replaced by the input-resolution tolerance above, which
does not depend on any implementation's error. The msd ratio is kept as a descriptive native-precision measure.

**Gates, per (implementation, workload), all three queries, every window:**
- finite distances;
- the track distance requirement;
- zero rank-non-equivalent query-windows;
- query-path top-k identical to the top-k of the full distance arrays;
- input unchanged.

Strict identity and near-tie counts are reported.

**Degeneracy suite (SYNTHETIC, seed 20261021, deciding; development used seed 777).** Evaluated per track:
- query identical to database rows;
- collinear, planar and mirror-image queries;
- exact duplicates (bitwise-equal distances);
- database permutation;
- rigid transform;
- near ties.

Failures are **numerical-limitation findings**. They do not change comparison eligibility, and they never justify
integration by themselves. V2b must pass F64 degeneracy for INTEGRATION_SUPPORTED.

## §5. Workloads (frozen: `results/cmp1_workloads.json`)

**Query frames:** drawn with seed 20261020 from positions never used by EB-1, BQ-1, BQ-2, BQ-3 or the development
checks.

| ID | Data | n_D × windows | Queries |
|---|---|---|---|
| R1_D1 | AdK open (real, L = 214) | 1015 × 207 | 570, 580, 736 |
| R2_D1plus | AdK open + closed (real) | 2693 × 207 | 932, 941, 942 |
| R3_D2 | mdCATH 1ha8A00 (real, L = 51) | 405 × 44 | 273, 319, 388 |
| S1_P2slice | SYNTHETIC_tiled, L = 1000 | 10,000 × 100 | tiled 20, 21, 22 |
| S2_P2D | SYNTHETIC_tiled, L = 1000 (P2) | 10,000 × 993 | tiled 23, 24, 25 |

- **Real workloads:** reused development data, non-duplicated MD frames. Duplicate database rows and query membership
  are counted and reported.
- **SYNTHETIC rows:** performance-only; they repeat frames.
- **Out of scope:** held-out cohort domains. No real P2-scale non-duplicated data exists locally; this is stated as a
  limitation.

## §6. Measurements and order

**Process:** each (round, workload, implementation) runs in a fresh process, one thread (`OMP_NUM_THREADS = 1`,
BLAS = 1, MDTraj `parallel=False`). Each process runs, in order:
1. timed import/setup;
2. an identical SYNTHETIC warm-up;
3. timed canonical load;
4. timed preparation;
5. Q1, Q2 and Q3, each timed for wall, CPU and kernel (L1).

The first-query span is measured end to end from the start of preparation to the end of Q1. Outputs are released
between queries.

**Order:** 5 rounds. In each round, all 20 cells (5 workloads × 4 implementations) run in a permutation from seed
20261022.

**Gating:** only cells whose F32 gate passed are timed.

**Recorded per process:** threads (threadpool and OS), `/proc/loadavg`, MemAvailable and SwapFree.

**Order of steps:** tests → references → gates → degeneracy suites → timing.

**Budget:** ≤ 2.0 CPU-h for CMP-1 in total (tests, references, gates, degeneracy and timing). Stop before any step
once the total exceeds 1.8 CPU-h; unrun cells are reported as incomplete.

**Memory:** a process starts only if MemAvailable − predicted peak ≥ 600 MB, with a wait of at most 10 min;
otherwise it is skipped. `RLIMIT_AS` is 6 GB. No application is closed.

## §7. Decision rules (frozen; `cmp1_criteria.py`; failure paths tested)

**Eligibility.** For workload w and track T, an implementation is eligible iff:
- its gate passes in T;
- all 5 timing repetitions completed;
- every timed output (indices + distances) equals the gate output.

**Strongest eligible comparator E\*(w, T).** The **eligible external** implementation (E-MDT or E-MDA) with the lowest
median repeated-query latency (Q2–Q3 pooled over repetitions). Ties go to the lower incremental peak RSS, then the
name. E\* is undefined if no external implementation is eligible. It is defined per workload and track; it is never
a pooled "fastest anywhere".

**Outcome per (w, T), mutually exclusive categories:**
- **V2b_advantage** (performance): V2b eligible; median_V2b ≤ 0.67 × median_E\*; every per-round paired ratio
  E\*/V2b ≥ 1.25; V2b absolute peak ≤ 8 GiB; and V2b incremental peak ≤ 1.25 × E\*'s.
- **comparator_preferable**: E\* is at least as fast as V2b (median_E\* ≤ median_V2b).
- **no_decisive_difference**: both defined, neither of the above.
- **insufficient_evidence**: V2b ineligible, E\* undefined, or unpaired rounds.

**Competitor-specific numerical limitations:** failures in a track, or in degeneracy sub-cases, are listed separately.
They make that competitor ineligible in that track and workload. They **never** yield V2b_advantage, and they never
justify integration on their own.

**Conclusion.** The primary track is F32, and all five workloads are required:
- **INTEGRATION_SUPPORTED** iff V2b_advantage on every required workload and V2b passes the F64 degeneracy suite.
- **EXISTING_LIBRARY_PREFERABLE** iff comparator_preferable on every required workload. The E\* per workload is named
  and recommended for this task.
- **NO_DECISIVE_ADVANTAGE** otherwise, when every workload has a defined outcome.
- **INSUFFICIENT_EVIDENCE** if any required workload lacks evidence.

**Secondary track:** the F64 result is reported as "among float64-accurate implementations". It cannot produce
INTEGRATION_SUPPORTED. A float64 requirement for the project would need its own reviewed decision, with evidence
that float32-class error matters for the application.

**Descriptive, not deciding:** kernel (L1) ratios, first-query spans, setup and preparation times.

## §8. Reporting

The report covers:
- accuracy per track, workload and implementation, with representation error;
- limitations;
- setup, preparation, first-query and repeated-query latency (medians and ranges);
- L1 versus L2;
- memory (absolute and incremental);
- per-workload E\* and outcome;
- coverage, including incomplete cells;
- the conclusion and a recommendation.

Results are not committed or pushed. No follow-on experiment starts automatically.
