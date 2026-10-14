> **Superseded 2026-10-13** by the frozen [`CMP1_COMPETITIVE_EVALUATION_PROTOCOL.md`](CMP1_COMPETITIVE_EVALUATION_PROTOCOL.md), which corrects this draft. The draft is kept as reviewed.

# CMP-1 Competitive Evaluation of Window Cα RMSD Retrieval: Protocol DRAFT (for review, not frozen)

**Status:** draft, 2026-10-13. Nothing is implemented, measured, installed or downloaded. No held-out data is used.
Production code is unchanged.

**Basis:**
- BQ-3 (`docs/B2_DIRECT_P2_QUALIFICATION_CLOSEOUT_2026-10-13.md`): the V2b harness qualified for prepared-database
  operation on SYNTHETIC P2-dimension rows; first-query operation did not meet 10 s.
- BQ-2 remains conclusion C.
- Neither integration option (`docs/B2_V2B_INTEGRATION_PROPOSAL.md`) proceeds before this evaluation.

**Question:** on an identical task, with equal information and equal opportunity to prepare the database, does the
V2b harness offer a predeclared advantage over strong existing implementations that justifies integrating it, rather
than using an existing library?

## §1. Task and outputs (identical for every implementation)

**Inputs:**
- one query frame and n_D database frames, given as Cα coordinates (Å) of the same residues in the same order;
- the window length W = 8;
- an explicit list of window start positions;
- k = 10.

**Distance:** for every window and database frame, the minimal RMSD over **proper rotations** (det = +1) after
removing each window's centroid. Unweighted, over the W Cα atoms. Reflections are not allowed.

**Search:** exhaustive; every database frame and every listed window, with no pruning or approximation.

**Outputs, per window:**
- (a) the distance vector (n_D values);
- (b) the top-k database positions under the canonical order: ascending (distance, database position).

Every implementation's distances pass through the **same** top-k routine (`topk_lex`), and its cost is included for
everyone.

**Equal input values:** all implementations receive the same float32-representable coordinate values.
- Real workloads use the S2 rounded arrays already used by EB-1/BQ.
- The SYNTHETIC tiled rows are built by float64 translation and so are **not** float32-representable. They would be
  cast once to float32 and back before being given to every implementation, a documented change from BQ-2/BQ-3.
- Values are passed in Å even to libraries whose native unit is nm. The RMSD is unit-agnostic, and converting to nm
  would alter the values by rounding. The adapter must document this.
- Internal arithmetic precision may differ, so each implementation is checked against its declared precision class
  (§5).

## §2. Implementations (verified primary sources)

| ID | Implementation | Role | Primary source (verified 2026-10-13) |
|---|---|---|---|
| T-V2b | V2b harness (`benchmarks/audit/b2opt/b2opt_lib.py`), unchanged | candidate | Horn 1987, JOSA A 4:629; Theobald 2005, Acta Cryst A 61:478, doi:10.1107/S0108767305015266; Liu, Agrafiotis & Theobald 2010, J Comput Chem 31:1561, doi:10.1002/jcc.21439 |
| T-V0 | EB-1 harness kernel (Horn key matrix + `eigvalsh`) | internal baseline | as above |
| E-MDT | MDTraj `md.rmsd(target, reference, frame, atom_indices, parallel, precentered)` | strong external, complete workflow and kernel | McGibbon et al. 2015, Biophys J 109:1528, doi:10.1016/j.bpj.2015.08.015. The paper reports a QCP engine about 3× faster than the original implementation, citing Haque, Beauchamp & Pande 2014, bioRxiv doi:10.1101/008631. Signature, precentering and OpenMP behaviour are from the MDTraj API documentation |
| E-MDA | MDAnalysis: `lib.qcprot.CalcRMSDRotationalMatrix` (kernel) and `analysis.rms.RMSD` / `rms.rmsd(..., superposition=True)` (workflow) | external, widely used | Michaud-Agrawal et al. 2011, J Comput Chem 32:2319, doi:10.1002/jcc.21787; Gowers et al. 2016, Proc. SciPy, doi:10.25080/Majora-629e541a-00e. The documentation states QCP superposition (Theobald). MDAnalysis 2.10.0 is already installed |
| E-CPP | CPPTRAJ (`rms` / `rms2d`) | strong external C++ workflow; **conditional** | Roe & Cheatham 2013, J Chem Theory Comput 9:3084, doi:10.1021/ct400341p. Needs an AmberTools installation (a download) |
| R-SVD | `b2opt_lib.ref_svd` (batched Kabsch SVD, float64) | correctness reference only, not timed competitively | — |

**Still to verify at the adapter review, not asserted here:**
- MDTraj's internal floating-point precision and whether `precentered=True` is valid together with `atom_indices`
  subsets. Neither is stated in the primary sources read.
- CPPTRAJ's exact command semantics for best-fit RMSD over atom masks.

**Not proposed:**
- Biopython QCPSuperimposer and `scipy.spatial.transform.Rotation.align_vectors`: per-pair Python loops, not strong
  baselines at this scale. SciPy may serve as an extra correctness spot check.
- Structure-search tools with a different task (e.g. Foldseek).

## §3. Two comparison levels (reported separately)

**L1, solver kernel.** The function from (one centred query window, n_D database windows) to n_D RMSD values.
- Each implementation receives inputs already arranged in its preferred prepared form.
- Excluded: data loading, per-window preparation and top-k.
- Measures the arithmetic engine only.

**L2, complete retrieval workflow.** The full single-query search over all listed windows, producing the §1 outputs.
- Includes per-window selection and slicing, the implementation's own centring, distance computation and the shared
  top-k.
- Database preparation is a separate, measured step (§4).
- This is the level relevant to integration.

## §4. Equal preparation opportunity and measurements

Every implementation may prepare the database once, in its documented best way, and the preparation is timed and its
memory counted:
- **T-V2b:** `prepare()`, centred windows plus G_B.
- **E-MDT:** per window, `atom_slice` + `center_coordinates()`, then `precentered=True` (if valid; §2).
- **E-MDA:** an in-memory Universe (MemoryReader) with per-window AtomGroups. Optional pre-centred copies are allowed
  if documented.
- **E-CPP:** loading and any per-mask preparation its documentation supports.

No implementation receives information about the query, the expected answer or near-tie structure.

**Measured per (implementation, workload)**, in fresh processes, ≥ 5 repetitions, interleaved seeded order:
- **setup:** imports and library initialization, separately;
- **database preparation:** wall and CPU;
- **first-query latency:** preparation followed immediately by query 1, measured as one span. It excludes input
  loading, which is recorded separately and identically for all implementations;
- **repeated-query latency:** queries 2…m on the prepared database;
- **peak RSS:** lifetime peak per process, plus exact prepared-array bytes where obtainable;
- **context:** CPU/wall ratio, OS thread counts, `/proc/loadavg`, MemAvailable, SwapFree.

**Threads:**
- **Primary comparison:** one thread for everyone (`OMP_NUM_THREADS = 1`, MDTraj `parallel=False`, BLAS = 1).
- **Secondary, descriptive:** all implementations at the same thread count T, where supported. Implementations without
  threading are reported as single-threaded and not penalized or credited.

**Memory and CPU:**
- The BQ-2/BQ-3 memory start rule (predicted peak + 600 MB headroom).
- A proposed cap of ≤ 2 CPU-h in total, including tests and verification.

## §5. Correctness and ranking (checked separately; before any timing; mandatory)

**Precision classes, declared per implementation before measurement:**

| Class | Distance tolerance vs R-SVD | Basis |
|---|---|---|
| float64 | 10⁻⁶ Å | as BQ-1–BQ-3 |
| float32 | proposed 10⁻³ Å | float32 relative rounding near 6·10⁻⁸ on coordinates up to about 10² Å, amplified by the square root near RMSD ≈ 0. To be justified by a numerical test before freezing |

**Checks:**
1. **Distances:** all finite, and max |d − d_ref| within the class tolerance, on every workload, query and window.
2. **Ranking** (separate from distances): rank equivalence as in BQ-1 §1.2, with τ_rank = 2 × the class tolerance.
   Strict-identity and near-tie counts are reported. Zero non-equivalent query-windows are required.
3. **Degeneracies (SYNTHETIC):**
   - collinear, planar and mirror-image windows (proper rotation only);
   - identical structures (RMSD 0);
   - exact duplicates (bitwise-equal distances, position order);
   - database permutation;
   - rigid transform of the query;
   - constructed near ties.

   The frozen-V2 failure (collinear identical windows) shows that QCP-type solvers can hit a precision floor here.
   External QCP implementations get the same tests.
4. **Output identity:** each timed run's top-k equals the verified top-k of the same implementation.

Any failure excludes that implementation from the speed criteria and is reported as a correctness result.

## §6. Workloads

**Real, non-duplicated development workloads** (reused development data):
- D1: AdK open state, distinct frames.
- D1+: AdK open + closed, 2,693 database frames.
- D2: mdCATH `1ha8A00`, the EB-1 development domain.
- Configurations span n_D and n_win as in the BQ-2 fit grid, with query frames not used before (seeded draw).

These are **reused development evidence**. They support engineering comparison only. They are not independent
scientific validation, and no held-out cohort domain is used.

**SYNTHETIC stress workloads (labelled separately, performance-only):**
- the BQ-2/BQ-3 tiled rows at S1, S4_P2slice and P2D dimensions;
- the degeneracy fixtures.

Tiled rows repeat frames, so they test scale, not ranking difficulty or data realism.

**Missing evidence:** no real, non-duplicated workload at P2 scale (L ≈ 1,000, 10,000 frames) exists locally.
Obtaining one would need a download authorization and a recorded selection rule. Until then, P2-scale statements
rest on SYNTHETIC rows.

## §7. Predeclared integration criterion (proposed values for review; frozen before any measurement)

Integration of V2b is justified only if **all** of the following hold:

1. **Correctness:** V2b passes §5 in the float64 class on every workload and degeneracy test.
2. **Advantage over the strongest external L2 workflow** E\*. E\* is the fastest external implementation that passes
   §5 in its own class, single-thread. V2b must meet either (a) or (b):
   - **(a) speed:** on every real development workload and on the SYNTHETIC P2-dimension workload, V2b's median
     repeated-query latency is ≤ 0.67 × E\*'s (≥ 1.5× faster), and every per-round paired ratio is ≥ 1.25; **or**
   - **(b) correctness:** E\* fails §5 ranking or degeneracy checks on a real workload or a degeneracy test while V2b
     passes, and V2b's median repeated-query latency is ≤ 1.25 × E\*'s.
3. **Memory:** at P2 dimensions, V2b's peak RSS is ≤ 8 GiB (R4) and ≤ 1.25 × E\*'s peak RSS, counting prepared
   arrays for both.
4. **First query:** reported for all implementations. It is not part of the criterion, because BQ-3 already showed
   V2b does not meet 10 s; any first-query claim needs its own criterion.

**If the criterion is not met,** the recommendation is **not** to integrate V2b. Instead, consider depending on
E\* (licence and maintenance reviewed separately) or leaving the task outside the library.

**If it is met,** an integration option is proposed for review. A successful comparison would still establish no
biological capability, clinical relevance or novelty over existing RMSD methods.

## §8. Proposed experimental adapter (not implemented)

A small module, `benchmarks/audit/cmp1/adapters.py` (harness only), with one function per implementation:
- `prepare(X_db, windows) -> handle`;
- `query(handle, x_query) -> (distances (n_win, n_D), topk)`.

Its rules:
- It calls only documented public APIs.
- It converts no units and changes no values beyond the documented layout or dtype each library requires.
- It is reviewed against §2's open questions before any timing.

**Environment.** MDTraj, and AmberTools/CPPTRAJ, would need installation, which means downloads. They are **not
authorized** in this task, and the exact versions would be pinned in the frozen protocol. Without that
authorization, CMP-1 could proceed only with T-V2b, T-V0 and E-MDA, and the strongest-baseline requirement (AGENTS.md
rule 6) would then be **unmet**. That limitation would have to be stated in any result.

## §9. Not in scope

- Production integration.
- Held-out data.
- Biological or oncology claims.
- GPU implementations.
- Approximate or pruned search.
