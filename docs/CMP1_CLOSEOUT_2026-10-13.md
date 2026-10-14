# CMP-1 Closeout and Clarifications (2026-10-13)

**Closes:** CMP-1 (`docs/CMP1_COMPETITIVE_EVALUATION_PROTOCOL.md`; report `benchmarks/audit/cmp1/REPORT.md`; analysis
`benchmarks/audit/cmp1/results/cmp1_analysis.json`).

**Unchanged:** the frozen protocol, code, measurements, analysis and report. Their hashes were re-verified against
`cmp1_protocol_freeze.json` and `cmp1_pre_measurement_record.json` before this addendum. This document only clarifies
and qualifies them.

## 1. Verdict (unchanged)

- **Primary conclusion (common F32 track): NO_DECISIVE_ADVANTAGE. V2b integration is unsupported.**
- **Per-workload outcomes (frozen rules):**
  - V2b_advantage on R1_D1, S1_P2slice and S2_P2D, where the strongest eligible comparator was MDAnalysis;
  - comparator_preferable on R2_D1plus and R3_D2, where it was MDTraj.
- **Retained, without novelty claims:**
  - V2b computes the same exhaustive proper-rotation RMSD retrieval at float64 accuracy;
  - it passed every gate and degeneracy sub-case;
  - it was faster than the eligible comparator on the three workloads above.

  These are engineering properties of a harness implementation of an established method (Horn/QCP). They are not a
  new algorithm, a new capability or a biological result.

## 2. Clarification: conflicting degeneracy-eligibility wording

**The conflict.** The frozen protocol contains two incompatible sentences:
- §4: degeneracy failures "do not change comparison eligibility";
- §7: "failures in a track, **or in degeneracy sub-cases** … make that competitor ineligible in that track and
  workload".

**The interpretation implemented before measurement.** `cmp1_criteria.eligible` (hashed at 2026-10-13T14:53:22Z,
before the run started at 14:53:28Z) uses only three conditions:
- the workload accuracy gate in the track;
- completion of all timing repetitions;
- identity of the timed outputs with the gate outputs.

Degeneracy results never enter eligibility; they are reported as numerical limitations. They enter the decision in one
place only: V2b must pass the F64 degeneracy suite for INTEGRATION_SUPPORTED. The §4 sentence is therefore the one
that was implemented, and the §7 phrase "or in degeneracy sub-cases" does not describe the executed rule.

**The other reading would not change the integration verdict.** Both external comparators failed at least one F32
degeneracy sub-case: MDTraj the identical, planar and mirror queries; MDAnalysis the collinear query. Under the §7
reading, no external comparator would be eligible on any workload, E\* would be undefined everywhere, and the
conclusion would be INSUFFICIENT_EVIDENCE. Integration would remain unsupported. The frozen protocol and verdict are
not edited.

## 3. Qualifications

- **MDTraj recommendation, scoped.** MDTraj 1.11.2 met the common F32 requirement only on **R2_D1plus and R3_D2**,
  and was the faster implementation there. It did **not** meet it on R1_D1, S1_P2slice or S2_P2D (distance errors
  above 1e-3 Å).
  - On R1_D1 this was one pair per query out of 210,105, at about 2 Å reference RMSD, and outside every top-10. This
    comes from a post-hoc diagnostic.
  - Rankings were F32-equivalent everywhere.
  - The recommendation therefore applies to workloads where its accuracy has been checked against the requirement;
    it is not a blanket substitution.
- **Absolute vs incremental memory** (medians; `cmp1_analysis.json`):

  | Workload | Implementation | Absolute peak RSS | Incremental peak RSS |
  |---|---|---|---|
  | R2_D1plus | V2b | 169 MB | 129 MB |
  | R2_D1plus | MDTraj | 135 MB | 79 MB |
  | R3_D2 | V2b | 45 MB | 4.9 MB |
  | R3_D2 | MDTraj | **60 MB** | 3.1 MB |

  On R3, MDTraj's **absolute** footprint is higher than V2b's because of its library import. Its incremental
  (data-dependent) footprint is lower. The frozen memory criterion used incremental peaks, so absolute-memory
  statements must be made separately.
- **MDAnalysis timings include per-pair Python dispatch.** The adapter calls `qcprot.CalcRMSDRotationalMatrix` once per
  (window, database row) from Python, and its kernel time is close to its total. V2b's advantage over E-MDA therefore
  largely reflects call granularity, not a better RMSD algorithm. No vectorized float64 external implementation was
  evaluated.
- **The 0.001 Å tolerance is an engineering acceptance criterion.** It is one tenth of the input coordinate
  quantization (0.01 Å). It shows numerical agreement with a float64 reference. It is **not** evidence that 0.001 Å
  (or any distance accuracy) is adequate or necessary for a biological question. That would need its own
  task-specific justification.
- **Bounded scope:**
  - one machine, one thread, four implementations;
  - CPPTRAJ not evaluated;
  - development and SYNTHETIC data;
  - not a field-wide best-tool claim.

## 4. CPU accounting

| Activity | CPU | Instrumented? |
|---|---|---|
| Deciding run: tests, references, gates, degeneracy suites, 85 timed cells | 0.521 h (`cmp1_measurements.json` `cpu_h_total`) | yes (parent + children `getrusage`) |
| Adapter-development checks: three runs of `cmp1_devchecks.py` before the freeze | not recorded; each run took a few minutes of wall time | **no** |
| Post-hoc MDTraj error diagnostic | 1.3 s (`cmp1_posthoc_mdtraj_error_diagnostic.json`) | yes (in-process) |
| Environment installation (download, offline install) | not recorded | no |

The 2 CPU-h cap was applied by the runner to the deciding run only. The uninstrumented activities were not counted
against it; their wall-time bound is given above.

## 5. Status

- V2b remains an **audit artifact** (`benchmarks/audit/b2opt/`). The integration options in
  `docs/B2_V2B_INTEGRATION_PROPOSAL.md` are not pursued.
- Future tasks reuse established numerical libraries where their checked accuracy is adequate
  (`docs/BIOLOGICAL_VALIDATION_DECISION_MEMO_2026-10-13.md` §6).
