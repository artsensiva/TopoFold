# BQ-3 Direct-P2 Qualification of Unchanged V2b: Memory Preflight (2026-10-13), Stopped

**Status:** stopped before any protocol freeze, test or measurement, as instructed for insufficient memory. Nothing
was run except a read-only memory preflight. No application was closed, no system setting was changed, and the
workload was not reduced.

**Preserved:** BQ-2 (conclusion C), BQ-1, EB-1 and all earlier experiments are unchanged. Their freeze hashes were
re-verified before the preflight.

## Requirement

The estimate comes from `benchmarks/audit/bq3/bq3_memory_preflight.py`, with values in
`benchmarks/audit/bq3/results/bq3_memory_preflight.json`.
- **Workload:** unchanged P2: n_D = 10,000, L = 1,000, W = 8, 993 windows, 3 queries. SYNTHETIC_tiled,
  performance-only.
- **Method:** process peaks are summed from the array shapes the existing code allocates:
  - fixed allowances: interpreter/NumPy baseline and per-window transient;
  - the D1+ array used to build the tiled rows;
  - resident rows;
  - prepared arrays (≈ 1.99 GB);
  - retained top-k outputs;
  - per-query distance arrays (gate and reference processes).
- **Headroom:** the BQ-2 headroom of 600 MB is added.
- **Sequencing:** processes run one at a time, so the requirement is the largest single process.

| Process | Required MemAvailable |
|---|---|
| Gate (V2b with distances and reference slice) — largest | ≈ 3.16 GB |
| Timing (V2b) | ≈ 3.00 GB |
| Reference | ≈ 1.26 GB |
| Baseline V0 timing | ≈ 1.02 GB |

**Available at the check:** MemAvailable ≈ 2.03 GB, with swap exhausted (SwapFree ≈ 0). The shortfall is ≈ 1.13 GB.
The direct P2 qualification therefore cannot run safely on this machine in its current state.

**To run it, one of the following is needed:**
- ≥ 3.2 GB MemAvailable at start, with this machine otherwise in a comparable state; or
- a separate host with that headroom.

A protocol and acceptance rule must then be frozen before execution, as specified in the task. That protocol was not
written here.

## Clarification (erratum, recorded without editing frozen files)

BQ-1's preparation model (p0 + p1·n_D·n_win) omitted a per-window term.
- The BQ-1 fit grid, n_D ∈ {250, 500, 1015} × n_win ∈ {20, 207}, varied n_win at fixed n_D. It **would have made that
  term identifiable**; the omission was a model-specification choice, not a design limitation.
- The sentence in `docs/B2_SCALING_CONFIRMATION_PROTOCOL.md` §3, "so the per-window term could not be separated", is
  therefore inaccurate.
- This does not change BQ-1's verdict, BQ-2's conclusion C, or either frozen document.

## Clarification (2026-10-13)

"Read-only memory preflight" above means that no input data was changed and no timing was performed. The preflight
inspected `/proc/meminfo` and `/proc/loadavg` and **wrote** one record file,
`benchmarks/audit/bq3/results/bq3_memory_preflight.json`. Writing this note and that record was file creation, not
system inspection.
