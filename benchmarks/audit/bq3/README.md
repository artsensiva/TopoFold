# BQ-3: direct-P2 qualification of unchanged V2b (stopped at preflight)

The only artifact here is a memory preflight. It loads no data and does no timing. It reads `/proc/meminfo` and
`/proc/loadavg` and **writes** `results/bq3_memory_preflight.json`.

The preflight showed that the run needs more available memory than the machine had, so no protocol was frozen and no
measurement was made. See [`docs/B2_DIRECT_P2_QUALIFICATION_PREFLIGHT.md`](../../../docs/B2_DIRECT_P2_QUALIFICATION_PREFLIGHT.md).

Run: `../wpd_local_retrieval/build/venv/bin/python -B bq3_memory_preflight.py`. Running it again overwrites the
recorded JSON with the current machine state.
