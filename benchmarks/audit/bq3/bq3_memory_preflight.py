#!/usr/bin/env python3
"""BQ-3 (direct-P2 qualification of unchanged V2b) memory preflight. Read-only: no data loaded, no timing.

Estimates the peak memory of each process the direct P2 run needs (from the array shapes used by the existing code:
bq2_run.load / b2opt_lib.prepare / search / ref_svd) plus the BQ-2 headroom (600 MB), compares the largest with
MemAvailable, and records swap state. Writes results/bq3_memory_preflight.json.
Run: ../wpd_local_retrieval/build/venv/bin/python -B bq3_memory_preflight.py
"""
import datetime
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'bq2'))
import bq2_lib as Q  # noqa: E402  (BASE_RSS, TRANSIENT, HEADROOM constants; unchanged)

N_D, L, W, N_WIN, N_Q, K = 10_000, 1000, 8, 993, 3, 10
F8 = 8


def meminfo():
    d = {}
    for line in open('/proc/meminfo'):
        k, v = line.split(':')
        d[k] = int(v.split()[0]) * 1024
    return d


def main():
    resident = (N_Q + N_D) * L * 3 * F8                   # query + database rows, float64
    source = 3706 * 214 * 3 * F8                          # D1+ array used to build tiled rows (freed before preparation)
    prepared = N_D * N_WIN * (W * 3 + 1) * F8             # centred windows + G_B
    retained = N_WIN * K * F8
    dist_one_query = N_WIN * N_D * F8                     # (n_win, n_D) float64 distance array
    fixed = Q.BASE_RSS + Q.TRANSIENT
    proc = {
        'timing_V2b': fixed + source + resident + prepared + retained,
        'gate_V2b': fixed + source + resident + prepared + retained + 2 * dist_one_query,
        'reference': fixed + source + resident + N_Q * dist_one_query,
        'timing_V0_baseline': fixed + source + resident + retained,
    }
    need = {k: v + Q.HEADROOM for k, v in proc.items()}
    mi = meminfo()
    worst = max(need, key=need.get)
    out = {'recorded_utc': datetime.datetime.now(datetime.UTC).isoformat(),
           'workload': {'n_D': N_D, 'L': L, 'W': W, 'n_win': N_WIN, 'queries': N_Q,
                        'label': 'SYNTHETIC_tiled, performance-only'},
           'components_bytes': {'baseline_interpreter_numpy': Q.BASE_RSS, 'transient_allowance': Q.TRANSIENT,
                                'tiling_source_D1p': source, 'resident_rows': resident, 'prepared_arrays': prepared,
                                'retained_topk': retained, 'distance_array_per_query': dist_one_query,
                                'headroom': Q.HEADROOM},
           'process_peak_estimate_bytes': proc, 'required_mem_available_bytes': need,
           'largest_requirement': {'process': worst, 'bytes': need[worst]},
           'processes_run_sequentially': True,
           'meminfo_bytes': {k: mi[k] for k in ('MemTotal', 'MemAvailable', 'SwapTotal', 'SwapFree')},
           'loadavg': open('/proc/loadavg').read().split()[:3]}
    out['sufficient'] = mi['MemAvailable'] >= need[worst]
    out['shortfall_bytes'] = max(0, need[worst] - mi['MemAvailable'])
    (HERE / 'results' / 'bq3_memory_preflight.json').write_text(json.dumps(out, indent=1) + '\n')
    for k, v in need.items():
        print(f'{k:20s} peak {proc[k] / 1e9:.3f} GB, required available {v / 1e9:.3f} GB')
    print(f"MemAvailable {mi['MemAvailable'] / 1e9:.3f} GB; SwapFree {mi['SwapFree'] / 1e6:.1f} MB; sufficient: "
          f"{out['sufficient']}; shortfall {out['shortfall_bytes'] / 1e9:.3f} GB")


if __name__ == '__main__':
    main()
