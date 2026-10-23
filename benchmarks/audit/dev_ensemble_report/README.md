# Descriptive ensemble comparison: demonstration on development data

Demonstration of the reusable command `python -m topofold_inputs.ensemble`. It takes named groups of trajectories and
produces residue-mapped tables, distributions, replica variability and temporal diagnostics.
- **Inputs:** two existing development data sets, run through unchanged components: the MDAnalysis adapter or the
  mdCATH bridge, then core κ/τ.
- **Evidence level:** descriptive only. No p-values, significance labels or uncertainty intervals. Development data
  (both data sets were used before; see `docs/DECISION_LOG.md`). No biological claim.

| Spec | Groups | Outputs |
|---|---|---|
| `spec_1ha8A00.json` | mdCATH `1ha8A00`, 320 K vs 348 K; 5 replicas each; source frames 43–429 | `results/1ha8A00/` |
| `spec_adk.json` | AdK open-start vs closed-start; one run each; source frames 0–1676, every 4th | `results/adk/` |

**Versions.**
- `results/` holds the v1 outputs (tool 0.1.0), kept as history.
- `results_v2/` holds the corrected outputs (tool 0.2.0): an origin-invariant circular temporal diagnostic, resultant
  lengths with a numerical conditioning rule, and explicit reasons for unavailable values. See
  `CORRECTION_NOTE_2026-10-23.md`, with verification in `results_v2/verification.json` (`verify_v2.py`).
- **Pooled distributions** weight members by their retained frame counts.

Each output directory holds `tables/*.csv` (distribution, replicas, between, temporal, exclusions), `report.md`
and `manifest.json`. `distributions.npz` and the per-member intermediate files are regenerated in the git-ignored
`build/`.

```bash
PYTHONPATH=python .venv/bin/python -m topofold_inputs.ensemble benchmarks/audit/dev_ensemble_report/spec_1ha8A00.json \
    --out OUTDIR --adapter-python benchmarks/audit/wpd_local_retrieval/build/venv/bin/python
```

## Checks

- **Reproducibility:** each spec was run twice. All tables, `distributions.npz` and `report.md` were byte-identical.
  `manifest.json` differs only in its creation time.
- **Mapping:** pooled JS distances agree with the values independently recorded by DEV-REP-1 (194 channel–window
  evaluations) and DEV-ENS-1 (846), within CSV rounding (max |Δ| 5 × 10⁻⁷).
- **Exclusions:** none on these data. Exclusion of windows touching gaps, missing Cα and breaks is tested on
  SYNTHETIC fixtures (`tests/test_ensemble_report.py`).
- **Resources per run:** about 30 CPU-s and 110 MiB (1ha8A00); about 19 CPU-s and 118 MiB (AdK).

## What the tables show (descriptive)

**mdCATH 1ha8A00, 320 K vs 348 K.**
- The largest pooled JS distances (0.21–0.35) are concentrated in residues 22–37. For example: τ 32–35 is 0.350,
  κ 22–24 is 0.259.
- In most of these windows the median distance between replicas of different temperatures is comparable to the
  median distance between replicas of the same temperature. For example, τ 32–35 is 0.614 between, against 0.564
  (320 K) and 0.393 (348 K) within.
- In κ 22–24, the 348 K replicas differ more among themselves (0.361) than the 320 K replicas (0.165).
- Median lag-1 autocorrelations are slightly higher at 320 K than at 348 K in every channel. For τ the v2 values are
  0.699 vs 0.633; the v1 values (0.646 vs 0.570) came from the superseded circular diagnostic.
- These are descriptions. Which windows differ, and why, is not established.

**AdK open-start vs closed-start.**
- The largest pooled JS distances occur at residues 119–124, 134–137, 142–146, 40–45 and 9–13 (up to 0.953 for
  τ 121–124), the same windows DEV-ENS-1 reported.
- There is one run per group, so replica variability is not available (n/a). Starting structure and sampling are
  confounded with any difference.
