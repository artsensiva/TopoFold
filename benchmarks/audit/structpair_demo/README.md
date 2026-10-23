# STRUCTPAIR-VERIFY-1: verification of the two-structure command on the KRAS inputs

`python -m topofold_inputs.structpair` (v0.1.0) is a thin command over existing components: the Cα adapter, κ/τ from
the core, `compare.window_distances`, MDAnalysis QCP RMSD and `ensemble.exclusion_causes`. This directory checks it
against the frozen KRAS-STRUCT-1 and KRAS-STRUCT-2 results. It also checks that several valid input specifications run
without code edits, and that invalid ones are refused.

```bash
.venv/bin/python benchmarks/audit/structpair_demo/fetch_inputs.py   # inputs are git-ignored; verifies SHA-256
.venv/bin/python benchmarks/audit/structpair_demo/verify.py          # fresh build/<UTC stamp>/ per run
```

## Specs

All in `specs/`:
- `kras1_9O0N_vs_9O0O.json`: same-numbering ranges with an altloc sensitivity policy.
- `kras2_5US4A_vs_5USJA.json`: ranges with a gap.
- `kras2_5US4A_vs_5USJA_explicit_exclude.json`: explicit pairs with `on_missing: exclude`.
- `demo_9O0NA_vs_5US4B_mismatch_flag.json`: an input-handling demonstration (G12 vs D12 flagged). Its values are not
  interpreted.
- `verify.py` also generates the 12 chain-pairing specs and three invalid specs in `build/`.

## Results (`results/verification.json`; compact outputs in `results/<run>/`)

| Check | Result |
|---|---|
| C1 KRAS-STRUCT-1 reproduction | 993 window values vs the frozen table: 0 mismatches, max \|difference\| 0. Global RMSD equal to the record (0.718 Å); numpy Kabsch differs by 1.1 × 10⁻¹³ Å. Spearman agreement and altloc-sensitivity maxima equal to the record |
| C2 KRAS-STRUCT-2 primary (5US4 A vs 5USJ A) | 918 values: 0 mismatches. Global RMSD equal. The same 5 excluded windows (numbering discontinuity 59→68). Altloc sensitivity "not applicable" (no altlocs) |
| C3 all 12 chain pairings | 12 of 12 runs exit 0. Global RMSD max \|difference\| 0. 2,754 recorded min/median/max aggregates: 0 mismatches |
| C4 explicit pairs, `on_missing: exclude` | Missing residues named exactly (5US4 A: 60–64; 5USJ A: 61–63; causes `missing_coordinates`, `missing_CA`). 903 values away from the gap identical to C2 |
| C5 residue-name mismatch | `flag`: A:12 GLY vs ASP flagged, with windows 11–13 marked. `reject`: exit 2 with the residue named |
| C6–C8 refusals | Non-empty output directory, missing chain and absent chain: each exits 2 with a message |

Unit tests: `tests/test_structpair.py`, 11 tests. They cover the spec rules, refusals, a Kabsch reference, equal
coverage, circular τ wrapping, and a SYNTHETIC end-to-end run with a missing residue, a residue-name mismatch and
byte-identical reruns.

Resources: 182 CPU-s for `verify.py` (16 command runs).
