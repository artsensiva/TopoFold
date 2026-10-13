# WP-B G0–G1 Feasibility Report: AdK (third-party MD, not SYNTHETIC)

Data: Orädd et al., *Sci. Adv.* 2021, Zenodo 10.5281/zenodo.5583119 (CC BY 4.0). References: RCSB 4AKE (open) and 1AKE (closed). The report covers data integrity (G0) and an independent, reference-based assessment (G1). **No TopoFold or candidate-feature scoring, retrieval, localization or threshold tuning was performed.** AdK is a development dataset, not held-out validation.

**Specification history.** Original freeze 2026-10-04T16:24:29.834471+00:00; Amendment 1 2026-10-04T16:36:08.777573+00:00; Amendment 2 2026-10-04T16:38:03.704000+00:00. All three were recorded before G0 or G1 ran and before any trajectory-derived structural quantity was computed. Amendment 1 answered review findings; Amendment 2 replaced a core algorithm that failed on a synthetic test (see `AMENDMENTS.md`). The original freeze record is preserved unchanged.

## G0: integrity and mapping

Result: **PASS** (finished 2026-10-04T16:53:52.752846+00:00). Stops: none.

| File | Size (B) | md5 and size match the source record | SHA-256 |
|---|---:|---|---|
| ORADD_open_state.gro | 6567118 | True | `25c90e78f549d97ed61814950efcdf9b859a0320734bd15bd09003677608a40f` |
| ORADD_closed_state.gro | 6568567 | True | `4d05b0d3767e1f5e948b2c038a42e9b1a2132e3085c5b50899f4bf797e9f16e4` |
| ORADD_open_state_traj.xtc | 780301896 | True | `39e77b58ded7279d27b9e331b2ff9e6d39a803183689612856df0aa800a40019` |
| ORADD_closed_state_traj.xtc | 580049352 | True | `f950ff0a7beadab0992d9c1cffaea540675b62d06ed9d1eadcc988f20f6f3cc9` |

**4AKE** (SHA-256 `ff798ee8791878eb58bac1c6bed32042f51b455512ac93b50bda8fa0ff0e7f78`; REVDAT dates ['28-FEB-24', '24-FEB-09']; hetero residues ['HOH']):
- chain A (primary): mapping valid **True**; missing none; mismatched none; residues with altlocs 0
- chain B (sensitivity): mapping valid **True**; missing none; mismatched none; residues with altlocs 0

**1AKE** (SHA-256 `651e952f55f1317f50f4e82b7f0e99053397032cd8da5b2896f79d0af9f619a9`; REVDAT dates ['25-DEC-24', '05-JUN-24']; hetero residues ['AP5', 'HOH']):
- chain A (primary): mapping valid **True**; missing none; mismatched none; residues with altlocs 0
- chain B (sensitivity): mapping valid **True**; missing none; mismatched none; residues with altlocs 0

| Trajectory | atoms gro / xtc | Cα | frames | nominal Δt (ps) | cadence irregularities | first–last (ps) | file duration (ns) | stated (ns) | box angle dev (°) | quantization (Å) |
|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|
| open_start | 95175 / 95175 | 214 | 2253 | 200.0 | 0 | 0–450400 | 450.4 | 450 | 0.0e+00 | 0.010 |
| closed_start | 95196 / 95196 | 214 | 1678 | 200.0 | 0 | 0–335400 | 335.4 | 335 | 0.0e+00 | 0.010 |

Verified facts are those measured from the files above. The stated durations come from the paper and are compared with the files, not assumed. Discrepancies and notes: open_start: 2253 frames vs 2251 expected from the stated 450 ns (recorded discrepancy, not a stop); closed_start: 1678 frames vs 1676 expected from the stated 335 ns (recorded discrepancy, not a stop). Coordinates are stored at the quantization shown; float64 arithmetic does not recover precision lost in the file.

## G1: independent reference-based assessment

**Binding to G0:** all checks passed (spec, raw-file and reference hashes match the G0 record; G0 ran after the Amendment-2 freeze).

**Alignment cores** (crystal Cα only; two-phase algorithm of Amendment 2):

| Core | size | converged | valid | reason if invalid |
|---|---:|---|---|---|
| A_1.0 | 60 | True | True |  |
| A_0.5 | 23 | True | False | size 23 < 30 |
| A_2.0 | 112 | True | True |  |
| B_1.0 | 68 | True | True |  |

Variants used: ['A_core1.0_f0.2', 'A_core1.0_f0.25', 'A_core1.0_f0.33', 'A_core2.0_f0.25', 'B_core1.0_f0.25']. Excluded (invalid cores): {'A_0.5': 'size 23 < 30'}.
d_ref (Å): A/A_1.0 8.16, A/A_2.0 8.24, B/B_1.0 8.07. Non-core segments (≥ 5 residues): [[10, 19], [30, 79], [88, 92], [114, 170], [172, 178], [200, 208], [210, 214]].

**Kabsch vs MDAnalysis (QCP):** 161 comparisons; max |diff| 2.0e-13 Å; predeclared tolerance 1e-06 Å → **PASS**.

### open_start (2253 frames)

Preprocessing: frames with the raw Cα chain split across the periodic box: 0. Sequential minimum-image reconstruction was applied to every frame; max consecutive Cα distance afterwards 4.11 Å. Raw files are unchanged.

| Observable (Å, primary core) | min | q05 | median | q95 | max | first 10 % mean | last 10 % mean | τ_int (frames, est.) | ESS (est.) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| rmsd_open | 0.40 | 2.30 | 4.21 | 10.98 | 14.68 | 2.58 | 11.09 | 307.3 | 7.3 |
| rmsd_closed | 4.72 | 6.47 | 8.65 | 12.72 | 16.68 | 7.89 | 12.39 | 231.6 | 9.7 |
| delta | -1.53 | -0.05 | 4.02 | 6.43 | 7.71 | 5.32 | 1.30 | 258.8 | 8.7 |
| core_rmsd_open | 0.29 | 0.76 | 0.98 | 1.29 | 1.56 | 0.85 | 1.03 | 204.2 | 11.0 |
| core_rmsd_closed | 0.64 | 0.93 | 1.15 | 1.42 | 1.70 | 1.05 | 1.09 | 79.5 | 28.3 |

### closed_start (1678 frames)

Preprocessing: frames with the raw Cα chain split across the periodic box: 0. Sequential minimum-image reconstruction was applied to every frame; max consecutive Cα distance afterwards 4.09 Å. Raw files are unchanged.

| Observable (Å, primary core) | min | q05 | median | q95 | max | first 10 % mean | last 10 % mean | τ_int (frames, est.) | ESS (est.) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| rmsd_open | 1.92 | 2.84 | 4.73 | 6.85 | 8.19 | 6.32 | 3.68 | 161.4 | 10.4 |
| rmsd_closed | 0.43 | 3.35 | 5.26 | 7.10 | 9.81 | 3.85 | 6.05 | 105.5 | 15.9 |
| delta | -7.75 | -3.27 | 0.45 | 4.20 | 7.08 | -2.47 | 2.36 | 137.7 | 12.2 |
| core_rmsd_open | 0.61 | 0.78 | 0.95 | 1.18 | 1.35 | 0.96 | 0.92 | 52.0 | 32.3 |
| core_rmsd_closed | 0.35 | 0.88 | 1.08 | 1.31 | 1.47 | 1.02 | 1.11 | 35.9 | 46.7 |

### Non-core segment RMSDs (pre-registered observable; from the per-frame CSVs)

Mean segment Cα RMSD (Å) after primary-core superposition, first vs last 10 % of frames.

**open_start**

| Segment | to open: first → last | to closed: first → last |
|---|---|---|
| 10–19 | 1.53 → 2.84 | 2.02 → 2.96 |
| 30–79 | 2.25 → 2.79 | 7.86 → 7.28 |
| 88–92 | 1.29 → 0.81 | 2.20 → 1.63 |
| 114–170 | 4.19 → 21.22 | 13.13 → 22.88 |
| 172–178 | 1.17 → 1.98 | 2.25 → 2.55 |
| 200–208 | 1.04 → 0.97 | 2.33 → 2.05 |
| 210–214 | 1.24 → 1.10 | 1.79 → 1.56 |

**closed_start**

| Segment | to open: first → last | to closed: first → last |
|---|---|---|
| 10–19 | 3.39 → 2.24 | 3.78 → 2.70 |
| 30–79 | 7.71 → 3.16 | 2.36 → 7.00 |
| 88–92 | 0.74 → 1.25 | 1.58 → 2.10 |
| 114–170 | 9.59 → 6.18 | 6.58 → 9.36 |
| 172–178 | 1.09 → 2.02 | 2.25 → 2.13 |
| 200–208 | 1.44 → 0.93 | 2.38 → 2.36 |
| 210–214 | 1.42 → 1.06 | 2.20 → 2.01 |

### Descriptive proximity classes by variant

A sustained run has an elapsed span ≥ 5.0 ns. Each cell gives frames / sustained runs / brief runs / longest span (ns). The classes are provisional and descriptive, not validated state definitions.

| Variant (radius Å) | category | trajectory | starts near own ref | near open | near closed | ambiguous |
|---|---|---|---|---|---|---|
| A_core1.0_f0.2 (1.63) | unresolved | open_start | False | 6 / 0 / 6 / 0.0 | 0 / 0 / 0 / 0.0 | 2247 / 2 / 4 / 391.8 |
| A_core1.0_f0.2 (1.63) | unresolved | closed_start | False | 0 / 0 / 0 / 0.0 | 2 / 0 / 1 / 0.2 | 1676 / 1 / 0 / 335.0 |
| A_core1.0_f0.25 (2.04) | unresolved | open_start | False | 48 / 0 / 35 / 0.6 | 0 / 0 / 0 / 0.0 | 2205 / 3 / 32 / 369.2 |
| A_core1.0_f0.25 (2.04) | unresolved | closed_start | False | 4 / 0 / 3 / 0.2 | 7 / 0 / 3 / 0.8 | 1667 / 5 / 1 / 155.0 |
| A_core1.0_f0.33 (2.69) | unresolved | open_start | False | 268 / 1 / 103 / 6.0 | 0 / 0 / 0 / 0.0 | 1985 / 9 / 95 / 240.4 |
| A_core1.0_f0.33 (2.69) | unresolved | closed_start | False | 57 / 0 / 34 / 1.2 | 33 / 0 / 13 / 1.4 | 1588 / 9 / 38 / 123.2 |
| A_core2.0_f0.25 (2.06) | unresolved | open_start | False | 56 / 0 / 38 / 0.8 | 0 / 0 / 0 / 0.0 | 2197 / 3 / 35 / 369.2 |
| A_core2.0_f0.25 (2.06) | unresolved | closed_start | False | 5 / 0 / 4 / 0.2 | 3 / 0 / 1 / 0.4 | 1670 / 4 / 1 / 155.0 |
| B_core1.0_f0.25 (2.02) | unresolved | open_start | False | 27 / 0 / 19 / 0.4 | 0 / 0 / 0 / 0.0 | 2226 / 3 / 16 / 391.8 |
| B_core1.0_f0.25 (2.02) | unresolved | closed_start | False | 2 / 0 / 2 / 0.0 | 2 / 0 / 1 / 0.2 | 1674 / 3 / 0 / 155.0 |

**Observed cross visits vs claimed transition evidence:**
- open_start: sustained visit near the other reference observed in 0 of 5 variants; transition evidence claimed (all variants, starting near own reference): **False**.
- closed_start: sustained visit near the other reference observed in 0 of 5 variants; transition evidence claimed (all variants, starting near own reference): **False**.

## Outcome

- Categories by variant: {'A_core1.0_f0.2': 'unresolved', 'A_core1.0_f0.25': 'unresolved', 'A_core1.0_f0.33': 'unresolved', 'A_core2.0_f0.25': 'unresolved', 'B_core1.0_f0.25': 'unresolved'}; variants agree: True.
- **Outcome: unresolved pending specified evidence.**

**Evidence that would resolve the outcome** (any revision needs a new freeze before re-analysis; the current criteria are not changed after inspection):
1. State criteria fixed independently of these trajectories, e.g. segment-level criteria derived from the crystal structures, or verified literature domain observables.
2. An independent scale for thermal fluctuation around each reference state, to justify a proximity radius; the current 0.25·d_ref is a pre-declared choice, not a calibrated one.
3. Replicate trajectories per starting state.
4. Clarification of the closed-start system: it was simulated without the Ap5A ligand present in 1AKE, so its stability near the closed reference is not expected a priori.

**Interpretation limits.**
- Proximity to a reference by RMSD is not proof of a conformational transition.
- With one trajectory per starting state and autocorrelated frames, no equilibrium populations, rates or free energies can be inferred. The τ_int and ESS values are exploratory diagnostics, not evidence of equilibrium sampling; the observables drift over the runs.
- Segment observations are descriptive changes in pre-registered observables, not claimed transitions or assigned conformational states.
- The synthetic core-algorithm test (Amendment 2) is evidence about that construction only, not proof about AdK. On AdK the primary core converged with the size shown above.
- Pooling the two trajectories mixes starting states by construction and cannot establish spontaneous kinetics.
- If no transition occurred, that is a limitation of this dataset for the proposed task, not a failure of any representation; none was evaluated here.
- G2 label definitions and the held-out selection remain prerequisites for any scoring protocol.
