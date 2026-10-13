# LE-1: Legacy κ–τ engine evaluation (development tier, AdK)

**THIRD-PARTY MD DERIVATIVE** for AdK results (Orädd et al. 2021 AdK MD, Zenodo 10.5281/zenodo.5583119, CC BY 4.0); SYNTHETIC for fixtures and probe. Protocol `docs/LEGACY_ENGINE_EVALUATION_PROTOCOL.md` (SHA-256 `1ee5848025956f54…`, frozen 2026-10-05T13:40:40.814176+00:00); Amendment 1 (`docs/LEGACY_ENGINE_EVALUATION_AMENDMENT_1.md`, recorded 2026-10-05T16:45:48.628156+00:00, no LE-1 result file existed; attempt 2 in progress; no outcome inspected). Synthetic verification before the run: `16 passed in 0.23s` (exit 0). A **geometric local-conformation retrieval benchmark on one protein**: not biological, not held-out.

**Run history.** Attempt 1 ended at `600/834 windows 1964s` with no result file and no exit line; cause not established; coincided with the end of an assistant session; no termination record exists; outcomes observed: none (log contains progress lines only; results are written only at completion). Attempt 2 (identical code and protocol) produced the results below in 2282.0 s.

## 1. Legacy method: claim–evidence matrix

| Mechanism / claim | Status | Evidence |
|---|---|---|
| Raw Cartesian Fréchet (`discrete_frechet_distance`, Python) | Å, no superposition: not SE(3)-invariant (KNOWN_ISSUES M5); not part of the index | code `distance.rs::discrete_frechet_distance_coords` |
| κ–τ Fréchet is SE(3)-invariant | Supported (features are invariant) | `test_independent_rigid_motions_…` |
| Torsion handled periodically | Supported | probe `wrap_pair_expected_0p02`; `test_torsion_wrap_boundary` |
| Residue correspondence preserved | **Not preserved**: flexible coupling maps distinct same-length windows to 0 | probe `misregistration_x_y_vs_x_x`; Rust-path `test_frechet_misregistration_numpy_and_rust` |
| All window features retained | **No**: pairs (κ at vertex s+1+j, τ about bond s+1+j), length min(#κ, #τ) = W−3; last interior κ dropped; bond lengths unused | probe `last_kappa_changed_by_1rad`; `test_rust_path_ignores_last_window_kappa` |
| Consistent units and weights | κ, τ both rad, weights 1:1 (fixed in `engine.rs`); channel ranges differ, so equal weights are a modeling choice; "zero hyperparameters" is false (weights, window, radius) | code |
| Empty input | **Distance 0** (perfect match) for empty κ–τ, DTW and writhe spectra | probe `empty_vs_nonempty`, `dtw_empty`, `writhe_empty` |
| Non-finite input | **NaN pairs silently ignored** by f64::max/min; NaN only if every pair is NaN | probe `nan_torsion`, `all_nan_torsions` |
| Degenerate geometry | Exact collinear/repeated points raise and abort the whole `fit` (m1); float32-rounded collinear triples pass the 1e-9 check silently with arbitrary τ (LE-D6) | `test_degenerate_inputs_in_binding_and_fit`, `test_float32_rounded_collinear_passes_binding_check_silently` |
| Writhe-filtered candidate search | Approximate (M6); writhe distance truncates unequal lengths | probe `writhe_unequal_len`; recall in §6 |
| Exhaustive local retrieval (`query_subcurve`) | Linear scan over all frames (correctly described in README) | code `engine.rs::query_subcurve_k_nearest` |
| Scalar circular-mean / BC scanning | Separate mechanism (bimodality.rs); characterized in WP-A/WP-D (C1) | WP-A, WP-D reports |
| Historical benchmark advantages | Unsupported: synthetic ensembles, circular labels, unequal information (KNOWN_ISSUES C1, C2) | KNOWN_ISSUES |

Rust probe outputs (SYNTHETIC, rustc 1.98.1 (48a229cea 2026-09-01)):

| Probe | Output | Meaning |
|---|---|---|
| `empty_vs_nonempty` | 0.000000 | empty invariants compare as identical (0.0) to any curve |
| `last_kappa_changed_by_1rad` | 0.000000 | last interior turning angle omitted by min(len) pairing |
| `first_kappa_changed_by_1rad` | 1.000000 | control: change in a paired κ is detected |
| `nan_torsion` | 0.000000 | NaN pointwise distance ignored by f64::max/min; result as if the pair matched |
| `far_torsion_same_position` | 0.905539 | control: change in a paired κ is detected |
| `all_nan_torsions` | NaN | NaN only when every pair is NaN |
| `wrap_pair_expected_0p02` | 0.020000 | periodic torsion difference handled correctly |
| `misregistration_x_y_vs_x_x` | 0.000000 | flexible coupling maps distinct same-length sequences to 0 |
| `writhe_unequal_len` | 0.000000 | writhe spectrum distance truncates to the shorter vector |
| `writhe_empty` | 0.000000 | empty writhe spectrum distance is 0 |
| `dtw_empty` | 0.000000 | empty DTW distance is 0 |
| `kappa_len5_tau_len2` | 0.000000 | only min(#kappa,#tau) pairs are compared |

**AdK verification gates (Amendment 1 §5):** {'reference_check': True, 'frechet_rust_numpy': True, 'binding_failures_zero': True, 'nonfinite_zero_best_set': True, 'regret_failures_zero': True, 'undefined_below_1pct': True}. Rust vs NumPy Fréchet max |Δ| 1.1e-15; non-finite distances (T_φψ) {'M1': 0, 'M2': 0, 'K': 0, 'A': 0, 'F': 0, 'F_lock': 0, 'B1': 0, 'B2': 0, 'G_RMSD': 0, 'G_F': 0}; regret failures {'M1': 0, 'M2': 0, 'K': 0, 'A': 0, 'F': 0, 'F_lock': 0, 'B1': 0, 'B2': 0, 'G_RMSD': 0, 'G_F': 0}.

## 2. Retrieval: equal-window-weighted mean NR@10 (0 = target-optimal, 1 = random; lower is better)

Intervals: 95 % two-way block bootstrap of the same statistic, query blocks of 100 | 250 frames, window blocks of 20 starts, 2,000 paired resamples. Database frames are not resampled.

### Primary target T_φψ

| Method | open_start/W5 | open_start/W8 | closed_start/W5 | closed_start/W8 |
|---|---|---|---|---|
| M1 τ chord-RMS | 0.777 [0.700, 0.838] [0.696, 0.839] | 0.691 [0.610, 0.757] [0.610, 0.756] | 0.799 [0.746, 0.846] [0.742, 0.845] | 0.705 [0.617, 0.772] [0.620, 0.772] |
| M2 τ circular-RMS | 0.778 [0.702, 0.839] [0.697, 0.840] | 0.692 [0.611, 0.757] [0.611, 0.758] | 0.799 [0.746, 0.846] [0.742, 0.845] | 0.705 [0.617, 0.772] [0.620, 0.772] |
| K κ-only RMS | 0.755 [0.683, 0.819] [0.684, 0.821] | 0.687 [0.611, 0.758] [0.607, 0.761] | 0.786 [0.723, 0.836] [0.717, 0.836] | 0.711 [0.620, 0.777] [0.621, 0.780] |
| A residue-aligned κ–τ | 0.673 [0.586, 0.745] [0.582, 0.747] | 0.607 [0.528, 0.676] [0.524, 0.679] | 0.696 [0.623, 0.753] [0.618, 0.754] | 0.621 [0.525, 0.690] [0.527, 0.691] |
| F legacy κ–τ Fréchet (Rust) | 0.707 [0.623, 0.775] [0.620, 0.776] | 0.638 [0.562, 0.703] [0.562, 0.706] | 0.727 [0.661, 0.781] [0.658, 0.780] | 0.652 [0.557, 0.721] [0.558, 0.723] |
| F_lock legacy pairs, lock-step max | 0.707 [0.623, 0.775] [0.620, 0.776] | 0.633 [0.554, 0.698] [0.554, 0.702] | 0.727 [0.661, 0.781] [0.658, 0.780] | 0.644 [0.549, 0.713] [0.549, 0.714] |
| B1 local Cα–Cα distances | 0.694 [0.611, 0.764] [0.609, 0.766] | 0.660 [0.599, 0.718] [0.598, 0.720] | 0.729 [0.657, 0.783] [0.653, 0.783] | 0.676 [0.580, 0.741] [0.580, 0.744] |
| B2 window Cα RMSD | 0.674 [0.588, 0.746] [0.585, 0.748] | 0.640 [0.574, 0.701] [0.573, 0.704] | 0.700 [0.627, 0.757] [0.623, 0.757] | 0.654 [0.558, 0.723] [0.559, 0.725] |
| G_RMSD global Cα RMSD | 0.959 [0.910, 1.004] [0.911, 1.001] | 0.959 [0.898, 1.014] [0.897, 1.012] | 0.955 [0.911, 0.992] [0.906, 0.989] | 0.937 [0.878, 0.982] [0.877, 0.983] |
| G_F legacy full-chain Fréchet | 0.949 [0.902, 0.986] [0.896, 0.988] | 0.940 [0.888, 0.984] [0.885, 0.986] | 0.942 [0.888, 0.990] [0.886, 0.990] | 0.918 [0.841, 0.982] [0.840, 0.983] |

Valid query counts per window (min–max, distinct values): open_start/W5 1013–1013 (1); open_start/W8 1013–1013 (1); closed_start/W5 755–755 (1); closed_start/W8 755–755 (1). Undefined NR units: open_start/W5 0; open_start/W8 0; closed_start/W5 0; closed_start/W8 0.

### Complementary target T_pp (peptide-plane atoms)

| Method | open_start/W5 | open_start/W8 | closed_start/W5 | closed_start/W8 |
|---|---|---|---|---|
| M1 τ chord-RMS | 0.652 [0.569, 0.717] [0.570, 0.720] | 0.442 [0.382, 0.499] [0.378, 0.499] | 0.684 [0.637, 0.729] [0.635, 0.726] | 0.481 [0.427, 0.536] [0.426, 0.533] |
| M2 τ circular-RMS | 0.652 [0.569, 0.718] [0.571, 0.720] | 0.442 [0.382, 0.500] [0.379, 0.500] | 0.684 [0.637, 0.729] [0.635, 0.726] | 0.481 [0.427, 0.536] [0.426, 0.533] |
| K κ-only RMS | 0.699 [0.643, 0.750] [0.641, 0.753] | 0.645 [0.578, 0.706] [0.577, 0.706] | 0.707 [0.659, 0.743] [0.656, 0.743] | 0.634 [0.570, 0.679] [0.566, 0.680] |
| A residue-aligned κ–τ | 0.472 [0.401, 0.532] [0.403, 0.532] | 0.331 [0.286, 0.377] [0.281, 0.377] | 0.500 [0.445, 0.541] [0.442, 0.541] | 0.360 [0.310, 0.401] [0.308, 0.402] |
| F legacy κ–τ Fréchet (Rust) | 0.506 [0.433, 0.564] [0.436, 0.564] | 0.406 [0.361, 0.448] [0.360, 0.446] | 0.533 [0.484, 0.573] [0.480, 0.572] | 0.429 [0.377, 0.471] [0.375, 0.471] |
| F_lock legacy pairs, lock-step max | 0.506 [0.433, 0.564] [0.436, 0.564] | 0.362 [0.317, 0.406] [0.314, 0.407] | 0.533 [0.484, 0.573] [0.480, 0.572] | 0.391 [0.342, 0.433] [0.339, 0.434] |
| B1 local Cα–Cα distances | 0.468 [0.406, 0.520] [0.407, 0.522] | 0.272 [0.238, 0.308] [0.238, 0.307] | 0.498 [0.450, 0.535] [0.446, 0.536] | 0.297 [0.258, 0.332] [0.255, 0.332] |
| B2 window Cα RMSD | 0.431 [0.364, 0.487] [0.367, 0.487] | 0.197 [0.162, 0.232] [0.160, 0.232] | 0.454 [0.404, 0.493] [0.400, 0.493] | 0.216 [0.187, 0.242] [0.186, 0.241] |
| G_RMSD global Cα RMSD | 0.945 [0.888, 0.997] [0.884, 0.997] | 0.908 [0.811, 0.996] [0.810, 0.989] | 0.945 [0.904, 0.978] [0.902, 0.978] | 0.906 [0.852, 0.951] [0.853, 0.950] |
| G_F legacy full-chain Fréchet | 0.951 [0.894, 0.997] [0.892, 0.997] | 0.939 [0.868, 0.999] [0.866, 0.997] | 0.936 [0.882, 0.981] [0.883, 0.981] | 0.910 [0.845, 0.969] [0.839, 0.972] |

Valid query counts per window (min–max, distinct values): open_start/W5 1013–1013 (1); open_start/W8 1013–1013 (1); closed_start/W5 755–755 (1); closed_start/W8 755–755 (1). Undefined NR units: open_start/W5 0; open_start/W8 0; closed_start/W5 0; closed_start/W8 0.

### Secondary: overlap@10 (median across windows of window means) vs random

| Cell | Target | random | M1 | M2 | K | A | F | F_lock | B1 | B2 | G_RMSD | G_F |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| open_start/W5 | phipsi | 0.0099 | 0.0259 | 0.0259 | 0.0263 | 0.0413 | 0.0376 | 0.0376 | 0.0319 | 0.0402 | 0.0096 | 0.0102 |
| open_start/W5 | pp | 0.0099 | 0.0463 | 0.0463 | 0.0420 | 0.0963 | 0.0830 | 0.0830 | 0.0987 | 0.1216 | 0.0096 | 0.0102 |
| open_start/W8 | phipsi | 0.0099 | 0.0432 | 0.0432 | 0.0378 | 0.0627 | 0.0555 | 0.0572 | 0.0410 | 0.0482 | 0.0095 | 0.0105 |
| open_start/W8 | pp | 0.0099 | 0.1008 | 0.1009 | 0.0622 | 0.1554 | 0.1270 | 0.1369 | 0.2066 | 0.3033 | 0.0094 | 0.0100 |
| closed_start/W5 | phipsi | 0.0132 | 0.0328 | 0.0328 | 0.0329 | 0.0515 | 0.0430 | 0.0430 | 0.0401 | 0.0462 | 0.0148 | 0.0146 |
| closed_start/W5 | pp | 0.0132 | 0.0597 | 0.0597 | 0.0506 | 0.1102 | 0.0986 | 0.0986 | 0.1133 | 0.1366 | 0.0150 | 0.0147 |
| closed_start/W8 | phipsi | 0.0132 | 0.0490 | 0.0490 | 0.0517 | 0.0728 | 0.0623 | 0.0646 | 0.0511 | 0.0592 | 0.0155 | 0.0150 |
| closed_start/W8 | pp | 0.0132 | 0.1185 | 0.1187 | 0.0742 | 0.1691 | 0.1395 | 0.1494 | 0.2371 | 0.3245 | 0.0179 | 0.0159 |

## 3. Paired contrasts (ΔNR@10 = X − Y on identical units; negative favours X)

Statement rule (Amendment 1 §4): τ_int ≤ 100 → 100-frame interval; ≤ 250 → 250-frame interval; > 250 or not estimable → **dependence_unresolved**. Both intervals shown. F − A is an overall comparison; only F − F_lock isolates flexible coupling; F_lock − A combines aggregation, the omitted κ and channel normalization.

**Target phipsi**

| Contrast | open_start/W5 | open_start/W8 | closed_start/W5 | closed_start/W8 |
|---|---|---|---|---|
| H1_A_minus_M2 | -0.105 [-0.132, -0.082] [-0.132, -0.082] τ_int 2 → **below_0** | -0.085 [-0.112, -0.062] [-0.112, -0.062] τ_int 3 → **below_0** | -0.103 [-0.134, -0.078] [-0.137, -0.077] τ_int 2 → **below_0** | -0.084 [-0.110, -0.062] [-0.110, -0.062] τ_int 2 → **below_0** |
| H1_A_minus_K | -0.083 [-0.106, -0.058] [-0.109, -0.058] τ_int 9 → **below_0** | -0.080 [-0.103, -0.055] [-0.107, -0.051] τ_int 76 → **below_0** | -0.090 [-0.111, -0.068] [-0.111, -0.067] τ_int 5 → **below_0** | -0.090 [-0.116, -0.060] [-0.115, -0.055] τ_int 11 → **below_0** |
| overall_F_minus_A | 0.034 [0.025, 0.048] [0.025, 0.047] τ_int 1 → **above_0** | 0.031 [0.020, 0.046] [0.019, 0.045] τ_int 5 → **above_0** | 0.031 [0.023, 0.040] [0.023, 0.040] τ_int 1 → **above_0** | 0.031 [0.028, 0.035] [0.027, 0.036] τ_int 4 → **above_0** |
| coupling_F_minus_Flock | 0.000 [0.000, 0.000] [0.000, 0.000] τ_int — → **dependence_unresolved** | 0.005 [-0.004, 0.013] [-0.005, 0.013] τ_int 3 → **unresolved** | 0.000 [0.000, 0.000] [0.000, 0.000] τ_int — → **dependence_unresolved** | 0.009 [0.006, 0.012] [0.006, 0.012] τ_int 1 → **above_0** |
| aggregation_boundary_Flock_minus_A | 0.034 [0.025, 0.048] [0.025, 0.047] τ_int 1 → **above_0** | 0.026 [0.018, 0.039] [0.018, 0.039] τ_int 5 → **above_0** | 0.031 [0.023, 0.040] [0.023, 0.040] τ_int 1 → **above_0** | 0.023 [0.019, 0.027] [0.019, 0.028] τ_int 3 → **above_0** |
| H4_GF_minus_F | 0.242 [0.198, 0.289] [0.193, 0.292] τ_int 90 → **above_0** | 0.302 [0.257, 0.351] [0.252, 0.351] τ_int 85 → **above_0** | 0.215 [0.193, 0.242] [0.195, 0.241] τ_int 8 → **above_0** | 0.266 [0.240, 0.303] [0.241, 0.302] τ_int 8 → **above_0** |
| H4_GRMSD_minus_B2 | 0.285 [0.226, 0.353] [0.225, 0.355] τ_int 4 → **above_0** | 0.319 [0.259, 0.390] [0.259, 0.388] τ_int 6 → **above_0** | 0.255 [0.218, 0.300] [0.215, 0.301] τ_int 20 → **above_0** | 0.283 [0.243, 0.338] [0.241, 0.341] τ_int 26 → **above_0** |
| descriptive_F_minus_B2 | 0.033 [0.021, 0.048] [0.020, 0.047] τ_int 1 → **above_0** | -0.002 [-0.019, 0.013] [-0.018, 0.013] τ_int 4 → **unresolved** | 0.027 [0.018, 0.037] [0.018, 0.038] τ_int 2 → **above_0** | -0.002 [-0.009, 0.007] [-0.009, 0.009] τ_int 6 → **unresolved** |
| H5_B2_minus_A | 0.001 [-0.007, 0.007] [-0.007, 0.007] τ_int 1 → **unresolved** | 0.034 [0.014, 0.052] [0.013, 0.052] τ_int 7 → **above_0** | 0.004 [-0.000, 0.007] [-0.001, 0.007] τ_int 1 → **unresolved** | 0.033 [0.027, 0.039] [0.026, 0.039] τ_int 3 → **above_0** |
| H5_B1_minus_A | 0.022 [0.003, 0.035] [0.002, 0.035] τ_int 2 → **above_0** | 0.054 [0.025, 0.079] [0.023, 0.081] τ_int 3 → **above_0** | 0.033 [0.024, 0.042] [0.023, 0.043] τ_int 2 → **above_0** | 0.055 [0.042, 0.069] [0.041, 0.068] τ_int 3 → **above_0** |
| M1_minus_M2 | -0.001 [-0.003, 0.000] [-0.003, 0.000] τ_int 6 → **unresolved** | -0.001 [-0.003, 0.000] [-0.003, 0.000] τ_int 4 → **unresolved** | -0.000 [-0.000, 0.000] [-0.000, 0.000] τ_int 1 → **unresolved** | -0.000 [-0.000, 0.000] [-0.000, 0.000] τ_int 1 → **unresolved** |

**Target pp**

| Contrast | open_start/W5 | open_start/W8 | closed_start/W5 | closed_start/W8 |
|---|---|---|---|---|
| H1_A_minus_M2 | -0.180 [-0.204, -0.159] [-0.204, -0.159] τ_int 2 → **below_0** | -0.111 [-0.137, -0.085] [-0.136, -0.086] τ_int 3 → **below_0** | -0.184 [-0.208, -0.164] [-0.210, -0.162] τ_int 17 → **below_0** | -0.121 [-0.143, -0.100] [-0.142, -0.101] τ_int 6 → **below_0** |
| H1_A_minus_K | -0.227 [-0.257, -0.198] [-0.256, -0.197] τ_int 2 → **below_0** | -0.314 [-0.356, -0.266] [-0.356, -0.267] τ_int 11 → **below_0** | -0.207 [-0.228, -0.184] [-0.229, -0.183] τ_int 5 → **below_0** | -0.274 [-0.310, -0.238] [-0.311, -0.234] τ_int 23 → **below_0** |
| overall_F_minus_A | 0.034 [0.029, 0.042] [0.028, 0.042] τ_int 1 → **above_0** | 0.075 [0.061, 0.090] [0.061, 0.092] τ_int 88 → **above_0** | 0.033 [0.028, 0.041] [0.028, 0.041] τ_int 1 → **above_0** | 0.070 [0.057, 0.084] [0.056, 0.085] τ_int 37 → **above_0** |
| coupling_F_minus_Flock | 0.000 [0.000, 0.000] [0.000, 0.000] τ_int — → **dependence_unresolved** | 0.044 [0.031, 0.058] [0.031, 0.059] τ_int 32 → **above_0** | 0.000 [0.000, 0.000] [0.000, 0.000] τ_int — → **dependence_unresolved** | 0.039 [0.026, 0.052] [0.026, 0.052] τ_int 9 → **above_0** |
| aggregation_boundary_Flock_minus_A | 0.034 [0.029, 0.042] [0.028, 0.042] τ_int 1 → **above_0** | 0.031 [0.027, 0.035] [0.027, 0.036] τ_int 54 → **above_0** | 0.033 [0.028, 0.041] [0.028, 0.041] τ_int 1 → **above_0** | 0.031 [0.029, 0.035] [0.028, 0.035] τ_int 17 → **above_0** |
| H4_GF_minus_F | 0.445 [0.414, 0.475] [0.412, 0.480] τ_int 102 → **above_0** | 0.533 [0.496, 0.572] [0.495, 0.576] τ_int 92 → **above_0** | 0.402 [0.378, 0.430] [0.380, 0.430] τ_int 10 → **above_0** | 0.480 [0.445, 0.519] [0.437, 0.522] τ_int 50 → **above_0** |
| H4_GRMSD_minus_B2 | 0.514 [0.477, 0.551] [0.473, 0.552] τ_int 56 → **above_0** | 0.711 [0.643, 0.779] [0.641, 0.775] τ_int 69 → **above_0** | 0.491 [0.463, 0.519] [0.463, 0.518] τ_int 22 → **above_0** | 0.690 [0.652, 0.726] [0.653, 0.725] τ_int 27 → **above_0** |
| descriptive_F_minus_B2 | 0.075 [0.066, 0.084] [0.066, 0.085] τ_int 2 → **above_0** | 0.209 [0.187, 0.237] [0.186, 0.237] τ_int 102 → **above_0** | 0.079 [0.073, 0.087] [0.073, 0.089] τ_int 40 → **above_0** | 0.213 [0.187, 0.239] [0.182, 0.247] τ_int 85 → **above_0** |
| H5_B2_minus_A | -0.041 [-0.049, -0.034] [-0.049, -0.034] τ_int 2 → **below_0** | -0.135 [-0.163, -0.114] [-0.162, -0.113] τ_int 27 → **below_0** | -0.046 [-0.054, -0.039] [-0.056, -0.039] τ_int 49 → **below_0** | -0.144 [-0.166, -0.122] [-0.169, -0.120] τ_int 78 → **below_0** |
| H5_B1_minus_A | -0.004 [-0.021, 0.012] [-0.021, 0.013] τ_int 2 → **unresolved** | -0.059 [-0.093, -0.033] [-0.093, -0.033] τ_int 12 → **below_0** | -0.003 [-0.018, 0.012] [-0.019, 0.013] τ_int 15 → **unresolved** | -0.063 [-0.077, -0.050] [-0.081, -0.048] τ_int 56 → **below_0** |
| M1_minus_M2 | -0.000 [-0.001, 0.000] [-0.001, 0.000] τ_int 4 → **unresolved** | -0.001 [-0.002, 0.000] [-0.002, 0.000] τ_int 15 → **unresolved** | 0.000 [-0.000, 0.000] [-0.000, 0.000] τ_int 1 → **unresolved** | 0.000 [-0.000, 0.000] [-0.000, 0.000] τ_int 1 → **unresolved** |

## 4. Development classification of A and F (primary target; Amendment-1 gates)

| Method | 100-frame cells → outcome | 250-frame cells → outcome | τ_int vs best (frames) | Final |
|---|---|---|---|---|
| A | within_margin, within_margin, within_margin, within_margin → retain_as_heldout_candidate | within_margin, within_margin, within_margin, within_margin → retain_as_heldout_candidate | —, —, —, — | **insufficient_evidence (dependence_unresolved)** |
| F | inconclusive, inconclusive, inconclusive, inconclusive → inconclusive | inconclusive, inconclusive, inconclusive, inconclusive → inconclusive | 1, 5, 1, 4 | **inconclusive** |

| Cell | Best observed | L_A [100] [250] | L_F [100] [250] |
|---|---|---|---|
| open_start/W5 | A | 0.000 [0.000, 0.022] [0.000, 0.022] | 0.105 [0.077, 0.150] [0.078, 0.149] |
| open_start/W8 | A | 0.000 [0.000, 0.000] [0.000, 0.000] | 0.080 [0.054, 0.109] [0.052, 0.110] |
| closed_start/W5 | A | 0.000 [0.000, 0.000] [0.000, 0.004] | 0.101 [0.081, 0.124] [0.082, 0.126] |
| closed_start/W8 | A | 0.000 [0.000, 0.000] [0.000, 0.000] | 0.083 [0.064, 0.105] [0.064, 0.106] |

## 5. Failure cases: windows where F − A (T_φψ) is largest

| Cell | i | residues | NR_F | NR_F_lock | NR_A | Δ(F − A) |
|---|---|---|---|---|---|---|
| open_start/W5 | 119 | 119–123 | 1.385 | 1.385 | 0.331 | 1.054 |
| open_start/W5 | 150 | 150–154 | 0.816 | 0.816 | 0.416 | 0.400 |
| open_start/W5 | 98 | 98–102 | 0.883 | 0.883 | 0.641 | 0.241 |
| open_start/W5 | 132 | 132–136 | 0.670 | 0.670 | 0.509 | 0.161 |
| open_start/W5 | 10 | 10–14 | 1.636 | 1.636 | 1.489 | 0.146 |
| open_start/W8 | 116 | 116–123 | 1.158 | 1.192 | 0.264 | 0.894 |
| open_start/W8 | 147 | 147–154 | 0.659 | 0.659 | 0.376 | 0.283 |
| open_start/W8 | 138 | 138–145 | 0.433 | 0.262 | 0.213 | 0.219 |
| open_start/W8 | 95 | 95–102 | 0.840 | 0.836 | 0.637 | 0.203 |
| open_start/W8 | 140 | 140–147 | 0.485 | 0.401 | 0.334 | 0.151 |
| closed_start/W5 | 131 | 131–135 | 0.783 | 0.783 | 0.478 | 0.304 |
| closed_start/W5 | 118 | 118–122 | 1.004 | 1.004 | 0.722 | 0.281 |
| closed_start/W5 | 134 | 134–138 | 0.479 | 0.479 | 0.216 | 0.263 |
| closed_start/W5 | 125 | 125–129 | 0.634 | 0.634 | 0.429 | 0.205 |
| closed_start/W5 | 133 | 133–137 | 0.491 | 0.491 | 0.297 | 0.195 |
| closed_start/W8 | 115 | 115–122 | 0.999 | 0.995 | 0.702 | 0.297 |
| closed_start/W8 | 123 | 123–130 | 0.342 | 0.342 | 0.198 | 0.144 |
| closed_start/W8 | 207 | 207–214 | 0.441 | 0.444 | 0.308 | 0.133 |
| closed_start/W8 | 169 | 169–176 | 0.593 | 0.562 | 0.466 | 0.127 |
| closed_start/W8 | 11 | 11–18 | 0.599 | 0.492 | 0.475 | 0.124 |

Windows where NR_F > NR_A and where NR_F > NR_F_lock: open_start/W5: 195, 0 of 210; open_start/W8: 200, 144 of 207; closed_start/W5: 195, 0 of 210; closed_start/W8: 198, 143 of 207.

## 6. Writhe cascade vs exhaustive full-chain Fréchet (same final distance)

- open_start: mean recall@10 1.000, minimum 1.000, queries with recall 1: 1.000; queries returning < 10 hits: 0.
- closed_start: mean recall@10 1.000, minimum 1.000, queries with recall 1: 1.000; queries returning < 10 hits: 0.

## 7. Computation (single thread; open_start, W = 8, windows 1–20 unless stated)

CPU: Intel(R) Core(TM) i5-8400T CPU @ 1.70GHz. Rust vs NumPy differences are implementation differences, not representation effects.

| Implementation | Steady state (s, median of 3) | Work unit | Peak RSS above baseline (MB) |
|---|---|---|---|
| np_features | 0.046 | 20 windows × all queries × all candidates | 0.0 |
| M1 | 4.534 | 20 windows × all queries × all candidates | 0.0 |
| M2 | 3.289 | 20 windows × all queries × all candidates | 0.0 |
| K | 0.946 | 20 windows × all queries × all candidates | 0.0 |
| A | 4.282 | 20 windows × all queries × all candidates | 0.0 |
| F_lock | 4.440 | 20 windows × all queries × all candidates | 0.0 |
| F_numpy | 31.867 | 20 windows × all queries × all candidates | 75.7 |
| B1 | 2.392 | 20 windows × all queries × all candidates | 0.0 |
| B2 | 43.354 | 20 windows × all queries × all candidates | 0.0 |
| F_rust_fit | 0.355 | one index fit over the database | 0.0 |
| F_rust_query | 6.031 | 20 windows × all queries × all candidates | 0.0 |
| GF_rust_query | 24.697 | 50 full-chain queries | 0.0 |
| cascade_query | 24.737 | 50 full-chain queries | 0.0 |

Setup: binding feature extraction for all 2253 frames 0.74 s; Rust index fit 0.355 s.

## 8. Prior work and novelty

Crossref verification establishes bibliographic metadata only. "Content evidence" states whether what is attributed to a work was checked against the primary source in this audit; novelty judgments resting on unchecked content are **provisional**.

| Key | Crossref-verified | DOI | Year |
|---|---|---|---|
| kabsch1976 | True | 10.1107/s0567739476001873 | 1976 |
| horn1987 | True | 10.1364/josaa.4.000629 | 1987 |
| theobald2005 | True | 10.1107/s0108767305015266 | 2005 |
| levitt1976 | True | 10.1016/0022-2836(76)90004-8 | 1976 |
| ramachandran1963 | True | 10.1016/s0022-2836(63)80023-6 | 1963 |
| oldfield1994 | False | — (not cited) | — |
| rackovsky1978 | True | 10.1021/ma60066a020 | 1978 |
| hu2011 | True | 10.1103/physreve.83.061908 | 2011 |
| altgodau1995 | True | 10.1142/s0218195995000064 | 1995 |
| jiang2008 | True | 10.1142/9781860947995_0016 | 2007 |
| sakoe1978 | True | 10.1109/tassp.1978.1163055 | 1978 |
| debrevern2000 | True | 10.1002/1097-0134(20001115)41:3<271::aid-prot10>3.0.co;2-z | 2000 |
| kolodny2002 | True | 10.1016/s0022-2836(02)00942-7 | 2002 |
| vankempen2024 | True | 10.1038/s41587-023-01773-0 | 2024 |
| mu2005 | True | 10.1002/prot.20310 | 2005 |
| altis2007 | True | 10.1063/1.2746330 | 2007 |
| rogen2003 | True | 10.1073/pnas.2636460100 | 2003 |
| michaudagrawal2011 | True | 10.1002/jcc.21787 | 2011 |
| mcgibbon2015 | True | 10.1016/j.bpj.2015.08.015 | 2015 |
| yianilos1993 | False | — (not cited) | — |
| mirarchi2024 | True | 10.1038/s41597-024-04140-z | 2024 |

| Topic | References (verified keys only) | Content attributed | Content evidence |
|---|---|---|---|
| Cα virtual-bond / virtual-dihedral representation | levitt1976 | Cα-only angles and dihedrals as a reduced backbone description | prior knowledge, not checked in this audit |
| Discrete curvature/torsion (Frenet) description of protein chains | rackovsky1978, hu2011 | Differential-geometry descriptors of protein backbones | prior knowledge, not checked in this audit |
| φ/ψ as local backbone description | ramachandran1963 | Backbone dihedral map | prior knowledge, not checked in this audit |
| Fréchet distance; DTW | altgodau1995, sakoe1978, jiang2008 | Coupling distances for curves/sequences, incl. discrete Fréchet for protein alignment | prior knowledge, not checked in this audit |
| Optimal superposition RMSD | kabsch1976, horn1987, theobald2005 | Closed-form rotation / quaternion RMSD | prior knowledge; Horn vs Kabsch agreement checked numerically here |
| Circular (sin, cos) dihedral embeddings | mu2005, altis2007 | dPCA on sin/cos of dihedrals | prior knowledge, not checked in this audit |
| Structural alphabets / local fragment libraries | debrevern2000, kolodny2002, vankempen2024 | Discretized local backbone states; fragment libraries; 3Di local descriptors for search | prior knowledge, not checked in this audit |
| Writhe / Gauss-integral descriptors | rogen2003 | Gauss integrals for structure classification | prior knowledge, not checked in this audit |
| MD analysis tools | michaudagrawal2011, mcgibbon2015 | Dihedral, RMSD and distance analyses | MDAnalysis used here (2.10.0); content otherwise prior knowledge |
| Held-out dataset | mirarchi2024 | mdCATH content and builder | checked against paper and code in WP-D |

| Category | Items |
|---|---|
| Established mathematics | discrete turning angles and dihedrals; circular distance; Fréchet/DTW coupling; quaternion/SVD superposition |
| Existing software functionality | dihedral, distance and RMSD analysis in MD toolkits; structure search with local descriptors (provisional: content not checked here) |
| Engineering integration (TopoFold) | Rust implementation of κ–τ extraction, exhaustive window scan, writhe VP-tree cascade; Python binding |
| New empirical findings (this audit, one protein) | §2–§6 of this report |
| Unsupported novelty claims | combined κ–τ as new; SE(3) invariance as distinctive; "zero hyperparameters"; Rust as a scientific contribution |

## 9. Post-hoc observations (after outcomes were seen; they do not change any pre-declared outcome)

- **Gate 7 and method A.** The pre-declared outcome for A is `insufficient_evidence (dependence_unresolved)`. τ_int of NR_A − NR_best was not estimable ([None, None, None, None]) because A is the observed best in every primary cell, so the series is identically zero. The gate therefore flags an exact zero, not unresolved dependence. For the nearest competitor the dependence diagnostic is short: τ_int of B2 − A on T_φψ = 1, 7, 1, 3 frames. A future amendment could take τ_int against the runner-up when X is best; it is not applied here.
- **F − F_lock at W = 5 is exactly 0 by construction.** With W − 3 = 2 legacy pairs every monotone coupling contains both diagonal points, so the discrete Fréchet equals the lock-step maximum; the "dependence unresolved" label there reflects a constant series, not uncertainty.
- **Writhe cascade.** Recall@10 vs exhaustive was 1.000 / 1.000, but its time per query (24.7 s per 50 queries) equals the exhaustive scan (24.7 s): with the default radius 0.5 after the writhe-normalization fix (M1), equal timing is consistent with the filter admitting most frames (candidate counts were not recorded). Exact recall here is therefore not evidence of an effective approximate index.
- Valid query counts were identical across windows and no NR unit was undefined, so the Amendment 1 §1 estimand change does not alter any number here; it matters for data with exclusions.

## 10. Conclusions (development tier, one protein, two trajectories)

**Implemented correctly (verified):** the legacy κ–τ extraction and exhaustive `query_subcurve` path match independent NumPy implementations (§1); torsion periodicity and SE(3) invariance hold.

**Legacy defects documented (not repaired):** last interior κ dropped by the pairing; flexible coupling maps misregistered windows to 0; empty inputs score 0; NaN pairs silently ignored; float32-rounded collinear triples pass the degeneracy check; one collinear frame aborts `fit`; writhe distance truncates unequal lengths; in the one measured workload (single thread, 50 full-chain queries, open_start, default radius 0.5) the writhe cascade took the same time as the exhaustive scan [corrected E1: scope].

**What preserved more information:** combining κ and τ (A) retrieved closer conformations than either channel alone under both targets (H1 `below_0` in 16 of 16 cell × contrast statements).

**What improved retrieval:** best method by cell mean NR@10 — T_φψ: open_start/W5 A, open_start/W8 A, closed_start/W5 A, closed_start/W8 A; T_pp: open_start/W5 B2, open_start/W8 B2, closed_start/W5 B2, closed_start/W8 B2. A − B2 (H5) on T_φψ: unresolved, above_0, unresolved, above_0 (B2 − A; above_0 favours A); on T_pp: below_0, below_0, below_0, below_0 (below_0 favours B2). The winner depends on the target family: the angular target favours the angular combination, the Cartesian peptide-plane target favours Cartesian RMSD (H6). Gains over random are moderate for every method (NR@10 well above 0).

**Legacy Fréchet (F):** worse than A in every cell for both targets (overall F − A above_0, above_0, above_0, above_0, above_0, above_0, above_0, above_0). Of that difference, flexible coupling (F − F_lock) contributes 0 at W = 5 by construction and unresolved, above_0, above_0, above_0 at W = 8 (T_φψ open/closed, T_pp open/closed). Which component is larger depends on target and window size (§13, corrected E1): F_lock − A (max aggregation, omitted κ and normalization together) is the whole difference at W = 5 and the larger part on T_φψ at W = 8, but flexible coupling is the larger part on T_pp at W = 8. Development outcome F (declared, unchanged): `inconclusive`.

**Local vs global (H4):** above_0/above_0, above_0/above_0, above_0/above_0, above_0/above_0, above_0/above_0, above_0/above_0, above_0/above_0, above_0/above_0 (G_F − F / G_RMSD − B2). Global rankings, including the legacy full-chain Fréchet, are close to random for local retrieval; historical advantages of the engine over global baselines are consistent with local-vs-global information rather than with the κ–τ Fréchet representation itself.

**Generalized beyond development data:** nothing. No held-out data were scored.

**Scientifically unvalidated:** any biological relevance (both targets are geometric reconstructions from the same trajectories); performance on other proteins, force fields or temperatures; database-frame variability; ρ as a usefulness threshold; novelty statements resting on unchecked prior-work content.

## 11. Recommendation (engine vs loader)

1. **Engineering decision [corrected E1]: deprioritize the legacy κ–τ Fréchet for same-sequence local window retrieval.** This is a prioritization for the tested task, not a scientific rejection: F's declared development outcome is formally `inconclusive` (its relative-loss interval straddles ρ) and stays so. The basis is that, on AdK, the residue-aligned combination of the same features (A) had lower mean NR@10 in every cell for both targets; F − B2 was never below 0 (above 0 in six cells, unresolved in the two T_φψ W = 8 cells, where F's point estimate was marginally lower); and flexible coupling was zero or adverse where resolved. Tasks needing unequal-length or cross-sequence matching were not tested.
2. **Carry A and B2 forward as competing candidates**, with B1 and M1/M2 as baselines; no single method is preferred because the ranking depends on the target family. Adding A to the held-out method set changes the frozen held-out protocol and requires a reviewed amendment **before** held-out scoring.
3. **Build the minimal mdCATH reader/adapter next** (established reader such as h5py, explicit residue mapping, Å units, box handling, provenance), once the schema-sample download and the h5py installation are authorized; neither is authorized at present. That is the remaining dependency for generalization.
4. Document the legacy defects in `docs/KNOWN_ISSUES.md` in a separate, reviewed change; no production code was modified here.

## 12. Limitations

- One protein (AdK), two trajectories, already used in WP-B/WP-D; M1, M2, K, B1, B2 results were known before this protocol. A, F, F_lock, the global methods and T_pp are new.
- Both targets are derived from the same trajectory; T_φψ favours angular and T_pp Cartesian methods.
- Intervals cover query-frame and residue-position resampling only; the database pool is fixed; 250-frame blocks give only 3–5 query blocks per trajectory; τ_int estimates are themselves uncertain; WP-B observed slow drift that blocks of this length may not capture.
- No multiplicity adjustment; statements are descriptive.
- Timing compares a NumPy research implementation with the Rust engine path; peak RSS is dominated by a shared baseline.

## 13. Editorial closeout E1 (interpretation corrections; results and declared outcomes unchanged)

Recorded in `results/le_editorial_log.json`. No experiment was rerun and no protocol or analysis changed.

**E1.1 Dependence flag for A was structural, not measured.** The declared outcome for A remains `insufficient_evidence (dependence_unresolved)`. Gate 7 computes τ_int of NR_A − NR_best with best = the observed best method. A was that best in every primary cell, so the series compared A with itself and was identically zero; the estimator returned no value. No long temporal dependence was measured for A. The flag is therefore an artefact of a self-comparison in the gate definition, not evidence of dependence; equally, it is not evidence that dependence is short. The diagnostic was simply undefined.

**E1.2 Attribution of the Fréchet gap, per target and window size.** F − A = (F − F_lock) + (F_lock − A) exactly on identical units. Values are paired equal-window means of ΔNR@10.

| Cell | F − A | F − F_lock (coupling) | F_lock − A (aggregation + omitted κ + normalization) | Larger component |
|---|---|---|---|---|
| open_start/W5/phipsi | 0.0343 (above_0) | 0.0000 (identity (exact 0)) | 0.0343 (above_0) | F_lock − A (coupling ≡ 0 by construction at W = 5) |
| open_start/W5/pp | 0.0340 (above_0) | 0.0000 (identity (exact 0)) | 0.0340 (above_0) | F_lock − A (coupling ≡ 0 by construction at W = 5) |
| open_start/W8/phipsi | 0.0315 (above_0) | 0.0050 (unresolved) | 0.0264 (above_0) | F_lock − A |
| open_start/W8/pp | 0.0747 (above_0) | 0.0441 (above_0) | 0.0306 (above_0) | coupling |
| closed_start/W5/phipsi | 0.0306 (above_0) | 0.0000 (identity (exact 0)) | 0.0306 (above_0) | F_lock − A (coupling ≡ 0 by construction at W = 5) |
| closed_start/W5/pp | 0.0330 (above_0) | 0.0000 (identity (exact 0)) | 0.0330 (above_0) | F_lock − A (coupling ≡ 0 by construction at W = 5) |
| closed_start/W8/phipsi | 0.0314 (above_0) | 0.0087 (above_0) | 0.0227 (above_0) | F_lock − A |
| closed_start/W8/pp | 0.0697 (above_0) | 0.0386 (above_0) | 0.0311 (above_0) | coupling |

The earlier sentence "the larger part is F_lock − A" was correct for T_φψ and for W = 5, and wrong for T_pp at W = 8, where coupling was larger. The components of F_lock − A are not separated.

**E1.3 Deprioritizing F is an engineering decision** for same-sequence local window retrieval (§11 item 1). F's formal development outcome remains `inconclusive`.

**E1.4 Timing statements are limited to the implementations and workload measured** (single thread, Intel(R) Core(TM) i5-8400T CPU @ 1.70GHz, open_start, W = 8, windows 1–20, 1013 queries × 1015 candidates; full-chain items: 50 queries). They compare particular NumPy research implementations with the Rust engine path through the Python binding, which includes per-call feature extraction and allocation. They are not statements about the representations, other implementations, thread counts, index radii or data sizes.

**E1.5 τ_int diagnostics above about 1.5 % of the series length are not resolvable (found after the results, while preparing held-out diagnostics).** A SYNTHETIC AR(1) calibration of the same estimator (`wpd_local_retrieval/results/SYNTHETIC_tau_int_calibration.json`) shows severe underestimation (median estimate < 0.8 × true) whenever τ̂/n ≥ 0.017; e.g. true τ_int 199 frames was estimated at a median of 66.6 from 1000 frames. The Amendment 1 §4 rule therefore cannot certify block adequacy for τ̂ above 0.015 n (15.2 frames for open_start, 11.3 for closed_start). 28 of 84 estimable contrast diagnostics exceed this; their statements (listed below) rest on block lengths whose adequacy is not established, and true dependence may exceed both tested block lengths. Declared statements and outcomes are unchanged; this is an interpretation caveat. On T_φψ, all F − A, F − F_lock, F_lock − A and A − M2 diagnostics are below this heuristic cutoff (see E1.6 for what that does and does not mean).

**E1.6 Correction to E1.5 wording.** The 0.015 n cutoff is a limited heuristic from an AR(1) simulation of this estimator, not a reliability guarantee. A τ̂ below it does not establish that dependence is short: processes with drift or several timescales (WP-B observed drift in these runs) can yield small τ̂, and in the AR(1) simulation underestimation already begins below the cutoff for short series. Diagnostics below the cutoff are therefore "not flagged by the heuristic", not "resolved"; the listed ones above it are flagged.

| Cell | Contrast | τ̂ (frames) | Declared statement |
|---|---|---|---|
| open_start/W5/phipsi | H4_GF_minus_F | 90 | above_0 |
| open_start/W5/pp | H4_GF_minus_F | 102 | above_0 |
| open_start/W5/pp | H4_GRMSD_minus_B2 | 56 | above_0 |
| open_start/W8/phipsi | H1_A_minus_K | 76 | below_0 |
| open_start/W8/phipsi | H4_GF_minus_F | 85 | above_0 |
| open_start/W8/pp | overall_F_minus_A | 88 | above_0 |
| open_start/W8/pp | coupling_F_minus_Flock | 32 | above_0 |
| open_start/W8/pp | aggregation_boundary_Flock_minus_A | 54 | above_0 |
| open_start/W8/pp | H4_GF_minus_F | 92 | above_0 |
| open_start/W8/pp | H4_GRMSD_minus_B2 | 69 | above_0 |
| open_start/W8/pp | descriptive_F_minus_B2 | 102 | above_0 |
| open_start/W8/pp | H5_B2_minus_A | 27 | below_0 |
| open_start/W8/pp | M1_minus_M2 | 15 | unresolved |
| closed_start/W5/phipsi | H4_GRMSD_minus_B2 | 20 | above_0 |
| closed_start/W5/pp | H1_A_minus_M2 | 17 | below_0 |
| closed_start/W5/pp | H4_GRMSD_minus_B2 | 22 | above_0 |
| closed_start/W5/pp | descriptive_F_minus_B2 | 40 | above_0 |
| closed_start/W5/pp | H5_B2_minus_A | 49 | below_0 |
| closed_start/W5/pp | H5_B1_minus_A | 15 | unresolved |
| closed_start/W8/phipsi | H4_GRMSD_minus_B2 | 26 | above_0 |
| closed_start/W8/pp | H1_A_minus_K | 23 | below_0 |
| closed_start/W8/pp | overall_F_minus_A | 37 | above_0 |
| closed_start/W8/pp | aggregation_boundary_Flock_minus_A | 17 | above_0 |
| closed_start/W8/pp | H4_GF_minus_F | 50 | above_0 |
| closed_start/W8/pp | H4_GRMSD_minus_B2 | 27 | above_0 |
| closed_start/W8/pp | descriptive_F_minus_B2 | 85 | above_0 |
| closed_start/W8/pp | H5_B2_minus_A | 78 | below_0 |
| closed_start/W8/pp | H5_B1_minus_A | 56 | below_0 |

