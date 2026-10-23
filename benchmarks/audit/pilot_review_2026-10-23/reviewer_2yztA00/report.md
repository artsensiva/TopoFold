# Descriptive ensemble comparison: Reviewer pilot: mdCATH 2yztA00, 320 K vs 348 K; exposed development data, descriptive only

Generated from the tables in `tables/`. **Descriptive only:** no p-values, significance labels or uncertainty intervals. Differences between groups can reflect geometry, dynamics, sampling or starting conditions; replica variability and temporal diagnostics are shown so that readers can judge this. Marginal per-window distributions can miss joint changes.

**Pooling:** each group's pooled distribution concatenates its members' retained frames, so members are weighted by their frame counts (`member_frames` in `tables/distribution.csv`).

**Circular summaries:** the resultant length R is reported; mean directions (and differences involving them) are shown as unavailable when R < 1e-06, a numerical conditioning rule, not a biological criterion.

## Inputs

- **320K**: 5 member(s): 320K_0 (117 frames, source frames 13–129, frame identity not_verifiable); 320K_1 (117 frames, source frames 13–129, frame identity not_verifiable); 320K_2 (117 frames, source frames 13–129, frame identity not_verifiable); 320K_3 (117 frames, source frames 13–129, frame identity not_verifiable); 320K_4 (117 frames, source frames 13–129, frame identity not_verifiable)
- **348K**: 5 member(s): 348K_0 (117 frames, source frames 13–129, frame identity not_verifiable); 348K_1 (117 frames, source frames 13–129, frame identity not_verifiable); 348K_2 (117 frames, source frames 13–129, frame identity not_verifiable); 348K_3 (117 frames, source frames 13–129, frame identity not_verifiable); 348K_4 (117 frames, source frames 13–129, frame identity not_verifiable)
- Core: native extension resolved (SHA-256 d90e7ed4f1d8036a…); analysis options {'bins': 36, 'max_lag': 50, 'pbc': 'none'} (per-member reading options, e.g. frame slice and PBC policy, are in `manifest.json`)

## Exclusions

0 channel–window combinations excluded (invalid frames in at least one member); listed in `tables/exclusions.csv`.

## 320K vs 348K: largest pooled JS distances per channel (descriptive)

| Channel | Window | Pooled JS | Between-member JS (median) | Within-member JS median 320K / 348K | Center difference |
|---|---|---|---|---|---|
| kappa | 0:47-49 | 0.399 | 0.429 | 0.253 / 0.457 | +0.224 |
| kappa | 0:12-14 | 0.255 | 0.162 | 0.193 / 0.125 | -0.075 |
| kappa | 0:46-48 | 0.238 | 0.373 | 0.218 / 0.247 | +0.065 |
| kappa | 0:63-65 | 0.221 | 0.447 | 0.483 / 0.484 | +0.112 |
| kappa | 0:61-63 | 0.214 | 0.437 | 0.476 / 0.390 | -0.064 |
| tau | 0:56-59 | 0.438 | 0.705 | 0.559 / 0.865 | +0.416 |
| tau | 0:3-6 | 0.401 | 0.413 | 0.404 / 0.271 | +0.289 |
| tau | 0:58-61 | 0.368 | 0.514 | 0.490 / 0.643 | +0.343 |
| tau | 0:46-49 | 0.358 | 0.474 | 0.485 / 0.519 | -0.190 |
| tau | 0:60-63 | 0.315 | 0.651 | 0.742 / 0.667 | +0.405 |
| d_i_i2 | 0:47-49 | 0.410 | 0.442 | 0.273 / 0.455 | -0.498 |
| d_i_i2 | 0:12-14 | 0.262 | 0.212 | 0.260 / 0.169 | +0.171 |
| d_i_i2 | 0:46-48 | 0.248 | 0.369 | 0.236 / 0.225 | -0.155 |
| d_i_i2 | 0:61-63 | 0.220 | 0.458 | 0.506 / 0.418 | +0.160 |
| d_i_i2 | 0:63-65 | 0.218 | 0.452 | 0.515 / 0.497 | -0.237 |
| d_i_i3 | 0:46-49 | 0.421 | 0.480 | 0.293 / 0.509 | -1.250 |
| d_i_i3 | 0:47-50 | 0.340 | 0.409 | 0.408 / 0.359 | -0.389 |
| d_i_i3 | 0:48-51 | 0.316 | 0.434 | 0.407 / 0.368 | +0.565 |
| d_i_i3 | 0:55-58 | 0.310 | 0.652 | 0.543 / 0.776 | -0.407 |
| d_i_i3 | 0:58-61 | 0.310 | 0.569 | 0.572 / 0.615 | -0.829 |

## Temporal diagnostics (median over member-windows with defined values; descriptive)

Autocorrelation C(k): linear channels use the centered series; the circular channel uses centered unit vectors, C(k) = Re Σ conj(u_t) u_{t+k} / Σ |u_t|², which does not depend on the angular origin. The decorrelation lag is the first lag with C(k) < 1/e.

| Channel | Group | Lag-1 autocorrelation | Decorrelation lag (frames) | Not reached within available lags | Unavailable (reason: count) |
|---|---|---|---|---|---|
| kappa | 320K | 0.181 | 1.0 | 0 | none |
| kappa | 348K | 0.138 | 1.0 | 0 | none |
| tau | 320K | 0.294 | 1.0 | 0 | none |
| tau | 348K | 0.228 | 1.0 | 0 | none |
| d_i_i2 | 320K | 0.171 | 1.0 | 0 | none |
| d_i_i2 | 348K | 0.120 | 1.0 | 0 | none |
| d_i_i3 | 320K | 0.207 | 1.0 | 0 | none |
| d_i_i3 | 348K | 0.147 | 1.0 | 0 | none |

Distribution rows with an unavailable summary: 0 (reasons in the `status` column). Single-member groups have no replica variability (shown as n/a).
