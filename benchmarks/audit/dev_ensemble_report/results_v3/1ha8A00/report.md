# Descriptive ensemble comparison: mdCATH 1ha8A00, 320 K vs 348 K (development data; exposed domain)

Generated from the tables in `tables/`. **Descriptive only:** no p-values, significance labels or uncertainty intervals. Differences between groups can reflect geometry, dynamics, sampling or starting conditions; replica variability and temporal diagnostics are shown so that readers can judge this. Marginal per-window distributions can miss joint changes.

**Pooling:** each group's pooled distribution concatenates its members' retained frames, so members are weighted by their frame counts (`member_frames` in `tables/distribution.csv`).

**Circular summaries:** the resultant length R is reported; mean directions (and differences involving them) are shown as unavailable when R < 1e-06, a numerical conditioning rule, not a biological criterion.

## Inputs

- **320K**: 5 member(s): 320K_0 (387 frames, source frames 43–429, frame identity not_verifiable); 320K_1 (387 frames, source frames 43–429, frame identity not_verifiable); 320K_2 (387 frames, source frames 43–429, frame identity not_verifiable); 320K_3 (387 frames, source frames 43–429, frame identity not_verifiable); 320K_4 (387 frames, source frames 43–429, frame identity not_verifiable)
- **348K**: 5 member(s): 348K_0 (387 frames, source frames 43–429, frame identity not_verifiable); 348K_1 (387 frames, source frames 43–429, frame identity not_verifiable); 348K_2 (387 frames, source frames 43–429, frame identity not_verifiable); 348K_3 (387 frames, source frames 43–429, frame identity not_verifiable); 348K_4 (387 frames, source frames 43–429, frame identity not_verifiable)
- Core: native extension resolved (SHA-256 d90e7ed4f1d8036a…); analysis options {'bins': 36, 'max_lag': 50, 'pbc': 'none'} (per-member reading options, e.g. frame slice and PBC policy, are in `manifest.json`)

## Exclusions

0 channel–window combinations excluded (invalid frames in at least one member). `tables/exclusions.csv` has one row per member and cause, with location and source frames.

## 320K vs 348K: largest pooled JS distances per channel (descriptive)

| Channel | Window | Pooled JS | Between-member JS (median) | Within-member JS median 320K / 348K | Center difference |
|---|---|---|---|---|---|
| kappa | 0:22-24 | 0.259 | 0.181 | 0.165 / 0.361 | +0.042 |
| kappa | 0:33-35 | 0.239 | 0.323 | 0.274 / 0.257 | -0.131 |
| kappa | 0:31-33 | 0.220 | 0.190 | 0.340 / 0.122 | -0.064 |
| kappa | 0:27-29 | 0.210 | 0.370 | 0.418 / 0.293 | -0.098 |
| kappa | 0:30-32 | 0.206 | 0.143 | 0.151 / 0.095 | -0.053 |
| tau | 0:32-35 | 0.350 | 0.614 | 0.564 / 0.393 | +1.632 |
| tau | 0:34-37 | 0.298 | 0.498 | 0.402 / 0.511 | -0.545 |
| tau | 0:30-33 | 0.288 | 0.254 | 0.252 / 0.157 | +0.029 |
| tau | 0:33-36 | 0.287 | 0.534 | 0.453 / 0.476 | +0.995 |
| tau | 0:31-34 | 0.286 | 0.420 | 0.406 / 0.235 | -0.233 |
| d_i_i2 | 0:22-24 | 0.261 | 0.179 | 0.157 / 0.385 | -0.088 |
| d_i_i2 | 0:33-35 | 0.236 | 0.318 | 0.269 / 0.249 | +0.272 |
| d_i_i2 | 0:31-33 | 0.221 | 0.226 | 0.345 / 0.146 | +0.153 |
| d_i_i2 | 0:27-29 | 0.219 | 0.381 | 0.411 / 0.292 | +0.229 |
| d_i_i2 | 0:28-30 | 0.205 | 0.358 | 0.357 / 0.307 | +0.245 |
| d_i_i3 | 0:30-33 | 0.292 | 0.328 | 0.402 / 0.212 | +0.610 |
| d_i_i3 | 0:32-35 | 0.287 | 0.427 | 0.412 / 0.393 | +0.345 |
| d_i_i3 | 0:27-30 | 0.284 | 0.466 | 0.532 / 0.362 | +1.001 |
| d_i_i3 | 0:29-32 | 0.256 | 0.371 | 0.339 / 0.297 | +0.543 |
| d_i_i3 | 0:47-50 | 0.246 | 0.477 | 0.470 / 0.540 | +0.223 |

## Temporal diagnostics (median over member-windows with defined values; descriptive)

Autocorrelation C(k): linear channels use the centered series; the circular channel uses centered unit vectors, C(k) = Re Σ conj(u_t) u_{t+k} / Σ |u_t|², which does not depend on the angular origin. The decorrelation lag is the first lag with C(k) < 1/e.

| Channel | Group | Lag-1 autocorrelation | Decorrelation lag (frames) | Not reached within available lags | Unavailable (reason: count) |
|---|---|---|---|---|---|
| kappa | 320K | 0.393 | 2.0 | 11 | none |
| kappa | 348K | 0.339 | 1.0 | 6 | none |
| tau | 320K | 0.699 | 10.5 | 24 | none |
| tau | 348K | 0.633 | 7.0 | 5 | none |
| d_i_i2 | 320K | 0.395 | 2.0 | 10 | none |
| d_i_i2 | 348K | 0.340 | 1.0 | 7 | none |
| d_i_i3 | 320K | 0.630 | 7.0 | 27 | none |
| d_i_i3 | 348K | 0.593 | 6.0 | 13 | none |

Distribution rows with an unavailable summary: 0 (reasons in the `status` column). Single-member groups have no replica variability (shown as n/a).
