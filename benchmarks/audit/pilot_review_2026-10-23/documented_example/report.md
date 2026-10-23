# Descriptive ensemble comparison: Pilot example: mdCATH 1ha8A00, 320 K vs 348 K, 3 replicas each

Generated from the tables in `tables/`. **Descriptive only:** no p-values, significance labels or uncertainty intervals. Differences between groups can reflect geometry, dynamics, sampling or starting conditions; replica variability and temporal diagnostics are shown so that readers can judge this. Marginal per-window distributions can miss joint changes.

**Pooling:** each group's pooled distribution concatenates its members' retained frames, so members are weighted by their frame counts (`member_frames` in `tables/distribution.csv`).

**Circular summaries:** the resultant length R is reported; mean directions (and differences involving them) are shown as unavailable when R < 1e-06, a numerical conditioning rule, not a biological criterion.

## Inputs

- **320K**: 3 member(s): 320K_0 (387 frames, source frames 43–429, frame identity not_verifiable); 320K_1 (387 frames, source frames 43–429, frame identity not_verifiable); 320K_2 (387 frames, source frames 43–429, frame identity not_verifiable)
- **348K**: 3 member(s): 348K_0 (387 frames, source frames 43–429, frame identity not_verifiable); 348K_1 (387 frames, source frames 43–429, frame identity not_verifiable); 348K_2 (387 frames, source frames 43–429, frame identity not_verifiable)
- Core: native extension resolved (SHA-256 d90e7ed4f1d8036a…); analysis options {'bins': 36, 'max_lag': 50, 'pbc': 'none'} (per-member reading options, e.g. frame slice and PBC policy, are in `manifest.json`)

## Exclusions

0 channel–window combinations excluded (invalid frames in at least one member); listed in `tables/exclusions.csv`.

## 320K vs 348K: largest pooled JS distances per channel (descriptive)

| Channel | Window | Pooled JS | Between-member JS (median) | Within-member JS median 320K / 348K | Center difference |
|---|---|---|---|---|---|
| kappa | 0:22-24 | 0.333 | 0.190 | 0.266 / 0.623 | +0.100 |
| kappa | 0:10-12 | 0.329 | 0.384 | 0.318 / 0.332 | +0.176 |
| kappa | 0:16-18 | 0.304 | 0.325 | 0.273 / 0.239 | +0.108 |
| kappa | 0:4-6 | 0.276 | 0.268 | 0.117 / 0.406 | +0.101 |
| kappa | 0:9-11 | 0.273 | 0.239 | 0.395 / 0.150 | -0.120 |
| tau | 0:9-12 | 0.471 | 0.540 | 0.425 / 0.353 | +0.719 |
| tau | 0:11-14 | 0.388 | 0.427 | 0.242 / 0.443 | -0.109 |
| tau | 0:23-26 | 0.384 | 0.520 | 0.446 / 0.614 | -0.287 |
| tau | 0:12-15 | 0.379 | 0.501 | 0.610 / 0.378 | +0.366 |
| tau | 0:1-4 | 0.369 | 0.424 | 0.404 / 0.349 | -0.987 |
| d_i_i2 | 0:22-24 | 0.337 | 0.243 | 0.301 / 0.631 | -0.225 |
| d_i_i2 | 0:10-12 | 0.325 | 0.373 | 0.322 / 0.357 | -0.370 |
| d_i_i2 | 0:16-18 | 0.308 | 0.330 | 0.287 / 0.246 | -0.248 |
| d_i_i2 | 0:9-11 | 0.282 | 0.250 | 0.395 / 0.164 | +0.248 |
| d_i_i2 | 0:4-6 | 0.276 | 0.274 | 0.131 / 0.420 | -0.213 |
| d_i_i3 | 0:23-26 | 0.352 | 0.447 | 0.426 / 0.406 | -0.611 |
| d_i_i3 | 0:24-27 | 0.351 | 0.290 | 0.167 / 0.547 | -0.659 |
| d_i_i3 | 0:9-12 | 0.343 | 0.418 | 0.325 / 0.366 | -0.238 |
| d_i_i3 | 0:12-15 | 0.341 | 0.438 | 0.580 / 0.364 | +0.737 |
| d_i_i3 | 0:19-22 | 0.298 | 0.327 | 0.533 / 0.379 | +0.399 |

## Temporal diagnostics (median over member-windows with defined values; descriptive)

Autocorrelation C(k): linear channels use the centered series; the circular channel uses centered unit vectors, C(k) = Re Σ conj(u_t) u_{t+k} / Σ |u_t|², which does not depend on the angular origin. The decorrelation lag is the first lag with C(k) < 1/e.

| Channel | Group | Lag-1 autocorrelation | Decorrelation lag (frames) | Not reached within available lags | Unavailable (reason: count) |
|---|---|---|---|---|---|
| kappa | 320K | 0.409 | 2.0 | 6 | none |
| kappa | 348K | 0.391 | 1.0 | 5 | none |
| tau | 320K | 0.698 | 10.0 | 19 | none |
| tau | 348K | 0.643 | 7.0 | 3 | none |
| d_i_i2 | 320K | 0.396 | 2.0 | 6 | none |
| d_i_i2 | 348K | 0.394 | 1.0 | 6 | none |
| d_i_i3 | 320K | 0.624 | 6.0 | 19 | none |
| d_i_i3 | 348K | 0.635 | 7.0 | 7 | none |

Distribution rows with an unavailable summary: 0 (reasons in the `status` column). Single-member groups have no replica variability (shown as n/a).
