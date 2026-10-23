# Descriptive ensemble comparison: Reviewer operational gap check: mdCATH 1s14B00; development data; exclusion handling only

Generated from the tables in `tables/`. **Descriptive only:** no p-values, significance labels or uncertainty intervals. Differences between groups can reflect geometry, dynamics, sampling or starting conditions; replica variability and temporal diagnostics are shown so that readers can judge this. Marginal per-window distributions can miss joint changes.

**Pooling:** each group's pooled distribution concatenates its members' retained frames, so members are weighted by their frame counts (`member_frames` in `tables/distribution.csv`).

**Circular summaries:** the resultant length R is reported; mean directions (and differences involving them) are shown as unavailable when R < 1e-06, a numerical conditioning rule, not a biological criterion.

## Inputs

- **320K**: 2 member(s): 320K_0 (40 frames, source frames 43–82, frame identity not_verifiable); 320K_1 (40 frames, source frames 43–82, frame identity not_verifiable)
- **348K**: 2 member(s): 348K_0 (40 frames, source frames 43–82, frame identity not_verifiable); 348K_1 (40 frames, source frames 43–82, frame identity not_verifiable)
- Core: native extension resolved (SHA-256 d90e7ed4f1d8036a…); analysis options {'bins': 36, 'max_lag': 50, 'pbc': 'none'} (per-member reading options, e.g. frame slice and PBC policy, are in `manifest.json`)

## Exclusions

10 channel–window combinations excluded (invalid frames in at least one member); listed in `tables/exclusions.csv`.

## 320K vs 348K: largest pooled JS distances per channel (descriptive)

| Channel | Window | Pooled JS | Between-member JS (median) | Within-member JS median 320K / 348K | Center difference |
|---|---|---|---|---|---|
| kappa | 0:2195-2197 | 0.581 | 0.642 | 0.393 / 0.880 | -0.164 |
| kappa | 0:2155-2157 | 0.570 | 0.602 | 1.000 / 0.303 | +0.343 |
| kappa | 0:2078-2080 | 0.533 | 0.596 | 0.156 / 0.605 | +0.170 |
| kappa | 0:2214-2216 | 0.530 | 0.613 | 0.352 / 0.737 | -0.228 |
| kappa | 0:2215-2217 | 0.513 | 0.586 | 0.432 / 0.634 | +0.154 |
| tau | 0:2213-2216 | 0.595 | 0.657 | 0.372 / 1.000 | -1.402 |
| tau | 0:2059-2062 | 0.593 | 0.654 | 0.316 / 1.000 | -0.745 |
| tau | 0:2195-2198 | 0.592 | 0.663 | 0.294 / 0.908 | -0.446 |
| tau | 0:2194-2197 | 0.577 | 0.615 | 0.197 / 0.939 | +0.382 |
| tau | 0:2062-2065 | 0.576 | 0.611 | 0.156 / 1.000 | +0.665 |
| d_i_i2 | 0:2195-2197 | 0.628 | 0.667 | 0.268 / 0.843 | +0.391 |
| d_i_i2 | 0:2155-2157 | 0.578 | 0.656 | 1.000 / 0.368 | -0.684 |
| d_i_i2 | 0:2079-2081 | 0.541 | 0.637 | 0.273 / 0.969 | -0.245 |
| d_i_i2 | 0:2214-2216 | 0.516 | 0.567 | 0.340 / 0.734 | +0.476 |
| d_i_i2 | 0:2078-2080 | 0.498 | 0.552 | 0.189 / 0.562 | -0.413 |
| d_i_i3 | 0:2196-2199 | 0.591 | 0.685 | 0.164 / 0.969 | +0.317 |
| d_i_i3 | 0:2213-2216 | 0.584 | 0.663 | 0.302 / 0.980 | +1.623 |
| d_i_i3 | 0:2062-2065 | 0.565 | 0.612 | 0.215 / 1.000 | +0.530 |
| d_i_i3 | 0:2059-2062 | 0.561 | 0.570 | 0.193 / 1.000 | +1.019 |
| d_i_i3 | 0:2154-2157 | 0.559 | 0.598 | 0.929 / 0.175 | -0.554 |

## Temporal diagnostics (median over member-windows with defined values; descriptive)

Autocorrelation C(k): linear channels use the centered series; the circular channel uses centered unit vectors, C(k) = Re Σ conj(u_t) u_{t+k} / Σ |u_t|², which does not depend on the angular origin. The decorrelation lag is the first lag with C(k) < 1/e.

| Channel | Group | Lag-1 autocorrelation | Decorrelation lag (frames) | Not reached within available lags | Unavailable (reason: count) |
|---|---|---|---|---|---|
| kappa | 320K | 0.059 | 1.0 | 0 | none |
| kappa | 348K | 0.044 | 1.0 | 0 | none |
| tau | 320K | 0.115 | 1.0 | 0 | none |
| tau | 348K | 0.100 | 1.0 | 0 | none |
| d_i_i2 | 320K | 0.055 | 1.0 | 0 | none |
| d_i_i2 | 348K | 0.040 | 1.0 | 0 | none |
| d_i_i3 | 320K | 0.073 | 1.0 | 0 | none |
| d_i_i3 | 348K | 0.078 | 1.0 | 0 | none |

Distribution rows with an unavailable summary: 0 (reasons in the `status` column). Single-member groups have no replica variability (shown as n/a).
