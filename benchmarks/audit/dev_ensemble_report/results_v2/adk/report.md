# Descriptive ensemble comparison: AdK open-start vs closed-start, single runs, 0-335.2 ns every 4th frame (development data)

Generated from the tables in `tables/`. **Descriptive only:** no p-values, significance labels or uncertainty intervals. Differences between groups can reflect geometry, dynamics, sampling or starting conditions; replica variability and temporal diagnostics are shown so that readers can judge this. Marginal per-window distributions can miss joint changes.

**Pooling:** each group's pooled distribution concatenates its members' retained frames, so members are weighted by their frame counts (`member_frames` in `tables/distribution.csv`).

**Circular summaries:** the resultant length R is reported; mean directions (and differences involving them) are shown as unavailable when R < 1e-06, a numerical conditioning rule, not a biological criterion.

## Inputs

- **open_start**: 1 member(s): open_start_0 (420 frames, source frames 0–1676, frame identity not_verifiable)
- **closed_start**: 1 member(s): closed_start_0 (420 frames, source frames 0–1676, frame identity not_verifiable)
- Core: native extension resolved (SHA-256 d90e7ed4f1d8036a…); analysis options {'bins': 36, 'max_lag': 50, 'pbc': 'none'} (per-member reading options, e.g. frame slice and PBC policy, are in `manifest.json`)

## Exclusions

0 channel–window combinations excluded (invalid frames in at least one member); listed in `tables/exclusions.csv`.

## open_start vs closed_start: largest pooled JS distances per channel (descriptive)

| Channel | Window | Pooled JS | Between-member JS (median) | Within-member JS median open_start / closed_start | Center difference |
|---|---|---|---|---|---|
| kappa | SYSTEM:120-122 | 0.767 | 0.767 | n/a / n/a | +0.381 |
| kappa | SYSTEM:134-136 | 0.740 | 0.740 | n/a / n/a | -0.373 |
| kappa | SYSTEM:41-43 | 0.737 | 0.737 | n/a / n/a | -0.407 |
| kappa | SYSTEM:143-145 | 0.731 | 0.731 | n/a / n/a | -0.457 |
| kappa | SYSTEM:9-11 | 0.685 | 0.685 | n/a / n/a | +0.379 |
| tau | SYSTEM:121-124 | 0.953 | 0.953 | n/a / n/a | -1.443 |
| tau | SYSTEM:134-137 | 0.830 | 0.830 | n/a / n/a | -2.937 |
| tau | SYSTEM:119-122 | 0.821 | 0.821 | n/a / n/a | -0.622 |
| tau | SYSTEM:10-13 | 0.788 | 0.788 | n/a / n/a | -2.074 |
| tau | SYSTEM:143-146 | 0.773 | 0.773 | n/a / n/a | +1.945 |
| d_i_i2 | SYSTEM:120-122 | 0.772 | 0.772 | n/a / n/a | -0.877 |
| d_i_i2 | SYSTEM:134-136 | 0.750 | 0.750 | n/a / n/a | +0.808 |
| d_i_i2 | SYSTEM:41-43 | 0.731 | 0.731 | n/a / n/a | +0.865 |
| d_i_i2 | SYSTEM:143-145 | 0.729 | 0.729 | n/a / n/a | +0.944 |
| d_i_i2 | SYSTEM:9-11 | 0.682 | 0.682 | n/a / n/a | -0.826 |
| d_i_i3 | SYSTEM:121-124 | 0.892 | 0.892 | n/a / n/a | -1.356 |
| d_i_i3 | SYSTEM:119-122 | 0.825 | 0.825 | n/a / n/a | -2.412 |
| d_i_i3 | SYSTEM:134-137 | 0.777 | 0.777 | n/a / n/a | +2.928 |
| d_i_i3 | SYSTEM:142-145 | 0.770 | 0.770 | n/a / n/a | +2.800 |
| d_i_i3 | SYSTEM:42-45 | 0.750 | 0.750 | n/a / n/a | -1.457 |

## Temporal diagnostics (median over member-windows with defined values; descriptive)

Autocorrelation C(k): linear channels use the centered series; the circular channel uses centered unit vectors, C(k) = Re Σ conj(u_t) u_{t+k} / Σ |u_t|², which does not depend on the angular origin. The decorrelation lag is the first lag with C(k) < 1/e.

| Channel | Group | Lag-1 autocorrelation | Decorrelation lag (frames) | Not reached within available lags | Unavailable (reason: count) |
|---|---|---|---|---|---|
| kappa | open_start | 0.131 | 1.0 | 2 | none |
| kappa | closed_start | 0.144 | 1.0 | 0 | none |
| tau | open_start | 0.226 | 1.0 | 7 | none |
| tau | closed_start | 0.264 | 1.0 | 0 | none |
| d_i_i2 | open_start | 0.126 | 1.0 | 2 | none |
| d_i_i2 | closed_start | 0.129 | 1.0 | 0 | none |
| d_i_i3 | open_start | 0.169 | 1.0 | 5 | none |
| d_i_i3 | closed_start | 0.200 | 1.0 | 0 | none |

Distribution rows with an unavailable summary: 0 (reasons in the `status` column). Single-member groups have no replica variability (shown as n/a).
