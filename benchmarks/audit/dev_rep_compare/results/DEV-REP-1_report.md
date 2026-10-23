# DEV-REP-1: replicated descriptive comparison, mdCATH 1ha8A00, 320 K vs 348 K

Generated from `DEV-REP-1_results.json` (SHA-256 e7532407fcc292f8…). Development evidence only. Five replicas per temperature, all started from one equilibrated configuration, so they measure reproducibility given that start. Pairwise replica distances share replicas and are dependent; categories are descriptive, not tests. 320 K vs 348 K is a simulation-condition contrast with no experimental endpoint.

## Inputs

- Source SHA-256 e06d8cdb5d8b6fef… (matches the WP-D record). Source frames 43–429 (387 per replica; frames 0–42 excluded as predeclared). Time: INFERRED: index x 1 ns (mdCATH paper); the file stores no timestamps
- Checks: source_hash_matches True, frames_per_replica True, source_frame_indices_identical True, inferred_time_labelled True, residue_identity_identical True, mapping_check_no_mismatch True, no_breaks True, native_extension_resolved True, total_frames_within_budget True
- Core: wrapper SHA-256 49780093095f27c8…, native extension resolved (SHA-256 d90e7ed4f1d8036a…)
- Excluded windows per channel: kappa 2, d_i_i2 2, tau 2, d_i_i3 2 (chain-end windows where the feature is undefined)

## Categories (primary, 36 bins)

| Channel | Consistent | Replica-dependent | Inconclusive |
|---|---|---|---|
| kappa | 0 | 0 | 49 |
| d_i_i2 | 0 | 0 | 49 |
| tau | 0 | 0 | 48 |
| d_i_i3 | 0 | 0 | 48 |

*Consistent*: every between-temperature replica-pair distance exceeds every within-temperature one. *Replica-dependent*: that holds only after leaving out one specific replica. *Inconclusive*: neither; this is not evidence of no difference.

Paired channels on the same windows (category κ|d(i,i+2) and τ|d(i,i+3), primary):

- kappa~d_i_i2: inconclusive|inconclusive 49
- tau~d_i_i3: inconclusive|inconclusive 48

Spearman of pooled 320-vs-348 distances between paired channels: kappa~d_i_i2 0.969, tau~d_i_i3 0.598

κ and d(i,i+2) are near-deterministic functions of each other; d(i,i+3) cannot distinguish τ from −τ. Agreement or disagreement does not show that one representation is better.

## Windows by category (primary; pooled JS distance with its leave-one-replica-out range)


## Sensitivity (predeclared; not used for selection)

| Configuration | Channel | Same category as primary | Spearman of pooled distance vs primary |
|---|---|---|---|
| bins18 | kappa | 49/49 | 0.986 |
| bins18 | d_i_i2 | 49/49 | 0.978 |
| bins18 | tau | 48/48 | 0.989 |
| bins18 | d_i_i3 | 48/48 | 0.996 |
| bins72 | kappa | 49/49 | 0.988 |
| bins72 | d_i_i2 | 49/49 | 0.970 |
| bins72 | tau | 48/48 | 0.992 |
| bins72 | d_i_i3 | 48/48 | 0.996 |
| tau_origin_half_bin | kappa | 49/49 | 1.000 |
| tau_origin_half_bin | d_i_i2 | 49/49 | 1.000 |
| tau_origin_half_bin | tau | 48/48 | 0.992 |
| tau_origin_half_bin | d_i_i3 | 48/48 | 1.000 |

## Resources

- Total CPU 40.2 s (budget 1800 s); wall 26.9 s; max RSS: bridge stage 110336 KiB, features and comparison 55832 KiB
