# structpair 0.2.0: explicitly selected multiple chains (specification, 2026-10-25)

Written **before** the implementation. It fixes the aggregation, the advisory comparison, ties, missing values and the output
schema. The implementation and its tests follow this text; a later change to the rules is a dated amendment of this file.

## 1. What the feature is, and what it is not

For two structures A and B the user names every chain of each entry that should take part. For every mapped residue the
command reports the **between-entry difference** (over the chain pairings A × B) and the **within-entry crystal-copy
variation** (over the chain pairs inside A and inside B).

**Mandatory statement (printed in the report, the summary, the manifest and the user guide):**
*Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of them. The
within-entry crystal-copy variation is a descriptive range of differences among the chains the user selected. It is not a
calibrated noise floor, not a statistical uncertainty, not a nucleotide or ligand effect and not an HDX replicate
estimate.* The advisory flag of section 6 is descriptive too: it is not a gate, not an adequacy verdict and not a
significance claim.

## 2. Spec changes (all explicit, nothing inferred)

- `structures.X.chain` (existing) names one chain. **Or** `structures.X.chains` lists every chain of entry X that takes
  part, together with `structures.X.primary_chain`, which must be one of them. `chain` and `chains` are exclusive; names must
  be unique; a named chain absent from the file is refused by the adapter. With more than one chain, `primary_chain` is
  required (no chain is chosen for the user).
- The correspondence applies **identically to every selected chain of an entry**: the same author residue numbers and
  insertion codes name the same position in each of its chains. A residue a chain lacks follows `on_missing`; alternate
  locations follow the entry's `altloc` policy per chain; residue names must agree between the chains of an entry
  (otherwise refused, or flagged when `resname_mismatch` is `flag`).
- Optional block `chain_variation`: `{"displacement_trim_fraction": 0.2, "advisory": {"statistic": "min", "within_factor": 1.0,
  "min_difference": {}}}`. Defaults are those shown.
- The primary pair (primary chain of A against primary chain of B) is the pair of the existing command. **Every legacy
  output file of a multi-chain run is the output of that primary pair**, and is identical to the output of a single-chain
  run with those two chains.

## 3. Channels and the value of one pair at one position

Seven channels. Window channels exist at the position of the window's centre, as today.

| Channel | Unit | Definition (existing component) |
|---|---|---|
| `disp` | Å | Cα displacement at the position after a trimmed fit of the two chains on the positions present in both: Kabsch (proper rotation, float64) on all of them, drop the floor(`displacement_trim_fraction` × n) most displaced (ties: lower position first), refit once; the displacement of every position under the refit. The K1 endpoint calculation (`benchmarks/audit/kras_k1/k1_endpoint.py`, frozen) |
| `d_kappa`, `d_i_i2`, `rmsd3`, `d_tau`, `d_i_i3`, `rmsd4` | rad, Å | Exactly the existing metrics (`window_metrics`, core κ/τ, MDAnalysis QCP), computed for the chain pair |

A pair value is **defined** only where both chains have the data the channel needs (existing exclusion rules for missing
residues, breaks and invalid windows; for `disp`, a finite Cα in both chains and at least 4 positions in the fit). Otherwise
it is undefined and is left out. An undefined value is never replaced and never counted as zero.

## 4. Aggregation

- **Pairs.** Cross pairings: every (chain of A, chain of B), `n_A × n_B`. Within A: every unordered pair of A's chains,
  `n_A (n_A − 1) / 2`; likewise within B. A pair with a single chain per entry contributes no within pair.
- **Statistics over the defined pair values at a position and channel:** count `n`, minimum, median, mean, maximum. The
  median of an even count is the mean of the two middle values. The mean is arithmetic on float64. All pair values are
  non-negative magnitudes (|Δκ|, circular |Δτ|, distances, RMSD, displacement), so the arithmetic mean never averages
  angles.
- **Within-entry values:** per entry, `n`, minimum and maximum over its within pairs. `within_max` is the maximum over the
  maxima of the entries that have a defined value. `within_scope` records which entries contributed: `both`, `A_only`,
  `B_only` or `none`.
- **Missing.** Each statistic uses only the defined values and states its `n`. With `n = 0` the statistic is empty (written
  as an empty cell), never zero.

## 5. Ties and equal values

- Rankings (top-n tables) keep the existing rule: value descending, then position ascending.
- Median of equal values is that value; min = max is reported as it is.
- Displacement trimming: equal displacements drop the lower position first.

## 6. Advisory comparison (descriptive)

At a position and channel, with `S` the chosen between-entry statistic (default `min`), `f` = `within_factor` (default 1.0)
and `floor` = `min_difference[channel]` (default 0):

`separated` = `S >= max(floor, f × within_max)`, inclusive within 1e-9 of the threshold.

- It is **not assessed** (empty, with a stated reason) when the between count is 0 or when `within_scope` is `none`. It is
  never written as false in those cases.
- `at_boundary` is true when `S` is within 1e-9 of the threshold.
- Default meaning: every defined cross pairing is at least as large as the largest within-entry value, the description
  KRAS-STRUCT-2 used for residues 33–36. **K1 settings**: statistic `mean`, `within_factor` 2.0, `min_difference` `{"disp": 1.0}`
  on the displacement channel reproduces the K1 endpoint-positive rule; it is a regression setting, not a default and not
  a verdict. K1's adequacy gate (number and spread of positive residues) is **not** part of the command.

## 7. Outputs and schema versions

- `STRUCTPAIR_VERSION` 0.2.0; manifest schema `topofold_inputs.structpair/2` (additive; results of schema `/1` stay valid and
  are read by the viewer as before).
- `summary.json` gains `within_entry_variation`: `status` (`available` if both entries have at least 2 chains, `partial` if
  one, `unavailable` if none), chains, pair counts, the advisory parameters, counts of separated and not-assessed
  positions per channel, and the mandatory statement. A single-chain run reports `unavailable` with the reason `one chain per
  entry`, and **its `windows.csv`, `residues.csv` and `exclusions.csv` are byte-identical to version 0.1.0**.
- A run with at least one entry of 2 or more chains also writes `chain_variation.csv` (long format, one row per position and
  channel), `chains.csv` (every selected chain with present, missing and altloc counts) and `chain_variation.png`, and adds
  a prominent section to `report.md`. Exclusions and provenance sections are unchanged and stay in the report.

## 8. Verification plan (engineering consistency, not biological validity)

1. SYNTHETIC unit tests for each rule above (aggregation, ties, missing, boundary, spec refusals).
2. A SYNTHETIC golden test: single-chain outputs byte-identical to fixtures generated by version 0.1.0 before the change.
3. The four committed single-chain demonstrations re-run with 0.2.0: `windows.csv`, `residues.csv`, `exclusions.csv`
   identical to the committed outputs.
4. KRAS-STRUCT-2 chains (5US4 A, B; 5USJ A–F): the window-metric aggregates (cross min, median, max over 12 pairings; copy
   maximum over 16 pairs) equal the frozen `KRAS-STRUCT-2_table.csv`.
5. K1 settings on the exact K1 chains and residues: `disp` between-mean, within-maximum, threshold and the flagged set equal
   `KRAS-K1_endpoint.json`.
6. Runtime and output size are reported.

## 9. Amendments (dated)

- **2026-10-25, after the first verification run:** the comparison with the frozen K1 record uses a declared tolerance of
  1e-4 Å (the comparisons that use the same reader as the frozen STRUCT-2 table keep 1e-9 Å). The cause was found before the
  tolerance was set: the frozen K1 code parses the PDB text in float64, while the adapter holds float32 positions
  (differing from the text parse by at most 3.8e-6 Å over the 8 chains; the K1 displacement values differ by at most
  6.6e-6 Å). No rule of sections 3–6 changed.

- **2026-10-25, amendment 2: advisory output semantics `advisory/2` (schema /3, structpair 0.3.0).** Sections 4 and 5 are unchanged;
  sections 1, 6 and 7 are amended as follows. Evidence for the change: `benchmarks/audit/structpair_advisory_correction_2026-10-25/`.
  - *Why.* In 0.2.0 the default minimum difference was 0 and the comparison was widened by 1e-9, so a zero between-entry value against
    a zero within-entry value (identical chains) was reported as `separated`, as was any difference below 1e-9; NaN and infinity were
    accepted as `within_factor` and as minimum differences. In the committed K1 regression record of 0.2.0 the flag was true at 78
    position-channel cells, 70 of them on channels for which no minimum difference had been configured.
  - *Rule.* The numeric between-entry and within-entry values are always reported. The flag is **assessed only for a channel with a
    configured minimum difference `m`**, a finite number above the numerical tolerance 1e-9 in that channel's unit. Then
    flag = `S >= max(m, within_factor × within_max)`, an exact float64 comparison (`S` the chosen between-entry statistic).
    **Equality at the positive boundary meets the rule** (the threshold is a declared "at least"); a value a hair below does not;
    `advisory_at_boundary` is true when |S − threshold| is within 1e-9, so a result that hinges on rounding is visible. It is left
    empty (never false) with a stated reason when no minimum difference is configured, when no cross pairing is defined, or when
    no within-entry pair is defined. It can never be true for values within the tolerance of zero, because the threshold exceeds it.
  - *Parsing.* Every number must be finite: `displacement_trim_fraction` in [0, 0.5), `within_factor` ≥ 0, each minimum difference
    above 1e-9 (NaN, infinity and non-positive values are refused with the key named).
  - *Names and versions.* The columns are `advisory_flag`, `advisory_at_boundary` and `advisory_not_assessed_reason` (were `separated`,
    `at_boundary`, `not_assessed_reason`); the summary counts are `advisory_configured`, `advisory_flag_true`, `advisory_flag_false`,
    `advisory_not_assessed`, `advisory_at_boundary`; the summary and manifest carry `advisory_semantics: advisory/2`; manifest schema
    `…/3`. Reports and figures call the result an advisory rule, never significance and never biological separation.
  - *Historical records.* Records of schema /2 keep the meaning `advisory/1` and are not rewritten; the viewer reads both and labels
    the historical one. Numeric columns are unchanged between the versions (verified on STRUCT-2 and K1).
