# WP-D Held-out Protocol (mdCATH): Amendment 3 (proposed final, for review before scoring)

**Status:** PROPOSED. It is bound by hash in `benchmarks/audit/wpd_local_retrieval/results/wpd_heldout_amendment3.json`.
- It was recorded after the authorized download of the 30 primary files and after structural eligibility.
- It was recorded **before any held-out target, feature, ranking or score was computed**.
- No held-out retrieval output exists.

**Preserved unchanged** (hashes re-verified in the record):
- the original protocol, Amendments 1 and 2;
- the selection (primary and reserve order);
- the frozen code files `mdcath_reader.py`, `heldout_dependence.py`, `heldout_decision.py`, `mdcath_schema_check.py`;
- every earlier result file.

**Unchanged substance:**
- method A (exact LE-1 definition) and the baselines;
- separate targets (T_φψ primary, T_pp reported separately);
- the Amendment 1 estimand;
- domain-level resampling, B = 2,000, seed 20261008, ρ = 0.10.

## §1. Mandatory vs descriptive checks (new)

| Class | Checks | Consequence |
|---|---|---|
| **Mandatory** | File size and SHA-256 vs the selection (`mdcath_download.py`, `heldout_runner.py`); binding extension identity (WP-A SHA-256); binding τ/κ vs NumPy reference ≤ 1e-9 on identical inputs (§3); reader invariants (unit attribute, `numFrames` = coords length, unique atom indices, diagonal box); target-oracle regret failure (AdK §10 A3) | Run stops; **exit status 2** |
| **Authorization** | Held-out scoring requested without `results/wpd_heldout_scoring_authorization.json` | Refused; **exit status 3** |
| **Eligibility** | Original §4.3 criteria (with the §5 correction) | Ledger entry; the domain is not scored. A replacement needs a documented failure and separate authorization; none is automatic |
| **Descriptive** | Frame-order plausibility, memory, O availability, dependence diagnostics, method-level regret failures (counted), cis-ω counts | Recorded only |

`mdcath_schema_check.py` (Amendment 2) is **descriptive only**: its verdict never gated anything and its exit status
is always 0. Mandatory gating is in `heldout_runner.py` and `mdcath_download.py`.

## §2. Optional backbone oxygen (implements Amendment 2 §3)

- N, CA and C are mandatory. O is optional per residue: `mdcath_reader_r1.backbone_r1` returns O = NaN and
  `o_present = False` for a residue without a unique O.
- A missing O never alters N/CA/C data, the Cα methods or T_φψ.
- A T_pp window (peptide atoms C_r, O_r, N_{r+1} for r = s … s+W−2) is defined only if O is present for all its
  r in **both** replicas. Otherwise it is recorded as `peptide-plane O unavailable` for T_pp alone.
- Tested: `test_optional_oxygen_affects_only_peptide_plane_windows`. With one O removed, T_φψ domain results are
  identical, and exactly the T_pp windows containing that residue are undefined.

## §3. Coordinate precision and the Rust/NumPy reference check (clarifies original §5)

Order of operations per trajectory:
1. Read the stored float32 values into float64.
2. Make whole along N → CA → C → N(+1) by sequential minimum image; O follows its own C.
3. Discard the first ⌊0.10 n⌋ frames.
4. Round the whole Cα coordinates to float32.

**The same float32-rounded values** (held in float64) feed every Cα method, the binding (as float32) and the NumPy
reference. The reference check therefore compares identical inputs at ≤ 1e-9 rad. The φ/ψ and T_pp targets use the
unrounded whole-made backbone (original §5).

Tested on an actual boundary-crossing SYNTHETIC fixture, in which every atom is wrapped into a 30–32 Å box:
- stored C–N distances exceed 10 Å;
- reconstruction restores the chain to within 1e-4 Å of the source, up to one lattice vector;
- the reference check passes;
- a 1e-6 rad perturbation of the binding output fails it.

## §4. Dependence diagnostics, revision 1 (replaces Amendment 2 §4 wording and handling)

**Corrected claims.**
- The τ̂ ≤ 0.015·n cutoff is a **limited heuristic from an AR(1) simulation** of the WP-B Sokal estimator
  (`results/SYNTHETIC_tau_int_calibration.json`). It is not a reliability guarantee.
- A τ̂ below the cutoff means only "not flagged by the heuristic". Drift or multi-timescale processes can yield a
  small τ̂, and AR(1) underestimation already begins below the cutoff for short series.
- Amendment 2's statements that the cutoff makes τ̂ "resolvable" or "calibrated" are withdrawn.

**Handling (`heldout_dependence_r1.classify_series`):**
- **Non-finite observations are not removed by compressing time.** The estimator uses only the longest contiguous
  finite segment, and the number of gaps and unavailable values is reported.
- **Declared identities:** a finite non-zero value is an `identity_violation` error. Non-finite values are counted
  as unavailable observations and are not a violation.

**Classes:**
- `structural_identity`
- `constant_observed`
- `insufficient_observations` (< 50 contiguous finite values)
- `estimation_failure`
- `above_heuristic_cutoff`
- `estimated_below_heuristic_cutoff`

**Other rules:**
- Runner-up rule (Amendment 2 §4) unchanged.
- Diagnostics are computed for every pair of best-set candidates, per domain × temperature × W × target.
- They remain descriptive and do not alter the domain bootstrap.

## §5. Terminal caps and reader memory (prospective corrections)

**Caps (corrects the wording of original §4.3).**
- In mdCATH, caps are CHARMM patch atoms inside the first and last residues (verified: CAY/CY/OY on residue 1;
  NT/CAT on the last residue). They are not separate ACE/NME residues.
- The eligibility rule "residues with exactly one N, CA and C, one consecutive run of ≥ 50" is unchanged. It
  selects every residue of the domain, terminals included.
- The sentence "caps lack CA and are therefore not window residues" is withdrawn.

**Reader memory (corrects the Amendment 2 schema-check note and an earlier chat statement).**
- `backbone_r1` returns float64 arrays linear in frames × residues of **one** trajectory: F × L × 3 × 3 plus
  F × L × 3.
- Trajectories are read one at a time.
- Chunked reads over the selected atoms bound only the transient read buffer (chunk × 4L × 3 float32), not the
  output.
- The runner holds two trajectories (query and database) and their derived arrays per temperature, plus per-window
  distance matrices of size n_q × n_d.
- Measured tracemalloc peaks are in `results/wpd_mdcath_schema_check.json`.

## §6. Runner specification (`heldout_runner.py`)

**Per domain** (only if eligible), for each temperature in {320, 348}:
1. Query = replica 0 and database = replica 1, after burn-in (§3).
2. Geometry and reference check (mandatory).
3. Per-(window, frame) exclusions: degenerate angles (sin < 1e-3), bonds > 2.0 Å, binding failures.
4. Windows W ∈ {5, 8}.
5. Methods M1, M2, B1, B2, B3, A and C1 (C1: R < 0.1 → +∞).
6. Targets T_φψ and T_pp (§2).
7. Tie-aware NR@10, NR@1, overlap and regret (AdK §6 / §10 A3).
8. Spearman with the original §6.1 seeds.
9. Per-window means, then the domain median over windows. The domain is defined in a cell only if every non-C1
   method has ≥ 10 defined windows.
10. C1 coverage-matched medians (Amendment 1 §3), with ≥ 10 matched windows.
11. Dependence diagnostics (§4).

**Cell stage.**
- `heldout_decision_a2.analyse_cell`: Amendment 1 estimand, candidate set including A, domain bootstrap.
- Run per temperature × W, separately for each target.
- Fewer than 24 defined domains in a cell → insufficient evidence.

**Authorization guard.** `score-heldout` refuses (exit 3) without an authorization record binding the protocol
documents by hash. It was not run in this task.

## §7. Verification evidence

- **Tests** (SYNTHETIC): `test_heldout_runner.py`, plus the earlier held-out test files. Counts and the command
  are in the record.
- **Development-sample validation** (`results/wpd_heldout_runner_dev_validation.json`, domain 1ha8A00, excluded
  superfamily):
  - eligibility passed;
  - reference check passed;
  - every cell was defined.
  - Its method values are development outputs and were not used to change anything.

## §8. Download manifest and eligibility ledger

- **Manifest:** `results/wpd_mdcath_download_manifest.json` lists the 30 primary files with URL, size, SHA-256,
  verification and times. No reserve was downloaded.
- **Ledger:** `results/wpd_mdcath_eligibility_ledger.json` records each file's integrity and its §4.3 criteria
  (structure only). Ineligible domains are listed with reasons. No replacement was made.

## §9. Effective specification (read in this order)

1. `WPD_LOCAL_RETRIEVAL_PROTOCOL.md` (AdK; target, methods M1–C1, metrics, tie rule, degeneracy rules)
2. `WPD_HELDOUT_MDCATH_PROTOCOL.md` (data, selection, trajectory handling, aggregation, decision)
3. `WPD_HELDOUT_MDCATH_AMENDMENT_1.md` (estimand, invalid denominators, C1 coverage, labels)
4. `WPD_HELDOUT_MDCATH_AMENDMENT_2.md` (method A, T_pp, dependence; §4 superseded here)
5. This amendment.

Scoring requires review of this specification and a separate authorization record.
