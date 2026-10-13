# WP-D Held-out Protocol (mdCATH): Amendment 2 (proposed)

**Status:** PROPOSED. It is bound by hash in `benchmarks/audit/wpd_local_retrieval/results/wpd_heldout_amendment2.json`.
It must be reviewed before any held-out scoring. Nothing held-out has been downloaded or scored, except the
schema-inspection sample, which is development data under the original protocol §3.1.

**Amends:** `docs/WPD_HELDOUT_MDCATH_PROTOCOL.md` (unchanged) and Amendment 1 (unchanged).

**Unchanged:**
- The domain selection: schema sample, the 30 primary domains, the 10 reserves and their order
  (`results/wpd_mdcath_selection.json`).
- Seeds, ρ = 0.10.
- The Amendment 1 estimand: L = (median NR_m − median NR_best) / (1 − median NR_best), with domain resampling,
  best re-selected per resample, invalid-denominator handling, labels and interpretation limits.

## §1. Add method A (exact, untuned)

**Definition.** For a window of W residues (start s):
- κ: the W − 2 production turning angles at interior vertices s+1 … s+W−2;
- τ: the W − 3 production Cα dihedrals s … s+W−4.

Both are from `topofold.compute_invariants` on the full-chain float32 frame (identity and reference checks as in
the original §5).

  D_A(a, b) = √( ½ · (1/(W−2)) Σ_v (κ_v^a − κ_v^b)² + ½ · (1/(W−3)) Σ_j wrap(τ_j^a − τ_j^b)² )

- wrap maps to (−π, π].
- Each channel is averaged over its own element count, then the channels are weighted equally (½ each).
  Equal weighting in rad² is a modeling choice.
- **No parameter of A was tuned.** Its definition was frozen in LE-1 (`docs/LEGACY_ENGINE_EVALUATION_PROTOCOL.md`
  §4) before any A outcome existed, and it is used here unchanged.
- Reference implementation: `benchmarks/audit/legacy_engine/legacy_lib.py::d_aligned`.

**Candidate set.**
- A joins the best-method candidates. The order becomes B2, B1, B3, M1, M2, A; it is used only to break exact
  ties.
- This changes the candidate set, **not the estimand**. If A is best, every other method's L is measured against
  A. That is the intended effect of adding a competitor.
- Implementation: `heldout_decision_a2.py` loads the unmodified Amendment 1 code and replaces only the candidate
  tuple.
- **Test:** with A not best, the results are identical to Amendment 1 (`test_heldout_dependence.py`).

**Not added:**
- Legacy κ–τ Fréchet F. It was deprioritized for this task as an engineering decision. Its LE-1 development
  outcome was formally inconclusive.
- No weighting or normalization variants of A.

**Development role disclosed.**
- A is added **because** it had the lowest mean NR@10 on the AdK φ/ψ target in LE-1. AdK is development data,
  and its superfamily (3.40.50.300) is already excluded from the held-out pool.
- That choice is made after seeing development results. It creates a selection effect that the held-out
  evaluation, not AdK, must test.
- On the AdK peptide-plane target, A was behind window Cα RMSD (B2) in all four cells.

## §2. Baselines retained

M1, M2, B1, B2 and B3 (= κ-only) remain candidates. C1 remains a coverage-matched comparator outside best
selection (Amendment 1 §3).

## §3. Targets kept separate

**Primary target: T_φψ** (original §3 of the WP-D AdK protocol). The overall outcomes of Amendment 1 §4 are
determined on T_φψ only.

**Separate target: T_pp** (new to the held-out tier).
- Definition: the optimal-superposition RMSD (Horn) of the peptide-plane atoms C_r, O_r and N_{r+1}, for
  r = i … i+W−2, with no Cα.
- It is computed and classified with exactly the same estimand, intervals and classes, but **reported
  separately** and never pooled with T_φψ.
- Its overall outcome is reported as its own result. It does not change the T_φψ outcome.
- Bias disclosure: T_φψ is angular and may favour angular methods; T_pp is Cartesian and may favour Cartesian
  superposition.

**Atom availability.**
- T_pp requires backbone O for every residue r in the window's peptide range. Window eligibility for T_pp adds
  that requirement.
- A window lacking an O atom is excluded **from T_pp only** and counted. The schema-sample verification below
  checks O naming, including the C-terminal residue.

## §4. Dependence diagnostics (descriptive; the inference unit stays the domain)

**Series.** For each domain × temperature × W × target and each reported contrast X − Y: the window-averaged
per-query difference of NR@10, in query-frame order. Each series is classified by
`heldout_dependence.classify_series`:

| Class | Rule |
|---|---|
| `structural_identity` | X and Y are declared identical for the configuration (in the held-out set, only X = Y). The series must be exactly zero; otherwise this is an error that stops the run |
| `constant_observed` | Not a declared identity, but the observed series is constant (range ≤ 1e-12·(1 + \|mean\|)); its value is reported |
| `estimation_failure` | One of: fewer than 50 finite values; estimator returned no value; or **τ̂ > 0.015·n** (`tau_not_resolvable_from_series_length`) |
| `estimated` | τ̂ reported in frames (ns) |

**Calibration.** The 0.015·n limit comes from a SYNTHETIC AR(1) calibration of the same Sokal estimator
(`results/SYNTHETIC_tau_int_calibration.json`). There, the median estimate fell below 0.8 × the true value
whenever τ̂/n ≥ 0.017. A larger τ̂ therefore cannot distinguish moderate from much longer dependence.

**Self-comparison.** When X is the best method, its diagnostic uses the **runner-up**, never X against itself.
This corrects the LE-1 gate-7 artefact.

**Role.** These diagnostics do not alter the domain-level bootstrap: domains are the independent unit, and
frames within a domain are never resampled. Per-domain values describe within-domain precision only. A high
share of `estimation_failure` is reported as a limitation, not hidden.

## §5. Reader (feasibility only)

A minimal reader `mdcath_reader.py` (h5py) is verified on the schema sample (`results/wpd_mdcath_schema_check.json`).
The held-out loader must use the conventions recorded there:
- atom names from `pdbProteinAtoms`;
- coordinates in Å; box in nm × 10;
- frame order as stored;
- explicit cap handling.

Any held-out-relevant convention that the schema check finds different from the original protocol §4.3/§5 is
listed in that file. Such a difference requires a further amendment before scoring.

## §6. Tests

`test_heldout_dependence.py` (SYNTHETIC) covers:
- the four diagnostic classes;
- identity violation;
- AR(1) τ recovery;
- runner-up selection;
- the candidate set (including A);
- equality with Amendment 1 when A is not best.
