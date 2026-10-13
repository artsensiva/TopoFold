# WP-B Specification Amendment 1 (review findings; recorded before any G0/G1 execution)

**Disclosure:** at the time of this amendment, no trajectory-derived structural quantity had
been computed or inspected:
- G0 and G1 had never run;
- `results/` contained only the original freeze record;
- during the ongoing download only file sizes were observed.

The amendments are therefore pre-inspection. The original `SPEC.md`, `wpb_spec.json` and
`results/adk_spec_freeze.json` are unchanged. Amendment hashes are recorded in
`results/adk_spec_amendment1_freeze.json`. Machine-readable values are in
`wpb_spec_amendment1.json`. Where this amendment and SPEC.md differ, this amendment governs.

## A1. Enforceable reference mapping (replaces SPEC §3 "common residues" wording)

- **Declared sequence:** the RCSB canonical sequence of entity 1 of 4AKE and of 1AKE. It is
  identical for both, 214 residues, and recorded in the JSON.
- **Correspondence:** trajectory Cα residue *i* (resid 1–214) maps to reference residue
  `resSeq` *i*. **Every** residue 1–214 must be present in each reference chain used, with a
  residue name consistent with the declared sequence (CHARMM `HSD`/`HSE`/`HSP` ≡ `H`).
- No "common residue" subset is used. If a residue is missing or mismatched:
  - in a **primary** chain (4AKE:A, 1AKE:A), G0 stops;
  - in a **sensitivity** chain (B), that chain is marked invalid, reported, and not used for
    fitting.
- **Alternate locations:** for each residue's Cα, use the record with blank altloc. If there
  is none, use altloc `A`. If neither exists, the residue is a mapping failure. Every residue
  with altlocs is recorded, with the altloc used.

## A2. Alignment-core convergence (replaces SPEC §4 "until stable")

- The iteration is **converged** only if the residue set at iteration *k* equals the set at
  *k−1* within the 50-iteration limit. Reaching the limit without that is **not converged**.
- **Primary core** (chain A, 1.0 Å): must be converged and have ≥ 30 residues; otherwise G1
  stops.
- **Sensitivity cores** (A 0.5 Å, A 2.0 Å, B 1.0 Å): the same validity test. An invalid core is
  reported with its reason and is **not** used for any fit or classification.

## A3. Observed visits, transition evidence and outcome precedence (replaces SPEC §7–§8 outcome rules)

### Definitions

| Term | Definition |
|------|------------|
| **Run** | A maximal sequence of consecutive frames with the same proximity class (near open / near closed / ambiguous). A run never spans a cadence irregularity (A5) |
| **Sustained visit** | A run with elapsed span t_last − t_first ≥ 5.0 ns (A5) |
| **Brief run** | A shorter run; reported separately |
| **Observed visit** (descriptive) | A sustained visit near a reference under one analysis variant |
| **Transition evidence** (claimed) | A trajectory whose first 5.0 ns is classified near its start reference contains a sustained visit near the other reference, **in every valid variant**. Otherwise any cross visit is reported as an observed visit only, not as transition evidence |
| **Ambiguity dominance** | In a given trajectory and variant, ambiguous frames are > 50 % of all frames |
| **Variant** | Proximity fractions {0.20, 0.25, 0.33} with the primary core and chain A; the primary fraction with each valid sensitivity core; the primary fraction with chain-B references if valid. The primary variant is fraction 0.25, chain A, core 1.0 Å |

### Per-variant category (first matching rule applies)

1. Ambiguity dominance in **either** trajectory → *unresolved*.
2. Sustained visits near **both** references within one trajectory → *suitable*.
3. Neither trajectory has a sustained visit near its own start reference → *unsuitable*.
4. No sustained visit near the other reference in either trajectory → *restricted*.
5. Otherwise → *unresolved*.

### Final outcome

- If every valid variant gives the same category, that category is the outcome.
- Otherwise the outcome is **unresolved pending specified evidence**, with the disagreement
  reported.
- *Suitable* therefore requires agreement across all variants. Ambiguity and disagreement take
  precedence over any visit.
- G0 failure means G1 is not run, and the outcome is *unsuitable*, pending repaired inputs
  (none may be substituted).

## A4. Numerical cross-check with a predeclared tolerance

Kabsch RMSD (SVD) is compared with MDAnalysis `rms.rmsd(..., center=True, superposition=True)`
(QCP algorithm) on core-fit RMSD to each primary reference. The comparison covers every 50th
frame of each trajectory and the crystal pair.

**Tolerance: |difference| ≤ 1e-6 Å.** Both use float64 on coordinates of at most ~100 Å, so
the expected discrepancy is round-off (≲ 1e-10 Å). 1e-6 Å is still more than three orders of
magnitude below the 0.005 Å xtc quantization. Exceeding it stops G1.

Focused tests (`test_g1.py`) cover:
- alignment: recovered rigid motion; reflection not allowed;
- periodic reconstruction: a split chain is made whole; a whole chain is unchanged;
- core convergence: converged, non-converged and too-small cases;
- classification edges: ties, boundary radius, 4.8 ns vs ≥ 5.0 ns runs, cadence gaps;
- outcome precedence.

## A5. Sustained duration from timestamps (replaces SPEC §7 "≥ 25 consecutive frames")

- Sustained means an elapsed span t_last − t_first ≥ 5.0 ns. SPEC.md's "25 consecutive frames
  (5 ns)" was internally inconsistent: 25 samples at 200 ps span 4.8 ns. The stated 5 ns intent
  governs; at 200 ps cadence it requires ≥ 26 frames.
- Observation counts (frames) and elapsed spans (ns) are both reported.
- **Nominal cadence:** the single Δt verified by G0. Any frame interval differing from it by
  more than 1e-3 ps is a cadence irregularity: runs are split there, and the irregularities are
  reported. A trajectory with more than one distinct Δt is reported; it is not resampled.

## A6. Binding G1 to G0

Before computing anything, G1 recomputes the SHA-256 of:
- `SPEC.md`, `wpb_spec.json`, `AMENDMENTS.md` and `wpb_spec_amendment1.json`;
- the four raw data files and the two reference PDB files.

It compares these with the values recorded in `results/adk_g0_provenance.json`, which must
report PASS under Amendment 1. Any mismatch, a missing amendment hash in G0, or a G0 record
older than the amendment freeze is rejected as stale, and G1 stops.

# Amendment 2 (recorded before any G0/G1 execution; no trajectory-derived quantity inspected)

## A2 revised: two-phase core derivation

**Finding.** A synthetic test written for A4 showed the A2 iteration failing when a large
segment moves rigidly. A test with 20 of 60 residues displaced by 6 Å collapsed to an empty
core on the first step. The reason: the initial fit over all residues leaves every deviation
above the cutoff. Open and closed AdK differ by large rigid movements of non-core segments, so
the pre-registered algorithm could fail for a reason unrelated to the data's quality. The
finding came from synthetic data only.

**Revised algorithm** (crystal structures only, deterministic):
- **Phase 1, trim.** Fit on the current set. While any deviation in the set is ≥ the cutoff,
  drop the worst ⌈0.1·n⌉ residues (ties broken by residue order), for at most 200 steps.
- **Phase 2, refine.** From the trimmed set, run the A2 fixed-point iteration (all residues with
  deviation < cutoff) for at most 50 steps.
- **Converged** only if phase 1 met its condition and phase 2 reached a fixed point.

Validity rules are unchanged: convergence and ≥ 30 residues for the primary core; sensitivity
cores are reported invalid rather than used.

The test case is unchanged and must recover exactly the 40 rigid residues. The non-convergence
test now reaches its limit through the trimming step limit, because a single refine iteration can
legitimately be stable after trimming; the contract (limit reached → not converged → invalid) is
unchanged.
