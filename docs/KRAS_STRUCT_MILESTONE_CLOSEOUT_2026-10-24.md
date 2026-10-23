# KRAS static-structure milestone: closeout (2026-10-24)

**Status: closed.** No further KRAS case is started. The frozen artifacts of both studies are unchanged:
- KRAS-STRUCT-1: `benchmarks/audit/kras_struct_cmp/`, freeze `results/KRAS-STRUCT-1_freeze.json`;
- KRAS-STRUCT-2: `benchmarks/audit/kras_struct_g12d/`, freeze `results/KRAS-STRUCT-2_freeze.json`.

Every number below is read from those results files. Values marked "measured" come directly from code output.
Explanations are labelled as such.

## What was compared

- **KRAS-STRUCT-1:** 9O0N (GDP) vs 9O0O (GMPPNP), wild-type KRAS(1–169), both with MRTX1133. One chain each; mapped
  residues 1–168.
- **KRAS-STRUCT-2:** 5US4 (GDP) vs 5USJ (GMPPNP), KRAS G12D, no pharmacological inhibitor modelled. The GDP crystal
  grew with a compound present, and no density was reported for it. Mapped residues 1–59 and 68–166, in all 8 chains;
  primary chain pair A/A, plus 11 other pairings and 16 same-entry copy pairs.
- **Metrics (both studies, identical coordinates and coverage per window):** |Δκ|, circular |Δτ|, window Cα distances,
  and 3- and 4-residue proper-rotation local RMSD.

## Recurring local differences and chain sensitivity (measured)

- **Recurs in both studies:** windows centred around residues **30–38** are among the ten largest for every metric in
  both studies.
- **Chain-robust in STRUCT-2:** for every metric, at least three windows within residues 29–40 are among the ten
  largest in **12 of 12** chain pairings. At window 33–36 the smallest value over the pairings (|Δτ| 2.89 rad; 4-residue RMSD 2.20 Å) far exceeds
  the largest same-entry copy difference (0.27 rad; 0.19 Å).
- **Structure- or chain-dependent:**
  - 102–110 is prominent only in STRUCT-1. That study has one chain per entry, so chain dependence could not be
    assessed. In 5USJ, copy-to-copy differences around 100–125 match or exceed the GDP–GMPPNP differences.
  - 72–76 ranks high in 2–11 of 12 pairings, depending on the metric.
  - 50–53 and 57–59 are prominent only in STRUCT-2.
- **Whole-profile agreement between the studies** is modest: Spearman 0.28–0.63 per metric.
- **Metric agreement** (both studies): the three-residue metrics are close to redundant (Spearman 0.79–0.97), and the
  four-residue metrics agree less (0.38–0.82).
- **Alternate conformations:**
  - STRUCT-1: switching from blank-then-A to highest occupancy changes each metric by at most 0.023–0.097 (radians or Å).
  - STRUCT-2: neither entry has alternate locations.

## Excluded region

**Switch II (residues 60–67) is excluded from STRUCT-2.** It is unobserved in at least one of the eight chains:
- 60–64 in 5US4 A;
- 60–66 in 5US4 B;
- 61–63 in 5USJ A;
- 61–67 in 5USJ D.

The comparison therefore says nothing about Switch II for G12D. STRUCT-1 observed it (59–66 ranked high there), so no
cross-study statement about Switch II is possible.

## Limitations

- **Crystal conditions:**
  - STRUCT-1: P 2₁2₁2₁ vs P 2 3, different crystallization conditions.
  - STRUCT-2: H 3 (two copies) vs P 1 (six copies), pH 9.1 with additives vs pH 5.5.
  - Packing and conditions can change loops.
- **Coordinate quality:**
  - Resolutions: 1.40/1.90 Å (STRUCT-1) and 1.83/1.94 Å (STRUCT-2).
  - STRUCT-1 has alternate locations (7 residues in 9O0N, 1 in 9O0O), including residue 106 inside a high-ranking
    window.
  - In 5USJ chains B and F, the Tyr32 side chain is incomplete.
- **Cross-study differences:** inhibitor presence, mutation (WT vs G12D), crystal form, pH and resolution all change
  together. None of the differences between the studies is attributed to MRTX1133 or G12D.
- **Interpretation:** the recurring block lies in the Switch-I region, which is known to differ between nucleotide
  states. That is not independent validation. No ensemble, population, kinetic, functional-state or
  treatment-response claim is made.

## Global RMSD: verified value and unresolved published value

- **Measured:** 0.718 Å over the 168 mapped Cα of 9O0N vs 9O0O. MDAnalysis QCP and an independent numpy SVD/Kabsch
  calculation agree to within 2 × 10⁻⁹ Å (`KRAS-STRUCT-1_rmsd_verification.json`). The two-structure command
  reproduces the value exactly.
- **The published value:** the paper (Bonsor et al. 2026) reports 0.3 Å without an atom set, residue range, software
  or outlier rejection. Its Supplementary Information (37.2 MB) was not inspected. **The definition remains
  unresolved**, and no further effort will be spent reproducing an undocumented alignment procedure.
- **Exploratory only:** an unverified emulation of PyMOL-style iterative outlier rejection gave 0.27 Å on 141 of
  168 Cα. This shows only that outlier rejection *could* produce a value of that size. It is not evidence of what the
  authors did, and it is not used anywhere else.

## Carried forward

- **Reusable two-structure command:** `python -m topofold_inputs.structpair`. It is verified against both studies
  (`benchmarks/audit/structpair_demo/`); the guide is `docs/STRUCTPAIR_USER_GUIDE.md`.
- **Possible KRAS follow-up** (not started; only if a specific question needs it): a GDP/GMPPNP pair in one crystal
  form, to reduce the packing confound.
