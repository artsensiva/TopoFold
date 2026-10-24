# KIN-ROB-1: robustness of named-atom extraction on public kinase structures (2026-10-24)

**Engineering and representation evidence only.** No state labels, thresholds or classification; no biological,
general-kinase or performance claim. KinCoRe labels were not used.

- **Protocol:** `run_spec.json` was frozen with its code (`results/KIN-ROB-1_freeze.json`, 2026-10-24T15:39:29Z)
  **before any panel coordinate file was downloaded** (`panel_files_present_at_freeze: []`). The run downloads the panel
  and records SHA-256 values in `results/KIN-ROB-1_results.json`.
- **Panel:** 9 chains from 7 entries of 5 kinases.

  | Entries | Kinase |
  |---|---|
  | 1ATP E | PKA (mouse) |
  | 1HCK A; 1FIN A and C | CDK2 |
  | 1A9U A; 1KV2 A | p38α |
  | 1MQ4 A | Aurora A |
  | 1UWH A and B | BRAF |

- **Checks C1–C10** (spec): role identity in UniProt; numbering mapping and insertion codes; role-residue identity in the
  structure; required-atom coverage; missing atoms vs REMARK 465/470; altloc handling and residue consistency; an
  independent fixed-column parser; Cα bit-identity with `extract_ca`; determinism; bond-length plausibility.

```bash
PYTHONPATH=python <python with MDAnalysis> benchmarks/audit/kin_rob_1/run_kin_rob.py run   # downloads ~2.5 MB
```

## Results (`results/KIN-ROB-1_summary.csv`, `results/KIN-ROB-1_results.json`)

**Outcomes:** 7 pass, 1 pass with exclusions, 1 fail.

| Chain | Outcome | Cause |
|---|---|---|
| 1HCK A, 1FIN A, 1FIN C, 1A9U A, 1MQ4 A, 1UWH A, 1UWH B | pass | All checks pass. Coverage 20/20 role atoms. The independent parser agrees within 1.73 × 10⁻⁵ (Å or rad) |
| 1KV2 A | pass with exclusions | DFG-Gly 170 and DFGxxX 173 are unobserved (REMARK 465). The adapter reports them as `residue_not_found`, consistent with REMARK 465 (C5 pass). ψ(Phe) and the A-loop N–O distance are NaN with the named cause; the other 9 descriptors are computed |
| 1ATP E | **fail** | **Spec error caught by C1.** The PKA role positions were entered in conventional literature numbering (K72, E91, DFG 184–186), which omits the initiator Met, but labelled as UniProt P05132 positions. UniProt has K73, E92, D185-F186-G187 and Y165-R166-D167. Every role was therefore one residue early: "DFG-Phe" was Asp184, which has no CZ. C5 then failed as well, because REMARK 470 rightly lists no missing CZ |

**Other observations:**
- **Mapping:**
  - Numbering identity was 1.000 for 7 chains and 0.996 for 1UWH A and B.
  - The 1UWH discrepancy is one residue (Ala at UniProt position 521, Lys in UniProt).
  - PDB numbering equals UniProt numbering in every chain except 1ATP E and 1UWH A/B, where UniProt = PDB + 1. The
    mapping found and applied that offset.
- **Cα identity and determinism:** Cα bit-identity with `extract_ca` held under both altloc policies, and repeat
  extractions were bit-identical, in 9 of 9 chains.
- **Geometry:** no bond-length warnings.
- **Resources:** 51.6 CPU-s, 62.8 s wall, 2.51 MB downloaded.

## Amendment KIN-ROB-1a (post hoc, `../kin_rob_1a/`)

- **What changed:** after the failure was observed, the PKA positions were corrected to UniProt numbering (+1 for every
  role) in a separate spec, frozen at 2026-10-24T15:41:26Z.
- **What ran:** only 1ATP E, with a byte-identical copy of the runner (same SHA-256 `e191213d…`).
- **Result:** pass. C1 passes (YRD motif), coverage is 20/20, C5 passes, and the independent parser agrees within
  8.8 × 10⁻⁷.
- The KIN-ROB-1 result above is unchanged. The amendment shows the extraction works on PKA once the mapping is right; it
  does not convert the frozen failure into a pass.

## What this panel did not exercise

- **Alternate locations:** no role residue in the 9 chains has them, so real-structure altloc handling was not tested
  here. It is covered by the committed SYNTHETIC fixture (`tests/fixtures/SYNTHETIC_kinase_like_altloc.pdb`) and by
  9O0N in the KRAS tests.
- **Insertion codes:** none in any panel chain, so likewise covered only by the SYNTHETIC fixtures.
- **Formats:** only PDB-format files; mmCIF-only entries are untested.

## Lessons for the capability

- **Role positions** must be supplied in the numbering they are labelled with. The C1 identity check is necessary,
  because literature numbering can differ from UniProt by the initiator Met.
- **Missing residues and atoms** are reported consistently with the deposition records. A descriptor that depends on
  them becomes NaN with a named cause, never zero.
