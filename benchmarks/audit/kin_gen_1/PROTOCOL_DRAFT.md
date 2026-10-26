> **SUPERSEDED 2026-10-26 by `benchmarks/audit/kin_gen_2/PROTOCOL.md`** (reference-source re-check, one scoring procedure, development/evaluation separation by kinase). This draft was never frozen; its text below is preserved unchanged.

# KIN-GEN-1: kinase-descriptor generalization (protocol DRAFT, 2026-10-24)

**Status: DRAFT. Not frozen, not run.** No KinCoRe label file and no new structure has been downloaded.

**Blocked on:** the owner's decision about KinCoRe label terms (`benchmarks/audit/kincore_intake/intake_record.json`,
decision *conditional-go*). Items marked *[fix at freeze]* are chosen before any deciding calculation.

## Question

On a stratified set of public kinase crystal structures labelled by KinCoRe (Modi & Dunbrack definitions), which
representations separate the expert-defined conformational classes? The classes are spatial DFG (DFGin / DFGinter /
DFGout), C-helix in/out, and the extended activation loop where it is labelled. The representations are:
1. TopoFold's side-chain-aware descriptors (`kinase_geometry`: D1, D2, αC Cβ–Cβ, A-loop N–O, Phe χ1, φ/ψ);
2. Cα-only TopoFold κ/τ windows;
3. window-matched Cα distances;
4. local Cα RMSD;
5. Cα proxies of the established distances.

## Evidence classes (fixed before any data)

| Comparison | Class | Why |
|---|---|---|
| (1) vs KinCoRe labels | **Geometry implementation check** | KinCoRe labels are computed from the same published geometric definitions (D1/D2, dihedrals, Cβ distance). Agreement shows that TopoFold computes them as published. It is **not independent biological validation** and not evidence of performance |
| (2)–(5) vs KinCoRe labels | **Representation capability** (development evidence on public structures) | The labels are expert-defined and are not derived from TopoFold's Cα outputs (AGENTS rule 5). Every Cα representation receives the same chains, residues and windows (rule 6) |
| Anything | Never biological validation, general kinase performance, or a claim about function, druggability or clinical relevance | Static crystal structures; KinCoRe labels describe geometry |

## Sample (fixed at freeze from the label file alone, before computing any descriptor)

- **One chain per PDB entry**, chosen by a fixed rule *[fix at freeze, e.g. KinCoRe's representative chain, else the first
  chain]*.
- **Stratified by spatial label**, with at most *[fix: e.g. 2]* entries per kinase per class, so that a few heavily
  studied kinases cannot dominate.
- **Target size:** *[fix: e.g. ≤ 60 entries per class where available]*. The DFGinter class may be small; its count is
  reported, not padded.
- **Inclusion:** X-ray, resolution ≤ *[fix: 2.5 Å]*, all role atoms present (`on_missing="reject"` for the side-chain
  representation; the same chains for the others), no unresolved altloc (`blank_then_A`, recorded).
- **Exclusions:** the KIN-DESC-1 chains (1OPJ, 1OPK, 2F4J), reported separately if sampled.
- **Inputs:** structures from RCSB (CC0), with SHA-256 per file. The download budget is set at freeze *[e.g. ≤ 50 MB]*.

## Residue roles and mapping

- **Roles** (`atom_adapter.KINASE_DESCRIPTOR_ROLES`): β3-Lys, αC-Glu, αC-Glu+4, X-DFG and its predecessor, DFG-Asp,
  DFG-Phe, DFG-Gly, DFGxxX and X-HRD.
- **Source of positions:** per kinase, from KinCoRe's alignment or annotated positions if the permitted files provide
  them. Otherwise by a fixed motif rule written at freeze. Each assignment is confirmed by exact sequence matching to
  UniProt.
- **Mapping failures** are exclusions with a recorded cause, never guesses.

## Metrics

- **Per scalar descriptor and binary contrast** (DFGin vs DFGout; C-helix in vs out; A-loop extended vs not): AUROC, with
  the direction fixed in advance from the published definitions where they give one. Otherwise the absolute deviation
  from 0.5 is reported, and labelled as such.
- **Window representations:** only predeclared regions (DFG ± 3 residues, the A-loop from DFG to APE, the αC helix).
  Results are reported as the count of windows with AUROC ≥ *[fix: 0.8]* and the best window. Multiplicity is handled by
  the permutation control below, not by p-values.
- **Group-aware uncertainty:** leave-one-kinase-out ranges of each AUROC. No p-values or confidence intervals
  (decision D33).

## Controls (predeclared)

- **Label permutation within kinase:** labels are shuffled among each kinase's entries (*[fix: 200 permutations]*,
  fixed seed). The best-window AUROC of each representation under permutation sets the level that window selection
  alone can reach.
- **Mapping shift:** every role is shifted by +1 residue. Established-descriptor agreement should collapse, which tests
  that it depends on correct mapping.

## Outcomes (ordered; the first match applies)

1. **Inconclusive:** fewer than *[fix: 10]* entries in either class of a contrast, or the leave-one-kinase-out range
   includes 0.5.
2. **Separates:** AUROC ≥ *[fix: 0.8]*, above the 95th percentile of the permutation control, and the
   leave-one-kinase-out minimum ≥ *[fix: 0.7]*.
3. **Does not separate:** AUROC < *[fix: 0.6]* or not above the permutation control.
4. **Partial:** every remaining case.

The implementation check (1) has its own rule. It passes if TopoFold's D1/D2 and αC Cβ–Cβ reproduce KinCoRe's spatial and
C-helix labels with agreement ≥ *[fix: 95 %]* under the published thresholds or cluster boundaries. Disagreements are
listed by entry. It is never reported as validation.

## Stop conditions

- The owner has not approved the label terms.
- The label file fails SHA-256 recording or validation.
- The role mapping fails for more than *[fix: 20 %]* of the sampled entries.
- The download budget would be exceeded.
- Any change after data inspection that is not recorded as a dated amendment before scoring.
