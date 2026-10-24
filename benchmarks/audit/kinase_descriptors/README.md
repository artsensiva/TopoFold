# KIN-DESC-1: kinase-state descriptor feasibility on four public Abl chains (WP2, 2026-10-24)

**Representation capability only.** These are public crystal structures (static). There is no biological, clinical or
novelty claim. Protocol: `run_spec.json`, frozen with the code in `results/KIN-DESC-1_freeze.json` before any deciding
value was computed.

```bash
.venv/bin/python benchmarks/audit/kinase_descriptors/run_kinase.py run   # needs data_pdb/ (git-ignored; SHA-256 in run_spec.json)
```

**Audit and provenance addendum (2026-10-24):** `AUDIT_2026-10-24.md` and
`results/KIN-DESC-1_provenance_addendum_2026-10-24.json`. They record output hashes and software versions; no frozen
file was changed.

## Design

- **Chains:**
  - 1OPJ A and B: DFG-out, A-loop closed; imatinib bound.
  - 1OPK A: DFG-out, A-loop open; PD166326 ("P16") in the entry; D363N.
  - 2F4J A: DFG-in, A-loop open; VX-680; H396P.
  - State labels are from the text of Ayaz et al. 2023, which is independent of TopoFold. The paper calls 1OPK "apo",
    but its HET records list P16; this is recorded, not resolved.
- **Numbering:** c-Abl isoform 1a. 1OPJ and 1OPK use isoform-1b numbering (offset −19). Offsets were found by exact
  sequence matching to UniProt P00519, with 281–446 identities per chain and the mismatches listed in the spec.
- **Representations:**
  - **Established** (Modi & Dunbrack, *PNAS* 116:6818, 2019, verified from PMC6452665):
    - D1 = Cα(E290)–Cζ(F382) and D2 = Cα(K271)–Cζ(F382);
    - αC Cβ(E286)–Cβ(K271), with ≤ 10 Å meaning αC-in;
    - the extended-A-loop backbone N(386)–O(360) distance;
    - Phe382 χ1;
    - φ/ψ of residues 380–382.
  - **Cα-only:**
    - TopoFold κ (3-residue windows) and τ (4-residue windows);
    - the window-matched Cα distances;
    - pairwise local Cα RMSD (SVD/Kabsch);
    - Cα proxies of the established distances.
- **Rule:** a representation separates a contrast if the smallest between-group difference exceeds the largest
  within-group difference and is at least the floor (0.5 Å or 0.2 rad, fixed in advance).
  - Window representations are evaluated only in predeclared regions: DFG 378–386, A-loop 381–409, αC 279–293.
  - Every Cα representation uses identical chains, residues and windows.
- **Reference check (rule 8):** the dihedral routine gives the analytic IUPAC values for the cis, trans and −90°
  reference cases. Over 2,000 SYNTHETIC quadruples it agrees with an independent formula to 3.6 × 10⁻¹⁵ rad, once that
  formula's opposite sign convention is accounted for. This was done before the freeze.

## Results (from `results/KIN-DESC-1_results.json` and `results/KIN-DESC-1_summary.csv`)

**DFG (in: 2F4J; out: 1OPJ A/B, 1OPK)**
- **Established:**
  - D1 separates: 6.83 Å vs 12.85–13.91 Å.
  - D2 separates: 13.63 Å vs 10.83–11.52 Å.
  - φ381 and φ382 separate.
  - Phe382 χ1 does not: 1OPK's χ1 (+53°) differs from all other chains (about −50°).
- **Cα proxies:**
  - proxy_D2 separates (12.31 Å vs 9.41–10.06 Å).
  - proxy_D1 does not (7.32 Å vs 7.99–12.79 Å).
- **Windows** (separating / tested):

  | Representation | DFG windows separating |
  |---|---|
  | τ | 2/6 (hinges 380 and 383) |
  | Local RMSD, 4 residues | 2/6 |
  | d(i−1, i+2) | 1/6 |
  | κ | 0/7 |
  | d(i−1, i+1) | 0/7 |
  | Local RMSD, 3 residues | 0/7 |

**A-loop (closed: 1OPJ A/B; open: 1OPK, 2F4J)**
- **Established:** N(386)–O(360) separates: 16.1–16.2 Å vs 2.87–5.56 Å. Of the open chains, 2F4J is within the 3.5 Å
  H-bond reference and 1OPK is not.
- **Cα proxies:** proxy_Aloop (Cα386–Cα360) separates (19.25–19.39 vs 5.27–6.49 Å), and so does proxy_D1.
- **Windows:**

  | Representation | A-loop windows separating |
  |---|---|
  | τ | 13/26 |
  | Local RMSD, 4 residues | 12/26 |
  | d(i−1, i+2) | 10/26 |
  | κ | 6/27 |
  | d(i−1, i+1) | 4/27 |
  | Local RMSD, 3 residues | 1/27 |

**αC: not assessable.** All four chains are αC-in by the Cβ–Cβ criterion (7.96–8.26 Å), so there is no contrast.

**Resources:** 1.1 CPU-s. The three PDB files total 1.06 MB.

## Interpretation (within the frozen limits)

- **Consistency with the labels:** the established side-chain-aware descriptors (D1 and D2 for DFG, the backbone H-bond
  distance for the A-loop) reproduce the paper's labels.
- **What Cα geometry captures on these chains:**
  - **A-loop open/closed** (a large backbone rearrangement) is separated by many Cα windows and by the Cα proxy.
  - **DFG in/out** is separated only at a few four-point windows (τ, 4-residue RMSD, d(i−1, i+2)) and by one Cα proxy.
  - Three-point windows (κ, d(i−1, i+1), 3-residue RMSD) do not separate DFG at all.
- **The Cα proxies don't stand in for the side-chain descriptors:** proxy_D1 separates the A-loop contrast rather than
  DFG, because the Phe ring position (Cζ) is not represented by Cα. This is consistent with the stated physics: χ1 and
  ring placement are invisible to Cα.
- **Practical implication for the ABL track, a hypothesis to test there:**
  - DFG and αC states need side-chain atoms (Cζ of F382, Cβ of K271 and E286). An evaluation that relies on them needs
    capability "B" (backbone/side-chain extraction).
  - Cα κ/τ alone is likely enough to see A-loop rearrangements.
- **Limits:**
  - 4 chains in groups of 1–3; two of them are crystal copies of one entry.
  - Ligands, mutations (D363N, H396P), constructs and crystal forms differ between groups.
  - Complete separation here is not a classification accuracy, and nothing generalizes beyond these chains.
