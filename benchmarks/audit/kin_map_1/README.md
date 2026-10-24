# KIN-MAP-1: kinase residue-role mapping helper on the frozen KIN-ROB-1 panel (2026-10-24)

**Engineering evidence only.** No state labels, thresholds or classification; no biological or general-kinase claim.
KinCoRe labels were not used. **Not blind:** the helper was written after KIN-ROB-1 and KIN-ROB-1a were run, and its
unit tests use a fixture derived from the same panel. This is an acceptance and regression check against recorded
outputs.

## What was built

`python/topofold_inputs/kinase_roles.py` (`map_kinase_roles`, `observed_chain`, `map_kinase_roles_in_universe`):

- **Inputs, all explicit:**
  - a structure chain;
  - a stated reference (accession, isoform, sequence);
  - a numbering convention: reference index = label + a stated offset, or author numbering of the chain.
- **What it resolves:** four anchors (β3 Lys, αC Glu, HRD/YRD, DFG Asp), each motif-checked in the reference. The ten
  descriptor roles are derived from them.
- **How it maps:** structure residues are placed by sequence alignment, never by author numbers. Unobserved reference
  residues are cheap to skip only at a geometric chain break (C(i−1)–N(i) > 2.0 Å) or where bonding is unknown.
  Skipping between peptide-bonded residues is scored as a construct deletion.
- **What it records:** mutations, insertions, unaligned termini, unobserved ranges, construct deletions, chain breaks,
  insertion codes, modified-residue names (e.g. TPO, SEP), and the author-numbering offset. The offset is a diagnostic
  only.
- **What it refuses, with named reasons:**
  - motif or anchor-order contradictions;
  - motif windows with zero or several hits;
  - identity < 0.90 or < 30 aligned residues;
  - ambiguous placement, found by recomputing the alignment with each role's mapping forbidden;
  - role or motif residues mutated in the structure (unless allowed explicitly);
  - roles deleted in the construct;
  - missing or unaligned author residues.
- **Offsets are never applied.** When the motifs fail, offsets under which they would pass are reported and not used.

## Protocol and run

- **Freeze:** `run_spec.json` was frozen with the code, the fixture and the hashes of the KIN-ROB-1/1a artifacts
  (`results/KIN-MAP-1_freeze.json`) before the run.
- **Inputs:** the git-ignored KIN-ROB-1 downloads were read with their SHA-256 verified. Nothing was downloaded.

```bash
PYTHONPATH=python <python with MDAnalysis> benchmarks/audit/kin_map_1/run_kin_map.py run
python3 benchmarks/audit/kin_map_1/build_panel_fixture.py   # rebuilds the committed test fixture from the downloads
```

## Results (`results/KIN-MAP-1_results.json`, `results/KIN-MAP-1_summary.csv`)

**21 of 21 cases met their frozen expectation:** 19 accepted, 2 refused.

| Case | Chains | Outcome |
|---|---|---|
| Frozen KIN-ROB-1 spec, labelled as UniProt (offset 0) | 8 non-PKA chains | Accepted. The ten reference positions equal the spec, the author residues equal those KIN-ROB-1 used, and the descriptors equal the recorded values exactly |
| The same, PKA 1ATP E | 1 | **Refused.** All four motif mismatches are named (M72, N91, I-Y-R at 164, T-D-F at 184). The diagnostic reports that offset +1 would satisfy every motif; it was **not applied** |
| KIN-ROB-1a positions (UniProt, offset 0) | 1ATP E | Accepted; equal to KIN-ROB-1a |
| KIN-ROB-1 positions declared as literature numbering (offset +1 stated) | 1ATP E | Accepted; equal to KIN-ROB-1a |
| Author numbering of the residues KIN-ROB-1/1a used | 9 | Accepted; equal to the above |
| Author numbering of the residues the failed run used (71, 90, 163, 183) | 1ATP E | **Refused;** four motif mismatches named |

**Other checks:**
- **Independent parser:** `observed_chain` (MDAnalysis) equals the independent fixed-column parse for 9 of 9 chains.
- **Frozen code and artifacts:** the SHA-256 of `ca_adapter.py`, `atom_adapter.py` and `kinase_geometry.py` equal the
  KIN-ROB-1 freeze, and the 12 KIN-ROB-1/1a spec, runner, README and result files are unchanged.
- **Recorded, not refused:**
  - 1KV2 A: DFG-Gly and DFGxxX are unobserved at a chain break (ranges 115–122 and 170–184), as in KIN-ROB-1.
  - 1UWH A/B: one mutation, reference K522 observed as Ala.
  - 1ATP E: TPO and SEP are mapped to Thr and Ser; they are in the alignment because MDAnalysis' `protein` selection
    omits them.
  - Author numbering equals UniProt − 1 for 1ATP and 1UWH, and equals UniProt for the other 7 chains.
- **Resources:** 24.7 CPU-s; 0 bytes downloaded.

**Erratum found (not edited, because KIN-ROB-1 is frozen):** `../kin_rob_1/README.md` describes the 1UWH discrepancy
as "Ala at UniProt position 521". Its own results record it as PDB residue 521 (`C2_numbering_mapping.mismatches`),
which is UniProt 522 at offset +1. The KIN-ROB-1 results are correct; the README sentence is not.

## Tests

- **`tests/test_kinase_roles.py`** (no MDAnalysis; runs in every CI Python job):
  - SYNTHETIC cases: missing residues at a break; construct deletions; insertions; role, motif and non-role mutations;
    insertion codes; ambiguous and unique motif windows; an ambiguous alignment (`DFGG` with one Gly unobserved) that
    peptide-bond information resolves; an order contradiction; a wrong reference; numbering validation; author
    numbering.
  - The committed panel fixture (`tests/fixtures/kin_rob_1_panel_observed_chains.json`; public data, not SYNTHETIC)
    for the 9 chains, including the original PKA refusal and three correct PKA declarations.
- **`tests/test_kinase_roles_structure.py`** (MDAnalysis): generated SYNTHETIC PDB files with real chain-break
  geometry, insertion codes and altlocs; mapping → `extract_atoms` → `kinase_descriptors`; and, locally only, the
  independent-parser comparison on the panel files.

## Limits

- **Coverage:** 9 chains of 5 kinases, PDB format only, and canonical isoforms only.
- **Motifs:** the anchor motifs are strict (K, E, H/Y-R-D, D-F-G). A kinase with a variant motif (e.g. DLG) is refused,
  and the helper has no override for that.
- **Ambiguity:** checked for role and motif positions only, not for every residue.
- **Scoring:** the alignment scores are fixed engineering choices, not fitted.
