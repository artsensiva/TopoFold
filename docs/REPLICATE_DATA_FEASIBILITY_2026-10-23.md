# Replicate-Data Feasibility (2026-10-23)

**Status:** feasibility only. Nothing was downloaded, scored or requested. C1 stays blocked and resumable; ABL stays
the preferred oncology-connected route if its data and permissions become available.

**Purpose:** find existing data with several simulation runs per condition, so that condition differences can be
compared with run-to-run variation (the limitation of DEV-ENS-1). Replicates strengthen the design; they do not by
themselves establish equilibration, biological relevance or clinical value.

## Candidates

| Dataset | Runs per condition; conditions | Construct / start / conditions | Length and format | Terms | Independent experimental endpoint | Use |
|---|---|---|---|---|---|---|
| **mdCATH, local non-cohort domains** (`1ha8A00`, six ineligible primaries); Mirarchi et al., *Sci. Data* 2024 | 5 replicas at each of 320, 348, 379, 413, 450 K | Same domain and force field in every run; **all replicas start from the same equilibrated configuration** (velocities differ) | 1 ns per frame, no timestamps; HDF5 (Å). `1ha8A00`: 51 residues, 430–450 frames per run at every temperature. `1s14B00`: 178 residues, 420–450 frames, backbone gap. The other five have some 120–210-frame runs | CC BY 4.0 (`CHECKPOINT_2026-10-13_PROVENANCE.md`); already local, pinned revision and SHA-256 recorded | None (temperature is a simulation condition, not a measured biological endpoint) | **Reproducibility test** |
| mdCATH, 30 cohort domains (local) | as above | as above | as above | CC BY 4.0 | None | **Excluded:** protected WP-D held-out set; not to be rescored |
| GPCRmd μOR apo vs BU72 (dynamics 151, 169) | 3 vs 3 | Same receptor (PDB 5C1M), both runs start from the agonist-bound active structure; CHARMM36m, ACEMD, 4 fs | 1.1 / 1.25 µs accumulated; 0.2 ns frames; XTC | Academic, non-commercial; commercial licence on request; large downloads by advance request | Pharmacology only; no residue-level measurement of the ensemble difference | Method calibration against PENSA; not first choice (terms, start-state confound) |
| ABL, Ayaz et al. 2023 (DESRES) | 3–10 runs per set (inventory inconsistent in the paper) | WT and resistance mutants from matched starts; a99SB-disp, 310 K | 10–20 µs runs; size unknown | Non-commercial, on request; request drafted, not sent | **HDX-MS (M472I vs WT), melting temperatures, kinetics** | **Biological validation** route, pending access |
| ABL MELD ensembles, Zenodo 10.5281/zenodo.21515187 (Brini lab) | One trajectory file per variant (WT, ligand, 4 single and 4 double mutants) | MELD: restraint-guided replica exchange, i.e. biased sampling; constructs and temperatures not stated | DCD, ≈ 0.3–1.5 GB per file | CC BY 4.0 | None stated; manuscript unpublished | **Not suitable:** no unbiased replicates; methods undocumented in the record |
| ATLAS (Vander Meersche et al., *NAR* 2024) | 3 replicates (project memo, not re-verified here) at one condition | Standardised CHARMM36m protocol | Per-protein downloads | CC BY-NC 4.0 (ATLAS site) | None for condition differences | Single condition only; not suitable for a condition comparison |

No openly licensed dataset was found that combines unbiased replicate runs per condition with an independent
residue-level experimental endpoint. The search was limited to existing project records, the ATLAS and GPCRmd
sites, and one Zenodo search.

## Recommendation

**Reproducibility test: mdCATH `1ha8A00`, 320 K vs 348 K, 5 replicas each.** Before any scoring, the run
specification would fix the frame range, burn-in rule (stated before seeing distances), bins and blocks.

- **Why:** already local, CC BY 4.0, never part of the WP-D held-out scoring, complete replicas at both temperatures,
  and small. It is exposed development data: WP-D schema check and held-out runner development validation
  (320/0, 320/1, 348/0, 348/1), EB-1 workload D2 and CMP-1 workload R3_D2 (corrected 2026-10-23 to include the WP-D
  validation and CMP-1). None of these compared per-window distributions between temperatures. `1s14B00` (backbone gap) is a possible
  second domain that also exercises break handling.
- **Question:** are between-temperature per-window distances larger than the replica-to-replica distances within a
  temperature, with the replica as the unit? Replicas share one starting configuration, so they measure
  reproducibility given that start, not independence of starting conditions.
- **Prerequisite (small engineering step):** a thin HDF5 → MDAnalysis in-memory Universe route using the existing
  `h5py` reader and the file's own topology, with units declared as Å and frame identity reported as not verifiable.
  No new file parser.
- **Not supported by it:** any biological or experimental claim; 320 K vs 348 K is a simulation-condition contrast.

**Biological validation stays separate:** ABL (DESRES trajectories with HDX-MS) remains the first choice if access
and terms are granted. No open alternative with an experimental endpoint was found in this search.
