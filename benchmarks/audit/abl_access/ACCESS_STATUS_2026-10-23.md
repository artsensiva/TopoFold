# ABL Biological Route: Access Status (2026-10-23)

ABL (Ayaz et al., *Nat Commun* 14:1885, 2023) remains the preferred biological-validation route. This records what was
obtained, what is blocked, and by whom it can be unblocked. No trajectories were requested or downloaded, and no
scoring was done.

| Item | Status | Evidence |
|---|---|---|
| Trajectory run inventory (apo WT / M472I, matched starts) | **Not obtained.** Only D. E. Shaw Research can provide it, on request; the repository does not record that the request was sent | Draft: `docs/ABL_DESRES_DATA_REQUEST_DRAFT.md` (owner to review and send) |
| Trajectory permissions | **Not obtained.** The paper states "non-commercial use … through contacting trajectories@deshawresearch.com" | Same request covers internal evaluation, publication, redistribution and later commercial use |
| HDX-MS raw data | **Public:** PRIDE PXD034008 (CC0); 43 raw Waters archives (WT, A344V, M472I; 10 s–4 h; duplicates), roughly 0.2–0.5 GB each. **Not downloaded.** The authors' processing pipeline (PLGS/DynamX) is not available locally; other processing routes were not evaluated | PRIDE project and file listing (API, read 2026-10-23) |
| HDX-MS processed peptide map | **Obtained** (CC0, 3,386 bytes): 80 peptides, construct residues 10–284, covering 252 of 275 positions. Columns: sequence, start, end, MH⁺, max uptake, retention time. **It contains no deuterium-uptake values.** The numbering is construct-relative, so mapping to c-Abl residues is needed | `data_pride_PXD034008/REI_ABL_KB_dynamX_export_18May22.csv` (git-ignored), SHA-256 `2fd6b624…` |
| Peptide-level deuterium uptake per state and time | **Not found in the inspected resources:** the PRIDE file listing (one processed file, the peptide map) and the paper's Source Data (MOESM6, 4.49 MB, downloaded; sheets for kinetics, FoldX, melting temperatures and IC50s, no HDX sheet). Processed values may exist elsewhere (for example with the authors) | `data_sourcedata_ayaz2023/41467_2023_36956_MOESM6_ESM.xlsx` (git-ignored), SHA-256 `1b32799c…` |

## To unblock (owner actions)

*Update 2026-10-24: the owner has sent both requests. Replies are pending, and whether apo M472I trajectories exist is
unknown until D. E. Shaw Research replies.*

1. Send the D. E. Shaw Research request (run inventory and permissions).
2. Ask the HDX authors (Engen laboratory, co-authors of the paper) for processed peptide-level uptake for apo WT and
   M472I with replicate values, or confirm a processing route for the PRIDE raw files.
   Drafts for both requests, not sent: `docs/ABL_DESRES_DATA_REQUEST_DRAFT.md` and `docs/ABL_HDX_DATA_REQUEST_DRAFT.md`.

## Alternative routes (kept open)

- A bounded search for openly licensed replicated MD with an independent experimental endpoint found no suitable
  data set. ATLAS and GPCRmd remain non-commercial options without such endpoints.
- Other routes stay resumable: C1 (Rascore terms), GPCRmd, and further mdCATH non-cohort domains for method work.
