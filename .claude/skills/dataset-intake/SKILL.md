---
name: dataset-intake
description: TopoFold dataset-intake checklist. Use before downloading, accepting or analysing any external dataset (MD trajectories, structures, HDX/NMR or other endpoints). Covers source, version, licence and permissions, checksums and raw-data inventory, large-download approval, independence and leakage, endpoint provenance, replicate and temporal structure, residue/chain/atom mapping, missing data, a machine-readable intake record and a go / conditional-go / no-go decision.
---

# TopoFold dataset intake

`AGENTS.md` and `CLAUDE.md` remain authoritative (especially rules 1, 2, 5, 6 and 9). This skill produces one intake
record per dataset before any analysis. It never scores data.

Usage: `/dataset-intake <dataset name, URL, DOI or local path>`. The record is written from
`intake_record.template.json` (this directory) to the location the task names, by default
`benchmarks/audit/<dataset_slug>_intake/intake_record.json`. Unknown values stay `null` with a note; never infer a value
from figures or guess it.

## 1. Source, version, licence and permitted use

- Primary source (repository, accession, DOI, paper) and the exact version, release date or commit.
- Licence text as stated by the source (quote it; record the URL and access date). Distinguish the data licence from the
  paper licence.
- Permissions needed for this project, each as granted / refused / conditional / unknown with evidence:
  internal evaluation, publication of derived results, redistribution of derived files, later commercial use.
- Non-commercial or "on request" terms are recorded as conditions, not assumed.
- Verify every identifier (DOI, PDB, PRIDE, UniProt) against the source (AGENTS rule 9).

## 2. Size and download approval

- Record total size and per-file sizes from listings or headers before downloading.
- Downloads above the task's stated budget (default: ask above 100 MB, or anything the task forbids) need explicit owner
  approval first. Metadata and listings are not trajectory downloads; state which was fetched.
- Downloaded files go into their own new directory, are treated as untrusted, and are read with `python -I` scripts kept
  elsewhere. Large raw data stay git-ignored.

## 3. Checksums and immutable raw inventory

- SHA-256 for every file at receipt, stored in the intake record; source checksums compared when published.
- Raw files are never edited; derived files go to a separate directory with their own hashes.
- A missing real input raises an error; never synthesise a substitute (AGENTS rule 2).

## 4. Independence, prior exposure and leakage

- Has TopoFold (or its developers) already seen or tuned on this dataset? Name prior use (e.g. development runs).
- Overlap with other datasets used for development or evaluation (same proteins, same trajectories, same structures).
- Are labels, states or endpoints derived from TopoFold output or from methods sharing its inputs? If so, not independent.
- Would any compared method receive extra information (target region, labels)? Record equal-information constraints.

## 5. Endpoint provenance (biological or expert-defined)

- What is measured, by whom, with which method; processed values available or only raw/figures.
- Experimental conditions vs simulation conditions (construct, ligand state, mutation, temperature, pH). Record every
  mismatch explicitly; a mismatch can make the comparison cross-state.
- Replicates and uncertainty supplied by the endpoint providers; significance rules used by them.
- Whether the endpoint is conformational (HDX, NMR, crystallographic state) or not (activity, potency).

## 6. Replicate structure and temporal dependence

- Number of runs per condition, lengths, frame interval, saved atoms, equilibration notes.
- Independence of runs: independent starting structures vs same start with new velocities (conditionally independent)
  vs segments of one trajectory. Record which.
- Expected autocorrelation and whether effective sample size can be estimated; incomplete or truncated runs.

## 7. Residue, chain and atom mapping

- Reference numbering (e.g. UniProt canonical isoform) and the offset or alignment for each source, derived by exact
  sequence matching, not by assumption. List every mismatch (mutations, tags, engineered residues).
- Chain selection rule fixed before looking at results; alternate locations and insertion codes recorded.
- Atoms available (Cα only, backbone, heavy atoms, hydrogens) vs atoms the planned descriptors need.

## 8. Missing data and exclusions

- Unobserved residues, missing atoms, gaps, failed or truncated runs, peptides without replicates.
- Exclusion rules written before analysis, each with its cause; counts recorded.

## 9. Decision

| Decision | When |
|---|---|
| **go** | licence and permissions cover the intended use; checksums recorded; mapping verified; endpoint independent; replicate structure adequate for the planned question |
| **conditional-go** | usable once named conditions are met (permission answer, approval of a large download, a reader, an owner decision on a confound); list each condition and who must act |
| **no-go** | licence or permissions prohibit the use, endpoint not independent, mapping impossible, or data unavailable |

Report the decision with its reasons, the open conditions, and what was and was not downloaded. State that nothing was
sent externally unless the owner sent it.

## Output of this skill

1. The completed intake record (JSON, from the template), with `null` plus a note for anything unknown.
2. A short summary in the required TopoFold report format (see the `scientific-audit` skill, section 10).
