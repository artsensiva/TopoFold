# TopoFold separate-agent computational-science usability review

Repository commit: `389f2d7d168cfa9b4b76af2366a4e2d1cce4289d`. Tool v0.2.0.

**Verdict: qualified operational pass.** Three CLI runs completed without code edits. Preparing input, running the workflow and reading the results succeeded. Identifying exact exclusion causes from the exclusions table alone was only partly successful.

## Scope

Installed-environment usability review on development data, not clean-install validation, independent human feedback or biological validation. The reviewer has prior familiarity from code reviews; this is not a blind novice trial. All created files are confined to this temporary review directory.

## Runs

| Run | Wall seconds | CPU seconds | Peak MiB | Exit |
|---|---:|---:|---:|---:|
| documented_example | 9.18 | 18.79 | 107.7 | 0 |
| reviewer_2yztA00 | 20.17 | 33.78 | 106.9 | 0 |
| gap_check | 17.07 | 21.47 | 110.2 | 0 |

## Scientific reading of the new input

Domain `2yztA00` was selected before scores were seen, from existing non-cohort development data. Five replicas at each temperature were compared over source frames 13-129 inclusive (117 frames each). Initial exclusion follows the existing 10% convention and establishes no equilibration. No independent experimental endpoint was used.

| Channel | Peak window | Pooled JS | Between-member median | Within 320 K | Within 348 K |
|---|---|---:|---:|---:|---:|
| kappa | 0:47-49 | 0.399 | 0.429 | 0.253 | 0.457 |
| tau | 0:56-59 | 0.438 | 0.705 | 0.559 | 0.865 |
| d_i_i2 | 0:47-49 | 0.410 | 0.442 | 0.273 | 0.455 |
| d_i_i3 | 0:46-49 | 0.421 | 0.480 | 0.293 | 0.509 |

The peak pooled differences coexist with substantial within-temperature variability. They identify locations for inspection, not established temperature effects. Kappa and its distance baseline highlight the same peak; that is not representation superiority. Torsion and the longer-distance baseline highlight different peaks, which also establishes no superiority.

## Checks

- All recorded output hashes for the three runs match. Independent entropy-form JS checks cover 254 new-domain rows; maximum difference from rounded CSV is 4.97e-07.
- C-alpha coordinate extraction, residue labels and frame indices match direct reads of the source HDF5 and embedded PDB. This is a same-decoder consistency check.
- Gap check: 10 channel-windows excluded. The exact cause is a numbering jump 2093 to 2117, found in the member JSON rather than the exclusions CSV.
- Repository preservation: 570 tracked files checked; changed files: []. Git status identical: True; HEAD unchanged: True.

## Practical obstacles

- Exact exclusion causes are not in exclusions.csv; intermediates are required. Deleting members/ removes that explanation.
- Environment setup is assumed: two existing interpreters and local MD data. A fresh-machine installation was not tested.
- Existing output directories can be overwritten; fresh output directories were used.
- The report displays a default global pbc value although per-member policies are the effective settings.
- The worked reading uses five replicas per group, whereas the handoff command uses three; identical numerical results should not be expected.

## Next step

Operationally usable on the documented installed environment; qualified pilot pass, with exclusion-cause criterion partly met. Next: inline exact exclusion reasons and first outside-user trial; no new inference campaign.

Machine-readable evidence: `REVIEW.json`. Exact specifications, commands, stdout/stderr, source provenance and per-run outputs are retained beside this report. No repository edits, builds, downloads, commits or pushes were performed. No fresh-environment install, independent human trial or biological validation was performed.
