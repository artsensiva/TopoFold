# Pilot Log

## Entry 1 — INTERNAL WALKTHROUGH by the implementer (2026-10-23)

**This is not independent user validation.** The walkthrough was done by the same agent that wrote the tool and the
handoff, so it can find mechanical obstacles but cannot show that the documentation is clear to a new user. No user
feedback has been collected yet.

**Procedure:** followed `docs/PILOT_HANDOFF_ENSEMBLE_TOOL.md` without editing the implementation, in a separate
scratch directory.

| Step | Action | Result |
|---|---|---|
| 1 | Ran the handoff's example command (`docs/pilot/spec_example_1ha8A00.json`, 3 vs 3 replicas) | Completed: 18.3 CPU-s, 9.7 s wall, 110 MiB peak; 0 exclusions; report and tables written |
| 2 | Prepared a new spec for another local mdCATH development domain, `1s14B00` (320 K vs 348 K, replicas 0–1, frames 43–399) | Completed: 24.4 CPU-s, 20.4 s wall, 124 MiB; 10 exclusions |
| 3 | Identified excluded regions from `tables/exclusions.csv` | Region identified: windows spanning residues 2091–2119. The labels show a numbering jump from 2093 to 2117 in all four members |
| 4 | Explained limitations from `report.md` | The report states descriptive scope, pooling weights, the circular conditioning rule and the temporal diagnostics; replica variability is visible in the between/within columns |

**Success criteria** (handoff §7), internal reading:
- Criteria 1, 2 and 4 were met without editing the implementation.
- Criterion 3 was **partly** met. The region was identifiable, but the exclusion reason is generic ("missing
  residue, break or degenerate geometry"). The specific cause (`numbering_discontinuity` after output index 76) is
  only in `members/<run_id>.json` under `breaks`.

**Obstacles found:**
1. **Discovering mdCATH contents:** the replica keys and frame counts per temperature are not shown by the tool. The
   walkthrough had to open the HDF5 file with h5py to choose `stop` (the shortest replica had 420 frames). The
   handoff now includes a one-line inspection command.
2. **Generic exclusion reason:** the specific break type requires opening intermediate files. The handoff now says
   where. Showing it in `exclusions.csv` would be a small tool change; it is not made in this pilot.
3. **Two interpreters:** the core extension lives in `.venv`, while MDAnalysis/h5py live in a separate audit
   interpreter. A missing `--adapter-python` gives a clear message and exit code 2.
4. **Output location:** the handoff's example wrote `pilot_out/` inside the repository checkout. The walkthrough used
   a scratch directory; the handoff now recommends an outside path.
5. **Intermediate size:** `members/` held about 12 MB for 4 members × 357 frames (178 residues), larger than the
   final outputs (about 2.4 MB). It can be deleted after a run.
6. **Silent overwrite:** rerunning into an existing `--out` replaces outputs without a warning.
7. **"Analysis options" display:** `report.md` lists analysis options, including a default `pbc` value. The PBC
   policy actually used is per member, in `manifest.json`; this could be misread.

**Next pilot step:** give the handoff to at least one computational scientist who did not build the tool, and record
their run, obstacles and explanation of limitations in a new entry. Do not paraphrase or invent their feedback.

## Entry 2 — SEPARATE-AGENT REVIEW, internal pilot closed (2026-10-23)

**This is not independent human validation and not a fresh-install test.** A separate agent with prior familiarity
with the project followed the handoff in the existing environment. Evidence is archived in
`benchmarks/audit/pilot_review_2026-10-23/`.

- **Runs:** three CLI runs completed without implementation changes (`1ha8A00` example; reviewer-chosen `2yztA00`;
  `1s14B00` gap check).
- **Criteria:** preparing input, completing the run, and explaining limitations were met. Identifying exact exclusion
  causes from `exclusions.csv` was partly met (the cause was only in the member JSON).
- **Repository preservation:** verified (570 tracked files unchanged).
- **Follow-up done:** tool v0.3.0 writes exact causes per member into `tables/exclusions.csv` and `report.md`
  (adapter break records, residue status and feature reason codes, with locations and source frames).
- **Internal pilot status: closed.**
- **Next:** an outside computational scientist, using `docs/pilot/OUTSIDE_USER_INVITATION.md` and
  `docs/pilot/FEEDBACK_FORM.md`. No feedback has been collected or invented.
