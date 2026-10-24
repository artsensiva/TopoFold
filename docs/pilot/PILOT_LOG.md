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

## Entry 3 — INTERNAL SMOKE TEST of `structpair` by the implementer (2026-10-24)

**Not independent user validation and not a clean-machine install.** The implementer (a Claude agent) followed
`docs/STRUCTPAIR_USER_GUIDE.md` §2 in a fresh git clone of commit `3b2fad2`, writing to fresh output directories.
Record: `benchmarks/audit/structpair_demo/results/smoke_test_2026-10-24.json`.

- **Deviations from a clean install:**
  - The existing interpreters were used: the prebuilt core extension in `.venv`, and MDAnalysis 2.10.0 in the
    git-ignored audit venv.
  - Inputs were copied with `fetch_inputs.py --source`; the RCSB download path was not exercised.
  - No Rust or maturin build was performed.
- **Results:**
  - Inputs: 4 of 4 hashes ok.
  - Runs: 2 of 2 spec runs exited 0 (10.6 and 14.9 CPU-s).
  - Outputs: 6 of 6 compared output files (`windows.csv`, `residues.csv` and `exclusions.csv` for both runs) were
    byte-identical to the reference results.
- **Corrected before release:**
  - The guide's expected stdout now says 10 exclusion rows: 5 windows, listed once per structure.
  - The guide's per-run CPU estimate was corrected to the measured values.
- **Next:** an outside computational scientist (owner nominates; invitation updated to mention both tools). No feedback
  has been collected or invented.

## Entry 4 — AI TESTER 1: external AI technical trial of checkpoint 6ce7b49 (2026-10-23/24)

**An external AI technical trial, not human validation, a human usability assessment or biological validation.**
Human consent fields (feedback form items 15 and 16) do not apply. Archived evidence:
`benchmarks/audit/ai_tester_1_trial_2026-10-24/`. The tester's report is reproduced unchanged as `REPORT.md` there.

- **Tester:** AI tester 1.
  - Model family as reported: GPT-6. Exact model identifier: unknown. Application: Codex desktop, version unknown.
  - **AI helper disclosed:** a parallel AI subagent reviewed the documentation, prespecified the additional mapping
    (5US4 A vs 5USJ A, residues 55–72) before outputs were seen, and wrote the tester's numerical validator.
  - No human or project-developer assistance.
- **Tested commit:** `6ce7b49922311e831fc09e701fcb43ae9e314625`, branch `scientific-rebuild`, resolved before checkout.
- **Environment:**
  - One platform only: Linux x86-64, system Python 3.14.4 and Rust 1.98.1.
  - A fresh clone, a fresh venv, and Cargo/pip caches, build products, inputs and outputs, all inside the trial
    directory.
- **Permission event:** setup was first refused by the tester application's own approval review. It succeeded after
  the owner granted a trial-only exception. This was not a TopoFold failure.

| Area | Outcome category | Record |
|---|---|---|
| Installation | **Success using a tester-devised workaround** | The guide had no concrete venv or viewer `pip` commands and assumed a `.venv` path; the tester supplied them (`evidence/setup-notes.json`). Each install command exited 0. NumPy resolved to 2.5.4 versus the guide's 2.5.3 |
| Computation | **Success following the documentation** | Inputs fetched and hash-verified. The documented kras2 command exited 0 in 7.0 s. Outputs were byte-identical to the shipped reference, with 2,939 numeric comparisons at maximum difference 0. The tester's additional explicit mapping (55–72, `on_missing: exclude`) exited 0, with 56 cause rows explaining 26 excluded window/structure combinations |
| Viewer interaction | **Success using a tester-devised workaround**, with one TopoFold defect | Launched with tester-supplied paths and versions. Six metric axes and units, sorting, filtering, row selection, dual highlights, exclusion explanations and provenance were checked in a real browser. **Defect:** the table selection was out of sync with the focus after a filter change (below). Not tested: 3D rotation, hover, residue-name-mismatch and altloc cases, other platforms, the repository test suite |
| Interpretation | **Success following the documentation** | Exclusions explained correctly, including the 59→68 correspondence gap and the residues missing in 5US4 (60–64) and 5USJ (61–63). The report states that static coordinate differences establish no dynamics, populations, function, ligand causation or treatment relevance. This is the tester's written interpretation, not a human comprehension test |

**Not an overall pass.** Installation and viewer use needed tester workarounds, and one viewer defect was found.

**Time box:**
- The trial started at 23:49:31 UTC on 2026-10-23 and its services stopped at 07:33:42 UTC on 2026-10-24, 464.2
  elapsed minutes. A usage-limit interruption came in between, and the owner then requested resumption.
- Active time is unknown: computation took about the first 10 minutes and the resumed browser work about 10 minutes.
- **The continuous 60-minute cap was not met.** See `evidence/16-timing-and-interruption.json`.

**Tester-tooling mistakes, not TopoFold defects** (`evidence/18-tester-tooling-events.json`,
`evidence/validator-corrections.json`):
- Playwright locator mismatches;
- an HTTP 400 from the tester's own evidence form;
- two read-only inventory mistakes;
- four corrections to the tester's validator. The first attempt is preserved.

**Findings and actions** (fix commit `f96c3a2`, separate from the tested commit `6ce7b49`):
1. **Viewer selection out of sync — confirmed and fixed.**
   - **Reproduced before editing:** a headless-Chrome replay of the reported sequence on `6ce7b49` gave table
     position 33 against focus and details 58 after "rows: all" (`fix_verification/before_checkpoint_6ce7b49.json`).
     It also showed stale marks after re-sorting and after switching comparisons, and focus carried across
     comparisons.
   - **Cause:** the table selection was a row index that persisted when the displayed rows changed.
   - **Fix:** the table is keyed by what it displays (comparison, displayed positions, focus) and seeded with the
     focus row. Focus is remembered per comparison.
   - **After the fix:** 7 of 7 browser steps agree (`fix_verification/after_fix.json`), and 2 regression tests were
     added.
2. **Setup unclear — confirmed and fixed.** `docs/STRUCTPAIR_USER_GUIDE.md` §2 now gives one explicit
   isolated-install sequence based on the commands that succeeded in the trial. It assumes no `.venv` path and records
   versions as used, not as a compatibility range.
3. **Open, not acted on:**
   - `fetch_inputs.py` has no output-directory option (a setup note, not a blocker);
   - the declared Rust minimum is inaccurate (already KNOWN_ISSUES m12).

**Next:** a short confirmation run by the tester of the selection behaviour and the setup text, against the fix
commit.

## Entry 5 — AI TESTER 1: bounded confirmation of commit 049a47f (2026-10-24)

**A bounded external AI confirmation, not human validation, a human usability assessment or biological validation.**
- **Tester:** the same as entry 4 (model family as reported GPT-6; identifier unknown; Codex desktop, version unknown).
- **Commit tested:** `049a47f196e39ec46cc8e2c9d0bd63f03600df71` (it includes the fix `f96c3a2`).
- **Environment:** a fresh checkout and an isolated sibling environment on the same Linux machine.
- No parallel agent and no human or developer assistance.
- **Evidence:** `benchmarks/audit/ai_tester_1_trial_2026-10-24/confirmation/`. The tester's manifest verifies 47 of 47
  files, and the original `REPORT.md` is kept byte-for-byte.

| Item | Outcome |
|---|---|
| Setup instructions (guide §2) | **Sufficient in the tested Linux environment** (Python 3.14.4, Rust 1.98.1). No missing or extra package, version change or `--adapter-python` was needed. The commit checkout, trial-local paths and launch environment settings were applied and recorded (`evidence/COMMANDS.md`) |
| Computation | **Succeeded.** kras2 exited 0 in 7.0 s with the expected stdout; the added kras1 run exited 0 in 6.1 s. All 4 core outputs of each run are byte-identical to the shipped references |
| Persistent selection mismatch | **Not observed after the fix.** After selecting 58 and switching to "all", row 33 was unchecked |
| Settled selection states | **Synchronized at every requested step:** sidebar, details and checked row agree, or no row is marked when the focus is filtered out or off-screen. The tester scrolled the virtualised table to inspect rows |
| Transient states during comparison switches | **Qualified, unresolved observation, low severity.** Immediately after each comparison switch, accessibility snapshots briefly paired the new sidebar focus with the old details and checked row (84/58/58 entering kras1; 58/84/84 returning). Settled snapshots and screenshots were synchronized. The evidence cannot separate a visible transient from accessibility-cache lag or Streamlit's display of the previous state during a rerun |
| Continuous time cap | **Not met:** a usage interruption and a computer restart came between the position-sort check and the comparison switches. Elapsed time 3 h 57 min 38 s; active work about 14 minutes, exact value unknown |

- **Not concluded:** an unconditional "the table never mismatches" result.
- **Reopen the transient observation only if either:**
  - a mismatch remains after the Streamlit rerun has settled; or
  - a human visibly observes a persistent wrong selection.
- No application code was changed in response.
- **Tester-tooling events, not TopoFold defects** (`evidence/18-errors-and-limitations.json`):
  - one canvas accessibility click failed ("no layout object"); clicking by coordinates worked;
  - sandbox stream-fd warnings appeared, with exit 0;
  - the pre-interruption step records were reconstructed from tool outputs and screenshots.
- **Status:** the AI trial is closed, and viewer work stops here. An independent human trial remains a future
  product-usability milestone.
