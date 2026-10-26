# Outside-user feedback form (descriptive comparison tools)

**Which tool did you run?** `topofold_inputs.ensemble` (groups of trajectories, `docs/PILOT_HANDOFF_ENSEMBLE_TOOL.md`)
and/or `topofold_inputs.structpair` (two experimental structures, `docs/STRUCTPAIR_USER_GUIDE.md`). For structpair, read
"spec" below as the two-structure spec, and "members" as structures A and B.

Please answer in your own words; short answers are fine. Record only what happened.

**About you (optional):** role / field; how often you analyse MD ensembles; tools you normally use for this.

## 0. Trial record (filled in with the facilitator; kept private unless §6 allows otherwise)
a. **Environment:** OS, Python version, how the core was built (or not), where MDAnalysis came from, and whether a
   viewer was used.
b. **Task:** the comparison you attempted and why it is relevant to your work (inputs can be described without naming
   them).
c. **Commands:** paste the commands as run, including failed ones.
d. **Assistance:** none / some. If some: at which step, what was asked, and what help was given (documentation
   pointer, explanation, workaround, code change).
e. **Failures:** errors or dead ends, and whether each one was resolved.
f. **Outcome** (tick one; for structpair, tick one per stage in §7):
   - completed **unaided**, following the documentation;
   - completed with a **self-devised workaround** (your own fix, no help; describe it in e);
   - completed **after assistance** (list the steps in d);
   - not completed.
g. **Understanding:** could you explain each exclusion from the outputs alone? Could you state what the result does
   **not** show (no populations, kinetics, functional state, or causal attribution to a ligand or mutation)?

## 1. Installation and environment
1. Did you use the provided environments or install your own? What did you have to install or configure?
2. Time from start to the first successful run (approximate).
3. Any errors? Please paste the message and the command.

## 2. Preparing input
4. Did you run the documented example? Did you write a spec for your own data? Which input type (MD files or mdCATH)?
5. What was unclear about the specification (groups, members, frame slices, PBC)?

## 3. Exclusions
6. Were any regions excluded? From `tables/exclusions.csv` / `report.md`, could you tell which residues, which
   members and why?
7. Was anything excluded that you expected to be included, or the reverse?

## 4. Reading the report
8. In your own words, what does the "largest pooled JS distances" table tell you, and what does it not tell you?
9. Could you tell whether a pooled difference was larger than the variation between replicas of the same group?
10. Were the temporal diagnostics understandable? Did they change how you read the differences?
11. Did you find any statement in the outputs misleading or overstated?

## 5. Practical usefulness
12. Would this help in your own work? For which question?
13. Compared with what you do now (e.g. RMSF, per-residue RMSD, PCA, PENSA or other tools), what is better, worse or
    missing?
14. What single change would most increase its usefulness to you?

## 6. Permission
15. May we quote your answers (anonymously / with attribution / not at all)?
16. May we record your name or affiliation in the repository's pilot log (yes / no)? Without a yes, the log
    describes you only as "outside tester 1" and stores no identifying details.

## 7. structpair human usability trial (`docs/pilot/STRUCTPAIR_HUMAN_TRIAL_PACKET.md`)

Fill in as you go. Times can be approximate; say how you measured them.

| Stage | Outcome: unaided / self-devised workaround / after assistance / not completed | Time | Notes (what happened, what you did) |
|---|---|---|---|
| Fresh checkout and environment, install (guide §2 steps 1–3) | | | |
| Installation check (guide §2 step 4) | | | |
| Own comparison: spec written and run (or fallback; say which) | | | |
| Manifest check (guide §3b) | | | |
| Viewer: opened, plots, table, exclusions, correspondence, provenance | | | |
| Interpretation written | | | |

17. **Setup time:** from the start of the clone to the first successful structpair run.
18. **Your comparison:** the question, the two structures (described, or named if you wish), the chains, the residue
    correspondence and the policies you chose, written **before** you looked at results. Own pair or fallback?
19. **3D view:** did it show the two Cα traces with the focus window highlighted? If not: what was shown (for example
    "WebGL is not supported"), which browser and operating system?
20. **Exclusions:** list each excluded window or group of windows and explain its cause in your own words. Which ones could you
    explain from the outputs alone?
21. **Interpretation (before reading anything else in this repository about these structures):** what does the comparison
    show? What does it not show? What could explain the differences besides the factor you are interested in?
22. **Anything misleading:** any label, number or sentence that you read one way and later found meant something else.
23. **Viewer usability:** what was hard to find or read (for example truncated columns, axes, table selection)?

