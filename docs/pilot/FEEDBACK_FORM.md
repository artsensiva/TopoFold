# Outside-user feedback form: ensemble tool (`topofold_inputs.ensemble`)

> **structpair testers: do not use this form.** The structpair human usability trial has its own, shorter form:
> `docs/pilot/STRUCTPAIR_FEEDBACK_FORM.md`. Sections 2–4 below ask about the ensemble tool (groups of trajectories, pooled
> JS distances, replicas) and do not apply to structpair.

**Tool:** `topofold_inputs.ensemble` (groups of trajectories, `docs/PILOT_HANDOFF_ENSEMBLE_TOOL.md`).

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
f. **Outcome** (tick one):
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

## 2. Preparing input (ensemble tool only)
4. Did you run the documented example? Did you write a spec for your own data? Which input type (MD files or mdCATH)?
5. What was unclear about the specification (groups, members, frame slices, PBC)?

## 3. Exclusions (ensemble tool only)
6. Were any regions excluded? From `tables/exclusions.csv` / `report.md`, could you tell which residues, which
   members and why?
7. Was anything excluded that you expected to be included, or the reverse?

## 4. Reading the report (ensemble tool only)
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
