# structpair feedback form (human usability trial)

**Use this form if you ran `topofold_inputs.structpair`.** It is the only form a structpair tester needs. (The general form
`docs/pilot/FEEDBACK_FORM.md` is for the ensemble tool, which this trial does not use; do not fill it in.)
This is a human usability trial of a research prototype, not a biological validation. Answer in your own words; short answers
are fine. Record only what happened, including failures and dead ends. Fill it in as you go; instructions:
`docs/pilot/STRUCTPAIR_HUMAN_TRIAL_PACKET.md`.

**About you (optional):** role / field; how often you compare protein structures computationally; tools you normally use.

## A. Trial record
a. **Environment:** OS, Python version, how the core was built, where MDAnalysis came from, the browser used for the viewer.
b. **Commit tested:** the commit you cloned.
c. **Commands:** paste the commands as run, including failed ones.
d. **Assistance:** none / some. If some: at which step, what you asked, and what help you were given (documentation pointer,
   explanation, workaround, code change). A documentation pointer counts as assistance.
e. **Failures and workarounds:** errors or dead ends, whether each was resolved, and any fix you devised yourself.
f. **Outcome per stage** (table below).

## B. Stages (fill in as you go; times can be approximate, say how you measured them)

Outcome classes: **unaided** (following the documentation) / **self-devised workaround** (your own fix, no help) /
**after assistance** / **not completed**.

| Stage | Outcome | Time | Notes (what happened, what you did) |
|---|---|---|---|
| Fresh checkout and environment, install (guide §2 steps 1–3) | | | |
| Installation check (guide §2 step 4) | | | |
| Own comparison: spec written and run (or fallback; say which) | | | |
| Manifest check (guide §3b) | | | |
| Viewer: opened, plots, table, exclusions, correspondence, provenance | | | |
| Interpretation written | | | |

## C. Questions

1. **Setup time:** from the start of the clone to the first successful structpair run. Anything you had to install or configure
   beyond the guide? Paste any error message and the command.
2. **Your comparison:** the question, the two structures (described, or named if you wish), the chains, the residue
   correspondence and the policies you chose, written **before** you looked at results. Own pair or fallback?
3. **Specification:** what was unclear when writing the spec (chains, correspondence, altloc and missing-residue policies)?
4. **3D view:** did it show the two Cα traces with the focus window highlighted? If not: what was shown (for example
   "WebGL is not supported"), which browser and operating system?
5. **Exclusions:** list each excluded window or group of windows and explain its cause in your own words. Which could you
   explain from the outputs alone? Was anything excluded that you expected to be included, or the reverse?
6. **Interpretation (before reading anything else in this repository about these structures):** what does the comparison show?
   What does it not show? What could explain the differences besides the factor you are interested in?
7. **Anything misleading:** any label, number or sentence that you read one way and later found meant something else.
8. **Viewer usability:** what was hard to find or read (for example truncated columns, axes, table selection)?
9. **Usefulness:** would this help in your own work, and for which question? Compared with what you do now, what is better,
   worse or missing? What single change would most increase its usefulness to you?

## D. Permission
10. May we quote your answers (anonymously / with attribution / not at all)?
11. May we record your name or affiliation in the repository's pilot log (yes / no)? Without a yes, the log describes you only
    as "outside tester 1" and stores no identifying details.
