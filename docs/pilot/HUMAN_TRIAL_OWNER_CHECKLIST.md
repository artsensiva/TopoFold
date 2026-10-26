# First independent human trial of structpair: owner checklist (not for the tester)

**A human usability trial, not a biological validation.** Nothing has been sent. The owner decides whom to invite and sends
the invitation personally; the project does not contact anyone.

## 1. Choose the tester
A computational scientist who has not worked on TopoFold, has not seen its audit records or AI-trial reports, and works with
protein structures. Record only what they permit (structpair form §D); otherwise the log calls them "outside tester 1".

## 2. What to send (exactly)
1. The invitation text in `docs/pilot/OUTSIDE_USER_INVITATION.md` (the part below the line), with `[name]`, `[Your name]`,
   `[repository link]` filled in and `[COMMIT]` replaced by the commit named in the handoff report of 2026-10-26.
2. A link to the repository at that commit (branch `scientific-rebuild`).
3. Three files, as links at that commit or as attachments: `docs/pilot/STRUCTPAIR_HUMAN_TRIAL_PACKET.md`,
   `docs/STRUCTPAIR_USER_GUIDE.md`, `docs/pilot/STRUCTPAIR_FEEDBACK_FORM.md`. (Do not send `docs/pilot/FEEDBACK_FORM.md`: it is the
   ensemble tool's form.)

Do **not** send: expected results, the AI-trial reports or their archive, audit READMEs, the Abl demonstration, or any hint
about which regions differ in the fallback example.

## 3. What to request back
- The completed `STRUCTPAIR_FEEDBACK_FORM.md`: A (trial record), B (per-stage outcome and time), C (setup time, comparison chosen
  before results, 3D view, exclusions, interpretation, anything misleading, viewer usability, usefulness) and D (permission).
- Optional: the spec, the stdout of each command, `manifest.json` of their run, screenshots.

## 4. During the trial
- Do not coach. If asked for help, answer only the question asked, and write down the time, the question and the answer;
  this is "assistance" in the record. A documentation pointer counts as assistance too.
- Do not fix code during the trial. If they hit a defect, note it; a fix comes after the trial, with its own commit.

## 5. After the trial
- Add an entry to `docs/pilot/PILOT_LOG.md` (next number) labelled HUMAN USABILITY TRIAL: tested commit, environment,
  per-stage outcome in four separate classes (unaided / self-devised workaround / after assistance / not completed), setup
  time, own pair or fallback, whether the 3D view worked, exclusions explained, the tester's interpretation as written, and
  every finding. Never roll partial completion into an overall pass. Archive the form (with permission) under
  `benchmarks/audit/` with SHA-256 hashes.
- Reproduce any reported defect on the tested commit before changing code.
