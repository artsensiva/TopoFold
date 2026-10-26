# Invitation for the first independent human trial of structpair (ready to send; not sent)

The project owner sends this text with the three attachments or links listed in
`docs/pilot/HUMAN_TRIAL_OWNER_CHECKLIST.md`. Replace the bracketed fields; nothing else needs editing. The earlier draft for the
ensemble tool is in the git history of this file (before 2026-10-26).

---

**Subject:** Would you try a small structure-comparison tool and tell me where it is confusing? (human usability trial, 1–2 hours)

Dear [name],

I am developing TopoFold, a research prototype. One of its commands, `structpair`, compares two protein structures residue by
residue on a correspondence you state: for every 3- and 4-residue window it reports local Cα backbone geometry, distance
baselines and local RMSD, lists every excluded window with its cause, and records hashes of inputs and outputs. A local,
read-only viewer shows the result. It is descriptive only: no statistics, no p-values and no claims about function.

Before going further I would value an outside computational scientist's view. **This is a usability trial of the software,
not a biological validation**, and there is no right answer to find.

**What I would ask (about 1–2 hours, on your own machine):**
- install it in a fresh checkout and environment, following the guide;
- run one comparison of two structures relevant to your own work (a documented fallback example is provided if you have none);
- look at the result in the viewer and write, in your own words, what you think it means and what it does not show;
- fill in the feedback form, noting where you got stuck, anything you had to work around, and any help you asked for.

Please work from the documentation alone at first. If you do ask me for help, that is fine; the form simply records it.

**Your data stay with you.** Nothing is uploaded. Send back only the feedback form, plus anything else you choose to share.
I will not quote you or name you anywhere without your permission (form §6).

**Materials:** the repository at commit `[COMMIT]` ([repository link]), the tester packet
(`docs/pilot/STRUCTPAIR_HUMAN_TRIAL_PACKET.md`), the guide (`docs/STRUCTPAIR_USER_GUIDE.md`) and the feedback form
(`docs/pilot/FEEDBACK_FORM.md`).

Any finding, including "this was not useful to me", helps. Thank you for considering it.

Best regards,
[Your name]
