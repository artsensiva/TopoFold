# structpair human usability trial: tester packet

**This is a human usability trial of a research prototype, not a biological validation.** It asks whether a computational
scientist who has not worked on TopoFold can install `topofold_inputs.structpair`, run a comparison that matters to them, read
the result in the viewer and explain what it means. Any outcome is useful, including "I could not get it to work" or "this
would not help me".

You need: a machine with Python 3 (the package declares 3.9 or later; only Linux x86-64 with Python 3.14 has been
tried so far, so any other platform is itself useful information), git, a Rust toolchain (`rustup`), a desktop web browser, and
about 1–2 hours. Everything runs locally. Nothing is uploaded, and nothing is sent to the project unless you choose to send it.

## What to do

1. **Start the clock and keep notes as you go** in your copy of `docs/pilot/FEEDBACK_FORM.md` (§0 and the structpair section
   at its end). Paste commands and error messages as they happened, including failed ones.
2. **Fresh checkout and environment.** Clone the repository at the commit named in the invitation into a new directory and
   create a new Python environment for it. Do not reuse an existing TopoFold checkout or environment.
3. **Install and check the installation** by following `docs/STRUCTPAIR_USER_GUIDE.md` §2, steps 1–4. That step runs a
   documented example whose expected terminal output is given in the guide; it is an installation check only.
4. **Run one comparison of your own.** Choose two experimental or model structures that matter to your work, write a spec
   (guide §3; several chains per entry: §3a) and run the command once. Decide your question, the chains and the residue
   correspondence **before** you look at any result, and write them in the form.
   - **If you have no suitable pair**, use the documented fallback instead: the second example of the guide's input set,
     ```bash
     PYTHONPATH=python python -m topofold_inputs.structpair \
         benchmarks/audit/structpair_demo/specs/kras1_9O0N_vs_9O0O.json --out ~/structpair_out/kras1
     ```
     (its inputs were downloaded by `fetch_inputs.py` in step 3). Record in the form that you used the fallback.
5. **Check the result against its manifest** (guide §3b) and **open it in the viewer** (guide §2 step 5 and §6). Look at the
   plots, the window table, the exclusions, the residue correspondence and the provenance. Note whether the 3D view showed the
   two Cα traces with the focus window highlighted; if it stayed empty, note the message and your browser.
6. **Write your interpretation before reading anything else about the structures in this repository**: what the comparison
   shows, which exclusions there are and why, and what it does not show. Use only the outputs, the guide and your own
   knowledge.
7. **Finish the form** (sections 0–6 and the structpair section) and stop the clock.

## Rules for a fair trial

- Work from the documentation. If you get stuck, try to solve it yourself first and record what you did (a **self-devised
  workaround**). If you then ask the project for help, record exactly what you asked and what you were told (**assistance**).
  These are recorded separately and both are useful results.
- There is no expected answer and no target finding. Please do not look for one in the repository's audit folders before you
  have written your interpretation.
- Do not send structures or results you are not allowed to share. You may describe your inputs without naming them.

## What to send back

The completed feedback form. Optionally: your spec, the stdout of each command, and `manifest.json` of your run (it contains
file paths and hashes, not coordinates). Screenshots of anything confusing are welcome.
