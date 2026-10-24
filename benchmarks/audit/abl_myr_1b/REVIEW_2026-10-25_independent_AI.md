# Independent AI scientific review of the frozen ABL-MYR-1b protocol (recorded 2026-10-25)

**Attribution and status.** This is an **independent AI scientific review**, not human peer review. The text below was supplied
to the repository on 2026-10-25 by the repository owner as the review's content; the reviewing model is not identified in the
supplied text. It reviewed `PROTOCOL_1B.md` against `benchmarks/audit/abl_myr_1/run_spec.json` and the frozen runner. It is
recorded unedited; only the parts that are review findings are included (the supplied message also contained a note on the
GitHub landing page and an instruction for the implementation step, which are not part of the review). Nothing in the frozen
protocol, the earlier freezes, the results or the owner's files was changed by recording it.

## Findings (verbatim)

The a2 check supports comparing these simulated Myr and noMyr ensembles descriptively. It does not establish a causal Myr effect
or biological validation. The source paper supports the chosen 1,000 ns exclusion and separate thermostat runs.

Answers to the protocol's six review questions:

1. **Comparability:** K1–K5 are sufficient to proceed with the restricted descriptive question, provided every remaining gate
   passes. They do not rule out unexamined confounds.
2. **Thermostats:** The sets provide a useful consistency check across simulation settings. They share starting structures, so
   agreement is not independent replication.
3. **Starting structures:** The small start-contrast value is worth reporting, but it uses Cα RMSD; the main window statistic
   uses histogram distances. The numbers cannot be compared as a calibrated control for the main result.
4. **Decision rule:** E > 0 in leave-one-run-out checks and the three-residue cutoff are reasonable descriptive filters. They
   provide no error-rate control, significance, or evidence of technological superiority.
5. **Time and features:** Excluding the first 1,000 ns follows the paper. Backbone frames every 2 ns can support the planned
   geometry description, but cannot establish kinetics; temporal dependence must remain visible in the report.
6. **αI helix:** Its absence would not prevent a survey of the kinase-domain windows that remain. It would prevent a direct claim
   about the αI helix or the Myr-binding region.

**Execution ambiguity to settle before scoring:** the original runner can continue with one thermostat set if the other falls
below eight usable runs, then label the outcome inconclusive. The 1b protocol also says a failed gate stops the run. Since the
question requires agreement across both sets, the review recommends a dated amendment saying that loss of either set stops the
bilateral comparison. Keep both frozen specifications unchanged.

## How the repository responds

- The G5 recommendation, a clarification that leave-one-run-out checks cover every retained run, the wording corrections for
  answers 2–5 and the αI rule are made in the dated, separately frozen amendment `protocol_amendment_a1.json`.
- `protocol_1b.json`, `PROTOCOL_1B.md`, the original specification and every earlier record stay unchanged.
