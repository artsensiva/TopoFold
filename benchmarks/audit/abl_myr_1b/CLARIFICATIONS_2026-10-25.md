# ABL-MYR-1b: reading guide for frozen wording (dated 2026-10-25)

The frozen records of ABL-MYR-1, 1a and 1b are unchanged. This guide states how some of their wording is to be read.

| Frozen wording (location) | Read as |
|---|---|
| a1 outcome `not_comparable`, K1 failed (a1 results and report, attempt records) | Production-start comparability was **unverified** because the reader did not expose the stored box. It does not demonstrate that the starts differ, and it does not demonstrate that they match. `not_comparable` is the frozen name for "criteria not met" |
| `potential_confound: []`, `potential_confounds: []` (generated reports and results JSON of the check) | No inspected item matched the predeclared confound rules. It means **no unexpected difference was detected among the inspected items**; it is not proof that no confound exists |
| `inherent_to_Myr_model` (classification lists) | Differences detected between conditions that follow from the presence of Myr (index-group names, the MYRP molecule, the water count) |
| K2 "thermostat_comparable", K3, K4 | Names of criteria with predeclared rules; passing means the rules were met for the inspected settings |
| `comparable_for_1b` (a2 results and report) | K1-K5 were met under the frozen rules and thresholds. It is not a statement that the production starts are identical (they are distinct, see the a2 report), nor that no confound exists |
| `PROTOCOL_1B.md`: "the thermostat targets, coupling settings, tc-group content and composition differ between conditions only in ways that follow from the presence of Myr" | for the inspected items, no unexpected difference was detected; the protocol's own list of what the evidence does not support applies in full |

The a2 production-start geometry (within- and between-condition Calpha RMSDs) is reported separately for each condition and set in the a2 report and in `PROTOCOL_1B.md`.
