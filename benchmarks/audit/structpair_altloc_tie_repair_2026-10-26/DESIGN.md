# structpair altloc-sensitivity tie handling (KNOWN_ISSUES m13): design, written before the implementation (2026-10-26)

**Engineering repair with regression tests. Not a new evaluation.** KIN-PAIR-1 keeps its frozen-rule outcome (IMPLEMENTATION FAILURE) and all its files;
its six pairs are replayed here only as an engineering regression, and the replay is stored in this directory.

## The defect
`altloc_sensitivity: highest_occupancy` is an optional second coordinate-selection policy. When a mapped Cα has alternate locations tied on
occupancy, the adapter raised and structpair aborted the whole run, including the valid primary comparison.

## Behaviour to implement
1. **Primary policy unchanged.** A tie under the PRIMARY `altloc: highest_occupancy` keeps its explicit refusal (same text, exit code 2, no outputs).
2. **Optional sensitivity only.** If only the optional calculation meets a tie:
   - the primary comparison is completed and written exactly as without the optional entry (same windows, residues, exclusions, summary values);
   - the sensitivity is reported as **unavailable** (never zero, never an alternate chosen silently, never a partial number):
     `summary.json` gets `altloc_sensitivity = {policy, status: "unavailable: alternate locations tied on occupancy", available: false, values_reported: false, reason, why_not_partial, ties: [...]}`
     with no `max_abs_change` key;
   - each tie lists structure (A or B), chain, residue number, insertion code, whether the residue is in the correspondence, the alternate-location identifiers and their occupancies;
   - the report states "UNAVAILABLE, not zero" and lists the ties; the viewer shows a warning with the same words; `manifest.json` records the outcome in `policies.altloc_sensitivity_outcome`.
3. **No partial sensitivity.** A partial value would need a new rule for the tied residue (choose or skip it) and would change windows around it and the global fit
   (the fit uses all positions), so the whole optional sensitivity is unavailable. Ties at residues outside the correspondence also make the optional policy refuse in the
   adapter (it resolves every residue of the selection); they are listed with `in_correspondence: false` and the sensitivity is unavailable for them too (no new selection rule).
4. **Reliable structured signal.** The adapter raises a subclass `AltlocTieError(AdapterError)` carrying the list of all ties of the structure (all residues, collected before raising; the first tie's
   message text is unchanged). The adapter command-line writes them to a JSON file given by the new optional `--error-json PATH` and still prints the same `rejected: ...` line and exits 2.
   structpair passes `--error-json` and converts only a present, well-formed `altloc_tie` record into its exception; any other failure (parsing, mapping, numerical, file, unrelated adapter refusals) is
   re-raised unchanged. Both structures are examined for the optional policy so ties in A and B are both reported.
5. **Versions.** `STRUCTPAIR_VERSION` 0.3.0 → 0.3.1 (the summary gains the unavailable record; nothing else changes). Manifest schema `/3` unchanged (additive fields only; schemas /1 to /3 are read as before).
   Adapter version unchanged (successful extraction outputs are identical; the change affects only the refusal record and a new optional flag).
6. **Not changed:** the CSV number format (6 significant figures); the BTK-2 precision question is diagnosed separately in this directory.

## Tests (before the implementation is accepted)
SYNTHETIC cases: no alternate locations; unique highest-occupancy alternate (computed); tie in the optional policy; ties in both structures and several residues; tie in the primary policy (still fails); an unrelated
error (not swallowed); the primary outputs equal the run without the optional entry; JSON, report and viewer distinguish unavailable from zero; historical schemas still load.
