# Explicit support for documented phosphorylated amino acids in the PDB-to-structpair workflow (KNOWN_ISSUES m14): design, written before the implementation (2026-10-26)

**Engineering change with regression tests; not a new evaluation.** KIN-PAIR-1 and its frozen outcome are unchanged; its PIM1 pairs are replayed only as an engineering regression
with new, derived specifications. Based on the local inventory (D84): the adapter's `protein` selection does not contain SEP, TPO or PTR, so observed, annotated, peptide-bonded
phosphoresidues were silently absent from the output.

## Scope
- Opt-in policy `modified_residues: "phospho_pdb"` (spec) / `--modified-residues phospho_pdb` (adapter); default `"none"` keeps every historical behaviour and output.
- Residue names **SEP, TPO, PTR only**, with parents **SER, THR, TYR** (checked against the wwPDB chemical component records: `results/component_parents_verified.json`). No other modified residue is admitted.
- **PDB-format files only** (topology file suffix .pdb or .ent, with MODRES records read from the file). Any other input format is refused under the policy. No claim for other formats or residues.
- Scope of the checks: residues in the chains the correspondence names (all selected chains if the adapter is run without a correspondence).

## What the current Cα adapter checks (documented, not extended for historical inputs)
Unique (chain, resid, icode) keys; one Cα per residue under the altloc policy; break records for a missing Cα, a chain/segid change, a numbering discontinuity
(next resid without insertion code, or next insertion code of the same resid), and a Cα–Cα distance above 4.2 Å in any frame. **It does not check peptide bonds** and never reads N or C atoms.
Unchanged by this task for every residue the policy does not touch.

## Admission rule for a residue named SEP, TPO or PTR (only under the policy)
1. Not admitted, recorded as excluded with its reason (never an error, never silently dropped): the residue lacks any of N, CA, C (a Cα-only input cannot verify polymer membership), or it is not peptide-bonded to a neighbour
   (C(i-1)–N(i) or C(i)–N(i+1) at most 2.0 Å, neighbours in the same chain in file order): a component such as a free ligand. A calcium ion (residue CA, atom CA) or any other name is never considered.
2. Refused with an explicit error naming the residue (the policy cannot document its identity): backbone and bonded, but **no MODRES record** for it; or a MODRES record whose residue name or **parent differs** from the table.
3. Otherwise admitted: its Cα enters the output with the same altloc policy as every residue, and the N and C atoms are chosen under that policy only for the bonding check.
The deposited residue name is **never renamed to the parent**; the parent, its source ("MODRES record, parent confirmed by the policy table"), the MODRES comment, the record type (ATOM or HETATM), chain, author number and insertion code are carried in the adapter report. Annotation is not evidence of functional state.

## Continuity
Admitted residues take part in the existing break rules unchanged (missing Cα, chain boundary, numbering adjacency, 4.2 Å). The N/C bonding check is an admission condition for these residues only; it is not applied to standard residues and does not alter any historical input.
Residue-name mismatch stays explicit: SEP versus SER is a name difference; `resname_mismatch: "flag"` marks every affected window and `"reject"` (the default) refuses; the documented parent supports the correspondence but does not hide the difference.

## Outputs when the policy is on (additive; nothing changes when it is off)
adapter report: per-residue `modification` and a `modified_residues` block (admitted, excluded with reasons); structpair: `residues.csv` columns `A_modification`, `B_modification`; `windows.csv` columns `modified_in_window3`, `modified_in_window4`;
`summary.json` block `modified_residues`; `manifest.json` `policies.modified_residues` and its outcome; the report section; `numeric_record.json` position fields; the viewer shows retained modifications and excluded cases.
Versions: structpair 0.4.0 → 0.5.0, adapter 0.1.0 → 0.2.0; manifest schema `/3` and numeric record schema `/1` unchanged (additive, opt-in); the policy is recorded in the manifest.

## Tests (before any real-data replay)
SYNTHETIC: correct SEP, TPO, PTR; missing Cα; a real chain break and a numbering gap; missing and contradictory annotation; insertion codes; alternate locations; a non-polymer SEP, a calcium ion; SEP versus SER handling; standard residues and default behaviour unchanged;
newly supported geometry against independent references and rigid-transform controls.

## Replay (engineering regression, after the tests)
New derived specifications for PIM1-1 and PIM1-2 add the omitted position; originals preserved. Expected window membership is computed independently from the revised correspondence and continuity rules and written **before** the structpair outputs are compared.
Unchanged windows are compared by complete residue identifiers; global fit and aggregates are explained; the four unaffected pairs are compared payload against payload (numeric record) separately from version and provenance changes. No tolerance changes.
