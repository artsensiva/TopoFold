# Figure: KRAS K1 design and decision branches (DRAFT, 2026-10-25)

![K1 design and decisions](K1_design_and_decisions_DRAFT_2026-10-25.png)

**Purpose.** To show, in one view, how the proposed KRAS K1 study would run and which outcomes would justify further
development or close the application, so the decision rule is visible before any data are scored.

**What it shows.** Four steps: an endpoint computed from two crystal structures only; a gate on that endpoint; the fetch
and gating of the simulation files; and a comparison of TopoFold with Cα baselines at an equal budget of residues. Four
outcomes follow, evaluated in the order C1, C2, then C3 or C4. Only C3 (green) justifies further development, and then
only as a separately frozen confirmation on an independent system. The grey outcomes and the two gates close the
application without a search for a favourable result.

**Main limitation.** It is a schematic of a proposed design, not a result: nothing was scored. The parameters (0.15, 0.5,
0.8, 90 %, the window and the run minimums) are design choices fixed before any data are read; they are not calibrated,
and the protocol claims no significance.

**Evidence class.** Proposed design (schematic). The only recorded facts in it are the run counts, the file count, the file
size and the atom count, read from the provenance checks.

**Source data.** `benchmarks/audit/kras_k1/protocol_params_DRAFT.json`, `benchmarks/audit/kras_k1/results/k1_provenance_checks.json`
and `benchmarks/audit/next_app_feasibility/results/feasibility_metadata.json`. SHA-256 of these inputs, of the generator and
of the PNG are in `K1_design_and_decisions_DRAFT_2026-10-25_provenance.json`.

**Reproduce.** `PYTHONPATH=python .venv/bin/python benchmarks/audit/kras_k1/plot_design_schematic.py` (about 1.6
CPU-seconds, matplotlib only). External material: none.
