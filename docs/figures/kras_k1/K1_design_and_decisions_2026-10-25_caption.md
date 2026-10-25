# Figure: KRAS K1 design and decision branches (frozen protocol v1, 2026-10-25)

![K1 design and decisions](K1_design_and_decisions_2026-10-25.png)

**Purpose.** To show, in one view, how the frozen KRAS K1 study runs and which outcomes would justify further development
or close the application, so the decision rule is visible before any simulation data are scored.

**What it shows.** Four steps: an endpoint computed from two published crystal structures only; a gate on the spread of
that endpoint; the fetch and gating of the simulation files (owner approval required); and a comparison of TopoFold with
Cα baselines at an equal budget of residues, where the TopoFold and baseline channels (T\*, B\*) are chosen once from the
primary analysis and stay fixed in every sensitivity check. Four outcomes follow, evaluated in the order C1, C2, C3, C4.
Only C3 (green) justifies further development, and then only a separately frozen, prospective evaluation on another
system; it is not functional validation. The grey outcomes and the gates close the application without a search for a
favourable result.

**Main limitation.** It is a schematic of a design, not a result: nothing was scored. The gate tests only whether the
crystal proxy has a spread-out comparison target; it does not show that the endpoint can discriminate methods. The
parameters (0.15, 0.5, 0.8, 90 %, the window and the run minimums) are design choices fixed before any data were read; they
are not calibrated, and no significance is claimed. The crystal pair was already analysed during development, so it is
not held out.

**Evidence class.** Frozen design (schematic). The only recorded facts in it are the run counts, the file count, the file
size and the atom count, read from the provenance checks.

**Source data.** `benchmarks/audit/kras_k1/protocol_params.json`, `benchmarks/audit/kras_k1/results/k1_provenance_checks.json`
and `benchmarks/audit/next_app_feasibility/results/feasibility_metadata.json`. SHA-256 of these inputs, of the generator and
of the PNG are in `K1_design_and_decisions_2026-10-25_provenance.json`.

**Reproduce.** `PYTHONPATH=python .venv/bin/python benchmarks/audit/kras_k1/plot_design_schematic.py` (about 1.7
CPU-seconds, matplotlib only). External material: none. This figure replaces the draft schematic of commit 2b75a0d,
which is kept in Git history.
