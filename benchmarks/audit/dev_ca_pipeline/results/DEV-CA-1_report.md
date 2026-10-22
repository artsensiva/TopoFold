# Cα feature inspection report: DEV-CA-1

Generated from the exported results files. Engineering output: descriptive geometry only; no state, pocket, allostery or biological inference.

## Input
- Source files: `ORADD_open_state.gro` (SHA-256 25c90e78f549d97e…), `ORADD_open_state_traj.xtc` (SHA-256 39e77b58ded7279d…)
- Reader: XTCReader; units: coordinates angstrom (source unit nm; converted to Å by MDAnalysis (convert_units=True)); frame identity: not_verifiable
- Frames: 751 (source frame indices 0–2250); residues: 214 (214 with Cα, 0 without)
- Breaks: 0 ; segments: 1
- Altloc decisions: 0; PBC policy: make_whole_sequential (frames changed: 0)

## Feature validity (cells = frames × positions)

| Reason | κ | τ |
|---|---|---|
| valid | 159212 | 158461 |
| missing_coordinates | 0 | 0 |
| terminal_or_break | 1502 | 1502 |
| degenerate_segment | 0 | 0 |
| near_collinear | 0 | 0 |
| run_too_short_for_core | 0 | 0 |
| core_error | 0 | 0 |
| non_finite_core_output | 0 | 0 |

- Non-finite values among valid cells: 0
- Core errors: 0; max float32 rounding of core input: 0 Å
- κ over valid cells (rad): min 0.2819, median 1.4761, max 1.8334
- τ over valid cells (rad): min -3.1416, max 3.1414; circular mean resultant length 0.2899

## Positions never valid in any frame

- κ (2): SYSTEM:1 [terminal_or_break], SYSTEM:214 [terminal_or_break]
- τ hinges (2): SYSTEM:1→SYSTEM:2 [terminal_or_break], SYSTEM:213→SYSTEM:214 [terminal_or_break]

Features file: `DEV-CA-1_features.npz` (SHA-256 9e91d69fdb591c0e…).
