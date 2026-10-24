# KIN-ROB-1a: post-hoc amendment of KIN-ROB-1 (PKA numbering), 2026-10-24

This is a **post-hoc amendment** written after KIN-ROB-1 chain 1ATP E failed check C1. The reason is recorded in
`run_spec.json` under `amendment`: PKA role positions had been given in literature numbering but labelled as UniProt
positions.
- **Changes:** only the P05132 positions (+1) and the panel, which is restricted to 1ATP.
- **Runner:** byte-identical to `../kin_rob_1/run_kin_rob.py`.
- **Result** (`results/KIN-ROB-1_results.json`): pass.
- The original KIN-ROB-1 failure stands. See `../kin_rob_1/README.md`.
