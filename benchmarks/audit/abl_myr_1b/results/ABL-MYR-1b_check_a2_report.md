> Amendment a2 of the ABL-MYR-1b check (check_amendment_a2.json): the stored box, positions and velocities of the 40 production md.tpr are read independently with GROMACS and compared with the earlier reader. One bounded attempt; criteria K1-K5 and their thresholds are those of the frozen check. 'not_comparable' means criteria not met, not that the systems are physically incomparable; 'comparable_for_1b' means criteria passed, not that no confound exists.

# ABL-MYR-1b comparability check (generated from ABL-MYR-1b_check_results.json)

Post-hoc check of run inputs and parameter files; no trajectory opened. ABL-MYR-1 (G1 stop) and ABL-MYR-1a (V4) remain valid. Not evidence about Ayaz WT/M472I, HDX, drug response or biological function; no causal Myr claim.

**Outcome: comparable_for_1b**

| Criterion | Pass |
|---|---|
| K1_usable_coordinates | True |
| K2_thermostat_comparable | True |
| K3_tc_groups | True |
| K4_composition | True |
| K5_replicates_not_copies | True |

R1: 40 of 40 run inputs usable and validated. Atom order equal to system_em.gro in 40 runs; protein bond lengths 0.9475–1.9623 Å; Cα RMSD to system_em.gro 0.265288–0.309325 Å (limit 1.0 Å); velocities present in 40 runs.

R2 production-start Cα RMSD after superposition (Å; min / median / max, pairs):

- Set1 within_Myr: 0.296 / 0.323 / 0.355 (45 pairs); heavy-atom median 0.472; distinct starts 10, distinct velocity sets 10
- Set1 within_noMyr: 0.309 / 0.331 / 0.357 (45 pairs); heavy-atom median 0.478; distinct starts 10, distinct velocity sets 10
- Set1 between_Myr_noMyr: 0.319 / 0.341 / 0.372 (100 pairs); heavy-atom median 0.502
- Set2 within_Myr: 0.298 / 0.316 / 0.338 (45 pairs); heavy-atom median 0.467; distinct starts 10, distinct velocity sets 10
- Set2 within_noMyr: 0.299 / 0.335 / 0.363 (45 pairs); heavy-atom median 0.485; distinct starts 10, distinct velocity sets 10
- Set2 between_Myr_noMyr: 0.302 / 0.335 / 0.358 (100 pairs); heavy-atom median 0.498
- Set1 vs Set2 Myr: Cα median 0.317 (100 pairs)
- Set1 vs Set2 noMyr: Cα median 0.332 (100 pairs)

R3 tc groups:

- Set1 Myr: ['Protein_MYRP_ATP_MG', 'Water_NA_CL'] ok=True; Protein_MYRP_ATP_MG: 4573 atoms, 280 protein residues + {'ATP': 1, 'MG': 1, 'MYRP': 1}; Water_NA_CL: 90034 atoms, {'CL': 87, 'NA': 100, 'SOL': 29949}
- Set1 noMyr: ['Protein_ATP_MG', 'Water_NA_CL'] ok=True; Protein_ATP_MG: 4529 atoms, 280 protein residues + {'ATP': 1, 'MG': 1}; Water_NA_CL: 90001 atoms, {'CL': 87, 'NA': 100, 'SOL': 29938}
- Set2 Myr: ['Protein_MYRP_ATP_MG', 'Water_NA_CL'] ok=True; Protein_MYRP_ATP_MG: 4573 atoms, 280 protein residues + {'ATP': 1, 'MG': 1, 'MYRP': 1}; Water_NA_CL: 90034 atoms, {'CL': 87, 'NA': 100, 'SOL': 29949}
- Set2 noMyr: ['Protein_ATP_MG', 'Water_NA_CL'] ok=True; Protein_ATP_MG: 4529 atoms, 280 protein residues + {'ATP': 1, 'MG': 1}; Water_NA_CL: 90001 atoms, {'CL': 87, 'NA': 100, 'SOL': 29938}

R4 thermostat and coupling (differences between conditions):

- Set1 nvt.mdp: tcoupl V-rescale, ref_t 300 300, tau_t 0.5 0.5, pcoupl no, gen_vel yes, gen_seed -1, define -DPOSRES; differences tc_grps (inherent_group_name)
- Set1 npt.mdp: tcoupl V-rescale, ref_t 300 300, tau_t 0.5 0.5, pcoupl berendsen, gen_vel no, gen_seed None, define -DPOSRES; differences tc_grps (inherent_group_name)
- Set1 production_md.mdp: tcoupl V-rescale, ref_t 300 300, tau_t 1.0 1.0, pcoupl Parrinello-Rahman, gen_vel no, gen_seed None, define None; differences compressed_x_grps (inherent_group_name), tc_grps (inherent_group_name)
- Set2 nvt.mdp: tcoupl V-rescale, ref_t 300 300, tau_t 0.5 0.5, pcoupl no, gen_vel yes, gen_seed -1, define -DPOSRES; differences tc_grps (inherent_group_name)
- Set2 npt.mdp: tcoupl V-rescale, ref_t 300 300, tau_t 0.5 0.5, pcoupl berendsen, gen_vel no, gen_seed None, define -DPOSRES; differences tc_grps (inherent_group_name)
- Set2 production_md.mdp: tcoupl Nose-Hoover, ref_t 300 300, tau_t 1.0 1.0, pcoupl Parrinello-Rahman, gen_vel no, gen_seed None, define None; differences compressed_x_grps (inherent_group_name), tc_grps (inherent_group_name)

R5 composition:

- Set1: differences {'MYRP': [1, 0], 'SOL': [29949, 29938]}; net charge {'Myr': 4e-06, 'noMyr': 4e-06}; inherent ['MYRP present only with Myr', 'SOL differs by 11 (0.0367% of noMyr)']; potential confounds []
- Set2: differences {'MYRP': [1, 0], 'SOL': [29949, 29938]}; net charge {'Myr': 4e-06, 'noMyr': 4e-06}; inherent ['MYRP present only with Myr', 'SOL differs by 11 (0.0367% of noMyr)']; potential confounds []

R6 classification:

- inherent_to_Myr_model: ['Set1 nvt.mdp tc_grps: Protein_MYRP_ATP_MG Water_NA_CL | Protein_ATP_MG Water_NA_CL', 'Set1 npt.mdp tc_grps: Protein_MYRP_ATP_MG Water_NA_CL | Protein_ATP_MG Water_NA_CL', 'Set1 production_md.mdp compressed_x_grps: Protein_MYRP_ATP_MG | Protein_ATP_MG', 'Set1 production_md.mdp tc_grps: Protein_MYRP_ATP_MG Water_NA_CL | Protein_ATP_MG Water_NA_CL', 'Set1 MYRP present only with Myr', 'Set1 SOL differs by 11 (0.0367% of noMyr)', 'Set2 nvt.mdp tc_grps: Protein_MYRP_ATP_MG Water_NA_CL | Protein_ATP_MG Water_NA_CL', 'Set2 npt.mdp tc_grps: Protein_MYRP_ATP_MG Water_NA_CL | Protein_ATP_MG Water_NA_CL', 'Set2 production_md.mdp compressed_x_grps: Protein_MYRP_ATP_MG | Protein_ATP_MG', 'Set2 production_md.mdp tc_grps: Protein_MYRP_ATP_MG Water_NA_CL | Protein_ATP_MG Water_NA_CL', 'Set2 MYRP present only with Myr', 'Set2 SOL differs by 11 (0.0367% of noMyr)']
- potential_confound: []
- not_separable_from_Myr_presence: ['Set1: production-start geometry differences between Myr and noMyr (see R2 between_Myr_noMyr)', 'Set2: production-start geometry differences between Myr and noMyr (see R2 between_Myr_noMyr)']


## a2: independent reader

- GROMACS 2024.5-conda_forge (mixed); matches the freeze: True
- reader agreement (validation d): 40 of 40 run inputs; max |dx| 9e-07 nm (tolerance 0.001), max |dv| 2e-07 nm/ps (tolerance 0.001); velocities in both readers in 40 runs; GROMACS box vs MDAnalysis extract_box_info max 4.9976e-05 nm (tolerance 0.0001)
- descriptive Myr vs noMyr differences read from the run inputs: {"Set1": {"ld-seed": {"Myr": ["\"-1324034442\"", "\"-16970875\"", "\"-555771257\"", "\"-727316524\"", "\"1237870838\"", "\"1587374816\"", "\"1978182650\"", "\"447128343\"", "\"450394376\"", "\"460642464\""], "noMyr": ["\"-1303091211\"", "\"-1663619742\"", "\"-30268644\"", "\"-528363102\"", "\"1463647475\"", "\"16510209\"", "\"1965847268\"", "\"37816963\"", "\"675246710\"", "\"897391451\""]}, "molblocks": {"Myr": ["[{\"molecules\": 1, \"moltype\": 0, \"name\": \"Protein\"}, {\"molecules\": 1, \"moltype\": 1, \"name\": \"MYR\"}, {\"molecules\": 1, \"moltype\": 2, \"name\": \"ATP\"}, {\"molecules\": 1, \"moltype\": 3, \"name\": \"MG_ion\"}, {\"molecules\": 29949, \"moltype\": 4, \"name\": \"SOL\"}, {\"molecules\": 100, \"moltype\": 5, \"name\": \"NA\"}, {\"molecules\": 87, \"moltype\": 6, \"name\": \"CL\"}]"], "noMyr": ["[{\"molecules\": 1, \"moltype\": 0, \"name\": \"Protein\"}, {\"molecules\": 1, \"moltype\": 1, \"name\": \"ATP\"}, {\"molecules\": 1, \"moltype\": 2, \"name\": \"MG_ion\"}, {\"molecules\": 29938, \"moltype\": 3, \"name\": \"SOL\"}, {\"molecules\": 100, \"moltype\": 4, \"name\": \"NA\"}, {\"molecules\": 87, \"moltype\": 5, \"name\": \"CL\"}]"]}, "nrdf": {"Myr": ["\"11457.8      180252\""], "noMyr": ["\"11353.8      180186\""]}}, "Set2": {"ld-seed": {"Myr": ["\"-1029611126\"", "\"-1112234001\"", "\"-1328640691\"", "\"-1412393019\"", "\"-446199120\"", "\"-585702785\"", "\"1713217542\"", "\"1772795873\"", "\"344649271\"", "\"671773923\""], "noMyr": ["\"-1994288323\"", "\"-2119429686\"", "\"-266825421\"", "\"-715025050\"", "\"1326919954\"", "\"1523986196\"", "\"579034685\"", "\"672904591\"", "\"813235772\"", "\"882845639\""]}, "molblocks": {"Myr": ["[{\"molecules\": 1, \"moltype\": 0, \"name\": \"Protein\"}, {\"molecules\": 1, \"moltype\": 1, \"name\": \"MYR\"}, {\"molecules\": 1, \"moltype\": 2, \"name\": \"ATP\"}, {\"molecules\": 1, \"moltype\": 3, \"name\": \"MG_ion\"}, {\"molecules\": 29949, \"moltype\": 4, \"name\": \"SOL\"}, {\"molecules\": 100, \"moltype\": 5, \"name\": \"NA\"}, {\"molecules\": 87, \"moltype\": 6, \"name\": \"CL\"}]"], "noMyr": ["[{\"molecules\": 1, \"moltype\": 0, \"name\": \"Protein\"}, {\"molecules\": 1, \"moltype\": 1, \"name\": \"ATP\"}, {\"molecules\": 1, \"moltype\": 2, \"name\": \"MG_ion\"}, {\"molecules\": 29938, \"moltype\": 3, \"name\": \"SOL\"}, {\"molecules\": 100, \"moltype\": 4, \"name\": \"NA\"}, {\"molecules\": 87, \"moltype\": 5, \"name\": \"CL\"}]"]}, "nrdf": {"Myr": ["\"11457.8      180252\""], "noMyr": ["\"11353.8      180186\""]}}}

