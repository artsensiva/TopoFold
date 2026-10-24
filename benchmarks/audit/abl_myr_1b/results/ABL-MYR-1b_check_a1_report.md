> Amendment a1 of the frozen check (check_amendment_a1.json): index files with repeated group names are accepted. The original attempt stopped at S3 and is kept as results/ABL-MYR-1b_check_results_attempt1.json.

# ABL-MYR-1b comparability check (generated from ABL-MYR-1b_check_results.json)

Post-hoc check of run inputs and parameter files; no trajectory opened. ABL-MYR-1 (G1 stop) and ABL-MYR-1a (V4) remain valid. Not evidence about Ayaz WT/M472I, HDX, drug response or biological function; no causal Myr claim.

**Outcome: not_comparable**

| Criterion | Pass |
|---|---|
| K1_usable_coordinates | False |
| K2_thermostat_comparable | True |
| K3_tc_groups | True |
| K4_composition | True |
| K5_replicates_not_copies | True |

R1: 0 of 40 run inputs usable and validated. Atom order equal to system_em.gro in 40 runs; protein bond lengths None–None Å; Cα RMSD to system_em.gro None–None Å (limit 1.0 Å); velocities present in 40 runs.

R2 production-start Cα RMSD after superposition (Å; min / median / max, pairs):

- Set1 within_Myr: 0.000 / 0.000 / 0.000 (0 pairs); heavy-atom median nan; distinct starts 10, distinct velocity sets 10
- Set1 within_noMyr: 0.000 / 0.000 / 0.000 (0 pairs); heavy-atom median nan; distinct starts 10, distinct velocity sets 10
- Set1 between_Myr_noMyr: 0.000 / 0.000 / 0.000 (0 pairs); heavy-atom median nan
- Set2 within_Myr: 0.000 / 0.000 / 0.000 (0 pairs); heavy-atom median nan; distinct starts 10, distinct velocity sets 10
- Set2 within_noMyr: 0.000 / 0.000 / 0.000 (0 pairs); heavy-atom median nan; distinct starts 10, distinct velocity sets 10
- Set2 between_Myr_noMyr: 0.000 / 0.000 / 0.000 (0 pairs); heavy-atom median nan
- Set1 vs Set2 Myr: Cα median nan (0 pairs)
- Set1 vs Set2 noMyr: Cα median nan (0 pairs)

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

