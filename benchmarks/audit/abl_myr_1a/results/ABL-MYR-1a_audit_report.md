# ABL-MYR-1a initialization audit (generated from ABL-MYR-1a_audit_results.json)

Post-hoc amendment; not the original frozen test. Engineering audit of starting files only: no trajectory, no descriptor, no condition difference. Not evidence about Ayaz WT/M472I, HDX, drug response or biological function.

**Verdict: V4_set1_comparability_failed**

| Comparison | Category | Protein atoms | max abs dx (nm) | RMSD (nm) |
|---|---|---|---|---|
| Q1 Set1 system.gro, Myr vs noMyr | identical | 4473 | 0.0 | 0.0 |
| Q2 Set2 system_em.gro, Myr vs noMyr | differs_large | 4473 | 11.087 | 0.288674293 |
| Q3 Set1 system_em.gro, Myr vs noMyr (calibration) | differs_large | 4473 | 11.087 | 0.288674293 |
| Q4 system_em.gro Set1 vs Set2, Myr | identical | 4473 | 0.0 | 0.0 |
| Q4 system_em.gro Set1 vs Set2, noMyr | identical | 4473 | 0.0 | 0.0 |
| Set1 Myr: system.gro vs system_em.gro | differs_large | 4473 | 7.489 | 10.096414298 |
| Set1 noMyr: system.gro vs system_em.gro | differs_large | 4473 | 7.513 | 10.10817115 |

- Set1: production_md.mdp identical between conditions: False; topol.top identical after removing Myr lines: False
  - differing: `Myr:compressed-x-grps=Protein_MYRP_ATP_MG | noMyr:compressed-x-grps=Protein_ATP_MG`
  - differing: `Myr:tc_grps=Protein_MYRP_ATP_MGWater_NA_CL | noMyr:tc_grps=Protein_ATP_MGWater_NA_CL`
  - differing: `Myr:SOL 29949 | noMyr:SOL 29938`
- Set2: production_md.mdp identical between conditions: False; topol.top identical after removing Myr lines: False
  - differing: `Myr:compressed-x-grps=Protein_MYRP_ATP_MG | noMyr:compressed-x-grps=Protein_ATP_MG`
  - differing: `Myr:tc_grps=Protein_MYRP_ATP_MGWater_NA_CL | noMyr:tc_grps=Protein_ATP_MGWater_NA_CL`
  - differing: `Myr:SOL 29949 | noMyr:SOL 29938`

Run inputs (md.tpr, MDAnalysis topology reader only):

- Set1 Myr: 10 of 10 readable; distinct protein digests 1; positions available: True
- Set1 noMyr: 10 of 10 readable; distinct protein digests 1; positions available: True
- Set2 Myr: 10 of 10 readable; distinct protein digests 1; positions available: True
- Set2 noMyr: 10 of 10 readable; distinct protein digests 1; positions available: True

Initialization stated in production_md.mdp (stated, not verified):

- Set1 Myr: integrator=md, tcoupl=V-rescale, ref_t=300	300, tau_t=1.0	1.0, pcoupl=Parrinello-Rahman, dt=0.002, nsteps=250000000, gen_vel=no, continuation=yes, constraints=h-bonds, nstxout_compressed=20000, compressed_x_grps=Protein_MYRP_ATP_MG
- Set1 noMyr: integrator=md, tcoupl=V-rescale, ref_t=300	300, tau_t=1.0	1.0, pcoupl=Parrinello-Rahman, dt=0.002, nsteps=250000000, gen_vel=no, continuation=yes, constraints=h-bonds, nstxout_compressed=20000, compressed_x_grps=Protein_ATP_MG
- Set2 Myr: integrator=md, tcoupl=Nose-Hoover, ref_t=300	300, tau_t=1.0	1.0, pcoupl=Parrinello-Rahman, dt=0.002, nsteps=250000000, gen_vel=no, continuation=yes, constraints=h-bonds, nstxout_compressed=20000, compressed_x_grps=Protein_MYRP_ATP_MG
- Set2 noMyr: integrator=md, tcoupl=Nose-Hoover, ref_t=300	300, tau_t=1.0	1.0, pcoupl=Parrinello-Rahman, dt=0.002, nsteps=250000000, gen_vel=no, continuation=yes, constraints=h-bonds, nstxout_compressed=20000, compressed_x_grps=Protein_ATP_MG

Derived: Set1 original checks pass = False; Set2 parameters pass = False; Set2 start direct evidence = True.

