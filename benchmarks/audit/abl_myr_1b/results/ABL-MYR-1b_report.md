# ABL-MYR-1b results (generated from ABL-MYR-1b_results.json)

**Post-hoc, descriptive evidence from simulated Abl kinase-domain ensembles (de Buhr and Graeter, Dryad).** The study is separate from Ayaz WT/M472I, HDX, drug response and biological function. The two thermostat sets are a consistency check, not independent replication. No causal, clinical or novelty claim is made.

**Outcome: baselines cover TopoFold (post-hoc, descriptive)** (decision-rule step 4)

## Gates

- G0 archives: PCA.zip 875084667 B, kinase_kinked_Myr.zip 132587277 B, kinase_kinked_noMyr.zip 132251259 B; lineage hashes verified: True
- G1: 107 members extracted (293374746 B)
- G2: topology ref, 753 atoms, 251 residues, φ/ψ comparator available True
- G3: identity 1.0000, reference span [254, 504], unobserved ranges []
- G4b: G4b verifies recorded evidence. It does not establish identical production starts and does not read any Set2 coordinate file
- G5: 0 runs excluded; sets analysed {'Set1': True, 'Set2': True}
- G6: pbc used none; worst excluded-window fraction 0.00398406374501992

## αI helix

Residues 518-531 represented: False. alphaI is not (fully) represented: the analysis is limited to the represented kinase-domain residues and no claim is made about the alphaI helix or the Myr-binding region.

## Set1

Runs retained: Myr 10, noMyr 10; leave-one-run-out recomputations: 20.

| Channel | Role | Windows | E > 0 | Robust (full and every leave-one-run-out) | Median E |
|---|---|---|---|---|---|
| kappa | TopoFold | 249 | 50 | 21 | -0.0174 |
| tau | TopoFold | 248 | 43 | 16 | -0.0243 |
| d_i_i2 | Calpha baseline | 249 | 41 | 20 | -0.0170 |
| d_i_i3 | Calpha baseline | 248 | 50 | 30 | -0.0199 |
| rmsf | Calpha baseline | 251 | 41 | 10 | -0.0212 |
| phi | stronger-information comparator | 250 | 37 | 20 | -0.0120 |
| psi | stronger-information comparator | 250 | 46 | 15 | -0.0137 |

Run-to-run variability (median over windows of Jensen–Shannon distances) and temporal diagnostics:

- kappa: within Myr 0.1435117867253486, within noMyr 0.16443284798605026, between 0.15700155946093255; lag-1 median 0.13838159522893975, max 0.9276051749743265; decorrelated within 50 frames in 0.9843373493975903
- tau: within Myr 0.1965892137047051, within noMyr 0.23218503049253236, between 0.2322550740195195; lag-1 median 0.24834368233639684, max 0.9625526503581678; decorrelated within 50 frames in 0.9689516129032258
- d_i_i2: within Myr 0.1479891096915743, within noMyr 0.16937794505815051, between 0.162635518052117; lag-1 median 0.12773247748172048, max 0.9251379263370793; decorrelated within 50 frames in 0.9835341365461847
- d_i_i3: within Myr 0.1644686732051121, within noMyr 0.1892414294729251, between 0.1930391519955858; lag-1 median 0.18323664627445319, max 0.9588443408888417; decorrelated within 50 frames in 0.980241935483871

Production starts (a2 results; descriptive, Calpha RMSD to system_em.gro of 20 retained runs): 0.271–0.309 Å.

## Set2

Runs retained: Myr 10, noMyr 10; leave-one-run-out recomputations: 20.

| Channel | Role | Windows | E > 0 | Robust (full and every leave-one-run-out) | Median E |
|---|---|---|---|---|---|
| kappa | TopoFold | 249 | 44 | 20 | -0.0190 |
| tau | TopoFold | 248 | 46 | 17 | -0.0258 |
| d_i_i2 | Calpha baseline | 249 | 42 | 21 | -0.0186 |
| d_i_i3 | Calpha baseline | 248 | 31 | 16 | -0.0268 |
| rmsf | Calpha baseline | 251 | 54 | 20 | -0.0208 |
| phi | stronger-information comparator | 250 | 37 | 14 | -0.0155 |
| psi | stronger-information comparator | 250 | 37 | 18 | -0.0156 |

Run-to-run variability (median over windows of Jensen–Shannon distances) and temporal diagnostics:

- kappa: within Myr 0.15200247170360826, within noMyr 0.1694459081153691, between 0.16066746300395182; lag-1 median 0.1464221543677709, max 0.9326572656258533; decorrelated within 50 frames in 0.9791164658634538
- tau: within Myr 0.22389598560324045, within noMyr 0.2262391085739188, between 0.22337304264204857; lag-1 median 0.25901929882046454, max 0.975102591816197; decorrelated within 50 frames in 0.9647177419354839
- d_i_i2: within Myr 0.15815561051806487, within noMyr 0.1701975107586093, between 0.16855030856832; lag-1 median 0.13587236762100693, max 0.9313926713728912; decorrelated within 50 frames in 0.978714859437751
- d_i_i3: within Myr 0.1830639543538114, within noMyr 0.2070085749559029, between 0.19786654349729027; lag-1 median 0.19143565048572525, max 0.9826527606474251; decorrelated within 50 frames in 0.9745967741935484

Production starts (a2 results; descriptive, Calpha RMSD to system_em.gro of 20 retained runs): 0.265–0.303 Å.

## Outcome detail

- R_T (TopoFold residues): [23, 24, 25, 92, 93, 94, 104, 105, 106, 107, 108, 109, 110, 111, 191, 192, 193, 198, 199, 200, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237]
- R_B (baseline residues): [23, 24, 25, 63, 92, 93, 94, 95, 102, 105, 106, 107, 108, 109, 110, 111, 112, 191, 192, 193, 195, 196, 197, 198, 199, 200, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 234, 235, 236, 237]
- TopoFold residues outside N2(R_B): []
- baseline residues outside N2(R_T): [63]
- φ/ψ comparator residues (reported separately): [105, 106, 107, 109, 110, 194, 199, 200, 225, 226, 227, 228, 229, 230, 234]

These are descriptive filters (E > 0 in the full set and in every leave-one-run-out recomputation; the three-residue rule). They give no error-rate control or significance, and no outcome is evidence of the superiority of a method.

