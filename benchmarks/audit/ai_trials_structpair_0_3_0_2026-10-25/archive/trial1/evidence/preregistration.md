# Pre-result analysis decisions
Frozen at 2026-10-25T09:58:01.557863+00:00. No structpair analysis has run; no shipped Abl spec, result, demonstration or internal trial report has been opened. Permitted README and guide contain background summaries, but were not used to choose regions.

Question: Where does the local C-alpha geometry of the shared Abl kinase sequence differ between these deposited structures, and how does that difference compare descriptively with the two crystal copies in 1OPJ? This comparison cannot isolate a drug or mutation effect.

Primary pair: first deposited protein chain A in each entry, selected by order, before results. Include both A and B of 1OPJ and A of 1OPK; retain every resulting pairing.

Correspondence: author numbers 248 through 534 inclusive, 287 positions. This is the full shared native sequence according to DBREF: both entries map to mouse ABL1 P00520 with author number minus 19. 1OPJ's 242-247 cloning tag is excluded by construct provenance, not by outcome. 1OPK-only N-terminal domains are outside shared coverage. SEQRES and observed ATOM CA names were cross-checked independently with audit_inputs.py. All 287 shared sequence positions agree except 382 (ASP in J, ASN in K), consistent with the SEQADV engineered mutation at reference 363. J A/B agree throughout. No insertion codes occur. See input-audit.json and pre-run-correspondence.csv.

Missing coordinates: include the full shared sequence and use on_missing=exclude so absent tail coordinates remain explicit missing rows. J A lacks 530-534, J B and K A lack 532-534. Do not impute coordinates or prune tail after results. Residue names: flag the known D382N difference and all affected windows, do not reject or silently ignore it. Altloc: blank_then_A, plus highest_occupancy sensitivity. No CA alternate locations were found; side-chain alternates do not define a CA choice.

Advisory: configure no positive minimum difference (empty dictionary), so all flags should remain unassessed. There is no independently justified practical-effect threshold for this descriptive comparison. Use the documented 0.2 displacement trim fraction, min statistic and within_factor=1.0; these do not establish significance.

Confounds known before results: 1OPJ is an isolated kinase-domain construct (293 deposited residues including tag); 1OPK is SH3-SH2-kinase (495 including tag), with D382N. Both contain MYR; J contains STI and chloride, K contains P16 and glycerol. Resolutions 1.75 and 1.80 A; space groups P1 and C2221. Crystallization differs (J pH 5.6, 277 K, PEG4000/MES/MgCl2; K pH 7.0, 293 K, PEG3350/KNO3). Thus construct, ligand, mutation and crystal environment co-vary. The two J chains are crystal copies, not independent experimental replicates; K has no within-entry pair.

Sources: https://www.rcsb.org/structure/1OPJ and https://www.rcsb.org/structure/1OPK; coordinate URLs and hashes in logs/download-abl.log. No external assistance.
