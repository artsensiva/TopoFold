# Pre-registration: 1OPJ (chains A, B) vs 1OPK (chain A), structpair 0.3.0

Written 2026-10-25 before any structpair run on these structures. Inspection evidence:
`notes/pre_run_inspection.txt` (script `scripts/inspect_ca.py`), PDB header records of the two files.

## Inputs (RCSB, downloaded 2026-10-25)
| Entry | URL | SHA-256 |
|---|---|---|
| 1OPJ | https://files.rcsb.org/download/1OPJ.pdb | 6d3444a1f081b1e3b5079babd515d4e848586ec7055647fad9d294b9ee8fb7fb |
| 1OPK | https://files.rcsb.org/download/1OPK.pdb | 867ed2b25b80e0e2ba52f731e48df25fb8e982b052b89440a123d567c4384974 |
(mmCIF copies also kept for reference; `inputs/SHA256SUMS`.)

## What the deposition records say (both: Nagar et al., Cell 112:859, 2003, PMID 12654251)
| | 1OPJ | 1OPK |
|---|---|---|
| Construct | mouse Abl1 kinase domain, author residues 248–534 (+ GAMDPS cloning tag 242–247) | mouse Abl1 SH3-SH2-kinase, 46–534 (+ GAMDPS tag 40–45); SH3/SH2 present from 83 |
| Mutation | none recorded | D382N (SEQADV: "ENGINEERED MUTATION", UNP Asp 363) |
| Ligands | STI (imatinib) + MYR (myristic acid) + Cl⁻, in each chain | P16 (a pyrido[2,3-d]pyrimidinone inhibitor) + MYR + glycerol |
| Resolution / R / Rfree | 1.75 Å / 0.210 / 0.242 | 1.80 Å / 0.196 / 0.221 |
| Space group, copies | P1, 2 chains | C2221, 1 chain |
| Crystallization | 22% PEG 4000, 100 mM MES pH 5.6, 200 mM MgCl2, 277 K | 20% PEG 3350, 200 mM KNO3, pH 7.0, 293 K |
| Models / altlocs | 1 model, no altloc atoms | 1 model, no altloc atoms |

Both entries use the same author numbering (DBREF: author 248 = UNP 229 in 1OPJ; author 46 = UNP 27 in 1OPK, i.e. the same offset of 19).

## Scientific question
Where, and by how much, does the local Cα backbone geometry of the Abl kinase domain (248–529) differ between
1OPJ chain A and 1OPK chain A, and how do those differences compare with the descriptive copy-to-copy range between
1OPJ chains A and B? This is a descriptive geometric question. The two entries differ simultaneously in construct
(isolated kinase domain vs SH3-SH2-kinase), ATP-site inhibitor (imatinib vs P16), the D382N mutation, crystal
form and crystallization conditions, so no difference can be attributed to any single one of these factors.

## Residue correspondence
`same_numbering: "248-529"` (282 positions).
- Checked by: identical DBREF offsets; per-residue Cα residue names compared over the common Cα positions 244–529.
- 243/244–247 excluded on purpose: in 1OPJ they are the GAMDPS cloning artefact; in 1OPK they are native linker
  residues (IYGV), so these are not corresponding residues.
- Upper limit 529: last Cα present in 1OPJ chain A (530–534 are REMARK 465 missing there).
- Only remaining residue-name difference in 248–529: 382 (Asp in 1OPJ A/B, Asn in 1OPK A), the recorded D382N
  mutation → `resname_mismatch: "flag"`.

## Primary pairing and chains
Primary: 1OPJ chain A vs 1OPK chain A (first deposited chain of each; chosen a priori). 1OPJ chain B is included as a
second chain of entry A for within-entry crystal-copy variation. 1OPK has one chain, so no within-entry range exists
for that entry.

## Missing coordinates, residue names, alternate conformations
- Missing: none expected in 248–529 for any of the three chains (Cα present for all). `on_missing: "reject"`, so an
  unexpected absence stops the run instead of silently shrinking coverage.
- Residue names: `"flag"` (for D382N only, as above).
- Altlocs: no altloc atoms in either file. Policy `blank_then_A` for both entries (the guide's example);
  `altloc_sensitivity: "highest_occupancy"` set anyway, expected to show no change.

## Advisory minimum difference
Configured for `disp` only: `{"disp": 1.0}` Å, defaults otherwise (statistic `min`, within_factor 1.0).
Reason: 1 Å Cα displacement is a conventional screening level, several times the coordinate error expected at
1.75–1.8 Å resolution; it gives the advisory column something to assess so the feature is exercised. It is not a
significance threshold and not calibrated. No minimum is set for κ/τ/distance/RMSD channels because I have no
independent basis for a threshold in those units. These settings will not be changed after seeing results.
