# Figure: 1OPJ chains A, B versus 1OPK chain A (descriptive; structpair 0.3.0)

![chain variation](chain_variation.png)

**Purpose.** To show the partial status of within-entry crystal-copy variation on a real pair: two chains in one entry, one in the other.

**What it shows.** Per mapped residue (282 positions: 248-529), the Cα displacement between 1OPJ and 1OPK after a trimmed fit:
the min–max band and median over the 2 cross pairings, and the largest difference between the two 1OPJ chains (dashed). **1OPK has a single
chain, so its within-entry variation is unavailable and appears as missing, not zero** (empty cells; hatched legend entry). No advisory rule is
configured, so no flag is assessed. For orientation only: the median over positions of the between-entry displacement is 0.46 Å and its largest
value is 33.1 Å; positions above 5 Å are 248-251,297,401-415 (PDB numbering); the largest
difference between the two 1OPJ chains is 5.3 Å, above 2 Å at 248-251,293-295.

**Main limitation: nothing here can be attributed to imatinib, nucleotide state, mutation or resistance.** From the deposited records (cross-checked
with RCSB): both entries are *Mus musculus* Abl1 (UniProt P00520, not human), expressed in insect cells, with myristate bound in each.
They differ in the ligands (1OPJ: CL (CHLORIDE ION), MYR (MYRISTIC ACID), STI (HETSYN STI-571; IMATINIB); 1OPK: MYR (MYRISTIC ACID), P16 (HETSYN PD166326), GOL (HETSYN GLYCERIN; PROPANE-1,2,3-TRIOL)), in the construct (kinase domain, 229–515 in
1OPJ; SH3-SH2-kinase domain, 27–515 in 1OPK, database numbering), in an engineered mutation (D363N, PDB residue 382, in 1OPK; flagged in
the output), in crystal form (1OPJ P 1, 2 chains; 1OPK C 2 2 21, 1 chain) and in resolution (1.75 and 1.8 Å). Neither entry
contains a nucleotide. Chains in one crystal are not independent experimental replicates, and crystal packing can affect each of them; the dashed
line is a descriptive range, not a noise floor, an uncertainty or an HDX replicate estimate.

**Missing residues (REMARK 465, inside the shared numbering 248–534).** 1OPJ chain A: 530-534; chain B: 532-534; 1OPK chain A: 532-534. The correspondence keeps only residues
observed in all three chains, so it ends at 529.

**Source data.** `chain_variation.csv`, `chains.csv`, `manifest.json` and `../abl_demo_record.json` (inputs and their SHA-256, header facts, checks).

**Reproduce.** `.venv/bin/python benchmarks/audit/structpair_abl_demo_2026-10-25/run_demo.py` (15 CPU-s, 7 s wall).
