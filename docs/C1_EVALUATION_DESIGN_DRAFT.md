# C1 Evaluation Design — DRAFT — not frozen; permission pending (2026-10-14, revision 2)

**Status:** unfrozen. C1 scoring is paused until the Rascore label-data terms are clarified
(`docs/C1_RASCORE_TERMS_REQUEST_DRAFT.md`). This draft uses only existing documents and the local metadata
`entry.tsv` (Rascore commit `8f77dc2`, SHA-256 `bac0fcfb…`).
- Counts: `benchmarks/audit/c1_feasibility/results/c1_isoform.json`, plus KRAS subsets computed from the same table.
- Nothing has been downloaded, implemented or scored.
- Revision 2: mutually exclusive outcomes including "materially inferior"; claim wording narrowed; weighting and
  bootstrap defined exactly; operational thresholds labelled as such; the meaning of "unliganded" clarified.

## 1. Candidate primary task

- **Population (candidate primary):** all KRAS X-ray chains labelled **Y32in.3P-ON** or **Y32out.2P-OFF**
  (switch-I, SW1).
- **Task:** for each query chain, retrieve the nearest eligible reference chain by SW1 local geometry and transfer its
  label (1-NN).
- **Claim boundary:** only recovery of an expert **structural** classification from Cα geometry. No claim about
  activity, signalling, drug response, pockets or clinical relevance. Nucleotide state, which co-defines the labels,
  is not analysed.
- **Excluded:**
  - GEF (absent in KRAS);
  - Outlier and Disordered chains;
  - non-X-ray structures;
  - chains with any missing Cα in residues 25–40.

| Subset | ON: chains / entries / crystal forms (largest CF) | OFF: chains / entries / crystal forms (largest CF) | Mixed CFs |
|---|---|---|---|
| All KRAS (candidate primary) | 133 / 52 / 21 (67) | 166 / 100 / 43 (19) | 3 |
| Excluding CF5 (sensitivity) | 66 / 38 / 20 (12) | 166 / 100 / 43 (19) | 3 |
| No pharmacological ligand (sensitivity) | 80 / 41 / 18 (31) | 68 / 36 / 18 (10) | 2 |

**CF5** is P 21 21 21: 67 ON chains in 14 entries (51 of them Q61H) and no OFF chains.

**"No pharmacological ligand"** means Rascore `pharmacological_class` = None. Nucleotides (GTP analogues, GDP), ions
and other non-pharmacological ligands remain present.

## 2. Query/reference construction (leakage control)

**Reference set per query:** all eligible chains except those sharing any of:
- (a) the query chain itself;
- (b) the same PDB entry;
- (c) the same Rascore crystal form;
- (d) the same named pharmacological ligand, when one is present. Four ligand names recur across crystal forms.

**Not controlled:**
- shared laboratories, constructs or deposition series. `entry.tsv` has no author or laboratory field.

Leave-one-crystal-form-out reduces leakage but does not establish independence.

**Use of labels:** labels only score transfer. They never build representations or select parameters.

## 3. Grouping, weighting and confounds

**Unit:** the crystal form (CF), kept for all KRAS, including CF5.
- **Within a CF,** chains are weighted equally. Entries are therefore **not** weighted equally: an entry with more
  chains in a CF contributes more.
- **Entry-equal weighting** (mean over chains within an entry, then entries within a CF) is reported as a sensitivity
  analysis.

**Predeclared sensitivity analyses (they never override the primary):**
- excluding CF5;
- no pharmacological ligand;
- entry-equal weighting;
- stratification by mutation group (WT, G12C, G12D, G12V, Q61H, other).

**Confounds visible in the metadata:**
- SP2-site ligands in 84/166 OFF chains and SP12-site ligands in 49/133 ON chains;
- a cysteine-light G12C construct in 59 OFF chains;
- Q61H dominating ON;
- CF5 holding half of the ON chains.

**Not resolvable from the metadata:**
- laboratory or deposition series;
- crystallisation conditions;
- crystal contacts at SW1;
- nucleotide state, which co-defines the labels (ON ≡ 3P, OFF ≡ 2P).

## 4. Representation and baselines (equal information: Cα of residues 25–40, UniProt P01116 numbering)

- **Tested — TopoFold "A":** √(½·mean Δκ² + ½·mean wrap(Δτ)²) over the 14 virtual bond angles (κ) and 13 virtual
  torsions (τ) of the 16-residue window. It is alignment-free, and it is the distance that met its criterion on the
  WP-D φ/ψ target.
- **Equal-information baselines:**
  - **B2:** Cα RMSD after optimal superposition, computed with an established library (CMP-1);
  - **B1:** Euclidean distance between internal Cα–Cα distance vectors (pairs |i−j| ≥ 2).
- **Stronger-information comparator, reported separately:** the Y32 CA–G12 CA distance, which adds residue 12 and is a
  proxy for the label rule.
- **Interpretation:**
  - A "superior" result means **better label recovery by A under this metric on this evaluation**. It does not mean
    that A contains information absent from the same coordinates.
  - "Non-inferior" means A is not materially worse than B2 here. It does **not** show useful absolute performance:
    B2 itself may be weak.
  - Absolute balanced accuracy with CIs is reported for A, B2, B1 and a label-shuffled reference, as context only.
  - No novelty claim: B2 and B1 are standard.

## 5. Metric, uncertainty and outcomes

**Primary metric — crystal-form-weighted balanced accuracy.**

For method m and class y ∈ {ON, OFF}:
- c_i^m ∈ {0, 1} is correct 1-NN transfer for query chain i.
- F_y is the set of CFs containing at least one class-y query chain.
- n_{f,y} is the number of class-y query chains in CF f.

The class recall is

  R_y^m = (1/|F_y|) · Σ_{f∈F_y} (1/n_{f,y}) · Σ_{i∈f, y_i=y} c_i^m

and BA^m = (R_ON^m + R_OFF^m)/2.

- Each CF contributes equally to the recall of every class it contains.
- A mixed CF contributes to both classes.
- **Ties:** smallest distance first, then lexicographic (PDB ID, chain).

**Estimand:** D = BA^A − BA^B2.

**Uncertainty — paired, stratified crystal-form bootstrap:**
- Strata: ON-only CFs, OFF-only CFs, and mixed CFs (3).
- Each stratum is resampled with replacement, keeping its size.
- A resampled CF brings **all** its query chains of both classes and its query results for A and B2 together, so
  mixed CFs are resampled jointly and the comparison stays paired.
- 10,000 resamples; percentile 95 % CI [L, U]; fixed seed.
- **Limitation:** queries are resampled with their realised nearest neighbours. Variability of the reference set itself
  is not captured.

**Margin:** δ = 0.05, **provisional**. In a two-class triage, a loss of 5 balanced-accuracy points is roughly one
additional misassigned query in twenty per class. That is judged material for a template-selection aid. It is a
judgement, not derived from the sample size or the precision below.

**Mutually exclusive outcomes**, applied in order:

| Outcome | Condition |
|---|---|
| **Materially inferior** | U < −δ |
| **Better label recovery** (under this metric on this evaluation) | L > 0 |
| **Non-inferior** | −δ < L ≤ 0 |
| **Inconclusive** | L ≤ −δ ≤ U |

**Precision (illustrative only).** Treating CFs as independent (ON 21, OFF 43) with a paired disagreement rate p_d, a
normal approximation 1.96·√((p_d/21 + p_d/43)/4) gives half-widths of 0.083, 0.117 and 0.143 for p_d = 0.1, 0.2 and
0.3.
- This is not a power calculation.
- Dependence among CFs, the CF5 concentration and the bootstrap's realised behaviour can all differ.
- Under those illustrative numbers, a non-inferiority outcome would require D to be observed above zero by roughly the
  half-width minus δ. A non-negligible chance of "inconclusive" should be expected.

## 6. Minimum adapter, proceed and stop

**Adapter** (proposed, not built): (PDB ID, chain, UniProt residues 25–40, plus 12 for the comparator) → owned Cα
array and a completeness flag.
- One established mmCIF parser, with SIFTS residue mapping.
- Altloc policy: highest occupancy, then "A".
- Records file hashes.
- No general loader and no new kernel.

**Proceed only if:**
1. Written permission covers internal evaluation and publication of metrics.
2. This design is frozen.
3. After mapping and exclusions, each class retains at least 15 CFs.
4. The residue mapping reproduces the `entry.tsv` sequence for every included chain.

**Stop if:**
- permission is denied or unclear;
- either class falls below 15 CFs;
- the mapping fails for more than 5 % of chains;
- label-defining information (Y32 OH, nucleotide) is found in the tested inputs.

The 15-CF minimum and the 5 % mapping-failure cutoff are **proposed operational choices** for feasibility. They are
not statistical guarantees: meeting them does not ensure a narrow or meaningful CI.

## Open items before freezing (after the permission response)

- Confirm or revise the provisional δ = 0.05.
- Confirm the bootstrap settings and the tie rule above.
