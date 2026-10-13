# WP-D Held-out Protocol (mdCATH): Selecting Cα-only Local-Conformation Components

**Status:** proposed final. It is frozen by `benchmarks/audit/wpd_local_retrieval/results/wpd_heldout_protocol_freeze.json`
before the domain draw. No trajectory has been downloaded and nothing has been scored. Downloads and scoring
require separate authorization.

**Relation to the AdK protocol.** `docs/WPD_LOCAL_RETRIEVAL_PROTOCOL.md` is frozen and unchanged. Its §7–§8
mdCATH text (marked provisional) is **superseded** by this document. In particular, this document replaces the
0.05-overlap margin and the 0.05 rad target floor. Target, methods, distances, tie rule and degeneracy rules
are inherited from the AdK protocol (§3, §4, §6, §10 A1–A3) unless stated otherwise here.

**Objective.** Select useful engine components honestly. M1/M2 (ordered τ), B1 (Cα–Cα distances), B2 (window
Cα RMSD) and B3 (κ block) are **competing candidates on equal terms**. C1 (scalar circular mean) is a
historical comparator. Nothing here assumes the ordered representation wins. Several methods may be selected,
or none.

---

## 1. Question and claim scope

For a residue window of W consecutive residues, how well does each Cα-only distance retrieve frames with
similar full-backbone local conformation (interior φ/ψ, AdK protocol §3)? This is tested **across replicas of
the same domain** at a fixed temperature, on CATH domains outside the AdK superfamily.

- **Can establish:** how far each method closes the gap between random and target-optimal retrieval, relative
  to the other methods, over a seeded sample of CATH superfamilies in mdCATH (CHARMM22*, TIP3P, 320 K and
  348 K).
- **Cannot establish:** performance on other force fields, other temperatures, multi-domain proteins, or for
  functional states; any claim about global conformational change; the value of a method for downstream
  tasks.

## 2. Data source (verified facts)

- **Dataset:** mdCATH (Mirarchi, Giorgino, De Fabritiis, *Sci. Data* 11, 1299 (2024),
  doi:10.1038/s41597-024-04140-z), CC BY 4.0.
  - Hugging Face dataset `compsciencelab/mdCATH`, pinned revision
    `5e3ed8aec62b689e01751db16275fdcdbc39e47f`.
  - 5,398 per-domain HDF5 files plus `mdCATH_domains.txt` and `mdcath_source.h5`.
- **Domain origin (paper):**
  - CATH release 4.2.0, S20 non-redundant set (14,433 domains);
  - restricted to 50–500 residues;
  - contiguous backbone only;
  - standard amino acids only;
  - termini capped (acetyl / N-methyl).
- **Simulation (paper):**
  - 20 ns NPT equilibration at 300 K;
  - the **final configuration of each system starts all 25 production runs** (5 temperatures × 5 replicas);
  - NVT production with a Langevin thermostat, 4 fs steps with hydrogen-mass repartitioning;
  - positions and forces every 1 ns.
- **File layout:** verified from the authors' public notebook `user/1_inspect_content.ipynb` and builder code
  at GitHub `compsciencelab/mdCATH` commit `71cbe93`; not from our own download.
  - per domain: groups `320`, `348`, `379`, `413`, `450`, each with replicas `0`–`4`;
  - per replica: `coords` (F × N × 3, Å), `forces`, `dssp`, `box` (3 × 3, nm, first frame), `rmsd`,
    `gyrationRadius`, `rmsf`, and attribute `numFrames`;
  - per domain: `resid`, `resname`, `element`, `z`, `pdbProteinAtoms`.
  - The atom set **includes hydrogens**: `z` contains 1 in the authors' example.
- **Builder behaviour that matters here:**
  - Coordinates are read from xtc files with no stated unwrapping, so periodic reconstruction is required.
  - Frames with heavy-atom RMSD to frame 0 ≥ 40 nm are dropped. This is practically never triggered, but
    the files store **no timestamps**, so a dropped frame cannot be detected; 1 ns spacing is assumed and
    disclosed.
- **Not used:** `mdcath_source.h5` holds per-replica dynamics summaries (radius-of-gyration range, DSSP
  composition). It is **not downloaded or consulted**, because selecting on it would inspect
  candidate-domain dynamics.
- **CATH metadata:** `cath-domain-list-v4_2_0.txt`, `cath-names-v4_2_0.txt` and the CLF 2.0 format README,
  from `download.cathdb.info/cath/releases/all-releases/v4_2_0/cath-classification-data/`. Hashes are
  recorded in the freeze file.
- **AdK in CATH 4.2.0:** the E. coli AdK chains 4AKE A/B and 1AKE A/B are single domains (`4akeA00`, …) in
  superfamily **3.40.50.300** ("P-loop containing nucleotide triphosphate hydrolases" in
  `cath-names-v4_2_0.txt`).

## 3. Development data and exclusions (before the draw)

The following domains are **development data**. Their superfamilies are removed from the held-out pool:

| Item | Superfamily | Reason |
|------|-------------|--------|
| AdK (WP-B, WP-D AdK tier) | 3.40.50.300 | Protocol developed and inspected on AdK |
| `1r9lA02` | 3.40.190.100 | Example domain in the authors' notebook; its frame counts were read during schema verification |
| Schema sample (§3.1) | its CATH superfamily | Downloaded first, used to develop and test the loader |

### 3.1 Schema-inspection sample

- **Choice rule (deterministic, metadata only):** the metadata-eligible domain (§4, criteria 1–3) with the
  smallest file size in bytes; ties are broken by domain name. It is fixed by rule before the draw.
- **Purpose:** to confirm cap residue names, atom names, coordinate dtype and quantization, the wrapping
  state and box layout; to develop the held-out loader; and to run the harness end to end as a dry run.
- **Status:** its results are development results, never held-out. Its superfamily is excluded from the
  pool.

### 3.2 Disclosure

- The AdK development results (`REPORT_ADK.md`) were known when this protocol was written, including when ρ
  (§7) was chosen.
- The AdK point-estimate illustration in §7.4 was computed after ρ was fixed.
- No mdCATH trajectory, DSSP, RMSD, radius-of-gyration or `mdcath_source.h5` content has been read. The
  only domain-specific trajectory metadata read is the `numFrames` table of `1r9lA02` in the authors'
  notebook.
- Before the freeze, aggregate metadata statistics were computed to size the budget: pool counts, number of
  superfamilies, and quantiles of file size and length. No individual domain was chosen or inspected.
- The paper names three example domains in its figures (`5j8eA00`, `2a06B02`, `2xryA01`). Only the caption
  text was read; the figures were not viewed. These domains stay eligible.

## 4. Eligibility and selection

### 4.1 Metadata eligibility, applied before the draw

1. Listed in `mdCATH_domains.txt` **and** present as `data/mdcath_dataset_<id>.h5` in the pinned tree, with a
   size and LFS SHA-256 in the tree listing.
2. Present in the CATH 4.2.0 domain list; superfamily = columns 2–5 (C.A.T.H).
3. CATH domain length (column 11) between 50 and 300 residues inclusive. The upper bound limits per-file
   size and runtime. It is a budget choice, not a scientific one.
4. Superfamily not in the exclusion set (§3).

### 4.2 Seeded two-stage draw (one domain per superfamily)

1. Sort eligible superfamilies by their numeric (C, A, T, H) tuple.
2. Draw a permutation with `numpy.random.default_rng(20261007).permutation(n_superfamilies)`.
3. Go through the superfamilies in permuted order. For each, choose one domain with `rng.integers(0, m)`
   from that same generator, over the superfamily's eligible domains sorted by name.
4. The first 30 are the **primary sample**. The next 10, in order, are **reserves**.

Superfamilies are drawn uniformly, so large superfamilies are not over-represented. Because there is one
domain per superfamily, every domain is its own cluster.

### 4.3 Download-time structural eligibility

This check is mechanical, uses no dynamics-derived quantity, and is applied before any scoring. A primary
domain that fails is replaced by the **next reserve in order**. Every replacement is logged with its reason.

- **File integrity:** SHA-256 equals the LFS object id in the pinned tree, and the size matches.
- **Layout:** root attribute `layout` = `mdcath-only-protein-v1.0`; one domain group; `numChains` = 1.
- **Required trajectories:** temperatures `320` and `348`, each with replicas `0` and `1`, each with
  `numFrames` ≥ 200 (≥ 200 ns at 1 ns per frame, so k = 10 is ≤ 5 % of any database after burn-in).
- **Shapes:** `coords` of shape (`numFrames`, `numProteinAtoms`, 3) and a diagonal `box`.
- **Backbone:** residues with exactly one N, one CA and one C. Caps lack CA and are therefore not window
  residues. These residues must form one run of consecutive `resid` values, and the run must have at least
  50 residues.
- **Coordinates:** finite in every frame of the four trajectories.
- **Abort:** if more than 10 reserves are exhausted, stop and report **Insufficient evidence**.

## 5. Trajectory handling (per domain, per temperature)

- **Query:** replica 0. **Database:** replica 1. 320 K and 348 K are **separate cells and never pooled**.
  379–450 K are excluded: the paper reports a shift toward coil at high temperature, which changes the
  question.
- **Shared start:**
  - Every replica at every temperature starts from the same equilibrated configuration, so replicas 0 and 1
    are **not independent initial conditions**. Their early frames are expected to be similar for that
    reason alone.
  - The first ⌊0.10 n⌋ frames of each replica (n = its `numFrames`) are excluded from query and database.
    This is the same 10 % fraction as the AdK gap rule, fixed in advance, not tuned.
  - A **sensitivity analysis without this exclusion** is reported descriptively and does not enter any
    decision.
- **Frame spacing:** assumed to be 1 ns, from the paper; there are no timestamps (§2).
- **Periodic boundaries:** box = diagonal of `box` × 10 Å. The Cα chain and the backbone (N, CA, C order) are
  made whole by sequential minimum image (WP-B `adk_lib.make_whole`). Changes are counted.
- **Precision:**
  - Every Cα method receives identical coordinates: values as stored, rounded once to float32 and held in
    float64. The number of values changed by the rounding is recorded.
  - The binding receives the same float32 values.
  - The target uses N, CA, C in float64 from the same whole-made frame, without rounding.
- **Windows:**
  - W ∈ {5, 8}, all start positions inside the backbone run.
  - Cells are temperature × W: **four cells** (320/W5, 320/W8, 348/W5, 348/W8).
- **Per-frame exclusions** (AdK §3 and §10 A1, unchanged): degenerate bond angles (sin < 1e-3) in the Cα
  torsions or in the target dihedrals, backbone bonds > 2.0 Å after reconstruction, or a binding error. Each
  excludes the (window, frame) pair for **all** methods and the target. Counts are reported.
- **Binding identity and reference check** (AdK §10 A2, unchanged):
  - Binding: WP-A wheel; extension SHA-256 and production source hashes must match the WP-A record.
  - Binding τ/κ vs NumPy on every frame, ≤ 1e-9 rad. A violation stops the run.
  - Horn vs Kabsch on sampled pairs, ≤ 1e-6 Å.

## 6. Metrics and aggregation

### 6.1 Per window (AdK §6 and §10 A3, unchanged)

- Tie-aware membership with tol = 1e-12·(1 + d_k): uniform selection within the boundary tie set.
- Overlap@k and raw regret@k (rad) for k ∈ {1, 10}.
- **Normalized regret@k** NR = regret / (mean_d T − oracle). Undefined when the denominator is ≤ 1e-12; such
  queries are excluded and counted.
- Regret in [−1e-10, 0) is set to 0 and counted. Values below −1e-10 are numerical failures; in the
  target-oracle computation itself they stop the run.
- Spearman ρ over 20,000 seeded (q, d) pairs. Seed = 20261005 + 1000·W + i + 100000·(domain index in the
  primary order) + 10000000·(0 for 320 K, 1 for 348 K).
- **C1 conditioning** (AdK A1): ill-conditioned (R < 0.1) queries have no C1 result; such database frames are
  +∞. Coverage is reported.

### 6.2 Primary quantity

NR@10. It combines the target-distance cost of a method's neighbours with the random reference (1) and the
oracle (0). Overlap@10 is secondary: with a continuous target and many near-equivalent frames it is small even
for good retrieval, as on AdK.

### 6.3 Aggregation (fixed order)

1. **Window:** mean NR@10 over the window's valid queries with defined NR. For C1, only queries with a C1
   result.
2. **Domain:** for domain d, cell c and method m, NR_m,c(d) = median over windows with a defined value. A
   domain is **undefined in a cell** if fewer than 10 windows are defined for any non-C1 method. It is then
   reported and left out of that cell.
3. **Cell:** medians across the defined domains of the cell.

### 6.4 C1 coverage-matched values

C1 is never compared on a different query set. For each window, every method's mean NR@10 is recomputed on
the queries where C1 has a result. Domain medians follow step 2. These **coverage-matched** values are used
for C1's classification only.

## 7. Decision criterion

### 7.1 Estimands (per cell c)

Best method:

- b(c) = argmin over {M1, M2, B1, B2, B3} of median_d NR_m,c(d).
- C1 is not eligible as best, because its coverage is conditional.
- An exact tie (|Δ| ≤ 1e-12) is broken by the order B2, B1, B3, M1, M2. The order is irrelevant to the
  decision, because tied methods have L = 0.

Gain over random, for method m:

- G_m,c = 1 − median_d NR_m,c(d).
- This is the fraction of the random-to-oracle regret gap that m closes (0 = random, 1 = oracle).

Relative loss vs the best method:

- L_m,c = median_d [NR_m,c(d) − NR_b,c(d)] / G_b,c.
- This is the fraction of the best method's gain that m gives up, using paired, per-domain differences.
- For C1, both numerator terms use coverage-matched values (§6.4). G_b is unchanged.

### 7.2 Uncertainty

Percentile bootstrap over domains:

- 2,000 resamples with `numpy.random.default_rng(20261008)`.
- Within a cell, domains are resampled with replacement from that cell's defined domains. The same resampled
  indices are used for every method, so the differences stay paired.
- **b is re-selected in every resample**, so the interval includes uncertainty about which method is best.
- The 95 % interval is the 2.5–97.5 % percentiles of L_m,c and G_m,c.
- Domains are the independent unit. Windows and frames within a domain are never resampled.

### 7.3 Practical margin ρ = 0.10 and per-cell classes

**Equivalence to the best method:**

| Class | Rule |
|-------|------|
| **E** (practically equivalent to the best) | Upper bound of L_m,c ≤ ρ |
| **I** (inferior) | Lower bound of L_m,c > ρ |
| **U** (inconclusive) | Otherwise |

**Informativeness:**

| Class | Rule |
|-------|------|
| **Informative** | Lower bound of G_m,c > ρ |
| **Uninformative** | Upper bound of G_m,c ≤ ρ |
| **Inconclusive** | Otherwise |

**Why ρ = 0.10.** It is a stated value judgment, not a calibrated quantity.

- Giving up at most a tenth of the best method's improvement over random retrieval is an acceptable cost for
  an engine component that may bring other properties: alignment-free comparison, invariance, a streaming
  update, or lower cost.
- The same fraction defines "informative": a method must close at least a tenth of the random-to-oracle gap.
- The margin is **relative**, so it is scale-free across domains, window sizes and temperatures. It
  automatically accounts for the random reference (G) and the target-distance cost (regret), which the
  replaced 0.05-overlap rule did not.

**Limitations of the criterion.**

- ρ is not derived from a downstream task. A different engine use (e.g., coarse clustering vs. fine
  retrieval) could justify a different margin.
- NR is an average target-distance cost. It does not show whether the errors concentrate in particular
  windows: terciles are reported for that.
- Percentile bootstrap intervals for medians of 24–30 domains are approximate. Re-selecting b makes them
  non-smooth.
- The 6 methods × 4 cells × 2 statistics are not adjusted for multiplicity. The intervals are decision
  thresholds, not significance tests.
- Replicas share a starting configuration (§5), so cross-replica retrieval may be easier than retrieval
  between independent runs.
- Overall, the result describes these domains and conditions only (§1).

### 7.4 AdK point-estimate illustration (computed after ρ was fixed; not a decision)

`results/wpd_adk_THIRD_PARTY_MD_DERIVATIVE_rule_illustration.json` applies §7.1 to the AdK windows, with
windows standing in for domains and no intervals.

- **B2:** best in every AdK cell.
- **M1/M2:** lose more than ρ of B2's gain at W = 5 but less at W = 8. The two W = 8 point estimates lie
  **just below** ρ.
- **C1:** closes less than ρ of the gap.

So the rule does not make the AdK result pass uniformly. Because the W = 8 values are near the boundary, ρ
was **not** adjusted after this was seen. Under §7.3, near-boundary point estimates are classed U unless the
held-out interval is narrow.

### 7.5 Overall outcome per method (precedence order)

| Order | Outcome | Rule |
|---|---|---|
| 0 | **Insufficient evidence** (whole study) | Fewer than 24 defined domains in any cell; reserves exhausted; or any stop-rule failure (hash, binding identity, reference check, target-oracle) |
| 1 | **Selected** | E and Informative in all four cells |
| 2 | **Not selected** | I in all four cells, or Uninformative in all four cells |
| 3 | **Cell-dependent** | At least one E-and-Informative cell and at least one I cell. The role is restricted to the E cells (e.g., W = 8 only) and stated as such |
| 4 | **Inconclusive** | Any other pattern, including any U that prevents rules 1–3 |

- C1 is classified with the same rules on coverage-matched values. Its coverage is reported beside the class.
- Several methods may be **Selected**. Choosing among them uses criteria outside this protocol (cost,
  invariance, interpretability), stated when that choice is made.

### 7.6 Reported but not deciding

- Overlap@1 and overlap@10 against the random reference k/n_D.
- Raw regret@10 in rad.
- Spearman ρ.
- NR@1.
- Terciles of target spread (median T per window) within each cell.
- Coverage and exclusion counts; C1 R distribution; cis-ω counts.
- Per-domain values.
- The no-burn-in sensitivity analysis.

## 8. Download, storage and runtime budget

- **Scope:** the exact files, sizes, URLs and LFS SHA-256 values are produced by
  `benchmarks/audit/wpd_local_retrieval/mdcath_select.py`. Its outputs are
  `results/wpd_mdcath_selection.json` and `HELDOUT_DOWNLOAD_REQUEST.md`, written after this protocol is
  frozen. The request covers:
  1. the schema sample (download first);
  2. the 30 primary files;
  3. reserves, **only** on a logged replacement, each from the listed set.
- **Whole files** are downloaded from the pinned revision so that each file can be checked against its LFS
  SHA-256. Range reads of only the needed replicas would transfer less, but cannot be verified that way.
- **Storage:** `benchmarks/audit/wpd_local_retrieval/data/mdcath/` (git-ignored). Raw files are kept until
  the held-out review closes. Derived Cα/backbone arrays go in `data/derived/`.
- **Runtime:** estimated in the selection output from the AdK throughput, with a nominal 450 frames per
  replica. This is an estimate, not an established benchmark.

## 9. Execution order (after authorization)

1. Download the schema sample and verify its SHA-256.
2. Develop and test the held-out loader on the schema sample; its results are development only.
3. Freeze the loader and harness by hash.
4. Download the primary files, verify them, and apply the §4.3 structural eligibility and replacements.
5. Score.
6. Report every method and cell, including Inconclusive and negative outcomes, without amending this protocol.

Any amendment after the freeze is recorded separately, with its reason, before any affected computation.
