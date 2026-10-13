# WP-B G0–G1 Specification (frozen before trajectory inspection)

**Scope:** G0 (data acquisition and integrity) and G1 (independent state assessment) for the
candidate AdK development dataset. **No TopoFold or candidate-feature scoring, no retrieval
comparison, no localization, no threshold tuning.** AdK is a development dataset, not held-out
validation.
**Freeze:** this file and `wpb_spec.json` are hashed in `results/adk_spec_freeze.json` before
any trajectory frame beyond the headers is read. Any later change is recorded as an amendment
with a reason; definitions are not revised to obtain a desired outcome.

**Prior knowledge disclosed.** On 2026-10-04 the following were inspected:
- the two `.gro` files;
- the headers of the first three frames of each `.xtc` (200 ps spacing, precision 1000/nm,
  95,175 / 95,196 atoms).

No conformations, RMSDs or any trajectory-derived structural quantity have been computed.

## 1. Sources (source statements vs facts to verify)

| Item | Source | Licence | Expected (Zenodo API) |
|------|--------|---------|------------------------|
| `ORADD_open_state.gro` | https://zenodo.org/records/5583119/files/ORADD_open_state.gro, DOI 10.5281/zenodo.5583119 | CC BY 4.0 | md5 a73df3d6a60e69515bb7a7758af59f17; 6,567,118 B |
| `ORADD_closed_state.gro` | same record | CC BY 4.0 | md5 541635b532a85a5e98ff14a7228c6a2b; 6,568,567 B |
| `ORADD_open_state_traj.xtc` | same record | CC BY 4.0 | md5 05de6fa777dd5c86ba3e4c32b2c9080f; 780,301,896 B |
| `ORADD_closed_state_traj.xtc` | same record | CC BY 4.0 | md5 014063067a9617dd89a39e35f7350d25; 580,049,352 B |
| 4AKE (open reference) | https://files.rcsb.org/download/4AKE.pdb | PDB data, CC0 | SHA-256 recorded on download (RCSB files may be remediated over time) |
| 1AKE (closed reference) | https://files.rcsb.org/download/1AKE.pdb | PDB data, CC0 | as above |

Source statements, from Orädd et al., *Sci. Adv.* 2021 (PMC8597995): unrestrained 450 ns (from
4AKE) and 335 ns (from 1AKE); CHARMM36, GROMACS 2016.4, 303.15 K. These are compared with the
files and never assumed.

## 2. Storage and environment

- **Raw data:** `benchmarks/audit/adk_feasibility/data/raw/`, never modified. Derived arrays
  go in `data/derived/`. Both are covered by a narrow ignore rule; ≈ 1.4 GB is needed
  (377 GB free at freeze).
- **Isolated environment:** `benchmarks/audit/adk_feasibility/build/venv`, created with
  `python3 -m venv` from CPython 3.14.4. Packages: `numpy==2.5.3`, `MDAnalysis==2.10.0`
  (cp314 wheel). The full `pip freeze` is recorded. Existing environments are not modified.

## 3. Atom/residue mapping and units

- Protein Cα atoms only: selection `protein and name CA`. Expected: 214 atoms, resid 1–214,
  one chain.
- The sequence must equal the 4AKE and 1AKE canonical sequences. Any mismatch stops the run.
- Crystal chains: the Cα of residues 1–214. For alternate locations, the first altloc (blank or
  `A`) is used and every altloc encountered is recorded. Missing residues are recorded; RMSDs then
  use the common residue set, which is recorded.
- Units: Å and ps, as reported by MDAnalysis, which converts xtc nm to Å. xtc stores coordinates
  quantized to 0.01 Å; float64 arithmetic does not recover precision lost in the file.

## 4. References and alignment core (independent of TopoFold and of trajectories)

- **Primary references:** chain A of 4AKE (open) and chain A of 1AKE (closed). Precedent: the
  same chains define the transition endpoints in Seyler et al. 2015 (verified, WP-B §2B). Chain B
  is evaluated as a sensitivity analysis only.
- **Core:** derived from the two crystal Cα sets only. Start from all common residues. Then
  repeatedly (at most 50 times) Kabsch-superpose 1AKE:A on 4AKE:A over the current set and keep
  the residues whose deviation is < 1.0 Å, until the set is stable. If fewer than 30 core
  residues remain, the run stops.
  - Sensitivity cutoffs 0.5 Å and 2.0 Å are reported (core sizes and the effect on the
    observables), never substituted for the primary core.
- Literature domain names (LID, NMP, CORE) are not used, because their residue ranges were not
  verified from a primary source. Non-core residues are reported as contiguous crystal-derived
  segments.
- **d_ref:** the Cα RMSD over all common residues between 4AKE:A and 1AKE:A after core
  superposition. This is the reference separation scale.

## 5. Periodic boundaries and molecular integrity

- The box must be orthorhombic in every frame (angles within 1e-3° of 90°); otherwise stop.
- Cα chains are made whole per frame by sequential minimum image: each Cα is shifted by box
  multiples so that it is nearest its predecessor.
- **Recorded:** the number of frames where the raw maximum consecutive Cα distance exceeded
  4.2 Å, before the repair.
- **Required after the repair:** maximum consecutive Cα distance ≤ 4.2 Å in every frame;
  otherwise stop.
- Original coordinates are never overwritten.

## 6. G1 observables (continuous, primary)

For each frame of each trajectory, kept separate (open-start and closed-start never pooled):
1. `rmsd_open`, `rmsd_closed`: Cα RMSD over all common residues to 4AKE:A and to 1AKE:A after
   Kabsch superposition on the core.
2. `core_rmsd_open`, `core_rmsd_closed`: fit quality on the core itself.
3. `delta = rmsd_closed − rmsd_open`.
4. Per non-core segment (≥ 5 residues): segment Cα RMSD to each reference after core
   superposition.

Summaries: time series, quantiles, per-trajectory first/last 10 % means. Autocorrelation:
integrated autocorrelation time via Sokal windowing (c = 5), with the effective sample size
reported as an **estimate**.

## 7. Descriptive proximity classification (provisional, not a validated state definition)

- **Primary rule.** A frame is *near open* if `rmsd_open ≤ 0.25·d_ref` and
  `rmsd_open < rmsd_closed`. *Near closed* is defined symmetrically. Every other frame is
  *ambiguous/intermediate*.
- **Sustained visit:** ≥ 25 consecutive frames (5 ns at the 200 ps spacing). Shorter runs are
  *brief* and are reported separately (count, longest run).
- **Sensitivity:** fractions 0.20 and 0.33 are reported alongside; the primary is not changed.
- **Transition evidence:** a trajectory starting near one reference contains a sustained visit
  near the other. Preference for a reference by RMSD is **not proof of a conformational
  transition**. Classification is descriptive. If the classes are dominated by ambiguous frames,
  transition classification is reported as unresolved.
- **Not inferred:** equilibrium populations, rates or free energies.

## 8. Proceed / stop rules and outcome

**Stop G0 on:**
- any md5 mismatch;
- an unreadable file;
- an atom-count mismatch between `.gro` and `.xtc`;
- missing or extra Cα;
- a sequence mismatch;
- non-monotonic timestamps;
- a non-orthorhombic box;
- failure to make a chain whole;
- a core of < 30 residues.

**Recorded discrepancy, not a stop:** a frame count or duration differing from the paper's
statement. md5 establishes file integrity; the paper's durations are source statements. *This
supersedes the earlier draft rule (feasibility doc §6) that made a frame-count disagreement a
stop condition. Reason: it conflated integrity with consistency of a source statement.*

**Outcome (one of):**

| Outcome | Condition |
|---------|-----------|
| **Suitable for the proposed development task** | Sustained visits near **both** references occur, with at least one trajectory containing sustained visits near both |
| **Suitable only for a restricted task** | Each trajectory stays near its own start reference, or ambiguous, with no sustained visit near the other. Only pooled starting-state discrimination is then feasible, and it cannot speak to spontaneous transitions |
| **Unsuitable** | G0 fails irreparably, or neither trajectory stays near its start reference |
| **Unresolved pending specified evidence** | Ambiguous frames dominate, or the sensitivity fractions disagree on the categories above |

## 9. Outputs

| File | Content |
|------|---------|
| `results/adk_spec_freeze.json` | Hashes and time of this spec |
| `results/adk_g0_provenance.json` | Sources, checksums, environment, topology/trajectory facts, commands, exit codes |
| `results/adk_g1_summary.json` | Core, d_ref, observables summaries, classification, autocorrelation, outcome |
| `results/adk_g1_frames_{open,closed}.csv` | Per-frame observables |
| `REPORT.md` | Concise feasibility report |

These are third-party MD data: no SYNTHETIC label. Every file states the source DOI.
