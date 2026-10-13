# Evaluation-Design Decision Memo: Next Real-MD Experiment

**Date:** 2026-10-04
**Status:** decision memo for review. No analysis was run for it; no TopoFold scoring.
**Inputs:**
- WP-A and WP-C reports;
- the WP-B G0–G1 report (`benchmarks/audit/adk_feasibility/REPORT.md`, outcome *unresolved
  pending specified evidence*, preserved unchanged);
- primary sources verified for this memo (listed in §6).

## 1. Three different phenomena, and which one the ordered representation addresses

| Phenomenon | What changes | What a window-local Cα torsion representation sees |
|------------|--------------|-----------------------------------------------------|
| **Local backbone-conformation change** | Backbone dihedrals of a few consecutive residues (loop flips, helix fraying or unwinding, hinge-residue torsion changes) | Yes, through Cα geometry. Each Cα virtual torsion τ_j (Cα j…j+3) is a function of the intervening backbone (ψ_{j+1}, ω_{j+1}, φ_{j+2} and bond geometry). The mapping is many-to-one, so changes that leave the Cα trace nearly unchanged are not seen |
| **Rigid-domain displacement** | Position and orientation of an internally rigid domain; its internal geometry is preserved | Only at the hinge residues, if their local torsions change. **Not** the amplitude or direction of the displacement: windows inside the moving domain do not change |
| **Interdomain rearrangement** | Relative distances and orientations between domains | Not observable from window-local features by construction (locality, tested to 1e-12 in WP-A). Needs interdomain distances or orientations |

**Conclusion.** The ordered (cos τ, sin τ) representation retained in WP-C addresses **local
backbone-conformation change, and its retrieval**. That matches the engine's first deliverable
(`ENGINE_MISSION_AND_ROADMAP.md` §7). It is not a measure of domain motion. The WP-B G1
observables (global and segment Cα RMSD after core superposition) measure rigid-domain and
interdomain change, which this representation is not designed to capture.

## 2. Can AdK support a development experiment for the local task?

**What G1 does and does not say.**
- G1 was unresolved because frames were rarely within 2.04 Å Cα RMSD of either crystal reference.
  That concerns **global reference proximity**.
- The large segment RMSDs (e.g. residues 114–170 ending 21–23 Å from both references after core
  superposition) are consistent with displacement and possibly internal change of that segment.
  They are **not** evidence of local torsion change, and must not be used as local labels.

**A full-backbone reference target exists in the data.** The trajectories are all-atom (G0: 95,175 /
95,196 atoms), so backbone φ/ψ can be computed from N, Cα and C. No compared method receives
those atoms; all of them get Cα only. A window-level φ/ψ distance therefore provides a
**continuous structural reference target** computed from information no method receives. It
is not independent in a physical sense: it is geometrically coupled to the Cα trace. It is not
experimental ground truth either.

**Caveat.** Because τ depends mainly on ψ and φ, this tests whether a Cα-only representation
tracks full-backbone local change on real dynamics: fidelity and retrieval on real data. It does
not test discovery of biologically meaningful sites.

### Blockers and limitations for an AdK *descriptive development* test (not held-out validation)

| Item | Kind | Status |
|------|------|--------|
| Backbone atoms N, Cα, C present and complete for all 214 residues in every frame | **Essential blocker** | **Resolved:** separate WP-D eligibility artifact (`benchmarks/audit/wpd_local_retrieval/results/wpd_adk_backbone_eligibility.json`). Both trajectories are eligible |
| Global open/closed state labels | Blocker **only for global-state tasks** | Unresolved (G1). The local task uses a continuous target instead, so it is not blocked |
| Domain definitions | Not needed for the local task | One primary-source definition verified: Formoso et al. 2015, LID 118–160, NMP 30–67, CORE 1–29, 68–117, 161–214, P-loop 7–13. Boundaries are convention-dependent and differ between papers. G1's crystal-derived non-core segments (30–79, 114–170) overlap these ranges, a descriptive consistency only |
| One trajectory per starting state | Limitation | Acceptable for a descriptive development test. Replicates are required before any inferential or population claim, not before descriptive retrieval |
| Simulation provenance | Limitation | Verified (Orädd et al., Methods): CHARMM-GUI systems from 4AKE/1AKE; CHARMM36; GROMACS 2016.4; 303.15 K; 25 ps restrained equilibration; semi-isotropic barostat, rationale not given. The MD served as a conformational pool for refining scattering models; frames every 1 ns were used |
| Ligand handling | Limitation | Not stated in the paper. Verified from the files (G0): both systems are apo, so 1AKE was simulated without Ap5A. Chain choice, water model and ion conditions are also not stated; the files show TIP3 and K⁺ |
| Drift and non-stationarity | Limitation | Observed in G1. τ_int/ESS were exploratory only. Frames are not independent, so a time-separated design is required |

**Verdict.** AdK is usable for an **exploratory, descriptive** development test of the local
task once the backbone-atom check passes. It is **not** a held-out evaluation and can support no
generalization claim.

## 3. Shortlist of independent real-MD datasets (not downloaded)

| Dataset | Content (verified from primary paper) | Access / licence | Independent local target | Cost | Fit to first deliverable |
|---------|----------------------------------------|------------------|---------------------------|------|--------------------------|
| **mdCATH** (Mirarchi, Giorgino, De Fabritiis, *Sci. Data* 2024, doi:10.1038/s41597-024-04140-z) | 5,398 CATH domains; all-atom systems; 5 temperatures (320, 348, 379, 413, 450 K) × 5 replicas, **all started from the same equilibrated configuration**. Per-domain HDF5 with coordinates (Å) over the N protein atoms, forces and per-frame DSSP; one frame per ns | Trajectories CC BY 4.0; HuggingFace download per domain | Backbone φ/ψ, provided the stored protein atoms include N, Cα, C (expected for protein atoms; to verify on a file). Per-frame DSSP computed by the dataset authors | Full set > 3 TB, but files are per-domain, so a few domains can be selected. Per-file size to be measured | **Good.** Many domains across superfamilies. Separate replicas give **temporal separation, but not independent starting conditions**. The low temperatures (320, 348 K) suit local change; the paper reports a shift toward coil at 413 and 450 K |
| **ATLAS** (Vander Meersche et al., *NAR* 2024, doi:10.1093/nar/gkad1084) | 3 × 100 ns replicates at 300 K; CHARMM36m; coordinates every 10 ps; protein-only download formats; DSSP and protein-block (Neq) analyses | Article CC BY 4.0. **Data licence not verified**, and must be before any use | Backbone φ/ψ; recomputed DSSP / protein blocks | ~10,000 frames per replicate in complete format; per-protein download | **Good.** A different force field and protocol from mdCATH make it a useful second held-out source |
| DESRES long trajectories (BPTI, ubiquitin, fast folders) | Long single trajectories | By request, non-commercial use according to secondary sources (unverified) | φ/ψ | Request | **Excluded for now** (access terms; commercial goals) |

Neither shortlisted dataset has been inspected. Domain or protein selection for held-out use
must follow the held-out procedure (feasibility doc §6): eligibility file written and hashed, then
a seeded draw, **before** any scoring.

## 4. Recommended next experiment (WP-D, to be specified as its own protocol)

**Question.** Do Cα-only window representations retrieve frames whose **full-backbone local
conformation** is closest, and how does the ordered representation compare with equal-information
Cα baselines?

| Element | Definition (to be frozen in the WP-D protocol) |
|---------|----------------------------------------------|
| Units | Residue windows (W ∈ {5, 8}, pre-registered) at all positions; frames as queries and database |
| **Continuous independent target** | Window φ/ψ distance: circular RMS over the φ/ψ of the residues spanned by the window's torsions, computed from N, Cα, C atoms that no method receives. Secondary target: DSSP agreement (H-bond based), where available |
| Methods (all Cα-only, same windows and frames) | Ordered D_ch and D_circ (WP-C); local Cα–Cα distance vector; window-superposed Cα RMSD; κ block; the historical scalar circular-mean distance as comparator |
| Metric | Overlap@k between each representation's k nearest frames and the k nearest by the φ/ψ target (k ∈ {1, 10}). Rank correlation (Spearman) between representation distance and target distance on pre-sampled frame pairs |
| **Temporal separation** | **mdCATH:** queries from one replica, database from another. These are different runs, not independent starting conditions. **AdK (single trajectories):** exact query, gap and database ranges in the WP-D protocol §5. The gap reduces temporal proximity but does not establish independence. Exploratory only |
| Statistical unit | Domain or protein, not frame. Uncertainty from a bootstrap over domains on the held-out set; no inference on AdK |
| **Success (retain the representation for the first deliverable)** | On held-out domains, median overlap@10 of D_ch is (a) within a pre-declared **descriptive** margin of 0.05 of, or above, the best equal-information baseline, and (b) above the historical scalar comparator. This is not a statistical non-inferiority test. Exact rules: `WPD_LOCAL_RETRIEVAL_PROTOCOL.md` §7 |
| **Failure** | D_ch is materially below the best Cα baseline on held-out domains. The representation is then restricted to its demonstrated advantages (chirality and branch continuity) or dropped in favour of distance features |
| Explicitly not claimed | Detection of states, sites, pockets or kinetics; superiority beyond the measured comparison |

**Two tiers, kept separate:**
1. **Exploratory development (AdK, data in hand).** Descriptive results only, used to debug the
   pipeline and the target definition. No generalization.
2. **Held-out evaluation (mdCATH).** Domains drawn by seeded random selection from CATH
   superfamilies disjoint from any development domain. Eligibility is fixed before download.
   Replica split, temperatures 320 and 348 K. ATLAS is a second held-out source once its data
   licence is verified.

## 5. Recommended route and remaining prerequisites

**Route:** *Cα-only local backbone-conformation retrieval against an independent full-backbone
φ/ψ target.* It runs as an exploratory AdK development tier, followed by a held-out mdCATH tier.
This route does not depend on the unresolved AdK global-state labels (G1/G2), which remain
blocked for global-state tasks.

**Prerequisites, in order:**
1. **WP-D protocol, frozen before any computation:**
   - target definition, including the exact torsion-to-φ/ψ residue mapping;
   - window sizes, metrics, margins and baselines;
   - temporal and replica splits;
   - statistical unit;
   - held-out eligibility and draw procedure.
2. **AdK backbone-atom check:** read-only, extends G0 to N/Cα/C presence and completeness. This
   is the essential blocker for the exploratory tier.
3. **mdCATH facts, from file metadata or a single sample file, under download authorization:**
   - which atoms the HDF5 `coords` contain (backbone atoms needed);
   - frame interval and length per temperature;
   - per-domain file size.
4. **Held-out selection:** eligibility file hashed, seeded draw, no scoring beforehand.
5. **Authorization** for the mdCATH subset download, and later for ATLAS once its data licence is
   verified.
6. **Implementation in harness code only,** extending the WP-C harness. No production API change.

## 6. Sources verified for this memo

- Orädd et al., *Sci. Adv.* 2021, Methods, via Europe PMC full text PMC8597995: simulation set-up
  as quoted in §2; ligand handling, chain and water model not stated.
- Formoso, Limongelli, Parrinello, *Sci. Rep.* 2015, doi:10.1038/srep08425 (PMC4325324): domain
  residue ranges in a figure legend.
- Mirarchi, Giorgino, De Fabritiis, *Sci. Data* 2024, doi:10.1038/s41597-024-04140-z
  (PMC11604666): dataset scope, licence, HDF5 fields, size.
- Vander Meersche et al., *Nucleic Acids Res.* 2024, doi:10.1093/nar/gkad1084 (PMC10767941):
  protocol, replicates, frame interval, formats. Article licence only.
- WP-B G0 provenance for file-level facts (apo systems, TIP3, K⁺, atom counts).

## 7. Corrections (2026-10-04, after review)

- τ vs φ/ψ: "depends mainly on ψ_{j+1}, φ_{j+2}" was replaced by a statement of geometric
  dependence on the intervening backbone, many-to-one.
- The φ/ψ target is described as a full-backbone **structural reference**, geometrically coupled
  to Cα and not experimental ground truth. The earlier text called it "independent".
- mdCATH replicas share their starting configuration. "Independent replicas" and "clean temporal
  separation" were overstatements.
- The 0.05 margin is a descriptive margin, not statistical non-inferiority.
- The AdK backbone-atom blocker is resolved by the WP-D eligibility artifact.
- The detailed design now lives in `WPD_LOCAL_RETRIEVAL_PROTOCOL.md`.
