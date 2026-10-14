# Decision Memo: First Oncology-Relevant Validation Task for Protein Local-Conformation Analysis and Retrieval

**Status:** decision memo for review, 2026-10-13.
- No task is frozen, no data has been downloaded or inspected for outcomes, and no TopoFold or baseline run has been
  made on any candidate.
- No threshold below was set by inspecting evaluation outcomes. Where a number is needed, the memo says how it is to
  be fixed before data access.

**Constraints applied:**
- `docs/ONCOLOGY_VALIDATION_REQUIREMENTS.md`: question first; independent ground truth; holdouts that match the claim;
  equal-information baselines distinct from stronger practical comparators; claims at the tier tested.
- `docs/ENGINE_MISSION_AND_ROADMAP.md` §6–§7: the first deliverable is protein local-conformation analysis and
  retrieval. Pocket detection, kinetics, free energies and allostery are outside it.
- AGENTS.md rules 5, 6 and 9.

## 1. Evidence reviewed

**In the project:**
- Roadmap §6 lists **ABL kinase conformational change** (track 1) and **KRAS conformational/pocket investigations**
  (track 4).
- Both record blockers:
  - Cα-only descriptors cannot see side-chain χ (C4/H5);
  - there are no inter-element distances and no cavity representation;
  - there is no real MD;
  - earlier ABL PDB IDs (2GQG, 1IEP) "must be re-verified before reuse".
- The only validated evidence is geometric: WP-D held-out retrieval within mdCATH. A and B2 each met the criterion on
  a different geometric target.
- CMP-1 shows that established RMSD libraries suffice computationally for retrieval at development scale.

**External label resources** (verified by literature search on 2026-10-13; nothing downloaded):

| Resource | What it labels | Source | Terms |
|---|---|---|---|
| **Kincore** | Each kinase chain in the PDB: activation-loop conformation (spatial DFG label + dihedral cluster label) and inhibitor type | Modi & Dunbrack, *PNAS* 116:6818 (2019), doi:10.1073/pnas.1814279116; Modi & Dunbrack, *NAR* 50:D654 (2022), doi:10.1093/nar/gkab920 | Download available; **licence not found** |
| **Rascore** | Human KRAS/NRAS/HRAS PDB structures: switch-I (SW1) and switch-II (SW2) conformations, from Y32/Y71 positions and backbone dihedrals, linked to nucleotide state and bound partners or inhibitor sites | Parker, Meyer, Golemis & Dunbrack, *Cancer Res* 82:2485 (2022), doi:10.1158/0008-5472.CAN-22-0804; 721 structures at publication | **Licence not found** |
| **PDB coordinates** | — | wwPDB usage policy | CC0 1.0 |

**Critical property of both label sets:** they are **derived from the same experimental coordinates**, by expert
clustering of specific atom positions and dihedrals.
- They are independent of TopoFold (rule 5), but they are **not** independent of geometry.
- Agreement with them supports a claim about *agreement with an expert geometric classification*, not about function.
- Kincore's spatial DFG label depends on the Phe side-chain position, which Cα-only descriptors cannot observe (C4).

## 2. What each kind of claim would require

| Claim | Example | Independent evidence that would justify it | Within reach of the first deliverable? |
|---|---|---|---|
| **Geometric similarity / retrieval** | "The method retrieves structures whose switch region has the same expert-defined conformation" | Expert structural labels (Rascore/Kincore), held-out structures, equal-information geometric baselines | **Yes**, if labels, terms and counts are verified |
| **Functional state** | "Local geometry indicates GTP- vs GDP-bound RAS" or "active vs inactive kinase" | Labels from non-geometric evidence: bound nucleotide identity from chemical components in the entry; activity assays. Held-out by protein or time. Must beat baselines that read the same coordinates | Partly. Nucleotide identity is a candidate secondary label (§3, C1), but its relation to geometry is already published, so it shows predictability, not new biology |
| **Drug response** | "Retrieval predicts which inhibitor class binds" or "predicts IC50 / cellular response" | Measured binding or activity data for the same proteins/variants, an assay-level split and practical comparators (docking, ligand-based models) | **No.** No such dataset has been identified or verified |
| **Pocket detection** | "Finds the SII-P / DFG-out pocket" | Cavity representation on all heavy atoms; pocket ground truth; pocket-detection baselines (e.g. fpocket, PocketMiner — to be verified) | **No.** No cavity representation exists (roadmap §6) |
| **Clinical benefit** | Any patient-outcome statement | Appropriate clinical studies | **Outside the project scope** (`ONCOLOGY_VALIDATION_REQUIREMENTS.md` §2) |

## 3. Candidate tasks (three compared)

### C1 — RAS switch-region local-conformation retrieval in experimental structures

**Intended user decision.** Given a new RAS G-domain structure or model (e.g. a design model or an MD snapshot of
KRAS), which expert-defined SW1/SW2 conformational class does its switch region most resemble? The use is triage:
selecting reference structures with a matching switch conformation for comparison or modelling.

**Ground truth:**
- Rascore SW1/SW2 conformation labels, which are geometric and expert-derived.
- Secondary, non-geometric: bound nucleotide state taken from ligand identity in the PDB entry.
- Neither comes from TopoFold.

**Data:** RCSB PDB coordinates (CC0); Rascore label tables (**terms unverified**).
- KRAS, NRAS and HRAS share G-domain residue numbering, so residue correspondence is simple. This must be verified per
  structure (insertions, tags, missing residues).
- Unresolved switch residues are common in RAS structures. The number of usable chains after exclusions is
  **unknown** and must be counted from metadata before any method is run.

**Leakage risks:**
- multiple chains or entries of the same complex;
- inhibitor series from one study;
- near-identical structures across releases;
- labels whose definition shares atoms with the descriptors (Y32/Y71 Cα positions lie inside the switch windows).

**Mitigation:**
- cluster by entry, ligand chemotype and sequence identity before splitting;
- report a baseline built directly from the label-defining coordinates as an upper reference.

**Equal-information baselines** (same chains and residues, Cα or backbone only, same reference set):
- switch-window Cα RMSD top-k, using an established library (MDTraj/MDAnalysis per CMP-1);
- backbone φ/ψ (sin, cos) distance, if backbone atoms are used;
- Cα–Cα distance vectors (B1).

**Stronger-information comparator, reported separately:** Rascore's own assignment rules, the label generator itself,
which uses side-chain positions.

**Held-out design:**
- **Temporal** split: entries released after a cutoff date fixed before data access are held out.
- **Across-isoform** variant: develop on KRAS, hold out NRAS/HRAS.
- The claim must match the split: within-family versus across-isoform.

**Success and failure (structure fixed now; numbers fixed before data access by a precision calculation):**
- **Primary metric:** top-1 label-transfer balanced accuracy over SW1 and SW2 classes, resampled by entry cluster.
- **Success:** the TopoFold descriptor (κ–τ window distance "A") is non-inferior to the strongest equal-information
  baseline within a margin set in advance.
- **Failure:** worse than that baseline beyond the margin.
- **Insufficient evidence:** too few held-out clusters per class.

**Minimum infrastructure:**
- read the declared residue range (e.g. switch-I and switch-II windows) of each chain from PDB/mmCIF with an
  established parser, handling altlocs and missing residues explicitly;
- join the label table;
- a provenance manifest.

No general loader and no new kernel are needed.

### C2 — Kinase activation-loop conformation retrieval (Kincore), with ABL1 as the held-out protein

**Intended user decision.** Which Kincore activation-loop class (e.g. DFG-in/out with a dihedral cluster) does a
kinase structure most resemble? The use is structure selection for kinase-inhibitor modelling (type I vs type II).

**Ground truth:** Kincore labels, geometric and expert-derived. The spatial DFG label depends on the Phe side-chain
position.

**Data:** PDB (CC0); Kincore tables (**terms unverified**). The number of ABL1 chains is **unknown**.

**Leakage risks:**
- homologous kinases;
- the same inhibitor across kinases;
- labels defined partly by side-chain atoms that Cα descriptors cannot see.

Cross-kinase residue correspondence needs a kinase alignment (from Kincore or a similar resource; to be verified).

**Baselines:** activation-loop Cα RMSD; backbone φ/ψ; Kincore's standalone assigner as a stronger practical
comparator.

**Held-out design:** develop on other kinases and hold out ABL1 (across-protein). Clusters by kinase family.

**Success and failure:** as in C1, per class. The expected weak point must be stated in advance: Cα-only descriptors
**cannot** resolve labels defined by Phe side-chain position (C4).

**Minimum infrastructure:** as C1, plus a residue-correspondence table across kinases. That is more work and more
leakage surface.

### C3 — Variant effects on RAS switch-state populations from MD

**Intended user decision.** Does an oncogenic variant shift the population of a switch conformation?

**Ground truth would need:**
- biophysical measurements of state populations per variant, independent of the simulations;
- MD with independent replicas per variant.

**Status:** neither has been identified or verified in this project, and the project has no RAS MD. Requirements §3
cannot be met now.

**Excluded as the first task:** essential evidence is missing. It is a candidate for a later tier only after a
verified experimental dataset is found.

## 4. Comparison

| | C1 RAS switch retrieval | C2 Kinase A-loop / ABL1 | C3 RAS variant MD |
|---|---|---|---|
| Oncology relevance | KRAS: oncogenic driver; inhibitor sites linked to switch conformations (Parker et al. 2022) | ABL1: imatinib-class inhibitors bind DFG-out kinases (relation to be cited from verified sources before use) | Direct variant question |
| Ground truth | Expert geometric labels + nucleotide identity | Expert geometric labels, partly side-chain-defined | Not identified |
| Fit to Cα/backbone descriptors | Good (backbone switch geometry), partial overlap with label atoms | Weak for the spatial DFG label (side chain) | — |
| Residue correspondence | Shared numbering (verify) | Needs cross-kinase alignment | Same protein |
| Data terms | PDB CC0; Rascore **unverified** | PDB CC0; Kincore **unverified** | — |
| Usable count | **Unknown** | **Unknown** | — |
| Infrastructure | Smallest | Larger | MD generation |

## 5. Recommendation

**Recommend C1** (RAS switch-region local-conformation retrieval in experimental structures) as the first validation
task. It is **conditional** on three verifications that are currently **missing**. Each is metadata-only and runs no
method.

1. **Rascore label access and terms:** the licence or usage terms for redistributing derived results, and an
   immutable versioned snapshot or retrieval date.
2. **Usable counts:** chains with resolved switch residues per SW1/SW2 class, after the exclusion and clustering rules
   are fixed in a protocol, to check that held-out clusters per class can meet a pre-set precision target.
3. **Label-definition overlap:** documenting which atoms define each label, so that the claim is stated as agreement
   with an expert geometric classification, with the label-defining baseline reported.

**If any verification fails, the conclusion is that essential evidence is missing.** C2 is not to be substituted
automatically: it has the side-chain-label blind spot and needs more infrastructure.

**What a positive C1 result would support:** "Backbone local-geometry retrieval recovers expert-defined RAS switch
conformations in held-out experimental structures at least as well as equal-information baselines."

**What it would not support:** functional activity, drug response, pocket detection or clinical benefit (§2).

**Unresolved choices for review:**
- **Split:** temporal or across-isoform as the primary held-out design.
- **Information level:** Cα-only or full backbone.
- **Metric:** the primary metric and its non-inferiority margin, with the precision target used to set it.
- **Secondary label:** whether nucleotide identity is included, and how it is extracted from PDB chemical components.
- **Licence fallback:** what to do if Rascore terms do not permit redistribution of derived labels, i.e. keep labels
  local and publish only IDs and hashes.

## 6. Engineering consequence

- **Reuse established libraries.** Computation uses established libraries where their checked accuracy is adequate:
  MDTraj/MDAnalysis RMSD for baselines, an established structure parser and NumPy.
  - CMP-1's 1e-3 Å acceptance is far below crystallographic coordinate uncertainty (typically tenths of an Å).
  - Accuracy at that level is therefore not the limiting factor here; this is an expectation to be stated, not a
    measured result.
- **V2b stays an audit artifact.** Performance is not a constraint at this scale (hundreds of chains, short switch
  windows). No kernel optimization or general loader is proposed.
- **Smallest adapter (proposed, not implemented):** a function mapping (PDB ID, chain, declared residue range) to an
  owned Cα (or backbone) array plus a residue mask. It would use an established parser (to be chosen and verified),
  record its provenance, and have explicit missing-residue/altloc policies. A label-join step reads the verified
  Rascore snapshot.
- **Order:** nothing is implemented before the §5 verifications and a frozen protocol.
