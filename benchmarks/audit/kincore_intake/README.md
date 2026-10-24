# KinCoRe intake (2026-10-24, `/dataset-intake`)

**Record:** `intake_record.json`. **Decision: conditional-go.** No data file was downloaded; only HTML pages and HTTP
headers were read.

- **Licence:**
  - The KinCoRe website states no licence or terms of use.
  - The DunbrackLab standalone GitHub repositories have no licence file.
  - The 2019 PNAS paper (and presumably its Datasets S1–S3) is CC BY-NC-ND 4.0.
- **Sizes:** `PK_labels_PDB.fasta` is 5.5 MB. The all-structures ZIP is 1.6 GB and is not needed, since structures come
  from RCSB, CC0.
- **Evidence class:** agreement between TopoFold's side-chain descriptors and KinCoRe labels is a **geometry
  implementation check**, not independent biological validation. For Cα-only representations, the labels are
  independent expert-defined targets.
- **Open condition (owner):** either approve internal, non-commercial use without a stated licence, or obtain written
  permission from the Dunbrack laboratory.

Protocol draft: `../kin_gen_1/PROTOCOL_DRAFT.md`.
