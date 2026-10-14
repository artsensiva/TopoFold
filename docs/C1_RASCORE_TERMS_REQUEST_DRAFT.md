# DRAFT — not sent. For the project owner to review, edit and send if appropriate.

**To:** Rascore maintainers (repository `mitch-parker/rascore`; Dunbrack laboratory, Fox Chase Cancer Center)
**Subject:** Terms for using Rascore RAS conformation labels in a structural-retrieval evaluation

---

Dear Rascore maintainers,

We are planning a small, preregistered evaluation of whether backbone (Cα-only) local geometry of the RAS switch-I
loop recovers the expert conformational labels described in Parker et al., *Cancer Research* 82:2485 (2022),
doi:10.1158/0008-5472.CAN-22-0804. The planned population is primarily KRAS, with other RAS isoforms considered only
in secondary analyses.

The data snapshot we would use is the label table in your GitHub repository:
- repository `https://github.com/mitch-parker/rascore`, commit `8f77dc201cf9b1eb838d8e0d4d15bc0397dc3769`
  (2024-03-19);
- file `src/rascore/util/data/entry.tsv`, SHA-256
  `bac0fcfbe9992fc816eb95d33e64a110b9dfd9baa16eb6352cc13f663ff8132a`.

The repository is licensed under Apache-2.0, and the article under CC BY-NC-ND 4.0. We could not determine whether the
article licence extends to Supplementary Datasets S1–S4, and we could not find a statement specific to the label data.
We would be grateful for clarification on the following:

1. **Internal evaluation:** may we use these labels to evaluate our methods internally?
2. **Publication:** may we publish results derived from them? This would cover summary metrics and, if permitted,
   per-chain labels or train/test splits referencing PDB chain IDs, with attribution to your paper and database.
3. **Commercial context:** our project is a research prototype that may later be used commercially. Would use of the
   labels in that context be permitted, and under which terms?
4. **Governing terms:** if the repository copy and the article's supplementary datasets are intended to carry
   different terms, which applies to the repository file above?

We will cite the paper and the Rascore database in any output. Thank you for making this resource available.

Kind regards,
[Your name]
[Affiliation, or "Independent researcher"]
[Contact details]
