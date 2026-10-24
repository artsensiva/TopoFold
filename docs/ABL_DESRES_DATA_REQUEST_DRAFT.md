# Request text (draft). Status 2026-10-24: the owner has sent this request. Not to be resent. Kept as the record of what was asked.

*Revised 2026-10-24: the requested sets now reflect the Methods text, which describes no apo M472I simulations (`docs/ABL_WT_M472I_PREREGISTRATION_DRAFT.md` §0).*

**To:** trajectories@deshawresearch.com (address given in the Data availability statement of Ayaz et al. 2023)
**Subject:** Request for Abl kinase trajectories from Ayaz et al., Nat Commun 2023 (WT and M472I)

---

Dear D. E. Shaw Research team,

We would like to request access to some of the Abl kinase trajectories from Ayaz et al., "Structural mechanism of a
drug-binding process involving a large conformational change of the protein target", *Nature Communications* 14:1885
(2023), doi:10.1038/s41467-023-36956-5.

**About us.** TopoFold is a research prototype for locating and comparing conformational changes in molecular
ensembles. It is not a commercial product, but it has future commercial ambitions, and we want our use of your data
to be clearly within your terms.

**Intended use.** A preregistered method evaluation: whether backbone-geometry comparisons of wild-type and M472I
ensembles localize differences that coincide with the HDX-MS differences reported in the paper. We would compare
against established methods and report negative results.

**Data requested.** Wild-type and M472I runs started from the same structure and in the same ligand state.
- **First choice:** apo wild type and apo M472I, because the HDX-MS comparison was made on the apo proteins. The Methods
  describe apo WT runs (from 1OPJ) but, as far as we can see, no apo M472I runs. Could you confirm whether apo M472I
  simulations exist?
- **Otherwise:** the imatinib-bound WT and M472I runs from 1OPJ (described as 5 × 20 µs plus 5 × 10 µs for WT, and
  10 × 10 µs for M472I). Optionally also the WT and M472I runs with imatinib and GNF-2 (3K5V), and the apo WT runs.

Could you tell us:
1. the actual run inventory for these two conditions (number of runs, length of each, and starting structure), since
   the Methods and Results sections give different counts for some sets;
2. the simulated construct and its residue numbering relative to human c-Abl;
3. simulation conditions (force field, temperature, protonation states, ions) for these runs;
4. the saved-frame interval and approximate file sizes;
5. whether trajectories containing protein heavy atoms only (no solvent) could be provided.

**Permissions.** We would be grateful for a separate answer to each:
1. **Internal method evaluation:** may we use the trajectories to evaluate our methods internally?
2. **Publication:** may we publish results derived from them, such as summary metrics and figures, with citation of
   the paper?
3. **Redistribution:** may we share derived files (e.g. extracted Cα coordinates) with reviewers or in a public
   repository, or must they stay private?
4. **Later commercial use:** if the project later becomes commercial, would continued use require a separate
   agreement, and with whom?

If you are not the right contact for the peptide-level HDX-MS data from the same paper, a pointer to the right person
would also be appreciated.

We will cite the paper and acknowledge D. E. Shaw Research in any output. Thank you for considering this request.

Kind regards,
[Your name]
[Affiliation, or "Independent researcher"]
[Contact details]
