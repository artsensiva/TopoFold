# Independent AI trial of TopoFold structpair 0.3.0: 1OPJ (A, B) vs 1OPK (A)

**This was an independent AI trial, not a human usability study.** The tester was Claude Opus 5.5
(`claude-opus-5-5`) running in the Claude desktop app (Code tab, built-in browser pane). The tester had no
earlier involvement with TopoFold and used no project memory, earlier trial reports or other agents' conclusions.
The shipped Abl demonstration and its conclusions were **not opened**. The only incidental exposure was one README
status line, "Kinase descriptors (KIN-DESC-1, four public Abl crystal chains…)", seen in a grep for install
instructions. No external assistance was received.

Tested commit: `19b5268180d533efedcbdef846e7b3f24250a767` (detached HEAD, contained in `origin/scientific-rebuild`).
Clone status after the trial: `git status --short` is empty. No files in the clone were changed. The guide's
`fetch_inputs.py` wrote the KRAS PDBs into the clone's git-ignored `benchmarks/audit/*/data_pdb/` directories.

## Plain-language summary

The tool installed and ran exactly as documented. The KRAS example reproduced the documented numbers. My
pre-registered Abl comparison (1OPJ chains A and B vs 1OPK chain A, residues 248–529) ran on the first attempt, and
all 9 output hashes matched the manifest.

The comparison shows one dominant local difference. Residues about 400–416 move by up to ~33 Å between the two entries
after a core fit; in Abl numbering this is the DFG motif and activation loop. The two copies of 1OPJ agree with each
other there to ≤1.4 Å. Smaller differences occur at the N-terminal end (248–252), the glycine-rich loop (264–273) and
residues 292–297 (β3–αC).

Several factors differ together between the entries: construct, ATP-site inhibitor, the D382N mutation, crystal form
and crystallization conditions. The comparison therefore cannot attribute any of these differences to a single cause.

The viewer worked, and it rejected a tampered copy. I found one reproducible viewer defect: after the results folder
is changed and changed back, the sidebar focus no longer matches the details panel and plots. The 3D view could not be
tested visually because the browser pane has no WebGL.

## Outcome by stage

| Stage | Outcome |
|---|---|
| Installation | **Completed following documentation** (guide §2 steps 1–5; no failures) |
| Command execution (KRAS smoke test + Abl) | **Completed following documentation**. Manifest re-checking used a self-written script, since no verifier is documented |
| Viewer interaction | **Completed using a self-devised workaround** (port 8501 was occupied, so I used 8517 and `--server.headless true`). The **3D highlight was NOT TESTED visually**: WebGL is unavailable in the browser pane, and only the plot's DOM data was inspected |
| Interpretation | **Completed** (own work, no assistance) |

## Time

- Wall-clock time: 2026-10-25 17:35:42 → about 17:55 +01:00, roughly 20 minutes.
- Active time: about the same, with no interruptions.
- AI tool latency is not comparable to human effort.

## 1. Environment and installation

- OS: Ubuntu 26.04.1 LTS, Linux 7.0.0-38-generic, x86-64.
- Toolchain: Python 3.14.4 (`/usr/bin/python3`), rustc 1.98.1, cargo 1.98.1.
- Starting state: a fresh `python3 -m venv` containing only pip 25.1.1. No relevant packages were preinstalled in the
  venv; MDAnalysis, Streamlit and maturin were not checked system-wide.
- Commands run, from the guide, in the `trial/clone` repository root:
  ```
  git clone https://github.com/artsensiva/TopoFold.git clone && git -C clone checkout 19b5268180d533efedcbdef846e7b3f24250a767
  python3 -m venv ../topofold-env && source ../topofold-env/bin/activate
  pip install maturin                      # maturin 1.15.0
  maturin develop --release                # 39 s wall; topofold 0.8.3 editable
  pip install MDAnalysis==2.10.0 matplotlib==3.11.2
  python benchmarks/audit/structpair_demo/fetch_inputs.py      # 4 × "ok" (SHA-256 checked)
  PYTHONPATH=python python -m topofold_inputs.structpair benchmarks/audit/structpair_demo/specs/kras2_5US4A_vs_5USJA.json --out ../outputs/kras2
  pip install streamlit==1.64.0 pandas==3.0.6 altair==6.3.0 plotly==7.1.0   # Pillow 13.0.0 → 12.3.0, as the guide notes
  ```
- Departures from the guide:
  - `--out` was `trial/outputs/kras2` instead of `~/structpair_out/kras2`, to keep artefacts in the trial directory.
  - The clone was checked out at the exact commit rather than the branch tip.
- KRAS stdout: `158 positions; global Cα RMSD 0.974 Å over 158 Cα; 10 exclusion rows; within-entry variation:
  unavailable; wrote ../outputs/kras2/report.md`. This matches the guide except for the added "within-entry
  variation" clause (documentation gap).
- KRAS manifest: 6 of 6 output hashes match (`logs/03_kras_manifest_check.log`). Comparison with the shipped reference:
  - `windows.csv`, `residues.csv` and `exclusions.csv` are byte-identical.
  - `summary.json` differs: it gains a `within_entry_variation` block, and the numbers are unchanged. The guide says
    the trial outputs were byte-identical, which is true only of the older reference (documentation gap).

## 2. Abl analysis (pre-registered)

- Pre-registration: `notes/PREREGISTRATION.md`, SHA-256 `d593bc64…`. It was hashed before the run, together with the
  spec (`notes/prereg_hashes_before_run.txt`).
- Inputs, downloaded 2026-10-25T17:39 from RCSB:

  | File | URL | SHA-256 |
  |---|---|---|
  | 1OPJ.pdb | https://files.rcsb.org/download/1OPJ.pdb | `6d3444a1f081b1e3b5079babd515d4e848586ec7055647fad9d294b9ee8fb7fb` |
  | 1OPK.pdb | https://files.rcsb.org/download/1OPK.pdb | `867ed2b25b80e0e2ba52f731e48df25fb8e982b052b89440a123d567c4384974` |

  mmCIF copies are also in `inputs/SHA256SUMS`.
- **Question:** where, and by how much, does the local Cα geometry of the Abl kinase domain (248–529) differ between
  1OPJ A and 1OPK A, compared with the copy-to-copy range between 1OPJ A and B?
- **Correspondence:** `same_numbering "248-529"`, 282 positions. It was checked in three ways:
  - Identical DBREF offsets: author = UniProt P00520 isoform I + 19, the "1b" numbering.
  - A per-residue name comparison over the common Cα positions.
  - Motif positions in author numbering: catalytic HRD/N at 380–382, DFG at 400–402, APE at 426–428, gatekeeper
    T334, and UniProt ACT_SITE 363 = author 382.
- **Excluded residues:**
  - 243–247: these are a GAMDPS cloning tag in 1OPJ but native linker residues in 1OPK.
  - 530–534: unobserved in 1OPJ A.
- **Primary pair:** 1OPJ A vs 1OPK A, chosen a priori as the first deposited chain of each entry. 1OPJ B was added as a
  second copy; 1OPK has one chain.
- **Policies:**
  - `on_missing: reject`: no gaps were expected, and none were found.
  - `resname_mismatch: flag`: for D382N only.
  - Altloc `blank_then_A`: neither file has any altloc atoms.
  - `altloc_sensitivity: highest_occupancy`: reported as not applicable.
- **Advisory:** `min_difference {"disp": 1.0}` Å, with defaults otherwise (statistic `min`, within_factor 1). It was
  chosen as a conventional screening level, well above the expected coordinate error at 1.75–1.8 Å resolution, and is
  not a significance threshold. No other channel was configured.
- **Run:** one attempt; the command was
  `PYTHONPATH=python python -m topofold_inputs.structpair ../specs/abl_1OPJ_AB_vs_1OPK_A.json --out ../outputs/abl_1OPJ_AB_vs_1OPK_A`.
  - Exit code 0, 6.7 s wall, 14.5 CPU-s user.
  - Stdout: `282 positions; global Cα RMSD 5.121 Å over 282 Cα; 0 exclusion rows; within-entry variation: partial`.
  - Nothing was changed after the run.
- **Manifest:** schema `topofold_inputs.structpair/3`. All 9 outputs match: `windows.csv`, `residues.csv`,
  `exclusions.csv`, `summary.json`, `figure.png`, `chain_variation.csv`, `chains.csv`, `chain_variation.png` and
  `report.md`.
- **Requested outputs:** all were produced. Some values are unavailable rather than missing outputs:
  - Within-B (1OPK) variation for every position and channel, because 1OPK has one chain.
  - The 3-residue metrics at 248 and 529, and the 4-residue metrics at 248, 528 and 529 (chain ends).
  - The altloc sensitivity (no altlocs).

  `exclusions.csv` has a header and no rows.

## 3. Viewer

- **Launch:** the same documented command, run through `.claude/launch.json` from `trial/clone`:
  `streamlit run apps/structpair_viewer/app.py --server.address localhost --server.port 8517 --server.headless true -- --results ../outputs`.
  - My first attempt, on port 8501, failed with "Port 8501 is not available", because another Streamlit process
    already held that port. I did not touch that process.
- **Browser:** the Claude desktop app's built-in browser, which is a real browser with interaction by
  click/type/scroll. It reports "WebGL is not supported", so the 3D Cα view never rendered.
- **Screenshots:** 28 captioned screenshots in `screenshots/` (`CAPTIONS.md`).

| Check | Result |
|---|---|
| Two metrics and units | d_kappa in rad and rmsd3 in Å, each on its own axis. The hover at window A:405..A:407 gives 0.6709 Å, which matches `windows.csv` (s03) |
| Sort and filter | Sorting by rmsd3 descending puts A:410..A:412 first, as in the report (s04). Filtering by "name mismatch in window" shows rows 380–383 and notes that the focus row is hidden (s09) |
| Select a window | Clicked row 161 (I aimed at 158; my mis-click). The immediate state still showed the old sidebar focus (s05). After ~3 s the row, sidebar (161 · A:409), details (408–410, rmsd3 0.5811 Å = CSV) and 3D trace data (A:408–410; DOM only) agreed (s06, s07). Second selection: row 134 (A:382) gave the same pattern (s10 → s11) |
| Switch to KRAS and back | KRAS restored its own focus, 79 (s12). On returning, the immediate state showed the KRAS body under the Abl combobox (s13). Settled, the Abl focus was 134, consistent across the row, sidebar, details and 3D data (s14) |
| Exclusion | Abl has none ("No window inside the mapped range was excluded", s20). KRAS row: 3-residue window at position 58 (A:58..A:68), cause `numbering_discontinuity`, location A:59→A:68 (s21, s22). In plain language: residues 60–67 are not in the correspondence (60–64 are also unobserved in 5US4), so A:59 and A:68 are not chain neighbours, and any window joining them would bridge a gap |
| Correspondence and chain variation | same_numbering over 282 positions, 0 missing, 1 name mismatch (s18). The chain-variation tab shows the pair counts (cross 2, withinA 1, withinB 0) and between and within statistics. Unavailable within-B values are shown as grey "None" (s15–s17) |
| Provenance | Shows the spec SHA-256 (= my pre-run hash), the input hashes and a per-output manifest tick, and states that "1OPJ.pdb/1OPK.pdb SHA-256 matches the manifest" (s19) |
| Tamper test | Copied the result to `tamper/abl_tampered` and changed rmsd3 at A:409 from 0.581117 to 0.981117 without updating the manifest. My first `sed` matched nothing and is logged. The viewer refused the copy: "files differ from the manifest hashes: windows.csv" (s23). The original still passes 9 of 9 (`logs/09_tamper.log`) |

## 4. Interpretation

This interpretation was written without opening any TopoFold Abl report. "disp" is the Cα displacement after the
tool's trimmed fit. Values are taken from `chain_variation.csv` (between-entry min / within-1OPJ max, in Å),
`windows.csv` and `report.md`.

1. **Largest local differences.**
   - **Location.** Residues 400–416 (DFG and activation loop; the region includes Y412) have the largest disp:
     - 405: 21.0 / 0.35
     - 409: 32.9 / 1.08
     - 413: 17.8 / 0.47
     - 416: 3.9 / 0.48
   - **Independent check.** A numpy fit on a core excluding 395–420 (core Cα RMSD 0.85 Å) gives the same profile, with
     a peak of 33.5 Å at 409 (`logs/06_indep_check.log`). Consecutive Cα–Cα distances are 3.77–3.84 Å in all three
     chains, so the displacement is not a numbering artefact.
   - **Window metrics.** All six metrics peak in the same region (`report.md` top-10 table): rmsd3 0.789 Å at
     A:410..412; d_i_i3 4.555 Å at A:403..406; d_tau 3.03 rad at A:410..413.
   - **How the metrics differ.** The window metrics are modest compared with disp. A ~33 Å displacement corresponds
     to only ~0.6–0.8 Å of 3-residue local RMSD. The loop is relocated mostly as a body with local re-folding, rather
     than being locally distorted everywhere. This is what local and global measures are expected to differ on, and
     it is not a contradiction.
   - **Metric agreement.** The three-residue metrics agree closely (Spearman 0.94–0.98). d_tau vs d_i_i3 agree less
     (0.53), because d_tau sees handedness and d_i_i3 does not.
   - **Smaller signals:**
     - β3–αC, 292–297: up to 5.97 / 1.23 at 297; d_tau 2.92 rad at A:296..299.
     - Glycine-rich loop, 264–273 (GGGQYG at 268–273): 1.0–2.6, with within-1OPJ values of 0.26–1.87 Å.
     - N-terminus, 248–252: 4.3–17.2 / 1.6–5.3.
2. **Chains within 1OPJ.**
   - Copy-to-copy differences are small in the activation loop (≤1.39 Å) and the core (my core RMSD between A and B
     is 0.60 Å). This makes the ~33 Å loop difference clearly larger than anything seen between the 1OPJ copies.
   - The copies differ by up to 5.3 Å at the N-terminus and up to 1.9 Å in the glycine-rich loop, so the advisory
     flags there are weak evidence.
   - Confidence is limited because there is only one within-entry pair. 1OPK's copy variation is unknown (one chain),
     so the "within" comparison uses 1OPJ only (`within_scope A_only`).
3. **Advisory status.**
   - For disp only, the rule is met at 45 positions and not met at 237; the other 6 channels are not assessed
     (`report.md` table).
   - Being met means only that the minimum over the two cross-pairs is ≥ max(1 Å, the 1OPJ A–B difference).
   - It is not significance, it gives no information about 1OPK copy variation, and some flags sit near the threshold
     (A:277 = 1.002 Å).
4. **Exclusions.**
   - None: all 282 positions are present in all three chains and there are no altlocs.
   - Values are unavailable only at the chain ends, where no complete window exists. These are not listed as
     exclusions.
5. **Justified conclusion and limits.**
   - **Justified:** in these two deposited models, the Abl activation loop (≈400–416) occupies a substantially
     different position (up to ~33 Å Cα displacement after a core fit), while the remainder of the kinase domain
     superposes to under 1 Å. The 1OPJ copies agree with each other in that loop.
   - **Confounds.** The entries differ simultaneously in:
     - construct: isolated kinase domain + tag vs SH3-SH2-kinase;
     - ATP-site ligand: imatinib (STI) vs P16;
     - catalytic-site mutation: D382N in 1OPK only (HRD→HRN);
     - space group: P1 with 2 copies vs C2221 with 1 copy;
     - crystallization: pH 5.6, MES, MgCl2, 277 K vs pH 7.0, KNO3, 293 K;
     - resolution: 1.75 vs 1.80 Å.

     Both entries contain myristic acid. The loop difference therefore cannot be attributed to the inhibitor, the
     mutation, the regulatory domains or packing.
   - **Not inferable:** populations, dynamics, functional activity, the effect of any drug, or any clinical or
     treatment relevance.
6. **Did the tool help?**
   - Yes for locating and quantifying the difference, and for showing it against copy variation, with clear caveats
     throughout.
   - The 5.12 Å global RMSD on the header card is dominated by about 15 loop residues. It is labelled "context only",
     but it is the most prominent number on the page.

## 5. Findings

### Reproducible defect

- **D1. Sidebar focus desynchronises after a results-folder change**, reproduced 2 of 2 times.
  - Steps (KRAS): set focus to 79 or 39 → change *Results folder* to `../tamper` → change it back to `../outputs`.
  - Settled after more than 8 s, and visible on screen: the sidebar reads "0 · A A:1 ↔ B A:1", while the details
    panel ("Focus window at position 79/39"), the plot marker and the 3D trace data (DOM check) stay on the old window
    (s24, s25, s28).
  - The next comparison switch resynchronises to position 0 and loses the old focus (s26).

### Confusing or potentially misleading output

- The exclusion `detail` says "break between adjacent output residues 58 and 59". These are row indices, while the
  location field gives A:59→A:68.
- The explanation column is truncated in the grid, and I found no way in the UI to read it in full.
- The plot x-axes are "window position" indices, not residue numbers; hovering is needed to map them.
- The `chain_variation.png` axis label says "one tick label per position", but the ticks are every ~35 positions.
- Unavailable values appear as grey "None" in the variation grid but as blank cells in the correspondence grid. The
  guide says "gaps or empty cells". Neither is zero.
- "No window inside the mapped range was excluded" is correct, but the chain-end positions without values are not
  mentioned there.
- When switching comparisons, the row filter carries over, which can produce an empty table (s12).
- Transient states during a rerun are visible for about 1–3 s (s05, s10, s13). This is expected per guide §6, and
  each one settled correctly.

### Documentation gaps

- The expected stdout lacks the "within-entry variation" clause.
- The "byte-identical summary.json" claim no longer holds against the shipped reference.
- §3a says the manifest schema is "…/2", but the outputs are `/3`.
- No manifest-verification command is documented.
- `fetch_inputs.py` writes into the clone.
- The guide does not mention that the viewer's 3D view needs WebGL, or what to do when port 8501 is busy.

### Environmental failures

- No WebGL in the browser pane, so the 3D view was not testable.
- Port 8501 was occupied.

### My own mistakes

- Two mis-clicks: one hit the wrong row, and one opened the get.webgl.org link, which reset the session.
- My first tamper `sed` was a no-op.
- I first put `launch.json` in the wrong directory.

## Artefacts

All artefacts are in the trial directory (`trial/`).

- `REPORT.md`
- `notes/`: pre-registration, inspection and hashes; `notes/output_hashes.txt` holds the SHA-256 of all 17 result
  files.
- `specs/abl_1OPJ_AB_vs_1OPK_A.json`, SHA-256 `0343bbf7…321c`.
- `inputs/` (PDB, mmCIF, `SHA256SUMS`).
- `outputs/kras2`, `outputs/abl_1OPJ_AB_vs_1OPK_A`.
- `tamper/abl_tampered`.
- `logs/01–09`.
- `screenshots/` + `CAPTIONS.md`.
- `scripts/` (self-written checks).
- `clone/` (unchanged).
- `topofold-env/`.
