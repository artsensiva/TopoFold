# Comparing two experimental structures with `topofold_inputs.structpair` — guide for outside users

**What it does.** It compares two single-model structures, one primary chain from each (and, if you name them, several
chains per entry: §3a), residue by residue on a correspondence that you state. For every 3- and 4-residue window it reports:
- Cα curvature κ and torsion τ for each structure (from the TopoFold core), and their differences;
- the Cα–Cα distances spanning the same windows, as equal-information baselines;
- the proper-rotation local Cα RMSD of the same windows (MDAnalysis QCP).

All of these use identical coordinates and residue coverage. It also reports:
- global Cα RMSD, with an independent numpy check;
- every excluded window, with the cause;
- alternate-conformation handling;
- a report and a figure;
- hashes of the inputs, software and outputs.

**What it does not do.**
- It does not align sequences, choose chains or rank chains.
- It does not find pockets, compute statistics or produce p-values.
- It makes no claim about ensembles, populations, kinetics, functional state or treatment response. Two static
  structures cannot support such claims.

Status: v0.3.0. Version 0.1.0 was verified against two recorded KRAS studies
(`benchmarks/audit/structpair_demo/README.md`); version 0.2.0 added explicitly named multiple chains and 0.3.0 corrected their advisory output (schema /3); every
single-chain output is unchanged (`benchmarks/audit/structpair_chains_demo/README.md`). One external AI technical trial (checkpoint `6ce7b49`) and one bounded
AI confirmation (`049a47f`) have been completed, and two further independent AI technical trials of 0.3.0 (`19b5268`)
(`docs/pilot/PILOT_LOG.md`, entries 4–6). **No independent human trial has been completed.** If you are trying it,
please use `docs/pilot/STRUCTPAIR_HUMAN_TRIAL_PACKET.md` and `docs/pilot/STRUCTPAIR_FEEDBACK_FORM.md`.

## 1. What a clean checkout does not contain

| Needed | Why it is absent | Where it comes from in §2 |
|---|---|---|
| The compiled `topofold` Python extension | Build product | `maturin develop --release` (needs a Rust toolchain) |
| **MDAnalysis** | Optional dependency, not in the core package | `pip install` into the same environment, which is the simplest option. A separate interpreter can be named with `--adapter-python` or `TOPOFOLD_ADAPTER_PYTHON` instead |
| **matplotlib** (optional) | Optional | `pip install`; without it the figure is skipped and the report says so |
| Viewer packages (optional) | Optional | `pip install` (§2, step 5) |
| The example structures | PDB files are git-ignored | `fetch_inputs.py` downloads 9O0N, 9O0O, 5US4 and 5USJ from RCSB (about 2.4 MB, CC0) and checks their SHA-256 |

## 2. Isolated install and the example (one environment)

This sequence matches the commands that succeeded in the AI technical trial of 2026-10-24 (AI tester 1, checkpoint
`6ce7b49`, Linux x86-64; `docs/pilot/PILOT_LOG.md` entry 4). The environment can live anywhere; nothing assumes a
`.venv` inside the checkout. Run everything from the repository root.

```bash
# 1. code
git clone --branch scientific-rebuild https://github.com/artsensiva/TopoFold.git
cd TopoFold

# 2. an isolated environment (any path; here a sibling directory)
python3 -m venv ../topofold-env
source ../topofold-env/bin/activate

# 3. the compiled core (Rust toolchain required; about 1 minute on the trial machine)
pip install maturin
maturin develop --release

# 4. computation dependencies, the example inputs and the example comparison
pip install MDAnalysis==2.10.0 matplotlib==3.11.2
python benchmarks/audit/structpair_demo/fetch_inputs.py
PYTHONPATH=python python -m topofold_inputs.structpair \
    benchmarks/audit/structpair_demo/specs/kras2_5US4A_vs_5USJA.json --out ~/structpair_out/kras2

# 5. optional: the results viewer (needs neither the core nor MDAnalysis)
pip install streamlit==1.64.0 pandas==3.0.6 altair==6.3.0 plotly==7.1.0   # same pins: apps/structpair_viewer/requirements.txt
streamlit run apps/structpair_viewer/app.py --server.address localhost -- --results ~/structpair_out
```

- **Expected stdout of step 4 (version 0.3.0):** `158 positions; global Cα RMSD 0.974 Å over 158 Cα; 10 exclusion rows; within-entry variation: unavailable; wrote …/report.md`.
  The `within-entry variation` clause is new since 0.2.0 (it is `unavailable` with one chain per entry).
  The 10 rows are 5 windows that span the 59→68 numbering gap, each listed once per structure.
- **Reference outputs:** `benchmarks/audit/structpair_demo/results/kras2/`. The three CSV files (`windows.csv`, `residues.csv`,
  `exclusions.csv`) are byte-identical to them. `summary.json` has the same content plus one added block,
  `within_entry_variation` (version 0.2.0 and later; the earlier 0.1.0 trial found it byte-identical). Check your own run with §3b.
- **`fetch_inputs.py` writes into the clone** (git-ignored `data_pdb/` folders next to each demonstration); it has no output-directory
  option. To keep downloads elsewhere, run a copy of it from another directory.
- `--out` must be a new or empty directory, so each run keeps its own record.
- `PYTHONPATH=python` makes the `topofold_inputs` package importable. It is not installed as a package.
- **Versions** (an exact record, not a compatibility range):

  | Component | Used in the trial of 2026-10-24 |
  |---|---|
  | Python | 3.14.4 |
  | Rust | 1.98.1 |
  | maturin | 1.15.0 |
  | NumPy | 2.5.4 (resolved by the build) |
  | MDAnalysis | 2.10.0 |
  | matplotlib | 3.11.2 |
  | Streamlit | 1.64.0 |
  | pandas | 3.0.6 |
  | Altair | 6.3.0 |
  | Plotly | 7.1.0 |

  The implementer's own checks used the same versions except NumPy 2.5.3. In the trial, installing the viewer packages
  changed Pillow from 13.0.0 to 12.3.0 after the computation had finished.
- **Untested:**
  - other Python, Rust and package versions, and platforms other than Linux x86-64. The README notes that the
    declared Rust minimum (1.80) is not honoured (KNOWN_ISSUES m12);
  - a viewer-only environment without the core;
  - mmCIF input, which depends on the MDAnalysis reader available to you.

## 3. Write your own spec

```json
{"name": "my comparison",
 "structures": {"A": {"label": "apo", "file": "apo.pdb", "chain": "A", "altloc": "blank_then_A"},
                "B": {"label": "holo", "file": "holo.pdb", "chain": "B", "altloc": "blank_then_A"}},
 "correspondence": {"same_numbering": "5-120,125-180"},
 "on_missing": "reject",
 "resname_mismatch": "reject",
 "altloc_sensitivity": "highest_occupancy"}
```

- **Paths:** relative paths are resolved against the spec file's directory.
- **Chain:** required for each structure. Choose it *before* looking at any result, for example the first deposited
  chain or the chain named in the paper.
  - To check chain dependence, run each additional pairing as its own spec and report all of them.
  - Never keep only the pairing that looks best.
- **Correspondence:** use `{"same_numbering": "ranges"}` when both entries use the same author numbering. Otherwise
  use `{"pairs": [[A_resid, B_resid], …]}`, with insertion codes written as `"52A"`. Derive the correspondence from
  the deposited sequences or a documented alignment, not from the results.
- **Missing residues:**
  - `on_missing: "reject"` stops at the first absent residue.
  - `"exclude"` keeps it as a missing row, and every window touching it is excluded with the cause named.
- **Alternate locations:** set a policy per structure. `reject` stops if any mapped Cα has an altloc. Set
  `altloc_sensitivity` to report how much each metric changes under a second policy.
  **Known limitation (KNOWN_ISSUES m13):** `highest_occupancy` stops the whole run if a mapped Cα has alternate locations with equal
  occupancy (for example A and B at 0.5). In 4 of the 6 pairs of KIN-PAIR-1 it did. If you meet "tie on occupancy", remove
  `altloc_sensitivity` from the spec (the primary policy is unaffected) or choose `blank_then_A` and compare by hand.
- **Residue names:** `resname_mismatch: "reject"` stops on any difference. `"flag"` allows it (for a point mutant,
  for example) and marks every affected window.

## 3a. Several chains per entry: within-entry crystal-copy variation (v0.2.0; advisory output corrected in v0.3.0)

> **Read this first.** Chains in one crystal are **not independent experimental replicates**, and crystal packing can affect
> each of them. The within-entry crystal-copy variation is a **descriptive range** of differences among the chains you
> selected. It is **not** a calibrated noise floor, **not** a statistical uncertainty, **not** a nucleotide or ligand effect and
> **not** an HDX replicate estimate. An advisory flag exists only if you configure one (below), and it is an advisory rule: not a gate,
> not an adequacy verdict, not significance and not biological separation.

Name every chain, and one primary chain, per entry (no chain is chosen for you):

```json
"structures": {"A": {"label": "5US4", "file": "5US4.pdb", "chains": ["A", "B"], "primary_chain": "A"},
               "B": {"label": "5USJ", "file": "5USJ.pdb", "chains": ["A", "B", "C", "D", "E", "F"], "primary_chain": "A"}},
"correspondence": {"same_numbering": "1-59,68-166"},
"chain_variation": {"displacement_trim_fraction": 0.2,
                    "advisory": {"statistic": "min", "within_factor": 1.0, "min_difference": {}}}   // the defaults: numbers only, no flag
```

- **The correspondence applies identically to every chain of an entry** (same author numbers and insertion codes). The
  missing-residue, alternate-location and residue-name rules apply per chain; chains of one entry must agree in residue
  name (or set `resname_mismatch: "flag"`).
- **The legacy tables are the primary pair** (primary chain of A against primary chain of B). They are identical to a
  run with only those two chains. With one chain per entry nothing changes, and the report says that within-entry
  variation is unavailable.
- **New files:** `chain_variation.csv` (long format, one row per position and channel), `chains.csv` (every chain with its
  missing residues and altloc decisions) and `chain_variation.png`. The summary, the manifest (schema `…/3` since version 0.3.0; `…/2` for 0.2.0 results) and the report
  carry the same statement; the viewer shows it on a tab and in a banner, next to the unchanged exclusion and provenance views.
- **Channels:** the six existing window metrics, computed for every chain pair, and `disp`, the Cα displacement after a trimmed fit
  (the K1 calculation). Per position and channel you get the between-entry values over all A × B pairings (count, min,
  median, mean, max) and, per entry, the minimum and maximum over the pairs inside that entry.
- **Missing values:** a pair contributes only where both chains have the data; the counts (`*_n`) say how many pairs did. An
  empty cell is unavailable, never zero.
- **Advisory rule (optional).** By default you get the numbers and **no flag**: it stays unassessed until you set a finite, positive
  `min_difference` for a channel, in that channel's unit (for example `{"disp": 1.0}` for 1 Å). Then `advisory_flag` is true where
  the between-entry statistic (default `min`) is at least max(`min_difference`, `within_factor` × the within-entry maximum). **Equality
  meets the rule**; a value a hair below does not; `advisory_at_boundary` marks values within 1e-9 of the threshold so you can see when a
  flag hinges on rounding. It is empty, never false, where it cannot be assessed (no minimum difference configured, no cross pairing,
  no within-entry pair). All numbers must be finite, and a minimum difference must exceed 1e-9. It is not significance and not
  biological separation. Older results (schema /2, version 0.2.0) used a default minimum difference of 0 and their column was
  called `separated`: they can show it true for identical chains, and the viewer labels them as historical. The full rules are in
  `docs/STRUCTPAIR_CHAINS_SPEC_2026-10-25.md` (amendment 2).

## 3b. Checking a result against its manifest

`manifest.json` records the SHA-256 of every output file, of both input structures and of the spec. To re-check them (standard
library only; no viewer needed):

```bash
PYTHONPATH=python python -m topofold_inputs.structpair_verify ~/structpair_out/kras2 [--json check.json]
```

It prints `ok`, `mismatch` or `missing` for each file and the counts (9 files for the KRAS example: 6 outputs, 2 inputs, 1 spec),
and exits 0 only if every listed file is present and matches. A file moved since the run (an input or the spec) shows as `missing`:
that is a failure of the check, not a pass. `manifest.json` does not hash itself, so editing the manifest together with a file is not
detected; keep a copy of the manifest's own hash if that matters. The viewer performs the same output check on every page update
(§6).

## 4. Resource use (measured on the implementer's 6-core Linux machine)

- **Per run:** about 10.6 CPU-s and 5.7 s wall for two ~160-residue chains when no second altloc policy is needed.
  About 14.9 CPU-s and 6.2 s wall when altloc sensitivity is computed, because both structures are extracted twice.
  These were measured in the internal smoke test; most of the time is interpreter and MDAnalysis start-up.
- **Disk:** about 0.6–0.9 MB per output directory, most of it in `work/` intermediates.
- **Memory:** small, well under 1 GB.
- **Several chains (v0.2.0):** the KRAS example with 2 + 6 chains (28 chain pairs, 158 positions, 7 channels) took about
  33.7 CPU-s and 20 s wall and wrote 0.75 MB without `work/` (1.7 MB with it). Time grows with the number of chains (one
  adapter run each) and chain pairs.

## 5. Reading the outputs

- **Files:**
  - `report.md`: the inputs, mapping, exclusions, largest windows, metric agreement and altloc sensitivity.
  - `windows.csv`: one row per window position; units are in the column names (`_rad`, `_angstrom`).
  - `residues.csv`: the mapping, residue status and altloc used.
  - `exclusions.csv`: one row per excluded window, structure and cause. In a break record the `detail` text says
    "break between adjacent output residues 58 and 59": those numbers are **row indices** (the `position` column of
    `residues.csv`), not residue numbers. The `location` column gives the author residues (here `A:59→A:68`). The viewer
    shows this translation in a `detail_reading` column, taken from `residues.csv`; the file itself is unchanged.
  - `manifest.json`: the hashes of inputs and outputs and the software versions.
  - With several chains (§3a): `chain_variation.csv`, `chains.csv`, `chain_variation.png`. Read the statement at the top of
    the report section before using them.
- **What a large value means:** the two coordinate sets differ locally there. Before reading biology into it, check:
  1. **Coordinate quality:** resolution, B-factors, altlocs (marked `*`), and incomplete or missing residues nearby.
  2. **Crystal context:** different space groups, packing contacts and crystallization conditions (pH, additives)
     can change loops.
  3. **Chain dependence:** if an entry has several copies, compare copy-to-copy differences within that entry. A
     difference smaller than the copy-to-copy variation is not distinctive.
  4. **Metric agreement:** κ, d(i−1, i+1) and 3-residue RMSD share one triangle and usually agree. τ adds
     handedness, which d(i−1, i+2) cannot see. Disagreement shows what each metric measures, not which one is right.
- **What not to conclude:**
  - That a region is functionally switched, populated or dynamic.
  - That a ligand or mutation *caused* a difference, unless the comparison isolates that factor and the confounds
    above are addressed.
  - That agreement with a known functional region validates the method: familiar regions are where people look first.
- **Global RMSD:** depends on the atom set, the residue range and any outlier rejection. Compare it with a published
  value only when the definition is stated. In this repository, for example, 0.718 Å over 168 Cα could not be
  reconciled with an undefined published 0.3 Å.

## 6. Browse results in the viewer (optional)

"TopoFold Structural Comparison" is a local, read-only viewer for structpair output directories. With the environment
from §2 active, run from the repository root:

```bash
streamlit run apps/structpair_viewer/app.py --server.address localhost -- --results ~/structpair_out
```

- `--results` takes one output directory or a folder containing several. Without it, the viewer opens the two KRAS
  demonstrations in `benchmarks/audit/structpair_demo/results/`.
- **Page updates:** each control change reruns the page once. While a rerun is in progress, Streamlit may briefly
  show elements from the previous state (greyed out). Read the selection, details and plots after the update has
  finished.
- **Port busy:** if port 8501 is taken ("Port 8501 is not available"), add `--server.port 8517` (any free port) and open that address.
- **3D view needs WebGL.** Browsers or panes without WebGL (some embedded browsers, remote desktops, software-rendering
  setups) show "WebGL is not supported" or an empty 3D panel. Nothing is wrong with the result: the plots, the window table and
  the details panel show the same focus window. Open the page in a desktop browser with WebGL enabled to see the 3D view.
- **Changing the Results folder** shows the first comparison of the new folder and keeps the focus remembered for each comparison.
  If a result fails its manifest check the page says so and stops; returning to a valid folder restores the remembered focus
  (fixed after the trial of 2026-10-25, D77; at commit 19b5268 the sidebar could show position 0 while the details kept the old window).
- `--server.address localhost` keeps the page on your own machine. Nothing is uploaded, and no external script is
  loaded.
- **What it shows:**
  - one plot per metric, each with its own axis and unit;
  - the window table, which can be filtered, sorted and selected. Selecting a row moves the focus window, and the
    table always marks the focus row, or no row when a filter hides it. The focus is kept per comparison. The table
    shows only the rows in view, so scroll to see a marked row that is off-screen;
  - the exclusions, each with its recorded cause and an explanation of that cause type;
  - the residue correspondence, keeping chain IDs and insertion codes exactly as in the files;
  - provenance, with every output re-checked against its manifest hash;
  - a 3D Cα view of both structures in which the focus window is highlighted.
- **Missing values:** values that are unavailable are shown as gaps or empty cells, never as zero.
- **Results only:**
  - The viewer computes nothing new. Its one calculation is a display-only superposition for the 3D view.
  - Result files are re-read and re-checked against the manifest on every page update, not only when the manifest
    changes. The manifest must hold a valid SHA-256 entry for `summary.json`, `windows.csv`, `residues.csv` and
    `exclusions.csv`. A result with a missing entry, or with a file that no longer matches, is refused with the file
    named. Optional outputs (`report.md`, `figure.png`) are checked when listed and only noted if absent. The manifest
    does not hash itself.
  - The input structures are re-hashed on every page update too. If one is missing or changed, the 3D view is
    disabled with the reason, and the rest of the viewer still works.
- **Dependencies:** Streamlit, pandas, Altair, Plotly and NumPy (versions used: §2). It does not import the compiled
  TopoFold core or MDAnalysis; a test checks this. A viewer-only environment has not been tested.
- **Older app:** `apps/streamlit_app.py` is a historical dashboard from before the re-validation. Its claims are not
  current; do not use it to read these results.

## 7. Feedback

Please record the following in `docs/pilot/STRUCTPAIR_FEEDBACK_FORM.md`, or send them to the project owner:
- whether you completed the setup and a run on your own structures;
- where you got stuck;
- whether you could explain every exclusion from the outputs alone;
- what you would need before using the output in your own work.
