# Comparing two experimental structures with `topofold_inputs.structpair` — guide for outside users

**What it does.** It compares two single-model structures, one chain from each, residue by residue on a correspondence
that you state. For every 3- and 4-residue window it reports:
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

Status: v0.1.0. It has been verified against two recorded KRAS studies
(`benchmarks/audit/structpair_demo/README.md`). One external AI technical trial (checkpoint `6ce7b49`) and one bounded
AI confirmation (`049a47f`) have been completed (`docs/pilot/PILOT_LOG.md`, entries 4 and 5). **No independent human
trial has been completed.** If you are trying it, please use `docs/pilot/FEEDBACK_FORM.md`.

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
pip install streamlit==1.64.0 pandas==3.0.6 altair==6.3.0 plotly==7.1.0
streamlit run apps/structpair_viewer/app.py --server.address localhost -- --results ~/structpair_out
```

- **Expected stdout of step 4:** `158 positions; global Cα RMSD 0.974 Å over 158 Cα; 10 exclusion rows; wrote …/report.md`.
  The 10 rows are 5 windows that span the 59→68 numbering gap, each listed once per structure.
- **Reference outputs:** `benchmarks/audit/structpair_demo/results/kras2/`. In the trial, the three CSV files and
  `summary.json` were byte-identical to them.
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
- **Residue names:** `resname_mismatch: "reject"` stops on any difference. `"flag"` allows it (for a point mutant,
  for example) and marks every affected window.

## 4. Resource use (measured on the implementer's 6-core Linux machine)

- **Per run:** about 10.6 CPU-s and 5.7 s wall for two ~160-residue chains when no second altloc policy is needed.
  About 14.9 CPU-s and 6.2 s wall when altloc sensitivity is computed, because both structures are extracted twice.
  These were measured in the internal smoke test; most of the time is interpreter and MDAnalysis start-up.
- **Disk:** about 0.6–0.9 MB per output directory, most of it in `work/` intermediates.
- **Memory:** small, well under 1 GB.

## 5. Reading the outputs

- **Files:**
  - `report.md`: the inputs, mapping, exclusions, largest windows, metric agreement and altloc sensitivity.
  - `windows.csv`: one row per window position; units are in the column names (`_rad`, `_angstrom`).
  - `residues.csv`: the mapping, residue status and altloc used.
  - `exclusions.csv`: one row per excluded window, structure and cause.
  - `manifest.json`: the hashes of inputs and outputs and the software versions.
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

Please record the following in `docs/pilot/FEEDBACK_FORM.md`, or send them to the project owner:
- whether you completed the setup and a run on your own structures;
- where you got stuck;
- whether you could explain every exclusion from the outputs alone;
- what you would need before using the output in your own work.
