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
(`benchmarks/audit/structpair_demo/README.md`). **No independent user has run it yet.** Your trial is that test; please
use `docs/pilot/FEEDBACK_FORM.md`.

## 1. What a clean checkout does not contain

| Needed | Why it is absent | How to get it |
|---|---|---|
| The compiled `topofold` Python extension | Build product | Rust toolchain + `pip install maturin`, then `maturin develop --release` in a virtual environment (README "Installation"). Expect a few minutes of compile time |
| An interpreter with **MDAnalysis** (≥ 2.x; 2.10.0 tested) | Optional dependency, not in the core package | Either install MDAnalysis into the same environment, or use a separate environment and pass it with `--adapter-python` / `TOPOFOLD_ADAPTER_PYTHON`. Inside this repository the implementer used `benchmarks/audit/wpd_local_retrieval/build/venv`, a local, git-ignored environment that you will not have |
| **matplotlib** (optional; 3.11.2 tested) | Optional | Without it the figure is skipped and the report says so |
| The example structures | PDB files are git-ignored | `python benchmarks/audit/structpair_demo/fetch_inputs.py` downloads 9O0N, 9O0O, 5US4 and 5USJ from RCSB (about 2.4 MB, CC0) and checks their SHA-256 |

**Tested only with:** Python 3.14.4, NumPy 2.5.3 and Linux x86-64, using PDB-format files.
- mmCIF input depends on the MDAnalysis reader available to you and has not been tested.
- Other Python versions are untested. The README states Python 3.9+ for the core.

## 2. Run the example

From the repository root, with `PY` = your interpreter that has the `topofold` extension and `ADAPTER` = your
interpreter with MDAnalysis (they can be the same):

```bash
$PY benchmarks/audit/structpair_demo/fetch_inputs.py
PYTHONPATH=python $PY -m topofold_inputs.structpair \
    benchmarks/audit/structpair_demo/specs/kras2_5US4A_vs_5USJA.json \
    --out ~/structpair_out/kras2 --adapter-python $ADAPTER
```

- **Expected output on stdout:** `158 positions; global Cα RMSD 0.974 Å over 158 Cα; 10 exclusion rows; wrote …/report.md`
  (5 windows span the 59→68 numbering gap; each is listed once per structure).
- **Reference outputs to compare against:** `benchmarks/audit/structpair_demo/results/kras2/`.
- `--out` must be a new or empty directory, so each run keeps its own record.

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

"TopoFold Structural Comparison" is a local, read-only viewer for structpair output directories:

```bash
.venv/bin/streamlit run apps/structpair_viewer/app.py --server.address localhost -- --results ~/structpair_out
```

- `--results` takes one output directory or a folder containing several. Without it, the viewer opens the two KRAS
  demonstrations in `benchmarks/audit/structpair_demo/results/`.
- `--server.address localhost` keeps the page on your own machine. Nothing is uploaded, and no external script is
  loaded.
- **What it shows:**
  - one plot per metric, each with its own axis and unit;
  - the window table, which can be filtered, sorted and selected (selecting a row moves the focus window);
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
- **Dependencies:** Streamlit, pandas, Altair, Plotly and NumPy (tested with Streamlit 1.64.0, pandas 3.0.6,
  Altair 6.3.0, Plotly 7.1.0 and NumPy 2.5.3). It needs neither the compiled TopoFold core nor MDAnalysis.
- **Older app:** `apps/streamlit_app.py` is a historical dashboard from before the re-validation. Its claims are not
  current; do not use it to read these results.

## 7. Feedback

Please record the following in `docs/pilot/FEEDBACK_FORM.md`, or send them to the project owner:
- whether you completed the setup and a run on your own structures;
- where you got stuck;
- whether you could explain every exclusion from the outputs alone;
- what you would need before using the output in your own work.
