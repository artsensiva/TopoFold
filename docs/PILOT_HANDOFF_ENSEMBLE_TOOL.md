# Pilot Handoff: Descriptive Ensemble Comparison (TopoFold v0.2.0 tool)

For a computational scientist who wants to inspect, without writing code, where local backbone geometry differs
between groups of MD trajectories. The tool is descriptive. It gives no p-values, significance labels or uncertainty
intervals, and it is not validated for inference (`benchmarks/audit/synth_daboot/CLOSEOUT_2026-10-23.md`).

## 1. The question it helps inspect

> Within one protein, which local backbone windows (three- and four-residue Cα windows) have different conformational
> distributions between two groups of trajectories, and are those differences larger than the differences between
> trajectories of the same group?

**Channels:**
- κ: Cα turning angle, a three-residue window;
- τ: Cα virtual torsion, a four-residue window;
- the equal-information Cα–Cα distances d(i, i+2) and d(i, i+3).

All four are computed from Cα only, so side-chain changes are invisible.

## 2. Requirements

- **Repository checkout**, with:
  - `.venv`: Python with NumPy and the `topofold` extension (built with `maturin develop --release`);
  - an interpreter with **MDAnalysis** (and **h5py** for mdCATH files). In this repository that is
    `benchmarks/audit/wpd_local_retrieval/build/venv/bin/python`. Pass it as `--adapter-python`.
- **Inputs:** either MDAnalysis-readable trajectory files (topology plus trajectories) or mdCATH HDF5 files.
  - Every member must have the same residues (chain, number, insertion code, residue name) in the same order.
  - Coordinates are converted to Å by MDAnalysis. For mdCATH, the file's unit attribute must be "Angstrom".

## 3. Minimal input and exact command

The input is a JSON file with named groups of members. Each member is one trajectory or one mdCATH replica, with an
optional frame slice (`start`, `stop`, `step`) and a PBC policy. Example
(`docs/pilot/spec_example_1ha8A00.json`, three replicas per temperature):

```json
{"name": "Pilot example: mdCATH 1ha8A00, 320 K vs 348 K, 3 replicas each",
 "groups": {
   "320K": [{"type": "mdcath", "h5": "benchmarks/audit/wpd_local_retrieval/data/mdcath/mdcath_dataset_1ha8A00.h5",
             "domain": "1ha8A00", "temperature": "320", "replica": "0", "start": 43, "stop": 430,
             "pbc": "make_whole_sequential"}, "... replicas 1 and 2 ..."],
   "348K": ["... the same with temperature 348 ..."]}}
```

For ordinary MD files, use instead:

```json
{"type": "mdanalysis", "topology": "system.gro", "trajectories": ["run1.xtc"], "start": 0, "stop": null, "step": 1,
 "pbc": "make_whole_sequential"}
```

Run from the repository root. Put `--out` outside the checkout; rerunning into an existing directory overwrites it
without warning.

```bash
PYTHONPATH=python .venv/bin/python -m topofold_inputs.ensemble docs/pilot/spec_example_1ha8A00.json \
    --out ~/pilot_out --adapter-python benchmarks/audit/wpd_local_retrieval/build/venv/bin/python
```

To list the replicas and frame counts in an mdCATH file before choosing `start`/`stop` (use the same `stop` for
every member if you want equal frame counts):

```bash
benchmarks/audit/wpd_local_retrieval/build/venv/bin/python -c "import h5py,sys; f=h5py.File(sys.argv[1]); d=sys.argv[2]; \
print({T: {r: int(f[d][T][r].attrs['numFrames']) for r in f[d][T]} for T in f[d] if T.isdigit()})" FILE.h5 DOMAIN
```

**Expected cost:** for this example, about 20 s and 110 MiB. Five replicas per group took about 30 CPU-s on the
development machine. Cost grows roughly linearly with frames × residues × members.

## 4. Outputs (`--out` directory)

| File | Content |
|---|---|
| `report.md` | Human-readable summary: inputs, exclusions, the largest differences per channel with replica context, temporal diagnostics |
| `tables/distribution.csv` | Per channel, window and group: frames, per-member frame counts, center, spread, quantiles (linear), resultant length (circular), status |
| `tables/replicas.csv` | Per group: spread of member centers, and median and maximum JS distance between members of the same group |
| `tables/between.csv` | Per group pair: pooled JS distance, median between-member JS distance, within-member medians for each group, center difference |
| `tables/temporal.csv` | Per member: lag-1 autocorrelation, decorrelation lag (first lag with C(k) < 1/e), status |
| `tables/exclusions.csv` | Windows excluded because any member had invalid frames (missing residue, chain or numbering break, degenerate geometry), with frame counts per member |
| `distributions.npz` | Pooled and per-member histograms for every channel, window and group, with bin edges |
| `manifest.json` | Inputs with SHA-256, member specifications, tool version, thresholds, native-extension hash |
| `members/` | Per-member intermediate coordinates and features (can be deleted; can be larger than the final outputs). `members/<run_id>.json` → `breaks` gives the specific reason for an exclusion (e.g. `numbering_discontinuity`, `chain_boundary`, `missing_CA`, `ca_ca_distance`) |

## 5. Reading one report: worked example

From `benchmarks/audit/dev_ensemble_report/results_v2/1ha8A00/report.md` (5 replicas per temperature):

- **Exclusions: 0.** Every window had valid frames in all members. If a residue is missing, or the numbering jumps,
  windows touching it appear in `tables/exclusions.csv` and nowhere else.
- **Largest differences, τ 32–35:**

  | Pooled JS | Between-member JS (median) | Within-member JS (median) | Center difference |
  |---|---|---|---|
  | 0.350 | 0.614 | 0.564 (320 K) / 0.393 (348 K) | +1.632 rad |

  - **Pooled JS** compares all 320 K frames with all 348 K frames. Members are weighted by their frame counts.
  - **Between-member JS** compares single replicas across temperatures. It is about as large as the within-member
    median at 320 K, so individual replicas differ from each other about as much as the two temperatures do. The
    pooled difference here cannot be separated from replica-to-replica variability by eye.
  - **Center difference** is the wrapped difference of circular means. It is shown only when both resultant lengths
    are numerically well conditioned. A large center difference with low resultant lengths describes broad,
    overlapping distributions; it does not mean a rigid shift.
- **Temporal diagnostics:**
  - τ has a median lag-1 autocorrelation of 0.699 (320 K) and 0.633 (348 K), with decorrelation lags of about 10 and
    7 frames.
  - 387 frames per replica therefore hold far fewer effectively independent samples than 387.
  - "Not reached within L lags" counts member-windows whose autocorrelation stayed above 1/e up to the last lag
    examined.

## 6. Pooled differences vs replica variability vs sampling

- **Pooled differences alone are not enough.** Compare the between-member medians with the within-member medians in
  `tables/between.csv`. If they are similar, the pooled difference is within the range of replica-to-replica
  variation.
- **Replica spread in `tables/replicas.csv`.** A single unusual member raises the within-member maximum.
- **Correlation in `tables/temporal.csv`.** Strong correlation means fewer effective samples. Different correlation
  between groups can create histogram differences even when the underlying distributions are the same (SYNTHETIC
  study NULLSTRESS-1).
- **Single-member groups** have no replica variability; their columns read n/a. A difference between single runs
  confounds condition, starting structure and sampling.
- **Description, not inference.** No p-values. Differences do not establish mechanism, equilibrium or biological
  function.

## 7. Pilot success, defined operationally

A pilot run counts as successful only if the user, **without editing the implementation**:
1. prepares an input specification for their own groups;
2. completes a run to a `report.md`;
3. identifies the excluded regions and the reasons for them from `tables/exclusions.csv`;
4. explains the output's limitations: descriptive only; pooled vs replica variability; sampling and correlation;
   Cα-only.

Record obstacles (commands, error messages, unclear outputs) in `docs/pilot/PILOT_LOG.md`. Entry 1 there is an
internal walkthrough by the implementer, not independent user validation.
