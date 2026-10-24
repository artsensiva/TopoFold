# Executed commands and deviations

All project commands ran from the fresh checkout root. Exact arguments, environment overrides, timestamps, exit codes and stdout/stderr are in the numbered JSON/log files. Only paths, isolation settings and requested fixed-checkpoint selection differ from the guide setup. The existing system Python and Rust toolchains were used.

## 01-clone

```sh
git clone https://github.com/artsensiva/TopoFold.git /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/confirmation-049a47f/checkout
```

Exit: 0. Evidence: 01-clone.json.

## 02-resolve

```sh
git rev-parse '049a47f^{commit}'
```

Exit: 0. Evidence: 02-resolve.json.

## 03-checkout

```sh
git checkout --detach 049a47f196e39ec46cc8e2c9d0bd63f03600df71
```

Exit: 0. Evidence: 03-checkout.json.

## 04-venv

```sh
python3 -m venv ../topofold-env
```

Exit: 0. Evidence: 04-venv.json.

## 05-maturin

```sh
bash -c 'source ../topofold-env/bin/activate && pip install maturin'
```

Exit: 0. Evidence: 05-maturin.json.

## 06-build

```sh
bash -c 'source ../topofold-env/bin/activate && maturin develop --release'
```

Exit: 0. Evidence: 06-build.json.

## 07-computation-deps

```sh
bash -c 'source ../topofold-env/bin/activate && pip install MDAnalysis==2.10.0 matplotlib==3.11.2'
```

Exit: 0. Evidence: 07-computation-deps.json.

## 08-fetch-inputs

```sh
bash -c 'source ../topofold-env/bin/activate && python benchmarks/audit/structpair_demo/fetch_inputs.py'
```

Exit: 0. Evidence: 08-fetch-inputs.json.

## 09-kras2

```sh
bash -c 'source ../topofold-env/bin/activate && PYTHONPATH=python python -m topofold_inputs.structpair benchmarks/audit/structpair_demo/specs/kras2_5US4A_vs_5USJA.json --out /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/confirmation-049a47f/results/kras2'
```

Exit: 0. Evidence: 09-kras2.json.

## 10-kras1-added

```sh
bash -c 'source ../topofold-env/bin/activate && PYTHONPATH=python python -m topofold_inputs.structpair benchmarks/audit/structpair_demo/specs/kras1_9O0N_vs_9O0O.json --out /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/confirmation-049a47f/results/kras1'
```

Exit: 0. Evidence: 10-kras1-added.json.

## 11-viewer-deps

```sh
bash -c 'source ../topofold-env/bin/activate && pip install streamlit==1.64.0 pandas==3.0.6 altair==6.3.0 plotly==7.1.0'
```

Exit: 0. Evidence: 11-viewer-deps.json.

## 12-viewer

```sh
bash -c 'source ../topofold-env/bin/activate && exec streamlit run apps/structpair_viewer/app.py --server.address localhost -- --results /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/confirmation-049a47f/results'
```

Exit: None. Evidence: 12-viewer.json.

## 13-viewer-resume

```sh
bash -c 'source ../topofold-env/bin/activate && exec streamlit run apps/structpair_viewer/app.py --server.address localhost -- --results /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/confirmation-049a47f/results'
```

Exit: see shutdown record. Evidence: 13-viewer-resume.json.

## browser-saver-resume

```sh
python3 -B /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/confirmation-049a47f/evidence/save_browser_evidence.py
```

Exit: see shutdown record. Evidence: browser-saver-resume.json.

## browser-saver

```sh
python3 -B /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/confirmation-049a47f/evidence/save_browser_evidence.py
```

Exit: None. Evidence: browser-saver.json.

## Every addition/change relative to §2

- Cloned into the explicit trial checkout directory, then added `git rev-parse '049a47f^{commit}'` and `git checkout --detach 049a47f196e39ec46cc8e2c9d0bd63f03600df71`; the guide uses a branch clone and `cd TopoFold`.
- Created the guide's sibling `../topofold-env` and ran its activation in each recorded shell, because tool calls do not share shell activation. Commands use `bash -c` and `&&` to preserve this requirement. No extra package, alternative version or adapter-python option was needed.
- Replaced `~/structpair_out` in the example and viewer command with the absolute trial `results` directory, respecting existing trial-only writes. The literal home-directory output location was not tested.
- Set CARGO_HOME, CARGO_TARGET_DIR, PIP_CACHE_DIR, XDG_CACHE_HOME, MPLCONFIGDIR and TMPDIR to new trial directories; set PYTHONDONTWRITEBYTECODE=1 and PIP_DISABLE_PIP_VERSION_CHECK=1 for isolation/logging. These settings are listed in each command JSON.
- Set STREAMLIT_SERVER_HEADLESS=true and STREAMLIT_BROWSER_GATHER_USAGE_STATS=false for the automated local launch. Invoked the documented `streamlit run apps/structpair_viewer/app.py --server.address localhost -- --results ...` through the activated shell using exec so its PID could be stopped. These are test-harness settings, not a demonstrated dependency workaround.
- Added the same structpair command with shipped `kras1_9O0N_vs_9O0O.json` and `results/kras1`, because §2 creates only kras2 but this requested selection test needs both comparisons.
- Added trial-only `run.py` command logging, copies of AGENTS.md/user guide, and a loopback HTTP screenshot saver. Reused only this tester utility from the previous trial; all TopoFold dependencies/builds/inputs/results were fresh.
- After the user-reported reboot, repeated only the viewer launch and screenshot-saver launch. Restored comparison=kras2 and focus=58 through browser controls; defaults all/position matched the last verified state. Did not rerun install or comparisons.
- Added final Git status/diff/rev-parse and Python/Rust/pip version capture; exact commands and exit codes are in 15-final-checks.json. Stopped only recorded trial-owned services, verifying their command lines and disappearance.
