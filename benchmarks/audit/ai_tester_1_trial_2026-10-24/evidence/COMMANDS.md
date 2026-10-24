# Executed commands

Commands below were executed, in fresh isolated trial dependencies. Most ran with cwd `/home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/checkout`; full per-command environment overrides are recorded in each JSON. Do not replay output commands into existing result directories. Viewer and evidence saver are stopped; detached service exit status is unknown, with shutdown verified separately.

### 01-clone

```sh
git clone --branch scientific-rebuild --single-branch https://github.com/artsensiva/TopoFold.git /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/checkout
```

Exit: 0; wall seconds: 34.31928296200931. Full environment, timestamps and stdout/stderr: `01-clone.json` and companion log (clone stdout is inside JSON).

### 02-checkout

```sh
git checkout --detach 6ce7b49922311e831fc09e701fcb43ae9e314625
```

Exit: 0; wall seconds: 0.1645919168367982. Full environment, timestamps and stdout/stderr: `02-checkout.json` and companion log (clone stdout is inside JSON).

### 03-initial-status

```sh
git status --porcelain=v1 --untracked-files=all
```

Exit: 0; wall seconds: 0.008328284602612257. Full environment, timestamps and stdout/stderr: `03-initial-status.json` and companion log (clone stdout is inside JSON).

### 04-venv

```sh
python3 -m venv /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/.venv
```

Exit: 0; wall seconds: 2.920644453726709. Full environment, timestamps and stdout/stderr: `04-venv.json` and companion log (clone stdout is inside JSON).

### 05-install-maturin

```sh
pip install maturin
```

Exit: 0; wall seconds: 5.729135917034. Full environment, timestamps and stdout/stderr: `05-install-maturin.json` and companion log (clone stdout is inside JSON).

### 06-maturin-build

```sh
maturin develop --release
```

Exit: 0; wall seconds: 54.6977568147704. Full environment, timestamps and stdout/stderr: `06-maturin-build.json` and companion log (clone stdout is inside JSON).

### 07-fetch-inputs

```sh
python benchmarks/audit/structpair_demo/fetch_inputs.py
```

Exit: 0; wall seconds: 3.1756560811772943. Full environment, timestamps and stdout/stderr: `07-fetch-inputs.json` and companion log (clone stdout is inside JSON).

### 08-install-computation-deps

```sh
pip install MDAnalysis==2.10.0 matplotlib==3.11.2
```

Exit: 0; wall seconds: 53.6415826450102. Full environment, timestamps and stdout/stderr: `08-install-computation-deps.json` and companion log (clone stdout is inside JSON).

### 09-documented-example

```sh
python -m topofold_inputs.structpair benchmarks/audit/structpair_demo/specs/kras2_5US4A_vs_5USJA.json --out /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/results/kras2 --adapter-python /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/.venv/bin/python
```

Exit: 0; wall seconds: 7.0452040187083185. Full environment, timestamps and stdout/stderr: `09-documented-example.json` and companion log (clone stdout is inside JSON).

### 10-additional-case

```sh
python -m topofold_inputs.structpair /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/specs/explicit_55_72.json --out /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/results/explicit_55_72 --adapter-python /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/.venv/bin/python
```

Exit: 0; wall seconds: 5.785191902890801. Full environment, timestamps and stdout/stderr: `10-additional-case.json` and companion log (clone stdout is inside JSON).

### 11-install-viewer-deps

```sh
pip install streamlit==1.64.0 pandas==3.0.6 altair==6.3.0 plotly==7.1.0
```

Exit: 0; wall seconds: 48.27822992997244. Full environment, timestamps and stdout/stderr: `11-install-viewer-deps.json` and companion log (clone stdout is inside JSON).

### 12-viewer

```sh
/home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/.venv/bin/streamlit run apps/structpair_viewer/app.py --server.address localhost -- --results /home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/results
```

Exit: None; wall seconds: unknown. Full environment, timestamps and stdout/stderr: `12-viewer.json` and companion log (clone stdout is inside JSON).

### 13-dependency-freeze

```sh
pip freeze
```

Exit: 0; wall seconds: 0.41524617467075586. Full environment, timestamps and stdout/stderr: `13-dependency-freeze.json` and companion log (clone stdout is inside JSON).

Independent validation command and exit0: `numerical-validation-command.json`. Final Git commands and exit0: `15-final-tracked-audit.json`. Browser operations and tester-only failures: `17-browser-observations.json`, `18-tester-tooling-events.json`.
