# Executed command ledger

Working directory for logged commands: /home/artem/Documents/Codex/2026-10-25/files-pasted-by-the-user-you/outputs/topofold-independent-trial/clone

### environment — self-devised
UTC 2026-10-25T09:56:16.412155+00:00; exit 0; wall 0.380 s.
```bash
uname -a; cat /etc/os-release; python3 --version; rustc --version; cargo --version; git --version; python3 -m pip list --format=json
```

### download-abl — self-devised
UTC 2026-10-25T09:56:40.267208+00:00; exit 0; wall 2.966 s.
```bash
for id in 1OPJ 1OPK; do curl -fL --connect-timeout 15 --max-time 60 https://files.rcsb.org/download/$id.pdb -o ../inputs/$id.pdb || exit; curl -fL --connect-timeout 15 --max-time 60 https://data.rcsb.org/rest/v1/core/entry/$id -o ../inputs/$id-entry.json || exit; done; sha256sum ../inputs/*
```

### audit-inputs — self-devised
UTC 2026-10-25T09:57:08.656817+00:00; exit 0; wall 0.048 s.
```bash
python ../audit_inputs.py
```

### install-core — guide section 2; environment and build paths relocated
UTC 2026-10-25T09:56:36.848038+00:00; exit 0; wall 75.307 s.
```bash
python3 -m venv ../env && source ../env/bin/activate && pip install maturin && maturin develop --release
```

### fetch-help — self-devised discovery; unchanged guide fetcher copied outside clone
UTC 2026-10-25T09:58:01.598319+00:00; exit 0; wall 0.091 s.
```bash
python ../smoke-assets/fetch_inputs.py --help
```

### fetch-kras — guide downloader; unchanged copy relocated outside clone
UTC 2026-10-25T09:58:20.407371+00:00; exit 0; wall 3.576 s.
```bash
python ../smoke-assets/fetch_inputs.py
```

### cli-help — self-devised help discovery
UTC 2026-10-25T09:58:42.159519+00:00; exit 0; wall 0.224 s.
```bash
PYTHONPATH=python python -m topofold_inputs.structpair --help
```

### install-computation — guide section 2 step 4 exact package pins
UTC 2026-10-25T09:58:17.179216+00:00; exit 0; wall 74.596 s.
```bash
pip install MDAnalysis==2.10.0 matplotlib==3.11.2
```

### kras-run — guide command; input spec and output paths relocated
UTC 2026-10-25T09:59:39.094524+00:00; exit 0; wall 7.531 s.
```bash
PYTHONPATH=python /usr/bin/time -v python -m topofold_inputs.structpair ../smoke-assets/specs/kras2_5US4A_vs_5USJA.json --out ../results/kras2
```

### build-record — self-devised
UTC 2026-10-25T09:59:49.802823+00:00; exit 0; wall 0.394 s.
```bash
git rev-parse HEAD; git status --short --ignored --untracked-files=all; sha256sum Cargo.lock docs/STRUCTPAIR_USER_GUIDE.md ../smoke-assets/fetch_inputs.py benchmarks/audit/structpair_demo/fetch_inputs.py ../smoke-assets/specs/kras2_5US4A_vs_5USJA.json benchmarks/audit/structpair_demo/specs/kras2_5US4A_vs_5USJA.json; pip freeze
```

### abl-run — new frozen spec following guide sections 3 and 3a
UTC 2026-10-25T10:00:02.510108+00:00; exit 0; wall 8.027 s.
```bash
PYTHONPATH=python /usr/bin/time -v python -m topofold_inputs.structpair ../specs/abl.json --out ../results/abl
```

### verify-manifests — self-devised
UTC 2026-10-25T10:00:33.598222+00:00; exit 0; wall 0.059 s.
```bash
python ../verify_manifests.py
```

### install-viewer — guide section 2 step 5 exact pins
UTC 2026-10-25T10:00:30.740236+00:00; exit 0; wall 52.664 s.
```bash
pip install streamlit==1.64.0 pandas==3.0.6 altair==6.3.0 plotly==7.1.0
```

### interpretation — self-devised
UTC 2026-10-25T10:02:21.068775+00:00; exit 0; wall 0.065 s.
```bash
python ../interpret_results.py
```

### verify-after-interruption — self-devised
UTC 2026-10-25T16:18:03.682202+00:00; exit 0; wall 0.063 s.
```bash
python ../verify_manifests.py
```

### final-verification — self-devised
UTC 2026-10-25T16:28:34.213383+00:00; exit 0; wall 0.446 s.
```bash
python ../verify_manifests.py && git rev-parse HEAD && git status --short --ignored --untracked-files=all && pip freeze
```

### Viewer (initial launch and restart)
Same documented command, both from the clone and isolated environment. See `viewer.log` and `viewer-resume.log` for exact start times. The interrupted original process did not record a normal exit.
```bash
streamlit run apps/structpair_viewer/app.py --server.address localhost -- --results ../results
```

Trial-only scripts, browser interactions and the corruption-test attempts are described in the report and tamper-test.json. No scientific command was rerun.