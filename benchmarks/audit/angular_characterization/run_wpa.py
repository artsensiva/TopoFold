#!/usr/bin/env python3
"""Run WP-A end to end and record provenance (SYNTHETIC fixtures only).

Steps: record sources -> generate fixtures -> Rust example (R1/R2) -> build wheel from the
current source -> install it into the isolated venv -> verify the loaded extension is the one
just built -> P1/REF -> compare -> report. Every command and exit code is recorded.

Usage (from the repository root or anywhere):
  .venv/bin/python benchmarks/audit/angular_characterization/run_wpa.py
The isolated venv (build/venv, with pinned numpy and pytest) must exist; see REPORT.md.
"""
import datetime
import subprocess
import sys
import zipfile
from pathlib import Path

from common import (BUILD, FIXTURES_TSV, FIXTURE_BIN, HERE, PROVENANCE_JSON, REPO, RUST_JSON,
                    dump_json, sha256_file)

VENV_PY = BUILD / 'venv' / 'bin' / 'python'
WHEELS = BUILD / 'wheels'
MATURIN = REPO / '.venv' / 'bin' / 'maturin'
SOURCES = ['crates/topofold-core/src/circular.rs', 'crates/topofold-core/src/bimodality.rs',
           'crates/topofold-core/src/geometry.rs', 'crates/topofold-core/src/types.rs',
           'crates/topofold-core/src/lib.rs', 'crates/topofold-python/src/lib.rs',
           'Cargo.lock', 'crates/topofold-core/examples/angular_characterization.rs',
           'docs/ANGULAR_REPRESENTATION_CHARACTERIZATION_PROTOCOL.md']
HARNESS = ['common.py', 'build_traces.py', 'reference.py', 'generate_fixtures.py',
           'run_python_path.py', 'compare.py', 'make_report.py', 'run_wpa.py', 'test_harness.py']

log = []


def run(cmd, cwd=REPO, check=True):
    p = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    log.append({'cmd': ' '.join(str(c) for c in cmd), 'cwd': str(Path(cwd).relative_to(REPO))
                if Path(cwd).is_relative_to(REPO) else str(cwd), 'exit': p.returncode,
                'stdout_tail': p.stdout.strip().splitlines()[-3:]})
    print(f'[exit {p.returncode}] {log[-1]["cmd"]}')
    if check and p.returncode != 0:
        print(p.stdout[-3000:], p.stderr[-3000:])
        finish(status='failed')
        sys.exit(1)
    return p


def out(cmd, cwd=REPO):
    return run(cmd, cwd).stdout.strip()


prov = {}


def finish(status):
    prov['status'] = status
    prov['commands'] = log
    dump_json(prov, PROVENANCE_JSON)


def main():
    prov['_meta'] = {'type': 'SYNTHETIC', 'producer': 'run_wpa.py',
                     'started_utc': datetime.datetime.now(datetime.UTC).isoformat()}
    prov['git'] = {'head': out(['git', 'rev-parse', 'HEAD']),
                   'status_porcelain': out(['git', 'status', '--porcelain']).splitlines()}
    prov['source_sha256'] = {s: sha256_file(REPO / s) for s in SOURCES}
    prov['harness_sha256'] = {h: sha256_file(HERE / h) for h in HARNESS if (HERE / h).exists()}
    prov['versions'] = {'rustc': out(['rustc', '--version']), 'cargo': out(['cargo', '--version']),
                        'maturin': out([MATURIN, '--version']),
                        'venv_python': out([VENV_PY, '-c', 'import sys;print(sys.version)']),
                        'venv_numpy': out([VENV_PY, '-c', 'import numpy;print(numpy.__version__)']),
                        'venv_pytest': out([VENV_PY, '-c', 'import pytest;print(pytest.__version__)'])}

    run([VENV_PY, '-B', HERE / 'generate_fixtures.py'], cwd=HERE)
    prov['fixtures_json_sha256'] = sha256_file(HERE / 'fixtures' / 'SYNTHETIC_fixtures.json')

    rust_arrays = (BUILD / 'arrays' / 'SYNTHETIC_rust_arrays.f64').relative_to(REPO)
    (BUILD / 'arrays').mkdir(parents=True, exist_ok=True)
    run(['cargo', 'run', '--release', '-q', '-p', 'topofold-core', '--example',
         'angular_characterization', '--', FIXTURES_TSV.relative_to(REPO),
         FIXTURE_BIN.relative_to(REPO), RUST_JSON.relative_to(REPO), rust_arrays])

    # Build the binding from the current source and install it into the isolated venv.
    for w in WHEELS.glob('*.whl'):
        w.unlink()
    run([MATURIN, 'build', '--release', '-i', VENV_PY, '--out', WHEELS])
    wheels = sorted(WHEELS.glob('*.whl'))
    assert len(wheels) == 1, wheels
    wheel = wheels[0]
    with zipfile.ZipFile(wheel) as z:
        so_names = [n for n in z.namelist() if n.endswith('.so')]
        import hashlib
        wheel_so_sha = {n: hashlib.sha256(z.read(n)).hexdigest() for n in so_names}
    run([VENV_PY, '-m', 'pip', 'install', '-q', '--force-reinstall', '--no-deps', wheel])
    ident = out([VENV_PY, '-B', '-c',
                 'import topofold, topofold.topofold as e, hashlib, pathlib;'
                 'p=pathlib.Path(e.__file__).resolve();'
                 'print(p); print(hashlib.sha256(p.read_bytes()).hexdigest())'], cwd=HERE)
    so_path, so_sha = ident.splitlines()
    prov['extension'] = {'wheel': wheel.name, 'wheel_sha256': sha256_file(wheel),
                         'wheel_so_sha256': wheel_so_sha, 'loaded_extension_path': so_path,
                         'loaded_extension_sha256': so_sha,
                         'loaded_equals_freshly_built': so_sha in wheel_so_sha.values(),
                         'loaded_from_isolated_venv': so_path.startswith(str(BUILD / 'venv'))}
    if not (prov['extension']['loaded_equals_freshly_built']
            and prov['extension']['loaded_from_isolated_venv']):
        finish('failed: loaded extension is not the freshly built one')
        sys.exit(1)

    run([VENV_PY, '-B', HERE / 'run_python_path.py'], cwd=HERE)
    cmp_ = run([VENV_PY, '-B', HERE / 'compare.py'], cwd=HERE, check=False)
    prov['compare_exit'] = cmp_.returncode
    finish('completed' if cmp_.returncode == 0 else 'completed: harness checks failed')
    rep = run([VENV_PY, '-B', HERE / 'make_report.py'], cwd=HERE, check=False)
    prov['report_exit'] = rep.returncode
    finish(prov['status'])
    return cmp_.returncode


if __name__ == '__main__':
    sys.exit(main())
