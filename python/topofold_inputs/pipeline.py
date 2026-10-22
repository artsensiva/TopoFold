"""One command: local input -> checked Cα coordinates -> κ/τ features -> JSON/NPZ + inspection report.

    PYTHONPATH=python .venv/bin/python -m topofold_inputs.pipeline TOPOLOGY [TRAJECTORY ...] --run-id ID --out DIR \\
        [--adapter-python PATH] [adapter options]

Stage 1 (adapter) needs MDAnalysis; stage 2 (features) needs the ``topofold`` extension. If the running interpreter
lacks MDAnalysis, stage 1 runs in ``--adapter-python`` (or $TOPOFOLD_ADAPTER_PYTHON) as a subprocess. The stages
exchange only the hashed adapter export files.
"""

import argparse
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

PKG_ROOT = Path(__file__).resolve().parents[1]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m topofold_inputs.pipeline", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="topology file, then trajectory files")
    ap.add_argument("--run-id", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--adapter-python", default=os.environ.get("TOPOFOLD_ADAPTER_PYTHON"))
    ap.add_argument("--min-segment", type=float, default=1e-6, help="degenerate segment threshold, Å")
    ap.add_argument("--collinear-sin-tol", type=float, default=1e-6, help="|sin κ| below this is near-collinear")
    args, adapter_args = ap.parse_known_args(argv)

    cmd = [*args.inputs, "--run-id", args.run_id, "--out", args.out, *adapter_args]
    if importlib.util.find_spec("MDAnalysis") is not None:
        from .__main__ import main as adapter_main
        rc = adapter_main(cmd)
    else:
        if not args.adapter_python:
            print("MDAnalysis is not importable here; pass --adapter-python (an interpreter with MDAnalysis)",
                  file=sys.stderr)
            return 2
        env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(PKG_ROOT), os.environ.get("PYTHONPATH", "")]))
        rc = subprocess.run([args.adapter_python, "-m", "topofold_inputs", *cmd], env=env).returncode
    if rc != 0:
        print(f"adapter stage failed (exit {rc}); no features computed", file=sys.stderr)
        return rc

    from .features import export_features, features_from_export
    out = Path(args.out)
    a_npz, a_json = out / f"{args.run_id}.npz", out / f"{args.run_id}.json"
    feats, _ = features_from_export(a_npz, a_json, min_segment_angstrom=args.min_segment,
                                    collinear_sin_tol=args.collinear_sin_tol)
    npz, js, md = export_features(feats, a_npz, a_json, out, args.run_id)
    k, t = feats.kappa_valid, feats.tau_valid
    print(f"features: kappa valid {int(k.sum())}/{k.size}, tau valid {int(t.sum())}/{t.size}")
    print(f"wrote {npz}\nwrote {js}\nwrote {md}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
