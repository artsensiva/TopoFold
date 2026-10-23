"""Obtain the four KRAS example structures and verify them against the recorded SHA-256 values (stdlib only).

    python benchmarks/audit/structpair_demo/fetch_inputs.py                 # download from files.rcsb.org (~2.4 MB)
    python benchmarks/audit/structpair_demo/fetch_inputs.py --source DIR    # copy from a local directory instead

Files are placed where the example specs expect them (git-ignored): kras_struct_cmp/data_pdb/ and
kras_struct_g12d/data_pdb/. A file whose hash differs from the record is refused (the archive copy may have been
remediated since; the recorded results then no longer apply to it). PDB archive data are CC0 1.0.
"""

import argparse
import hashlib
import shutil
import sys
import urllib.request
from pathlib import Path

AUDIT = Path(__file__).resolve().parents[1]
FILES = {
    "9O0N.pdb": ("kras_struct_cmp/data_pdb", "b8ab08fedf0c9c330def32232e984112f6e497225d43489ebd2febcd4c50a92e"),
    "9O0O.pdb": ("kras_struct_cmp/data_pdb", "8ac5d4be7827bbc4c938cc733f09c2a1f1292b3e33bee575ae84a49fbd0e0d35"),
    "5US4.pdb": ("kras_struct_g12d/data_pdb", "a1d0c9105887a5b64deacf92d17ab5e62a2cdb1d3913f48b8203e9c84f36d22d"),
    "5USJ.pdb": ("kras_struct_g12d/data_pdb", "dab6b79531b0c046050958b914e380d8ba613b8e995522acb5ea323f974dbefe"),
}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", help="local directory holding the four .pdb files (no network access)")
    a = ap.parse_args(argv)
    bad = 0
    for name, (sub, want) in FILES.items():
        dst = AUDIT / sub / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        if not (dst.is_file() and hashlib.sha256(dst.read_bytes()).hexdigest() == want):
            if a.source:
                shutil.copyfile(Path(a.source) / name, dst)
            else:
                urllib.request.urlretrieve(f"https://files.rcsb.org/download/{name}", dst)
        got = hashlib.sha256(dst.read_bytes()).hexdigest()
        ok = got == want
        bad += not ok
        print(f"{name}: {'ok' if ok else 'HASH DIFFERS (' + got[:12] + ')'}  {dst}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
