"""Metadata-only inventory of the Folding@home COVID Moonshot Mpro dynamics data (AWS Open Data, CC0).

    python3 benchmarks/audit/mpro_dataset_search/list_moonshot_metadata.py

Lists S3 object keys and sizes (no trajectory is downloaded). It fetches two small metadata files (saved as fah-README.md and
fah-metadata.csv, about 45 kB) into this directory, and writes listing_summary.json with per-project RUN/CLONE counts and
sizes. Projects (from the dataset README): 13430/13431 apo monomer, 13432/13433 holo monomer, 13434/13435 apo dimer,
13436/13437 holo dimer; even = His41(0) Cys145(0), odd = His41(+) Cys145(-).
"""

import collections
import datetime as dt
import hashlib
import json
import re
import statistics
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUCKET = "https://fah-public-data-covid19-moonshot-dynamics.s3.us-east-2.amazonaws.com"
PREFIX = "SARS-CoV-2-Mpro-Diamond-structures/"
PROJECTS = [str(p) for p in range(13430, 13438)]
KEY_RE = re.compile(r"<Key>([^<]*)</Key><LastModified>[^<]*</LastModified><ETag>[^<]*</ETag><Size>(\d+)</Size>")


def list_keys(prefix):
    tok, out = None, []
    while True:
        q = {"list-type": "2", "max-keys": "1000", "prefix": prefix}
        if tok:
            q["continuation-token"] = tok
        t = urllib.request.urlopen(f"{BUCKET}/?{urllib.parse.urlencode(q)}", timeout=60).read().decode()
        out += [(k, int(s)) for k, s in KEY_RE.findall(t)]
        m = re.search(r"<NextContinuationToken>([^<]*)", t)
        if not m:
            return out
        tok = m.group(1)


def main():
    meta = {}
    for name, local in (("README.md", "fah-README.md"), ("fah-metadata.csv", "fah-metadata.csv")):
        data = urllib.request.urlopen(f"{BUCKET}/{PREFIX}setup/{name}", timeout=60).read()
        (HERE / local).write_bytes(data)
        meta[local] = {"source_key": f"{PREFIX}setup/{name}", "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
    projects = {}
    for p in PROJECTS:
        per = collections.defaultdict(list)
        for k, s in list_keys(f"{PREFIX}solute-hdf5-trajectories/{p}/"):
            m = re.search(r"run(\d+)-clone(\d+)\.h5$", k)
            if m:
                per[int(m.group(1))].append(s)
        sizes = [s for v in per.values() for s in v]
        projects[p] = {"runs_with_files": len(per), "files": len(sizes), "total_GB": round(sum(sizes) / 1e9, 2),
                       "clones_per_run_distribution": dict(sorted(collections.Counter(len(v) for v in per.values()).items())),
                       "file_MB_median": round(statistics.median(sizes) / 1e6, 1) if sizes else None,
                       "file_MB_min": round(min(sizes) / 1e6, 2) if sizes else None,
                       "file_MB_max": round(max(sizes) / 1e6, 1) if sizes else None}
    out = {"id": "MPRO-DATASET-SEARCH-1", "listed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "bucket": BUCKET, "prefix": PREFIX, "license": "CC0 (AWS Open Data registry entry 'foldingathome-covid19')",
           "metadata_files": meta, "projects": projects,
           "statement": "Object listing only; no trajectory downloaded. Trajectory lengths are not verified from file contents."}
    (HERE / "listing_summary.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
