"""Post-run diagnostic for the ABL-MYR-1 G1 stop: which required archive members exist, from zip directory listings only.

    PYTHONPATH=python .venv/bin/python benchmarks/audit/abl_myr_1/g1_member_audit.py

Not part of the frozen protocol and not a re-run: it extracts and opens no member, applies no threshold and computes no
scientific quantity. It reuses the frozen runner's member list and matching rule, and reads the dataset README's file tree
(a documentation file) to say whether a missing member is missing in the published dataset itself.
"""

import json
import re
import sys
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(HERE))
import run_abl_myr as R  # noqa: E402


def audit_archives(data_dir: Path, cfg: dict) -> dict:
    out = {}
    for arch, keys in R.wanted_members(cfg).items():
        names = zipfile.ZipFile(Path(data_dir) / arch).namelist()
        missing = []
        for k, suffix in keys.items():
            try:
                R.resolve_member(names, suffix)
            except R.GateError as e:
                missing.append({"key": k, "suffix": suffix, "matches": len(e.detail["hits"])})
        out[arch] = {"entries": len(names), "required": len(keys), "present": len(keys) - len(missing), "missing": missing}
    return out


def readme_set_gro_files(text: str) -> dict:
    """{model: {set: [.gro file names]}} from the dataset README's directory tree."""
    out, model, cur = {}, None, None
    for ln in text.splitlines():
        m = re.match(r"^\+---(kinase_\w+|SH3-SH2-KD_\w+)\s*$", ln)
        if m:
            model, cur = m.group(1), None
            out.setdefault(model, {})
            continue
        m = re.search(r"[+\\]---(Set\d)\s*$", ln)
        if m and model:
            cur = m.group(1)
            out[model].setdefault(cur, [])
            continue
        m = re.search(r"\|\s+([\w.\-]+\.gro)\s*$", ln)
        if m and model and cur:
            out[model][cur].append(m.group(1))
    return out


def main() -> int:
    spec = json.loads(R.SPEC.read_text())
    data = HERE / "data_dryad"
    readme = data / "README.md"
    intake = json.loads((ROOT / "benchmarks/audit/debuhr_abl_intake/intake_record.json").read_text())
    dryad = next(f["sha256_source"] for f in intake["raw_inventory"]["files"] if f["path"] == "README.md")
    tree = readme_set_gro_files(readme.read_text(encoding="utf-8", errors="replace"))
    rec = {"purpose": "diagnostic after the original ABL-MYR-1 run stopped at G1; not part of the frozen protocol",
           "spec_sha256": R.sha256_file(R.SPEC),
           "archives_listing_only": audit_archives(data, R.CFG),
           "dataset_readme": {"sha256": R.sha256_file(readme), "dryad_listed_sha256": dryad,
                              "matches_dryad_listing": R.sha256_file(readme) == dryad,
                              "set_gro_files_in_tree": {m: tree[m] for m in ("kinase_kinked_Myr", "kinase_kinked_noMyr")}}}
    (HERE / "results/ABL-MYR-1_g1_member_audit.json").write_text(json.dumps(rec, indent=1) + "\n")
    print(json.dumps({a: (v["present"], v["required"]) for a, v in rec["archives_listing_only"].items()}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
