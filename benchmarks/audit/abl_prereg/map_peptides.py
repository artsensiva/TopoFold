"""Map the HDX-MS peptide list of PRIDE PXD034008 (construct numbering) to c-Abl isoform 1a numbering.

    python3 benchmarks/audit/abl_prereg/map_peptides.py

Inputs:
- the peptide map `REI_ABL_KB_dynamX_export_18May22.csv` (CC0, git-ignored; SHA-256 checked);
- UniProt P00519 canonical sequence (`reference/UniProt_P00519.fasta`, CC BY 4.0).

Each peptide sequence is located in P00519 by exact matching. The offset (c-Abl = construct + offset) is accepted only if
every uniquely matching peptide gives the same offset. Peptides that do not match exactly are compared at that offset
and their differences are listed (they are expected to be the mutant peptides). Writes
`hdx_peptide_mapping.json`. This uses sequences only: no uptake values exist in the input, and none are inferred.
"""

import csv
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PEPTIDES = ROOT / "benchmarks/audit/abl_access/data_pride_PXD034008/REI_ABL_KB_dynamX_export_18May22.csv"
PEPTIDES_SHA = "2fd6b624a52fd5261c838f85a9c637f29b94e098297beb37f58c76b71705ca4a"
FASTA = HERE / "reference/UniProt_P00519.fasta"
OUT = HERE / "hdx_peptide_mapping.json"


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main() -> int:
    if not PEPTIDES.is_file():
        print(f"missing input: {PEPTIDES} (git-ignored; see benchmarks/audit/abl_access/ACCESS_STATUS_2026-10-23.md)",
              file=sys.stderr)
        return 2
    if sha256(PEPTIDES) != PEPTIDES_SHA:
        print("peptide map SHA-256 differs from the recorded value", file=sys.stderr)
        return 2
    seq = "".join(l.strip() for l in FASTA.read_text().splitlines() if not l.startswith(">"))
    rows = list(csv.reader(PEPTIDES.open()))
    header_i = next(i for i, r in enumerate(rows) if r and r[0] == "Sequence")
    peps = [r for r in rows[header_i + 1:] if r and r[0]]
    offsets, unmatched = set(), []
    for s, st, en, *_ in peps:
        hits = [i for i in range(len(seq)) if seq.startswith(s, i)]
        if len(hits) == 1:
            offsets.add(hits[0] + 1 - int(st))
        else:
            unmatched.append((s, int(st), int(en), len(hits)))
    if len(offsets) != 1:
        print(f"no single consistent offset: {sorted(offsets)}", file=sys.stderr)
        return 1
    off = offsets.pop()
    out_peps = []
    for s, st, en, mhp, maxu, rt in peps:
        st, en = int(st), int(en)
        ref = seq[st + off - 1: en + off]
        diffs = [{"cAbl_residue": st + off + k, "reference": ref[k], "peptide": s[k]}
                 for k in range(len(s)) if k < len(ref) and ref[k] != s[k]]
        out_peps.append({"sequence": s, "construct_start": st, "construct_end": en,
                         "cAbl_start": st + off, "cAbl_end": en + off, "max_uptake": int(maxu),
                         "exact_match_to_P00519": not diffs and len(ref) == len(s), "differences": diffs})
    cover = sorted({r for p in out_peps for r in range(p["cAbl_start"], p["cAbl_end"] + 1)})
    res = {"id": "ABL-HDX-PEPTIDE-MAP-1", "inputs": {"peptide_map": {"path": str(PEPTIDES.relative_to(ROOT)), "sha256": PEPTIDES_SHA,
                                                                    "source": "PRIDE PXD034008 (CC0)"},
                                                     "uniprot": {"path": str(FASTA.relative_to(ROOT)), "sha256": sha256(FASTA),
                                                                 "accession": "P00519 canonical (isoform 1a)", "license": "CC BY 4.0"}},
           "offset_construct_to_cAbl_1a": off, "n_peptides": len(out_peps),
           "n_exact_unique_matches": len(out_peps) - len(unmatched),
           "non_matching_peptides": [{"sequence": p["sequence"], "cAbl_start": p["cAbl_start"], "differences": p["differences"]}
                                     for p in out_peps if not p["exact_match_to_P00519"]],
           "coverage_cAbl": {"first": cover[0], "last": cover[-1], "n_residues_covered": len(cover),
                             "span_length": cover[-1] - cover[0] + 1},
           "landmarks_in_P00519": {str(p): seq[p - 1] for p in (271, 286, 315, 344, 381, 382, 383, 472)},
           "peptides": out_peps,
           "statement": "Sequence-based numbering map only. The input contains no deuterium-uptake values; none are inferred."}
    OUT.write_text(json.dumps(res, indent=1) + "\n")
    print(json.dumps({k: res[k] for k in ("offset_construct_to_cAbl_1a", "n_peptides", "n_exact_unique_matches",
                                           "non_matching_peptides", "coverage_cAbl")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
