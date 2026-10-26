"""KIN-PAIR-1 pre-protocol bookkeeping (metadata only; no coordinates, no TopoFold score).

    python3 benchmarks/audit/kin_pair_1/prep_protocol_inputs.py

1. Verifies the declared candidate kinase accessions against UniProt (gene, organism, name, length).
2. Collects every PDB-ID-like token in the repository's development records, keeps those that RCSB confirms, and records the
   UniProt accessions of their polymer entities: the exposure history that the protocol excludes.
Writes results/candidate_targets_verified.json and results/exposure_history.json. Standard library only.
"""

import json
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CANDIDATES = {   # gene -> UniProt accession, declared before any RCSB structure query (oncology-relevant human protein kinases)
    "EGFR": "P00533", "ERBB2": "P04626", "MET": "P08581", "ALK": "Q9UM73", "KIT": "P10721", "FLT3": "P36888",
    "FGFR1": "P11362", "KDR": "P35968", "RET": "P07949", "SRC": "P12931", "LCK": "P06239", "JAK2": "O60674",
    "BTK": "Q06187", "MAP2K1": "Q02750", "CHEK1": "O14757", "PLK1": "P53350", "AKT1": "P31749", "CDK6": "Q00534",
    "MAPK1": "P28482", "SYK": "P43405", "IRAK4": "Q9NWZ3", "PIM1": "P11309", "PDPK1": "O15530"}
SCAN = ["docs/*.md", "docs/pilot/*.md", "benchmarks/audit/*/README.md", "benchmarks/audit/*/run_spec.json",
        "benchmarks/audit/*/PROTOCOL*.md", "benchmarks/audit/*/results/*.json", "benchmarks/README.md", "README.md"]
SKIP_PARTS = ("REPORT-ChatGPT", "Научно")


def get(url):
    for k in range(4):
        try:
            return json.load(urllib.request.urlopen(url, timeout=30))
        except Exception:
            if k == 3:
                return None


def main():
    tsv = urllib.request.urlopen("https://rest.uniprot.org/uniprotkb/search?query=" + "+OR+".join(
        f"accession:{a}" for a in CANDIDATES.values()) + "&fields=accession,gene_primary,organism_id,length,protein_name&format=tsv&size=50",
        timeout=30).read().decode().splitlines()[1:]
    got = {r.split("\t")[0]: r.split("\t") for r in tsv}
    ver = {}
    for gene, acc in CANDIDATES.items():
        r = got.get(acc)
        ver[gene] = {"accession": acc, "uniprot_gene": r[1] if r else None, "organism_id": r[2] if r else None,
                     "length": int(r[3]) if r else None, "name": r[4] if r else None,
                     "verified": bool(r and r[1].split()[0] == gene and r[2] == "9606")}
    (HERE / "results/candidate_targets_verified.json").write_text(json.dumps(ver, indent=1) + "\n")
    toks = set()
    for pat in SCAN:
        for p in sorted(ROOT.glob(pat)):
            if any(s in str(p) for s in SKIP_PARTS) or not p.is_file():
                continue
            toks |= set(re.findall(r"(?<![A-Za-z0-9_])([1-9][A-Za-z0-9]{3})(?![A-Za-z0-9_])", p.read_text(errors="ignore").upper()))
    tracked = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True).stdout.split()
    for f in tracked:
        m = re.fullmatch(r"(?:.*/)?([1-9][A-Za-z0-9]{3})\.(?:pdb|cif)", f)
        if m:
            toks.add(m.group(1).upper())
    toks = sorted(t for t in toks if re.search(r"[A-Z]", t) and re.search(r"[0-9]", t))
    exposed, rejected = {}, []
    for t in toks:
        e = get(f"https://data.rcsb.org/rest/v1/core/entry/{t}")
        if not e:
            rejected.append(t)
            continue
        accs = set()
        for i in e.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", []):
            pe = get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{t}/{i}") or {}
            for r in pe.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers", []) or []:
                if r.get("database_name") == "UniProt":
                    accs.add(r["database_accession"])
        exposed[t] = sorted(accs)
    (HERE / "results/exposure_history.json").write_text(json.dumps(
        {"method": "PDB-ID-like tokens in the listed repository records and tracked .pdb/.cif file names; kept if RCSB confirms the entry",
         "scan_globs": SCAN, "n_tokens": len(toks), "n_rejected_not_in_rcsb": len(rejected), "entries": exposed,
         "exposed_uniprot_accessions": sorted({a for v in exposed.values() for a in v})}, indent=1) + "\n")
    print(sum(v["verified"] for v in ver.values()), "of", len(ver), "candidate accessions verified;", len(exposed), "exposed entries;",
          len({a for v in exposed.values() for a in v}), "exposed accessions")


if __name__ == "__main__":
    sys.exit(main())
