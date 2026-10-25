"""Data inventory of the small public Abl HDX-MS files (2026-10-25). Not an analysis and not a study.

    python3 -I benchmarks/audit/abl_hdx_lead_2026-10-25/inventory_hdx_lead.py --fetch-dir <new empty directory outside the repo>

Inspects only: the 3,386-byte DynamX export of PRIDE PXD034008 (provider record and SHA-1 verified), the article's directly
linked Source Data workbook (structure only: sheet names, dimensions, header strings) and the Supplementary Information PDF
(text layer only). No raw mass-spectrometry file and no trajectory is fetched. Reports which fields are present in each file,
which are absent, and which are only stated in the paper. Writes results/hdx_inventory.json. Needs `pdftotext` (poppler).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import html
import io
import json
import re
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / "results/hdx_inventory.json"
REF = ROOT / "benchmarks/audit/abl_prereg/reference/UniProt_P00519.fasta"
PRIDE = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD034008"
CSV_NAME = "REI_ABL_KB_dynamX_export_18May22.csv"
ESM = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-36956-5/MediaObjects/41467_2023_36956_MOESM{n}_ESM.{ext}"
MONO = {"G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841, "T": 101.04768, "C": 103.00919, "L": 113.08406,
        "I": 113.08406, "N": 114.04293, "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259, "M": 131.04049,
        "H": 137.05891, "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931}
H2O, PROTON = 18.010565, 1.007276          # monoisotopic water mass and proton mass (constants)


def get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "TopoFold-hdx-inventory/1"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()


def digests(b: bytes) -> dict:
    return {"bytes": len(b), "sha1": hashlib.sha1(b).hexdigest(), "sha256": hashlib.sha256(b).hexdigest()}


def read_fasta(p: Path) -> str:
    return "".join(ln.strip() for ln in p.read_text().splitlines() if not ln.startswith(">"))


def parse_export(text: str) -> dict:
    """Rows of a DynamX peptide export: the preamble rows before the header are kept as given."""
    rows = list(csv.reader(io.StringIO(text, newline="")))
    h = next(i for i, r in enumerate(rows) if r and r[0] == "Sequence")
    return {"preamble": rows[:h], "header": rows[h], "rows": [r for r in rows[h + 1:] if any(r)]}


def file_shape(raw: bytes) -> dict:
    """Physical shape of the export: lines, line endings, trailing newline, the preamble row and the first and last rows."""
    text = raw.decode("ascii")
    lines = text.split("\n")
    trailing = lines[-1] == ""
    body = lines[:-1] if trailing else lines
    crlf = all(ln.endswith("\r") for ln in body)
    rows = [ln.rstrip("\r") for ln in body]
    return {"n_lines": len(body), "line_ending": "CRLF" if crlf else "mixed or LF", "trailing_newline": trailing, "first_line": rows[0],
            "header_line": rows[1], "first_data_row": rows[2], "last_data_row": rows[-1]}


def mhp(seq: str) -> float:
    """Monoisotopic [M+H]+ of a peptide."""
    return sum(MONO[c] for c in seq) + H2O + PROTON


def locate(seq: str, ref: str, start: int, variants: dict) -> dict:
    """Offset (reference index minus file Start) of the exact sequence match, in the wild-type reference or in a variant."""
    for name, r in (("wild type", ref), *variants.items()):
        i = r.find(seq)
        if i >= 0:
            return {"matched_in": name, "offset": i + 1 - start}
    return {"matched_in": None, "offset": None}


def describe_export(parsed: dict, ref: str) -> dict:
    rows = parsed["rows"]
    n = len(rows)
    seqs = [r[0] for r in rows]
    starts, ends = [int(r[1]) for r in rows], [int(r[2]) for r in rows]
    variants = {}
    for name, (pos, wt, mut) in {"A344V": (344, "A", "V"), "M472I": (472, "M", "I")}.items():
        if ref[pos - 1] == wt:
            variants[name] = ref[:pos - 1] + mut + ref[pos:]
    loc = [locate(s, ref, st, variants) for s, st in zip(seqs, starts)]
    covered = sorted({i for st, en in zip(starts, ends) for i in range(st, en + 1)})
    mh_err = max(abs(mhp(s) - float(r[3])) for s, r in zip(seqs, rows))
    return {
        "header": parsed["header"], "preamble_rows": [",".join(r) for r in parsed["preamble"]], "n_peptides": n,
        "unique_sequence_range": len({(r[0], r[1], r[2]) for r in rows}), "unique_sequences": len(set(seqs)),
        "file_numbering": {"first_start": min(starts), "last_end": max(ends), "covered_residues": len(covered),
                           "span_residues": max(ends) - min(starts) + 1},
        "length_equals_end_minus_start_plus_1": sum(len(s) == e - st + 1 for s, st, e in zip(seqs, starts, ends)),
        "reference_check": {"reference": "UniProt P00519 (benchmarks/audit/abl_prereg/reference/UniProt_P00519.fasta)",
                            "residue_344_in_reference": ref[343], "residue_472_in_reference": ref[471],
                            "matched_in_wild_type": sum(x["matched_in"] == "wild type" for x in loc),
                            "matched_only_in_a_variant": {k: sum(x["matched_in"] == k for x in loc) for k in variants},
                            "unmatched": sum(x["matched_in"] is None for x in loc),
                            "offsets_seen": sorted({x["offset"] for x in loc if x["offset"] is not None}),
                            "variant_peptides": [(s, st, e, x["matched_in"]) for s, st, e, x in zip(seqs, starts, ends, loc) if x["matched_in"] not in ("wild type",)]},
        "mhp_is_monoisotopic_MH_plus": {"max_abs_difference_Da": mh_err, "rows_within_0.001": sum(abs(mhp(s) - float(r[3])) < 1e-3 for s, r in zip(seqs, rows))},
        "maxuptake_equals_len_minus_1_minus_prolines": sum(int(r[4]) == len(r[0]) - 1 - r[0].count("P") for r in rows),
        "RT_min_max": [min(float(r[5]) for r in rows), max(float(r[5]) for r in rows)],
    }


def xlsx_inventory(data: bytes) -> dict:
    z = zipfile.ZipFile(io.BytesIO(data))
    wb = z.read("xl/workbook.xml").decode()
    rel = z.read("xl/_rels/workbook.xml.rels").decode()
    rid = {}
    for m in re.finditer(r"<Relationship ([^>]*)/?>", rel):
        a = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        rid[a["Id"]] = a["Target"]
    ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        ss = [html.unescape(re.sub(r"<[^>]+>", "", m)) for m in re.findall(r"<si>(.*?)</si>", z.read("xl/sharedStrings.xml").decode(), flags=re.S)]
    sheets = []
    for m in re.finditer(r"<sheet ([^>]*)/?>", wb):
        a = dict(re.findall(r'([\w:]+)="([^"]*)"', m.group(1)))
        x = z.read("xl/" + rid[a["r:id"]].replace("/xl/", "").lstrip("/")).decode()
        dim = re.search(r'<dimension ref="([^"]+)"', x)
        strings = []
        for row in re.findall(r"<row [^>]*>(.*?)</row>", x, flags=re.S)[:3]:
            for attrs, v in re.findall(r"<c ([^>]*)>(?:<f>.*?</f>)?<v>(.*?)</v>", row):
                if 't="s"' in attrs and int(v) < len(ss):
                    strings.append(ss[int(v)])
        sheets.append({"name": a["name"], "dimension": dim.group(1) if dim else None, "header_strings": strings[:8]})
    return {"n_sheets": len(sheets), "sheets": sheets}


def hdx_terms(inv: dict) -> list:
    pat = re.compile(r"hdx|deuter|uptake|peptide|exchange", re.I)
    return [s["name"] for s in inv["sheets"] if pat.search(s["name"]) or any(pat.search(h) for h in s["header_strings"])]


def pdf_facts(pdf: Path) -> dict:
    text = subprocess.run(["pdftotext", "-layout", str(pdf), "-"], capture_output=True, text=True, check=True).stdout
    flat = re.sub(r"\s+", " ", text)
    cap = re.search(r"\(B\) Full list of peptides with corresponding relative deuterium uptake at different exposure time points\.", flat)
    imgs = subprocess.run(["pdfimages", "-list", str(pdf)], capture_output=True, text=True).stdout.splitlines()[2:]
    tab = {k: re.search(p, flat) for k, p in {
        "time_course_min": r"HDX time course ([\d., ]+) minutes", "controls": r"HDX controls (\d+ undeuterated) (\d+ undeuterated) (\d+ undeuterated)",
        "identified": r"(\d+) identified (\d+) identified (\d+) identified", "followed": r"(\d+) followed c (\d+) followed c (\d+) followed c",
        "coverage": r"Sequence coverage (\S+) (\S+) (\S+)", "replicates": r"Replicates (\d+ technical per state)",
        "repeatability": r"Repeatability b (±[\d.]+ relative Da)", "meaningful": r"Meaningful (>[\d.]+ Da)"}.items()}
    return {"supplementary_table_1": {k: (list(m.groups()) if m and m.lastindex and m.lastindex > 1 else (m.group(1) if m else None)) for k, m in tab.items()},
            "supp_fig_4B_caption_found": bool(cap), "raster_images_in_pdf": len(imgs),
            "text_layer_lines_that_look_like_peptide_rows": sum(bool(re.match(r"^\s*[A-Z]{6,}\b", ln)) for ln in text.splitlines()),
            "no_statistical_tests_statement": "No statistical tests were applied to the HDX-MS measurements" in flat}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fetch-dir", required=True)
    args = ap.parse_args(argv)
    d = Path(args.fetch_dir)
    d.mkdir(parents=True, exist_ok=False)
    proj = json.loads(get(PRIDE))
    files = json.loads(get(PRIDE + "/files?pageSize=200"))
    rec = next(f for f in files if f["fileName"] == CSV_NAME)
    url = next(x["value"] for x in rec["publicFileLocations"] if "FTP" in x["name"]).replace("ftp://", "https://")
    raw = get(url)
    (d / CSV_NAME).write_bytes(raw)
    dg = digests(raw)
    provider = {"project": "PXD034008", "licence": proj.get("license"), "submission_type": proj.get("submissionType"),
                "publication_date": proj.get("publicationDate"), "file_record": {"fileName": rec["fileName"], "category": rec["fileCategory"]["value"],
                "fileSizeBytes": rec["fileSizeBytes"], "checksum_sha1": rec["checksum"], "total_downloads_at_check": rec.get("totalDownloads")},
                "download_url": url, "size_matches": dg["bytes"] == rec["fileSizeBytes"], "sha1_matches": dg["sha1"] == rec["checksum"]}
    if not (provider["size_matches"] and provider["sha1_matches"]):
        raise SystemExit("STOP: the downloaded file does not match the provider record")
    desc = describe_export(parse_export(raw.decode("ascii")), read_fasta(REF))
    desc["file_shape"] = file_shape(raw)
    sd, si = get(ESM.format(n=6, ext="xlsx")), get(ESM.format(n=1, ext="pdf"))
    (d / "MOESM6.xlsx").write_bytes(sd)
    (d / "MOESM1.pdf").write_bytes(si)
    wbinv = xlsx_inventory(sd)
    present = {"Sequence": True, "Start": True, "End": True, "MHP": True, "MaxUptake": True, "RT": True}
    absent = ["protein state or condition label (WT, A344V, M472I)", "exchange time", "deuterium uptake value", "replicate identifier or replicate count",
              "uptake error or standard deviation", "undeuterated control labels", "charge state", "any per-peptide value that differs between states"]
    from_paper = {"states": "WT, A344V and M472I Abl kinase domain (Ayaz et al. 2023 Methods)", "exchange times": "0.167, 1, 10, 60 and 240 minutes (Supplementary Table 1)",
                  "replicates": "2 technical per state (Supplementary Table 1)", "repeatability": "±0.12 relative Da (Supplementary Table 1)",
                  "meaningful difference": "> 0.5 Da, no statistical test (Supplementary Table 1)", "uptake": "relative, not corrected for back-exchange (Methods)",
                  "peptide counts": "59, 66 and 73 identified; 76, 74 and 76 followed (Supplementary Table 1)"}
    out = {"id": "ABL-HDX-LEAD-INVENTORY-2026-10-25", "inspected_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "scope": "data inventory of small public files; no raw mass-spectrometry file or trajectory fetched; no analysis of any endpoint",
           "provider_record": provider, "csv": {"file": CSV_NAME, **dg, "inspection": desc, "fields_present": present, "fields_absent": absent},
           "fields_stated_only_in_the_paper": from_paper,
           "source_data_workbook": {"url": ESM.format(n=6, ext="xlsx"), **digests(sd), **wbinv, "sheets_with_hdx_terms": hdx_terms(wbinv),
                                    "note": "structure only; sheets are binding kinetics, FoldX, melting temperatures and IC50 values"},
           "supplementary_information_pdf": {"url": ESM.format(n=1, ext="pdf"), **digests(si), **pdf_facts(d / "MOESM1.pdf")},
           "local_copies": "outside the repository (--fetch-dir); not committed"}
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    c = desc
    print(json.dumps({"provider_match": [provider["size_matches"], provider["sha1_matches"]], "n_peptides": c["n_peptides"],
                      "coverage": c["file_numbering"], "reference": {k: c["reference_check"][k] for k in ("matched_in_wild_type", "matched_only_in_a_variant", "unmatched", "offsets_seen")},
                      "mhp_max_err": c["mhp_is_monoisotopic_MH_plus"]["max_abs_difference_Da"], "workbook_sheets": wbinv["n_sheets"],
                      "workbook_hdx_sheets": out["source_data_workbook"]["sheets_with_hdx_terms"], "pdf": {k: out["supplementary_information_pdf"][k] for k in ("supp_fig_4B_caption_found", "raster_images_in_pdf", "text_layer_lines_that_look_like_peptide_rows")}}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
