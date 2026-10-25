"""Source verification for the oncology evaluation search (2026-10-25): metadata and short text excerpts only.

    python3 -I benchmarks/audit/oncology_search_2026-10-25/collect_sources.py

Re-queries public APIs (Crossref, PRIDE, Europe PMC) and records, in results/sources_verified.json, what each cited source
actually states: bibliographic data and licence of the cited articles, the data-availability statements, the PRIDE project
metadata and file listing of the Abl HDX-MS data, and the sentences of the KRAS simulation paper that describe the
biochemical comparison data. Nothing is downloaded beyond API metadata and article text; no trajectory, raw mass-spectrometry
file or coordinate is fetched. Submitter e-mail addresses returned by PRIDE are not stored.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "results/sources_verified.json"
ROOT = HERE.parents[2]
C1_SOURCES = ROOT / "benchmarks/audit/c1_feasibility/results/c1_sources.json"
DOIS = {"ayaz_2023": "10.1038/s41467-023-36956-5", "pantsar_2018": "10.1371/journal.pcbi.1006458",
        "hunter_2015": "10.1158/1541-7786.MCR-15-0203", "rennella_2025": "10.1039/d4cb00233d",
        "parker_2022": "10.1158/0008-5472.CAN-22-0804"}
PMC = {"ayaz_2023": "PMC10076256", "pantsar_2018": "PMC6147662", "rennella_2025": "PMC11629925"}
EXCERPT_MAX = 700


def get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "TopoFold-source-check/1", "Accept": "application/json, */*"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def crossref(doi: str) -> dict:
    m = json.loads(get("https://api.crossref.org/works/" + urllib.parse.quote(doi, safe="/")))["message"]
    return {"doi": doi, "title": m["title"][0], "journal": (m.get("container-title") or [""])[0], "volume": m.get("volume"),
            "pages_or_article": m.get("page") or m.get("article-number"),
            "first_authors": [a.get("family") for a in m.get("author", [])][:3], "year": (m.get("issued", {}).get("date-parts") or [[None]])[0][0],
            "licence_urls": [lic.get("URL") for lic in m.get("license", [])]}


def plain_text(xml: str) -> str:
    import html
    return html.unescape(re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", xml))).strip()


def excerpt(text: str, pattern: str, before: int = 120, after: int = 520) -> str | None:
    m = re.search(pattern, text)
    return text[max(0, m.start() - before):m.start() + after][:EXCERPT_MAX] if m else None


def article(pmcid: str) -> tuple:
    raw = get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML")
    return plain_text(raw.decode("utf8", "replace")), hashlib.sha256(raw).hexdigest(), len(raw)


def hdx_file_summary(names_sizes: list) -> dict:
    """Counts of raw files per protein and exchange time from the PRIDE file names (e.g. ABL_WT_10minD2O)."""
    per, tot = {}, 0
    for name, size in names_sizes:
        tot += size
        m = re.search(r"ABL_(WT|A344V|M472I)_(ctrl\d*|\d+(?:sec|min|h))", name)
        if m and name.endswith(".raw.zip"):
            per.setdefault(m.group(1), []).append(m.group(2))
    return {"raw_files_per_protein": {k: len(v) for k, v in sorted(per.items())},
            "labels_per_protein": {k: sorted(v) for k, v in sorted(per.items())}, "total_bytes": tot}


def pride(acc: str, with_files: bool) -> dict:
    d = json.loads(get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}"))
    rec = {"accession": acc, "title": d.get("title"), "publication_date": d.get("publicationDate"),
           "submitter": f"{d['submitters'][0].get('name')} ({d['submitters'][0].get('affiliation')})" if d.get("submitters") else None,
           "instruments": [i.get("name") for i in d.get("instruments", [])]}
    if with_files:
        f = json.loads(get(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/files?pageSize=200"))
        files = f if isinstance(f, list) else f.get("_embedded", {}).get("files", [])
        rec["files_listed"] = len(files)
        rec["files_by_category"] = {}
        for x in files:
            c = (x.get("fileCategory") or {}).get("value")
            rec["files_by_category"][c] = rec["files_by_category"].get(c, 0) + 1
        rec["non_raw_files"] = [(x["fileName"], x.get("fileSizeBytes")) for x in files if (x.get("fileCategory") or {}).get("value") != "RAW"]
        rec.update(hdx_file_summary([(x["fileName"], x.get("fileSizeBytes") or 0) for x in files]))
    return rec


def main() -> int:
    out = {"id": "ONCOLOGY-SEARCH-2026-10-25-sources", "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "scope": "API metadata and short article excerpts only; no trajectory, raw mass-spectrometry file or coordinate fetched",
           "crossref": {k: crossref(v) for k, v in DOIS.items()}}
    texts = {k: article(v) for k, v in PMC.items()}
    out["articles"] = {k: {"pmcid": PMC[k], "xml_sha256": t[1], "xml_bytes": t[2]} for k, t in texts.items()}
    ay, pa, re_ = texts["ayaz_2023"][0], texts["pantsar_2018"][0], texts["rennella_2025"][0]
    out["statements"] = {
        "ayaz_data_availability": excerpt(ay, r"Data availability", 0, 650),
        "ayaz_hdx_method": excerpt(ay, r"Hydrogen exchange was performed", 0, 700),
        "ayaz_hdx_duplicate_and_relative": excerpt(ay, r"Deuterium levels were not corrected", 0, 330),
        "ayaz_hdx_result": excerpt(ay, r"The results confirmed that M472I increases HDX", 0, 560),
        "ayaz_2f4j_simulations": excerpt(ay, r"Simulations starting from the DFG-in/AL-open conformation", 0, 520),
        "pantsar_hunter_sentences": excerpt(pa, r"An interesting discrepancy among KRAS G12X mutants", 0, 520),
        "pantsar_hunter_raf": excerpt(pa, r"exhibit reduced affinity to Raf", 60, 420),
        "pantsar_replica_variation": excerpt(pa, r"individual replicas displayed variation", 200, 260),
        "rennella_data_availability": excerpt(re_, r"Data availability", 0, 480),
        "rennella_trajectory_statement": excerpt(re_, r"Trajectory frames were written", 0, 150)}
    out["pride"] = {"PXD034008": pride("PXD034008", True), "PXD054924": pride("PXD054924", False)}
    hit = json.loads(get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:26037647%20AND%20SRC:MED&format=json&resultType=core"))["resultList"]["result"][0]
    out["hunter_2015_europepmc"] = {"pmid": "26037647", "pmcid": hit.get("pmcid"), "isOpenAccess": hit.get("isOpenAccess"),
                                    "note": "no open-access full text in Europe PMC: the per-mutant values were not verified against the primary text"}
    c1 = json.loads(C1_SOURCES.read_text())
    out["rascore_project_evidence"] = {"file": "benchmarks/audit/c1_feasibility/results/c1_sources.json",
                                       "sha256": hashlib.sha256(C1_SOURCES.read_bytes()).hexdigest(),
                                       "retrieved_utc": c1["retrieved_utc"], "entries": [f["what"][:160] for f in c1["files"]][:6]}
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    p = out["pride"]["PXD034008"]
    print(json.dumps({"crossref": len(out["crossref"]), "excerpts_found": sum(v is not None for v in out["statements"].values()),
                      "excerpts_total": len(out["statements"]), "pxd034008_files": p["files_listed"], "raw_per_protein": p["raw_files_per_protein"],
                      "total_GB": round(p["total_bytes"] / 1e9, 2)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
