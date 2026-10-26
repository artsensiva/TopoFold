"""KIN-GEN-2 step 1: re-check the official KinCoRe sources (definitions, version history, data terms) and record what is stated and what is not.

    python3 benchmarks/audit/kin_gen_2/reference_source_check.py

Reads only HTML text, HTTP headers, the GitHub licence API and the open-access paper XML (small text retrievals; no label file, no code, no supplementary data is
downloaded). Each resource is recorded separately; a missing statement is recorded as unresolved, never as permission. Output: results/reference_source_check.json.
"""

import hashlib
import json
import re
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent


def get(url, headers=None):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0", **(headers or {})})
    try:
        r = urllib.request.urlopen(req, timeout=90)
        return r.status, dict(r.headers), r.read()
    except Exception as e:     # noqa: BLE001
        return None, {}, str(e).encode()


def text(b):
    t = re.sub(r"<script.*?</script>|<style.*?</style>", "", b.decode("utf8", "ignore"), flags=re.S)
    t = re.sub(r"<[^>]+>", " ", t)
    return re.sub(r"\s+", " ", re.sub(r"&[a-z]+;", " ", t))


out = {"id": "KIN-GEN-2-REFERENCE-SOURCE-CHECK", "pages": {}, "github": {}, "paper": {}}
terms_words = ("licen", "terms of use", "copyright", "creative commons", "permission", "commercial", "redistribut")
for u in ("https://dunbrack.fccc.edu/kincore/", "https://dunbrack.fccc.edu/kincore/download", "https://dunbrack.fccc.edu/kincore/about", "https://dunbrack.fccc.edu/kincore/help"):
    st, h, b = get(u)
    t = text(b)
    m = re.search(r"Last updated (\d{4}-\d{2}-\d{2})", t)
    out["pages"][u] = {"status": st, "bytes": len(b), "sha256": hashlib.sha256(b).hexdigest(), "last_updated_text": m.group(1) if m else None,
                       "http_last_modified": h.get("Last-Modified"), "terms_words_present": {w: (w in t.lower()) for w in terms_words}}
home = text(get("https://dunbrack.fccc.edu/kincore/")[2])
out["site_statements"] = {"news_latest": re.findall(r"News - ([A-Za-z]+ \d+, \d{4})\. ([^.]*\.)", home)[:2],
                          "spatial_definition": ("D1 = dist(αC-Glu(+4) Cα, DFG-Phe Cζ); D2 = dist(β3-Lys Cα, DFG-Phe Cζ); three groups DFGin, DFGout, DFGinter; the home page states no numeric boundaries"
                                                 if "D1 = dist" in home else "not found"),
                          "numeric_thresholds_on_site": bool(re.search(r"D1\s*(≤|≥|<|>)\s*\d", home))}
for repo in ("DunbrackLab/Kincore-standalone2", "DunbrackLab/Kincore-standalone3"):
    st, h, b = get(f"https://api.github.com/repos/{repo}")
    d = json.loads(b) if st == 200 else {}
    out["github"][repo] = {"status": st, "license": d.get("license"), "pushed_at": d.get("pushed_at"), "default_branch": d.get("default_branch")}
st, h, b = get("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6452665/fullTextXML")
x = b.decode("utf8", "ignore")
lic = re.findall(r"<license.*?</license>", x, flags=re.S)
tx = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", x))
out["paper"] = {"citation": "Modi V, Dunbrack RL Jr. Defining a new nomenclature for the structures of active and inactive kinases. PNAS 116:6818-6827 (2019), doi:10.1073/pnas.1814279116, PMC6452665",
                "status": st, "licence_statement": re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", lic[0])).strip() if lic else None,
                "sha256": hashlib.sha256(b).hexdigest(),
                "spatial_definition_in_text": "average-linkage hierarchical clustering of the (D1, D2) distances into three groups" if "average linkage hierarchical clustering" in tx else "not found",
                "numeric_boundaries_in_text": bool(re.search(r"D[12][^.]{0,60}(≤|≥|<|>)\s*\d", tx)),
                "supplementary_datasets_named": sorted(set(re.findall(r"Dataset S\d", tx))), "supplementary_downloaded": False}
out["assessment"] = {
    "route_A_official_label_snapshot": {"usable": False, "reasons": [
        "website (data files, 'Last updated 2026-05-25'): no licence or terms statement on the home, download, about or help pages: unresolved",
        "GitHub standalone software repositories: licence field is null (no licence file): unresolved; software terms do not govern the label data in any case",
        "paper text and supplementary datasets: CC BY-NC-ND 4.0 as stated for the article (adaptations not permitted); the licence of the website data is not thereby established",
        "no pinned, versioned label snapshot with stated terms was found (the site is a live resource updated in place)"]},
    "route_B_independent_labels": {"usable": True, "basis": [
        "the published definition gives the two distances (D1, D2) and the method (average-linkage clustering into three groups) but no numeric boundaries, so official boundaries cannot be reproduced and are not claimed",
        "labels are generated here from public PDB coordinates with an independently written implementation; the cluster boundaries are fitted on development data only; no official label, code or table is copied",
        "the groups are named 'D1D2 groups' (in, inter, out) and are never called KinCoRe labels or biological measurements"]},
    "selected_route": "B", "unresolved": ["terms of use of the KinCoRe website data", "terms of use of the KinCoRe GitHub software", "whether the CC BY-NC-ND article licence extends to the supplementary datasets"]}
(HERE / "results/reference_source_check.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps({"pages": {k: v["last_updated_text"] for k, v in out["pages"].items()}, "github": {k: v["license"] for k, v in out["github"].items()}, "paper_licence": out["paper"]["licence_statement"],
                  "numeric_in_paper": out["paper"]["numeric_boundaries_in_text"], "site_numeric": out["site_statements"]["numeric_thresholds_on_site"], "route": out["assessment"]["selected_route"]}, indent=1))
