"""Minimal stdlib reader for .xlsx cell values (no formulas evaluated, no dependency)."""
import re
import zipfile
import xml.etree.ElementTree as ET

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main", "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def _col(ref):
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group():
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_xlsx(path):
    """Returns {sheet name: list of rows (lists; None for empty cells)}."""
    z = zipfile.ZipFile(path)
    sst = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
            sst.append("".join(t.text or "" for t in si.iter("{%s}t" % NS["m"])))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    out = {}
    for sh in wb.find("m:sheets", NS):
        tgt = rels[sh.get("{%s}id" % NS["r"])]
        tgt = tgt.lstrip("/")
        root = ET.fromstring(z.read(tgt if tgt.startswith("xl/") else "xl/" + tgt))
        rows = []
        for row in root.iter("{%s}row" % NS["m"]):
            cells = {}
            for c in row.findall("m:c", NS):
                v = c.find("m:v", NS)
                t = c.get("t")
                if t == "inlineStr":
                    val = "".join(x.text or "" for x in c.iter("{%s}t" % NS["m"]))
                elif v is None:
                    val = None
                elif t == "s":
                    val = sst[int(v.text)]
                elif t in ("str", "e"):
                    val = v.text
                elif t == "b":
                    val = bool(int(v.text))
                else:
                    val = float(v.text)
                cells[_col(c.get("r"))] = val
            if cells:
                rows.append([cells.get(i) for i in range(max(cells) + 1)])
        out[sh.get("name")] = rows
    return out
