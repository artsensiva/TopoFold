"""KIN-PAIR-1 panel selection (protocol.json; frozen before this runs).

    python3 benchmarks/audit/kin_pair_1/select_panel.py freeze        # records spec and code hashes
    .venv/bin/python benchmarks/audit/kin_pair_1/select_panel.py run  # needs the freeze record; network: RCSB, UniProt

Uses primary deposition metadata (RCSB) and input-quality checks on the downloaded PDB files (fixed-column parser,
sequence mapping to UniProt). No TopoFold score is computed and no plot is drawn. Every candidate is recorded with the
reason for inclusion or exclusion in results/KIN-PAIR-1_selection.json. Downloaded files go to data_pdb/ (git-ignored).
"""

import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import kinpair as kp  # noqa: E402

FREEZE = HERE / "results/KIN-PAIR-1_freeze_1.json"
OUT = HERE / "results/KIN-PAIR-1_selection.json"
DATA = HERE / "data_pdb"
FROZEN_FILES = ["protocol.json", "PROTOCOL.md", "select_panel.py", "kinpair.py", "results/candidate_targets_verified.json",
                "results/exposure_history.json"]


def sha256_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# --- network (replaceable in tests) -----------------------------------------------------------------------------------

class Net:
    def __init__(self):
        self.bytes = 0
        self.files = 0

    def _get(self, url, data=None, headers=None):
        for k in range(4):
            try:
                req = urllib.request.Request(url, data=data, headers=headers or {})
                return urllib.request.urlopen(req, timeout=60).read()
            except Exception as e:                          # noqa: BLE001
                if k == 3:
                    raise RuntimeError(f"{url}: {e}") from e
                time.sleep(1 + k)
        raise RuntimeError(f"{url}: unreachable")

    def json(self, url):
        try:
            return json.loads(self._get(url))
        except RuntimeError:
            return None

    def search(self, accession, protocol):
        ec = protocol["entry_criteria"]
        q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match", "value": accession}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                "operator": "exact_match", "value": "UniProt"}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "exptl.method", "operator": "exact_match", "value": ec["method"]}},
            {"type": "terminal", "service": "text", "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined", "operator": "less_or_equal",
                "value": ec["max_resolution_angstrom"]}}]},
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": 10000},
                                "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"},
                                         {"sort_by": "rcsb_id", "direction": "asc"}]}}
        d = json.loads(self._get("https://search.rcsb.org/rcsbsearch/v2/query", json.dumps(q).encode(),
                                 {"Content-Type": "application/json"}))
        return [x["identifier"] for x in d["result_set"]], d["total_count"]

    def pdb_text(self, entry):
        """PDB-format file text, or None if RCSB has no PDB-format file for the entry."""
        try:
            raw = self._get(f"https://files.rcsb.org/download/{entry}.pdb")
        except RuntimeError:
            return None
        self.bytes += len(raw)
        self.files += 1
        DATA.mkdir(exist_ok=True)
        (DATA / f"{entry}.pdb").write_bytes(raw)
        return raw.decode("latin-1")

    def uniprot_fasta(self, acc):
        txt = self._get(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta").decode()
        return "".join(txt.splitlines()[1:]), hashlib.sha256(txt.encode()).hexdigest()


# --- entry metadata ---------------------------------------------------------------------------------------------------

def entry_detail(net, entry, accession, protocol):
    """Summary of an entry from RCSB metadata (the target chains, the primary ligand, resolution, space group)."""
    e = net.json(f"https://data.rcsb.org/rest/v1/core/entry/{entry}")
    if not e:
        return {"entry": entry, "usable": False, "reason": "RCSB entry record unavailable"}
    ids = e["rcsb_entry_container_identifiers"]
    chains, lengths, mut = [], [], 0
    for i in ids.get("polymer_entity_ids", []):
        p = net.json(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entry}/{i}") or {}
        refs = p.get("rcsb_polymer_entity_container_identifiers", {}).get("reference_sequence_identifiers", []) or []
        if any(r.get("database_name") == "UniProt" and r.get("database_accession") == accession for r in refs):
            chains += p["rcsb_polymer_entity_container_identifiers"].get("auth_asym_ids", [])
            lengths.append(p.get("entity_poly", {}).get("rcsb_sample_sequence_length"))
            mut += p.get("entity_poly", {}).get("rcsb_mutation_count") or 0
    lig = []
    for i in ids.get("non_polymer_entity_ids", []) or []:
        n = net.json(f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{entry}/{i}") or {}
        cid = (n.get("pdbx_entity_nonpoly") or {}).get("comp_id")
        if cid and cid not in protocol["additive_comp_ids_not_primary_ligand"]:
            lig.append((n.get("rcsb_nonpolymer_entity", {}).get("formula_weight") or 0.0, cid))
    lig.sort(key=lambda t: (-t[0], t[1]))
    res = e["rcsb_entry_info"].get("resolution_combined") or [None]
    return {"entry": entry, "usable": bool(chains), "reason": "" if chains else "no chain of the target entity",
            "target_chains": sorted(set(chains)), "resolution": res[0], "release_date": e["rcsb_accession_info"]["initial_release_date"][:10],
            "space_group": (e.get("symmetry") or {}).get("space_group_name_H_M"), "primary_ligand": lig[0][1] if lig else None,
            "other_ligands": [c for _, c in lig[1:]], "n_mutations_rcsb": mut, "entity_sequence_lengths": lengths}


def entry_quality(net, d, accession, ref_seq, protocol, texts):
    """Download the PDB file, parse it, choose the chain by the chain rule, map it to UniProt, apply the quality rules."""
    q = protocol["quality"]
    txt = net.pdb_text(d["entry"])
    if txt is None:
        return {"ok": False, "reason": "no PDB-format file at RCSB"}
    texts[d["entry"]] = txt
    P = kp.parse_pdb(txt)
    if not P.residues:
        return {"ok": False, "reason": "no polymer residues parsed"}
    cand = [c for c in d["target_chains"] if c in P.residues]
    if not cand:
        return {"ok": False, "reason": f"none of the target chains {d['target_chains']} is in the file"}
    counts = {c: sum(r.read for r in P.residues[c]) for c in cand}
    chain = sorted(cand, key=lambda c: (-counts[c], c))[0]
    m = kp.map_chain(P, chain, ref_seq, accession)
    n_read = sum(P.residues[chain][j].read for j in range(len(P.residues[chain])) if P.residues[chain][j].key in m["index"])
    problems = list(m["problems"])
    if n_read < q["min_read_mapped_residues_per_chain"]:
        problems.append(f"only {n_read} mapped residues are read by the adapter (< {q['min_read_mapped_residues_per_chain']})")
    return {"ok": not problems, "reason": "; ".join(problems), "chain": chain, "chain_read_counts": counts, "n_read_mapped": n_read,
            "identity": m.get("identity"), "coverage": m.get("coverage"), "n_aligned": m.get("n_aligned"),
            "mismatches": m.get("mismatches"), "dbref_segments": m.get("dbref_segments"), "sha256": P.text_sha256,
            "n_models": P.n_models, "mapping": m, "P": P}


# --- selection --------------------------------------------------------------------------------------------------------

def select(net, protocol, targets, exposed_ids, exposed_acc, log):
    rec = {"targets": {}, "ranking": [], "panel": [], "stop": None}
    ranked = []
    for gene, acc in targets:
        if acc in exposed_acc:
            rec["targets"][gene] = {"accession": acc, "status": "excluded", "reason": "target exposed during TopoFold development"}
            continue
        ids, total = net.search(acc, protocol)
        elig = [i for i in ids if i.upper() not in exposed_ids]
        rec["targets"][gene] = {"accession": acc, "search_total": total, "eligible_by_metadata": len(elig),
                                "excluded_exposed_entries": sorted(i for i in ids if i.upper() in exposed_ids)}
        ranked.append((-len(elig), acc, gene, elig))
    ranked.sort(key=lambda t: (t[0], t[1]))
    rec["ranking"] = [{"rank": k + 1, "gene": g, "accession": a, "eligible_by_metadata": -n} for k, (n, a, g, _) in enumerate(ranked)]
    need = protocol["pairs_per_target"] * 2
    for k, (n, acc, gene, elig) in enumerate(ranked):
        t = rec["targets"][gene]
        if len(rec["panel"]) >= protocol["targets_needed"]:
            t["status"] = "not examined (panel complete)"
            continue
        if k >= protocol["max_targets_examined"]:
            t["status"] = "not examined (max_targets_examined reached)"
            continue
        if -n < need:
            t["status"] = "excluded"
            t["reason"] = f"fewer than {need} entries eligible by metadata"
            continue
        ref_seq, ref_sha = net.uniprot_fasta(acc)
        pairs, examined, texts = build_pairs(net, protocol, gene, acc, elig, ref_seq, log)
        t.update({"rank": k + 1, "uniprot_sequence_length": len(ref_seq), "uniprot_fasta_sha256": ref_sha,
                  "examined_entries": [{kk: v for kk, v in x.items() if kk not in ("mapping", "P")} for x in examined]})
        if net.files > protocol["max_files_downloaded_total"] or net.bytes > protocol["max_scientific_download_mb"] * 1e6:
            t["status"] = "stopped"
            t["reason"] = "download limit reached"
            rec["stop"] = "download limit reached"
            break
        if len(pairs) < protocol["pairs_per_target"]:
            t["status"] = "excluded"
            t["reason"] = f"{len(pairs)} usable disjoint pair(s) found (< {protocol['pairs_per_target']})"
            continue
        t["status"] = "selected"
        for pr in pairs:
            rec["panel"].append({"gene": gene, "accession": acc, "uniprot_sequence_length": len(ref_seq), **pr})
    rec["files_downloaded"], rec["bytes_downloaded"] = net.files, net.bytes
    return rec


def build_pairs(net, protocol, gene, acc, elig, ref_seq, log):
    q = protocol["quality"]
    examined, texts, used, pairs = [], {}, set(), []
    info = {}

    def get(entry):
        if entry in info:
            return info[entry]
        if len(examined) >= protocol["max_entries_examined_per_target"] or net.files >= protocol["max_files_downloaded_total"]:
            return None
        d = entry_detail(net, entry, acc, protocol)
        row = {"entry": entry, "resolution": d.get("resolution"), "primary_ligand": d.get("primary_ligand"), "status": "", "reason": ""}
        if not d["usable"]:
            row.update(status="excluded", reason=d["reason"])
        else:
            qq = entry_quality(net, d, acc, ref_seq, protocol, texts)
            row.update({k: v for k, v in qq.items() if k not in ("ok", "reason", "mapping", "P")})
            row.update({k: d[k] for k in ("target_chains", "release_date", "space_group", "other_ligands", "n_mutations_rcsb")})
            row["status"], row["reason"] = ("usable", "") if qq["ok"] else ("excluded", qq["reason"])
            d["quality"] = qq
        examined.append(row)
        info[entry] = (d, row)
        log(f"  {gene} {entry}: {row['status']} {row['reason'][:70]}")
        return info[entry]

    order = list(elig)
    for a_pos, e1 in enumerate(order):
        if len(pairs) >= protocol["pairs_per_target"]:
            break
        if e1 in used:
            continue
        g1 = get(e1)
        if g1 is None:
            break
        d1, r1 = g1
        if r1["status"] != "usable" or not d1["primary_ligand"]:
            if r1["status"] == "usable":
                r1["status"], r1["reason"] = "not paired", "no primary ligand"
            continue
        for e2 in order[a_pos + 1:]:
            if e2 in used:
                continue
            g2 = get(e2)
            if g2 is None:
                break
            d2, r2 = g2
            if r2["status"] != "usable" or not d2["primary_ligand"] or d2["primary_ligand"] == d1["primary_ligand"]:
                continue
            c = kp.correspondence(d1["quality"]["P"], d1["quality"]["mapping"], d1["quality"]["chain"],
                                  d2["quality"]["P"], d2["quality"]["mapping"], d2["quality"]["chain"])
            if c["n_shared"] < q["min_shared_positions"]:
                r2["reason"] = (r2["reason"] + f"; only {c['n_shared']} positions shared with {e1}").strip("; ")
                continue
            used |= {e1, e2}
            pairs.append({"pair_id": f"{gene}-{len(pairs) + 1}", "A": e1, "B": e2, "chain_A": d1["quality"]["chain"],
                          "chain_B": d2["quality"]["chain"], "sha256_A": d1["quality"]["sha256"], "sha256_B": d2["quality"]["sha256"],
                          "n_read_mapped_A": d1["quality"]["n_read_mapped"], "n_read_mapped_B": d2["quality"]["n_read_mapped"],
                          "n_shared_positions": c["n_shared"], "span_not_shared": c["n_span_not_shared"],
                          "meta_A": {k: d1[k] for k in ("resolution", "release_date", "space_group", "primary_ligand", "other_ligands", "n_mutations_rcsb", "target_chains")},
                          "meta_B": {k: d2[k] for k in ("resolution", "release_date", "space_group", "primary_ligand", "other_ligands", "n_mutations_rcsb", "target_chains")}})
            r1["status"] = r2["status"] = "paired"
            break
    return pairs, examined, texts


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "freeze":
        import datetime as dt
        FREEZE.write_text(json.dumps({"study": "KIN-PAIR-1", "stage": "selection protocol and code frozen", "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                                      "sha256": {f: sha256_file(HERE / f) for f in FROZEN_FILES}, "files_present_at_freeze": sorted(p.name for p in DATA.glob("*.pdb")) if DATA.exists() else []}, indent=1) + "\n")
        print("frozen", FREEZE)
        return 0
    if cmd != "run":
        print(__doc__)
        return 2
    fr = json.loads(FREEZE.read_text())
    for f, h in fr["sha256"].items():
        if sha256_file(HERE / f) != h:
            print(f"STOP: {f} differs from the freeze")
            return 1
    protocol = json.loads((HERE / "protocol.json").read_text())
    cands = json.loads((HERE / protocol["candidate_targets_file"]).read_text())
    exp = json.loads((HERE / protocol["exposure_file"]).read_text())
    targets = [(g, v["accession"]) for g, v in cands.items() if v["verified"]]
    import resource
    t0 = time.time()
    rec = select(Net(), protocol, targets, {k.upper() for k in exp["entries"]}, set(exp["exposed_uniprot_accessions"]), lambda s: print(s, flush=True))
    ru = resource.getrusage(resource.RUSAGE_SELF)
    rec["resources"] = {"cpu_seconds": ru.ru_utime + ru.ru_stime, "wall_seconds": time.time() - t0}
    rec["freeze_utc"] = fr["utc"]
    OUT.write_text(json.dumps(rec, indent=1, default=str) + "\n")
    print("panel:", [(p["pair_id"], p["A"], p["B"]) for p in rec["panel"]], "| files", rec["files_downloaded"], "bytes", rec["bytes_downloaded"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
