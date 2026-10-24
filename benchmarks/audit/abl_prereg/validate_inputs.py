"""Validate a delivered ABL trajectory inventory and/or processed HDX table against the preregistration schemas.

    .venv/bin/python benchmarks/audit/abl_prereg/validate_inputs.py --inventory INVENTORY.json
    .venv/bin/python benchmarks/audit/abl_prereg/validate_inputs.py --hdx UPTAKE.csv

Checks only structure and consistency; it computes no comparison and no score. Exit code 0 if valid, 1 otherwise.
Requires `jsonschema`.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import sys
from pathlib import Path

import jsonschema

HERE = Path(__file__).resolve().parent
INV_SCHEMA = json.loads((HERE / "schemas/trajectory_inventory.schema.json").read_text())
HDX_SCHEMA = json.loads((HERE / "schemas/hdx_uptake.schema.json").read_text())
MAPPING = HERE / "hdx_peptide_mapping.json"
FASTA = HERE / "reference/UniProt_P00519.fasta"
MUTATIONS = {"M472I": {472: ("M", "I")}, "A344V": {344: ("A", "V")}, "WT": {}}


def _errors(validator, doc):
    return [f"{'/'.join(map(str, e.absolute_path)) or '(root)'}: {e.message}" for e in validator.iter_errors(doc)]


def validate_inventory(doc: dict) -> dict:
    errs = _errors(jsonschema.Draft202012Validator(INV_SCHEMA), doc)
    if errs:
        return {"valid": False, "errors": errs}
    systems = {s["system_id"]: s for s in doc["systems"]}
    if len(systems) != len(doc["systems"]):
        errs.append("duplicate system_id")
    ids = [r["run_id"] for r in doc["runs"]]
    if len(set(ids)) != len(ids):
        errs.append("duplicate run_id")
    for r in doc["runs"]:
        if r["system_id"] not in systems:
            errs.append(f"run {r['run_id']}: unknown system_id {r['system_id']}")
            continue
        expect = r["length_us"] * 1000.0 / r["frame_interval_ns"]
        if abs(r["n_frames"] - expect) > max(2.0, 0.01 * expect):
            errs.append(f"run {r['run_id']}: n_frames {r['n_frames']} inconsistent with length/interval ({expect:.0f})")
    for s in systems.values():
        if s["construct_last_cAbl"] < s["construct_first_cAbl"]:
            errs.append(f"system {s['system_id']}: construct range reversed")
        if s["numbering"] != "cAbl_1a" and "numbering_offset_to_cAbl_1a" not in s:
            errs.append(f"system {s['system_id']}: numbering is {s['numbering']} but no offset to c-Abl 1a is given")
    # comparable WT/M472I pairs: same start structure, same ligands, same force field and temperature
    key = lambda s: (s["starting_structure_pdb"].upper(), tuple(sorted(s["ligands"])), s["force_field"], s["temperature_K"])  # noqa: E731
    groups = collections.defaultdict(dict)
    for s in systems.values():
        groups[key(s)].setdefault(s["variant"], []).append(s["system_id"])
    nruns = collections.Counter(r["system_id"] for r in doc["runs"])
    pairs = [{"start": k[0], "ligands": list(k[1]) or ["none (apo)"], "force_field": k[2], "temperature_K": k[3],
              "WT_runs": sum(nruns[x] for x in v["WT"]), "M472I_runs": sum(nruns[x] for x in v["M472I"])}
             for k, v in groups.items() if "WT" in v and "M472I" in v]
    blocking = [p for p in doc["permissions"] if doc["permissions"][p]["status"] in ("refused", "unknown")
                and p in ("internal_evaluation", "publication_of_derived_results")]
    return {"valid": not errs, "errors": errs, "comparable_WT_M472I_pairs": pairs,
            "apo_pair_available": any(p["ligands"] == ["none (apo)"] for p in pairs),
            "permissions_blocking_evaluation": blocking}


def _coerce(row: dict) -> dict:
    out = {}
    for k, v in row.items():
        if v is None or v == "":
            continue
        if k in ("construct_start", "construct_end", "replicate", "max_uptake_theoretical"):
            out[k] = int(v)
        elif k in ("labeling_time_s", "uptake_Da", "uptake_sd_Da"):
            out[k] = float(v)
        elif k in ("charge_states_combined", "back_exchange_corrected"):
            out[k] = {"true": True, "false": False}.get(str(v).strip().lower(), v)
        else:
            out[k] = v
    return out


def validate_hdx(rows: list[dict]) -> dict:
    v = jsonschema.Draft202012Validator(HDX_SCHEMA)
    errs = []
    data = []
    for i, raw in enumerate(rows):
        row = _coerce(raw)
        e = _errors(v, row)
        errs += [f"row {i + 1}: {x}" for x in e]
        if not e:
            data.append(row)
    if errs:
        return {"valid": False, "errors": errs[:50], "n_errors": len(errs)}
    off = json.loads(MAPPING.read_text())["offset_construct_to_cAbl_1a"] if MAPPING.exists() else None
    seq = "".join(l.strip() for l in FASTA.read_text().splitlines() if not l.startswith(">"))
    keys = collections.Counter()
    reps = collections.defaultdict(set)
    sd = collections.defaultdict(bool)
    for r in data:
        if r["construct_end"] < r["construct_start"] or len(r["peptide_sequence"]) != r["construct_end"] - r["construct_start"] + 1:
            errs.append(f"{r['peptide_sequence']}: length does not match construct_start..construct_end")
        if off is not None:
            ref = list(seq[r["construct_start"] + off - 1: r["construct_end"] + off])
            for pos, (wt, mut) in MUTATIONS[r["state"]].items():
                if r["construct_start"] + off <= pos <= r["construct_end"] + off:
                    ref[pos - (r["construct_start"] + off)] = mut
            if "".join(ref) != r["peptide_sequence"]:
                errs.append(f"{r['peptide_sequence']} ({r['state']}): sequence does not match P00519 at offset {off}")
        k = (r["peptide_sequence"], r["construct_start"], r["state"], r["ligands"], r["labeling_time_s"])
        keys[k + (r["replicate"],)] += 1
        reps[k].add(r["replicate"])
        sd[k] |= "uptake_sd_Da" in r
    errs += [f"duplicate row {k}" for k, n in keys.items() if n > 1]
    errs += [f"fewer than 2 replicates and no uptake_sd_Da: {k}" for k, s in reps.items() if len(s) < 2 and not sd[k]]
    # every compared state must cover the same (peptide position, ligands, time) set
    cover = collections.defaultdict(set)
    for (pep, start, state, lig, t) in reps:
        cover[(state, lig)].add((start, len(pep), t))
    states = sorted(cover)
    for a in states:
        for b in states:
            if a < b and a[1] == b[1] and cover[a] != cover[b]:
                errs.append(f"states {a} and {b} differ in peptide/time coverage ({len(cover[a] ^ cover[b])} items)")
    return {"valid": not errs, "errors": errs[:50], "n_errors": len(errs), "n_rows": len(data),
            "states": [f"{s}|{l}" for s, l in states],
            "labeling_times_s": sorted({r["labeling_time_s"] for r in data})}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inventory")
    ap.add_argument("--hdx")
    a = ap.parse_args(argv)
    if not (a.inventory or a.hdx):
        ap.error("give --inventory and/or --hdx")
    ok = True
    if a.inventory:
        r = validate_inventory(json.loads(Path(a.inventory).read_text()))
        print(json.dumps({"inventory": r}, indent=1))
        ok &= r["valid"]
    if a.hdx:
        r = validate_hdx(list(csv.DictReader(open(a.hdx))))
        print(json.dumps({"hdx": r}, indent=1))
        ok &= r["valid"]
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
