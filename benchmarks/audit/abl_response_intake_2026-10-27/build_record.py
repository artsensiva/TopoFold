"""Builds results/abl_response_intake_record.json (dataset-intake record + inventory) from results/inputs_inspection.json,
the raw files' hashes and the provenance facts recorded below (each with its source). Run: python3 -I build_record.py
"""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ins = json.loads((HERE / "results/inputs_inspection.json").read_text())
zen = json.loads((HERE / "work/zenodo.json").read_text())
t = ins["endpoint_table"]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


raw = []
for p in sorted((HERE / "data_raw").rglob("*")):
    if p.is_file() and "hauser_s5_extracted" not in p.parts:
        raw.append({"path": str(p.relative_to(HERE)), "bytes": p.stat().st_size, "sha256": sha(p)})
ledger = {
    "kept_raw_bytes": sum(r["bytes"] for r in raw),
    "transient_not_kept_bytes": {"bioRxiv_239012_v5_full_pdf (read once, overwritten)": 7101532, "Nature_pdf_second_fetch_duplicate_of_kept_file": 1647637, "PMC_HTML_challenge_page": 1816},
    "metadata_only_fetches": ["Crossref JSON x2", "Europe PMC search JSON", "Zenodo record JSON (3,975 bytes)", "RCSB entry JSON x8", "GitHub API (tree, commits, README)", "HTTP headers of ACS supplementary URLs (HTML challenge pages, 5.6 kB each)"],
    "not_downloaded": [{"item": "Zenodo 10.5281/zenodo.3350897 Data_S3.zip", "bytes": zen["files"][0]["size"], "reason": "831,795,523 bytes, 16.6 times the 50 MB ceiling"},
                       {"item": "Aldeghi et al. Data S1 (XLSX) and Data S2 (ZIP) on pubs.acs.org", "bytes": None, "reason": "served an HTML challenge page to a non-browser client; size unknown; owner browser download possible"}],
}
ledger["scientific_download_total_bytes_upper_bound"] = ledger["kept_raw_bytes"] + sum(ledger["transient_not_kept_bytes"].values())
rec = {
    "schema": "topofold.dataset-intake/1",
    "dataset": {"name": "Abl kinase inhibitor-response benchmark of Hauser et al. 2018 (144 variant-drug records), reused by Aldeghi et al. 2019", "slug": "abl_response_intake_2026-10-27", "type": "compiled literature endpoint (cell-viability IC50 and binding constants) plus wild-type structure models",
                "description": "144 Abl mutant-inhibitor ΔΔG values derived from published IC50 (and, for two drugs, Kd) values, with Prime and FEP+ calculations; eight wild-type Abl:inhibitor models. Aldeghi et al. use the same 144 records."},
    "source": {"primary_reference": ["Hauser K. et al., Commun. Biol. 1, 70 (2018)", "Aldeghi M., Gapsys V., de Groot B. L., ACS Cent. Sci. 5, 1468-1474 (2019)"],
               "doi_or_accession": ["10.1038/s42003-018-0075-x (Crossref: published 2018-06-13)", "10.1021/acscentsci.9b00590 (Crossref: published 2019-08-13)", "10.5281/zenodo.3350897 (Zenodo record, Data S3)", "PMC6716344 (Europe PMC full text of Aldeghi)", "PMC6110136 (Hauser; page not retrievable, challenge page)"],
               "url": "https://github.com/kehauser/Predicting-resistance-of-clinical-Abl-mutations-to-targeted-kinase-inhibitors-using-FEP",
               "version_or_release": "GitHub master at commit d5c6c508acbfe6fe8a47969c7cf93cbe069c065b (2018-04-25T22:45:22Z), last commit; no releases", "accessed_utc": "2026-10-27"},
    "licence": {"data_licence_as_stated": None, "licence_url": None,
                "quoted_text": "GitHub repository: licence field null and no licence file in the 6-file tree. Hauser article: CC BY 4.0 (Crossref licence URL https://creativecommons.org/licenses/by/4.0; article text 'licensed under a Creative Commons Attribution 4.0 International License'). Aldeghi article: 'open access article published under a Creative Commons Attribution (CC-BY) License' (Europe PMC XML). Zenodo Data S3: cc-by-4.0 (record metadata).",
                "paper_licence": "CC BY 4.0 (both articles); the licence of the Hauser supplementary files as deposited on GitHub, of the ACS Data S1/S2 and of the underlying third-party tables is not established by these statements"},
    "permissions": {
        "internal_evaluation": {"status": "conditional", "evidence": "Hauser Data availability: 'All relevant data are publicly available ... at GitHub'; no licence text for the repository files."},
        "publication_of_derived_results": {"status": "conditional", "evidence": "Derived numbers with citation are covered by the CC BY articles; the repository files' own terms are unstated."},
        "redistribution_of_derived_files": {"status": "unknown", "evidence": "No licence for the GitHub files; third-party source values (Schrock 2013 abstract, Pemovska 2015, Davis 2011) have their own terms."},
        "later_commercial_use": {"status": "unknown", "evidence": "Not stated for the repository files or the compiled third-party values; CC BY articles permit commercial reuse of the article."}},
    "size_and_download": {"total_bytes": ledger["scientific_download_total_bytes_upper_bound"], "files_listed": 6, "downloaded": raw, "approval_required": False, "approval_obtained": None,
                          "budget_note": "ceiling 50 MB set by the task; ledger in `download_ledger`"},
    "download_ledger": ledger,
    "raw_inventory": {"location": "benchmarks/audit/abl_response_intake_2026-10-27/data_raw/ (git-ignored)", "git_ignored": True, "files": raw,
                      "extracted_from_Supplementary_Data_5": {"files": 63, "note": "extracted to data_raw/hauser_s5_extracted/ (git-ignored); archive sha256 in files list", "structure_file_sha256": {k: v["sha256"] for k, v in ins["structures"].items()}}},
    "independence": {"prior_topofold_use": "none before this intake (no Abl response data was used by TopoFold studies; the Abl structures of KIN-DESC-1, ABL-MYR and structpair demonstrations are different files)",
                     "overlap_with_other_datasets": "Aldeghi et al. 2019 reuse the 144 records compiled by Hauser et al. 2018 (same endpoint, same wild-type crystal structure inputs); the two papers are not independent measurements. Both Abl sets overlap in protein with the de Buhr & Gräter and Ayaz Abl material only by gene.",
                     "labels_from_topofold_output": False, "equal_information_constraints": "every compared method receives the same 144 records, the same eight wild-type templates and the same site list; no method receives a resistant/susceptible label of the held-out site"},
    "endpoint": {"measurement": "ΔΔG = RT ln(IC50_mut / IC50_WT) at 300 K (RT = 0.5961 kcal/mol per the sheets; positive = resistance), derived by the authors from published IC50 (six inhibitors: cell-viability assays per the paper) and published binding constants (erlotinib, gefitinib; Davis et al. 2011, phosphorylated Abl block per this intake's value match). Assumption of the authors: Michaelis-Menten competitive inhibition with K_M >> [ATP], IC50 change reflecting binding change only.",
                 "provider": "Hauser et al. compile three sources: Schrock et al. Blood 122, 3992 (2013) (a conference abstract; five inhibitors, 105 records), Pemovska et al. Nature 519, 102-105 (2015) (axitinib, 26 records), Davis et al. Nat. Biotechnol. 29, 1046-1051 (2011) (erlotinib and gefitinib, 13 records), per Table 1 of the paper; the underlying tables were not downloaded or verified here.",
                 "processed_values_available": True, "replicates": "none supplied per record; the authors estimate experimental variability from inter-laboratory comparison (RMSE 0.36 kcal/mol, 97 ΔΔGs) and from ΔIC50 versus ΔKd (0.81 kcal/mol, 27 mutations)", "significance_rule": "none per record",
                 "conformational": False, "condition_mismatches_with_simulation": ["cell viability or binding endpoint versus static models", "models are non-phosphorylated, the Davis erlotinib/gefitinib values are from the phosphorylated block", "construct: kinase domain S229-K512 of isoform IA models, assays in cells (full-length BCR-ABL or kinase construct, not verified)", "wild-type IC50 not in the sheet; implied per inhibitor"]},
    "replicates": {"conditions": [], "runs_per_condition": {}, "run_lengths": None, "frame_interval": None, "run_independence": "no trajectories exist in the Hauser deposit; the calculations are alchemical free-energy and Prime estimates (outputs in Supplementary Data 3, 4)", "temporal_dependence_notes": None, "incomplete_runs": None},
    "mapping": {"reference_numbering": "UniProt P00519 canonical sequence (ABL1 isoform IA numbering, gatekeeper Thr315), verified against the fetched sequence", "method": "exact residue identity: the wild-type residue of all 31 mutations equals the P00519 residue at that position; every residue of the eight renumbered models equals the P00519 residue at its number",
                "offsets": {"renumbered files": 0}, "mismatches": ["models contain ACE/NME caps; erlotinib and gefitinib models cover 233-502 only, the others 229-512"], "chain_rule": "protein chain = the chain with the most standard residues (B for axitinib and imatinib, A otherwise), as the authors chose", "atoms_available": "all heavy atoms and hydrogens of the protein; ligand coordinates in mol2/sdf in the same frame (minimum ligand-protein heavy-atom distance 2.8-3.2 Å)", "atoms_needed": "Cα for TopoFold windows; side-chain heavy atoms for contact features"},
    "missing_data_and_exclusions": {"unobserved": "no gaps in numbering within the modelled ranges (loops and termini were built by the authors); built regions are model coordinates, not observed", "exclusion_rules": ["records without IC50 and ΔΔG values in Supplementary Data 3 (135 of 279 rows) are not endpoint records", "H396P (axitinib, gefitinib, erlotinib) excluded by the authors (proline unsupported by FEP+)", "double mutants excluded by the authors"],
                                "counts": {"rows_in_sheet": t["rows_in_sheet"], "rows_with_values": t["rows_with_full_values"]}},
    "decision": {"value": "conditional-go", "allowed_values": ["go", "conditional-go", "no-go"],
                 "reasons": ["endpoint is independent of TopoFold and publicly derived; mapping verified", "terms of the repository files unresolved (internal use only; no redistribution until resolved)", "no mutant structures and no trajectories in the deposit"],
                 "open_conditions": ["owner decision on the terms of the GitHub files for anything beyond internal analysis", "optional owner browser download of Aldeghi Data S1 and S2 (not retrievable by script) to check whether mutant models exist"], "decided_utc": "2026-10-27"},
    "sent_externally": "nothing",
    "exposure": "All 144 records, the 31 mutations, the 20 sites and the eight templates were inspected in this intake; they are development material and cannot serve as an untouched test set.",
    "notes": ["Study decision (B) and design are in REPORT.md."],
}
(HERE / "results/abl_response_intake_record.json").write_text(json.dumps(rec, indent=1))
print("wrote record; raw bytes kept", ledger["kept_raw_bytes"], "upper bound", ledger["scientific_download_total_bytes_upper_bound"])
