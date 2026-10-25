# Next-application feasibility (2026-10-25): metadata only

The decision record is [`docs/NEXT_APPLICATION_FEASIBILITY_2026-10-25.md`](../../../docs/NEXT_APPLICATION_FEASIBILITY_2026-10-25.md).
Nothing here scores a method or freezes a protocol.

| File | Content |
|---|---|
| `collect_metadata.py` | Collects public metadata (stdlib only; run with `python3 -I`):<br>• the Zenodo record of the KRAS G12X MD, with file sizes and MD5;<br>• two XTC frame headers from each of 4 trajectories, read with HTTP range requests;<br>• one 38 kB topology file, fetched into `--fetch-dir` outside the repository, for the sequence mapping to UniProt P01116 isoforms;<br>• RCSB metadata for the candidate endpoint entries and the 120 Moonshot PDB entries;<br>• the UniProt nsp5 range;<br>• S3 subset sizes for the Moonshot projects 13430 and 13432 |
| `probe_moonshot_h5.py` | Run with the MDAnalysis environment's `python -I`, which has h5py. Reads the HDF5 metadata, the time axis and the topology of two Moonshot files through HTTP range requests (hard cap 8 MB per file; coordinates are never read). Maps the sequence to UniProt P0DTD1 nsp5. |
| `results/feasibility_metadata.json`, `results/moonshot_h5_probe.json` | Generated records |

```bash
python3 -I benchmarks/audit/next_app_feasibility/collect_metadata.py --fetch-dir <new empty directory outside the repo>
benchmarks/audit/wpd_local_retrieval/build/venv/bin/python -I benchmarks/audit/next_app_feasibility/probe_moonshot_h5.py
```

Licences are those stated by the sources:
- Folding@home AWS data: CC0;
- PDB: CC0;
- the Zenodo KRAS record: CC BY 4.0;
- UniProt: CC BY 4.0.

Nothing was sent externally.
