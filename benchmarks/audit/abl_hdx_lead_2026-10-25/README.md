# Abl HDX-MS public data lead (2026-10-25): inventory only

Decision record: [`docs/ABL_HDX_PUBLIC_DATA_INVENTORY_2026-10-25.md`](../../../docs/ABL_HDX_PUBLIC_DATA_INVENTORY_2026-10-25.md).
No analysis, no raw mass-spectrometry file, no trajectory.

| File | Content |
|---|---|
| `inventory_hdx_lead.py` | Fetches the 3,386-byte DynamX export (size and SHA-1 checked against the PRIDE record), the article's Source Data workbook and Supplementary Information PDF (into a directory outside the repository), and records which fields are present, absent or only stated in the paper. Needs `pdftotext`/`pdfimages` (poppler) |
| `results/hdx_inventory.json` | The generated record |

```bash
python3 -I benchmarks/audit/abl_hdx_lead_2026-10-25/inventory_hdx_lead.py --fetch-dir <new empty directory outside the repo>
```

Terms: PRIDE project licence CC0; the article is CC BY 4.0. Nothing was sent externally.
