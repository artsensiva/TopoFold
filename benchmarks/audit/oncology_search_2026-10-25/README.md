# Oncology evaluation search (2026-10-25): metadata only

The decision record is [`docs/ONCOLOGY_EVALUATION_SEARCH_2026-10-25.md`](../../../docs/ONCOLOGY_EVALUATION_SEARCH_2026-10-25.md).
Nothing here scores a method, freezes a protocol or downloads a data file.

| File | Content |
|---|---|
| `collect_sources.py` | Re-queries Crossref, PRIDE and Europe PMC and records what the cited sources state (bibliographic data and licences, data-availability statements, the PRIDE file listing of the Abl HDX-MS data, the sentences describing the biochemical comparison data). Standard library only; PRIDE submitter e-mail addresses are not stored |
| `results/sources_verified.json` | The generated record every number in the decision record is read from |

```bash
python3 -I benchmarks/audit/oncology_search_2026-10-25/collect_sources.py
```

Terms are those of the sources: Ayaz et al. (CC BY 4.0), Pantsar et al. (CC BY 4.0), Rennella et al. (CC BY-NC); the Hunter et
al. article has no licence in Crossref and no open text in Europe PMC. Nothing was sent externally.
