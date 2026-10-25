# KRAS K1 crystal endpoint report (step S1, 2026-10-25)

**Decision: K1_ENDPOINT_INADEQUATE** — application closed.

Generated from [`KRAS-K1_endpoint.json`](KRAS-K1_endpoint.json) by `plot_k1_endpoint.py`; nothing is typed by hand. Figure:
[`../figures/KRAS-K1_endpoint_residues_2026-10-25.png`](../figures/KRAS-K1_endpoint_residues_2026-10-25.png). No MD trajectory was downloaded, opened or scored.

## Result
- **Endpoint E:** 8 residues of 153 scored: 31–33, 35–37, 59, 74.
- **Merged segments** (gap of up to 2 residues): 31–33, 35–37 (6 positives); 59 (1 positive); 74 (1 positive).
- **Gate (frozen, inclusive):** |E| between 8 and 38: met;
  at least 2 segments of at least 3 residues: **not met** (1 found).
- **Boundary residues** (within 1e-9 Å of the threshold): none.

## Why the endpoint is not adequate
Only 1 of the 3 segments has at least 3 residues, and the frozen gate needs 2. The proxy therefore does not contain a sufficiently distributed comparison target for a method comparison. K1 is closed with this reason; no alternative chain, threshold, endpoint or system is tried (protocol section 4.2).
The gate tests only the spread of the target; it never showed that the endpoint could discriminate methods.

## Displacement context (Å)
- Between-state mean: median 0.30, maximum 4.58. Within-state (copy-to-copy) maximum: median 0.36, maximum 4.13.
- 22 residues have a between-state mean of at least 1 Å; 11 have at least 2 × the within-state value.
- The within-state value exceeds the between-state mean at 93 of 153 residues; copy-to-copy variation is therefore a limit on what the two entries can say about the nucleotide, although both values are small at most residues (medians above).
- Residues with a mean of at least 1 Å that are not positive: 30, 47–48, 68–72, 103–107, 122.
- Positives whose threshold is set by the within-state value rather than by 1 Å: 31–32, 74.

## Exclusions and missing residues
- Scored: 153 residues (4–59 and 68–164); the observed set within 1–166 equals the frozen KRAS-STRUCT-2 set.
- Residues not observed (of 1–169) per chain: GDP:A 60–64; GDP:B 60–66; GNP:A 61–63, 169; GNP:B 168–169; GNP:C 168–169; GNP:D 61–67, 168–169; GNP:E 167–169; GNP:F 167–169.
- Not scored because unobserved in at least one chain: 60–67.

## Sensitivity to the predefined chain handling (reported only; never changes E or the gate)
- primary pair only: |E| = 10 (27, 31–33, 35–37, 51, 59, 74); Jaccard with E 0.80; segments of at least 3 residues: 1
- median over pairs: |E| = 8 (31–33, 35–37, 59, 74); Jaccard with E 1.00; segments of at least 3 residues: 1

## Not applicable
No MD download is proposed: the gate failed and K1 is closed.

## Independent check
A separate implementation (Horn quaternion superposition and an independent PDB reader) reproduces the recorded values: maximum absolute difference 5.2e-14 Å over 12 between-state and 16 within-state pairs; endpoint set identical (`KRAS-K1_endpoint_verification.json`).

## Provenance
- **Frozen commit:** `fc6a3fa1e6c8e29246e22e55585e8d4f13fbaf75`; HEAD at run `fc6a3fa1e6c8e29246e22e55585e8d4f13fbaf75`.
- **SHA-256:** protocol `6bb17ade8a906fa440cc48c6cadf8353c6ce43776779a2e911ebc47ae3817686`; parameters `a1b6576fa1ae54b0f6dd9d3e5f3dcedd34f09e356f910792dc7506cfd30284fb`; endpoint code `7051dc351d6b1b9602bcd0d7bfd1fba32473d75bb5e9f403751cdf9ff65959ab`;
  tests `test_kras_k1_endpoint.py 25560189d6bd7987778f7ba5063ae31f984b6a792cd8ba8e30f4b94f85628899, test_kras_k1_protocol.py cb5094bbed529803f46b8e7639ac9ec01ed91014467bdc348802af87415c59e7`; freeze record `24bd7b96afe77b695ada63e8ecd5e5d2c5303b02be67b1d7e6987177a0e18382`.
- **Inputs (equal to the KRAS-STRUCT-2 spec):** 5US4 `a1d0c9105887a5b64deacf92d17ab5e62a2cdb1d3913f48b8203e9c84f36d22d`; 5USJ `dab6b79531b0c046050958b914e380d8ba613b8e995522acb5ea323f974dbefe`.
- **Endpoint record:** `245407ffaf8e82ea3721c20a54c39d0eea905c1fb0548b7ebf666a400875976c`.
- **Resources:** 0.22 CPU-s, 0.05 s wall, 42 MB peak memory; numpy 2.5.3, Python 3.14.4.
- **Evidence class:** external static structural proxy, already exposed during KRAS-STRUCT-2 development; not held out, not prospective, not independent validation, not functional ground truth.
