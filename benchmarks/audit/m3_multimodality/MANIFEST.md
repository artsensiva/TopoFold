# M3 Multimodality Detection: Artifact Manifest

**Generated:** 2026-10-03
**Repository Branch:** scientific-rebuild
**Repository HEAD:** 95244fd486be1a96e29cc0cb37bb5cb244fa9407

This manifest records the provenance, checksums, and regenerability status of all promoted M3 scientific audit artifacts.

---

## STATUS AND NAVIGATION DOCUMENTS

| Path | Purpose | SHA-256 |
|------|---------|---------|
| `M3_STATUS.md` | Current M3 status and classification | `491491f...b35f3a` |
| `README.md` | Navigation and reproducibility | `cf0703d...c1c9a` |
| `MANIFEST.md` | This provenance manifest | (self) |

---

## CHECKPOINT REPORTS

| Path | Checkpoint | Status | SHA-256 |
|------|-----------|--------|---------|
| `reports/M3_A2f_CHECKPOINT_REPORT.md` | **M3-A2f** | **CURRENT** | `7e2fb3b...67bb07` |
| `reports/M3_A2f_CLOSURE.md` | **M3-A2f** | **CURRENT** | `7122cb6...659ee9` |
| `reports/M3_A2f_LIKELIHOOD_DEGENERACY.md` | **M3-A2f** | **CURRENT** | `dc249a3...cf5dac` |
| `reports/M3_A2e_CHECKPOINT_REPORT.md` | M3-A2e | Superseded (synthesized) | `9b1d488...a5d21` |
| `reports/M3_A2d_CHECKPOINT_REPORT.md` | M3-A2d | Supporting | `0f08aa0...5975e7` |
| `reports/M3_A2c_CHECKPOINT_REPORT.md` | M3-A2c | Supporting | `23892a3...33471c` |
| `reports/M3_A2b_CHECKPOINT_REPORT.md` | M3-A2b | Superseded | `252a7f6...88ef67` |
| `reports/M3_A2_CHECKPOINT_REPORT.md` | M3-A2 | Superseded | `155167d...140072` |
| `reports/M3_A_CHECKPOINT_REPORT.md` | M3-A | Superseded | `dbe6054...7fbdfac69` |
| `reports/M3_A2_vs_M3_A2b_CORRECTIONS.md` | M3-A2→M3-A2b | Supporting | `dd2c922...55170c` |

---

## DIAGNOSTIC SCRIPTS

| Path | Purpose | Regenerable | SHA-256 |
|------|---------|-------------|---------|
| `scripts/product_vm_fast_accurate.py` | **Definitive** corrected fitter (U estimator) | Yes | `bb1a12f...260f47` |
| `scripts/product_vm_three_estimators.py` | **Definitive** U/C/P framework | Yes | `aee65e4...ddfac69` |
| `scripts/m3_a2f_analysis_framework.py` | Corrected analysis tools | Yes | `626795...cb0a0181` |
| `scripts/investigate_wrapped_normal_anomaly.py` | Convergence investigation | Yes (stochastic) | `7a02cc7...914471` |
| `scripts/m3_a2e_full_benchmark.py` | M3-A2e benchmark orchestration | Yes (stochastic) | `182d7c1...d777adc` |
| `scripts/m3_a2d_corrected_fitter.py` | Historical corrected fitter | Yes | `7cb6e66...cfba891` |
| `scripts/m3_a2c_sanity_audit.py` | Historical sanity audit | Yes | `e60995b...7341814` |
| `scripts/m3_a2b_streamlined.py` | Historical validation | Yes | `898d343...271ec58f` |

---

## RESULT FILES

| Path | Purpose | SHA-256 |
|------|---------|---------|
| `results/m3_a2e_checkpoint6_exact_family.csv` | Exact Product-vM K=1 nulls (3/2250 = 0.13%) | `4b4e89e...aba20fa31` |
| `results/m3_a2e_checkpoint9_k2_recovery.csv` | K2 power (1620/1620 = 100% for ≥60°) | `28c73b6...b7dd1821` |
| `results/m3_a2e_checkpoint10_wrapped_normal.csv` | Wrapped-normal σ=1.0 anomaly | `73eaabc...d36767b267` |
| `results/m3_a2c_exact_null_results.csv` | Early exact-null validation | `4d28b21...dc4622c1` |

---

## REPRODUCIBILITY NOTES

### Python Environment
- Required: NumPy (≥1.20), SciPy (≥1.7)
- **Limitation:** Versions not pinned; platform differences may affect exact numerical results
- **Standard:** Distributional reproducibility (not bitwise-identical)

### Invocation from Repository Root
```bash
# Exact-family validation
python benchmarks/audit/m3_multimodality/scripts/m3_a2e_full_benchmark.py

# Convergence investigation
python benchmarks/audit/m3_multimodality/scripts/investigate_wrapped_normal_anomaly.py
```

### Missing Provenance
- Repository HEAD at analysis: Most historical checkpoints conducted in /tmp without recording repository state
- Current HEAD `95244fd...` is promotion HEAD, not analysis HEAD

---

## KEY CORRECTIONS APPLIED DURING PROMOTION

**M3_STATUS.md** (SHA-256: `491491f...`):
- Corrected K2 power: 1620/1620 (was incorrectly 750/750)
- Corrected U terminology: "numerical-ceiling local-EM heuristic" (not "unrestricted")
- Corrected P description: "not yet selected/validated" (removed unverified Gamma claim)
- Corrected initialization: "split + 2 random" (not "3 random")
- Corrected dependent-toroidal examples: specific verified families

**M3_A2e_CHECKPOINT_REPORT.md** (SHA-256: `9b1d488...`):
- Added SYNTHESIZED provenance metadata
- Corrected exact-family design: κ={1,5,20}, F={100,250,500,1000,2500}
- Corrected K2 recovery design and denominator: 1620/1620 (not 750/750)
- Removed "false discovery" terminology (replaced with "K2-selection frequency")

---

**Promotion Date:** 2026-10-03
**Promotion HEAD:** 95244fd486be1a96e29cc0cb37bb5cb244fa9407
