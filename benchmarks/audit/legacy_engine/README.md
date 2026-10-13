# LE-1: legacy κ–τ engine evaluation (audit harness, not production)

Protocol `docs/LEGACY_ENGINE_EVALUATION_PROTOCOL.md` (frozen: `results/le_protocol_freeze.json`), Amendment 1
`docs/LEGACY_ENGINE_EVALUATION_AMENDMENT_1.md` (`results/le_amendment1.json`). Report: `REPORT.md`.
AdK outputs are THIRD-PARTY MD derivatives (Zenodo 10.5281/zenodo.5583119); fixtures and probe inputs are SYNTHETIC.

| File | Role |
|------|------|
| `legacy_lib.py` | Research distances (K, A, legacy Fréchet mirror, F_lock, peptide-plane atoms) |
| `test_legacy.py` | SYNTHETIC verification incl. the actual Rust path via the binding |
| `rust_probe/` | Separate cargo workspace: legacy edge cases (empty, NaN, unequal lengths) |
| `run_dev.py` | AdK development run (reuses WP-D loader and rules) |
| `analyse_dev.py`, `test_analyse.py` | Original analysis (preserved; never run on AdK) |
| `analyse_dev_a1.py`, `test_analyse_a1.py` | Amended analysis (Amendment 1) |
| `timing_dev.py` | Single-thread timing and peak RSS |
| `verify_prior_work.py` | Crossref metadata check |
| `make_report.py` | Writes `REPORT.md` from `results/*.json` |

Run order: `test_legacy.py` → `run_dev.py` → `analyse_dev_a1.py` → `timing_dev.py` → `make_report.py`
(venv: `../wpd_local_retrieval/build/venv`). `results/logs/le_run_attempt1.log`: interrupted first attempt
(cause not established; see `results/le_run_attempts.json`).
