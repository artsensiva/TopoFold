# PROJECT STATE — read this first (checkpoint 2026-10-28)

**Status: PAUSED BY OWNER, effective 2026-10-28.** No development task is authorized. A new session reports this status, answers questions and **waits for an explicit owner task**. Nothing listed under "possible future work" is authorized by appearing here.

Latest scientific checkpoint: commit `495298f` (MET exploratory evaluation, decision log D93). The pause documentation is a later documentation-only commit on the same branch (`git log -3`).

## 1. Mission and what is supported
- **Mission (owner hypothesis, untested commercially):** a reproducible kinase structural-evidence and variant-response research workflow built on TopoFold's SE(3)-invariant local geometry of Cα backbones (curvature κ, torsion τ, writhe).
- **Supported capabilities (engineering, with tests; not validated science):** Rust core `topofold` with Python bindings; input adapters (`python/topofold_inputs`: Cα and atom adapters, features, `kinase_roles`, phospho support); `structpair` 0.5.0 (adapter 0.2.0) two-structure comparison with equal-information baselines, hashes and numeric record; `ensemble` descriptive group comparison; local read-only viewer (`apps/structpair_viewer`); frozen audit harnesses in `benchmarks/audit/` (index: its README). Two independent AI technical trials and one confirmation exist; **no independent human trial has been completed**.
- Commands: see §7.

## 2. Principal results (details: `docs/DECISION_LOG.md`, `benchmarks/audit/README.md`)
- **Positive (engineering/mathematical):** writhe normalisation and circular statistics repaired (KNOWN_ISSUES M1, M2); reference checks and tests for geometry; structpair verified on recorded KRAS and Abl comparisons; deterministic MET template-feature exporter (D92).
- **Negative / not met:** v0.8.3 scientific claims withdrawn (synthetic data had been presented as MD; estimator artefacts as free energies). KIN-PAIR-1: frozen-rule outcome IMPLEMENTATION FAILURE (D81). **MET evaluation (D93, one frozen run):** angular κ/τ features and Cα distances did **not** meet the pre-declared improvement criterion beyond ordinary ligand/site context in either model family; the pre-declared rule says do not continue geometric-feature development on that evidence. Conditional on this dataset, the wild-type-template features and the blocked procedure; not equivalence, not proof of no effect, **not a statement about all uses of geometry**.
- **Inconclusive / descriptive:** KIN-GEN-2 insufficient evidence (8 DFG-out entries from 2 kinases; D87/D88); ABL-MYR-1b descriptive result (baselines covered the flagged residues); K1 stopped at an inadequate crystal endpoint; local geometry largely redundant with Cα-distance and local-RMSD baselines in 10 of 12 KIN-PAIR cases; ESM-1b log-likelihood column modestly improved the MET comparison (2–4%, below the 5% threshold). **No ESM campaign is authorized.**
- **Unsupported claims (do not make):** molecular-dynamics or free-energy results from this code; pocket detection; superiority over other methods (including the MET authors' model); drug-response or resistance prediction; clinical or commercial usefulness; human usability validation; commercial or redistribution rights to any third-party data.

## 3. Parked directions (not authorized)
Parked: wild-type-template geometric features for kinase variant-response prediction (D93); Hauser/Aldeghi Abl response table as a pilot (D90/D91). Completed and closed, not to be reopened without an owner decision: KIN-GEN-2, KIN-PAIR-1, K1, ABL-MYR-1/1a/1b, the MET export and evaluation, the PIM1 annotation, the modified-residue and phospho work.
Possible future discussion, subject to a fresh owner decision (not a task): practical structural-analysis value of existing tools (structpair, ensemble, viewer); mutant-specific structural information plus a prospective measurement if kinase response is ever reopened; further KNOWN_ISSUES items.

## 4. Repository state
- Active branch: `scientific-rebuild` (unreleased research work; remote `origin/scientific-rebuild`). `main` is historical: it carries the v0.8.3 text with a README notice pointing to the active branch; local `main` may be behind `origin/main`. **Nothing has been merged; no release exists; tags `v0.8.3`, `v0.8.3-pre-audit` and older are unchanged.** `Cargo.toml` and `pyproject.toml` still read 0.8.3 (historical strings, not a claim). No force-push, no history rewriting, no default-branch change.
- CI: `.github/workflows/ci.yml` (Rust, Python 3.9/3.12/3.14, viewer, adapter tests with MDAnalysis, audit reproductions with the MET evaluator tests). Last scientific commit was green.

## 5. Unresolved issues and external requests (actual status)
- `docs/KNOWN_ISSUES.md` is the backlog; C1–C5 and M3–M8 remain open as written there; M1, M2 resolved; m13, m13b, m14 addressed (D82–D86).
- **DESRES/Ayaz trajectories:** request was SENT by the owner; the reply offers non-commercial access through an eligible institution and leaves publication, derived-feature redistribution and commercial use unanswered; treated as **unavailable for present planning** (D88). **Never resend.**
- **Processed HDX request:** SENT by the owner; **no reply recorded; pending.** Never resend.
- **Data terms:** Hauser GitHub/supplement terms unresolved; MET data used internally on the documented CC BY article / MIT repository basis; no commercial or redistribution rights established (D90–D93).
- Reviewer role: the GPT/Codex reviewer is read-only unless the owner explicitly authorizes an exception.

## 6. Protected material
- **Owner files (never modify, move, stage or delete):** the four root files named "Научно-техническая оценка TopoFold …" (`.md`, `.pdf` ×2, `.md` variant), `docs/REPORT-ChatGPT-Independent-tester-24-10-2026.md`, `…-chat.md`, `…_2.md` (untracked), and the ignored `SomeHistory.md`.
- **Frozen artifacts (do not edit; correct by dated addenda):** KIN-PAIR-1, K1, ABL-MYR-1b, BC detector, statistical-calibration work, AI-trial archives, KIN-GEN-2 (and its addendum), `met_feature_export_2026-10-27` (frozen exporter 6094572), `met_eval_2026-10-28` (frozen at 84215b8, results 495298f), inventory, phospho and PIM1 records.

## 7. Environment, setup and commands (paths relative to the repository root `~/TopoFold`)
- **Main environment:** `.venv` (Python 3.14.4; list in `docs/ENVIRONMENT_2026-10-28_main_venv_pip_freeze.txt`: numpy 2.5.3, scipy 1.18.1, scikit-learn 1.9.1, pandas 3.0.6, matplotlib 3.11.2, pytest 9.1.1, maturin 1.15.0, streamlit 1.64.0). Rust 1.98.1 / cargo 1.98.1.
- **Adapter environment (MDAnalysis 2.10.0):** `benchmarks/audit/wpd_local_retrieval/build/venv` (regenerable; package list `benchmarks/audit/ENVIRONMENT_audit_venv_pip_freeze.txt`).
- **Rebuild:** `python3 -m venv .venv && .venv/bin/pip install maturin numpy pytest scipy scikit-learn matplotlib pandas streamlit && .venv/bin/maturin develop --release`; add MDAnalysis only for adapter tests. `Cargo.lock` is git-ignored (a copy is in the private backup).
- **Lightweight integrity checks (no science rerun):** `git status -sb; git log -3 --oneline; git fsck --no-dangling`; `.venv/bin/python benchmarks/audit/kin_pair_1/verify_freeze.py` (documents known amendment differences); `.venv/bin/python -m pytest tests -q` (476 passed, 17 skipped, about 10 minutes locally on 2026-10-28; it launches a temporary Streamlit server for the viewer browser check, which exits with the test); `cargo test --workspace` after Rust changes.
- **Launch tools:** `PYTHONPATH=python .venv/bin/python -m topofold_inputs.structpair SPEC.json --out DIR`; `PYTHONPATH=python .venv/bin/python -m topofold_inputs.ensemble SPEC.json --out DIR`; viewer `.venv/bin/streamlit run apps/structpair_viewer/app.py --server.address localhost -- --results DIR` (details `docs/STRUCTPAIR_USER_GUIDE.md`). No server or job runs while paused.
- **Local-only inputs:** git-ignored raw inputs (downloaded structures, owner-downloaded Dryad archives, MD datasets, assay tables, derived `work/` folders) are not in Git. Per-study provenance is in each study's intake record and `benchmarks/audit/checkpoint_2026-10-13_manifest.json`. A private backup with hashes was made on 2026-10-28 outside the repository (`~/backups/topofold-2026-10-28/`; its manifest is deliberately not committed). **It is on the same computer; no off-device copy was verified.** See the final checkpoint report for exact status.

## 8. Reading order
`AGENTS.md` (rules) → this file → `README.md` banner → `docs/DECISION_LOG.md` (D87–D93 newest) → `docs/ENGINE_MISSION_AND_ROADMAP.md` → `docs/KNOWN_ISSUES.md` → `benchmarks/audit/README.md` → the relevant study's `REPORT.md`.

## 9. Restart prompt (Claude or GPT/Codex)
> Read `AGENTS.md` and `docs/PROJECT_STATE.md` in the TopoFold repository, then inspect `git status` and `git log -3`. The project is PAUSED BY OWNER: report the checkpoint status, the open questions you see and the files you would need, and wait for an explicit development task from me. Do not start science, features, downloads, training or reruns; do not touch owner files or frozen artifacts; do not resend the DESRES or HDX requests; do not merge, tag, release or force-push. The GPT/Codex reviewer stays read-only unless I explicitly authorize an exception. Follow the scientific-integrity rules and the report format in AGENTS.md.
