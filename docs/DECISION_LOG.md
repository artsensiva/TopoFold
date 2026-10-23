# Decision Log

One line per recorded decision, newest last. Details stay in the linked documents, which are not restated here. The
programme and current queue are in [`ENGINE_MISSION_AND_ROADMAP.md`](ENGINE_MISSION_AND_ROADMAP.md) §12; the
gate procedure for new entries is §12.6. Source paths are relative to `docs/` unless they begin with `benchmarks/`
(repository root).

| ID | Date | Decision | Evidence level supported | Source |
|---|---|---|---|---|
| D01 | 2026-10-04 | Scalar circular-mean/BC scan component: **CHANGE** (restrict; not a state detector) | A | `benchmarks/audit/angular_characterization/REPORT.md` |
| D02 | 2026-10-04 | AdK WP-B G0–G1: unresolved for global-state tasks; recommend WP-D local retrieval | — | `EVALUATION_DESIGN_DECISION_MEMO.md` |
| D03 | 2026-10-12 | Legacy Fréchet coupling deprioritized for same-sequence windows (LE-1; formal outcome inconclusive) | B (development tier) | `benchmarks/audit/legacy_engine/REPORT.md` |
| D04 | 2026-10-12 | WP-D held-out: A retained for the φ/ψ target, B2 for the peptide-plane target, within mdCATH only | B (scoped) | `WPD_HELDOUT_MDCATH_CLOSEOUT_ADDENDUM_2026-10-12.md` |
| D05 | 2026-10-12 | EB-1 Amendment 2: decision **engine** from a projected B2 single-query time | Engineering | `EB1_AMENDMENT_2.md` |
| D06 | 2026-10-13 | BQ-1 acceptance not met; BQ-2 conclusion C (insufficient evidence) | A + engineering | `benchmarks/audit/b2opt/REPORT.md`, `benchmarks/audit/bq2/REPORT.md` |
| D07 | 2026-10-13 | BQ-3: prepared-database operation qualified on SYNTHETIC P2-dimension rows; preparation + first query did not meet 10 s | Engineering | `B2_DIRECT_P2_QUALIFICATION_CLOSEOUT_2026-10-13.md` |
| D08 | 2026-10-13 | CMP-1: no decisive advantage; V2b integration not supported; V2b stays an audit artifact; integration proposal historical | A + engineering | `CMP1_CLOSEOUT_2026-10-13.md`, `B2_V2B_INTEGRATION_PROPOSAL.md` |
| D09 | 2026-10-13 | First biological-validation candidate: RAS switch-region structural label recovery (C1) | — | `BIOLOGICAL_VALIDATION_DECISION_MEMO_2026-10-13.md` |
| D10 | 2026-10-14 | C1 blocked and paused pending Rascore label-data terms; KRAS-only preferred candidate population; design unfrozen; no margin selected | — | `C1_DECISION_NOTE_2026-10-13.md`, `C1_EVALUATION_DESIGN_DRAFT.md` |
| D11 | 2026-10-22 | Programme adopted as a proposal: near-term mission, evidence levels A–E, Phases 0–5, portfolio, verification and adaptation rules, concurrency limits | — | Roadmap §1, §12 |
| D12 | 2026-10-22 | MI double regularization recorded as a bounded engineering issue (M7); fixing it would not validate allostery | A (pending) | `KNOWN_ISSUES.md` M7 |
| D13 | 2026-10-22 | ABL recorded as a candidate subsequent application for Phase 3, not a replacement for KRAS; data access not assumed | — | Roadmap §12.3 |
| D14 | 2026-10-22 | **Owner accepted D11–D13 with corrections:** C1 stays blocked (not closed) and resumable while permission is unresolved; any alternative scientific evaluation must pass its own data-access, feasibility, protocol and budget gates; Mpro backbone analysis is separate from pocket analysis, and only the latter needs cavity observables; ensemble comparison may use any justified system; state-detection research starts only when a selected application needs discrete states | — | Roadmap §12.2, §12.3, §12.7 |
| D15 | 2026-10-22 | Owner authorized the bounded MI regularization repair (production code, tests, documentation). C1 stays paused; ABL/EN feasibility not started | — | `KNOWN_ISSUES.md` M7 |
| D16 | 2026-10-22 | MI regularization repaired: ε added once to marginal and joint covariances in both paths; regression tests against closed-form and eigenvalue references; outputs can change. M7 otherwise open | A (numerical consistency of the regularized estimator only) | `KNOWN_ISSUES.md` M7 |
| D17 | 2026-10-22 | MI documentation qualified (plug-in estimate of the fitted regularized Gaussian model; ridge need not reduce MI; determinant-fallback limitation disclosed); estimator unchanged | — | `KNOWN_ISSUES.md` M7 |
| D18 | 2026-10-22 | EN feasibility: neither candidate ready (non-commercial-only terms; ABL also needs a data request). ABL WT vs M472I recommended; GPCRmd μOR retained as an alternative. C1 stays blocked and resumable | — | `EN_FEASIBILITY_2026-10-22.md` |
| D19 | 2026-10-22 | EN note corrected: permissions separated by provider (Rascore unresolved, not non-commercial-only; DESRES; GPCRmd); owner authorization ≠ provider permission; GPCRmd frame counts and exact-permutation limit corrected. ABL request drafted, not sent; evaluation checklist added. ABL conditional; C1 blocked and resumable | — | `EN_FEASIBILITY_2026-10-22.md`, `ABL_DESRES_DATA_REQUEST_DRAFT.md` |
| D20 | 2026-10-22 | Owner authorized implementation of a reusable Cα input adapter over MDAnalysis (`python/topofold_inputs`), optional and outside the core package dependencies. Engineering only; C1 and ABL keep their access and protocol gates | — | `python/topofold_inputs/README.md` |
| D21 | 2026-10-23 | Accepted sequence: **A** reliable Cα pipeline, then **C** one bounded development exercise on permitted local data, then **B** full backbone only when a selected comparison requires it. A completed (frame-identity contract, conservative insertion-code continuity, unit enforcement, core κ/τ with validity reasons). C completed as DEV-CA-1 on local AdK (CC BY 4.0): all pre-specified checks pass. Engineering readiness only; no scientific validation. C1 and ABL keep their separate gates | Engineering; development evidence | `python/topofold_inputs/README.md`, `benchmarks/audit/dev_ca_pipeline/` |
| D22 | 2026-10-23 | Provenance correction: DEV-CA-1 `core.binary_sha256` hashed the Python wrapper, not the native extension; historical native binary recorded as unknown (overwritten, no retained wheel). Future exports record wrapper and native extension separately, or `unresolved`. Measurements unchanged | — | `benchmarks/audit/dev_ca_pipeline/PROVENANCE_CORRECTION_2026-10-23.md` |

## Rejected or withdrawn claims (kept as history)

These claims are not current. Their records are kept unchanged where they are:
- v0.8.3 benchmarks presented as real MD, estimator artefacts reported as free energies, hard-coded statistics and
  withdrawn benchmarks (PROTAC): `KNOWN_ISSUES.md` C1–C5 and M8.
- Earlier speed figures and Foldseek/PocketMiner comparisons that were never measured: `KNOWN_ISSUES.md` m9.
- Architecture records that predate the audit: `docs/adr/`, `docs/whitepaper.md`, `docs/DEVELOPMENT_HISTORY.md`. The
  roadmap and `KNOWN_ISSUES.md` take precedence over them.
- The V2b integration options: `B2_V2B_INTEGRATION_PROPOSAL.md`, marked historical and not approved (D08).
