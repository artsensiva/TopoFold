# TopoFold Scientific Audit Artifacts

This directory contains scientific validation audit trails for TopoFold core methods.

**Purpose:** Preserve the evidence chain for methodological decisions, including failed hypotheses, implementation bugs, and superseded conclusions.

---

## Workstream index (reconciled 2026-10-13)

Each directory holds harness code, frozen-protocol references, results files and a report. **Harness code is not production code** (see below). Frozen protocols live in [`../../docs/`](../../docs/).

| Directory | Workstream | Data | Entry point | Status (short) |
|---|---|---|---|---|
| [`m3_multimodality/`](m3_multimodality/) | M3 Product-vM multimodality detector audit | SYNTHETIC | [`M3_STATUS.md`](m3_multimodality/M3_STATUS.md) | Authoritative status in `M3_STATUS.md`; not production-ready |
| [`angular_characterization/`](angular_characterization/) | WP-A angular characterization of the BC torsion scan | SYNTHETIC | [`REPORT.md`](angular_characterization/REPORT.md) | Characterization of one component |
| [`ordered_representation/`](ordered_representation/) | WP-C ordered angular representation | SYNTHETIC | [`REPORT.md`](ordered_representation/REPORT.md) | Candidate representation, fixtures only |
| [`adk_feasibility/`](adk_feasibility/) | WP-B AdK real-MD feasibility (G0–G1) | Third-party MD (CC BY 4.0) | [`REPORT.md`](adk_feasibility/REPORT.md) | Data integrity and reference assessment; development data |
| [`wpd_local_retrieval/`](wpd_local_retrieval/) | WP-D local retrieval: AdK development tier and mdCATH held-out evaluation | Third-party MD derivatives (CC BY 4.0) | [`REPORT_HELDOUT.md`](wpd_local_retrieval/REPORT_HELDOUT.md), [`REPORT_ADK.md`](wpd_local_retrieval/REPORT_ADK.md) | Held-out evaluation completed; see closeout addendum |
| [`legacy_engine/`](legacy_engine/) | LE-1 legacy κ–τ engine evaluation | AdK derivative + SYNTHETIC | [`REPORT.md`](legacy_engine/REPORT.md) | Legacy Fréchet coupling deprioritized (engineering decision) |
| [`eb1/`](eb1/) | EB-1 engineering benchmark, Amendments 1–2 | Development data | [`README.md`](eb1/README.md) | EB-1/A1 insufficient evidence; A2 decision engine (projection) |
| [`b2opt/`](b2opt/) | BQ-1 B2 single-query optimization (V0/V1/V2/V2b) | Development data | [`REPORT.md`](b2opt/REPORT.md) | Acceptance not met; V2b correct and faster on tested workloads |
| [`bq2/`](bq2/) | BQ-2 scaling confirmation of unchanged V2b | Development data + SYNTHETIC_tiled (performance-only) | [`REPORT.md`](bq2/REPORT.md) | Conclusion C (insufficient evidence) |
| [`bq3/`](bq3/) | BQ-3 direct-P2 qualification | — | [`README.md`](bq3/README.md) | Stopped at memory preflight; no measurement |

Large or third-party inputs that are not in Git, and how to retrieve or regenerate them: [`../../docs/CHECKPOINT_2026-10-13_PROVENANCE.md`](../../docs/CHECKPOINT_2026-10-13_PROVENANCE.md).

---

## Audit Principles

### Why Version Scientific Audit Artifacts?

**Transparency and Reproducibility:**
- Document how conclusions were reached, including errors discovered
- Preserve evidence for methodological decisions
- Enable independent verification of validation claims
- Support publication reproducibility requirements

**NOT a Clean Success Narrative:**
- Failed hypotheses are scientifically valuable (document what doesn't work)
- Implementation bugs and their discovery are part of the evidence
- Superseded conclusions remain in the historical record with supersession markers
- The messy scientific process is intentionally visible

### Distinction: Audit vs. Production Code

| Audit Artifacts | Production Code |
|----------------|-----------------|
| **Location:** `benchmarks/audit/` | **Location:** `src/`, `python/topofold/` |
| **Purpose:** Diagnostic, validation, evidence | **Purpose:** Production implementation |
| **Quality:** Research/experimental quality | **Quality:** Production quality |
| **Promotion:** Explicit approval required | **Promotion:** Not automatic from audit |

**IMPORTANT:** Audit artifact inclusion does NOT imply production readiness.

---

## Current Workstreams

### M3: Toroidal Multimodality Detection

> **Superseded summary (marker added 2026-10-13).** The status lines below date from M3-A2f (2026-10-03). Later checkpoints (M3-A2g) changed the classification. The authoritative current status is [`m3_multimodality/M3_STATUS.md`](m3_multimodality/M3_STATUS.md). The text below is kept for the record.

**Location:** [m3_multimodality/](m3_multimodality/)

**Current Status:** **M3-A2f CLOSED**

**Classification:**
```
VALID ALGORITHMIC HEURISTIC /
OPTIMIZATION ROBUSTNESS UNRESOLVED
```

**Product von Mises detector:**
- Exact-family validation: 0.13% (3/2250)
- K2 power (q=2, ≥60°): 100% (1620/1620)
- **PRIMARY BLOCKER:** 97% EM convergence failure
- **NOT production-ready**

See [m3_multimodality/M3_STATUS.md](m3_multimodality/M3_STATUS.md)

---

## Historical Checkpoints and Supersession

### What is a "Superseded" Checkpoint?

A checkpoint whose **conclusions were invalidated** by later work, but **retained for provenance**.

**Supersession Header Format:**
```markdown
SUPERSEDED SCIENTIFIC CONCLUSION

This report is retained for provenance.

Invalidated by: [checkpoint]
Reason: [one-line explanation]
Current status: see [STATUS.md]
```

### Why Preserve Failed Hypotheses?

1. **Document negative results**
2. **Prevent repeated mistakes**
3. **Methodological transparency**
4. **Complete evidence chain**

---

## Promotion to Production

**Audit artifacts are NOT automatically production code.**

Before promoting any audit implementation to `src/` or `python/topofold/`:

1. **Scientific validation complete** (all MANDATORY requirements met)
2. **External review and approval** (documented in ADR)
3. **Code quality standards** (production-quality implementation, tests, documentation)
4. **API design** (public interface, backwards compatibility)
5. **Integration testing** (works with existing TopoFold pipeline)

---

**Repository:** TopoFold
**Branch:** scientific-rebuild
**Last Updated:** 2026-10-03
