# TopoFold Scientific Audit Artifacts

This directory contains scientific validation audit trails for TopoFold core methods.

**Purpose:** Preserve the evidence chain for methodological decisions, including failed hypotheses, implementation bugs, and superseded conclusions.

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
