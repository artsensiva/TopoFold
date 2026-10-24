---
name: scientific-audit
description: TopoFold scientific-audit checklist. Use before freezing or running any deciding analysis, when reviewing a result, protocol or claim, or before writing a final report. Covers freeze-before-run, evidence classes, preregistered metrics and stop rules, provenance hashes, corruption tests, frozen-artifact preservation, negative outcomes, claim boundaries and the required report format.
---

# TopoFold scientific audit

`AGENTS.md` and `CLAUDE.md` are the authoritative rules and are always loaded; this skill does not repeat or override
them. It turns them into an ordered checklist for one analysis or claim. Work through the sections that apply and write
the outcome of each check into the analysis record (spec, results JSON, README), not only into chat.

Usage: `/scientific-audit <path to a run_spec.json, results directory, protocol draft or claim text>`.
If no target is given, ask which artifact to audit. Never start a deciding run from inside this skill.

## 1. Identify the target and its status

- What is being audited: protocol draft, frozen protocol, run, result, report text, or claim.
- Find the freeze record (`results/*_freeze.json` or equivalent). Record its UTC time, spec SHA-256 and code SHA-256s.
- Status is one of: draft (not frozen) · frozen, not run · run · reported. A draft must say "not frozen" wherever it is
  cited.

## 2. Freeze before any deciding calculation

- Before the deciding run: the question, inputs (with SHA-256), mapping, exclusions, metrics, thresholds, outcome
  categories, stop rules and permitted interpretation are written in a spec file, and the spec and code hashes are
  recorded (pattern: `benchmarks/audit/*/run_*.py freeze` → `results/*_freeze.json`).
- The runner refuses to run if any hash differs (`STOP: spec or code differs from the freeze`).
- Checks done before the freeze are limited to metadata, availability, mapping and reference checks of helper functions
  on SYNTHETIC inputs. Nothing that computes the deciding quantity on real data.
- Verify now: recompute the spec and code hashes and compare with the freeze record. Report the match explicitly.

## 3. Evidence class (label every number)

| Class | Meaning | May support |
|---|---|---|
| SYNTHETIC | generated data with known truth | method behaviour, calibration; must say SYNTHETIC in names and titles |
| Development | real data used while building or tuning | engineering readiness; never a performance claim |
| Held-out | real data untouched until a frozen evaluation | the frozen question only |
| Public-data observation | facts read from primary records (papers, PDB, repositories) | context; cite the exact source |
| Biological validation | agreement with an independent experimental endpoint under a frozen protocol | only the stated endpoint and conditions |

Check that the report never upgrades a class (development described as validation, SYNTHETIC as real, a draft as
frozen, an AI trial as human validation).

## 4. Preregistered question, metrics, thresholds and stop rules

- One primary question; secondary analyses labelled as such.
- Metrics defined with units; circular quantities handled circularly.
- Thresholds and outcome categories fixed before data are seen, including "not supportive" and "inconclusive".
- Stop rules listed (missing permission, failed mapping, insufficient replicates, failed validation).
- Equal information for every compared method (AGENTS rule 6); the strongest reasonable baseline is named.
- No p-value, interval or "significance" unless a calibrated test is part of the frozen protocol (decision D33 conditions;
  AGENTS rule 4). Otherwise report descriptive ranges and say why.
- Labels or states never come from TopoFold's own output (AGENTS rule 5).

## 5. Provenance and hashes

- Inputs: source, version or accession, licence, SHA-256 of every file used.
- Code: SHA-256 of runner and library modules in the freeze record; software versions (Python, NumPy, MDAnalysis, core
  native-extension hash where relevant).
- Outputs: SHA-256 in a manifest; results numbers live in machine-readable files (`results/*.json`, CSV) and text reads
  from them (AGENTS rule 3).
- Downloads and CPU time recorded against the stated budget.

## 6. Corruption or mutation tests (when justified)

Use when a result could arise from a pipeline artefact, mapping error or leakage. Choose before running, for example:
shuffled residue mapping, swapped condition labels, a deliberately broken or truncated input that must be refused,
tampered output hashes that must be detected, a mirror image for superposition code. State what each test should show.
Skip with a reason when not justified.

## 7. Preserve frozen artifacts

- Never edit a frozen spec, runner, result or archived evidence. A needed change is a dated amendment or a new versioned
  protocol, recorded before any re-scoring.
- Before committing, check `git status` / `git diff --stat` shows no modification to frozen paths; re-verify recorded
  hashes when an artifact is copied or archived.

## 8. Outcomes, including negative and inconclusive

- Report the outcome category exactly as defined in the frozen spec, with the measured numbers, sample sizes and
  uncertainty (or the reason none is given).
- Report negative and inconclusive results with the same prominence; no re-scoring on the same data after seeing them.
- Separate post-hoc observations from preregistered results.

## 9. Claim boundaries

- State what the result does not show (populations, kinetics, function, mechanism, causation, clinical relevance,
  generality beyond the tested systems) as applicable.
- Name confounders (ligands, mutations, constructs, crystal forms, replicate dependence, starting structures).
- Language rules: AGENTS rule 10.

## 10. Report format (required for every TopoFold task)

1. **Plain-language overview** – what was done, learned, why it matters, what is unresolved; measured results separated
   from hypotheses and interpretations.
2. **Technical and scientific report** – work per package; question, methods, data, decision criteria; results with
   units, n and uncertainty; limitations and claim boundaries; files changed; checks with exact counts (never "all
   checks passed"); resources and downloads; Git commit, working tree and remote status; provenance and unexpected
   findings.
3. **Owner decisions and actions** – only what needs the owner, why, and what stays blocked.
4. **Recommended next step** – one task, with alternatives (value, cost, risk, prerequisites) and the reason.

Also state explicitly whether anything was sent externally.

## Output of this skill

Write an audit note (in the analysis directory or the reply) with one line per section above:
`<section>: pass | fail | not applicable – <evidence: path, hash or quoted line>`.
A single fail in sections 2, 3, 5 or 7 blocks reporting the result as frozen evidence until resolved.
