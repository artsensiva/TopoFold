# AI tester 1: external AI technical trial of checkpoint 6ce7b49 (2026-10-23/24)

Archived evidence of an external AI technical trial. **This is not human validation, not a human usability assessment
and not biological validation.** Pilot log entry: `docs/pilot/PILOT_LOG.md`, entry 4.

## What is here

The files are unchanged copies of the tester's trial directory, verified as identical after copying (85 of 85).
Relative links in `REPORT.md` work from this folder.

| Path | Content |
|---|---|
| `REPORT.md` | The tester's report. Identical to `docs/REPORT-ChatGPT-Independent-tester-24-10-2026.md`, as supplied by the owner |
| `run.py` | The tester's command wrapper: isolated environment variables, per-command logs and timing |
| `evidence/` | Commands (`COMMANDS.md`), per-command logs and JSON records (exit codes, timing, peak RSS), versions (`trial-metadata.json`, `13-dependency-freeze.log`), setup notes, validator scripts and results (including the preserved first attempt), browser observations and 8 screenshots, timing and interruption record, service shutdown, final tracked-file audit, and the tester's evidence manifest `evidence-sha256.json` |
| `specs/explicit_55_72.json` | The additional mapping prespecified by the tester |
| `results/kras2/`, `results/explicit_55_72/` | The tester's structpair outputs, without the `work/` intermediates |
| `fix_verification/` | **Not trial evidence.** The implementer's real-browser checks of the selection finding: on checkpoint `6ce7b49` (`before_…`) and after the fix (`after_…`), made with `apps/structpair_viewer/browser_check.py` |

## Verification performed on receipt (2026-10-24)

- **Tester's evidence manifest** (`evidence/evidence-sha256.json`, created 07:36:20 UTC): 68 of 68 listed files match;
  none is missing or modified.
- **Not covered by that manifest:** `run.py`, `specs/` and `results/`.
  - Each result directory was checked against its own structpair `manifest.json` instead: 6 of 6 outputs match in
    each.
  - `results/kras2` `windows.csv`, `residues.csv`, `exclusions.csv` and `summary.json` are byte-identical to
    `benchmarks/audit/structpair_demo/results/kras2/`.

## Not archived

These stay in the tester's trial directory:
`/home/artem/Documents/Codex/2026-10-24/act-as-an-independent-computational-software/outputs/topofold-ai-tester-1-trial/`

| Item | Size (du) |
|---|---|
| `.venv/` (trial Python environment) | 829 MB |
| `target/` (Cargo build products) | 195 MB |
| `cache/` (pip and matplotlib caches) | 213 MB |
| `checkout/` (fresh clone, including the downloaded structures) | 102 MB |
| `cargo/` (downloaded crates) | 52 MB |
| `tmp/` | 4 KB |
| `results/*/work/` (structpair intermediates) | 308 KB |

The owner's chat transcript is `docs/REPORT-ChatGPT-Independent-tester-24-10-2026-chat.md` (SHA-256
`01bb4c50…ba1d7`, 276 KB with embedded images). It is kept untracked and unchanged, and was consulted only selectively.

## Confirmation of commit 049a47f (`confirmation/`)

Bounded AI confirmation of the selection fix and the setup text. Pilot log entry 5.

| Path | Content |
|---|---|
| `confirmation/REPORT.md` | The tester's confirmation report, unchanged. The owner's formatted copy, `docs/REPORT-ChatGPT-Independent-tester-24-10-2026_2.md`, differs only in table-separator formatting and 7 appended image links to absolute local paths (its words are otherwise identical); it is kept untracked and unchanged. Use the relative links in the archived copy |
| `confirmation/run.py`, `confirmation/evidence/` | Commands (`COMMANDS.md`), per-command JSON and logs, step outcomes (`16-step-outcomes.json`), timing (`17-timing.json`), errors and limitations, final checks, the tester's manifest (`evidence-sha256.json`), and 7 screenshots with accessibility records (`evidence/browser/`) |
| `confirmation/results/kras1/`, `confirmation/results/kras2/` | Outputs without `work/` |

**Verification on receipt:**
- the tester's manifest: 47 of 47 entries match;
- 63 of 63 archived copies are identical to the originals;
- each result directory matches its own `manifest.json` (6 of 6 files each), and its 4 core outputs are byte-identical
  to `benchmarks/audit/structpair_demo/results/`;
- all 34 relative links in the archived reports and in this README resolve.

**Not archived:** they remain in the confirmation directory,
`…/topofold-ai-tester-1-trial/confirmation-049a47f/`:
- `topofold-env/` (829 MB);
- `target/` (195 MB);
- `cache/` (213 MB);
- `checkout/` (106 MB);
- `cargo/` (52 MB);
- `tmp/`;
- `results/*/work/`.
