# Claude Code setup for TopoFold (optional tooling)

`AGENTS.md` and `CLAUDE.md` remain the authoritative project rules and are always loaded. Everything below is optional:
nothing here changes code, data, protocols or results.

## Project skills (version-controlled)

| Skill | Purpose | Invoke |
|---|---|---|
| `scientific-audit` (`.claude/skills/scientific-audit/SKILL.md`) | Checklist for freezing protocols, evidence classes, preregistered metrics and stop rules, provenance hashes, corruption tests, frozen-artifact preservation, negative outcomes, claim boundaries and the required report format | `/scientific-audit <spec, results dir, draft or claim>` |
| `dataset-intake` (`.claude/skills/dataset-intake/`) | Intake checklist and machine-readable record (`intake_record.template.json`) ending in go / conditional-go / no-go | `/dataset-intake <dataset name, URL, DOI or path>` |

- **Loading:** Claude Code loads a skill's short description at session start and the full text only when the skill is
  invoked. Claude may also apply a skill when its description matches the task.
- **Context cost:** measured here as about 96 and 117 tokens of description per session, and about 1,700 and 1,400
  tokens when invoked.
- **Make a skill manual-only:** add `disable-model-invocation: true` to its frontmatter. It then leaves the listing and
  runs only through `/name`.
- **Turn a skill off locally without editing it:** add `"skillOverrides": {"scientific-audit": "off"}` to
  `.claude/settings.local.json`. Deleting the directory removes it.
- **Git:** only these two skill directories are tracked. `.gitignore` keeps the rest of `.claude/` local, including
  `settings.local.json` and any personal skills.

## Code intelligence (per machine, not version-controlled)

These are Anthropic's official LSP plugins from the `claude-plugins-official` marketplace. Each only starts a language
server binary that you install yourself; neither adds context at session start (`claude plugin details` reports ~0
always-on tokens). Claude gets an `LSP` tool for navigation (definitions, references, symbols) and diagnostics after it
edits a file.

| Plugin | Binary | Install the binary | Enable the plugin for this repository only |
|---|---|---|---|
| `rust-analyzer-lsp` | `rust-analyzer` | `rustup component add rust-analyzer` | `claude plugin install rust-analyzer-lsp@claude-plugins-official --scope local` |
| `pyright-lsp` | `pyright-langserver` | `npm install -g pyright` (Microsoft, MIT) | `claude plugin install pyright-lsp@claude-plugins-official --scope local` |

- **Scope:** local scope writes `enabledPlugins` to `.claude/settings.local.json`, which is git-ignored, so the change
  applies to this repository on this machine only. Language-server binaries are machine-specific, so the plugins are not
  enabled for collaborators through a committed file.
- **Activation:** a newly installed plugin becomes active after `/reload-plugins` or in a new session.
- **Disable:** `claude plugin disable <plugin>@claude-plugins-official`.
- **Remove:** `claude plugin uninstall <plugin>@claude-plugins-official --scope local`.
- **Check:** `claude plugin list` shows scope and status; the `/plugin` **Errors** tab shows a missing binary.
- **Diagnostics are advisory.** Pyright reports type-level issues; for example `python/topofold_inputs/compare.py:174`
  passes `sum(...)` of arrays where the type checker infers `int | ndarray`. Do not change scientific code only to
  silence a checker. Treat a reported issue as a normal code change with its own tests.

## Verification record (2026-10-24)

- **Rust**, in the setup session:
  - `documentSymbol` on `crates/topofold-core/src/circular.rs` returned 15 symbols.
  - `goToDefinition` from `bimodality.rs:37` resolved to `circular.rs:73:8`.
  - `rust-analyzer diagnostics` (read-only) on `topofold-core` reported no errors.
- **Python**, in a fresh non-interactive session limited to the LSP tool:
  - `goToDefinition` from `structpair.py:191:17` resolved to `compare.py:65:5`.
  - `findReferences` returned 9 references in 3 files.
  - The Pyright CLI flagged both deliberate errors in a scratch file outside the repository.
- **Skills:** both were invoked explicitly on read-only targets. No repository file was changed by the trials.

## Not configured here

No MCP servers, GitHub or Playwright integrations, or plugin bundles are added by this setup.
