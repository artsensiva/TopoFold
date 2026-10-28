# CLAUDE.md

Project instructions for Claude Code. The rules shared by all AI assistants are in AGENTS.md and are imported here:

@AGENTS.md

## Claude-specific notes
- Start each session by reading `docs/PROJECT_STATE.md` (the canonical restart document), then `git status` and `git log -3`. While the project is **PAUSED BY OWNER**, report the status and wait for an explicit owner task; do not ask the owner to pick a backlog issue and do not begin work. After the owner lifts the pause, read `docs/KNOWN_ISSUES.md` and ask which issue to work on.
- Application memory is not the record: everything needed to restart is in the repository and `docs/PROJECT_STATE.md`. Never commit private chat histories, account information or credentials.
- Use plan mode for any change that touches more than one crate or changes a public Python API.
- After changing Rust code, run `cargo test --workspace`; after changing bindings, rebuild with `maturin develop --release` and run `pytest tests/`.
