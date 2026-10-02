# CLAUDE.md

Project instructions for Claude Code. The rules shared by all AI assistants are in AGENTS.md and are imported here:

@AGENTS.md

## Claude-specific notes
- Start each session by reading `docs/KNOWN_ISSUES.md` and asking which issue to work on.
- Use plan mode for any change that touches more than one crate or changes a public Python API.
- After changing Rust code, run `cargo test --workspace`; after changing bindings, rebuild with `maturin develop --release` and run `pytest tests/`.
