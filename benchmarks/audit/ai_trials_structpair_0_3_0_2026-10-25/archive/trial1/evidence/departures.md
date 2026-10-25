# Documentation departures and attempts
- All files reside in the owner-authorized separate trial directory. Clone command additionally used --single-branch and checkout --detach at the exact requested commit.
- Early cat AGENTS.md failed because clone was still running. This was a tester timing mistake; no scientific run attempted.
- curl in the default sandbox returned exit 6 (Could not resolve host: github.com). Authorized elevated networking succeeded. Environmental constraint, not a TopoFold failure.
- Python environment is ../env, rather than the example ../topofold-env; paths are explicitly permitted to vary in the guide.
- CARGO_TARGET_DIR points outside clone to ../build; pip cache is inside trial. Existing shared Rust registry cache was used (fresh Python environment and freshly compiled target, not a fresh Rust installation). System Python initially lacked TopoFold, numpy, MDAnalysis, matplotlib, pandas, Streamlit, Altair, Plotly and maturin. System pip inventory retained.
- Build automatically created ignored Cargo.lock; no tracked source was changed. Temporary wheel was made by maturin in /tmp; this was automatic build behavior.
- To satisfy the owner's separate-input requirement, copied the guide KRAS fetch_inputs.py and its KRAS spec unchanged to smoke-assets/ and smoke-assets/specs/. Help shows no output-directory flag. The unchanged downloader uses its parent directory and saved all four PDBs under trial/kras_struct_*/data_pdb; spec relative paths still resolve correctly. This is self-devised path relocation, not external assistance. No shipped results read or overwritten.
- run.py is a trial-only command/log wrapper. It supplies PATH/VIRTUAL_ENV, PIP_CACHE_DIR, CARGO_TARGET_DIR, PYTHONDONTWRITEBYTECODE and records stdout/stderr, exit code, command, working directory and wall duration. Commands and origins are in logs/commands.jsonl.
- audit_inputs.py is a trial-only independent PDB record/sequence audit, not a TopoFold source modification.
