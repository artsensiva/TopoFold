//! Multi-model PDB trajectory streaming parser.

#![forbid(unsafe_code)]

use std::io::BufRead;
use topofold_core::pdb::parse_pdb_ca;
use topofold_core::types::BackboneTrace;

use crate::error::TrajectoryError;

/// Reads a multi-model PDB trajectory into a sequence of C-alpha [`BackboneTrace`] frames.
///
/// Handles files with `MODEL` / `ENDMDL` blocks as well as single-model PDBs.
pub fn read_pdb_trajectory<R: BufRead>(
    mut reader: R,
    chain: Option<char>,
) -> Result<Vec<BackboneTrace>, TrajectoryError> {
    let mut lines = Vec::new();
    let mut traces = Vec::new();
    let mut in_model = false;
    let mut model_found = false;

    let mut line_buf = String::new();
    while reader.read_line(&mut line_buf)? > 0 {
        let trimmed = line_buf.trim();
        if trimmed.starts_with("MODEL") {
            model_found = true;
            in_model = true;
            lines.clear();
        } else if trimmed.starts_with("ENDMDL") {
            if in_model && !lines.is_empty() {
                let model_str = lines.join("");
                let (trace, _) = parse_pdb_ca(model_str.as_bytes(), chain)?;
                traces.push(trace);
                lines.clear();
            }
            in_model = false;
        } else {
            lines.push(line_buf.clone());
        }
        line_buf.clear();
    }

    if !model_found && !lines.is_empty() {
        let model_str = lines.join("");
        let (trace, _) = parse_pdb_ca(model_str.as_bytes(), chain)?;
        traces.push(trace);
    }

    if traces.is_empty() {
        return Err(TrajectoryError::EmptyTrajectory);
    }

    Ok(traces)
}
