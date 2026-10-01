//! In-memory streaming parser for PDB C-alpha protein backbones.
//!
//! Extracts only the C-alpha coordinates and residue metadata without the overhead
//! of parsing all side-chain heavy atoms and solvent molecules.

#![forbid(unsafe_code)]

use nalgebra::Point3;
use std::io::BufRead;
use thiserror::Error;

use crate::types::BackboneTrace;

/// Errors arising during PDB format parsing.
#[derive(Debug, Error, PartialEq, Clone)]
pub enum PdbError {
    /// An I/O error occurred while reading the stream.
    #[error("I/O error at line {line}: {message}")]
    IoError {
        /// Line number where the error occurred.
        line: usize,
        /// Detail error message.
        message: String,
    },

    /// An invalid coordinate format was encountered.
    #[error("Invalid coordinate at line {line}, column {column}: {value}")]
    InvalidCoordinate {
        /// Line number where the error occurred.
        line: usize,
        /// Column name ("X", "Y", or "Z").
        column: &'static str,
        /// Raw value that failed parsing.
        value: String,
    },

    /// An invalid residue sequence number was encountered.
    #[error("Invalid residue sequence number at line {line}: {value}")]
    InvalidResidueNumber {
        /// Line number where the error occurred.
        line: usize,
        /// Raw string that failed parsing.
        value: String,
    },

    /// The PDB file contained no C-alpha atoms for the specified chain.
    #[error("No C-alpha atoms found for chain '{0:?}'")]
    EmptyTrace(Option<char>),
}

/// Metadata associated with an individual C-alpha residue in the backbone trace.
#[derive(Debug, Clone, PartialEq)]
pub struct ResidueMeta {
    /// 3-letter amino acid code (e.g. "ALA", "GLY", "TRP").
    pub name: String,
    /// Residue sequence identifier (PDB resSeq).
    pub seq_id: i32,
    /// Chain identifier (PDB chainID).
    pub chain_id: char,
    /// Isotropic temperature factor (B-factor in Å^2).
    pub b_factor: f64,
}

/// Parses C-alpha atoms from any reader implementing `BufRead`.
///
/// # Arguments
/// - `reader`: Buffered input stream of PDB file content.
/// - `target_chain`: Optional chain ID to filter (e.g. `Some('A')`). If `None`, parses all chains sequentially.
///
/// # Errors
/// Returns [`PdbError`] if parsing fails or no C-alpha atoms are found.
pub fn parse_pdb_ca<R: BufRead>(
    reader: R,
    target_chain: Option<char>,
) -> Result<(BackboneTrace, Vec<ResidueMeta>), PdbError> {
    let mut coords = Vec::new();
    let mut metadata = Vec::new();

    for (line_idx, line_res) in reader.lines().enumerate() {
        let line_num = line_idx + 1;
        let line = line_res.map_err(|e| PdbError::IoError {
            line: line_num,
            message: e.to_string(),
        })?;

        // Standard PDB coordinate records require at least 54 columns
        if line.len() < 54 {
            continue;
        }

        let record = &line[0..6];
        if record != "ATOM  " && record != "HETATM" {
            continue;
        }

        // Atom name: columns 13-16
        let atom_name = line[12..16].trim();
        if atom_name != "CA" {
            continue;
        }

        // Alternate location indicator: column 17 (0-indexed 16)
        let alt_loc = line.chars().nth(16).unwrap_or(' ');
        if alt_loc != ' ' && alt_loc != 'A' && alt_loc != '1' {
            continue;
        }

        // Chain identifier: column 22 (0-indexed 21)
        let chain_id = line.chars().nth(21).unwrap_or(' ');
        if let Some(target) = target_chain {
            if chain_id != target {
                continue;
            }
        }

        // Residue name: columns 18-20
        let res_name = line[17..20].trim().to_string();

        // Residue sequence number: columns 23-26
        let seq_id_str = line[22..26].trim();
        let seq_id: i32 = seq_id_str.parse().map_err(|_| PdbError::InvalidResidueNumber {
            line: line_num,
            value: seq_id_str.to_string(),
        })?;

        // Orthogonal coordinates X, Y, Z: columns 31-38, 39-46, 47-54
        let x_str = line[30..38].trim();
        let x: f64 = x_str.parse().map_err(|_| PdbError::InvalidCoordinate {
            line: line_num,
            column: "X",
            value: x_str.to_string(),
        })?;

        let y_str = line[38..46].trim();
        let y: f64 = y_str.parse().map_err(|_| PdbError::InvalidCoordinate {
            line: line_num,
            column: "Y",
            value: y_str.to_string(),
        })?;

        let z_str = line[46..54].trim();
        let z: f64 = z_str.parse().map_err(|_| PdbError::InvalidCoordinate {
            line: line_num,
            column: "Z",
            value: z_str.to_string(),
        })?;

        // Temperature factor: columns 61-66
        let b_factor = if line.len() >= 66 {
            line[60..66].trim().parse().unwrap_or(0.0)
        } else {
            0.0
        };

        coords.push(Point3::new(x, y, z));
        metadata.push(ResidueMeta {
            name: res_name,
            seq_id,
            chain_id,
            b_factor,
        });
    }

    if coords.is_empty() {
        return Err(PdbError::EmptyTrace(target_chain));
    }

    Ok((BackboneTrace::new(coords), metadata))
}

/// Convenience function to parse a PDB string in memory.
pub fn parse_pdb_str(
    pdb_content: &str,
    target_chain: Option<char>,
) -> Result<(BackboneTrace, Vec<ResidueMeta>), PdbError> {
    parse_pdb_ca(pdb_content.as_bytes(), target_chain)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_PDB: &str = "\
HEADER    SAMPLE PROTEIN
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N  
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C  
ATOM      3  C   THR A   1      15.685  12.755   5.133  1.00  9.19           C  
ATOM      4  O   THR A   1      15.268  13.825   5.594  1.00  9.85           O  
ATOM      8  CA  CYS A   2      17.200  10.165   7.382  1.00  7.50           C  
ATOM     15  CA  PRO A   3      14.150   8.220   8.910  1.00  8.20           C  
ATOM     22  CA  ILE B   4      11.520   6.100  11.230  1.00  6.90           C  
";

    #[test]
    fn test_parse_pdb_chain_a() {
        let (trace, metas) = parse_pdb_str(SAMPLE_PDB, Some('A')).expect("Parse chain A");
        assert_eq!(trace.len(), 3);
        assert_eq!(metas.len(), 3);

        assert_eq!(metas[0].name, "THR");
        assert_eq!(metas[0].seq_id, 1);
        assert_eq!(metas[0].chain_id, 'A');
        assert_eq!(metas[0].b_factor, 10.80);

        assert_eq!(metas[1].name, "CYS");
        assert_eq!(metas[2].name, "PRO");

        assert_eq!(trace.coordinates()[0], Point3::new(16.967, 12.784, 4.338));
    }

    #[test]
    fn test_parse_pdb_all_chains() {
        let (trace, metas) = parse_pdb_str(SAMPLE_PDB, None).expect("Parse all chains");
        assert_eq!(trace.len(), 4);
        assert_eq!(metas.len(), 4);
        assert_eq!(metas[3].chain_id, 'B');
    }
}
