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

/// Internal accumulator for multi-atom residue parsing.
struct ResidueAccumulator {
    chain_id: char,
    seq_id: i32,
    i_code: char,
    ca_coord: Option<Point3<f64>>,
    cb_coord: Option<Point3<f64>>,
    n_coord: Option<Point3<f64>>,
    c_coord: Option<Point3<f64>>,
    meta: Option<ResidueMeta>,
}

/// Parses C-alpha and C-beta atoms from any reader implementing `BufRead`.
///
/// For residues with a C-beta atom, extracts both coordinates.
/// For Glycine (which lacks C-beta), calculates the deterministic pseudo-Cbeta coordinate
/// from the backbone bisector:
/// $$\mathbf{v}_{\text{bisect}} = \operatorname{normalized}((\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_N) + (\mathbf{r}_{\text{C}\alpha} - \mathbf{r}_C))$$
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
    let mut residues: Vec<ResidueAccumulator> = Vec::new();
    let mut res_map: std::collections::HashMap<(char, i32, char), usize> = std::collections::HashMap::new();

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
        if atom_name != "CA" && atom_name != "CB" && atom_name != "N" && atom_name != "C" {
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

        // Insertion code: column 27 (0-indexed 26)
        let i_code = line.chars().nth(26).unwrap_or(' ');

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

        let key = (chain_id, seq_id, i_code);
        let current_match = residues
            .last()
            .map(|r| (r.chain_id, r.seq_id, r.i_code) == key)
            .unwrap_or(false);

        let res_idx = if current_match {
            residues.len() - 1
        } else if let Some(&idx) = res_map.get(&key) {
            idx
        } else {
            let idx = residues.len();
            res_map.insert(key, idx);
            residues.push(ResidueAccumulator {
                chain_id,
                seq_id,
                i_code,
                ca_coord: None,
                cb_coord: None,
                n_coord: None,
                c_coord: None,
                meta: None,
            });
            idx
        };

        let pt = Point3::new(x, y, z);
        match atom_name {
            "CA" => {
                residues[res_idx].ca_coord = Some(pt);
                residues[res_idx].meta = Some(ResidueMeta {
                    name: res_name,
                    seq_id,
                    chain_id,
                    b_factor,
                });
            }
            "CB" => {
                residues[res_idx].cb_coord = Some(pt);
            }
            "N" => {
                residues[res_idx].n_coord = Some(pt);
            }
            "C" => {
                residues[res_idx].c_coord = Some(pt);
            }
            _ => {}
        }
    }

    let mut ca_coords = Vec::new();
    let mut cb_coords = Vec::new();
    let mut metadata = Vec::new();

    // First pass: collect residues that have CA
    for r in residues {
        if let (Some(ca), Some(meta)) = (r.ca_coord, r.meta) {
            let cb = if let Some(cb_pt) = r.cb_coord {
                cb_pt
            } else if let (Some(n_pt), Some(c_pt)) = (r.n_coord, r.c_coord) {
                crate::ribbon::compute_glycine_pseudo_cbeta(ca, n_pt, c_pt)
            } else {
                // Placeholder to be filled in second pass using backbone bisector
                ca
            };
            ca_coords.push(ca);
            cb_coords.push(cb);
            metadata.push(meta);
        }
    }

    if ca_coords.is_empty() {
        return Err(PdbError::EmptyTrace(target_chain));
    }

    // Second pass: fill in any remaining CB placeholders using backbone bisectors
    for i in 0..cb_coords.len() {
        if cb_coords[i] == ca_coords[i] {
            let prev_ca = if i > 0 {
                ca_coords[i - 1]
            } else if ca_coords.len() > 1 {
                ca_coords[i] + (ca_coords[i] - ca_coords[i + 1])
            } else {
                ca_coords[i] + nalgebra::Vector3::new(0.0, 0.0, 1.0)
            };
            let next_ca = if i + 1 < ca_coords.len() {
                ca_coords[i + 1]
            } else if ca_coords.len() > 1 {
                ca_coords[i] + (ca_coords[i] - ca_coords[i - 1])
            } else {
                ca_coords[i] + nalgebra::Vector3::new(0.0, 0.0, -1.0)
            };
            cb_coords[i] = crate::ribbon::compute_pseudo_cbeta_from_ca(prev_ca, ca_coords[i], next_ca);
        }
    }

    let trace = BackboneTrace::with_cbeta(ca_coords, cb_coords);
    Ok((trace, metadata))
}

/// Parses a PDB stream into an oriented [`crate::ribbon::RibbonTrace`].
pub fn parse_pdb_ribbon<R: BufRead>(
    reader: R,
    target_chain: Option<char>,
) -> Result<(crate::ribbon::RibbonTrace, Vec<ResidueMeta>), PdbError> {
    let (trace, metas) = parse_pdb_ca(reader, target_chain)?;
    let ribbon = crate::ribbon::RibbonTrace::new(
        trace.coordinates().to_vec(),
        trace.cb_coordinates().unwrap().to_vec(),
    )
    .map_err(|e| PdbError::IoError {
        line: 0,
        message: e.to_string(),
    })?;
    Ok((ribbon, metas))
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

    #[test]
    fn test_parse_pdb_ribbon_with_glycine() {
        let pdb_text = "\
ATOM      1  N   ALA A   1       0.000   1.000   0.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       0.000   0.000   0.000  1.00 10.00           C
ATOM      3  C   ALA A   1       1.000   0.000   0.000  1.00 10.00           C
ATOM      4  CB  ALA A   1       0.000   0.000   1.520  1.00 10.00           C
ATOM      5  N   GLY A   2       2.000   0.000   0.000  1.00 10.00           N
ATOM      6  CA  GLY A   2       3.000   0.000   0.000  1.00 10.00           C
ATOM      7  C   GLY A   2       4.000   0.000   0.000  1.00 10.00           C
";
        let (ribbon, metas) = parse_pdb_ribbon(pdb_text.as_bytes(), Some('A')).expect("Parse ribbon");
        assert_eq!(ribbon.len(), 2);
        assert_eq!(metas[0].name, "ALA");
        assert_eq!(metas[1].name, "GLY");

        // Residue 1 has real CB at z = 1.52
        assert_eq!(ribbon.cb_coords()[0], Point3::new(0.0, 0.0, 1.52));

        // Residue 2 (GLY) lacks CB, pseudo-CB calculated deterministically without NaN
        let cb_gly = ribbon.cb_coords()[1];
        assert!(!cb_gly.x.is_nan() && !cb_gly.y.is_nan() && !cb_gly.z.is_nan());
        let dist = (cb_gly - ribbon.ca_coords()[1]).norm();
        assert!((dist - crate::ribbon::STANDARD_CA_CB_BOND_LENGTH).abs() < 1e-9);
    }
}
