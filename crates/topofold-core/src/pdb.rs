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
    let mut res_map: std::collections::HashMap<(char, i32, char), usize> =
        std::collections::HashMap::new();

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
        let seq_id: i32 = seq_id_str
            .parse()
            .map_err(|_| PdbError::InvalidResidueNumber {
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
            cb_coords[i] =
                crate::ribbon::compute_pseudo_cbeta_from_ca(prev_ca, ca_coords[i], next_ca);
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

/// Internal accumulator for multi-atom RNA nucleotide parsing.
struct RnaNucleotideAccumulator {
    chain_id: char,
    seq_id: i32,
    i_code: char,
    res_name: String,
    p_coord: Option<Point3<f64>>,
    c4_coord: Option<Point3<f64>>,
    c1_coord: Option<Point3<f64>>,
    n_coord: Option<Point3<f64>>,
    o5_coord: Option<Point3<f64>>,
    c5_coord: Option<Point3<f64>>,
    c3_coord: Option<Point3<f64>>,
    b_factor: f64,
}

/// Parses an RNA ribbon trace from any reader implementing `BufRead`.
///
/// Recognizes RNA residues (A, G, C, U, ADE, GUA, CYT, URA) and extracts
/// Phosphorus (P), Ribose C4', Ribose C1', and the glycosidic nitrogen (N9 for purines, N1 for pyrimidines).
/// For 5'-terminal nucleotides lacking a Phosphorus atom, a deterministic pseudo-Phosphorus coordinate
/// is reconstructed from available terminal ribose atoms.
///
/// # Arguments
/// - `reader`: Buffered input stream of PDB file content.
/// - `target_chain`: Optional chain ID to filter (e.g. `Some('X')`). If `None`, parses all chains.
///
/// # Errors
/// Returns [`PdbError`] if parsing fails or no RNA nucleotides are found.
pub fn parse_pdb_rna<R: BufRead>(
    reader: R,
    target_chain: Option<char>,
) -> Result<(crate::rna::RnaRibbonTrace, Vec<ResidueMeta>), PdbError> {
    let mut nucleotides: Vec<RnaNucleotideAccumulator> = Vec::new();
    let mut nt_map: std::collections::HashMap<(char, i32, char), usize> =
        std::collections::HashMap::new();

    let rna_names = [
        "A", "G", "C", "U", "ADE", "GUA", "CYT", "URA", "DA", "DG", "DC", "DT",
    ];

    for (line_idx, line_res) in reader.lines().enumerate() {
        let line_num = line_idx + 1;
        let line = line_res.map_err(|e| PdbError::IoError {
            line: line_num,
            message: e.to_string(),
        })?;

        if line.len() < 54 {
            continue;
        }

        let record = &line[0..6];
        if record != "ATOM  " && record != "HETATM" {
            continue;
        }

        let res_name = line[17..20].trim().to_uppercase();
        if !rna_names.contains(&res_name.as_str()) {
            continue;
        }

        let chain_id = line.chars().nth(21).unwrap_or(' ');
        if let Some(target) = target_chain {
            if chain_id != target {
                continue;
            }
        }

        let alt_loc = line.chars().nth(16).unwrap_or(' ');
        if alt_loc != ' ' && alt_loc != 'A' && alt_loc != '1' {
            continue;
        }

        let seq_id_str = line[22..26].trim();
        let seq_id: i32 = seq_id_str
            .parse()
            .map_err(|_| PdbError::InvalidResidueNumber {
                line: line_num,
                value: seq_id_str.to_string(),
            })?;

        let i_code = line.chars().nth(26).unwrap_or(' ');

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

        let b_factor = if line.len() >= 66 {
            line[60..66].trim().parse().unwrap_or(0.0)
        } else {
            0.0
        };

        let key = (chain_id, seq_id, i_code);
        let current_match = nucleotides
            .last()
            .map(|r| (r.chain_id, r.seq_id, r.i_code) == key)
            .unwrap_or(false);

        let nt_idx = if current_match {
            nucleotides.len() - 1
        } else if let Some(&idx) = nt_map.get(&key) {
            idx
        } else {
            let idx = nucleotides.len();
            nt_map.insert(key, idx);
            nucleotides.push(RnaNucleotideAccumulator {
                chain_id,
                seq_id,
                i_code,
                res_name: res_name.clone(),
                p_coord: None,
                c4_coord: None,
                c1_coord: None,
                n_coord: None,
                o5_coord: None,
                c5_coord: None,
                c3_coord: None,
                b_factor,
            });
            idx
        };

        let pt = Point3::new(x, y, z);
        let atom_name = line[12..16].trim().to_uppercase();

        let is_purine = ["A", "G", "ADE", "GUA", "DA", "DG"].contains(&res_name.as_str());

        match atom_name.as_str() {
            "P" => nucleotides[nt_idx].p_coord = Some(pt),
            "C4'" | "C4*" => nucleotides[nt_idx].c4_coord = Some(pt),
            "C1'" | "C1*" => nucleotides[nt_idx].c1_coord = Some(pt),
            "O5'" | "O5*" => nucleotides[nt_idx].o5_coord = Some(pt),
            "C5'" | "C5*" => nucleotides[nt_idx].c5_coord = Some(pt),
            "C3'" | "C3*" => nucleotides[nt_idx].c3_coord = Some(pt),
            "N9" if is_purine => nucleotides[nt_idx].n_coord = Some(pt),
            "N1" => {
                if !is_purine || nucleotides[nt_idx].n_coord.is_none() {
                    nucleotides[nt_idx].n_coord = Some(pt);
                }
            }
            "N3" if !is_purine && nucleotides[nt_idx].n_coord.is_none() => {
                nucleotides[nt_idx].n_coord = Some(pt);
            }
            _ => {}
        }
    }

    let mut p_coords = Vec::new();
    let mut c4_coords = Vec::new();
    let mut c1_coords = Vec::new();
    let mut n_coords = Vec::new();
    let mut names = Vec::new();
    let mut seq_ids = Vec::new();
    let mut chain_ids = Vec::new();
    let mut metadata = Vec::new();

    for nt in nucleotides {
        if let (Some(c1), Some(n_pt)) = (nt.c1_coord, nt.n_coord) {
            let c4 = nt.c4_coord.unwrap_or(c1);

            let p = if let Some(p_pt) = nt.p_coord {
                p_pt
            } else if let (Some(o5), Some(c5)) = (nt.o5_coord, nt.c5_coord) {
                // Reconstruct pseudo-P from 5'-terminal O5' and C5'
                let diff = o5 - c5;
                let norm = diff.norm();
                if norm >= crate::geometry::GEOMETRY_EPSILON {
                    o5 + 1.60 * (diff / norm)
                } else {
                    c4 + nalgebra::Vector3::new(0.0, 0.0, 3.5)
                }
            } else if let Some(c3) = nt.c3_coord {
                let diff = c4 - c3;
                let norm = diff.norm();
                if norm >= crate::geometry::GEOMETRY_EPSILON {
                    c4 + 2.50 * (diff / norm)
                } else {
                    c4 + nalgebra::Vector3::new(0.0, 0.0, 3.5)
                }
            } else {
                c4 + nalgebra::Vector3::new(0.0, 0.0, 3.5)
            };

            p_coords.push(p);
            c4_coords.push(c4);
            c1_coords.push(c1);
            n_coords.push(n_pt);
            names.push(nt.res_name.clone());
            seq_ids.push(nt.seq_id);
            chain_ids.push(nt.chain_id);
            metadata.push(ResidueMeta {
                name: nt.res_name,
                seq_id: nt.seq_id,
                chain_id: nt.chain_id,
                b_factor: nt.b_factor,
            });
        }
    }

    if p_coords.is_empty() {
        return Err(PdbError::EmptyTrace(target_chain));
    }

    let trace = crate::rna::RnaRibbonTrace::with_meta(
        p_coords, c4_coords, c1_coords, n_coords, names, seq_ids, chain_ids,
    )
    .map_err(|e| PdbError::IoError {
        line: 0,
        message: e.to_string(),
    })?;

    Ok((trace, metadata))
}

/// Convenience function to parse an RNA PDB string in memory.
pub fn parse_pdb_rna_str(
    pdb_content: &str,
    target_chain: Option<char>,
) -> Result<(crate::rna::RnaRibbonTrace, Vec<ResidueMeta>), PdbError> {
    parse_pdb_rna(pdb_content.as_bytes(), target_chain)
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
        let (ribbon, metas) =
            parse_pdb_ribbon(pdb_text.as_bytes(), Some('A')).expect("Parse ribbon");
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

    #[test]
    fn test_parse_pdb_rna() {
        let rna_pdb = "\
ATOM      1  C4'   C X  13       0.659  -6.175  -0.764  1.00 66.28           C
ATOM      2  C1'   C X  13       0.134  -3.945  -1.350  1.00 65.21           C
ATOM      3  N1    C X  13       0.664  -2.744  -0.554  1.00 64.39           N
ATOM      4  O5'   C X  13       0.047  -6.959   1.489  1.00 65.39           O
ATOM      5  C5'   C X  13       1.106  -6.846   0.528  1.00 66.02           C
ATOM      6  P     G X  14       2.936  -7.140  -3.702  1.00 66.85           P
ATOM      7  C4'   G X  14       2.484  -6.666  -7.391  1.00 63.85           C
ATOM      8  C1'   G X  14       2.302  -4.321  -7.106  1.00 61.99           C
ATOM      9  N9    G X  14       3.125  -3.376  -6.331  1.00 59.83           N
";
        let (rna_trace, metas) = parse_pdb_rna(rna_pdb.as_bytes(), Some('X')).expect("Parse RNA");
        assert_eq!(rna_trace.len(), 2);
        assert_eq!(metas.len(), 2);
        assert_eq!(metas[0].name, "C");
        assert_eq!(metas[0].seq_id, 13);
        assert_eq!(metas[1].name, "G");
        assert_eq!(metas[1].seq_id, 14);

        // Residue 13 lacked P, pseudo-P reconstructed from O5' and C5' without NaN
        let p13 = rna_trace.p_coords()[0];
        assert!(!p13.x.is_nan() && !p13.y.is_nan() && !p13.z.is_nan());
        // Residue 14 has real P
        let p14 = rna_trace.p_coords()[1];
        assert_eq!(p14, Point3::new(2.936, -7.140, -3.702));
    }
}
