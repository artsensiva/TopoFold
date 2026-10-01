//! Trajectory I/O error definitions.

#![forbid(unsafe_code)]

use thiserror::Error;
use topofold_core::pdb::PdbError;

/// Errors arising during trajectory I/O (DCD, PDB multi-model).
#[derive(Debug, Error)]
pub enum TrajectoryError {
    /// Underlying I/O error.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// Unexpected end of file.
    #[error("Unexpected end of file while reading trajectory stream")]
    UnexpectedEof,

    /// Invalid header format or magic number.
    #[error("Invalid trajectory header: {0}")]
    InvalidHeader(String),

    /// Record size marker does not match expected length.
    #[error("Fortran record size mismatch: expected {expected} bytes, found {found} bytes")]
    RecordSizeMismatch {
        /// Expected byte length.
        expected: usize,
        /// Actual byte length.
        found: usize,
    },

    /// Atom index out of bounds.
    #[error("Atom index {index} out of bounds (frame contains {n_atoms} atoms)")]
    IndexOutOfBounds {
        /// Requested atom index.
        index: usize,
        /// Total atoms available in frame.
        n_atoms: usize,
    },

    /// Frame atom count does not match header specification.
    #[error("Atom count mismatch: expected {expected} atoms, found {found} atoms")]
    AtomCountMismatch {
        /// Expected atom count.
        expected: usize,
        /// Actual atom count.
        found: usize,
    },

    /// PDB format parsing error.
    #[error("PDB trajectory error: {0}")]
    Pdb(#[from] PdbError),

    /// Trajectory contains no valid coordinate frames.
    #[error("Trajectory contains no valid coordinate frames")]
    EmptyTrajectory,
}
