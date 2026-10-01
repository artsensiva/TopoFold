//! Error types for geometric operations on protein backbone traces.

use thiserror::Error;

/// Errors arising during discrete differential geometry computations.
#[derive(Debug, Error, PartialEq, Clone)]
pub enum GeometryError {
    /// The trace does not have enough points to compute the requested differential primitive.
    #[error("Chain length too short: requires at least {required} C-alpha atoms, found {actual}")]
    InsufficientPoints {
        /// Minimum number of C-alpha atoms required.
        required: usize,
        /// Actual number of atoms provided.
        actual: usize,
    },

    /// Two consecutive C-alpha atoms are virtually coincident.
    #[error("Degenerate segment between C-alpha {index} and {next_index}: distance {length:.6} Å is below threshold")]
    DegenerateSegment {
        /// Index of the first atom.
        index: usize,
        /// Index of the second atom.
        next_index: usize,
        /// Measured Euclidean distance in Ångströms.
        length: f64,
    },

    /// Three consecutive C-alpha atoms are collinear, making the osculating plane ill-defined.
    #[error("Collinear segments at C-alpha vertex {vertex_index}: osculating plane is ill-defined (turning angle: {turning_angle:.6} rad)")]
    CollinearSegments {
        /// Index of the internal vertex where collinearity occurs.
        vertex_index: usize,
        /// Turning angle in radians (near 0 or near pi).
        turning_angle: f64,
    },
}
