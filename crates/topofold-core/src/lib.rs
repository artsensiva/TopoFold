//! # TopoFold Core Engine
//!
//! High-performance discrete differential geometry engine for protein backbones,
//! extracting rotation- and translation-invariant (SE(3)-invariant) curve signatures.
//!
//! ## Overview
//!
//! Traditional conformational analysis relies on Cartesian coordinates and linear PCA,
//! which requires expensive pairwise structural alignment (Kabsch algorithm) and is
//! sensitive to rotational drift.
//!
//! TopoFold models the C-alpha protein backbone as a discrete 3D space curve, representing
//! conformations via intrinsic local differential invariants:
//! - Segment lengths: `l_i`
//! - Discrete curvature (turning angles): `kappa_i \in [0, \pi]`
//! - Discrete torsion (osculating plane dihedrals): `tau_i \in (-\pi, \pi]`
//!
//! These invariants are strictly SE(3)-invariant and chiral-sensitive.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod error;
pub mod geometry;
pub mod invariants;
pub mod types;

// Re-export primary types and functions for ergonomic top-level use
pub use error::GeometryError;
pub use geometry::{
    compute_binormals, compute_curvatures, compute_tangents, compute_torsions,
    extract_curve_invariants, GEOMETRY_EPSILON,
};
pub use invariants::CurveInvariants;
pub use types::BackboneTrace;
