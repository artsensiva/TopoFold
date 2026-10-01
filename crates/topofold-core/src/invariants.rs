//! Discrete geometric invariants of 3D polygonal curves.

/// Complete intrinsic differential geometry invariants of a backbone trace.
///
/// Under the Fundamental Theorem of Discrete Space Curves, a polygonal curve
/// with strictly positive segment lengths and non-zero turning angles is uniquely
/// determined up to a global rigid-body transformation (SE(3)) by:
/// 1. Segment bond lengths: `l_i` (N - 1 values)
/// 2. Discrete turning angles (curvature): `kappa_i` (N - 2 values)
/// 3. Discrete signed dihedral angles (torsion): `tau_i` (N - 3 values)
#[derive(Debug, Clone, PartialEq)]
pub struct CurveInvariants {
    /// Segment lengths `l_i = ||r_{i+1} - r_i||` for `i = 0 .. N-2`. Length = N - 1.
    pub segment_lengths: Vec<f64>,

    /// Discrete curvature (turning angles) `kappa_i` in radians `[0, pi]`, for `i = 1 .. N-2`.
    /// Length = N - 2.
    pub curvatures: Vec<f64>,

    /// Discrete torsion (signed dihedral angles) `tau_i` in radians `(-pi, pi]`, for `i = 1 .. N-3`.
    /// Length = N - 3.
    pub torsions: Vec<f64>,
}

impl CurveInvariants {
    /// Total integrated curvature along the backbone trace (sum of turning angles).
    #[must_use]
    pub fn total_curvature(&self) -> f64 {
        self.curvatures.iter().sum()
    }

    /// Average bond length along the backbone trace.
    #[must_use]
    pub fn mean_bond_length(&self) -> f64 {
        if self.segment_lengths.is_empty() {
            0.0
        } else {
            self.segment_lengths.iter().sum::<f64>() / self.segment_lengths.len() as f64
        }
    }
}
