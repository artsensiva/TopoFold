//! Discrete Differential Geometry for Protein Backbone C-alpha Curves.
//!
//! Treats a protein backbone as a 3D piecewise-linear space curve and extracts
//! intrinsic, SE(3)-invariant geometric invariants:
//! - Discrete segment bond lengths: `l_i`
//! - Discrete tangent unit vectors: `T_i`
//! - Discrete curvature (turning angles): `kappa_i` in `[0, pi]`
//! - Discrete torsion (dihedral angle between consecutive osculating planes): `tau_i` in `(-pi, pi]`
//!
//! References:
//! - Hu, S., & Niemi, A. J. (2010). Discrete Frenet frames and knot theory of proteins.
//!   *Physical Review E*, 82(4), 041912.
//! - Bobenko, A. I., & Suris, Y. B. (2008). *Discrete Differential Geometry: Integrable Structure*.

#![forbid(unsafe_code)]

use nalgebra::{UnitVector3, Vector3};

use crate::error::GeometryError;
use crate::invariants::CurveInvariants;
use crate::types::BackboneTrace;

/// Numerical epsilon for singular configurations (degenerate lengths and collinear vectors).
pub const GEOMETRY_EPSILON: f64 = 1e-9;

/// Computes the unit tangent vectors `T_i` and Euclidean lengths for each segment
/// between consecutive C-alpha atoms.
///
/// For `i = 0 .. N-2`:
/// - `\Delta r_i = r_{i+1} - r_i`
/// - `l_i = ||\Delta r_i||`
/// - `T_i = \Delta r_i / l_i`
///
/// Output:
/// - Segment lengths count: `N - 1`
/// - Tangent vectors count: `N - 1`
///
/// # Errors
/// Returns [`GeometryError::InsufficientPoints`] if `trace.len() < 2`.
/// Returns [`GeometryError::DegenerateSegment`] if any consecutive distance is below [`GEOMETRY_EPSILON`].
pub fn compute_tangents(
    trace: &BackboneTrace,
) -> Result<(Vec<f64>, Vec<UnitVector3<f64>>), GeometryError> {
    let coords = trace.coordinates();
    let n = coords.len();

    if n < 2 {
        return Err(GeometryError::InsufficientPoints {
            required: 2,
            actual: n,
        });
    }

    let mut lengths = Vec::with_capacity(n - 1);
    let mut tangents = Vec::with_capacity(n - 1);

    for i in 0..n - 1 {
        let diff: Vector3<f64> = coords[i + 1] - coords[i];
        let len = diff.norm();

        if len < GEOMETRY_EPSILON {
            return Err(GeometryError::DegenerateSegment {
                index: i,
                next_index: i + 1,
                length: len,
            });
        }

        lengths.push(len);
        tangents.push(UnitVector3::new_unchecked(diff / len));
    }

    Ok((lengths, tangents))
}

/// Computes the discrete curvature `kappa_i` (turning angle in radians) at each internal vertex.
///
/// For internal vertex `i` (joining incoming segment `T_{i-1}` and outgoing segment `T_i`):
/// - `cos(kappa_i) = T_{i-1} . T_i`
/// - `kappa_i = acos(clamp(T_{i-1} . T_i, -1.0, 1.0))`
///
/// Output length: `N - 2` (where `N` is the number of atoms, and `tangents.len() == N - 1`).
///
/// # Errors
/// Returns [`GeometryError::InsufficientPoints`] if `tangents.len() < 2`.
pub fn compute_curvatures(tangents: &[UnitVector3<f64>]) -> Result<Vec<f64>, GeometryError> {
    if tangents.len() < 2 {
        return Err(GeometryError::InsufficientPoints {
            required: 2,
            actual: tangents.len(),
        });
    }

    let mut curvatures = Vec::with_capacity(tangents.len() - 1);

    for i in 0..tangents.len() - 1 {
        let dot = tangents[i].dot(&tangents[i + 1]).clamp(-1.0, 1.0);
        let kappa = dot.acos();
        curvatures.push(kappa);
    }

    Ok(curvatures)
}

/// Computes the unit binormal vectors `B_i` normal to each internal vertex osculating plane.
///
/// At vertex `i` (spanned by consecutive tangents `T_{i-1}` and `T_i`):
/// `B_i = (T_{i-1} x T_i) / ||T_{i-1} x T_i||`
///
/// Output length: `N - 2` (matching the number of curvature vertices).
///
/// # Errors
/// Returns [`GeometryError::InsufficientPoints`] if `tangents.len() < 2`.
/// Returns [`GeometryError::CollinearSegments`] if two consecutive tangents are parallel or antiparallel.
pub fn compute_binormals(
    tangents: &[UnitVector3<f64>],
) -> Result<Vec<UnitVector3<f64>>, GeometryError> {
    if tangents.len() < 2 {
        return Err(GeometryError::InsufficientPoints {
            required: 2,
            actual: tangents.len(),
        });
    }

    let mut binormals = Vec::with_capacity(tangents.len() - 1);

    for i in 0..tangents.len() - 1 {
        let cross = tangents[i].cross(&tangents[i + 1]);
        let cross_norm = cross.norm();

        if cross_norm < GEOMETRY_EPSILON {
            let dot = tangents[i].dot(&tangents[i + 1]).clamp(-1.0, 1.0);
            return Err(GeometryError::CollinearSegments {
                vertex_index: i + 1,
                turning_angle: dot.acos(),
            });
        }

        binormals.push(UnitVector3::new_unchecked(cross / cross_norm));
    }

    Ok(binormals)
}

/// Computes the discrete torsion `tau_i` (signed dihedral angle) between consecutive osculating planes.
///
/// The torsion between vertex `i` and vertex `i+1` is defined around the mutual segment axis `T_i`:
/// - `cos(tau_i) = B_i . B_{i+1}`
/// - `sin(tau_i) = (B_i x B_{i+1}) . T_i`
/// - `tau_i = atan2(sin(tau_i), cos(tau_i)) in (-pi, pi]`
///
/// Output length: `N - 3`.
///
/// # Errors
/// Returns [`GeometryError::InsufficientPoints`] if `tangents.len() < 3`.
/// Propagates [`GeometryError::CollinearSegments`] if osculating planes are ill-defined.
pub fn compute_torsions(tangents: &[UnitVector3<f64>]) -> Result<Vec<f64>, GeometryError> {
    if tangents.len() < 3 {
        return Err(GeometryError::InsufficientPoints {
            required: 3,
            actual: tangents.len(),
        });
    }

    let binormals = compute_binormals(tangents)?;
    let mut torsions = Vec::with_capacity(binormals.len() - 1);

    for i in 0..binormals.len() - 1 {
        let b_curr = &binormals[i];
        let b_next = &binormals[i + 1];
        let t_hinge = &tangents[i + 1]; // Mutual segment axis between vertex i and i+1

        let cos_tau = b_curr.dot(b_next).clamp(-1.0, 1.0);
        let sin_tau_vec = b_curr.cross(b_next);
        let sin_tau = sin_tau_vec.dot(t_hinge);

        let tau = sin_tau.atan2(cos_tau);
        torsions.push(tau);
    }

    Ok(torsions)
}

/// Extracts the complete SE(3)-invariant curve signature `(l_i, kappa_i, tau_i)` for a C-alpha backbone trace.
///
/// # Errors
/// Returns an error if the trace has fewer than 4 points or contains degenerate/collinear geometry.
pub fn extract_curve_invariants(trace: &BackboneTrace) -> Result<CurveInvariants, GeometryError> {
    if trace.len() < 4 {
        return Err(GeometryError::InsufficientPoints {
            required: 4,
            actual: trace.len(),
        });
    }

    let (segment_lengths, tangents) = compute_tangents(trace)?;
    let curvatures = compute_curvatures(&tangents)?;
    let torsions = compute_torsions(&tangents)?;

    Ok(CurveInvariants {
        segment_lengths,
        curvatures,
        torsions,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Isometry3, Point3, Translation3, UnitQuaternion};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use std::f64::consts::PI;

    /// Generates a mock C-alpha backbone resembling an ideal right-handed alpha-helix.
    ///
    /// Alpha-helix canonical parameters:
    /// - Radius: r approx 2.3 Å
    /// - Rise per residue: delta_z approx 1.5 Å
    /// - Angular pitch: delta_theta = 100 deg (approx 1.7453 rad)
    /// - Virtual C_alpha - C_alpha distance: approx 3.8 Å
    fn generate_mock_alpha_helix(n_residues: usize) -> BackboneTrace {
        let r = 2.3;
        let delta_z = 1.5;
        let delta_theta = 100.0 * PI / 180.0;

        let coords: Vec<Point3<f64>> = (0..n_residues)
            .map(|i| {
                let theta = (i as f64) * delta_theta;
                let z = (i as f64) * delta_z;
                Point3::new(r * theta.cos(), r * theta.sin(), z)
            })
            .collect();

        BackboneTrace::new(coords)
    }

    /// Computes angular distance on S^1 considering +/- pi branch cut wrapping.
    fn angular_difference(a: f64, b: f64) -> f64 {
        let mut diff = (a - b).abs();
        if diff > PI {
            diff = (2.0 * PI - diff).abs();
        }
        diff
    }

    #[test]
    fn test_se3_strict_invariance() {
        let n_atoms = 60;
        let original_trace = generate_mock_alpha_helix(n_atoms);
        let invariants_orig = extract_curve_invariants(&original_trace)
            .expect("Mock alpha-helix should compute invariants successfully");

        assert_eq!(invariants_orig.segment_lengths.len(), n_atoms - 1);
        assert_eq!(invariants_orig.curvatures.len(), n_atoms - 2);
        assert_eq!(invariants_orig.torsions.len(), n_atoms - 3);

        // Deterministic pseudo-random number generator for reproducible verification
        let mut rng = StdRng::seed_from_u64(0x42_DEAD_BEEF);
        let n_trials = 25;

        for trial in 0..n_trials {
            // Generate arbitrary random rotation R in SO(3) via random unit quaternion
            let q_rot = UnitQuaternion::from_euler_angles(
                rng.gen_range(-PI..PI),
                rng.gen_range(-PI..PI),
                rng.gen_range(-PI..PI),
            );

            // Generate arbitrary translation vector t in R^3 (large dynamic range)
            let translation = Translation3::new(
                rng.gen_range(-10_000.0..10_000.0),
                rng.gen_range(-10_000.0..10_000.0),
                rng.gen_range(-10_000.0..10_000.0),
            );

            // Construct rigid-body isometry g in SE(3)
            let isometry = Isometry3::from_parts(translation, q_rot);

            // Transform all C-alpha coordinates: r' = R * r + t
            let transformed_coords: Vec<Point3<f64>> = original_trace
                .coordinates()
                .iter()
                .map(|p| isometry.transform_point(p))
                .collect();

            let transformed_trace = BackboneTrace::new(transformed_coords);
            let invariants_trans = extract_curve_invariants(&transformed_trace)
                .expect("Transformed trace must produce valid invariants");

            // 1. Validate bond lengths invariance
            for i in 0..invariants_orig.segment_lengths.len() {
                let diff_l = (invariants_orig.segment_lengths[i] - invariants_trans.segment_lengths[i]).abs();
                assert!(
                    diff_l < 1e-6,
                    "Trial {trial}: Bond length invariance failed at segment {i}: orig={}, trans={}, diff={diff_l}",
                    invariants_orig.segment_lengths[i], invariants_trans.segment_lengths[i]
                );
            }

            // 2. Validate discrete curvature invariance
            for i in 0..invariants_orig.curvatures.len() {
                let diff_k = (invariants_orig.curvatures[i] - invariants_trans.curvatures[i]).abs();
                assert!(
                    diff_k < 1e-6,
                    "Trial {trial}: Curvature invariance failed at vertex {i}: orig={}, trans={}, diff={diff_k}",
                    invariants_orig.curvatures[i], invariants_trans.curvatures[i]
                );
            }

            // 3. Validate discrete torsion invariance
            for i in 0..invariants_orig.torsions.len() {
                let diff_tau = angular_difference(invariants_orig.torsions[i], invariants_trans.torsions[i]);
                assert!(
                    diff_tau < 1e-6,
                    "Trial {trial}: Torsion invariance failed at hinge {i}: orig={}, trans={}, diff={diff_tau}",
                    invariants_orig.torsions[i], invariants_trans.torsions[i]
                );
            }
        }
    }

    #[test]
    fn test_chiral_reflection_inversion() {
        // Under spatial reflection (improper rotation in O(3) with det = -1):
        // Curvature kappa must be preserved: kappa' = kappa
        // Torsion tau must flip sign: tau' = -tau
        let original_trace = generate_mock_alpha_helix(30);
        let orig = extract_curve_invariants(&original_trace).unwrap();

        // Parity reflection across xy-plane: (x, y, z) -> (x, y, -z)
        let reflected_coords: Vec<Point3<f64>> = original_trace
            .coordinates()
            .iter()
            .map(|p| Point3::new(p.x, p.y, -p.z))
            .collect();

        let reflected_trace = BackboneTrace::new(reflected_coords);
        let reflected = extract_curve_invariants(&reflected_trace).unwrap();

        for i in 0..orig.curvatures.len() {
            let diff_k = (orig.curvatures[i] - reflected.curvatures[i]).abs();
            assert!(diff_k < 1e-6, "Curvature must be parity-invariant");
        }

        for i in 0..orig.torsions.len() {
            let sum_tau = angular_difference(orig.torsions[i], -reflected.torsions[i]);
            assert!(sum_tau < 1e-6, "Torsion must flip sign under chiral reflection");
        }
    }

    #[test]
    fn test_degenerate_segment_detection() {
        let coords = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0), // Duplicate point
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ];
        let trace = BackboneTrace::new(coords);
        let result = extract_curve_invariants(&trace);

        match result {
            Err(GeometryError::DegenerateSegment { index, next_index, .. }) => {
                assert_eq!(index, 0);
                assert_eq!(next_index, 1);
            }
            other => panic!("Expected DegenerateSegment error, got {other:?}"),
        }
    }

    #[test]
    fn test_collinear_segments_detection() {
        let coords = vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0), // Collinear: 0, 1, 2 lie along x-axis
            Point3::new(2.0, 1.0, 0.0),
        ];
        let trace = BackboneTrace::new(coords);
        let result = extract_curve_invariants(&trace);

        match result {
            Err(GeometryError::CollinearSegments { vertex_index, .. }) => {
                assert_eq!(vertex_index, 1);
            }
            other => panic!("Expected CollinearSegments error, got {other:?}"),
        }
    }
}
