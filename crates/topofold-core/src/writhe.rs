//! Discrete Writhe and Gauss Linking Integrals for Space Curves.
//!
//! Writhe ($\text{Wr}$) is a topological and geometric invariant that quantifies
//! the degree of spatial non-planarity, 3D supercoiling, and helical handedness of a curve.
//! Unlike local curvature, writhe is sensitive to non-local self-wrapping.
//!
//! References:
//! - van Oosterom, A., & Strackee, J. (1983). The solid angle of a plane triangle.
//!   *IEEE Transactions on Biomedical Engineering*, BME-30(2), 125–126.
//! - Klenin, K., & Langowski, J. (2000). Computation of writhe in modeling of supercoiled DNA.
//!   *Biopolymers*, 54(5), 307–317.
//! - Røgen, P., & Fain, B. (2003). Automatic classification of protein structure by using Gauss integrals.
//!   *PNAS*, 100(1), 119–124.

#![forbid(unsafe_code)]

use nalgebra::{Point3, Vector3};
use std::f64::consts::PI;

use crate::error::GeometryError;
use crate::types::BackboneTrace;

/// Computes the signed solid angle subtended by a spherical triangle formed by three unit vectors
/// $\mathbf{a}, \mathbf{b}, \mathbf{c}$ from the origin.
///
/// Implements the unconditionally stable, branch-cut-free formula of van Oosterom & Strackee (1983):
/// $$\tan\left(\frac{1}{2}\Omega\right) = \frac{\mathbf{a} \cdot (\mathbf{b} \times \mathbf{c})}{1 + \mathbf{a}\cdot\mathbf{b} + \mathbf{b}\cdot\mathbf{c} + \mathbf{c}\cdot\mathbf{a}}$$
#[inline]
pub fn spherical_triangle_solid_angle(
    a: &Vector3<f64>,
    b: &Vector3<f64>,
    c: &Vector3<f64>,
) -> f64 {
    let det = a.dot(&b.cross(c));
    let denom = 1.0 + a.dot(b) + b.dot(c) + c.dot(a);
    2.0 * det.atan2(denom)
}

/// Computes the mutual writhe contribution between two directed straight line segments
/// $S_1 = [\mathbf{r}_1 \to \mathbf{r}_2]$ and $S_2 = [\mathbf{r}_3 \to \mathbf{r}_4]$.
///
/// Decomposes the spherical quadrangle formed by the direction vectors from $S_1$ to $S_2$
/// into two oriented spherical triangles, avoiding branch-cut discontinuities.
#[inline]
pub fn segment_pair_writhe(
    r1: &Point3<f64>,
    r2: &Point3<f64>,
    r3: &Point3<f64>,
    r4: &Point3<f64>,
) -> f64 {
    let r13: Vector3<f64> = r3 - r1;
    let r14: Vector3<f64> = r4 - r1;
    let r24: Vector3<f64> = r4 - r2;
    let r23: Vector3<f64> = r3 - r2;

    let norm_r13 = r13.norm();
    let norm_r14 = r14.norm();
    let norm_r24 = r24.norm();
    let norm_r23 = r23.norm();

    if norm_r13 < 1e-12 || norm_r14 < 1e-12 || norm_r24 < 1e-12 || norm_r23 < 1e-12 {
        return 0.0;
    }

    let u13 = r13 / norm_r13;
    let u14 = r14 / norm_r14;
    let u24 = r24 / norm_r24;
    let u23 = r23 / norm_r23;

    // Decompose spherical quadrangle (u13, u14, u24, u23) into two spherical triangles:
    // T1 = (u13, u14, u24), T2 = (u13, u24, u23)
    let omega1 = spherical_triangle_solid_angle(&u13, &u14, &u24);
    let omega2 = spherical_triangle_solid_angle(&u13, &u24, &u23);

    (omega1 + omega2) / (4.0 * PI)
}

/// Computes the total discrete writhe of a C-alpha backbone trace.
///
/// Evaluates the double sum over all pairs of non-adjacent segments:
/// $$\text{Wr} = \sum_{i=0}^{N-3} \sum_{j=i+2}^{N-2} w(i, j)$$
///
/// # Errors
/// Returns [`GeometryError::InsufficientPoints`] if `trace.len() < 4`.
pub fn compute_total_writhe(trace: &BackboneTrace) -> Result<f64, GeometryError> {
    let coords = trace.coordinates();
    let n = coords.len();

    if n < 4 {
        return Err(GeometryError::InsufficientPoints {
            required: 4,
            actual: n,
        });
    }

    let mut wr = 0.0;
    for i in 0..n - 2 {
        let p_i = &coords[i];
        let p_i1 = &coords[i + 1];
        for j in (i + 2)..n - 1 {
            let p_j = &coords[j];
            let p_j1 = &coords[j + 1];
            wr += segment_pair_writhe(p_i, p_i1, p_j, p_j1);
        }
    }

    Ok(wr)
}

/// Computes the localized writhe spectrum along the backbone using a sliding window.
///
/// For each internal residue $c \in [w, N - 1 - w]$:
/// Computes the discrete writhe of the subchain spanning $[c - w, c + w]$.
///
/// Output length: `N - 2 * window_radius`.
///
/// # Arguments
/// - `trace`: The C-alpha backbone trace.
/// - `window_radius`: Half-width of the window in residues (e.g. 5 for an 11-residue window).
///
/// # Errors
/// Returns [`GeometryError::InsufficientPoints`] if `trace.len() < 2 * window_radius + 2`.
pub fn compute_local_writhe(
    trace: &BackboneTrace,
    window_radius: usize,
) -> Result<Vec<f64>, GeometryError> {
    let coords = trace.coordinates();
    let n = coords.len();
    let min_required = 2 * window_radius + 2;

    if n < min_required {
        return Err(GeometryError::InsufficientPoints {
            required: min_required,
            actual: n,
        });
    }

    let mut spectrum = Vec::with_capacity(n - 2 * window_radius);

    for center in window_radius..(n - window_radius) {
        let start = center - window_radius;
        let end = center + window_radius; // inclusive

        let mut local_wr = 0.0;
        for i in start..end - 1 {
            let p_i = &coords[i];
            let p_i1 = &coords[i + 1];
            for j in (i + 2)..end {
                let p_j = &coords[j];
                let p_j1 = &coords[j + 1];
                local_wr += segment_pair_writhe(p_i, p_i1, p_j, p_j1);
            }
        }
        spectrum.push(local_wr);
    }

    Ok(spectrum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Isometry3, Translation3, UnitQuaternion};
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    fn generate_mock_helix(n_atoms: usize, right_handed: bool) -> BackboneTrace {
        let r = 2.3;
        let delta_z = 1.5;
        let sign = if right_handed { 1.0 } else { -1.0 };
        let delta_theta = sign * 100.0 * PI / 180.0;

        let coords: Vec<Point3<f64>> = (0..n_atoms)
            .map(|i| {
                let theta = (i as f64) * delta_theta;
                let z = (i as f64) * delta_z;
                Point3::new(r * theta.cos(), r * theta.sin(), z)
            })
            .collect();

        BackboneTrace::new(coords)
    }

    #[test]
    fn test_writhe_chirality_and_se3_invariance() {
        let n_atoms = 25;
        let right_helix = generate_mock_helix(n_atoms, true);
        let left_helix = generate_mock_helix(n_atoms, false);

        let wr_right = compute_total_writhe(&right_helix).expect("Writhe computation");
        let wr_left = compute_total_writhe(&left_helix).expect("Writhe computation");

        // Chirality check: opposite sign under mirror reflection
        assert!(
            (wr_right + wr_left).abs() < 1e-12,
            "Writhe must be chiral-antisymmetric: right={wr_right}, left={wr_left}"
        );
        assert!(
            wr_right.abs() > 0.1,
            "Helical curve must have significant non-zero writhe"
        );

        // SE(3) invariance check under random rotation and translation
        let mut rng = StdRng::seed_from_u64(0xFEED_FACE);
        for trial in 0..50 {
            let q = UnitQuaternion::from_euler_angles(
                rng.gen_range(-PI..PI),
                rng.gen_range(-PI..PI),
                rng.gen_range(-PI..PI),
            );
            let t = Translation3::new(
                rng.gen_range(-10_000.0..10_000.0),
                rng.gen_range(-10_000.0..10_000.0),
                rng.gen_range(-10_000.0..10_000.0),
            );
            let iso = Isometry3::from_parts(t, q);

            let trans_coords: Vec<Point3<f64>> = right_helix
                .coordinates()
                .iter()
                .map(|p| iso.transform_point(p))
                .collect();
            let trans_trace = BackboneTrace::new(trans_coords);

            let wr_trans = compute_total_writhe(&trans_trace).expect("Transformed writhe");
            let diff = (wr_right - wr_trans).abs();
            assert!(
                diff < 1e-6,
                "Trial {trial}: Writhe SE(3) invariance violation: diff={diff}"
            );
        }
    }

    #[test]
    fn test_local_writhe_spectrum() {
        let n_atoms = 30;
        let helix = generate_mock_helix(n_atoms, true);
        let window_radius = 4; // 9-residue window

        let spectrum = compute_local_writhe(&helix, window_radius).expect("Local writhe");
        assert_eq!(spectrum.len(), n_atoms - 2 * window_radius);

        // In a uniform helix, the local writhe in the interior should be constant
        let mean = spectrum.iter().sum::<f64>() / spectrum.len() as f64;
        for (i, &val) in spectrum.iter().enumerate() {
            assert!(
                (val - mean).abs() < 1e-4,
                "Interior local writhe in uniform helix should be constant: index={i}, val={val}, mean={mean}"
            );
        }
    }
}
